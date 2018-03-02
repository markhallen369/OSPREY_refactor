/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug.triage;

import cern.colt.matrix.DoubleFactory1D;
import edu.duke.cs.osprey.bbfree.BBFreeBlock;
import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.plug.FunnyVoxBuilder;
import edu.duke.cs.osprey.plug.LPChecks;
import edu.duke.cs.osprey.plug.RCTuplePolytope;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.NoFeasibleSolutionException;
import org.apache.commons.math3.optim.linear.Relationship;

/**
 *
 * Prune consecutive pairs of residues based on if there is a valid and
 * Ramachandran-allowed conf
 * In practice expect dihedrals to vary slow enough that need not search comprehensively
 * DEBUG!!! might check this
 * 
 * @author mhall44
 */
public class RamaPruningStep extends GeomPruningStep {
    
    ArrayList<LinearConstraint> BFBBoxConstr = null;//main box constraints on BFB
    
    HashMap<Residue,ArrayList<ArrayList<LinearConstraint> > > bbVoxConstr = new HashMap<>();
    //single-voxel constraints for each residue and bbVox
    
    
    RamaPruningStep(ArrayList<RCClassification> classesAtRes, PruningLookup prooningLookup, 
            RCClass.ClassType classType, RCClass.ClassType classType2, ArrayList<Residue> shellResidues, 
            CentralStericCache stericCache, ConfSpace confSpace) {
        super(classesAtRes, prooningLookup, classType, classType2, shellResidues, stericCache);
        
        BBFreeBlock bfb = FunnyVoxBuilder.getBBFreeBlock(confSpace);
        
        for(int pos=0; pos<confSpace.numPos; pos++){
            PositionConfSpace pcSpace = confSpace.posFlex.get(pos);
            if(pcSpace.RCs.get(0).bbVoxNum!=-1){//there is FVB flexibility
                
                //collect BFBBoxConstr if haven't yet
                if(BFBBoxConstr==null)
                    collectBFBBoxConstr(pcSpace.RCs.get(0), bfb.getDOFs().size());
                
                ArrayList<ArrayList<LinearConstraint> > resBBVoxConstr = new ArrayList<>();
                for(int rcNum=0; rcNum<pcSpace.RCs.size(); rcNum++){
                    RC rc = pcSpace.RCs.get(rcNum);
                    if(rcNum>0 && rc.bbVoxNum==0)//all bb vox's have been collected
                        break;
                    if(rcNum!=rc.bbVoxNum)
                        throw new RuntimeException("ERROR: expected FVB ordering");
                    
                    ArrayList<LinearConstraint> voxConstr = RCTuplePolytope.transferConstraints(
                            rc.DOFs, bfb.getDOFs(), rc.linConstr);
                    resBBVoxConstr.add(voxConstr);
                }
                
                bbVoxConstr.put(pcSpace.res, resBBVoxConstr);
            }
        }
    }

    
    private void collectBFBBoxConstr(RC rc, int numBBDOFs){
        //collect the main box constr on the bfb based on this rc
        BFBBoxConstr = new ArrayList<>();
        int bbDOFCount = 0;
        for(int dof=0; dof<rc.DOFs.size(); dof++){
            if(rc.DOFs.get(dof) instanceof BBFreeDOF){
                double[] lbCoeffs = new double[numBBDOFs];
                lbCoeffs[bbDOFCount] = 1;
                LinearConstraint lb = new LinearConstraint(lbCoeffs, Relationship.GEQ, rc.DOFmin.get(dof));
                BFBBoxConstr.add(lb);
                
                double[] ubCoeffs = new double[numBBDOFs];
                ubCoeffs[bbDOFCount] = 1;
                LinearConstraint ub = new LinearConstraint(ubCoeffs, Relationship.LEQ, rc.DOFmax.get(dof));
                BFBBoxConstr.add(ub);
                bbDOFCount++;
            }
        }
        
        if(BFBBoxConstr.size()!=2*numBBDOFs)
            throw new RuntimeException("ERROR wrong number of BBFreeDOFs");
    }
    
    
            
    private ArrayList<LinearConstraint> bbVoxConstr(RCClass cl){
        if(cl.origBBEquivalent==null){
            if(bbVoxConstr.containsKey(cl.res))//there are constr for this res
                return bbVoxConstr.get(cl.res).get(cl.bbVoxNum);
            else
                return new ArrayList<>();
        }
        else
            return bbVoxConstr.get(cl.origBBEquivalent.res).get(cl.bbVoxNum);
    }
    
    
    @Override
    boolean canPrune(RCClassTuple tup) {
        if(tup.size()!=2)
            throw new RuntimeException("ERROR: Rama pruning is for pairs");
        
        //we can actually do the feasibility check for non-consecutive pairs,
        //but only consecutive pairs can really have a Ramachandran check
        RCClass cl1 = tup.get(0);
        RCClass cl2 = tup.get(1);
        
        //OK we only really care about BB vox nums
        ArrayList<LinearConstraint> constr = new ArrayList(BFBBoxConstr);
        ArrayList<LinearConstraint> constr1 = bbVoxConstr(cl1);
        ArrayList<LinearConstraint> constr2 = bbVoxConstr(cl2);
        if(constr1.isEmpty() || constr2.isEmpty())//no bb flex for one or other, won't prune
            return false;
        
        constr.addAll( constr1 );
        constr.addAll( constr2 );
        
        //this could be slow but let's see
        double[] feasiblePt;//feasible values for BFB free DOFs
        try {
            feasiblePt = LPChecks.getFeasiblePt(constr);
        }
        catch(NoFeasibleSolutionException e){//not feasible
            return true;
        }
        
        return false;//DEBUG!!!  ignoring rama here, there are vox variability concerns
        //I think Rama is best posed as a series of lin constr in free DOF space.  
        //That u can combine with these other constr if u want 2.  (Even 1-wise!!)
        //bfb.setDOFs(DoubleFactory1D.dense.make(feasiblePt));
        //double CCoords[] = res1.getCoordsByAtomName("N");
        
        //SO actually true rama is a three-res thing
        //but we can do a two-res equivalent where we see if the given phi and psi both make
        //sense on their own...
        //this could still get a lot of the pruning (sort of creates separable envelope of Rama plot)
        //hmm also variability across voxel could be a significant issue indeed...
    }
    
}
