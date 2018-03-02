/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tupexp;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.epic.EPICEnergyFunction;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.tools.ObjectIO;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Implementation of a tuple expander that expands conformational energies--
 * including continuous flexibility, non-pairwise terms, etc.
 * 
 * @author mhall44
 */
public class ConfETupleExpander extends TupleExpander {
    
    SearchProblem sp;
    
    boolean fullPLUGPruning;//only allow sampled conformations that PLUG considers feasible
    
    public ConfETupleExpander(SearchProblem sp, LUTECompStatus compStatus){
        super(sp.confSpace.numPos, sp.confSpace.getNumRCsAtPos(), sp.pruneMat.getPruningInterval(), sp.luteSettings, compStatus);
        this.sp = sp;
        fullPLUGPruning = sp.plugMat!=null;//DEBUG!!!
        if(fullPLUGPruning)
            canCheckPartialPruning = false;
    }

    double worstELBDiff = 0;
        
    @Override
    public double scoreAssignmentList(int[] assignmentList) {
        if(sp.useEPIC){//Faster if we can score by EPIC
            
            double E = sp.EPICMinimizedEnergy(assignmentList);
           
            RCTuple tup = new RCTuple(assignmentList);
            if(E==Double.POSITIVE_INFINITY){//this is going to be a problem if used as a true value
                if(isPruned(tup))
                    throw new RuntimeException("ERROR: Scoring pruned conformation: "+tup.stringListing());
                else if(!TESampleSet.rejectInf)//this might happen if we're doing rejectInf
                    throw new RuntimeException("ERROR: Infinite E for unpruned conf: "+tup.stringListing());
            }
            
            if(fullPLUGPruning){
                if( ! sp.plugMat.isTupleFeasible(tup) ){
                    throw new RuntimeException("ERROR: Scoring PLUG-infeasible conformation: "+tup.stringListing());
                }
            }
            
            
            /*
            //code for debugging inaccurate pairwise lower bounds
            double LB = sp.lowerBound(assignmentList);
            if(E - LB  < -10){
                System.out.println("ERROR: E="+E+" LB="+LB);
                System.out.print( "Assignments: ");
                for(int a : assignmentList)
                    System.out.print(a+", ");
                
                RCTuple tup = new RCTuple(assignmentList);
                EPICEnergyFunction efunc = sp.epicMat.internalEnergyFunction(tup);
                MoleculeModifierAndScorer objFcn = new MoleculeModifierAndScorer(efunc,sp.epicMat.getConfSpace(),tup);
                
                Minimizer minim = new CCDMinimizer(objFcn,false);
                DoubleMatrix1D bestDOFVals = minim.minimize();
                System.out.println("Best DOF vals check: "+bestDOFVals);
                double Echeck = objFcn.getValue(bestDOFVals);
                System.out.println("E check: "+Echeck);
                
                System.out.println("EPIC term values: ");
                efunc.printAllTermValues();
                System.out.println("End EPIC term values");
                
                if(E-LB < worstELBDiff){
                    worstELBDiff = E-LB;
                    System.out.println("Outputting obj fcn to problem_obj_fcn"+worstELBDiff+".dat");
                    ObjectIO.writeObject(objFcn, "problem_obj_fcn"+worstELBDiff+".dat");
                }
            }*/
            
            return E;
        }
        else//Score by minimizing energy function directly
            return sp.minimizedEnergy(assignmentList);
    }

    
    @Override
    public boolean isPruned(RCTuple tup) {
        if (sp.pruneMat.isPruned(tup))
            return true;
        if(fullPLUGPruning){//pruning without marking--we only need expansion to be right on PLUG-feasible voxels
            if(!sp.plugMat.isTupleFeasible(tup))
                return true;
        }
        return false;
    }

    
    @Override
    public void pruneTuple(RCTuple tup) {
        sp.pruneMat.markAsPruned(tup);
    }

    
    @Override
    public ArrayList<RCTuple> higherOrderPrunedTuples(RCTuple tup) {
        //list higher-order pruned tuples containing the pair tup
        
        if(tup.pos.size() != 2)
            throw new RuntimeException("ERROR: higherOrderPrunedTuples is meant to take an RC pair as argument");
                
        HigherTupleFinder<Boolean> htf = 
                sp.pruneMat.getHigherOrderTerms(tup.pos.get(0), tup.RCs.get(0), tup.pos.get(1), tup.RCs.get(1));
        
        if(htf!=null){
            ArrayList<RCTuple> otherTups = htf.listInteractionsWithValue(true);
            //otherTups are recorded as what tup interacts with...add in tup to get the whole pruned tuple
            for(RCTuple otherTup: otherTups){
                for(int i=0; i<2; i++){
                    otherTup.pos.add(tup.pos.get(i));
                    otherTup.RCs.add(tup.RCs.get(i));
                }
            }
            
            return otherTups;
        }
        else//no interactions
            return new ArrayList<>();
    }

}
