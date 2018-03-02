/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.aibis;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import com.joptimizer.optimizers.JOptimizer;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.ematrix.epic.SAPE;
import edu.duke.cs.osprey.ematrix.epic.SeriesFitter;
import static edu.duke.cs.osprey.ibis.IntegrableDOF.RT;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.plug.LPChecks;
import static edu.duke.cs.osprey.plug.LPChecks.constrAsFunction;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.NoFeasibleSolutionException;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.linear.UnboundedSolutionException;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;

/**
 *
 * An integral DOF that is one of the continuous conformational degrees of freedom
 * 
 * @author mhall44
 */
public class ContIntegrableDOF implements IntegrableDOF {
    
    DegreeOfFreedom dof;
    PolytopeMatrix plugMtx;
    int dofIndex;
    
    public ContIntegrableDOF(DegreeOfFreedom dof, PolytopeMatrix plugMat, int dofIndex){
        this.dof = dof;
        this.dofIndex = dofIndex;
        plugMtx = plugMat;
    }
    
    
    
    DofPLUGBounds nonClashingBounds(ConfSample samp, SparseLinearEnergy baseE){
        ArrayList<DegreeOfFreedom> fixedDOFs = new ArrayList<>();
        for(DegreeOfFreedom confDOF : baseE.getContDOFs()){
            if( ! confDOF.getName().equalsIgnoreCase(dof.getName()) ){
                fixedDOFs.add(confDOF);
            }
        }
        
        return nonClashingBounds(samp, fixedDOFs);
    }
    
    
    @Override
    public double boltzmannIntegratedEnergy(SparseLinearEnergy baseE, ConfSample samp) {
        DofPLUGBounds dpb = nonClashingBounds(samp, baseE);
        switch(dpb.status){
            case CLASHING:
                System.out.println("Warning: Could not find non-clashing bounds; sample w/i numerical error of clashing");
                return 0;
            case INAPPLICABLE://no need to integrate
                double E = baseE.evalEnergy(samp);
                if(Double.isInfinite(E)){
                    int aaa = 5;
                }
                return E;
            case BOUNDED:
                double[] bounds = dpb.bounds;
                
                
                //DEBUG!!!!!!!
                //bounds = calcFullInterval(samp.vox);
                
                
                MoleculeModifierAndScorer mms = baseE.makeMMS(samp);//same voxel will apply across the range of integration
                int numIntegBins = 20;//DEBUG!!! will ultimately probably want 
                //some kind of adaptive integration, but for Simpson's just need a large even # of bins

                double q=0;
                int wtSum = 0;
                //let's use Simpson's rule.  
                //Normalization: we just average over the interval (width of voxel not non-clashing region).  Shifts due to normalization
                //cancel between same-seq voxels anyway as long as vox widths consistent between confs.  
                
                //System.out.println("INTEGRATING");
                
                double minE = Double.POSITIVE_INFINITY;
                
                for(int step=0; step<=numIntegBins; step++){
                    double dofVal = bounds[0] + (bounds[1]-bounds[0])*step/numIntegBins;
                    
                    //HashMap<String,Double> altDOFVals = new HashMap<>();
                    //altDOFVals.putAll(samp.contDOFVals);
                    //altDOFVals.put(dof.getName(), dofVal);
                    double[] altDOFValsArr = samp.contDOFValsArr.clone();
                    altDOFValsArr[dofIndex] = dofVal;
                    
                    ConfSample altSamp = new ConfSample(samp.vox, altDOFValsArr);
                    int wt = simpsonCoeff(step,numIntegBins);
                    double stepE = baseE.evalEnergy(altSamp,mms);
                    q += wt * Math.exp(-stepE/RT);
                    wtSum += wt;
                    
                    
                    minE = Math.min(minE, stepE);
                    //System.out.println(dofVal+" "+stepE);
                }

                
                //if(true)
                //    return minE;//DEBUG!!!!
                
                q /= wtSum;

                //currently normalized by width of non-clashing interval, get real interval
                double fullIntervalWidth = calcFullIntervalWidth(samp.vox);
                if(Math.abs(fullIntervalWidth)<1e-10)
                    throw new RuntimeException("ERROR: Trying to integrate wrt discrete DOF");
                q *= (bounds[1]-bounds[0]) / fullIntervalWidth;
                //DEBUG!!
                
                
                double ans = -RT*Math.log(q);
                if(Double.isInfinite(ans)){
                    throw new RuntimeException("ERROR: Infinite free energy");//can be caused by overflow of exp
                }
                
                //System.out.println("INTEGRATED.  FULL INTERVAL WIDTH: "+fullIntervalWidth+" G: "+ans);

                
                return ans;
            default:
                throw new RuntimeException("ERROR: Unexpected status");
        }
    }
    
    
    public void scanIntegratedEnergies(SparseLinearEnergy baseE, ConfSample samp, int scanDOFIndex){
        //scan integrated energies wrt an additional dof (numbered among all DOFs)
        System.out.println("SCANNING INTEGRATED E.  x, reg E, integ E, plugwidth, avg E");
        
        //figure out voxel bounds on scanDOF
        DegreeOfFreedom scanDOF = plugMtx.cSpace.confDOFs.get(scanDOFIndex);
        double scanStart=0, scanEnd=0;
        for(int posCount=0; posCount<samp.vox.size(); posCount++){
            RC rc = plugMtx.cSpace.posFlex.get(samp.vox.pos.get(posCount)).RCs.get(samp.vox.RCs.get(posCount));
            for(int dofIndex=0; dofIndex<rc.DOFs.size(); dofIndex++){
                if(rc.DOFs.get(dofIndex).getName().equalsIgnoreCase(scanDOF.getName())){//found bounds on our dof
                    scanEnd = rc.DOFmax.get(dofIndex);
                    scanStart = rc.DOFmin.get(dofIndex);
                    break;
                }
            }
        }
        
        ArrayList<LinearConstraint> polytope = plugMtx.getFullStericPolytope(samp.vox, plugMtx.cSpace.confDOFs);
        
        //now scan
        int numBins = 20;
        for(int step=0; step<=numBins; step++){
            double dofVal = scanStart + (scanEnd-scanStart)*step/numBins;

            double[] altDOFValsArr = samp.contDOFValsArr.clone();
            altDOFValsArr[scanDOFIndex] = dofVal;
            if(LPChecks.isPointInPolytope(polytope, altDOFValsArr)){
                ConfSample altSamp = new ConfSample(samp.vox, altDOFValsArr);
                double regE = baseE.evalEnergy(altSamp);
                double integE = boltzmannIntegratedEnergy(baseE,altSamp);
                DofPLUGBounds dpb = nonClashingBounds(altSamp,baseE);
                double plugWidth = dpb.status==DPBStatus.BOUNDED ? (dpb.bounds[1]-dpb.bounds[0]) : 0;
                double fullWidth = calcFullIntervalWidth(samp.vox);
                double avgE = integE + RT*Math.log(plugWidth/fullWidth);
                System.out.println(dofVal+", "+regE+", "+integE+", "+plugWidth+", "+avgE);
            }
            else {
                System.out.println(dofVal+", X");
            }
        }
        
        
        System.out.println("DONE SCANNING INTEGRATED E");
    }
    
    
    double[] calcFullInterval(RCTuple vox){
        for(int posCount=0; posCount<vox.size(); posCount++){
            RC rc = plugMtx.cSpace.posFlex.get(vox.pos.get(posCount)).RCs.get(vox.RCs.get(posCount));
            for(int dofIndex=0; dofIndex<rc.DOFs.size(); dofIndex++){
                if(rc.DOFs.get(dofIndex).getName().equalsIgnoreCase(dof.getName())){//found bounds on our dof
                    return new double[] {rc.DOFmin.get(dofIndex), rc.DOFmax.get(dofIndex)};
                }
            }
        }
        throw new RuntimeException("ERROR: DOF "+dof.getName()+" not defined for voxel: "+vox.stringListing());
    }
    
    double calcFullIntervalWidth(RCTuple vox){
        for(int posCount=0; posCount<vox.size(); posCount++){
            RC rc = plugMtx.cSpace.posFlex.get(vox.pos.get(posCount)).RCs.get(vox.RCs.get(posCount));
            for(int dofIndex=0; dofIndex<rc.DOFs.size(); dofIndex++){
                if(rc.DOFs.get(dofIndex).getName().equalsIgnoreCase(dof.getName())){//found bounds on our dof
                    return rc.DOFmax.get(dofIndex) - rc.DOFmin.get(dofIndex);
                }
            }
        }
        throw new RuntimeException("ERROR: DOF "+dof.getName()+" not defined for voxel: "+vox.stringListing());
    }
    
    static private int simpsonCoeff(int step, int numBins){
        //For Simpson integration with numBins bins (and thus numBins+1 func evals),
        //get the coefficient for the step'th func val
        if(step==0 || step==numBins)
            return 1;
        else return step%2==0 ? 2 : 4;
    }
    
    
    private class DofPLUGBounds {
        double[] bounds;
        DPBStatus status;
        
        private DofPLUGBounds(double[] bounds, DPBStatus status) {
            this.bounds = bounds;
            this.status = status;
        }
    }
    
    enum DPBStatus {
        BOUNDED,//PLUG imposes non-empty bounds on the DOF
        CLASHING,///no bounds (could arise if sample is in a corner of the voxel,
            //and after equality constraints that define the other DOFs in an inequality, sample
            //is left outside teh voxel within numerical error)
        INAPPLICABLE//DOF doesn't apply to this voxel
    };

    DofPLUGBounds nonClashingBounds(ConfSample samp, Iterable<DegreeOfFreedom> fixedDOFs){
        //find the bounds on this DOF indicated by samp, with degrees of freedom fixed as indicated in baseE
        
        //DEBUG!!  this relies on the thing where "steric" polytope actually contains voxel bounds (replicated excessively)
        ArrayList<DegreeOfFreedom> allDOFs = plugMtx.cSpace.confDOFs;
        ArrayList<LinearConstraint> constrList = plugMtx.getFullStericPolytope(samp.vox, allDOFs);
        
        ArrayList<String> allDOFNames = new ArrayList<>();
        for(DegreeOfFreedom curDOF : allDOFs)
            allDOFNames.add(curDOF.getName());
        
        //now impose the add'l equality constraints implied by samp
        for(DegreeOfFreedom fixedDOF : fixedDOFs){
            double[] q = new double[allDOFs.size()];
            int fixedDOFIndex = allDOFNames.indexOf(fixedDOF.getName());
            q[fixedDOFIndex] = 1;
            LinearConstraint eqConstr = new LinearConstraint( q, Relationship.EQ, samp.contDOFValsArr[fixedDOFIndex] );//samp.contDOFVals.get(fixedDOF.getName()) );
            constrList.add(eqConstr);
        }
        
        int dofIndex = allDOFNames.indexOf(dof.getName());
        boolean dofActive = false;//dof constrained and thus active in this conformation
        for(LinearConstraint constr : constrList){
            if( Math.abs(constr.getCoefficients().getEntry(dofIndex)) > 1e-10 ){
                dofActive = true;
                break;
            }
        }
        
        if(dofActive){
            //finally, maximize and minimize this DOF by simplex w/i the constraints
            double[] bounds = boundDOFInPolytope(constrList, dofIndex, false);
            if(bounds==null)
                return new DofPLUGBounds(bounds,DPBStatus.CLASHING);
            return new DofPLUGBounds(bounds,DPBStatus.BOUNDED);
        }
        else {
            //dof doesn't apply to this voxel
            return new DofPLUGBounds(null,DPBStatus.INAPPLICABLE);
        }
    }
    
    static double[] boundDOFInPolytope(ArrayList<LinearConstraint> polytope, int dofNum, boolean errorOnInfeasible){
        double[] q = new double[polytope.get(0).getCoefficients().getDimension()];
        q[dofNum] = 1;
        LinearObjectiveFunction objFcn = new LinearObjectiveFunction(q,0);
        SimplexSolver ss = new SimplexSolver();
        LinearConstraintSet constrSet = new LinearConstraintSet(polytope);
        double bounds[] = new double[2];
        
        try{
            bounds[0] = ss.optimize(objFcn, constrSet, GoalType.MINIMIZE).getValue();
            bounds[1] = ss.optimize(objFcn, constrSet, GoalType.MAXIMIZE).getValue();
        }
        catch(NoFeasibleSolutionException e){
            if(errorOnInfeasible)
                throw e;
            else
                return null;//there are no interior points
        }
        /*catch(UnboundedSolutionException e){
            //DEBUG!!!
            //and probs reduce # of constraints!!  redundant!!!
            for(int d=0; d<polytope.get(0).getCoefficients().getDimension(); d++){
                double[] dub = new double[polytope.get(0).getCoefficients().getDimension()];
                dub[d] = 1;
                polytope.add(new LinearConstraint(dub,Relationship.GEQ,-1000000000.));
                polytope.add(new LinearConstraint(dub,Relationship.LEQ,1000000000.));
            }
            constrSet = new LinearConstraintSet(polytope);
            bounds[0] = ss.optimize(objFcn, constrSet, GoalType.MINIMIZE).getValue();
            bounds[1] = ss.optimize(objFcn, constrSet, GoalType.MAXIMIZE).getValue();
            int trex = 9;
        }*/
        
        return bounds;
    }


    @Override
    public FeatureSet featureSetForInteg(FeatureSet baseFeat) {
        //basically go through all the DenseFeatureSet stripping out polynomial & SAPE terms that depend on this DOF
        //ultimately this will strip out everything but the constant 
        //(and backbone SAPEs that are effectively constant--so remove SAPE when poly down to const)
        TupleMatrixGeneric<DenseFeatureSet> oldMtx = baseFeat.featMatrix;
        FeatureSet ans = new FeatureSet(baseFeat.confSpace);
        
        int checkOldNumFeat = 0;
        
        for(RCTuple tup : baseFeat.unprunedFeatTuples()){
            DenseFeatureSet oldFS = oldMtx.getTupleValue(tup);
            ans.addDenseFeatureSet(denseFeatureSetForInteg(oldFS), tup);
            checkOldNumFeat += oldFS.getNumFeatures();
            if(baseFeat.freezeSAPE && baseFeat.extractSAPE(oldFS)!=null)
                checkOldNumFeat--;
        }
        
        if(checkOldNumFeat != baseFeat.numFeatures)
            throw new RuntimeException("ERROR: Unexpected number of features in unintegrated FeatureSet--might be param sharing (unsupported here)");
        //currently not supporting param sharing, e.g. integrating over RCs, 
        //and then going back to integrating over more cont DOFs.  May not ever need to...other order seems better
        
        return ans;
    }
    
    
    private DenseFeatureSet denseFeatureSetForInteg(DenseFeatureSet oldFS){
        if(oldFS instanceof DenseFeatureSet.Const)
            return oldFS;
        else if(oldFS instanceof DenseFeatureSet.EPIC){
            DenseFeatureSet.EPIC efs = (DenseFeatureSet.EPIC)oldFS;
            
            ArrayList<DegreeOfFreedom> redDOFs = new ArrayList<>();//all oldFS' dofs expect this one
            ArrayList<Double> centerDOFVals = new ArrayList<>();//center of reduced EPIC term.  Arbitrary but could help numerically
            ArrayList<Integer> redDOFIndexList = new ArrayList<>();
            for(int oldDOFNum=0; oldDOFNum<efs.dofs.size(); oldDOFNum++){
                if(!efs.dofs.get(oldDOFNum).getName().equalsIgnoreCase(dof.getName())){
                    redDOFs.add(efs.dofs.get(oldDOFNum));
                    centerDOFVals.add(efs.center.get(oldDOFNum));
                    redDOFIndexList.add(efs.dofIndices[oldDOFNum]);
                }
            }
            
            if(redDOFs.isEmpty())//no DOFs left: we now have just a constant
                return new DenseFeatureSet.Const();
            
            SAPE redSAPE = null;
            if(efs.sapeTerm!=null)
                redSAPE = sapeForDOFs(efs.sapeTerm, redDOFs, efs.moveableResTypes);
            
            
            int[] redDOFIndices = new int[redDOFIndexList.size()];
            for(int q=0; q<redDOFIndexList.size(); q++)
                redDOFIndices[q] = redDOFIndexList.get(q);
            
            return new DenseFeatureSet.EPIC(redSAPE, efs.order, redDOFs, centerDOFVals, efs.moveableResTypes, redDOFIndices);
        }
        else
            throw new RuntimeException("ERROR: Unexpected DenseFeatureSet type");
    }
    
    //DEBUG!!!  Do high-RT comparison to minim
    
    SAPE sapeForDOFs(SAPE oldSAPE, ArrayList<DegreeOfFreedom> dofs, HashMap<Integer,ResidueTemplate> moveableResTypes){
        //grab the terms in oldSAPE that depend on dofs but not on any other DOFs in the system
        //hence these terms are fully defined by dofs
        //and should include all the close-range interactions that are fully defined by dofs
        HashSet<String> dofNames = new HashSet<>();
        for(DegreeOfFreedom curDOF : dofs)
            dofNames.add(curDOF.getName());
        
        HashMap<Integer,HashSet<Integer> > atomsToInclude = new HashMap<>();
        for(int resIndex : moveableResTypes.keySet()){
            HashSet<Integer> resAtomsToInclude = new HashSet<>();
            for(DegreeOfFreedom includedDOF : dofs)//get all the atoms depending on dofs
                resAtomsToInclude.addAll(dofAffectedAtoms(includedDOF,resIndex,moveableResTypes.get(resIndex)));
            for(DegreeOfFreedom excludedDOF : plugMtx.cSpace.confDOFs){//remove all the atoms depending on DOFs not in dofs
                if(!dofNames.contains(excludedDOF.getName())){//check it's really an excluded DOF
                    resAtomsToInclude.removeAll(dofAffectedAtoms(excludedDOF,resIndex,moveableResTypes.get(resIndex)));
                }
            }
            atomsToInclude.put(resIndex, resAtomsToInclude);
        }
        return new SAPE(oldSAPE, atomsToInclude);
    }
    
    static ArrayList<Integer> dofAffectedAtoms(DegreeOfFreedom dof, int resIndex, ResidueTemplate resType){
        //specify indices (among atoms in the specified template) of the atoms that dof can move
        //in the residue with the specified index in the molecule
        if(dof instanceof FreeDihedral){
            if(((FreeDihedral)dof).getResidue().indexInMolecule==resIndex){
                int dihNum = ((FreeDihedral)dof).getDihedralNumber();
                if(dihNum < resType.numDihedrals)//This dihedral can move this residue type
                    return resType.dihedralMovingAtoms.get(dihNum);
            }
            
            //dof doesn't apply to residue
            return new ArrayList<>();
        }
        else{//DEBUG!!!  this maintains correctness because we make sure not to make SAPEs depend on undefined DOFs 
            //Good for CATS & strand motions, and most DEEPer stuff
            //However if atoms can move but not with this DOF, e.g. stacked perturbations, SAPE could miss some useful interactions
            boolean resMoving = false;
            for(Residue movingRes : dof.getBlock().listResidues()){
                if(movingRes.indexInMolecule==resIndex)
                    resMoving = true;
            }
            ArrayList<Integer> ans = new ArrayList<>();
            if(resMoving){
                for(int a=0; a<resType.templateRes.atoms.size(); a++)
                    ans.add(a);
            }
            return ans;
        }
    }

    @Override
    public String description() {
        return "ContIntegrableDOF: "+dof.getName();
    }
    
}
