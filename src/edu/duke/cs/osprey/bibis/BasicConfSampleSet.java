/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.bibis;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.plug.LPChecks;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tupexp.BasicPruningTupleExpander;
import edu.duke.cs.osprey.tupexp.TESampleSet;
import edu.duke.cs.osprey.tupexp.TupleExpander;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.Relationship;

/**
 *
 * Just a list of samples and (optionally) their energies
 * (can be used w/ or w/o the energies--w/o=regular ConfSampleSet)
 * 
 * @author mhall44
 */
public class BasicConfSampleSet implements EnergiedConfSampleSet {
    
    ArrayList<ConfSample> samples;
    ArrayList<Double> energies = null;//can be null if using this as a ConfSampleSet


    
    
    protected static class VoxDrawTupleExpander extends BasicPruningTupleExpander {//for drawing samples

        FeatureSet featSet;
        int numSampPerFeach = 30;//DEBUG!!!
        
        VoxDrawTupleExpander(PruningMatrix pruneMat, PolytopeMatrix plugMat, FeatureSet featSet){
            super(pruneMat,plugMat);
            this.featSet = featSet;
        }
        
        @Override
        public int numSamplesNeeded(int tup){
            //for a tuple (index in tuples), how many samples are needed?
            return numSampPerFeach*featSet.featMatrix.getTupleValue(getTuples().get(tup)).getNumFeatures();
        }
    }
    
    public BasicConfSampleSet(FeatureSet featSet, PolytopeMatrix plugMtx, PruningMatrix pruneMat){
        //The main issue in sampling is since we have a lot of terms with small support,
        //make sure that each RC tuple in featSet is represented by enough samples
        //and also that the total number of samples is > number of params
        //DEBUG!!! might be convenient to test if this really matters much
        //Anyway constructor will just draw samples, energies can be handled separately
        ArrayList<int[]> voxSamples = drawVoxSamples(featSet, pruneMat, plugMtx);
        samples = new ArrayList<>();
        ArrayList<DegreeOfFreedom> DOFs = plugMtx.cSpace.confDOFs;
        for(int[] voxSamp : voxSamples){
            RCTuple tup = new RCTuple(voxSamp);
            ArrayList<LinearConstraint> tope = plugMtx.getFullStericPolytope(tup,DOFs);//DEBUG!! again, uses funny def'n of steric
            ArrayList<Double> dofVals = sampleFeasPt(tope);
            if(dofVals.isEmpty())//no constraints or DOFs...just add all 0s
                dofVals.addAll(Collections.nCopies(DOFs.size(),0.));
            
            HashMap<String,Double> contDOFVals = new HashMap<>();
            for(int dofNum=0; dofNum<DOFs.size(); dofNum++)
                contDOFVals.put(DOFs.get(dofNum).getName(), dofVals.get(dofNum));
            samples.add(new ConfSample(tup,contDOFVals));
        }
    }
    
    
    public static ArrayList<Double> sampleFeasPt(ArrayList<LinearConstraint> polytope){
        double[] pt = LPChecks.getFeasiblePt(polytope);//.getInteriorPt(polytope);
        RealVector vec = new ArrayRealVector(pt);//.toArray());
        //Let's go through a round or 2 of Gibbs sampling within the polytope
        int gibbsRounds = 2;
        for(int round=0; round<gibbsRounds; round++){
            for(int d=0; d<vec.getDimension(); d++){
                //figure out the bounds on this dimension
                double[] dofBounds = boundSingleDOF(vec, polytope, d);
                double drawVal = dofBounds[0] + Math.random()*(dofBounds[1]-dofBounds[0]);
                vec.setEntry(d, drawVal);
            }
        }
        
        ArrayList<Double> ans = new ArrayList<>();
        for(int d=0; d<vec.getDimension(); d++)
            ans.add(vec.getEntry(d));
        return ans;
    }
    
    static double[] boundSingleDOF(RealVector vec, ArrayList<LinearConstraint> polytope, int dofNum){
        //If we want to change just DOF #dofNum in vec, what bounds do
        //the constraints in polytope impose on it?
        double maxVal = Double.POSITIVE_INFINITY;
        double minVal = Double.NEGATIVE_INFINITY;
        
        for(LinearConstraint constr : polytope){
            //boundary is at q*x+qi*dxi=r --> dxi = (r-q*x)/qi
            double qi = constr.getCoefficients().getEntry(dofNum);
            if(Math.abs(qi)>1e-8){//this constr actually involves this DOF
                double curConstrVal = constr.getCoefficients().dotProduct(vec);//q*x
                double xi = vec.getEntry(dofNum);
                double dx = (constr.getValue()-curConstrVal)/qi;//how much we must move x to get to the boundary

                if( (constr.getRelationship()==Relationship.GEQ) == (qi>0) )//constraint is xi >= cur_xi + dx
                    minVal = Math.max(minVal, xi+dx);
                else//xi <= cur_xi+dx
                    maxVal = Math.min(maxVal, xi+dx);

                //check feasibility
                if( (curConstrVal>=constr.getValue()) != (constr.getRelationship()==Relationship.GEQ) ){
                    if(Math.abs(dx)>0.001){//this should be true within numerical error...
                        ObjectIO.writeObject(polytope, "POLYTOPE.dat");
                        System.out.println("VEC: "+vec);
                        System.out.println("DOF NUM: "+dofNum);
                        System.out.println("CONSTR: "+constr);
                        throw new RuntimeException("ERROR: infeasible pt");
                    }
                }
                if(constr.getRelationship()==Relationship.EQ)
                    throw new RuntimeException("ERROR: Equality constraints not supported here");
            }
        }
        
        if(Double.isInfinite(minVal)||Double.isInfinite(maxVal)){//no constraints on this DOF...must not be active in this voxel
            if(Double.isFinite(minVal)||Double.isFinite(maxVal)){
                ObjectIO.writeObject(polytope, "POLYTOPE.dat");
                System.out.println("VEC: "+vec);
                System.out.println("DOF NUM: "+dofNum);
                throw new RuntimeException("ERROR: One-sided bounds on DOF");
            }
            maxVal = minVal = 0;
        }
        
        return new double[] {minVal,maxVal};
    }
     
    
    private ArrayList<int[]> drawVoxSamples(FeatureSet featSet, PruningMatrix pruneMat, PolytopeMatrix plugMat){
        //this assumes all the tuples with terms in featSet are feasible!  
        //For example the constructor of SparseLinearEnergy from an EPIC matrix ensures this
        VoxDrawTupleExpander te = new VoxDrawTupleExpander(pruneMat, plugMat, featSet);
        TESampleSet tss = new TESampleSet(te);
        ArrayList<RCTuple> tupleList = featSet.unprunedFeatTuples();
        
        for(int t=0; t<tupleList.size(); t++){
            te.getTuples().add(tupleList.get(t));
            tss.addTuple(t);
            tss.updateSamples(t);
        }
        
        if(tss.getNumSamples() < 2*tupleList.size()){
             te.numSampPerFeach *= 2*te.getTuples().size()/tss.getNumSamples()+1;//+1 to round up
            
            for(int t=0; t<tupleList.size(); t++)//get the additional samples we need
                tss.updateSamples(t);
        }
        
        return tss.getSamples();
    }
    
    public BasicConfSampleSet(ArrayList<ConfSample> samples){
        this.samples = samples;//sometimes shallow copy is good.  Leave out energies for now
    }

    @Override
    public void checkError(EnergyModel f) {
        double ssd = 0;
        for(int index=0; index<samples.size(); index++){
            double dev = f.evalEnergy(samples.get(index)) - energies.get(index);
            ssd += dev*dev;
        }
        System.out.println("RMS error: "+Math.sqrt(ssd/samples.size()));
    }

    @Override
    public EnergiedConfSampleSet integrateEnergies(IntegrableDOF dof, EnergyModel E) {
        BasicConfSampleSet ans = new BasicConfSampleSet(samples);
        ans.energies = new ArrayList<>();
        for(ConfSample samp : samples){
            ans.energies.add(dof.boltzmannIntegratedEnergy(E, samp));
        }
        return ans;
    }

    @Override
    public EnergiedConfSampleSet evalEnergies(EnergyModel E) {
        BasicConfSampleSet ans = new BasicConfSampleSet(samples);
        ans.energies = new ArrayList<>();
        for(ConfSample samp : samples){
            ans.energies.add(E.evalEnergy(samp));
        }
        return ans;
    }
    
    
    @Override
    public void checkPLUGViolations(PolytopeMatrix plugMat){
        int numViolations = 0;
        for(ConfSample samp : samples){
            ArrayList<DegreeOfFreedom> allDOFs = plugMat.cSpace.confDOFs;
            ArrayList<LinearConstraint> tope = plugMat.getFullStericPolytope(samp.vox, allDOFs);
            double[] pt = new double[allDOFs.size()];
            for(int d=0; d<allDOFs.size(); d++){
                String dofName = allDOFs.get(d).getName();
                pt[d] = samp.contDOFVals.containsKey(dofName) ? samp.contDOFVals.get(dofName) : 0;
            }
            if(!LPChecks.isPointInPolytope(tope, pt))
                numViolations++;
        }
        System.out.println("BasicConfSampleSet has "+numViolations+" PLUG violations among "+samples.size()+" samples");
    }
    
    @Override
    public int getNumSamples() {
        return samples.size();
    }
    
    
}
