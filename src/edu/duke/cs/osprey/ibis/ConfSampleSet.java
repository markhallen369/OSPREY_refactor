/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ibis;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.plug.LPChecks;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tupexp.BasicPruningTupleExpander;
import edu.duke.cs.osprey.tupexp.TESampleSet;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import static java.util.stream.Collectors.toList;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.Relationship;

/**
 *
 * This sample set wants to avoid spending a lot of time initializing samples
 * So it draws only as many voxels as needed to get the discrete parameters
 * and draw a bunch of continuous voxels for each
 * 
 * 
 * @author mhall44
 */
public class ConfSampleSet {
    
    static double MODLSQ_THRESH_BUFFER = 10;//DEBUG!!!
    
    protected class SampleBatch {
        RCTuple vox;
        //ArrayList<HashMap<String,Double> > contDOFVals = new ArrayList<>();
        ArrayList<Double> energies = new ArrayList<>();
        ArrayList<double[]> contDOFValsArr = new ArrayList<>();//DEBUG!!
        
        
        //DEBUG!!!  For full interval integration we will need sample outside the PLUG zone
        ArrayList<Double> energiesUnplugged = new ArrayList<>();
        ArrayList<double[]> contDOFValsArrUnplugged = new ArrayList<>();
        
        protected SampleBatch(RCTuple vox){this.vox = vox;}
    }
    
    ArrayList<SampleBatch> batches;
    int numSamples;
    
    public ConfSampleSet(){
        //starting empty
        batches = new ArrayList<>();
        numSamples = 0;
    }
    
    public ConfSampleSet(FeatureSet featSet, PolytopeMatrix plugMtx, PruningMatrix pruneMat){
        //The main issue in sampling is since we have a lot of terms with small support,
        //make sure that each RC tuple in featSet is represented by enough samples
        //and also that the total number of samples is > number of params
        //DEBUG!!! might be convenient to test if this really matters much
        //Anyway constructor will just draw samples, energies can be handled separately

        List<RCTuple> voxSamples = drawVoxSamples(featSet,pruneMat,plugMtx).stream().map(q->new RCTuple(q)).collect(toList());
        int[] batchSizes = chooseBatchSizes(voxSamples, featSet);
        
        batches = new ArrayList<>();
        numSamples = 0;
        ArrayList<DegreeOfFreedom> DOFs = plugMtx.cSpace.confDOFs;
        for(int voxSampNum=0; voxSampNum<voxSamples.size(); voxSampNum++){
            RCTuple tup = voxSamples.get(voxSampNum);
            ArrayList<LinearConstraint> tope = plugMtx.getFullStericPolytope(tup,DOFs);//DEBUG!! again, uses funny def'n of steric
            SampleBatch batch = new SampleBatch(tup);
            for(int sampNum=0; sampNum<batchSizes[voxSampNum]; sampNum++){
                ArrayList<Double> dofVals = sampleFeasPt(tope);
                if(dofVals.isEmpty())//no constraints or DOFs...just add all 0s
                    dofVals.addAll(Collections.nCopies(DOFs.size(),0.));

                /*HashMap<String,Double> contDOFVals = new HashMap<>();
                for(int dofNum=0; dofNum<DOFs.size(); dofNum++)
                    contDOFVals.put(DOFs.get(dofNum).getName(), dofVals.get(dofNum));
                //batch.contDOFVals.add(contDOFVals);
                batch.contDOFValsArr.add(makeContDOFValsArr(plugMtx.cSpace, contDOFVals));
                */
                double contDOFValsArr[] = new double[DOFs.size()];
                for(int d=0; d<DOFs.size(); d++)
                    contDOFValsArr[d] = dofVals.get(d);
                batch.contDOFValsArr.add(contDOFValsArr);
                
                numSamples++;
            }
            
            if(!ContIntegrableDOF.usePLUGIntegral){
                //we want a somewhat sterically feasible sample,
                //but allow some sampling everywhere
                //so create a truly random point and randomly sample along the line
                //between that point and our PLUG-constrained sample
                
                for(double cdv[] : batch.contDOFValsArr){
                    double alpha = Math.random();
                    
                    double cdvUnplugged[] = new double[cdv.length];
                    for(int dof=0; dof<cdv.length; dof++){
                        double[] boxConstr = ContIntegrableDOF.calcFullInterval(batch.vox, plugMtx, 
                                plugMtx.cSpace.confDOFs.get(dof).getName(), false);
                        if(cdv[dof]<boxConstr[0]-1e-6||cdv[dof]>boxConstr[1]+1e-6)
                            throw new RuntimeException("ERROR: PLUG value inconsistent with box constr");
                        double randPt = boxConstr[0] + Math.random()*(boxConstr[1]-boxConstr[0]);
                        cdvUnplugged[dof] = randPt*alpha + cdv[dof]*(1-alpha);
                    }
                    
                    batch.contDOFValsArrUnplugged.add(cdvUnplugged);
                }
            }
            
            
            
            batches.add(batch);
        }
    }
    
    
    double[] makeContDOFValsArr(ConfSpace cSpace, HashMap<String,Double> contDOFVals){
        List<String> allDOFNames = cSpace.confDOFs.stream().map(s->s.getName()).collect(toList());
        double[] contDOFValsArr = new double[allDOFNames.size()];
        for(String name : contDOFVals.keySet())
            contDOFValsArr[allDOFNames.indexOf(name)] = contDOFVals.get(name);
        return contDOFValsArr;
    }
    
    
    private int[] chooseBatchSizes(List<RCTuple> voxSamples, FeatureSet featSet){
        //given the voxel samples defining each batch,
        //decide how many samples to put in each batch
        int[] ans = new int[voxSamples.size()];
        int numSampPerFeach = 10;//DEBUG!!
        int totTupSampToGo = featSet.numFeatures*numSampPerFeach;//how many (tuple,sample) pairs still need to be drawn
        TupleMatrixGeneric<Integer> numSampForTups = new TupleMatrixGeneric<Integer>(
                featSet.featMatrix.getNumPos(), 
                featSet.featMatrix.getNumConfAtPos(), Double.POSITIVE_INFINITY, 0);
        
        while(totTupSampToGo>0){
            //loop through the batches, adding to them if it helps
            //adding one at a time in order to distribute more evenly across batches
            boolean gotSomething = false;
            for(int b=0; b<voxSamples.size(); b++){
                boolean sampleNeeded = false;
                for( RCTuple tup : featSet.tuplesForSamp(new ConfSample(voxSamples.get(b),null)) ){
                    int tupNumFeat = featSet.featMatrix.getTupleValue(tup).getNumFeatures();
                    if(tupNumFeat*numSampPerFeach > getNumSampForTup(numSampForTups,tup)){
                        sampleNeeded = true;
                        break;
                    }
                }
                if(sampleNeeded){//add the sample
                    ans[b]++;
                    gotSomething = true;
                    for( RCTuple tup : featSet.tuplesForSamp(new ConfSample(voxSamples.get(b),null)) ){
                        DenseFeatureSet fs = featSet.featMatrix.getTupleValue(tup);
                        if(fs==null)
                            throw new RuntimeException("ERROR: Drew pruned voxel");
                        int tupNumFeat = fs.getNumFeatures();
                        if(tupNumFeat*numSampPerFeach > getNumSampForTup(numSampForTups,tup)){
                            totTupSampToGo--;
                        }
                        numSampForTups.setTupleValue(tup, getNumSampForTup(numSampForTups,tup)+1);
                    }
                }
                if(totTupSampToGo<=0)
                    break;
            }
            
            if(!gotSomething)//would lead to infinite loop...complain
                throw new RuntimeException("ERROR: Still need more samples for some tuples but couldn't get them for any batch!");
        }
        
        return ans;
    }
    
    
    private static int getNumSampForTup(TupleMatrixGeneric<Integer> numSampForTups, RCTuple tup){
        //look up with 0 by default
        Integer val = numSampForTups.getTupleValue(tup);
        if(val==null)
            return 0;
        else
            return val;
    }
    

     
    
    private ArrayList<int[]> drawVoxSamples(FeatureSet featSet, PruningMatrix pruneMat, PolytopeMatrix plugMat){
        //this assumes all the tuples with terms in featSet are feasible!  
        //For example the constructor of SparseLinearEnergy from an EPIC matrix ensures this
        BasicPruningTupleExpander te = new BasicPruningTupleExpander(pruneMat, plugMat);
        TESampleSet tss = new TESampleSet(te);
        ArrayList<RCTuple> tupleList = featSet.unprunedFeatTuples();
        
        for(int t=0; t<tupleList.size(); t++){
            te.getTuples().add(tupleList.get(t));
            tss.addTuple(t);
        }
        
        if(tss.getNumSamples() < 2*tupleList.size()){//overall number of voxels must be big enough to not give an underdetermined problem
             te.scaleUpNumSampsPerTuple( 2*te.getTuples().size()/tss.getNumSamples()+1 );//+1 to round up
            
            for(int t=0; t<tupleList.size(); t++)//get the additional samples we need
                tss.updateSamples(t);
        }
        
        return tss.getSamples();
    }
    
    private ConfSampleSet(ArrayList<SampleBatch> batches, int numSamples){
        this.batches = batches;//sometimes shallow copy is good. 
        this.numSamples = numSamples;
    }

    public double checkError(SparseLinearEnergy f) {
        double ssd = 0;
        int numUnpluggedSamples = 0;//how many of these did we use (only use sub-threshold ones)
        for(SampleBatch batch : batches){
            MoleculeModifierAndScorer mms = f.makeMMS(new ConfSample(batch.vox,null));
            double thresh = Double.POSITIVE_INFINITY;//DEBUG!!!  thresh based on plug samples
            for(int index=0; index<batch.contDOFValsArr.size(); index++){
                double dev = f.evalEnergy(new ConfSample(batch.vox,batch.contDOFValsArr.get(index)), mms)
                        - batch.energies.get(index);
                ssd += dev*dev;
                
                if(Double.isInfinite(dev)){
                    int aaa = 0;
                }
                
                thresh = Math.min(thresh,batch.energies.get(index)+MODLSQ_THRESH_BUFFER);
            }
            
            for(int index=0; index<batch.contDOFValsArrUnplugged.size(); index++){
                double fval = f.evalEnergy(new ConfSample(batch.vox,batch.contDOFValsArrUnplugged.get(index)), mms);
                double dev = fval - batch.energiesUnplugged.get(index);
                if(fval<thresh || batch.energiesUnplugged.get(index)<thresh){
                    ssd += dev*dev;
                    numUnpluggedSamples++;
                }
            }
        }
        double rmsError = Math.sqrt(ssd/(numSamples+numUnpluggedSamples));
        System.out.println("RMS error: "+rmsError);
        
        if(Double.isInfinite(rmsError)){
            int bbb = 4;
        }
        
        return rmsError;
    }

    public ConfSampleSet integrateEnergies(IntegrableDOF dof, SparseLinearEnergy E) {
        ConfSampleSet ans = new ConfSampleSet(batches, numSamples);
        for(SampleBatch batch : batches){
            batch.energies.clear();
            for(int index=0; index<batch.contDOFValsArr.size(); index++){
                ConfSample confSamp = new ConfSample(batch.vox,batch.contDOFValsArr.get(index));
                
                
                /*if(dof.description().contains("DIH41.1")){//DEBUG!!!
                    if(index==0 && batch==batches.get(0)){//This one is interesting
                        RCTuple interestingVox = new RCTuple(new int[] {1,17,3,8,2,4,4});
                        double[] interestingDOFVals = new double[] {177.53009750733216,0,0,0,0,
                            -64.31559509756225,90.5543885241102,5.652702575046551,
                            -60.47593826638652,0,0};
                        ConfSample interestingSamp = new ConfSample(interestingVox, interestingDOFVals);
                        ((ContIntegrableDOF)dof).scanIntegratedEnergies(E, interestingSamp, 5);
                    }
                    ((ContIntegrableDOF)dof).scanIntegratedEnergies(E, confSamp, 5);
                }*/
                
                
                double sampleE = dof.boltzmannIntegratedEnergy(E, confSamp);
                batch.energies.add(sampleE);
            }
            
            batch.energiesUnplugged.clear();
            for(int index=0; index<batch.contDOFValsArrUnplugged.size(); index++){
                ConfSample confSamp = new ConfSample(batch.vox,batch.contDOFValsArrUnplugged.get(index));
                double sampleE = dof.boltzmannIntegratedEnergy(E, confSamp);
                batch.energiesUnplugged.add(sampleE);
            }
        }
        return ans;
    }

    public ConfSampleSet evalEnergies(SparseLinearEnergy E) {
        ConfSampleSet ans = new ConfSampleSet(batches, numSamples);
        for(SampleBatch batch : batches){
            batch.energies.clear();
            MoleculeModifierAndScorer mms = E.makeMMS(new ConfSample(batch.vox,null));
            for(int index=0; index<batch.contDOFValsArr.size(); index++){
                double sampleE = E.evalEnergy(new ConfSample(batch.vox,batch.contDOFValsArr.get(index)),mms);
                batch.energies.add(sampleE);
            }
            
            batch.energiesUnplugged.clear();
            for(int index=0; index<batch.contDOFValsArrUnplugged.size(); index++){
                double sampleE = E.evalEnergy(new ConfSample(batch.vox,batch.contDOFValsArrUnplugged.get(index)),mms);
                batch.energiesUnplugged.add(sampleE);
            }
        }
        return ans;
    }
    
    
    public void checkPLUGViolations(PolytopeMatrix plugMat){
        int numViolations = 0;
        for(SampleBatch batch : batches){
            for(int cdfIndex=0; cdfIndex<batch.contDOFValsArr.size(); cdfIndex++){
                //for(HashMap<String,Double> cdf : batch.contDOFVals){
                ConfSample samp = new ConfSample(batch.vox, batch.contDOFValsArr.get(cdfIndex));
                ArrayList<DegreeOfFreedom> allDOFs = plugMat.cSpace.confDOFs;
                ArrayList<LinearConstraint> tope = plugMat.getFullStericPolytope(samp.vox, allDOFs);
                /*double[] pt = new double[allDOFs.size()];
                for(int d=0; d<allDOFs.size(); d++){
                    String dofName = allDOFs.get(d).getName();
                    pt[d] = samp.contDOFVals.containsKey(dofName) ? samp.contDOFVals.get(dofName) : 0;
                }*/
                if(!LPChecks.isPointInPolytope(tope, samp.contDOFValsArr))//pt))
                    numViolations++;
            }
        }
        System.out.println("ContBatchSampleSet has "+numViolations+" PLUG violations among "+numSamples+" samples");
    }
    
    
    
    
    public int getNumSamples() {
        return numSamples;
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
    
}
