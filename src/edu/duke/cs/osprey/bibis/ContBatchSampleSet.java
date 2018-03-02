/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.bibis;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.plug.LPChecks;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
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
public class ContBatchSampleSet implements EnergiedConfSampleSet {
    
    protected class SampleBatch {
        RCTuple vox;
        ArrayList<HashMap<String,Double> > contDOFVals = new ArrayList<>();
        ArrayList<Double> energies = new ArrayList<>();
        ArrayList<double[]> contDOFValsArr = new ArrayList<>();//DEBUG!!
        
        protected SampleBatch(RCTuple vox){this.vox = vox;}
    }
    
    ArrayList<SampleBatch> batches;
    int numSamples;
    
    public ContBatchSampleSet(FeatureSet featSet, PolytopeMatrix plugMtx, PruningMatrix pruneMat){
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
                ArrayList<Double> dofVals = BasicConfSampleSet.sampleFeasPt(tope);
                if(dofVals.isEmpty())//no constraints or DOFs...just add all 0s
                    dofVals.addAll(Collections.nCopies(DOFs.size(),0.));

                HashMap<String,Double> contDOFVals = new HashMap<>();
                for(int dofNum=0; dofNum<DOFs.size(); dofNum++)
                    contDOFVals.put(DOFs.get(dofNum).getName(), dofVals.get(dofNum));
                batch.contDOFVals.add(contDOFVals);
                batch.contDOFValsArr.add(makeContDOFValsArr(plugMtx.cSpace, contDOFVals));
                numSamples++;
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
    
    private ContBatchSampleSet(ArrayList<SampleBatch> batches, int numSamples){
        this.batches = batches;//sometimes shallow copy is good. 
        this.numSamples = numSamples;
    }

    @Override
    public void checkError(EnergyModel f) {
        double ssd = 0;
        for(SampleBatch batch : batches){
            MoleculeModifierAndScorer mms = f.makeMMS(new ConfSample(batch.vox,null));
            for(int index=0; index<batch.contDOFVals.size(); index++){
                double dev = f.evalEnergy(new ConfSample(batch.vox,batch.contDOFVals.get(index),batch.contDOFValsArr.get(index)), mms)
                        - batch.energies.get(index);
                ssd += dev*dev;
            }
        }
        System.out.println("RMS error: "+Math.sqrt(ssd/numSamples));
    }

    @Override
    public EnergiedConfSampleSet integrateEnergies(IntegrableDOF dof, EnergyModel E) {
        ContBatchSampleSet ans = new ContBatchSampleSet(batches, numSamples);
        for(SampleBatch batch : batches){
            batch.energies.clear();
            for(int index=0; index<batch.contDOFVals.size(); index++){
                double sampleE = dof.boltzmannIntegratedEnergy(E, new ConfSample(batch.vox,batch.contDOFVals.get(index),batch.contDOFValsArr.get(index)));
                batch.energies.add(sampleE);
            }
        }
        return ans;
    }

    @Override
    public EnergiedConfSampleSet evalEnergies(EnergyModel E) {
        ContBatchSampleSet ans = new ContBatchSampleSet(batches, numSamples);
        for(SampleBatch batch : batches){
            batch.energies.clear();
            MoleculeModifierAndScorer mms = E.makeMMS(new ConfSample(batch.vox,null));
            for(int index=0; index<batch.contDOFVals.size(); index++){
                double sampleE = E.evalEnergy(new ConfSample(batch.vox,batch.contDOFVals.get(index),batch.contDOFValsArr.get(index)),mms);
                batch.energies.add(sampleE);
            }
        }
        return ans;
    }
    
    
    @Override
    public void checkPLUGViolations(PolytopeMatrix plugMat){
        int numViolations = 0;
        for(SampleBatch batch : batches){
            for(int cdfIndex=0; cdfIndex<batch.contDOFVals.size(); cdfIndex++){
                //for(HashMap<String,Double> cdf : batch.contDOFVals){
                ConfSample samp = new ConfSample(batch.vox, batch.contDOFVals.get(cdfIndex), batch.contDOFValsArr.get(cdfIndex));
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
        }
        System.out.println("ContBatchSampleSet has "+numViolations+" PLUG violations among "+numSamples+" samples");
    }
    
    
    @Override
    public int getNumSamples() {
        return numSamples;
    }
    
}
