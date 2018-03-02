/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ibis;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import java.util.ArrayList;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.ConjugateGradient;
import org.apache.commons.math3.linear.RealLinearOperator;
import org.apache.commons.math3.linear.RealVector;

/**
 *
 * Fit a SparseLinearEnergy by iteratively fitting manageable numbers of parameters
 * 
 * @author mhall44
 */
public class TupIterFitter implements SLEFitter {

    static double trainingResidImprovementThresh = 0.01;//DEBUG!!! 
    
    
    FeatureSet featSet;
    ConfSampleSet sampleSet;
    double curResid;
    
    
    public TupIterFitter(FeatureSet featSet, ConfSampleSet trainingSamples) {
        this.featSet = featSet;
        sampleSet = trainingSamples;
        if(!featSet.freezeSAPE)
            throw new RuntimeException("ERROR: Frozen SAPE expected");
    }
    
    @Override
    public double[] fitCoeffs(double[] initCoeffs) {
        System.out.println("STARTING ITERFITTER.  Training RMS errors will be shown.");
        long startTime = System.currentTimeMillis();
        double[] coeffs = initCoeffs;
        curResid = Double.POSITIVE_INFINITY;
                
        while(needToKeepFitting(coeffs)){
            if(ContIntegrableDOF.usePLUGIntegral){
                refitParams(new DiscreteParamsFittingPartition(featSet,sampleSet), coeffs);
                for(RCTuple tup : featSet.unprunedFeatTuples()){
                    if(featSet.featMatrix.getTupleValue(tup) instanceof DenseFeatureSet.EPIC)
                        refitParams(new RCTupleFittingPartition(featSet,sampleSet,tup), coeffs);
                }
            }
            else {
                refitParamsModLSQ(null, coeffs);
                for(RCTuple tup : featSet.unprunedFeatTuples()){
                    if(featSet.featMatrix.getTupleValue(tup) instanceof DenseFeatureSet.EPIC)
                        refitParamsModLSQ(tup, coeffs);
                }
            }
        }
        
        System.out.println("FITTING DONE.  Time (ms): "+(System.currentTimeMillis()-startTime));
        return coeffs;
    }
    
    boolean needToKeepFitting(double[] coeffs){
        double newResid = sampleSet.checkError(new SparseLinearEnergy(featSet, coeffs));
        boolean keepFitting = (newResid < curResid - trainingResidImprovementThresh);
        curResid = newResid;
        return keepFitting;
    }
    
    
    private void refitParamsModLSQ(RCTuple tup, double[] coeffs){
        int oldNumSamples = -1;//DEBUG!!  iterated til number of samples used stabilizes, for convenience
        
        while(true) {
            //we'll use non-PLUG samples for this,
            //but will exclude them if they are correctly placed above threshold
            ConfSampleSet iterSamps = iterSampleSet(tup,coeffs);
            int sampleStabilizationThresh = Math.max(iterSamps.getNumSamples()/1000, 1);//DEBUG!!
            
            FittingPartition fp;
            if(tup==null)//fit discrete params
                fp = new DiscreteParamsFittingPartition(featSet,iterSamps);
            else
                fp = new RCTupleFittingPartition(featSet,sampleSet,tup);
            refitParams(fp,coeffs);
            
            if(Math.abs(iterSamps.numSamples-oldNumSamples)<sampleStabilizationThresh)
                break;
            oldNumSamples = iterSamps.numSamples;
        }
    }
    
    ConfSampleSet iterSampleSet(RCTuple fittingTup, double[] coeffs){
        //figure out which samples to use for an iteration of modified least squares
        //just exclude those unplugged samples whose energy is correctly above threshold
        //threshold is minimum PLUG energy value + MODLSQ_THRESH_BUFFER
        ConfSampleSet ans = new ConfSampleSet();
        SparseLinearEnergy f = new SparseLinearEnergy(featSet, coeffs);//current model of the energy
        
        for(ConfSampleSet.SampleBatch oldBatch : sampleSet.batches){
            
            if(fittingTup!=null)//skip useless samples for this tuple
                if(!RCTupleFittingPartition.matches(oldBatch.vox, fittingTup))
                    continue;
            
            ConfSampleSet.SampleBatch newBatch = ans.new SampleBatch(oldBatch.vox);

            newBatch.contDOFValsArr.addAll(oldBatch.contDOFValsArr);
            newBatch.energies.addAll(oldBatch.energies);
            ans.numSamples += oldBatch.energies.size();

            MoleculeModifierAndScorer mms = f.makeMMS(new ConfSample(oldBatch.vox,null));

            double minPluggedE = Double.POSITIVE_INFINITY;
            for(double[] cdv : oldBatch.contDOFValsArr){
                minPluggedE = Math.min(minPluggedE, f.evalEnergy(new ConfSample(oldBatch.vox,cdv), mms));
            }

            double thresh = minPluggedE + ConfSampleSet.MODLSQ_THRESH_BUFFER;

            for(int index=0; index<oldBatch.contDOFValsArrUnplugged.size(); index++){
                double fval = f.evalEnergy(new ConfSample(oldBatch.vox,oldBatch.contDOFValsArrUnplugged.get(index)), mms);
                if(fval<thresh || oldBatch.energiesUnplugged.get(index)<thresh){
                    //use this unplugged samples
                    newBatch.contDOFValsArr.add(oldBatch.contDOFValsArrUnplugged.get(index));
                    newBatch.energies.add(oldBatch.energiesUnplugged.get(index));
                    ans.numSamples++;
                }
            }

            ans.batches.add(newBatch);
        }
        
        return ans;
    }
    
    private void refitParams(FittingPartition fp, double[] coeffs){
        //before this, prepare the partition, including everything needed to apply the oeprator
        
        RealLinearOperator AtA = fp.makeNormalOperator();
        RealVector Atb = fp.makeNormalRHS(coeffs);
        ArrayList<Integer> refitParamIndices = fp.refitParamIndices;
        
        double initCoeffs[] = new double[fp.numParams];
        for(int p=0; p<fp.numParams; p++)
            initCoeffs[p] = coeffs[refitParamIndices.get(p)];
        
        double refitParamVals[] = conjugateGradientSoln(AtA, Atb, initCoeffs);
        
        for(int q=0; q<refitParamVals.length; q++){//DEBUG!!!
            if(Math.abs(refitParamVals[q])>10000){
                int fff = 456;
            }
        }
               
        
        for(int p=0; p<fp.numParams; p++)
            coeffs[refitParamIndices.get(p)] = refitParamVals[p];
    }
 
    
    
    
    
    static double[] conjugateGradientSoln(RealLinearOperator M, RealVector v, double[] initCoeffs){
        ConjugateGradient cg = new ConjugateGradient(100000,1e-8,false);//max_iter; delta; whether to check pos def
        //delta is target ratio of residual norm to true vals norm
        //long startTime = System.currentTimeMillis();
        RealVector ans = cg.solve(M, v, new ArrayRealVector(initCoeffs));
        //System.out.println( "Conjugate gradient fitting time (ms): " + (System.currentTimeMillis()-startTime) );
        return ans.toArray();    
    }
}
