/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tupexp;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * Specifies the objective function used for fitting: basic least squares
 * 
 * 
 * @author mhall44
 */
class FittingObjFcn implements Serializable {
    
    //OK let's do basic least squares for true values up to thresh
    double thresh = Double.POSITIVE_INFINITY;
    double k;//and for the portion of true value above thresh, we allow relative error k
    
    
    enum WeightingMode {
        REGULAR,
        ONEWT,
        THRESHWT
    }
    
    WeightingMode wtMode = WeightingMode.REGULAR;
    
    FittingObjFcn(){}//basic least squares, regular weighting
    
    FittingObjFcn(double thresh, double k, boolean useRelWt, boolean useThreshWt){//modified least squares
        this.thresh = thresh;
        this.k = k;
        
        if(useThreshWt)
            wtMode = WeightingMode.THRESHWT;
        else if(useRelWt)
            wtMode = WeightingMode.ONEWT;//this seems to be the more effective weighting
        else
            wtMode = WeightingMode.REGULAR;
        
        if(k<0 || k>1)
            throw new RuntimeException("ERROR: Invalid modified least squares k value");
    }
      
    
    double computeResid(double fitVal, double trueVal){
        //compute the residual given the true value of a sample and the value of 
        //the fit at that sample
        double penalizedError;
        //if(trueVal<=thresh)//DEBUG!!!!  TURNING OFF MOD LSQ
            penalizedError = fitVal-trueVal;
        /*else {//DEBUG!!!!!
            double excess = trueVal - thresh;
            double fitExcess  = fitVal-thresh;
            if(fitExcess < (1-k)*excess)
                penalizedError = fitExcess-(1-k)*excess;
            else if(fitExcess > (1+k)*excess)
                penalizedError = fitExcess-(1+k)*excess;
            else
                penalizedError = 0;
        }*/
        
        return penalizedError*penalizedError;
    }
    
    /*double[] goodRegionBoundsForSample(double trueVal){
        //given the true value for a sample, return the bounds on the "good region"
        //for the fit value for that sample
        double bounds[] = new double[2];
        if(trueVal<=thresh){//only trueVal is good
            bounds[0] = trueVal;
            bounds[1] = trueVal;
        }
        else {//penalize only larger deviations
            bounds[0] = thresh + (1-k)*(trueVal-thresh);
            bounds[1] = thresh + (1+k)*(trueVal-thresh);
        }
        
        return bounds;
    }*/
    
    
    CGTupleFitter makeTupleFitter(TupleIndexMatrix tim, ArrayList<int[]> samp, int numTuples, double[] trueVals){
       
        ArrayList<Double> tv = new ArrayList<>();//trueVals as ArrayList for use by computeWeights
        for(double tvv : trueVals)
            tv.add(tvv);
        ArrayList<Double> weights = computeWeights(tv, 0);
        
        
        //if(thresh==Double.POSITIVE_INFINITY)//basic least squares//DEBUG!!!
            return new CGTupleFitter(tim, samp, numTuples, trueVals, weights);
        
            //DEBUG!!!!
        /*ArrayList<double[]> goodRegionBounds = new ArrayList<>();
        int numSamp = samp.size();
        for(int s=0; s<numSamp; s++)
            goodRegionBounds.add(goodRegionBoundsForSample(trueVals[s]));
        
        return new IterativeCGTupleFitter(tim, samp, numTuples, goodRegionBounds);*/
    }
    
    
    
    
    //weighted processing of sample set
    double unnormalizedWeight(double trueVal){//trueVal has 0 as constTerm
        double relCutoff = 1;//cutoff above which to use relative weighting
        switch(wtMode){
            case REGULAR:
                return 1;
            case THRESHWT:
                relCutoff = thresh;
            case ONEWT:
                break;
                //onewt will keep 1
        }
        //we want relative weighting for large trueVal (errors here relatively less important)
        //but not for small (would overweight) or negative (doesn't make sense)
        //this scheme transitions continuously between the two
        if(trueVal<relCutoff)
            return 1;
        else{
            //relative weighting for square residual,
            //scaled to be continuous
            double x = relCutoff/trueVal;
            return x*x;
        }
    }
    
    
    ArrayList<Double> computeWeights(ArrayList<Double> offsetTrueVals, double energyOffset){
        //weights for the whole sample set
        int numSamp = offsetTrueVals.size();
        ArrayList<Double> weights = new ArrayList<>();
        double wtsum = 0;
        
        for(int s=0; s<numSamp; s++){
            double targetVal = offsetTrueVals.get(s) - energyOffset;
            double unnormWt = unnormalizedWeight(targetVal);
            weights.add(unnormWt);
            wtsum += unnormWt;
        }
        
        for(int s=0; s<numSamp; s++)
            weights.set( s, weights.get(s) * (numSamp/wtsum) );
        
        return weights;
    }
    
    ArrayList<Double> computeAllResids(ArrayList<Double> offsetTrueVals, ArrayList<Double> offsetFitVals, 
            double energyOffset){
        
        ArrayList<Double> weights = computeWeights(offsetTrueVals, energyOffset);
        ArrayList<Double> sampleResids = new ArrayList<>();
        int numSamp = weights.size();
        
        for(int s=0; s<numSamp; s++){
            double targetVal = offsetTrueVals.get(s) - energyOffset;
            double sampResid = computeResid(offsetFitVals.get(s)-energyOffset, targetVal);
            
            sampleResids.add(weights.get(s)*sampResid);
        }
        
        return sampleResids;
    }
    
}
