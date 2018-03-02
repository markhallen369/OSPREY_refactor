/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.aibis;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.epic.SeriesFitter;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.RealLinearOperator;

/**
 *
 * Refits the parameters for a particular RCTuple
 * 
 * @author mhall44
 */
public class RCTupleFittingPartition extends FittingPartition {
    
    RCTuple fittingTup;
    DenseFeatureSet.EPIC efs;
    
    ArrayList<ConfSampleSet.SampleBatch> samples;
    
    static double regularizationLambda = 1e-4;
    
    public RCTupleFittingPartition(FeatureSet featSet, ConfSampleSet sampleSet, RCTuple fittingTup){
        this.featSet = featSet;
        this.fittingTup = fittingTup;
        
        DenseFeatureSet fs = featSet.featMatrix.getTupleValue(fittingTup);
        if( ! (fs instanceof DenseFeatureSet.EPIC))
            throw new RuntimeException("ERROR: RCTupleFittingPartition currently just for EPIC features");
        efs = (DenseFeatureSet.EPIC)fs;
        numParams = efs.getNumFeatures();
        
        refitParamIndices = new ArrayList<>();
        int offset = featSet.featOffsets.getTupleValue(fittingTup);
        for(int p=0; p<numParams; p++)
            refitParamIndices.add(offset+p);
        
        
        //DEBUG!!!  may not need so many samples...  (though if sample set isn't excessive probably not a huge issue on average)
        samples = new ArrayList<>();
        numSamples = 0;
        for(ConfSampleSet.SampleBatch batch : sampleSet.batches){
            //check if batch.vox matches fittingTup
            if(matches(batch.vox, fittingTup)){
                samples.add(batch);
                numSamples += batch.contDOFValsArr.size();
            }
        }
        
        prepareOperators();
    }
    
    
    private static boolean matches(RCTuple a, RCTuple b){
        //is b a part of a?
        for(int index=0; index<b.size(); index++){
            if(a.RCAtPos(b.pos.get(index)) != b.RCs.get(index))
                return false;
        }
        return true;
    }
    
    
    private void prepareOperators() {
        Array2DRowRealMatrix directOp = new Array2DRowRealMatrix(numSamples, numParams);
        
        int numDOFs = efs.dofIndices.length;
        
        DoubleMatrix1D coeffVals = DoubleFactory1D.dense.make(numParams);
        DoubleMatrix1D dofVals = DoubleFactory1D.dense.make(numDOFs);
        
        int s = 0;
        for(ConfSampleSet.SampleBatch batch : samples){
            for(double[] contDOFValsArr : batch.contDOFValsArr){
                for(int d=0; d<numDOFs; d++)
                    dofVals.set(d, contDOFValsArr[efs.dofIndices[d]]);
                dofVals.assign(efs.center, Functions.minus);
                SeriesFitter.calcSampParamCoeffs(coeffVals, dofVals, numDOFs, true, efs.order, efs.order, null);
                for(int p=0; p<numParams; p++)
                    directOp.setEntry(s, p, coeffVals.get(p));
                s++;
            }
        }
        
        directOperator = directOp;
        transposeOperator = (RealLinearOperator) directOp.transpose();
    }
    
    
    @Override
    double[] computeTargetVals(double coeffs[]) {
        double[] targetVals = new double[numSamples];
        
        double nonFitCoeffs[] = coeffs.clone();
        
        for(int c : refitParamIndices){
            nonFitCoeffs[c] = 0;
        }
        
        int s=0;
        for(ConfSampleSet.SampleBatch batch : samples){
            MoleculeModifierAndScorer mms = featSet.makeMMS(new ConfSample(batch.vox,null));
            for(int index=0; index<batch.contDOFValsArr.size(); index++){
                double otherTermsE = featSet.evalEnergy(new ConfSample(batch.vox,batch.contDOFValsArr.get(index)), nonFitCoeffs, mms);
                targetVals[s] = batch.energies.get(index) - otherTermsE;
                s++;
            }
        }
        
        return targetVals;
    }
    
    
    
    @Override
    RealLinearOperator makeNormalOperator(){
        //Some coefficients are growing very big.  They seem to mostly cancel each other
        //but this may be contributing to error eventually
        //and definitely causes infinite energies once exp overflows.  
        //Fighting this with a little regularization
        Array2DRowRealMatrix AtA = ((Array2DRowRealMatrix)transposeOperator).multiply((Array2DRowRealMatrix)directOperator);
        for(int p=0; p<numParams; p++)
            AtA.addToEntry(p, p, regularizationLambda);//parameter sizes penalized much less than energy dev of course
        return AtA;
    }
    
}
