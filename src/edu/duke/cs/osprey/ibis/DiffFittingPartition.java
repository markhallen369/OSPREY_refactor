/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


package edu.duke.cs.osprey.ibis;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import java.util.ArrayList;
import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.RealLinearOperator;

/**
 *
 * Fit differences between samples in which only the movingDOFs move, and only within the voxel
 * We fit only parameters whose features depend on the movingDOFs, 
 * and assuming that the terms depending on both movingDOFs and other stuff have been fit properly using second derivs,
 * we will get valid values for those parameters
 * 
 * 
 * @author mhall44
 */

/*
public class DiffFittingPartition extends FittingPartition {

    int[] movingDOFs;//in each diff, only the specified DOFs move.  
    
    //we'll only do fits where these dofs apply, within the voxel where we are trying to apply them.  
    
    
    
    
    private class SamplePair {
        //A pair of samples.  
        batchThese;
        //just create diffs for all moving dofs in sample set on boolean flags, shouldn't be super expensive
        int[] vox;
        double[] contDOFVals1, contDOFVals2;
        double targetE1, targetE2;
    }
    
    ArrayList<SamplePair> samples;
    
    
    
    
    
    //so what DenseFeatureSets do we care about??
    //When evaluating, we will cancel out any DenseFeatureSet (or SAPE) contribs not depending on the movingDOFs
    //if a term does depend on the movingDOFs then it is expected to have fittable params though
    //(possibly along with non-fittable ones)
    
    //ASSUMPTIONS: cancellation of terms w/o moving dofs; fittable coeffs are <=quadratic
    //what params are we actually fitting:
    //answer: linear and quadratic terms limited to the movingDOFs.  
    //so our params, then, are those applicable terms for each tuple.  
    //hmm we should probably explicitly build a record of them.  
    private class FittablePolynomial {//the portion of the polynomial in a DenseFeatureSet.EPIC
        //that we can actually fit.  All the fittable params are part of one of these
        FittablePolynomial(DenseFeatureSet.EPIC fs, int offset){
            this.offset = offset;
            if(fs.order!=2)
                throw new RuntimeException("ERROR: Expected quadratic EPIC here");//may need more later...
            
            for(int i=0; i<fs.dofIndices.length;  i++){
                if(contains(movingDOFs,fs.dofIndices[i])){
                    //fittable linear monomial
                    linTermDOFs.add(fs.dofIndices[i]);
                    refitParamIndices.add(offset+i+1);
                }
            }
            int coeffCount = fs.dofIndices.length+1;//counting coeffs within EPIC term
            for(int i=0; i<fs.dofIndices.length; i++){
                for(int j=0; j<=i; j++){
                    if(contains(movingDOFs,fs.dofIndices[i]) && contains(movingDOFs,fs.dofIndices[j])){
                        //fittable quadratic monomial
                        quadTermDOFs.add(new int[] {fs.dofIndices[i],fs.dofIndices[j]});
                        refitParamIndices.add(offset+coeffCount);
                    }
                    coeffCount++;
                }
            }
                
            //also update refitParamIndices to account for the new FittablePolynomial
        }
        
        int offset;//offset of this FittablePolynomial's coeffs in the overall fittable coeffs
        ArrayList<Integer> linTermDOFs = new ArrayList<>();
        ArrayList<int[]> quadTermDOFs = new ArrayList<>();
    }
    
    TupleMatrixGeneric<FittablePolynomial> fittablePolynomials;
    
    
    DiffFittingPartition(FeatureSet featSet, ConfSampleSet sampleSet, int movingDOFs[], RCTuple restrTup){
        //for fitting diffs in the moving DOFs, restricted to confs in restrTup
        this.featSet = featSet;
        this.movingDOFs = movingDOFs;
                
        //figure out parameters 
        int curParamIndex = 0;
        fittablePolynomials = new TupleMatrixGeneric<>(featSet.confSpace,Double.POSITIVE_INFINITY,null);
        refitParamIndices = new ArrayList<>();
        samples = new ArrayList<>();
        for(RCTuple tup : featSet.unprunedFeatTuples()){
            DenseFeatureSet tupFS = featSet.featMatrix.getTupleValue(tup);
            if(tupFS instanceof DenseFeatureSet.EPIC){
                DenseFeatureSet.EPIC eTupFS = (DenseFeatureSet.EPIC)tupFS;
                if(hasOverlap(eTupFS.dofIndices, movingDOFs)){//tup has non-cancelling terms
                    FittablePolynomial fpoly = new FittablePolynomial(eTupFS, curParamIndex);
                    fittablePolynomials.setTupleValue(tup, fpoly);
                    int fpolyNumParams = fpoly.linTermDOFs.size() + fpoly.quadTermDOFs.size();
                    curParamIndex += fpolyNumParams;
                    addSamplesForTup(tup, sampleSet, fpolyNumParams);
                }
            }
        }
        
        numParams = curParamIndex;
        numSamples = samples.size();
        
        prepareOperators();
    }
    
    
    private void addSamplesForTup(RCTuple tup, ConfSampleSet sampleSet, int numParamsForTup){
        samples.add(asNeeded);//one option is just to make sure everything with the specified tuple
        //goes into samples...
    }
    
    
    private void prepareOperators() {
        OpenMapRealMatrix directOp = new OpenMapRealMatrix(numSamples, numParams);
        //DEBUG!!  May be faster to implement a RealLinearOperator with ArrayLists of indices, entries for each row?
        //since we just want this sparse matrix for matrix-vec mult
        //same goes for CrossDerivFittingPartition
        
        for(int s=0; s<numSamples; s++){
            for(RCTuple tup : featSet.tuplesForSamp(new ConfSample(samples.get(s).vox,null))){
                //if(hasMovingDOFs(tup)){
                FittablePolynomial fpoly = fittablePolynomials.getTupleValue(tup);
                if(fpoly!=null){//tup has moving DOFs
                    //put in the difference in feature values for each movable feature
                    int numLinTerms = fpoly.linTermDOFs.size();
                    int numQuadTerms = fpoly.quadTermDOFs.size();
                    for(int linTermNum=0; linTermNum<numLinTerms; linTermNum++){
                        int dof1 = fpoly.linTermDOFs.get(linTermNum);
                        double x = samples.get(s).contDOFVals1[dof1];
                        double y = samples.get(s).contDOFVals2[dof1];
                        directOp.setEntry(s, fpoly.offset+linTermNum, x-y);
                    }
                    for(int quadTermNum=0; quadTermNum<numQuadTerms; quadTermNum++){
                        int qd[] = fpoly.quadTermDOFs.get(quadTermNum);
                        double x1 = samples.get(s).contDOFVals1[qd[0]];
                        double x2 = samples.get(s).contDOFVals1[qd[1]];
                        double y1 = samples.get(s).contDOFVals2[qd[0]];
                        double y2 = samples.get(s).contDOFVals2[qd[1]];
                        directOp.setEntry(s, fpoly.offset+numLinTerms+quadTermNum, x1*x2-y1*y2);
                    }
                }
            }
        }
        
        directOperator = directOp;
        transposeOperator = (RealLinearOperator) directOp.transpose();//DEBUG!! hopefully this works
    }

    @Override
    double[] computeTargetVals(double coeffs[]) {
        double[] targetVals = new double[numSamples];
        for(int s=0; s<numSamples; s++){
            targetVals[s] = samples.get(s).targetE1 - samples.get(s).targetE2;
            ConfSample samp1 = new ConfSample(samples.get(s).vox, samples.get(s).contDOFVals1);
            ConfSample samp2 = new ConfSample(samples.get(s).vox, samples.get(s).contDOFVals1);
            targetVals[s] -= featSet.evalEnergy(samp1, coeffs, movingDOFs)//feed in particular DOFs to use...
                    - featSet.evalEnergy(samp2, coeffs, movingDOFs);
        }
        
        return targetVals;
    }

    
    static boolean hasOverlap(int[] a, int[] b){
        for(int c : a){
            for(int d : b){
                if(d==c)
                    return true;
            }
        }
        return false;
    }
    
    static boolean contains(int[] a, int b){
        for(int c : a){
            if(c==b)
                return true;
        }
        return false;
    }
    
}

*/
