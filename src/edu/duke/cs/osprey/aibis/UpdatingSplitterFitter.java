/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.aibis;

import edu.duke.cs.osprey.confspace.RCTuple;
import java.util.ArrayList;
import java.util.Iterator;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.ConjugateGradient;
import org.apache.commons.math3.linear.RealLinearOperator;
import org.apache.commons.math3.linear.RealVector;

/**
 *
 * This fitter has two advantages over typical least-squares fitting:
 * (1) Only updating parameters close enough to the integrated DOF to change significantly
 * (2) Splitting up the fit by using derivatives and differences
 * 
 * 
 * @author mhall44
 */

/*
public class UpdatingSplitterFitter implements SLEFitter {

    ConfSampleSet samples;
    FeatureSet featSet;


    public UpdatingSplitterFitter(ConfSampleSet samples, FeatureSet featSet) {
        this.samples = samples;
        this.featSet = featSet;
    }
    
    
    @Override
    public double[] fitCoeffs(double[] initCoeffs) {
        double[] coeffs = initCoeffs;
        
        //this ordering of re-fits will give valid fits for each stage given results of previous stages
        //but fits only manageable numbers of parameters at once
        fitSecondDerivs(coeffs);//this fits parameters that affect the second derivatives of the energy
        fitDiffs(coeffs);//this fits parameters that affect intra-voxel differences
        refitDiscreteParams(coeffs);//and now we just need to fit the discrete parameters
        
        //validate();//DEBUG!!  and repeat above steps with lower thresh if needed? 
        //or even refit bits of the energy one term at a time?
        //There is the regular validation, but validating here could be used to improve the fit
        
        return coeffs;
        //splitting into sub-methods will be helpful for profiling...
    }
    
    
    //MAJOR STEPS IN FITTING
    
    private void fitSecondDerivs(double[] coeffs){
        //first compute cross-terms between residues
        
        //DEBUG!!  this so far is intended for sidechain dihedrals (different DOF's for each res)
        //do need cross terms btw sc and bb!   try even between all dof's??  see what's fast & accurate
        
        for(int pos1=0; pos1<featSet.confSpace.numPos; pos1++){
            for(int pos2=0; pos2<pos1; pos2++){
                for(FittingPartition fp : CrossDerivFittingPartition.makeSidechainPairCrossDerivs(featSet, samples, pos1, pos2)){
                    refitParams(fp, coeffs);
                }
            }
        } //DEBUG!! EASY VERSION NO UPDATE
        
        /*Iterator<int[]> pairEnumerator = resPairEnumerator();//enumerates in order away from integrated DOF
        while(pairEnumerator.hasNext()){
            int[] resPair = pairEnumerator.next();
            double oldCoeffs[] = coeffs.clone();
            for(FittingPartition fp : CrossDerivFittingPartition.makeSidechainPairCrossDerivs(featSet, samples, resPair[0], resPair[1]))
                refitParams(fp, coeffs);
            double dev = coeffs - oldCoeffs;
            if(dev<pairEnumerator.devThresh)
                pairEnumerator.reportPairLowDev();//this will affect what other pairs it might want to return
        }*/
   /* }
    
    
    private void fitDiffs(double[] coeffs){
        //now fit differences between confs only differing for a single residue's sidechain dihedrals
        //DEBUG!! will need to adjust this for bb stuff.
        //in general need a partition of DOFs such that intra-partition DOF terms are handled here,
        //and inter-partition cross terms handled by the cross derivs func
        
        int numRCAtPos[] = featSet.confSpace.getNumRCsAtPos();
        
        for(int pos=0; pos<featSet.confSpace.numPos; pos++){
            fitDiffsForPos(pos, coeffs, numRCAtPos);
        } //DEBUG!! EASY VERSION NO UPDATE
        
        
        /*Iterator<Integer> singleEnumerator = resEnumerator();//enumerates in order away from integrated DOF
        while(singleEnumerator.hasNext()){
            int pos = singleEnumerator.next();
            double oldCoeffs[] = coeffs.clone();
            fitDiffsForPos(pos, coeffs, numRCAtPos);
            double dev = coeffs - oldCoeffs;
            if(dev<singleEnumerator.devThresh)
                singleEnumerator.reportLowDev();//this will affect what other pairs it might want to return
        }*/
  /*  }
    
    
    private void fitDiffsForPos(int pos, double[] coeffs, int numRCAtPos[]){
        for(int rc=0; rc<numRCAtPos[pos]; rc++){
            RCTuple single = new RCTuple(pos,rc);
            DenseFeatureSet fs = featSet.featMatrix.getOneBody(pos,rc);
            if(fs instanceof DenseFeatureSet.EPIC){//we can refit
                int dofIndices[] = ((DenseFeatureSet.EPIC)fs).dofIndices;
                FittingPartition fp = new DiffFittingPartition(featSet, samples, dofIndices, single);
                refitParams(fp, coeffs);
            }
        }
    }
    
    
    private void refitDiscreteParams(double[] coeffs){
        //DEBUG!!  could use a way to fit less of the discrete params for update purposes
        
        //OK at this point all continuous terms should be fit properly.  
        //So refit the discrete terms assuming the continuous terms are right.
        DiscreteParamsFittingPartition dpfp = new DiscreteParamsFittingPartition(featSet, samples);
        refitParams(dpfp, coeffs);
    }
        

        //NOTE: energy caching within samples (e.g. what is the energy and what model was it for)
        //will also be key
        //but hopefully can figure out a little more about energy models first, decide how to cache
                
        //Cross-validate and if needed can refit with lower thresh?  Or even refit 
        //bits of the energy one at a time?
    
    
    
    
    
    private void refitParams(FittingPartition fp, double[] coeffs){
        //before this, prepare the partition, including everything needed to apply the oeprator
        
        RealLinearOperator AtA = fp.makeNormalOperator();
        RealVector Atb = fp.makeNormalRHS(coeffs);
        ArrayList<Integer> refitParamIndices = fp.refitParamIndices;
        
        double initCoeffs[] = new double[fp.numParams];
        for(int p=0; p<fp.numParams; p++)
            initCoeffs[p] = coeffs[refitParamIndices.get(p)];
        
        double refitParamVals[] = conjugateGradientSoln(AtA, Atb, initCoeffs);
        for(int p=0; p<fp.numParams; p++)
            coeffs[refitParamIndices.get(p)] = refitParamVals[p];
    }
    
    
    static double[] conjugateGradientSoln(RealLinearOperator M, RealVector v, double[] initCoeffs){
        ConjugateGradient cg = new ConjugateGradient(100000,1e-8,false);//max_iter; delta; whether to check pos def
        //delta is target ratio of residual norm to true vals norm
        long startTime = System.currentTimeMillis();
        RealVector ans = cg.solve(M, v, new ArrayRealVector(initCoeffs));
        System.out.println( "Conjugate gradient fitting time (ms): " + (System.currentTimeMillis()-startTime) );
        return ans.toArray();    
    }
    
}
*/