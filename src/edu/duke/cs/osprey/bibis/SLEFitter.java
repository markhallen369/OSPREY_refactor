/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.bibis;

import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.ConjugateGradient;
import org.apache.commons.math3.linear.RealLinearOperator;
import org.apache.commons.math3.linear.RealVector;

/**
 *
 * fits coefficients for a SparseLinearEnergy
 * Thus implementations are linear fitters
 * 
 * @author mhall44
 */

    

public interface SLEFitter extends Serializable {
    
    double[] fitCoeffs(double[] initCoeffs);


    class BasicCGFitter implements SLEFitter {
        //most basic SLEFitter
        //Intended for a BasicConfSampleSet

        FeatureSet featSet;
        BasicConfSampleSet trainingSamples;

        public BasicCGFitter(FeatureSet featSet, BasicConfSampleSet trainingSamples) {
            this.featSet = featSet;
            this.trainingSamples = trainingSamples;
            if(featSet.freezeSAPE)
                throw new RuntimeException("ERROR: Melted SAPE expected");
        }
        
        
        @Override
        public double[] fitCoeffs(double[] initCoeffs) {
            ConjugateGradient cg = new ConjugateGradient(100000,1e-8,false);//max_iter; delta; whether to check pos def
            //delta is target ratio of residual norm to true vals norm

            int numSamp = trainingSamples.samples.size();
            
            //set up the normal equations
            RealLinearOperator AtA = new RealLinearOperator() {
                @Override
                public int getRowDimension() {return featSet.numFeatures;}
                
                @Override
                public int getColumnDimension() {return featSet.numFeatures;}

                @Override
                public RealVector operate(RealVector rv) throws DimensionMismatchException {
                    //first apply A
                    double Arv[] = new double[numSamp];
                    for (int s = 0; s < numSamp; s++)
                        Arv[s] = featSet.evalEnergy(trainingSamples.samples.get(s), rv.toArray());

                    //then apply A^T to Arv
                    double ans[] = new double[featSet.numFeatures];
                    for (int s = 0; s < numSamp; s++)
                        featSet.updateATProduct(trainingSamples.samples.get(s), Arv[s], ans);
                    
                    return new ArrayRealVector(ans, false);//make RealVector without copying ans
                }
            };
            
            double atb[] = new double[featSet.numFeatures];
            //apply A^T to true vals
            for (int s = 0; s < numSamp; s++) {
                featSet.updateATProduct(trainingSamples.samples.get(s), trainingSamples.energies.get(s), atb);
            }
            RealVector Atb = new ArrayRealVector(atb);
            
            
            long startTime = System.currentTimeMillis();
            RealVector ans = cg.solve(AtA, Atb, new ArrayRealVector(initCoeffs));

            System.out.println( "Conjugate gradient fitting time (ms): " + (System.currentTimeMillis()-startTime) );

            return ans.toArray();
        }
    }

    public static class ContBatchCGFitter implements SLEFitter {
        //intended for a ContBatchSampleSet
        
        FeatureSet featSet;
        ContBatchSampleSet trainingSamples;

        public ContBatchCGFitter(FeatureSet featSet, ContBatchSampleSet trainingSamples) {
            this.featSet = featSet;
            this.trainingSamples = trainingSamples;
            if(featSet.freezeSAPE)
                throw new RuntimeException("ERROR: Melted SAPE expected");
        }
        
        @Override
        public double[] fitCoeffs(double[] initCoeffs) {
            ConjugateGradient cg = new ConjugateGradient(100000,1e-8,false);//max_iter; delta; whether to check pos def
            //delta is target ratio of residual norm to true vals norm

            int numSamp = trainingSamples.numSamples;
            
            //set up the normal equations
            RealLinearOperator AtA = new RealLinearOperator() {
                @Override
                public int getRowDimension() {return featSet.numFeatures;}
                
                @Override
                public int getColumnDimension() {return featSet.numFeatures;}

                @Override
                public RealVector operate(RealVector rv) throws DimensionMismatchException {
                    //first apply A
                    double Arv[] = new double[numSamp];
                    
                    int s=0;
                    for(ContBatchSampleSet.SampleBatch batch : trainingSamples.batches){
                        MoleculeModifierAndScorer mms = featSet.makeMMS(new ConfSample(batch.vox,null));
                        for(HashMap<String,Double> contDOFVals : batch.contDOFVals){
                            Arv[s] = featSet.evalEnergy(new ConfSample(batch.vox,contDOFVals), rv.toArray(), mms);
                            s++;
                        }
                    }
                    
                    //then apply A^T to Arv
                    double ans[] = new double[featSet.numFeatures];
                    s=0;
                    for(ContBatchSampleSet.SampleBatch batch : trainingSamples.batches){
                        MoleculeModifierAndScorer mms = featSet.makeMMS(new ConfSample(batch.vox,null));
                        for(HashMap<String,Double> contDOFVals : batch.contDOFVals){
                            featSet.updateATProduct(new ConfSample(batch.vox,contDOFVals), Arv[s], ans, mms);
                            s++;
                        }
                    }
                    
                    return new ArrayRealVector(ans, false);//make RealVector without copying ans
                }
            };
            
            double atb[] = new double[featSet.numFeatures];
            //apply A^T to true vals
            for(ContBatchSampleSet.SampleBatch batch : trainingSamples.batches){
                MoleculeModifierAndScorer mms = featSet.makeMMS(new ConfSample(batch.vox,null));
                for(int index=0; index<batch.contDOFVals.size(); index++){
                    HashMap<String,Double> contDOFVals = batch.contDOFVals.get(index);
                    double sampE = batch.energies.get(index);
                    featSet.updateATProduct(new ConfSample(batch.vox,contDOFVals), sampE, atb, mms);
                }
            }

            RealVector Atb = new ArrayRealVector(atb);
            
            
            long startTime = System.currentTimeMillis();
            RealVector ans = cg.solve(AtA, Atb, new ArrayRealVector(initCoeffs));

            System.out.println( "Conjugate gradient fitting time (ms): " + (System.currentTimeMillis()-startTime) );

            return ans.toArray();
        }
    }
    
    
    public static class FrozenSAPECGFitter implements SLEFitter {
        //intended for a ContBatchSampleSet
        //with frozen SAPE
        
        FeatureSet featSet;
        ContBatchSampleSet trainingSamples;

        public FrozenSAPECGFitter(FeatureSet featSet, ContBatchSampleSet trainingSamples) {
            this.featSet = featSet;
            this.trainingSamples = trainingSamples;
            if(!featSet.freezeSAPE)
                throw new RuntimeException("ERROR: Frozen SAPE expected");
        }
        
        @Override
        public double[] fitCoeffs(double[] initCoeffs) {
            ConjugateGradient cg = new ConjugateGradient(100000,1e-8,false);//max_iter; delta; whether to check pos def
            //delta is target ratio of residual norm to true vals norm

            int numSamp = trainingSamples.numSamples;
            
            //precalculate coeff-independent enegies
            double frozenE[] = new double[numSamp];
            if(featSet.hasSAPE){
                int s = 0;
                for(ContBatchSampleSet.SampleBatch batch : trainingSamples.batches){
                    MoleculeModifierAndScorer mms = featSet.makeMMS(new ConfSample(batch.vox,null));
                    for(int index=0; index<batch.contDOFVals.size(); index++){
                        frozenE[s] = featSet.evalFrozenEnergy(new ConfSample(batch.vox,batch.contDOFVals.get(index),batch.contDOFValsArr.get(index)), mms);
                        s++;
                    }
                }
            }
            
            //set up the normal equations
            RealLinearOperator AtA = new RealLinearOperator() {
                @Override
                public int getRowDimension() {return featSet.numFeatures;}
                
                @Override
                public int getColumnDimension() {return featSet.numFeatures;}

                @Override
                public RealVector operate(RealVector rv) throws DimensionMismatchException {
                    //first apply A
                    double Arv[] = new double[numSamp];
                                        
                    int s=0;
                    for(ContBatchSampleSet.SampleBatch batch : trainingSamples.batches){
                        //MoleculeModifierAndScorer mms = featSet.makeMMS(new ConfSample(batch.vox,null));
                        for(int index=0; index<batch.contDOFVals.size(); index++){
                            Arv[s] = featSet.evalNoSAPEEnergy(new ConfSample(batch.vox,batch.contDOFVals.get(index),batch.contDOFValsArr.get(index)), rv.toArray());
                            s++;
                        }
                    }
                    
                    //then apply A^T to Arv
                    double ans[] = new double[featSet.numFeatures];
                    s=0;
                    for(ContBatchSampleSet.SampleBatch batch : trainingSamples.batches){
                        //MoleculeModifierAndScorer mms = featSet.makeMMS(new ConfSample(batch.vox,null));
                        for(int index=0; index<batch.contDOFVals.size(); index++){
                        //for(HashMap<String,Double> contDOFVals : batch.contDOFVals){
                            featSet.updateATProduct(new ConfSample(batch.vox,batch.contDOFVals.get(index),batch.contDOFValsArr.get(index)), Arv[s], ans);
                            s++;
                        }
                    }
                    
                    return new ArrayRealVector(ans, false);//make RealVector without copying ans
                }
            };
            
            double atb[] = new double[featSet.numFeatures];
            //apply A^T to true vals
            int s=0;
            for(ContBatchSampleSet.SampleBatch batch : trainingSamples.batches){
                //MoleculeModifierAndScorer mms = featSet.makeMMS(new ConfSample(batch.vox,null));
                for(int index=0; index<batch.contDOFVals.size(); index++){
                    HashMap<String,Double> contDOFVals = batch.contDOFVals.get(index);
                    double sampE = batch.energies.get(index);
                    featSet.updateATProduct(new ConfSample(batch.vox,contDOFVals,batch.contDOFValsArr.get(index)), sampE-frozenE[s], atb);
                    s++;
                }
            }

            RealVector Atb = new ArrayRealVector(atb);
            
            
            long startTime = System.currentTimeMillis();
            
            RealVector ans;
            if(initCoeffs==null)
                ans = cg.solve(AtA, Atb);
            else
                ans = cg.solve(AtA, Atb, new ArrayRealVector(initCoeffs));

            System.out.println( "Conjugate gradient fitting time (ms): " + (System.currentTimeMillis()-startTime) );

            return ans.toArray();
        }
    }
    
}
