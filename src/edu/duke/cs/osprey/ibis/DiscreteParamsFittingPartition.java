/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ibis;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealLinearOperator;
import org.apache.commons.math3.linear.RealVector;

/**
 *
 * A FittingPartition for fitting all the discrete parameters, assuming the continuous ones are correct
 * //DEBUG!!  Will probably want to support leaving out some (those farther from integ dof) in the case of larger designs
 * 
 * Speedup compared to fitting the same samples in a big fit comes from the lack of continuous evaluations within the fitting iterations
 * and possibly fast computation of the integrated energies to fit to using caching?
 * 
 * @author mhall44
 */
public class DiscreteParamsFittingPartition extends FittingPartition {
    //DEBUG!!  This is based on the idea (valid for Const and EPIC DenseFeatureSet)
    //that features are discrete (i.e. const over a voxel) iff they are the first in their DenseFeatureSet
    //should probably keep this up even with other features--it's cheap to include a constant offset
    
    ArrayList<ConfSample> samples;//since the goal is to cover different discrete parameters, batching is not helpful here
    ArrayList<Double> targetEnergies;
    
    TupleMatrixGeneric<Integer> tupleIndexMatrix;
    //maps tuples to their index in featSet.unprunedFeatTuples,
    //which is the index of their discrete coeff among the params being fitted
    //and of course numParams will be the size of unprunedFeatTuples
    
    static double regularizationLambda = 1e-6;
    
    
    public DiscreteParamsFittingPartition(FeatureSet featSet, ConfSampleSet sampleSet){
        this.featSet = featSet;
        
        //first figure out the parameters
        tupleIndexMatrix = new TupleMatrixGeneric<>(featSet.confSpace,  Double.POSITIVE_INFINITY, null);
        ArrayList<RCTuple> tupList = featSet.unprunedFeatTuples();
        numParams = tupList.size();
        refitParamIndices = new ArrayList<>();
        for(int p=0; p<numParams; p++){
            RCTuple tup = tupList.get(p);
            tupleIndexMatrix.setTupleValue(tup, p);
            refitParamIndices.add( featSet.featOffsets.getTupleValue(tup) );
        }
        
        //now grab the samples.  Batches in ConfSampleSet are meant to cover the discrete space
        //so just use those
        samples = new ArrayList<>();
        Random rand = new Random();
        targetEnergies = new ArrayList<>();
        for(ConfSampleSet.SampleBatch batch : sampleSet.batches){
            if( ! batch.contDOFValsArr.isEmpty()){
                int indexInBatch = rand.nextInt(batch.energies.size());
                samples.add(new ConfSample(batch.vox, batch.contDOFValsArr.get(indexInBatch)));
                targetEnergies.add(batch.energies.get(indexInBatch));
            }
        }
        numSamples = samples.size();
        
        directOperator = makeDirectOperator();
        transposeOperator = makeTransposeOperator();
    }
    
    private RealLinearOperator makeDirectOperator() {//for the relevant parameters, this is just like in traditional CG fitting
        return new RealLinearOperator() {
            @Override
            public int getRowDimension() {
                return numSamples;
            }

            @Override
            public int getColumnDimension() {
                return numParams;
            }

            @Override
            public RealVector operate(RealVector rv) throws DimensionMismatchException {
                double Arv[] = new double[numSamples];
                    
                for(int s=0; s<numSamples; s++){
                    for(RCTuple tup : featSet.tuplesForSamp(samples.get(s))){//DEBUG!! might need to speed this up by explicit looping?
                        Arv[s] += rv.getEntry(tupleIndexMatrix.getTupleValue(tup));
                    }
                }
                
                return new ArrayRealVector(Arv);
            }
        };
    }
    
    private RealLinearOperator makeTransposeOperator(){
        return new RealLinearOperator() {
            @Override
            public int getRowDimension() {
                return numParams;
            }

            @Override
            public int getColumnDimension() {
                return numSamples;
            }

            @Override
            public RealVector operate(RealVector rv) throws DimensionMismatchException {
                double[] ans = new double[numParams];
                
                for(int s=0; s<numSamples; s++){
                    double sampVal = rv.getEntry(s);
                    for(RCTuple tup : featSet.tuplesForSamp(samples.get(s))){//might need to speed this up by explicit looping?
                        ans[tupleIndexMatrix.getTupleValue(tup)] += sampVal;
                    }
                }
                
                return new ArrayRealVector(ans);
            }
        };
    }

    
    
    @Override
    double[] computeTargetVals(double coeffs[]) {
        double[] targetVals = new double[numSamples];
        
        double contOnlyCoeffs[] = coeffs.clone();
        
        for(RCTuple tup : featSet.unprunedFeatTuples()){
            contOnlyCoeffs[featSet.featOffsets.getTupleValue(tup)] = 0;
        }
        
        for(int s=0; s<numSamples; s++)
            targetVals[s] = targetEnergies.get(s) - featSet.evalEnergy(samples.get(s), contOnlyCoeffs);
        
        return targetVals;
        //DEBUG!! We can make a class similar to this for fitting any set of params!!  but the operator will be different (explicit Sparse Matrix)
        //less params
        //used to refine fit perhaps
    }

    
    
    
    @Override
    RealLinearOperator makeNormalOperator(){
        if(transposeOperator.getColumnDimension()!=directOperator.getRowDimension())
            throw new RuntimeException("ERROR: Dimension mismatch");
                
        return new RealLinearOperator(){
            @Override
            public int getRowDimension() {
                return transposeOperator.getRowDimension();
            }

            @Override
            public int getColumnDimension() {
                return directOperator.getColumnDimension();
            }

            @Override
            public RealVector operate(RealVector rv) throws DimensionMismatchException {
                RealVector ans = transposeOperator.operate(directOperator.operate(rv));
                for(int p=0; p<ans.getDimension(); p++)
                    ans.addToEntry(p, regularizationLambda*rv.getEntry(p));
                return ans;
            }
        };
    }
    
    
    
}
