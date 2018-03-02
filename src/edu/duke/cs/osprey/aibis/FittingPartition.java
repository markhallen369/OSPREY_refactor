/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.aibis;

import java.util.ArrayList;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealLinearOperator;
import org.apache.commons.math3.linear.RealVector;

/**
 *
 * This is a subset of the parameters to be fit together (rather than fitting all parameters at once)
 * along with some stuff to do the fitting
 * 
 * @author mhall44
 */
public abstract class FittingPartition {
    
    
    int numParams=-1, numSamples=-1;
    FeatureSet featSet=null;//the feature set whose parameters we're partitioning
    
    
    //The following will be computed as needed by subclass methods (null til then)
    RealLinearOperator transposeOperator=null, directOperator=null;//A^t and A respectively, where we fit Ax to b
    ArrayList<Integer> refitParamIndices = null;//indices in general coeffs of the parameters being refit

    
    //NO CONSTRUCTOR HERE...SUBCLASS CONSTRUCTORS SHOULD CONSTRUCT ALL FIELDS
    
    
      
    RealLinearOperator makeNormalOperator(){
        return compose(transposeOperator, directOperator);
    }
    
    //IMPORTANT: although the normal operator operates in the space of just the coefficients being fit,
    //the coeffs taken as argument here are the full coeffs (since other param values may be needed to compute the RHS)
    RealVector makeNormalRHS(double coeffs[]){
        double targetVals[] = computeTargetVals(coeffs);
        return transposeOperator.operate(new ArrayRealVector(targetVals));
    }
    
    

    
    
    
    public static RealLinearOperator compose(RealLinearOperator A, RealLinearOperator B){
        if(A.getColumnDimension()!=B.getRowDimension())
            throw new RuntimeException("ERROR: Dimension mismatch");
                
        return new RealLinearOperator(){
            @Override
            public int getRowDimension() {
                return A.getRowDimension();
            }

            @Override
            public int getColumnDimension() {
                return B.getColumnDimension();
            }

            @Override
            public RealVector operate(RealVector rv) throws DimensionMismatchException {
                return A.operate(B.operate(rv));
            }
            
        };
    }
    
    
    abstract double[] computeTargetVals(double coeffs[]);//prepare the target vals: (energy being fit) - (contribution from non-fit parameters and SAPE/constant stuff)            


}



