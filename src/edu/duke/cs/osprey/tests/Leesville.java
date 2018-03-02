/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import cern.colt.matrix.DoubleFactory1D;
import com.joptimizer.functions.LinearMultivariateRealFunction;
import com.joptimizer.functions.PSDQuadraticMultivariateRealFunction;
import com.joptimizer.optimizers.JOptimizer;
import com.joptimizer.optimizers.OptimizationRequest;

/**
 *
 * Quadratic programming test
 * 
 * 
 * @author mhall44
 */
public class Leesville {
    
    public static void main(String[] args){
        double P[][] = new double[][] { {1,0}, {0,1} };
        double q[] = new double[] {3,5};
        PSDQuadraticMultivariateRealFunction objFcn = new PSDQuadraticMultivariateRealFunction(P,q,0);
        
        OptimizationRequest or = new OptimizationRequest();
        or.setF0(objFcn);
        
        or.setInitialPoint(new double[] {4,6});//this is funny because it's really 3,5!!!  LOL

        //hoping to improve efficiency
        or.setToleranceFeas(0.0001);
        or.setToleranceInnerStep(0.0001);
        or.setToleranceKKT(0.0001);
        
        double qpTol = 1e-5;
        or.setTolerance(qpTol);
        LinearMultivariateRealFunction[] box = new LinearMultivariateRealFunction[] {
            new LinearMultivariateRealFunction(new double[] {-1,0},-1000),
            new LinearMultivariateRealFunction(new double[] {1,0},-1000),
            new LinearMultivariateRealFunction(new double[] {0,-1},-1000),
            new LinearMultivariateRealFunction(new double[] {0,1},-1000)
        };

        or.setFi(box);
        
        JOptimizer opt = new JOptimizer();
	opt.setOptimizationRequest(or);
        try{
            opt.optimize();
        }
        catch(Exception e){
            if(e.getMessage().contains("initial point must be strictly feasible")){
                //return null;
                throw new RuntimeException("ERROR: PLUGEnhancedMinimizer asked to minimize in "
                        + "a PLUG-infeasible voxel");
            }
            e.printStackTrace();
            throw new RuntimeException(e.getMessage());
        }

        double[] ans = opt.getOptimizationResponse().getSolution();
        System.out.println("Best x: "+DoubleFactory1D.dense.make(ans));
    }
    
}
