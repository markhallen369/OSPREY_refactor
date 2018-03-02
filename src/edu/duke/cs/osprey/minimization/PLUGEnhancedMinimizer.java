/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import com.joptimizer.functions.LinearMultivariateRealFunction;
import com.joptimizer.functions.PSDQuadraticMultivariateRealFunction;
import com.joptimizer.optimizers.JOptimizer;
import com.joptimizer.optimizers.OptimizationRequest;
import edu.duke.cs.osprey.confspace.RCTuple;
import static edu.duke.cs.osprey.minimization.SQPMinimizer.qpTol;
import edu.duke.cs.osprey.plug.LPChecks;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import java.util.ArrayList;
import java.util.Random;

/**
 *
 * This minimizer is enhanced in two ways:
 * -It uses PLUG to get a declashed initial point (combine this with initFixableDOFs),
 * ensuring that all minimized energies are downhill from something clash-free
 * -It can use multiple draws with perturbed starting points to ensure the minimum
 * is robust to small perturbations
 * 
 * 
 * @author mhall44
 */
public class PLUGEnhancedMinimizer implements Minimizer {
    
    MoleculeModifierAndScorer mms;
    LinearMultivariateRealFunction[] polytope;

    //To keep the scoring from depending too much on small perturbations
    //to the minimization problem,
    //we can rerun it with slightly perturbed starting points
    //until we come up with numConsistent draws within consistencyThresh
    public static double consistencyThresh = 0.01;
    int numConsistent = 1;
    public static double perturbationSigma = 0.1;//this is per-DOF and relative to voxel width
    
    static double polytopeTol = 0.001;//tolerance on the polytope constr
    //since we have to be able to minimize anything that passed PLUG pruning
    
    public PLUGEnhancedMinimizer(MoleculeModifierAndScorer objFcn, 
            PolytopeMatrix plugMat, RCTuple RCTup, int numConsistent){
        mms = objFcn;
        this.numConsistent = numConsistent;
        polytope = plugMat.getJoptPolytope(RCTup, objFcn.getDOFs(), polytopeTol);//DEBUG!!  assuming this has box constr
    }
    
    
    @Override
    public Result minimize() {
        CCDMinimizer baseMinimizer = new CCDMinimizer(mms,false);
        DoubleMatrix1D baseInitVals = plugComputeInitVals();
        if(numConsistent==1){
            baseMinimizer.singleInitVal = keepInRange(baseInitVals);
            return baseMinimizer.minimize();
        }
        else {
            ArrayList<Result> baseResults = new ArrayList<>();
            while(true){
                baseMinimizer.singleInitVal = keepInRange(randomlyPerturb(baseInitVals));
                baseResults.add(baseMinimizer.minimize());
                for(int s=0; s<baseResults.size(); s++){
                    int numAgreeing=0;//see how many other scores agree with score s
                    for(int s2=0; s2<s; s2++){
                        if( Math.abs(baseResults.get(s).energy-baseResults.get(s2).energy) 
                                < consistencyThresh ){
                            numAgreeing++;
                        }
                    }
                    if(numAgreeing+1>=numConsistent)
                        return baseResults.get(s);
                }
                
                if(baseResults.size()%(3*numConsistent) == 0){
                    System.out.println("Warning: minimized energy with "+baseResults.size()
                        +" different random initial values but no "+numConsistent+
                        " minima agreed within "+consistencyThresh);
                }
            }
        }
    }
    
    
    private DoubleMatrix1D keepInRange(DoubleMatrix1D x){
        //Move any elements of x that are outside the box constraints back in
        //this is for making initial values
        x = x.copy();
        double tol = 1e-10;//let's try to get them in range by at least tol
        //won't affect starting point quality noticeably, but can help with in-range checks, etc.
        DoubleMatrix1D[] constr = mms.getConstraints();
        for(int dof=0; dof<constr[0].size(); dof++){
            double ub = constr[1].get(dof)-tol;
            double lb = constr[0].get(dof)+tol;
            if(x.get(dof)>ub)
                x.set(dof,ub);
            else if(x.get(dof)<lb)
                x.set(dof,lb);
        }
        return x;
    }
    
    
    private DoubleMatrix1D plugComputeInitVals(){
        //minimize distance from center of mms constr, subject to plug constr
        DoubleMatrix1D top = mms.getConstraints()[1];
        DoubleMatrix1D bottom = mms.getConstraints()[0];
        DoubleMatrix1D center = top.copy().assign(bottom,Functions.plus).assign(Functions.mult(0.5));
        
        //Need to do a convex quadratic optimization: minimize ||x-center||^2
        //subject to keeping x in polytope
        int numDOFs = mms.getNumDOFs();
        double[][] P = new double[numDOFs][numDOFs];
        double[] q = new double[numDOFs];
        for(int i=0; i<numDOFs; i++){
            P[i][i] = 1.;
            q[i] = -center.get(i);
        }
        PSDQuadraticMultivariateRealFunction objFcn = new PSDQuadraticMultivariateRealFunction(P,q,0);
        
        OptimizationRequest or = new OptimizationRequest();
        or.setF0(objFcn);
        
        //or.setInitialPoint(center.toArray());
        DoubleMatrix1D interiorPt = LPChecks.getInteriorPt(polytope);
        if(interiorPt==null){
            throw new RuntimeException("ERROR: PLUGEnhancedMinimizer asked to minimize in "
                        + "a PLUG-infeasible voxel");
        }
        or.setInitialPoint(interiorPt.toArray());
        
        //hoping to improve efficiency
        or.setToleranceFeas(0.0001);
        or.setToleranceInnerStep(0.0001);
        or.setToleranceKKT(0.0001);
        
        or.setTolerance(qpTol);
        or.setFi(polytope);
        
        JOptimizer opt = new JOptimizer();
	opt.setOptimizationRequest(or);
        try{
            opt.optimize();
        }
        catch(Exception e){
            e.printStackTrace();
            throw new RuntimeException(e.getMessage());
        }

        double[] ans = opt.getOptimizationResponse().getSolution();
        return DoubleFactory1D.dense.make(ans);
    }
    
    
    private DoubleMatrix1D randomlyPerturb(DoubleMatrix1D x){
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(x.size());
        Random rand = new Random();
        DoubleMatrix1D voxWidths = mms.getConstraints()[1].copy().assign(mms.getConstraints()[0],Functions.minus);
        for(int dof=0; dof<x.size(); dof++){
            ans.set(dof, x.get(dof) + voxWidths.get(dof)*perturbationSigma*rand.nextGaussian() );
        }
        return ans;
    }
    
    
    
    
    
}
