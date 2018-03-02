/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.jet.math.Functions;
import com.joptimizer.functions.LinearMultivariateRealFunction;
import com.joptimizer.functions.PSDQuadraticMultivariateRealFunction;
import com.joptimizer.optimizers.JOptimizer;
import com.joptimizer.optimizers.OptimizationRequest;
import edu.duke.cs.osprey.energy.derivs.EnergyDerivs;
import edu.duke.cs.osprey.plug.LPChecks;
import edu.duke.cs.osprey.tools.ObjectIO;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.BitSet;

/**
 *
 * @author mhall44
 * Trying sequential quadratic programming, including the SQP+ algorithm of Morales et al
 * Constraints will be linear, which simplifies the algorithm
 * The biggest issue is that the Hessian may often be nonconvex, and SQP+ is supposed to help in such cases
 * 
 */
public class SQPMinimizer implements Minimizer, Serializable {

    
    /*private static class SerializableLinearMultivariateRealFunction extends LinearMultivariateRealFunction implements Serializable {
        SerializableLinearMultivariateRealFunction(LinearMultivariateRealFunction lmrf){
            super(lmrf.getQ().toArray(), lmrf.getR());
        }
        SerializableLinearMultivariateRealFunction(double[] q, double r){
            super(q,r);
        }
    }//DEBUG!!!  just want to be able to serialize this minimizer for debugging*/
    
    
    transient LinearMultivariateRealFunction[] constr;//linear constraints defining the feasible region (in f(x)<=0 form)
    int numConstr;
    ObjectiveFunction objFcn;
        
    
    public DoubleMatrix1D x;//current estimate of best DOF vals
    //DEBUG!!!  (on publicness)
    
    
    static final boolean SQPPlus = true;//use SQP+ (w/ equality constraint step+line search)
    static final double backtrackingTau = 0.8;//for line search
    static final double sigma = 0.5;//for sufficent decrease condition in line search
    
    //tolerances for computing things 
    static final double numDiffInterval = 1e-3;//interval for computing finite differences
    static final double qpTol = 1e-5;//convergence tolerance for quadratic programming
    //tolerances for evaluating things--higher so that can tolerate the errors from computing things
    static final double constrViolTol = 1e-3;
    static final double distConvTol = 0.01;//evaluating convergence
    static final double improvementTol = 1e-7;//if convex optimization step only 
    //gives this much improvement, gradient in directions not blocked by constraints must be negligible
    static final double redundancyTol = 1e-6;//for checking redundancy of constraints
    
    static final int maxNumIter=100;
    
    double minTime;
    
    
    public ArrayList<DoubleMatrix1D> xTrack = new ArrayList<>();//DEBUG!!!
    
    
    double[] boxCenter;//DEBUG!!!
    
    public SQPMinimizer(){}
    
    public SQPMinimizer(ObjectiveFunction of, LinearMultivariateRealFunction[] extraConstr){
        //This minimizer is intended for an ObjectiveFunction whose constraints won't change!
        //For example a MoleculeModifierAndScorer
        
        objFcn = of;
        
        
        //build constraints
        if(extraConstr==null)//no constraints
            extraConstr = new LinearMultivariateRealFunction[0];
        int n = of.getNumDOFs();
        
        
        //DEBUG!!!!  removing redundant constr
        ArrayList<LinearMultivariateRealFunction> ajlnv = new ArrayList<>();
        for(LinearMultivariateRealFunction f : extraConstr){
            boolean redundant = false;
            for(LinearMultivariateRealFunction g : ajlnv){
                if(Math.abs(f.getR()-g.getR())<redundancyTol){
                    if(Algebra.DEFAULT.norm2(f.getQ().copy().assign(g.getQ(),Functions.minus))<redundancyTol*redundancyTol){
                        redundant = true;
                        break;
                    }
                }
            }
            if(!redundant)
                ajlnv.add(f);
        }
        constr = ajlnv.toArray(new LinearMultivariateRealFunction[ajlnv.size()]);
        
        
        //DEBUG!!!  box constr already provided if extraConstr nonempty...
        if(constr.length==0){
            constr = new LinearMultivariateRealFunction[2*n+extraConstr.length];
            DoubleMatrix1D DOFBounds[] = of.getConstraints();
            for(int dim=0; dim<n; dim++){
                double coeffsLB[] = new double[n];
                coeffsLB[dim] = -1;
                constr[2*dim] = new LinearMultivariateRealFunction(coeffsLB, DOFBounds[0].get(dim));
                double coeffsUB[] = new double[n];
                coeffsUB[dim] = 1;
                constr[2*dim+1] = new LinearMultivariateRealFunction(coeffsUB, -DOFBounds[1].get(dim));
            }
            System.arraycopy(extraConstr, 0, constr, 2*n, extraConstr.length);
        }
        
        numConstr = constr.length;
    }
    
    
    @Override
    public Minimizer.Result minimize() {
        //First figure out the constraints and initial values
        //(since the minimizer might be used for many rotameric states, etc.,
        //this can't be done in the constructor)
        if(!compInitVals())//Initialize x
            return null;//No initial values (constraints are in conflict)
        
        minimizeFromCurPoint();
        if(x==null)//constraints too close to conflict to minimize
            return null;
        
        return new Minimizer.Result(x, objFcn.getValue(x));
    }
    
    
    

    public void minimizeFromCurPoint(){//minimize from current starting vector
        
        long minStartTime = System.currentTimeMillis();
        
        //double oldE = Double.POSITIVE_INFINITY;                           
        objFcn.setDOFs(x);
                
        for(int iter=0; iter<maxNumIter; iter++) {
            
            xTrack.add(x.copy());//DEBUG!!!
            if(iter==90){//put at 90 so last 10 steps can be compared to what happens here
                ObjectIO.writeObject(this, "SQP_MINIMIZER_STEP90.dat");
            }
            
            
            double E = objFcn.getValue(x);
            //System.out.println("Iteration: "+iter+" Energy: "+E);
            /*if(isConverged(E,oldE))
                break;
            oldE = E;*/
            DoubleMatrix1D oldX = x.copy();

            
            //OK now make a quadratic approximation and optimize it
            //DoubleMatrix1D grad = calcGradient(x);
            //DoubleMatrix2D hess = calcHessian(x);
            EnergyDerivs ed = new EnergyDerivs((MoleculeModifierAndScorer)objFcn, x);//DEBUG!! this version of ed is for mms
            DoubleMatrix1D grad = ed.getGrad();
            DoubleMatrix2D hess = ed.getHess();
            
            
            DoubleMatrix2D B = posSemidefApprox(hess);//may also try BFGS here, in which case may not need hess at all            
            
            
            //DEBUG!!!
            /*DoubleMatrix1D altGrad = calcGradient(x);//altCalcGradient(x);
            DoubleMatrix2D altHess = calcHessian(x);//altCalcHessian(x);
            DoubleMatrix1D gradRelErr = altGrad.copy().assign(grad,Functions.minus).assign(grad,Functions.div);
            DoubleMatrix2D hessErr = altHess.copy().assign(hess,Functions.minus).assign(hess,Functions.div);
            */
            
            //the first step is to solve the convex approximate problem: minimize 0.5*(x-x0)B(x-x0) + grad*(x-x0)
            //equivalently minimize grad*x + 0.5*xBx - (Bx0)x
            DoubleMatrix1D qConvex = grad.copy().assign( Algebra.DEFAULT.mult(B,x), Functions.minus );
            double[] qpSoln = solveQP(B.toArray(), qConvex.toArray(), x.toArray());
            if(qpSoln==null){//feasible region too small
                x=null;
                return;
            }
            DoubleMatrix1D convProblemSoln = DoubleFactory1D.dense.make( qpSoln );
            DoubleMatrix1D convProblemStep = convProblemSoln.copy().assign(x, Functions.minus);
            
            //calculate approximate improvement from ineq step
            double ineqApproxImprovement = - ( convProblemStep.zDotProduct(grad)
                    + convProblemStep.zDotProduct(Algebra.DEFAULT.mult(B,convProblemStep))/2 );

            if(ineqApproxImprovement<improvementTol){
                //this only happens if gradient (in allowed directions) is within numerical error of 0
                //thus indicates a local minimum at x
                //System.out.println("SQP converged by negligible improvement in "+iter+" iterations");
                break;
            }

            if(SQPPlus){
                BitSet activeConstr = pickActiveConstr(convProblemSoln);
                DoubleMatrix1D adjGrad = grad.copy().assign( Algebra.DEFAULT.mult(hess,x), Functions.minus );
                DoubleMatrix1D eqProblemSoln = minQuadraticOnManifold(hess, adjGrad, activeConstr);
                
                boolean eqUseful;
                DoubleMatrix1D eqProblemStep = null;
                if(eqProblemSoln==null){//couldn't minimize quadratic on manifold
                    //probably means no unique solution.  Won't use it.  (Should be rare, gives a warning).
                    eqUseful = false;
                }
                else {//got a solution, see if it's useful in improving energy
                    eqProblemStep = eqProblemSoln.copy().assign(convProblemSoln, Functions.minus);

                    //In regions not convex within the hyperplane of active constraints,
                    //eqProblemStep might go to a higher-energy saddle point
                    //in this case the other extreme could be a better estimate
                    //this is mainly useful for following gently sloping valleys with deep walls
                    //(the one negative eigenvalue is along the valley)
                    //note: this complicates entropy computation.  I guess we'd need to integrate along the valley?
                    //The arbitrary cutoff at the voxel edge is a problem
                    //This would be best investigated using more accurate energy functions
                    //VDW provides an approx Rama potential but not sure how realistic this behavior is
                    //as it stands we can't count on finding convex minima within voxels
                    //so the arbitrary voxel borders are a big problem for within-voxel entropy calculations
                    double eqApproxImprovement = - ( eqProblemStep.zDotProduct(grad)
                            + eqProblemStep.zDotProduct(Algebra.DEFAULT.mult(hess,eqProblemStep))/2 );
                    if(eqApproxImprovement<0)
                        eqProblemStep.assign(Functions.mult(-1));

                    scaleWithinConstr(eqProblemStep, convProblemSoln);



                    double eqTarget = E - sigma*ineqApproxImprovement;
                    eqUseful = lineSearch(convProblemSoln, eqProblemStep, 0.25, eqTarget, 0);//line search from eqProblemStep back to adjSoln,
                    //going back to 0.25, stopping if function value less than ineqApproxImprovement
                    //scale eqProblemStep if successful, else return false
                
                    //DEBUG!!!
                    /*eqProblemSoln.assign(x);
                    eqProblemSoln.assign(eqProblemStep,Functions.plus);
                    double eqSolnE = objFcn.getValue(eqProblemSoln);
                    double eqSolnApprox = approxObjFcn(eqProblemSoln, hess, grad);
                    double convSolnE = objFcn.getValue(convProblemSoln);
                    double convSolnApprox = approxObjFcn(convProblemSoln, hess, grad);
                    double convSolnConvApprox = approxObjFcn(convProblemSoln, B, grad);*/
                    //if(iter==28){
                        //System.out.println("ConvSoln for iter "+iter+" : "+convProblemSoln);
                        //scanDirection(convProblemStep, hess, grad);*/
                    //    System.exit(0);
                    //}
                }
                                
                
                if(eqUseful){
                    x.assign(convProblemStep, Functions.plus);
                    x.assign(eqProblemStep, Functions.plus);
                }
                else {
                    //eq problem not useful.  Line search for convex step
                    lineSearch(x, convProblemStep, 0, E, sigma*ineqApproxImprovement);
                    //moving target now. If scaling by alpha, target should be E - sigma*alpha*ineqApproxImprovement
                    x.assign(convProblemStep, Functions.plus);
                }
            }
            else//basic SQP
                x.assign(convProblemStep, Functions.plus);
            //an intermediate possibility is basic + line search
            
            if(isConverged(oldX)){
                //System.out.println("SQP converged by small step in "+iter+" iterations");
                break;
            }
            
            if(iter==maxNumIter-1){
                System.out.println("Warning: SQP did not converge!");
                System.exit(0);//DEBUG!!!
            }
        }

        minTime = System.currentTimeMillis() - minStartTime;
    }
    
    
    void scanDirection(DoubleMatrix1D step, DoubleMatrix2D hess, DoubleMatrix1D grad){
        System.out.println("SCANNING DIRECTION.  fac,E,approx(2nd order)");
        for(double fac=1e-3; fac<1; fac*=1.2){
            DoubleMatrix1D z = x.copy().assign(step, Functions.plusMult(fac));
            double E = objFcn.getValue(z);
            double approx = approxObjFcn(z,hess,grad);
            System.out.println(fac+" "+E+" "+approx);
        }
        System.out.println("DONE SCANNING DIRECTION.  ");
    }
    
    double approxObjFcn(DoubleMatrix1D z, DoubleMatrix2D H, DoubleMatrix1D g){
        //evaluate the approximate energy E(x) + g*(z-x) + 0.5*(z-x)'*H*(z-x)
        DoubleMatrix1D relZ = z.copy().assign(x,Functions.minus);
        double ans = objFcn.getValue(x);
        ans += g.zDotProduct(relZ);
        ans += 0.5*relZ.zDotProduct( Algebra.DEFAULT.mult(H,relZ) );
        return ans;
    }
    
    BitSet pickActiveConstr(DoubleMatrix1D z){
        //figure out what constraints are active at z
        BitSet ans = new BitSet();
        for(int constrNum=0; constrNum<constr.length; constrNum++){
            if( constr[constrNum].value(z.toArray()) > -constrViolTol )//not significantly in interior of allowed region
                ans.set(constrNum);
        }
        return ans;
    }
        
    
    void scaleWithinConstr(DoubleMatrix1D change, DoubleMatrix1D base){
        //base is a vector satisfying the constr
        //Find maximal a (0<=a<=1) such that base+a*change satisfies the constr
        //Then make change *= a
        DoubleMatrix1D sum = base.copy().assign(change, Functions.plus);
        
        for(LinearMultivariateRealFunction constraint : constr){
            double sumConstrVal = constraint.value(sum.toArray());
            if(sumConstrVal>constrViolTol){//significant violation
                double baseConstrVal = constraint.value(base.toArray());
                if(baseConstrVal>constrViolTol)
                    throw new RuntimeException("ERROR: Base vector violates constraints");
                //Want baseConstrVal + a*(sumConstrVal-baseConstrVal) = 0
                double a = baseConstrVal / (baseConstrVal-sumConstrVal);
                change.assign(Functions.mult(a));
                sum = base.copy().assign(change, Functions.plus);
            }
        }
    }
    
    
    /*boolean isConverged(double E, double oldE){
        if(Double.isInfinite(E)&&Double.isInfinite(oldE)){
            //Infinite energy two steps in a row implies we're not getting out
            //Could also mean we started at infinity, but due to the initial value computation setup,
            //this should only happen for a uniformly infinite objective function
            return true;
        }
        else if ( oldE - E < EConvTol ) {
            //Considered to be converged if energy improvement is small enough
            //or if numerical error is causing the energy to rise
            //(for EConvTol small enough, or for very rugged regions like clashes, the rise may exceed EConvTol)
            
                        //now that there's a gradient and Hessian maybe can use those for a better convergence criterion
            return true;
        }
        else
            return false;
    }*/
    
    //let's try monitoring convergence by how far consecutive steps go.  
    //This should work regardless of energy ranges involved,
    //and if can get into convex well as expected, will be a good sign of convergence.  
    //Can perhaps compare active set to last iter as well.  
    boolean isConverged(DoubleMatrix1D oldX){
        double lastStepDistSq = Algebra.DEFAULT.norm2( x.copy().assign(oldX,Functions.minus) );
        //System.out.println("Squared distance travelled in last step: "+lastStepDistSq);
        return lastStepDistSq <= distConvTol*distConvTol;
    }
    
    
    boolean lineSearch(DoubleMatrix1D startVec, DoubleMatrix1D fullStep, double minFrac,
            double baseE, double targetImprovementSlope){
        //do a backtracking line search from startVec+fullStep to startVec
        //goal is to find z = startVec+alpha*fullStep such that E(z) <= baseE - alpha*targetImprovementSlope
        //the function will scale fullStep *= alpha to indicate what it found
        //if get to alpha<minFrac without sufficient energy improvement, bail and return false        
        for(double alpha=1; alpha>=minFrac; alpha*=backtrackingTau){
            DoubleMatrix1D z = startVec.copy().assign(fullStep, Functions.plusMult(alpha));
            double curE = objFcn.getValue(z);
            if(curE <= baseE - alpha*targetImprovementSlope){//found good improvement
                fullStep.assign(Functions.mult(alpha));
                return true;
            }
            if(alpha<improvementTol){//kind of weird...if Taylor series is good, at least
                //the target improvements should happen at small steps (on the scale of improvements in conv step, if small)
                //may be numerical error.  
                System.out.println("Warning: line search not converging at alpha="+alpha+" E-baseE"+(curE-baseE)
                        +" targetImprovementSlope: "+targetImprovementSlope);
                fullStep.assign(Functions.mult(0));
                return false;
            }
        }
        
        return false;//didn't find good improvement
    }
    
    public static boolean lowerTol = true;//DEBUG!!!
    
    
    double[] solveQP(double[][] P, double[] q, double[] initVal){
        //Minimize qx+x^tPx/2 over the problem constraints
        //Keep a set of JOptimizerConstr for this purpose

        PSDQuadraticMultivariateRealFunction approxObjFcn = new PSDQuadraticMultivariateRealFunction(P,q,0);
        
        // Inequality constraints
        //ConvexMultivariateRealFunction[] inequalities = new ConvexMultivariateRealFunction[1];

        OptimizationRequest or = new OptimizationRequest();
        or.setF0(approxObjFcn);
        
        //or.setInitialPoint(initVal);
        //DEBUG!!!!
        or.setInitialPoint(boxCenter);

        if(lowerTol){
            or.setToleranceFeas(0.0001);
            or.setToleranceInnerStep(0.0001);
            or.setToleranceKKT(0.0001);
        }
        
        or.setTolerance(qpTol);
        or.setFi(constr);
        
        JOptimizer opt = new JOptimizer();
	opt.setOptimizationRequest(or);
        try{
            opt.optimize();
        }
        catch(Exception e){
            //if averaging of opposite corners doesn't yield interior point, feasible region is tiny
            //consider it negligible
            if(e.getMessage().contains("initial point must be strictly feasible"))
                return null;
            e.printStackTrace();
            throw new RuntimeException(e.getMessage());
        }

        return opt.getOptimizationResponse().getSolution();
    }
    
    
    DoubleMatrix1D minQuadraticOnManifold(DoubleMatrix2D P, DoubleMatrix1D q, BitSet activeConstr){
        //If want Mx=v, minimize qx + 0.5*x*P*x - x*M*lambda (lambda are Lagrange multipliers)
        //taking grad wrt x, q + P*x - M^T*lambda = 0, and also of course Mx = v
        //solve this system for x and lambda.  It has the right number of eqs.  
        //Colt will be convenient here.
        int numConstr = activeConstr.cardinality();
        int n = objFcn.getNumDOFs();
        
        //prepare active constr in (Mx=v) form
        DoubleMatrix2D M = DoubleFactory2D.dense.make(numConstr,n);
        DoubleMatrix1D v = DoubleFactory1D.dense.make(numConstr);
        int overallConstrNum = -1;
        for(int constrCount=0; constrCount<numConstr; constrCount++){
            overallConstrNum = activeConstr.nextSetBit(overallConstrNum+1);
            //put in M, v
            M.viewRow(constrCount).assign(constr[overallConstrNum].getQ());
            v.set(constrCount, -constr[overallConstrNum].getR());
        }
        
        DoubleMatrix2D lhs = DoubleFactory2D.dense.make(numConstr+n, numConstr+n);
        //lhs consists of blocks for equations corresponding to x dimensions and constr,
        //columns also correspond to x dimensions or Lagrange multipliers (one for each constr)
        DoubleMatrix2D lhsxxPart = lhs.viewPart(0, 0, n, n);
        lhsxxPart.assign(P);
        
        DoubleMatrix2D lhsxcPart = lhs.viewPart(0, n, n, numConstr);
        lhsxcPart.assign(Algebra.DEFAULT.transpose(M));
        lhsxcPart.assign(Functions.mult(-1));
        
        //constraints' values do not depend on Lagrange multipliers
        lhs.viewPart(n, 0, numConstr, n).assign(M);
        
        DoubleMatrix2D rhs = DoubleFactory2D.dense.make(numConstr+n,1);//column vector
        rhs.viewColumn(0).viewPart(0,n).assign(q);
        rhs.viewColumn(0).viewPart(0,n).assign(Functions.mult(-1));
        rhs.viewColumn(0).viewPart(n, numConstr).assign(v);
        
        //DEBUG!! If obj fcn is not even convex on the hyperplane defined by the active constraints,
        //this method will get stuck at a saddle point, 
        //maybe even maximum if it is concave on this manifold
        
        DoubleMatrix2D ans = null;
        try{
            ans = Algebra.DEFAULT.solve(lhs, rhs);
        }
        catch(Exception e){
            System.out.println("Warning: minimizing quadratic on manifold failed: "+e.getMessage());
            return null;
        }
        
        return ans.viewColumn(0).viewPart(0, n);
    }
    
    
    //The following may be replaced at some point by analytical derivatives
    DoubleMatrix1D calcGradient(DoubleMatrix1D z){
        int n = objFcn.getNumDOFs();
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(n);
        objFcn.setDOFs(z);
        for(int dim=0; dim<n; dim++){
            double diff = objFcn.getValForDOF(dim, z.get(dim)+numDiffInterval) - 
                    objFcn.getValForDOF(dim, z.get(dim)-numDiffInterval);
            ans.set(dim, diff/(2*numDiffInterval));
            objFcn.setDOF(dim, z.get(dim));
        }
        return ans;
    }
    
    DoubleMatrix2D calcHessian(DoubleMatrix1D z){
        int n = objFcn.getNumDOFs();
        DoubleMatrix2D ans = DoubleFactory2D.dense.make(n,n);
        objFcn.setDOFs(z);
        for(int dim=0; dim<n; dim++){
            for(int dim2=0; dim2<dim; dim2++){//cross-second derivatives
                objFcn.setDOF(dim, z.get(dim)+numDiffInterval);
                double diff1 = objFcn.getValForDOF(dim2, z.get(dim2)+numDiffInterval) - 
                        objFcn.getValForDOF(dim2, z.get(dim2)-numDiffInterval);
                objFcn.setDOF(dim, z.get(dim)-numDiffInterval);
                double diff2 = objFcn.getValForDOF(dim2, z.get(dim2)+numDiffInterval) - 
                        objFcn.getValForDOF(dim2, z.get(dim2)-numDiffInterval);
                double deriv2 = (diff1-diff2)/(4*numDiffInterval*numDiffInterval);
                objFcn.setDOF(dim, z.get(dim));
                objFcn.setDOF(dim2, z.get(dim2));
                ans.set(dim, dim2, deriv2);
                ans.set(dim2, dim, deriv2);
            }
            
            //now do the same-dimension second derivative
            double diff = objFcn.getValForDOF(dim, z.get(dim)+numDiffInterval) + 
                    objFcn.getValForDOF(dim, z.get(dim)-numDiffInterval)
                    - 2*objFcn.getValForDOF(dim, z.get(dim));
            ans.set(dim, dim, diff/(numDiffInterval*numDiffInterval));
        }
        return ans;
    }
    
    DoubleMatrix1D altCalcGradient(DoubleMatrix1D z){
        int n = objFcn.getNumDOFs();
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(n);
        for(int dim=0; dim<n; dim++){
            DoubleMatrix1D zPert = z.copy();
            zPert.set(dim, z.get(dim)+10*numDiffInterval);//mess with interval a little to check stability
            double Eup = objFcn.getValue(zPert);//slower but want to check consistency
            zPert.set(dim, z.get(dim)-10*numDiffInterval);
            double Edown = objFcn.getValue(zPert);
            ans.set(dim, (Eup-Edown)/(20*numDiffInterval));
        }
        return ans;
    }
    
    DoubleMatrix2D altCalcHessian(DoubleMatrix1D z){
        int n = objFcn.getNumDOFs();
        DoubleMatrix2D ans = DoubleFactory2D.dense.make(n,n);
        for(int dim=0; dim<n; dim++){
            DoubleMatrix1D zPert = z.copy();
            zPert.set(dim, z.get(dim)+10*numDiffInterval);//mess with interval a little to check stability
            DoubleMatrix1D gradUp = altCalcGradient(zPert);//slower but want to check consistency
            zPert.set(dim, z.get(dim)-10*numDiffInterval);
            DoubleMatrix1D gradDown = altCalcGradient(zPert);
            
            DoubleMatrix1D ansRow = ans.viewRow(dim);
            ansRow.assign(gradUp);
            ansRow.assign(gradDown, Functions.minus);
            ansRow.assign(Functions.mult(0.05/numDiffInterval));
        }
        return ans;
    }
    
    
    static DoubleMatrix2D posSemidefApprox(DoubleMatrix2D M){
        //Make M positive semidefinite by zeroing out its negative eigenvalues
        EigenvalueDecomposition ed = new EigenvalueDecomposition(M);
        DoubleMatrix2D D = ed.getD();
        
        for(int dim=0; dim<M.rows(); dim++){
            if(D.get(dim,dim)<0)
                D.set(dim, dim, constrViolTol);
        }
        
        return Algebra.DEFAULT.mult(ed.getV(), Algebra.DEFAULT.mult(D, Algebra.DEFAULT.transpose(ed.getV())));
    }
    
    
    boolean compInitVals(){
        //compute good initial values, set to x
        
        x = LPChecks.getInteriorPt(constr);
        if(x==null)//constraints inconsistent
            return false;
        boxCenter = x.copy().toArray();//this should do for box center as well (used to initialize QP subproblems)
        return true;
        
        //DEBUG!!!  this is not very efficient, 
        //also assumes have all box constraints (in (LB,UB) pairs in order of dim) followed by
        //all other constraints.  
        //going to just use JOptimizer to find the spot closest to the center of the box
        //that satisfies all the constraints
        //by adding one non-box constraint at a time
        
        /*int n = objFcn.getNumDOFs();
        boxCenter = new double[n];
        for(int dim=0; dim<n; dim++)
            boxCenter[dim] = 0.5 * (constr[2*dim].getR()-constr[2*dim+1].getR());
            
        //create (quadratic) function: squared distance to center of box
        ConvexMultivariateRealFunction distSq = FunctionsUtils.createCircle(n, 0, boxCenter);
        
        int numNonBoxConstr = constr.length-2*n;
        
        double curCenter[] = boxCenter;
        for(int nbc=0; nbc<numNonBoxConstr; nbc++){
            LinearMultivariateRealFunction[] ineq = new LinearMultivariateRealFunction[nbc+1];
            System.arraycopy(constr, 2*n, ineq, 0, nbc+1);
            
            //do the optimization
            OptimizationRequest or = new OptimizationRequest();
            or.setF0(distSq);
            or.setInitialPoint(curCenter);
            or.setTolerance(qpTol);
            or.setFi(ineq);

            JOptimizer opt = new JOptimizer();
            opt.setOptimizationRequest(or);
            try{
                opt.optimize();
            }
            catch(Exception e){
                e.printStackTrace();
                throw new RuntimeException(e.getMessage());
            }

            curCenter = opt.getOptimizationResponse().getSolution();
        }
        
        x = DoubleFactory1D.dense.make(curCenter);
        return true;*/
    }
    
    
    //for serialization
    //the issue is that LinearMultivariateRealFunction isn't serializable
    private void readObject(java.io.ObjectInputStream stream)
         throws IOException, ClassNotFoundException {
         stream.defaultReadObject();
         constr = new LinearMultivariateRealFunction[numConstr];
         for(int c=0; c<numConstr; c++){
             double q[] = (double[])stream.readObject();
             double r = stream.readDouble();
             constr[c] = new LinearMultivariateRealFunction(q,r);
         }
    }
    
    private void writeObject(java.io.ObjectOutputStream stream)
         throws IOException {
         stream.defaultWriteObject();
         for(int c=0; c<numConstr; c++){
             stream.writeObject(constr[c].getQ().toArray());
             stream.writeDouble(constr[c].getR());
         }
    }
    
}
