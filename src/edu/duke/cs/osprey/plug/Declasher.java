/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import com.joptimizer.functions.LinearMultivariateRealFunction;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.epic.EPICEnergyFunction;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.minimization.SQPMinimizer;
import static edu.duke.cs.osprey.plug.LPChecks.constrAsFunction;
import java.util.ArrayList;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.NoFeasibleSolutionException;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;

/**
 *
 * Finds the point in a voxel that maximally buffers the linearized clashes
 * 
 * @author mhall44
 */
public class Declasher {
    
    public static double declashedE(SearchProblem sp, RCTuple tup){
        return declash(sp,tup).energy;
    }
    
    public static Minimizer.Result declash(SearchProblem sp, RCTuple tup){
        EPICEnergyFunction efunc = sp.epicMat.internalEnergyFunction(tup,false);
        if(efunc==null){//indicates bad conf
            return new Minimizer.Result(null, Double.POSITIVE_INFINITY);
        }

        //efunc.includeMinE = true;//might was well just difference the cont parts...
        MoleculeModifierAndScorer mms = new MoleculeModifierAndScorer(efunc,sp.epicMat.getConfSpace(),tup); 
        ArrayList<LinearConstraint> constrList = sp.plugMat.getFullStericPolytope(tup, mms.getDOFs());

        LinearMultivariateRealFunction constr[] = new LinearMultivariateRealFunction[constrList.size()];
        for(int c=0; c<constrList.size(); c++){
            constr[c] = LPChecks.toLinearMultivariateRealFunction(constrList.get(c));
        }
        
        DoubleMatrix1D pt = getDeclashedPt(constr);
        if(pt==null)
            return new Minimizer.Result(null, Double.POSITIVE_INFINITY);
        else
            return new Minimizer.Result(pt, mms.getValue(pt));
    }
    
    
    public static DoubleMatrix1D getDeclashedPt(LinearMultivariateRealFunction[] constr){
        //Identify clash-based constraints (not parallel to an axis)
        //Find the point satisfying constr that maximizes the buffer a
        //such that the clash-based constraints are satisfied within buffer a
        //return null if can't get a positive buffer (>minBuffer)
        double minBuffer = 1e-6;
        
        //Add a as last variable
        //Then constr: regular or 
        //f(x)+a<0
        //maximize a
        
        if(constr.length==0)
            return DoubleFactory1D.dense.make(0);
        
        ArrayList<LinearConstraint> polytope = new ArrayList<>();
        boolean hasClashConstr = false;
        for(LinearMultivariateRealFunction f : constr){
            double fCoeffs[] = f.getQ().toArray();
            double[] augCoeffs = new double[fCoeffs.length+1];
            System.arraycopy(fCoeffs, 0, augCoeffs, 0, fCoeffs.length);
            if(!isOneHot(f)){//add buffer
                augCoeffs[fCoeffs.length] = 1;
                hasClashConstr = true;
            }
            LinearConstraint augConstr = new LinearConstraint(augCoeffs, Relationship.LEQ, -f.getR());
            polytope.add(augConstr);
            
            
            //DEBUG!!!!!!  checking gets center if just do vox constr w/ buffer
            /*if(isOneHot(f)){
                augCoeffs[fCoeffs.length] = 1;
                hasClashConstr = true;
                LinearConstraint augConstr = new LinearConstraint(augCoeffs, Relationship.LEQ, -f.getR());
            polytope.add(augConstr);
            }*/
        }
        
        if(!hasClashConstr)//LP here will be unbounded, just use voxel center
            return LPChecks.getInteriorPt(constr);
        
        double bufferCoeffs[] = new double[polytope.get(0).getCoefficients().getDimension()];
        bufferCoeffs[bufferCoeffs.length-1] = 1;
        LinearObjectiveFunction bufferFunc = new LinearObjectiveFunction(bufferCoeffs, 0);
        
        SimplexSolver ss = new SimplexSolver();
        LinearConstraintSet polytopeConstr = new LinearConstraintSet(polytope);
        PointValuePair ansPair = ss.optimize(bufferFunc, polytopeConstr, GoalType.MAXIMIZE);
        
        double[] augAns = ansPair.getPoint();
        if(augAns[augAns.length-1]>minBuffer){//buffer OK
            double[] ans = new double[augAns.length-1];
            System.arraycopy(augAns, 0, ans, 0, ans.length);
            System.out.println("Buffer: "+augAns[augAns.length-1]);//DEBUG!!
            return DoubleFactory1D.dense.make(ans);
        }
        else//constraints inconsistent, no declashed pt possible
            return null;
    }

    
    private static boolean isOneHot(LinearMultivariateRealFunction f){
        double tol = 1e-10;
        int nonZeroCount = 0;
        for(double d : f.getQ().toArray()){
            if(Math.abs(d)>tol){
                nonZeroCount++;
                if(Math.abs(Math.abs(d)-1)>tol)
                    return false;
            }
        }
        return nonZeroCount==1;
    }
    
    
    //declash is generally several kcal/mol better worse than min, way more than regular/plug min diff
    //even for sc
    //worse than vox center
    //this probably makes it a poor choice to use instead of min.  Maybe to define well for G?
}
