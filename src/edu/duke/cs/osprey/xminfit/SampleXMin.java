/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.xminfit;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.derivs.EnergyDerivs;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import java.util.ArrayList;

/**
 *
 * Tracking the conformational minimization (within an xminfit)
 * for a single voxel sample
 * 
 * @author mhall44
 */
public class SampleXMin {
    
    SearchProblem searchSpace;//sample defined within this search space
    int[] conf;
    
    
    public SampleXMin(SearchProblem searchSpace, int[] conf) {
        this.searchSpace = searchSpace;
        this.conf = conf;
    }
    
    double[] compareXMinToMatrix(XMFMatrix XMinMatrix){
        //do a full minimization, get (this energy) - (energy from xminmatrix),
        //report {this energy difference, distance in CATS space between minima}
        MoleculeModifierAndScorer mmsFrozenCATS = makeObjFcn(XMinMatrix);
        MoleculeModifierAndScorer mmsFull = makeObjFcn(null);
        Minimizer.Result frozenCATSMin = minDOFVals(mmsFrozenCATS);
        Minimizer.Result fullMin = minDOFVals(mmsFull);
        DoubleMatrix1D xCATSFrozen = getCATSDOFs(frozenCATSMin.dofValues, mmsFrozenCATS);
        DoubleMatrix1D xCATSFull = getCATSDOFs(fullMin.dofValues, mmsFull);
        double catsDist = dist(xCATSFrozen,xCATSFull);
        return new double[] {fullMin.energy-frozenCATSMin.energy, catsDist};
    }
    
    double dist(DoubleMatrix1D a, DoubleMatrix1D b){
        DoubleMatrix1D diff = a.copy().assign(b,Functions.minus);
        return Math.sqrt(diff.zDotProduct(diff));
    }
    
    static Minimizer.Result minDOFVals(MoleculeModifierAndScorer mms){
        return new CCDMinimizer(mms,false).minimize();
    }
    
    
    double curEnergy(XMFMatrix XMinMatrix) {
        MoleculeModifierAndScorer mms = makeObjFcn(XMinMatrix);
        return minDOFVals(mms).energy;
    }
    
    DoubleMatrix1D curCATSGrad(XMFMatrix XMinMatrix){
        MoleculeModifierAndScorer mms = makeObjFcn(XMinMatrix);
        DoubleMatrix1D fullGrad = new EnergyDerivs(mms, minDOFVals(mms).dofValues).getGrad();
        //DEBUG!!  not accounting for motion in sc min.  Purely 1st order
        return getCATSDOFs(fullGrad,mms);
    }
    
    MoleculeModifierAndScorer makeObjFcn(XMFMatrix XMinMatrix){
        //if XMinMatrix not null, freeze CATS vals to the values it implies
        RCTuple confTup = new RCTuple(conf);
        EnergyFunction ef = searchSpace.epicMat.internalEnergyFunction(confTup,true);
        MoleculeModifierAndScorer ans = new MoleculeModifierAndScorer(ef, searchSpace.epicMat.getConfSpace(), confTup);
        if(XMinMatrix!=null){
            DoubleMatrix1D catsDOFVals = XMinMatrix.getXMin(conf);
            DoubleMatrix1D[] constr = ans.getConstraints();
            for(int a=0; a<2; a++)
                setCATSDOFs(constr[a], catsDOFVals, ans);
        }
        return ans;
    }
    
    private void setCATSDOFs(DoubleMatrix1D x, DoubleMatrix1D xCats, MoleculeModifierAndScorer mms){
        //given DOF definitions in mms, set the CATS dofs in x to xCats
        ArrayList<Integer> catsIndices = catsDOFIndices(mms);
        if(catsIndices.size()!=xCats.size())
            throw new RuntimeException("ERROR: wrong number of CATS DOFs");
        for(int c=0; c<catsIndices.size(); c++)
            x.set(catsIndices.get(c), xCats.get(c));
    }
    
    private DoubleMatrix1D getCATSDOFs(DoubleMatrix1D x, MoleculeModifierAndScorer mms){
        ArrayList<Integer> catsIndices = catsDOFIndices(mms);
        DoubleMatrix1D xCats = DoubleFactory1D.dense.make(catsIndices.size());
        for(int c=0; c<catsIndices.size(); c++)
            xCats.set(c, x.get(catsIndices.get(c)));
        return xCats;
    }
    
    static ArrayList<Integer> catsDOFIndices(MoleculeModifierAndScorer mms) {
        //which dofs are cats?
        ArrayList<Integer> ans = new ArrayList<>();
        ArrayList<DegreeOfFreedom> DOFs = mms.getDOFs();
        for(int dof=0; dof<DOFs.size(); dof++){
            if(DOFs.get(dof) instanceof BBFreeDOF)
                ans.add(dof);
        }
        return ans;
    }
    
    
    
    static double curESum(ArrayList<SampleXMin> samples, XMFMatrix XMinMatrix){
        double ESum = 0;
        for(SampleXMin samp : samples){
            ESum += samp.curEnergy(XMinMatrix);
        }
        return ESum;
    }
}
