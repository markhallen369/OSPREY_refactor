/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.bbfree;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.ematrix.epic.SeriesFitter;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * This function checks if a particular series for full DOFs as a function of free DOFs
 * is accurate over a particular voxel
 * in terms of constraints being nearly satisfied
 * 
 * @author mhall44
 */
public class VoxelSeriesChecker {
    
    
    static int numSamples = 50;
    
    int numRes;//number of residues (including anchor residues) in the loop
    int numFreeDOFs, numFullDOFs;
    double[][] NCoord, CACoord, CCoord;//their main BB atom coords in the current sample
    
    double[][] fullDOFPolys;
    
    
    ArrayList<Double> targetConstraintVals;
    
    PepPlaneLinModel pepPlanes[];
    
    DoubleMatrix2D freeDOFMatrix;//rows of this are coefficients of free DOFs as lin comb of full

    DoubleMatrix1D freeDOFCenter;
    
    boolean fixAnchors;
    boolean useDistSqrt;
    
    public VoxelSeriesChecker(List<Residue> residues, int numFreeDOFs, int numFullDOFs, 
            double[][] fullDOFPolys, PepPlaneLinModel[] pepPlanes, DoubleMatrix2D freeDOFMatrix,
            DoubleMatrix1D freeDOFCenter, boolean fixAnchors, boolean useDistSqrt
    ) {
        //to be called in init conf!!
        numRes = residues.size();
        this.numFreeDOFs = numFreeDOFs;
        this.numFullDOFs = numFullDOFs;
        this.fullDOFPolys = fullDOFPolys;
        this.pepPlanes = pepPlanes;
        this.freeDOFMatrix = freeDOFMatrix;
        this.freeDOFCenter = freeDOFCenter;
        this.fixAnchors = fixAnchors;
        this.useDistSqrt = useDistSqrt;
        
        NCoord = new double[numRes][];
        CACoord = new double[numRes][];
        CCoord = new double[numRes][];
        
        for(int resNum=0; resNum<numRes; resNum++){
            NCoord[resNum] = residues.get(resNum).getCoordsByAtomName("N");
            CACoord[resNum] = residues.get(resNum).getCoordsByAtomName("CA");
            if(resNum==numRes-1)//use actual, instead of projected, C'
                CCoord[resNum] = residues.get(resNum).getCoordsByAtomName("C");
        }
        
        for(int resNum=0; resNum<numRes-1; resNum++){//make projected C' from N, CA
            CCoord[resNum]= pepPlanes[resNum].calcCCoords(CACoord[resNum], NCoord[resNum+1],
                    CACoord[resNum+1], true);
        }
        
        targetConstraintVals = calcConstraintVals();
    }
    
    
    
    
    
    double[] sampleResid(double[][] freeDOFVoxel){
        //draw a sample from the voxel and measure its constr resid and free-DOF resid
        
        //draw free DOFs for sample...
        DoubleMatrix1D sampFreeDOFs = DoubleFactory1D.dense.make(numFreeDOFs);
        
        for(int freeDOF=0; freeDOF<numFreeDOFs; freeDOF++){
            double voxWidth = freeDOFVoxel[1][freeDOF] - freeDOFVoxel[0][freeDOF];
            sampFreeDOFs.set( freeDOF, freeDOFVoxel[0][freeDOF] + Math.random()*voxWidth );
        }
        
        //OK now do full DOFs, placing them in NCoord and CACoord
        int fullDOFCount = 0;
        DoubleMatrix1D fullDOFVals = DoubleFactory1D.dense.make(numFullDOFs);
        
        for(int resNum=0; resNum<numRes; resNum++){
            
            if(resNum>0){//N not a free atom in first res
                for(int dim=0; dim<3; dim++){
                    NCoord[resNum][dim] = evalFullDOF(fullDOFCount,sampFreeDOFs);
                    fullDOFVals.set(fullDOFCount, NCoord[resNum][dim]);
                    fullDOFCount++;
                }
            }
            
            //CA's at end residues can't move if fixAnchors
            if( (!fixAnchors) || (resNum>0&&resNum<numRes-1) ){
                for(int dim=0; dim<3; dim++){
                    CACoord[resNum][dim] = evalFullDOF(fullDOFCount,sampFreeDOFs);
                    fullDOFVals.set(fullDOFCount, CACoord[resNum][dim]);
                    fullDOFCount++;
                }
            }
        }
        
        //Once N and CA in place, can calc C'.  Use plane projection, to match constr in jac
        for(int resNum=0; resNum<numRes-1; resNum++){
            CCoord[resNum]= pepPlanes[resNum].calcCCoords(CACoord[resNum], NCoord[resNum+1],
                    CACoord[resNum+1], true);
        }
        
        //OK now handle add up constraint resids!
        ArrayList<Double> sampConstraintVals = calcConstraintVals();
        List<Double> constrTarget = targetConstraintVals;
        //DEBUG!!!!!!!
        //sampConstraintVals = calcPepPlaneConstraintVals(2);
        //constrTarget = targetConstraintVals.subList(8, 11);
        
        int numConstr = sampConstraintVals.size();
        
        double constrResid = 0;
        for(int c=0; c<numConstr; c++){
            double dev = sampConstraintVals.get(c) - constrTarget.get(c);
            constrResid += dev*dev;
        }
        
        constrResid /= numConstr;//normalize resid by # of constraints
        
        
        DoubleMatrix1D freeDOFsCheck = Algebra.DEFAULT.mult(freeDOFMatrix, fullDOFVals);
        freeDOFsCheck.assign(freeDOFCenter, Functions.minus);
        freeDOFsCheck.assign(sampFreeDOFs, Functions.minus);//calc deviation in free DOFs
        double freeDOFResid = freeDOFsCheck.zDotProduct(freeDOFsCheck) / numFreeDOFs;
        
        return new double[] {constrResid, freeDOFResid};
    }
    
    
    
    private ArrayList<Double> calcConstraintVals(){
        //Calculate the values of the constraint variables, based on the current NCoord,
        //CACoord, and CCoord
        ArrayList<Double> ans = new ArrayList<>();
        
        for(int resNum=0; resNum<numRes; resNum++){
            if(resNum>0){//peptide plane distance constrs
                if(useDistSqrt){
                    ans.add( VectorAlgebra.distance(NCoord[resNum], CACoord[resNum]) );
                    ans.add( VectorAlgebra.distance(CACoord[resNum-1], CACoord[resNum]) );
                    ans.add( VectorAlgebra.distance(NCoord[resNum], CACoord[resNum-1]) );
                }
                else {
                    ans.add( VectorAlgebra.distsq(NCoord[resNum], CACoord[resNum]) );
                    ans.add( VectorAlgebra.distsq(CACoord[resNum-1], CACoord[resNum]) );
                    ans.add( VectorAlgebra.distsq(NCoord[resNum], CACoord[resNum-1]) );
                }
            }
            
            if(fixAnchors||(resNum>0&&resNum<numRes-1)){
                if(useDistSqrt && resNum>0 && resNum<numRes-1){
                    ans.add( VectorAlgebra.distance(NCoord[resNum], CCoord[resNum]) );
                }
                else {
                    ans.add( VectorAlgebra.dot( VectorAlgebra.subtract(NCoord[resNum], CACoord[resNum]),
                            VectorAlgebra.subtract(CCoord[resNum], CACoord[resNum]) ) );
                }
            }
        }
        
        return ans;
    }
    
    
    /*private ArrayList<Double> calcPepPlaneConstraintVals(int pepPlaneNum){
        ArrayList<Double> ans = new ArrayList<>();
        
        int resNum = pepPlaneNum+1;
        ans.add( VectorAlgebra.distance(NCoord[resNum], CACoord[resNum]) );
        ans.add( VectorAlgebra.distance(CACoord[resNum-1], CACoord[resNum]) );
        ans.add( VectorAlgebra.distance(NCoord[resNum], CACoord[resNum-1]) );
       
        return ans;
    }*/
    
    
    private double evalFullDOF(int fullDOF, DoubleMatrix1D sampFreeDOFs){
        return SeriesFitter.evalSeries(fullDOFPolys[fullDOF], sampFreeDOFs, 
                    numFreeDOFs, true, BBFreeBlock.polyOrder);
    }
    
    public double getConstraintsResid(double[][] freeDOFVoxel){
        //Get the residual of the constraint variables within the voxel
        //Also print the free-DOF resid so we can make sure it's not too big
        double constrResid = 0;
        double freeDOFResid = 0;
        
        for(int samp=0; samp<numSamples; samp++){
            double resids[] = sampleResid(freeDOFVoxel);
            constrResid += resids[0];
            freeDOFResid += resids[1];
        }
        
        constrResid /= numSamples;
        freeDOFResid /= numSamples;
        System.out.println("Free DOF resid: "+freeDOFResid);
        
        return constrResid;
    }
    
}
