/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.energy.derivs;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.util.ArrayList;

/**
 *
 * Coordinates of a residue, and derivatives thereof wrt internal coordinates
 * 
 * @author mhall44
 */
public class ResCoordDerivs {
    
    //these are flattened multiple-dimensional arrays
    //ordering of dimensions based roughly on order of access in ResPairEnergyDerivs, for cache efficiency
    double x[];//coordinates (indices: atom, 3D dimension)
    double grad[];//gradients (indices: atoms, DOF, 3D dimension)
    double hess[];//Hessians (indices: atoms, first DOF, second DOF, 3D dimension)
    
    ArrayList<Integer> dofs;//DOFs (as indexed in the MoleculeModifierAndScorer we're working with)
    int numDOFs;
    int numAtoms;
    
    static final double numDiffInterval = 0.001;//going to rely on numerical derivatives for now
    
    
    
    public ResCoordDerivs(Residue res, ArrayList<Integer> dofs, MoleculeModifierAndScorer mms){
        this.dofs = dofs;
        numDOFs = dofs.size();
        numAtoms = res.atoms.size();
        
        x = res.coords.clone();
        calcGrad(mms,res);
        calcHess(mms,res);
    }
    
    public ResCoordDerivs(Residue res){
        //constructor for an immobile residue
        dofs = new ArrayList<>();
        numDOFs = 0;
        numAtoms = res.atoms.size();
        x = res.coords.clone();
        grad = new double[0];
        hess = new double[0];
    }
    
    private void calcGrad(MoleculeModifierAndScorer mms, Residue res){
        double derivs[][] = calcFirstDerivs(mms, res, 0);
        grad = new double[numAtoms*numDOFs*3];
        for(int at=0; at<numAtoms; at++){
            for(int dof=0; dof<numDOFs; dof++){
                for(int dim=0; dim<3; dim++){
                    setGrad(at, dof, dim, derivs[dof][3*at+dim]);
                }
            }
        }
    }
    
    private double[][] calcFirstDerivs(MoleculeModifierAndScorer mms, Residue res, int firstDOF){
        //calculate derivatives with respect to degrees of freedom firstDOF through numDOFs
        //answer[i] is derivative of res coords wrt DOF i
        double ans[][] = new double[numDOFs][numAtoms*3];
        
        for(int resDOFIndex=firstDOF; resDOFIndex<numDOFs; resDOFIndex++){
            int overallDOFIndex = dofs.get(resDOFIndex);
            double curDOFVal = mms.getCurValueOfDOF(overallDOFIndex);
            mms.setDOF(overallDOFIndex, curDOFVal+numDiffInterval);
            double coordsUp[] = res.coords.clone();
            mms.setDOF(overallDOFIndex, curDOFVal-numDiffInterval);
            double coordsDown[] = res.coords.clone();
            mms.setDOF(overallDOFIndex, curDOFVal);
            for(int k=0; k<numAtoms*3; k++)
                ans[resDOFIndex][k] = (coordsUp[k]-coordsDown[k])/(2*numDiffInterval);
        }
        
        return ans;
    }
    
    private void calcHess(MoleculeModifierAndScorer mms, Residue res){
        hess = new double[numAtoms*numDOFs*numDOFs*3];
        
        for(int resDOFIndex=0; resDOFIndex<numDOFs; resDOFIndex++){
            //in general, earlier DOFs are more likely to be expensive to set, hence this scheme
            int overallDOFIndex = dofs.get(resDOFIndex);
            double curDOFVal = mms.getCurValueOfDOF(overallDOFIndex);
            mms.setDOF(overallDOFIndex, curDOFVal+numDiffInterval);
            double coordsUp[] = res.coords.clone();//coordinates with only dof raised, not the others
            double derivsUp[][] = calcFirstDerivs(mms, res, resDOFIndex+1);
            
            mms.setDOF(overallDOFIndex, curDOFVal-numDiffInterval);
            double coordsDown[] = res.coords.clone();
            double derivsDown[][] = calcFirstDerivs(mms, res, resDOFIndex+1);
            
            mms.setDOF(overallDOFIndex, curDOFVal);
            
            //handle same-DOF second derivs
            for(int at=0; at<numAtoms; at++){
                for(int dim=0; dim<3; dim++){
                    double deriv2 = (coordsUp[3*at+dim]-2*getCoord(at,dim)+coordsDown[3*at+dim])/(numDiffInterval*numDiffInterval);
                    setHess(at, resDOFIndex, resDOFIndex, dim, deriv2);
                }
                for(int resDOFIndex2=resDOFIndex+1; resDOFIndex2<numDOFs; resDOFIndex2++){
                    for(int dim=0; dim<3; dim++){
                        double deriv2 = (derivsUp[resDOFIndex2][3*at+dim]-derivsDown[resDOFIndex2][3*at+dim])/(2*numDiffInterval);
                        setHess(at, resDOFIndex, resDOFIndex2, dim, deriv2);
                        setHess(at, resDOFIndex2, resDOFIndex, dim, deriv2);
                    }
                }
            }
        }
    }
        
    
    public double getCoord(int atomNum, int dim){
        return x[3*atomNum+dim];
    }
    
    public double getGrad(int atomNum, int dof, int dim){
        return grad[dim+3*(dof+numDOFs*atomNum)];
    }
    
    public double getHess(int atomNum, int dof, int dof2, int dim){
        return hess[dim+3*(dof2+numDOFs*(dof+numDOFs*atomNum))];
    }
    
    
    private void setGrad(int atomNum, int dof, int dim, double val){
        grad[dim+3*(dof+numDOFs*atomNum)] = val;
    }
    
    private void setHess(int atomNum, int dof, int dof2, int dim, double val){
        hess[dim+3*(dof2+numDOFs*(dof+numDOFs*atomNum))] = val;
    }
        
    
    
}
