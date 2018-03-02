/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.energy.derivs;

import java.util.HashMap;

/**
 *
 * Uses already-computed coordinate derivatives to get the gradient and Hessian
 * of a ForcefieldEnergy with respect to the (internal) coordinates
 * 
 * 
 * @author mhall44
 */
public class ResPairEnergyDerivs {
    
    //stuff used in computation
    double[][][] ffTermsAndDerivs;//from ForcefieldEnergy
    double dist[][];//distances between pairs of atoms
    ResCoordDerivs cd1, cd2;
    int numAtoms1, numAtoms2, numDOFs;

    //OK let's say the first res1NumDOFs dofs are just those of res 1
    //then we'll create a mapping from overall DOFs to res 2 DOFs, marking -1 if none
    int res1NumDOFs;
    int[] res2DOFIndices;
    
    
    //the energies that we will want
    double E;
    double[] grad;
    double[][] hess;
    
    
    public ResPairEnergyDerivs(double[][][] ffTermsAndDerivs, ResCoordDerivs cd1, ResCoordDerivs cd2){
        this.cd1 = cd1;
        this.cd2 = cd2;
        this.ffTermsAndDerivs = ffTermsAndDerivs;
        
        numAtoms1 = cd1.numAtoms;
        numAtoms2 = cd2.numAtoms;

        setupDOFMapping();
        
        precomputeDistances();//between atoms
        
        calcE();//maybe?
        calcGrad();
        calcHess();
    }
    

    
    private void precomputeDistances(){
        dist = new double[numAtoms1][numAtoms2];
        for(int i=0; i<numAtoms1; i++){
            for(int j=0; j<numAtoms2; j++){
                double distsq = 0;
                for(int dim=0; dim<3; dim++){
                    double diff = cd1.getCoord(i,dim) - cd2.getCoord(j,dim);
                    distsq += diff*diff;
                }
                dist[i][j] = Math.sqrt(distsq);
                
                //If i and j are the same atom this messes up some formulas, but energy should be 0
                //so set distance to 1 arbitrarily, while checking that energy=0 so it doesn't matter
                if(distsq<1e-10){
                    dist[i][j] = 1;
                    if(Math.abs(ffTermsAndDerivs[i][j][0])>1e-10)
                        throw new RuntimeException("ERROR: Unexpected self terms...");
                }
            }
        }
    }
    
    private double coordDiff(int i, int j, int k){
        //differences in coord k between atoms i and j
        return cd1.getCoord(i,k) - cd2.getCoord(j,k);
    }
    
    
    private void calcE(){
        E = 0;
        for(int i=0; i<numAtoms1; i++){
            for(int j=0; j<numAtoms2; j++){
                E += ffTermsAndDerivs[i][j][0];
            }
        }
    }
    
    private void calcGrad(){
        grad = new double[numDOFs];
        
        for(int i=0; i<numAtoms1; i++){
            double[] q = new double[3];
            for(int j=0; j<numAtoms2; j++){
                for(int k=0; k<3; k++)
                    q[k] += ffTermsAndDerivs[i][j][1] * coordDiff(i,j,k) / dist[i][j];
            }
            for(int a=0; a<numDOFs; a++){
                if(dofAffectsRes1(a)){
                    for(int k=0; k<3; k++)
                        grad[a] += q[k]*cd1.getGrad(i, res1DOFIndex(a), k);
                }
            }
        }
        
        for(int j=0; j<numAtoms2; j++){
            double[] q = new double[3];
            for(int i=0; i<numAtoms1; i++){
                for(int k=0; k<3; k++)
                    q[k] += ffTermsAndDerivs[i][j][1] * (-coordDiff(i,j,k)) / dist[i][j];
            }
            for(int a=0; a<numDOFs; a++){
                if(dofAffectsRes2(a)){
                    for(int k=0; k<3; k++)
                        grad[a] += q[k]*cd2.getGrad(j, res2DOFIndex(a), k);
                }
            }
        }
    }
    
    
    private void calcHess(){
        //OK there are a bunch of terms
        hess = new double[numDOFs][numDOFs];
        
        //term 1
        for(int i=0; i<numAtoms1; i++){
            double[] q = new double[3];
            for(int j=0; j<numAtoms2; j++){
                for(int k=0; k<3; k++)
                    q[k] += ffTermsAndDerivs[i][j][1] * coordDiff(i,j,k) / dist[i][j];
            }
            for(int a=0; a<numDOFs; a++){
                if(dofAffectsRes1(a)){
                    for(int b=0; b<numDOFs; b++){
                        if(dofAffectsRes1(b)){
                            for(int k=0; k<3; k++)
                                hess[a][b] += q[k]*cd1.getHess(i, res1DOFIndex(a), res1DOFIndex(b), k);
                        }
                    }
                }
            }
        }
       
        //Term 1 with i,j switched...
        for(int j=0; j<numAtoms2; j++){
            double[] q = new double[3];
            for(int i=0; i<numAtoms1; i++){
                for(int k=0; k<3; k++)
                    q[k] += ffTermsAndDerivs[i][j][1] * (-coordDiff(i,j,k)) / dist[i][j];
            }
            for(int a=0; a<numDOFs; a++){
                if(dofAffectsRes2(a)){
                    for(int b=0; b<numDOFs; b++){
                        if(dofAffectsRes2(b)){
                            for(int k=0; k<3; k++)
                                hess[a][b] += q[k]*cd2.getHess(j, res2DOFIndex(a), res2DOFIndex(b), k);
                        }
                    }
                }
            }
        }
        
        
        //term 2
        for(int i=0; i<numAtoms1; i++){
            double[][] q = new double[3][3];
            for(int j=0; j<numAtoms2; j++){
                double dist3 = dist[i][j]*dist[i][j]*dist[i][j];
                for(int k=0; k<3; k++){
                    for(int L=0; L<3; L++){
                        q[k][L] -= ffTermsAndDerivs[i][j][1] * coordDiff(i,j,k) * coordDiff(i,j,L) / dist3;
                    }
                    q[k][k] += ffTermsAndDerivs[i][j][1] / dist[i][j];
                }
            }
            for(int a=0; a<numDOFs; a++){
                if(dofAffectsRes1(a)){
                    double r[] = new double[3];
                    for(int k=0; k<3; k++){
                        for(int L=0; L<3; L++){
                            r[L] += cd1.getGrad(i, res1DOFIndex(a), k) * q[k][L];
                        }
                    }
                    for(int b=0; b<numDOFs; b++){
                        if(dofAffectsRes1(b)){
                            for(int L=0; L<3; L++)
                                hess[a][b] += cd1.getGrad(i, res1DOFIndex(b), L) * r[L];
                        }
                    }
                }
            }
        }
        
        
        //term 2, switched
        for(int j=0; j<numAtoms2; j++){
            double[][] q = new double[3][3];
            for(int i=0; i<numAtoms1; i++){
                double dist3 = dist[i][j]*dist[i][j]*dist[i][j];
                for(int k=0; k<3; k++){
                    for(int L=0; L<3; L++){
                        q[k][L] -= ffTermsAndDerivs[i][j][1] * coordDiff(i,j,k) * coordDiff(i,j,L) / dist3;
                    }
                    q[k][k] += ffTermsAndDerivs[i][j][1] / dist[i][j];
                }
            }
            for(int a=0; a<numDOFs; a++){
                if(dofAffectsRes2(a)){
                    double r[] = new double[3];
                    for(int k=0; k<3; k++){
                        for(int L=0; L<3; L++){
                            r[L] += cd2.getGrad(j, res2DOFIndex(a), k) * q[k][L];
                        }
                    }
                    for(int b=0; b<numDOFs; b++){
                        if(dofAffectsRes2(b)){
                            for(int L=0; L<3; L++)
                                hess[a][b] += cd2.getGrad(j, res2DOFIndex(b), L) * r[L];
                        }
                    }
                }
            }
        }
        
        
        //term 3 
        for(int j=0; j<numAtoms2; j++){
            for(int a=0; a<numDOFs; a++){
                if(dofAffectsRes1(a)){
                    double q[] = new double[3];
                    for(int i=0; i<numAtoms1; i++){
                        double dist3 = dist[i][j]*dist[i][j]*dist[i][j];
                        double r=0;
                        for(int k=0; k<3; k++)
                            r += coordDiff(i,j,k) * cd1.getGrad(i, res1DOFIndex(a), k);
                        for(int k=0; k<3; k++)
                            q[k] -= ffTermsAndDerivs[i][j][1] * (-r*coordDiff(i,j,k)/dist3 + cd1.getGrad(i,res1DOFIndex(a),k)/dist[i][j]);
                    }
                    for(int b=0; b<numDOFs; b++){
                        if(dofAffectsRes2(b)){
                            for(int k=0; k<3; k++)
                                hess[a][b] += q[k] * cd2.getGrad(j, res2DOFIndex(b), k);
                        }
                    }
                }
            }
        }
        
        
        //term 3, switched
        for(int i=0; i<numAtoms1; i++){
            for(int a=0; a<numDOFs; a++){
                if(dofAffectsRes2(a)){
                    double q[] = new double[3];
                    for(int j=0; j<numAtoms2; j++){
                        double dist3 = dist[i][j]*dist[i][j]*dist[i][j];
                        double r=0;
                        for(int k=0; k<3; k++)
                            r += (-coordDiff(i,j,k)) * cd2.getGrad(j, res2DOFIndex(a), k);
                        for(int k=0; k<3; k++)
                            q[k] -= ffTermsAndDerivs[i][j][1] * (r*coordDiff(i,j,k)/dist3 + cd2.getGrad(j,res2DOFIndex(a),k)/dist[i][j]);
                    }
                    for(int b=0; b<numDOFs; b++){
                        if(dofAffectsRes1(b)){
                            for(int k=0; k<3; k++)
                                hess[a][b] += q[k] * cd1.getGrad(i, res1DOFIndex(b), k);
                        }
                    }
                }
            }
        }
        
        
        //term 4, ii
        for(int i=0; i<numAtoms1; i++){
            double[][] q = new double[3][3];
            for(int j=0; j<numAtoms2; j++){
                double dist2 = dist[i][j]*dist[i][j];
                for(int k=0; k<3; k++){
                    for(int L=0; L<3; L++){
                        q[k][L] += ffTermsAndDerivs[i][j][2] * coordDiff(i,j,k) * coordDiff(i,j,L) / dist2;
                    }
                }
            }
            for(int a=0; a<numDOFs; a++){
                if(dofAffectsRes1(a)){
                    double r[] = new double[3];
                    for(int k=0; k<3; k++){
                        for(int L=0; L<3; L++){
                            r[L] += cd1.getGrad(i, res1DOFIndex(a), k) * q[k][L];
                        }
                    }
                    for(int b=0; b<numDOFs; b++){
                        if(dofAffectsRes1(b)){
                            for(int L=0; L<3; L++)
                                hess[a][b] += cd1.getGrad(i, res1DOFIndex(b), L) * r[L];
                        }
                    }
                }
            }
        }
        
        //term 4, ij
        for(int i=0; i<numAtoms1; i++){
            for(int b=0; b<numDOFs; b++){
                if(dofAffectsRes2(b)){
                    double[] q = new double[3];
                    for(int j=0; j<numAtoms2; j++){
                        double r = 0;
                        for(int k=0; k<3; k++)
                            r += (-coordDiff(i,j,k)) * cd2.getGrad(j, res2DOFIndex(b), k);
                        r *= ffTermsAndDerivs[i][j][2] / dist[i][j];
                        for(int k=0; k<3; k++){
                            q[k] += r * coordDiff(i,j,k) / dist[i][j];
                        }
                    }
                    for(int a=0; a<numDOFs; a++){
                        if(dofAffectsRes1(a)){
                            for(int k=0; k<3; k++){
                                hess[a][b] += cd1.getGrad(i, res1DOFIndex(a), k) * q[k];
                            }
                        }
                    }
                }
            }
        }
        
        //term 4, ji
        for(int j=0; j<numAtoms2; j++){
            for(int b=0; b<numDOFs; b++){
                if(dofAffectsRes1(b)){
                    double[] q = new double[3];
                    for(int i=0; i<numAtoms1; i++){
                        double r = 0;
                        for(int k=0; k<3; k++)
                            r += coordDiff(i,j,k) * cd1.getGrad(i, res1DOFIndex(b), k) / dist[i][j];
                        r *= ffTermsAndDerivs[i][j][2];
                        for(int k=0; k<3; k++){
                            q[k] += r * (-coordDiff(i,j,k)) / dist[i][j];
                        }
                    }
                    for(int a=0; a<numDOFs; a++){
                        if(dofAffectsRes2(a)){
                            for(int k=0; k<3; k++){
                                hess[a][b] += cd2.getGrad(j, res2DOFIndex(a), k) * q[k];
                            }
                        }
                    }
                }
            }
        }
                  
        //term 4, jj
        for(int j=0; j<numAtoms2; j++){
            double[][] q = new double[3][3];
            for(int i=0; i<numAtoms1; i++){
                double dist2 = dist[i][j]*dist[i][j];
                for(int k=0; k<3; k++){
                    for(int L=0; L<3; L++){
                        q[k][L] += ffTermsAndDerivs[i][j][2] * (-coordDiff(i,j,k)) * (-coordDiff(i,j,L)) / dist2;
                    }
                }
            }
            for(int a=0; a<numDOFs; a++){
                if(dofAffectsRes2(a)){
                    double r[] = new double[3];
                    for(int k=0; k<3; k++){
                        for(int L=0; L<3; L++){
                            r[L] += cd2.getGrad(j, res2DOFIndex(a), k) * q[k][L];
                        }
                    }
                    for(int b=0; b<numDOFs; b++){
                        if(dofAffectsRes2(b)){
                            for(int L=0; L<3; L++)
                                hess[a][b] += cd2.getGrad(j, res2DOFIndex(b), L) * r[L];
                        }
                    }
                }
            }
        }
    }
    
    
   
    
    
    private void setupDOFMapping(){
        res1NumDOFs = cd1.dofs.size();
        
        HashMap<Integer,Integer> res1Res2DOFMapping = new HashMap<>();//map indices of DOFs within res1
        //to indices within res2 (if they are shared DOFs)
        for(int res2DOFIndex=0; res2DOFIndex<cd2.dofs.size(); res2DOFIndex++){
            int res1DOFIndex = cd1.dofs.indexOf(cd2.dofs.get(res2DOFIndex));
            if(res1DOFIndex!=-1)
                res1Res2DOFMapping.put(res1DOFIndex,res2DOFIndex);
        }
        
        numDOFs = res1NumDOFs + cd2.dofs.size() - res1Res2DOFMapping.size();
        //number of DOFs is sum of each res' number of DOFs, minus number shared
        
        res2DOFIndices = new int[numDOFs];
        for(int dof=0; dof<res1NumDOFs; dof++){//first handle res2DOFIndices for res 1 DOFs
            if(res1Res2DOFMapping.containsKey(dof))
                res2DOFIndices[dof] = res1Res2DOFMapping.get(dof);
            else
                res2DOFIndices[dof] = -1;
        }
        //now handle res2-only DOFs
        int res2DOFIndex = 0;
        for(int dof=res1NumDOFs; dof<numDOFs; dof++){
            while(res1Res2DOFMapping.containsValue(res2DOFIndex))//skip shared DOFs
                res2DOFIndex++;
            res2DOFIndices[dof] = res2DOFIndex;
            res2DOFIndex++;
        }
    }
    
    //handling residue mappings
    public int mmsDOFIndex(int a){
        //map DOF index for this ResPairEnergyDerivs
        //to index in the MoleculeModifierAndScorer used to make it
        if(a<res1NumDOFs)
            return cd1.dofs.get(a);
        else
            return cd2.dofs.get(res2DOFIndex(a));
    }
    
    private boolean dofAffectsRes1(int a){
        return a<res1NumDOFs;
    }
    
    private int res1DOFIndex(int a){//map DOF index for this ResPairEnergyDerivs
        //to DOF index within first residue's DOFs
        if(a>=res1NumDOFs)
            throw new RuntimeException("ERROR: DOF doesn't affect res1");
        return a;
    }
    
    private boolean dofAffectsRes2(int b){
        return res2DOFIndices[b] != -1;
    }
    
    private int res2DOFIndex(int b){
        return res2DOFIndices[b];//will return -1 if not in res2, but that shouldn't happen
    }
    
}
