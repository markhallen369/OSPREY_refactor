/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.bbfree.BBFreeBlock;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.DihedralRotation;
import edu.duke.cs.osprey.dof.deeper.RamachandranChecker;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.RigidBodyMotion;
import edu.duke.cs.osprey.tools.RotationMatrix;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.util.Arrays;
import java.util.HashSet;

/**
 *
 * Trying out big CATS voxels, hopefully big enough to cover all motions
 * accessible continuously from central conf w/o Ramachandran violation
 * 
 * 
 * @author mhall44
 */
public class LionPark {
    
    //to figure out: 1. Can paths across the voxel keep grad f nonsingular (based on central conf at
    //middle)
    //If not then need to consider splitting
    //2. If so, construct some "partial conf" voxels and see how Taylor series quality is
    //If not then need to create some kind of system for quick eval over multiple "Taylor voxels"
    
    
    String pdbFileName;
    Molecule origMolec;
    
    DoubleMatrix2D origMz;
    double center[];
    
    double maxMotion = 3;
        
    public LionPark(String pdbFileName){
        this.pdbFileName = pdbFileName;
        origMolec = PDBFileReader.readPDBFile(pdbFileName);
        int numRes = origMolec.residues.size();
        center = origMolec.residues.get(numRes/2).getCoordsByAtomName("CA");
        origMz = BBFreeBlock.getOrthogVectors( new BBFreeBlock(origMolec.residues,false).getConstrJac() );
    }
    
    public static void main(String[] args){
        ConfigFileParser cfp = new ConfigFileParser(args);//args 1, 3+ are configuration files
	cfp.loadData();
        
        //new LionPark("PEP.pdb").exploreFDet();
        Molecule molec = PDBFileReader.readPDBFile("PEP.pdb");
        double voxLim[] = new double[2*molec.residues.size()+4];
        Arrays.fill(voxLim, 1);
        /*for(int a=0; a<6; a++){
            voxLim[a]/=4;
        }*/
        BBFreeBlock bfb = new BBFreeBlock(molec.residues, false, voxLim, true);

        DoubleMatrix1D x = DoubleFactory1D.dense.random(bfb.getDOFs().size()).assign(Functions.mult(1));
        if(!bfb.setDOFs(x))
            throw new RuntimeException("ERROR: Failed to set DOFs");
        PDBFileWriter.writePDBFile(molec, "LION.pdb");
    }
    
    
    DoubleMatrix1D randomStep(){
        DoubleMatrix1D ans = DoubleFactory1D.dense.random(2*origMolec.residues.size()+4);
        ans.viewPart(0,6).assign(Functions.mult(0.03));
        return ans;
    }
    
    DoubleMatrix1D oneDirStep(int DOFNum){
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(2*origMolec.residues.size()+4);
        if(DOFNum<6)
            ans.set(DOFNum, 0.03);
        else
            ans.set(DOFNum, 1);
        return ans;
    }
            
            
    void exploreFDet(){
        
        //So in design I think I want to choose where the backbone will go within ~2 A
        //and then there are Rama constraints within that specification
        //both are roughly convex
        //so let's randomly pick a direction in trans/rot/dih space
        //and graph det grad f walking out til we hit Rama or 2 A trouble
        int numPaths = 2*origMolec.residues.size()+4;
        int numRes = origMolec.residues.size();
        int middleRes = numRes/2;
                
        
        for(int path=0; path<numPaths; path++){
            System.out.println("PATH "+path+": ");
            Molecule molec = PDBFileReader.readPDBFile(pdbFileName);
            
            DoubleMatrix1D step = randomStep();
            //DoubleMatrix1D step = oneDirStep(path+3);//will step this way til either Rama or 
            //distance gone too far
            
            
            while(true){
                
                boolean ramaBad = false;
                boolean movedTooMuch = false;
                
                for(int res=0; res<numRes; res++){
                    Residue curRes = molec.residues.get(res);
                    if(res>0 && res<numRes-1){
                        if(!RamachandranChecker.getInstance().checkByAAType(curRes)[2])
                            ramaBad = true;
                    }
                    
                    Residue origRes = origMolec.residues.get(res);
                    for(String atName : new String[] {"N","CA"}){
                        int atIndex = curRes.getAtomIndexByName(atName);
                        double dist = VectorAlgebra.distance(curRes.coords,atIndex,origRes.coords,atIndex);
                        if(dist>maxMotion)
                            movedTooMuch = true;
                    }
                }
                
                
                if(ramaBad || movedTooMuch){
                    System.out.println("HIT END OF DESIRED AREA!  RAMABAD: "+ramaBad+" MOVEDTOOM"
                            + "UCH: "+movedTooMuch);
                    break;
                }
                else {
                    DoubleMatrix2D constrJac = new BBFreeBlock(molec.residues,false).getConstrJac();
                    DoubleMatrix2D jac = DoubleFactory2D.dense.compose(
                        new DoubleMatrix2D[][] { new DoubleMatrix2D[] {Algebra.DEFAULT.transpose(origMz)},
                        new DoubleMatrix2D[] {constrJac}
                    } );
                    System.out.println(Algebra.DEFAULT.det(jac));
                }
                
                
                //apply DOF step
                applyRigidBodyMotion(molec, step.viewPart(0,6));
                int dofCount = 6;
                for(int res=middleRes-1; res>=0; res--){
                    applyPhiBackward(molec, res+1, step.get(dofCount++));
                    applyPsiBackward(molec, res, step.get(dofCount++));
                }
                for(int res=middleRes+1; res<numRes; res++){
                    applyPsiForward(molec, res-1, step.get(dofCount++));
                    applyPhiForward(molec, res, step.get(dofCount++));
                }
            }
            
            System.out.println("PATH DONE");
            PDBFileWriter.writePDBFile(molec, "PATH"+path+".pdb");
        }
    }
    
    
    void applyRigidBodyMotion(Molecule molec, DoubleMatrix1D params){
        //apply a rigid body motion to the whole molecule
        //parameters 3-5 are for rotation
        RotationMatrix mate = new RotationMatrix(1,0,0,params.get(3),true);
        mate = mate.multiply(new RotationMatrix(0,1,0,params.get(4),true));
        mate = mate.multiply(new RotationMatrix(0,0,1,params.get(5),true));
        
        double cen2[] = center.clone();
        for(int dim=0; dim<3; dim++)
            cen2[dim] += params.get(dim);
        
        RigidBodyMotion mosh = new RigidBodyMotion(center,mate,cen2);
        for(Residue res : molec.residues)
            mosh.transform(res.coords);
    }
    
    void applyPhiBackward(Molecule molec, int res, double angle){
        //apply a phi backbone dihedral, propagating backward down the chain.  Angle in degrees.  
        Residue curRes = molec.residues.get(res);
        DihedralRotation dih = new DihedralRotation(curRes.getCoordsByAtomName("CA"), curRes.getCoordsByAtomName("N"), angle);
        dih.transform(curRes.coords, curRes.getAtomIndexByName("H"));
        for(int res2=0; res2<res; res2++)
            dih.transform(molec.residues.get(res2).coords);
    }
    
    void applyPhiForward(Molecule molec, int res, double angle){
        Residue curRes = molec.residues.get(res);
        DihedralRotation dih = new DihedralRotation(curRes.getCoordsByAtomName("N"), curRes.getCoordsByAtomName("CA"), angle);
        HashSet fixedAtoms = new HashSet(Arrays.asList("N", "H", "CA"));
        for(int a=0; a<curRes.atoms.size(); a++){
            if( ! fixedAtoms.contains(curRes.atoms.get(a).name) )
                dih.transform(curRes.coords,a);
        }
        for(int res2=res+1; res2<molec.residues.size(); res2++)
            dih.transform(molec.residues.get(res2).coords);
    }
    
    void applyPsiBackward(Molecule molec, int res, double angle){
        Residue curRes = molec.residues.get(res);
        DihedralRotation dih = new DihedralRotation(curRes.getCoordsByAtomName("C"), curRes.getCoordsByAtomName("CA"), angle);
        HashSet fixedAtoms = new HashSet(Arrays.asList("C", "O", "CA"));
        for(int a=0; a<curRes.atoms.size(); a++){
            if( ! fixedAtoms.contains(curRes.atoms.get(a).name) )
                dih.transform(curRes.coords,a);
        }
        for(int res2=0; res2<res; res2++)
            dih.transform(molec.residues.get(res2).coords);
    }
    
    void applyPsiForward(Molecule molec, int res, double angle){
        Residue curRes = molec.residues.get(res);
        DihedralRotation dih = new DihedralRotation(curRes.getCoordsByAtomName("CA"), curRes.getCoordsByAtomName("C"), angle);
        dih.transform(curRes.coords, curRes.getAtomIndexByName("O"));
        for(int res2=res+1; res2<molec.residues.size(); res2++)
            dih.transform(molec.residues.get(res2).coords);
    }
    
    
    
}
