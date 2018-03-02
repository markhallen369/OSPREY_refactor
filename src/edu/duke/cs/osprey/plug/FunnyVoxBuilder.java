/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug;

import cern.colt.matrix.DoubleMatrix2D;
import edu.duke.cs.osprey.bbfree.BBFreeBlock;
import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.Relationship;

/**
 *
 * build funny voxels for the free backbone block
 * 
 * 
 * @author mhall44
 */
public class FunnyVoxBuilder {
    
    SearchProblem sp;
    
    int splitNum;//how many pieces to cut each dimension of voxel into
    //I think about 4 could be reasonable.  Note the number of BB vox created per res
    //is splitNum^6, so be careful: e.g. 4^6 = 4096

    double resPoseVoxWidth;//maybe about an angstrom?  As in CATS paper
    
    
    
    //i would like to cover a grid of 2 angstroms in each direction for each N and CA coord.
    //Albeit, obviously most pairs for consecutive res will not be tractable.  
    
    public FunnyVoxBuilder(SearchProblem sp, int splitNum, double resPoseVoxWidth) {
        this.sp = sp;
        this.splitNum = splitNum;
        this.resPoseVoxWidth = resPoseVoxWidth;
    }
    
    
    
    public void buildFunnyVox(){
        BBFreeBlock bfb = getBBFreeBlock();
        ArrayList<BBFreeDOF> bbVoxDOFs = bfb.getDOFs();
        
        for(int resNumInBlock=0; resNumInBlock<bfb.getResidues().size(); resNumInBlock++){
            DoubleMatrix2D resDOFCombos = bfb.selectResFreeDOFVectors(resNumInBlock);//6 linear combinations of the free DOFs
            //approximately representing the pose of this res
            Residue res = bfb.getResidues().get(resNumInBlock);
            
            ArrayList<ArrayList<LinearConstraint>> BBVoxList = generateBBVoxels(resDOFCombos);//list of lin constr
            
            int pos = -1;
            for(int p=0; p<sp.confSpace.numPos; p++){
                if(sp.confSpace.posFlex.get(p).res == res){
                    pos = p;
                    break;
                }
            }
            if(pos==-1)
                throw new RuntimeException("ERROR res in BB block not in confspace!!");
            
            ArrayList<RC> RCs = sp.confSpace.posFlex.get(pos).RCs;
            
            int numSidechainRCs = RCs.size();
            int newRCCounter = 0;
            for(int rc=0; rc<numSidechainRCs; rc++){
                RC sidechainRC = RCs.remove(0);
                for(int bbVoxNum=0; bbVoxNum<BBVoxList.size(); bbVoxNum++){
                    //make newRC = sidechainRC \cup BBVox
                    ArrayList<LinearConstraint> BBVox = BBVoxList.get(bbVoxNum);
                    RC newRC = new RC(sidechainRC);//should currently not have any lin constr
                    newRC.linConstr.addAll( constraintsForNewDOFs(BBVox,newRC.DOFs,bbVoxDOFs) );
                    newRC.bbVoxNum = bbVoxNum;
                    newRC.RCIndex = newRCCounter;
                    newRCCounter++;
                    RCs.add(newRC);
                }
            }
        }
    }

    
    
    ArrayList<ArrayList<LinearConstraint>> generateBBVoxels(DoubleMatrix2D resDOFCombos){
        //voxel up the combos
        ArrayList<ArrayList<LinearConstraint>> baseVox = new ArrayList<>();
        baseVox.add(new ArrayList<>());
        return generateBBVoxelsHelper(baseVox, resDOFCombos, 0);
    }
    
    ArrayList<ArrayList<LinearConstraint>> generateBBVoxelsHelper(ArrayList<ArrayList<LinearConstraint>> curVoxels,
            DoubleMatrix2D resDOFCombos, int nextDOF){
        //add gridding-up of nextDOF (column of resDOFCombos) to the voxels in curVoxels
        ArrayList<ArrayList<LinearConstraint>> ans = new ArrayList<>();
        double dofStart = -0.5*resPoseVoxWidth*splitNum;
        double coeffs[] = resDOFCombos.viewColumn(nextDOF).toArray();
        
        for(ArrayList<LinearConstraint> oldVox : curVoxels){
            for(int gridpt=0; gridpt<splitNum; gridpt++){
                ArrayList<LinearConstraint> newVox = new ArrayList<>();
                newVox.addAll(oldVox);
                double lbDOFVal = dofStart+gridpt*resPoseVoxWidth;
                newVox.add( new LinearConstraint(coeffs,Relationship.GEQ,lbDOFVal) );
                newVox.add( new LinearConstraint(coeffs,Relationship.LEQ,lbDOFVal+resPoseVoxWidth) );
                ans.add(newVox);
            }
        }
        
        if(nextDOF+1==resDOFCombos.columns())
            return ans;
        else//need to deal with more DOFs
            return generateBBVoxelsHelper(ans, resDOFCombos, nextDOF+1);
    }
    
    
    
    static ArrayList<LinearConstraint> constraintsForNewDOFs(ArrayList<LinearConstraint> oldConstr,
            ArrayList<DegreeOfFreedom> newDOFs, ArrayList<? extends DegreeOfFreedom> oldDOFs){
        //given a list of constraints on an "old" set of DOFs,
        //return these same constraints acting on a "new" set (must be superset of the old set)
        
        //figure out the mapping between the DOF lists
        int[] DOFIndexMapping = new int[oldDOFs.size()];
        Arrays.fill(DOFIndexMapping,-1);
        for(int oldDOF=0; oldDOF<oldDOFs.size(); oldDOF++){
            for(int newDOF=0; newDOF<newDOFs.size(); newDOF++){
                if(oldDOFs.get(oldDOF)==newDOFs.get(newDOF)){
                    DOFIndexMapping[oldDOF] = newDOF;
                }
            }
        }
        
        ArrayList<LinearConstraint> newConstrs = new ArrayList<>();
        
        for(LinearConstraint constr : oldConstr){
            double[] newCoeffs = new double[newDOFs.size()];
            for(int oldDOF=0; oldDOF<oldDOFs.size(); oldDOF++)
                newCoeffs[DOFIndexMapping[oldDOF]] = constr.getCoefficients().getEntry(oldDOF);
            LinearConstraint newConstr = new LinearConstraint(newCoeffs,constr.getRelationship(),constr.getValue());
            newConstrs.add(newConstr);
        }
        
        return newConstrs;
    }
    
    BBFreeBlock getBBFreeBlock(){
        return getBBFreeBlock(sp.confSpace);
    }
    
    public static BBFreeBlock getBBFreeBlock(ConfSpace confSpace){
        BBFreeBlock ans = null;
        for(DegreeOfFreedom dof : confSpace.confDOFs){
            if(dof instanceof BBFreeDOF){
                BBFreeBlock block = (BBFreeBlock)dof.getBlock();
                if(block!=ans && ans!=null){
                    throw new RuntimeException("ERROR: Can't currently build funny vox if there's >1 free-BB block");
                }
                ans = block;
            }
        }
        if(ans==null)
            throw new RuntimeException("ERROR: Can't build funny vox, no free-BB block");
        return ans;
    }
    
    
    
    
    
}
