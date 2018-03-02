/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 *
 * This is for running on clusters with a time limit
 * Precomputes matrices including LUTE; if cut off just rerun and can even resume
 * during EPIC or LUTE precomp
 * when finished, the matrices will be both loaded into sp and saved to disk
 * 
 * @author mhall44
 */
public class ResumableMatrixPrecomputer {
    
    SearchProblem sp;
    
    //parameters controlling what's precomputed
    double pruningInterval;
    boolean useEMat;
    boolean useEPIC;
    boolean usePLUG;
    boolean useTupExp;
    
    public ResumableMatrixPrecomputer(ConfigFileParser cfp){
        sp = cfp.getSearchProblem();
    }
    
    public void precomputeMatrices(){
        //the actual precomputation
        
        
        
    }
    
    
    
    static boolean load(SearchProblem sp, double pruningInterval, boolean useEMat,
            boolean useEPIC, boolean usePLUG, boolean useTupExp){
        //load precomputed matrices into a search problem  
        //in case of I-value issues we will allow it to return false
        
        //there will always be a pruning matrix
        sp.pruneMat = (PruningMatrix) ObjectIO.readObject(sp.name+".PRUNEMAT.dat", true);
        if(sp.pruneMat==null){
            System.out.println("FAILED TO LOAD RESUMABLY PRECOMPUTED MATRICES: No pruning matrix");
            return false;
        }
        
        double pruneMatInterval = sp.pruneMat.getPruningInterval();
        if(pruneMatInterval<pruningInterval-0.00001){
            System.out.println("FAILED TO LOAD RESUMABLY PRECOMPUTED MATRICES: Pruning matrix interval="
                    +pruneMatInterval+" but should be at least "+pruningInterval);
            return false;
        }
        
        if(usePLUG){//let's have an error if PLUG expected but missing or bad
            sp.plugMat = (PolytopeMatrix) ObjectIO.readObject(sp.name+".PLUGMAT.dat", false);
        }
        
        //OK load the rest as usual
        if(useEMat)
            sp.loadEnergyMatrix();
        if(useEPIC)
            sp.loadEPICMatrix();
        if(useTupExp)
            sp.loadTupExpEMatrix();
        
        return true;
    }
}
