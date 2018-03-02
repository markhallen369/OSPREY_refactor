/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import java.util.ArrayList;
import org.apache.commons.math3.optim.linear.LinearConstraint;

/**
 *
 * This is kind of like a ConfSearch, but (for now?) doing it separate from findGMEC
 * It requires no energy modeling...hopefully fast, could be good for prescreen
 * or even library design
 * 
 * @author mhall44
 */
public class PureGeomDFS {
    
    PolytopeMatrix mat;
    int numRCsAtPos[];
    boolean printStruct = true;//make a PDB file from the returned struckchur.  
    
    //OK to start we'll just go in regular order
    //but for big systems better ordering could be key
    //both in terms of what variables to try (least options first? most? dynamic?),
    //and in what order to try them (once RC rejected try further ones from it?)
    
    public PureGeomDFS(PolytopeMatrix mat){
        this.mat = mat;
        numRCsAtPos = mat.cSpace.getNumRCsAtPos();
    }
    
    public RCTuple findConf(){
        //find a good conf
        RCTuple curTup = new RCTuple();
        RCTuple conf = doDFS(curTup,0);
        
        if(printStruct && conf!=null){
            for(int pos=0; pos<mat.cSpace.numPos; pos++){//conf will be all pos in order
                String AAType = mat.cSpace.posFlex.get(pos).RCs.get(conf.RCs.get(pos)).AAType;
                mat.cSpace.mutDOFs.get(pos).mutateTo(AAType);
            }
            ArrayList<DegreeOfFreedom> DOFs = new ArrayList(mat.calcDOFBounds(conf).keySet());
            ArrayList<LinearConstraint> polytope = mat.getFullPolytope(conf);
            double[] feasPt = LPChecks.getFeasiblePt(polytope);
            for(int dof=0; dof<DOFs.size(); dof++){
                DOFs.get(dof).apply(feasPt[dof]);
            }
            PDBFileWriter.writePDBFile(mat.cSpace.m, "PURE_GEOM_DFS.pdb", 0);
        }
            
        return conf;
    }
    
    RCTuple doDFS(RCTuple startTuple, int nextPos){
        if(nextPos==mat.cSpace.numPos)//all pos done (ASSUMING FORWARD ORDERING)
            return startTuple;
        
        for(int rc=0; rc<numRCsAtPos[nextPos]; rc++){
            RCTuple biggerTup = startTuple.addRC(nextPos, rc);
            if(isTupleOK(biggerTup)){
                RCTuple fullTup = doDFS(biggerTup,nextPos+1);
                if(fullTup!=null)
                    return fullTup;
            }
        }
        
        return null;//no way to get a legit full conf containing startTuple
    }
    
    
    boolean isTupleOK(RCTuple tup){
        return mat.isTupleFeasible(tup);
        //DEBUG!!  Also check uc feasibility here (probs based on
        //if the required atoms can have uc contacts w/i vox, rather than
        //trying to figure out which contacts they are (combinatorially hard) & impose lin constr)
    }

    
}
