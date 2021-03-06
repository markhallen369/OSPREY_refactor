/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.aibis;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import java.io.Serializable;
import java.util.HashMap;
import java.util.List;
import static java.util.stream.Collectors.toList;

/**
 *
 * @author mhall44
 */
public class ConfSample implements Serializable {
    
    RCTuple vox;
    double[] contDOFValsArr;//values for each of the DOFs in the IBIS system (can be 0 if inapplicable)
    
    //If doing explicit bulk solvent will want a list of coords for the solvent molecules too
    
        
    //later will cache energy too
    
    public ConfSample(RCTuple vox, double[] contDOFValsArr){
        this.vox = vox;
        this.contDOFValsArr = contDOFValsArr;
    }
    /*public ConfSample(RCTuple vox, HashMap<String,Double> contDOFVals){
        if(contDOFVals!=null)
            throw new RuntimeException("ERROR DEPRECATED I THINK");
        this.vox = vox;
        this.contDOFVals = contDOFVals;
    }
    
    public ConfSample(RCTuple vox, HashMap<String,Double> contDOFVals, double[] contDOFValsArr){
        this.vox = vox;
        this.contDOFVals = contDOFVals;
        this.contDOFValsArr = contDOFValsArr;
    }*/
    
    public int getRC(int pos){
        return vox.RCAtPos(pos);
    }
    
    /*public void applyContDOFs(MoleculeModifierAndScorer mms){
        for(int dofNum=0; dofNum<mms.getNumDOFs(); dofNum++){
            //DEBUG!!! may benefit from settings DOFs at once
            //Look into blocks for MMS!!!
            DegreeOfFreedom dof = mms.getDOFs().get(dofNum);
            if(contDOFVals.containsKey(dof.getName()))//we have a value for this dof
                mms.setDOF(dofNum, contDOFVals.get(dof.getName()));
        }
    }*/
    
}
