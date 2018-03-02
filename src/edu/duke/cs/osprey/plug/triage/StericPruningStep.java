/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug.triage;

import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author mhall44
 */
public class StericPruningStep extends GeomPruningStep {

    static double stericPruningOverlap = 1;//how much overlap at central conf is considered
    //to suffice for steric pruning (i.e., w/o even checking if clash can be escaped)
    
    
    StericPruningStep(ArrayList<RCClassification> classesAtRes, PruningLookup prooningLookup, 
            RCClass.ClassType classType, ArrayList<Residue> shellResidues, CentralStericCache stericCache) {
        super(classesAtRes, prooningLookup, classType, shellResidues, stericCache);
    }
    
    StericPruningStep(ArrayList<RCClassification> classesAtRes, PruningLookup prooningLookup, 
            RCClass.ClassType classType, RCClass.ClassType classType2, ArrayList<Residue> shellResidues, CentralStericCache stericCache) {
        super(classesAtRes, prooningLookup, classType, classType2, shellResidues, stericCache);
    }
    
    
    boolean canPruneClass(RCClass cl){
        //Can prune class based on intra+shell sterics
        ResSphereModel rcModel = stericCache.getSphereModel(cl);
        if(rcModel==null)//if steric cache returns a non-null model, no intra clashes
            return true;
        HashMap<Residue,ResPairClashExemptions> exemptionMap = stericCache.getShellExemptions(cl);//what shell residue may have clash exemptions
        
        for(int shellPos=0; shellPos<shellResidues.size(); shellPos++){
            Residue res = shellResidues.get(shellPos);
            ResSphereModel sm = stericCache.getShellSphereModel(shellPos);
            ResPairClashExemptions exemptions = null;
            if(exemptionMap.containsKey(res))
                exemptions = exemptionMap.get(res);
            if(rcModel.clashesWith(sm, stericPruningOverlap, exemptions))
                return true;
        }
        
        return false;
    }
    
    
    boolean canPrunePair(RCClass cl1, RCClass cl2){
        //prune pair based on sterics
        ResSphereModel rcModel1 = stericCache.getSphereModel(cl1);
        ResSphereModel rcModel2 = stericCache.getSphereModel(cl2);
        if(rcModel1==null || rcModel2==null)
            throw new RuntimeException("ERROR: Trying to prune pair but at least one RC bad by itself");
        
        ResPairClashExemptions exemptions = stericCache.getRCPairExemptions(cl1, cl2);
        return rcModel1.clashesWith(rcModel2, stericPruningOverlap, exemptions);
    }
    
    @Override
    boolean canPrune(RCClassTuple tup) {
        //will prune if the cores of the RC classes clash, indicating clash unavoidable for all 
        //RC bindings of class(es)
        switch(tup.size()){
            case 1:
                return canPruneClass(tup.get(0));
            case 2:
                return canPrunePair(tup.get(0),tup.get(1));
            default:
                throw new RuntimeException("ERROR: Steric pruning only makes sense for singles and pairs");
        }
    }
    
    
}
