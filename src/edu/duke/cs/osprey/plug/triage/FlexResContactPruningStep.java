/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug.triage;

import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * @author mhall44
 */
public class FlexResContactPruningStep extends GeomPruningStep {
    
    
    static double contactBuffer = 1.5;//Mutable residue RCs are required to be
    //within this distance of non-mutable (i.e. target) residues, measured between VDW spheres, to pass this
    //pruning step.  
    //If at this distance they could probably get into contact by within-voxel minimization

    HashSet<Residue> mutableRes = null;
    HashSet<Integer> nonBBFlexPos = new HashSet<>();
    
    FlexResContactPruningStep(ArrayList<RCClassification> classesAtRes, PruningLookup prooningLookup, 
            RCClass.ClassType classType, ArrayList<Residue> shellResidues,
            HashSet<Residue> mutableRes, CentralStericCache stericCache) {
        super(classesAtRes, prooningLookup, classType, shellResidues, stericCache);
        this.mutableRes = mutableRes;
        findNonBBFlexPos();
    }
    
    private void findNonBBFlexPos(){
        for(int pos=0; pos<classesAtRes.size(); pos++){
            boolean posNonBB = true;
            for( RCClass cl : classesAtRes.get(pos).singletonClassLookup.values() ){
                if(cl.bbVoxNum!=-1){//bb flex detected
                    posNonBB = false;
                    break;
                }
            }
            if(posNonBB)
                nonBBFlexPos.add(pos);
        }
    }
    
    @Override
    boolean canPrune(RCClassTuple tup) {
        //will prune if the first residue is mutable
        //but can't touch any non-mutable residues (DOING NON-FLEX-BB INSTEAD)
        //yeah I know this is kind of a weird criterion
        //but it's a simple example of a packing criterion that might sometimes be useful
        
        //only works on singles, which may not touch any RC of the non-bb-flex res
        //(whose RCs we thus must iterate over)
        
        if(tup.size()>1)
            throw new RuntimeException("ERROR: Flex-res contact pruning only makes sense for singles");
        
        RCClass cl = tup.get(0);
        if(!mutableRes.contains(cl.res))
            return false;
        
        ResSphereModel rcModel = stericCache.getSphereModel(cl, true);
        if(rcModel==null)//let invalid conf be pruned
            return true;
        HashMap<Residue,ResPairClashExemptions> exemptionMap = stericCache.getShellExemptions(cl);//what shell residue may have clash exemptions
        
        //first, the shell is nonbbflex so touching it rules out pruning
        for(int shellPos=0; shellPos<shellResidues.size(); shellPos++){
            Residue res = shellResidues.get(shellPos);
            ResSphereModel sm = stericCache.getShellSphereModel(shellPos);
            ResPairClashExemptions exemptions = null;
            if(exemptionMap.containsKey(res))
                exemptions = exemptionMap.get(res);
            if(rcModel.clashesWith(sm, -contactBuffer, exemptions))//negative buffer...don't want to miss touching
                return false;
        }
        
        //now check all RC's at non-BB-flex positions
        for(int pos : nonBBFlexPos){
            for( RCClass cl2 : classesAtRes.get(pos).singletonClassLookup.values() ){
                ResSphereModel rcModel2 = stericCache.getSphereModel(cl2, true);
                if(rcModel2!=null){
                    ResPairClashExemptions exemptions = stericCache.getRCPairExemptions(cl, cl2);
                    if( rcModel.clashesWith(rcModel2, -contactBuffer, exemptions) )
                        return false;
                }
            }
        }
        
        
        return true;
    }
    

}
