/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug.triage;

import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * for recording/looking up what RCClassTuples are pruned
 * 
 * @author mhall44
 */
public class PruningLookup {
    
    //DEBUG!!!  This may be slo
    //but we can just record what all classes are pruned
    //lookup can go up hierarchy (parent pruned implies child pruned)
    
    HashSet<RCClass> proondClasses = new HashSet<>();
    HashMap<RCClass,HashSet<RCClass>> proondPairz = new HashMap<>();
    
    void markAsPruned(RCClassTuple tup){
        switch(tup.size()){
            case 1:
                proondClasses.add(tup.get(0));
                break;
            case 2:
                RCClass key = tup.get(0);
                if( ! proondPairz.containsKey(key) )
                    proondPairz.put(key, new HashSet<>());
                proondPairz.get(key).add(tup.get(1));
                break;
            default:
                throw new RuntimeException("ERROR TUPLE TOO BIG");     
        }
    }
    
    boolean checkIfPruned(RCClassTuple tup){
        //go through pair reordering, checking parent classes, etc.
        if(checkClassPruned(tup.get(0)))
            return true;
        if(tup.size()>1){
            if(checkClassPruned(tup.get(1)))
                return true;
            if(checkPairPruned(tup.get(0),tup.get(1)))
                return true;
            if(checkPairPruned(tup.get(1),tup.get(0)))//may be recorded in either order...
                return true;
            
            if(tup.size()>2)
                throw new RuntimeException("ERROR TUPLE TOO BIG");     
        }
        return false;
    }
    
    
    boolean checkClassPruned(RCClass cl){
        if(proondClasses.contains(cl))
            return true;
        if(cl.parent!=null)
            if(checkClassPruned(cl.parent))
                return true;
        return false;
    }
    
    private boolean checkPairPruned(RCClass cl1, RCClass cl2){
        //check if the given pair, in the given order, is pruned 
        //DON'T USE TO CHECK IF THE PAIR IS PRUNED IN ANY WAY
        if(checkPairPrunedHelper(cl1,cl2))
            return true;
        if(cl1.parent!=null)
            if(checkPairPruned(cl1.parent,cl2))
                return true;
        return false;
    }
    
    
    private boolean checkPairPrunedHelper(RCClass cl1, RCClass cl2){
        //Checks if (cl1, cl2) is pruned, with cl1 pruned explicitly (not superclass)
        //this is really only meaningful as a sub-checkin checkPairPruned
        if(proondPairz.containsKey(cl1)){
            if(proondPairz.get(cl1).contains(cl2))
                return true;
            if(cl2.parent!=null){
                if(checkPairPrunedHelper(cl1,cl2.parent))
                    return true;
            }
        }
        return false;
    }
    
   
    
}
