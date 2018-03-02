/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.aibis;

import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * This discrete DOF represents the different RC choices
 * for a particular AA type at a particular position
 * 
 * @author mhall44
 */
public class RCIntegrableDOF implements IntegrableDOF {

    int pos;
    HashMap<String,ArrayList<Integer> > RCsByAAType;
    HashMap<Integer,String> AATypeByRC;
    
    
    public RCIntegrableDOF(int pos, ArrayList<RC> rcListAtPos){
        this.pos = pos;
        RCsByAAType = new HashMap<>();
        AATypeByRC = new HashMap<>();
        for(int rcIndex=0; rcIndex<rcListAtPos.size(); rcIndex++){
            String AAType = rcListAtPos.get(rcIndex).AAType;
            AATypeByRC.put(rcIndex, AAType);
            if( ! RCsByAAType.containsKey(AAType) )
                RCsByAAType.put(AAType, new ArrayList<>());
            RCsByAAType.get(AAType).add(rcIndex);
        }
    }
    
    
    @Override
    public double boltzmannIntegratedEnergy(SparseLinearEnergy baseE, ConfSample samp) {
        //integral is just a sum in this case
        String AAType = AATypeByRC.get(samp.getRC(pos));
        double q = 0;
        for(int rc : RCsByAAType.get(AAType)){
            RCTuple altTup = samp.vox.subtractPos(pos).addRC(pos,rc);
            ConfSample altSamp = new ConfSample(altTup,samp.contDOFValsArr);//likely won't be any contDOFVals but maybe a bb motion or something...
            q += Math.exp(-baseE.evalEnergy(altSamp)/RT);
        }
        return -RT*Math.log(q);
    }

    @Override
    public FeatureSet featureSetForInteg(FeatureSet baseFeat) {
        //assumes these RCs already have just constant terms...
        //set all features for RCs of this AAType other than the first to be the same as the first feature
        FeatureSet ans = new FeatureSet(baseFeat.confSpace);
        
        HashMap<String,Integer> mergedOffsets = new HashMap<>();//merge offsets for tuples that differ
        //only in having different RCs at pos (from the specified set)
        
        HashMap<Integer,RCTuple> previouslySeenOffsets = new HashMap<>();//for detecting redundant offsets in baseFeat
        //and replicating them here.  Maps the offset to baseFeat to the RCTuple for which it was first seen
                
        for(RCTuple tup : baseFeat.unprunedFeatTuples()){
            DenseFeatureSet tupFS = baseFeat.featMatrix.getTupleValue(tup);
            int baseOffset = baseFeat.featOffsets.getTupleValue(tup);
            boolean redundant = previouslySeenOffsets.containsKey(baseOffset);
            if(!redundant)
                previouslySeenOffsets.put(baseOffset,tup);
            
            if(tup.pos.contains(pos)){//need to merge
                String key = makeKey(tup);
                if(mergedOffsets.containsKey(key)){//already have a tuple that is like tup except different RC at same AA type
                    ans.addRedundantDenseFeatureSet(tupFS, tup, mergedOffsets.get(key));
                }
                else {//first tuple with this AA type: add as usual, use this offset for the others
                    if(redundant)
                        ans.addRedundantDenseFeatureSet(tupFS, tup, ans.featOffsets.getTupleValue(previouslySeenOffsets.get(baseOffset)));
                    else
                        ans.addDenseFeatureSet(tupFS, tup);
                    mergedOffsets.put(key, ans.featOffsets.getTupleValue(tup));
                }
            }
            else{//no need to merge
                if(redundant)//make tup redundant with the same stuff it was redundant with in baseFeat
                    ans.addRedundantDenseFeatureSet(tupFS, tup, ans.featOffsets.getTupleValue(previouslySeenOffsets.get(baseOffset)));
                else
                    ans.addDenseFeatureSet(tupFS, tup);
            }
        }
        
        return ans;
    }
    
    
    String makeKey(RCTuple tup){
        //A key for a tuple, shared by tuples with the same AA type at pos,
        //but otherwise unique for tuples
        RCTuple redTup = tup.subtractPos(pos);
        String AAType = AATypeByRC.get(tup.RCAtPos(pos));
        return AAType+";"+redTup.stringListing();
    }

    @Override
    public String description() {
        return "RCIntegrableDOF "+pos;
    }
}
