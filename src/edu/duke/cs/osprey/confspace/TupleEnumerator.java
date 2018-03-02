/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;
import java.util.stream.Collectors;

/**
 *
 * @author mhall44
 */
public class TupleEnumerator {
    //Makes lists of residue or RC tuples for pruning, tuple expansion, etc.
    
    PruningMatrix pruneMat;//for figuring out which are unpruned
    EnergyMatrix emat;//for figuring out which tuples interact strongly
    PolytopeMatrix plugMat;//ditto.  Can be null if needed
    int numPosTotal;//total number of positions

    
    public TupleEnumerator(PruningMatrix pruneMat, EnergyMatrix emat, PolytopeMatrix plugMat, int numPosTotal) {
        this.pruneMat = pruneMat;
        this.emat = emat;
        this.plugMat = plugMat;
        this.numPosTotal = numPosTotal;
    }
    
    
    
    public ArrayList<RCTuple> enumerateUnprunedTuples(int numPosInTuple){
        //enumerate all unpruned tuples of the specified size, at any position
        //e.g. if numPosInTuple==2 then enumerate all unpruned pairs
        ArrayList<ArrayList<Integer>> posTupleCand = allPositionTuples(numPosInTuple);
        return enumerateUnprunedTuples( posTupleCand );
    }
    
    
    public ArrayList<RCTuple> enumerateUnprunedTuples( ArrayList<ArrayList<Integer>> posTupleCand  ){
        //Enumerate all unpruned RC tuples at the specified position tuples
        
        ArrayList<RCTuple> allCandidates = new ArrayList<>();
        
        for(ArrayList<Integer> posTuple : posTupleCand){
            ArrayList<RCTuple> candidatesAtPos = pruneMat.unprunedRCTuplesAtPos(posTuple);
            allCandidates.addAll(candidatesAtPos);
        }
        
        return allCandidates;
    }
    
    
    public ArrayList<ArrayList<Integer>> allPositionTuples(int numPosInTuple){
        //get all possible sets of positions, of the specified size
        //(sets represented as tuples in ascending order, and must have all distinct positions)

        ArrayList<ArrayList<Integer>> ans = new ArrayList<>();
        
        if(numPosInTuple==1){
            for(int pos=0; pos<numPosTotal; pos++){
                ArrayList<Integer> singleton = new ArrayList<>();
                singleton.add(pos);
                ans.add(singleton);
            }
        }
        else {
            ArrayList<ArrayList<Integer>> reducedTups = allPositionTuples(numPosInTuple-1);
            
            for(ArrayList<Integer> redTup : reducedTups){
                int lastPosInTup = redTup.get(numPosInTuple-2);
                
                for(int newPos=lastPosInTup+1; newPos<numPosTotal; newPos++){
                    ArrayList<Integer> fullTup = (ArrayList<Integer>)redTup.clone();
                    fullTup.add(newPos);
                    ans.add(fullTup);
                }
            }
        }
        
        return ans;
    }
    
    
    public ArrayList<ArrayList<Integer>> topPositionTriples(int numPartners){
        //Get only the strongest-interacting position triples
        
        //enumerate all interactions first...  (pretty quick because not looping over RCs)
        ArrayList<ArrayList<Integer>> allPositionTriples = allPositionTuples(3);
        ArrayList<ArrayList<Integer>> ans = new ArrayList<>();
        
        boolean strongInteraction[][] = strongInteractingPairs(numPartners);
        
        //pick out the strong ones
        for(ArrayList<Integer> posTriple : allPositionTriples) {
            int strongInteractionCount = 0;
            for(int i=0; i<3; i++){
                int p1 = posTriple.get(i);
                int p2 = posTriple.get((i+1)%3);
                if(strongInteraction[p1][p2])
                    strongInteractionCount++;
           }
           
           if(strongInteractionCount>=2)
               ans.add(posTriple);
        }
        
        return ans;
    }
    
    
    public boolean[][] strongInteractingPairs(int numPartners){
        //For each position pos, the numPartners other positions that interact most strongly with pos
        //are counted as "strong" interactions
        //return a matrix of which pairwise interactions are "strong"
        //Using lower-bound energy matrix should capture these trends without possible noise from initial bad tup-exp fit...
        
        int numPos = emat.getNumPos();
        double strongestPairE[][] = emat.topPairwiseInteractions();
        
        //Next, use these to figure out what are the numPartners top interaction partners for each residue
        //make a matrix of "strong" interactions based on this
        boolean strongInteraction[][] = new boolean[numPos][numPos];
        
        for(int pos=0; pos<numPos; pos++){
            
            PriorityQueue<TupE> posTop = new PriorityQueue<>();
            
            for(int pos2=0; pos2<numPos; pos2++){
                if(pos!=pos2){
                    posTop.add( new TupE(new RCTuple(pos2,-1), strongestPairE[pos][pos2]) );
                    if(posTop.size()>numPartners)//throw out subpar interactions
                        posTop.poll();
                }
            }
            
            for(TupE top : posTop){
                int pos2 = top.tup.pos.get(0);
                strongInteraction[pos][pos2] = true;
                strongInteraction[pos2][pos] = true;
            }
        }
        
        return strongInteraction;
    }
    
    
    public ArrayList<RCTuple> clique2PairTriples(int numTriplesPerType){
        //enumerate RC triples based on containing either two or three strong pairwise interactions in emat
        //we take the best triples for each measure (numTriplesPerType of each)
        
        ArrayList<RCTuple> tupList = new ArrayList<>();
        
        PriorityQueue<TupE> topCliques = new PriorityQueue<>();//top triples in terms of all pairs strong
        PriorityQueue<TupE> top2Pair = new PriorityQueue<>();//top triples in terms of having 2 strong pairs
        
        ArrayList<RCTuple> possibleTriples = enumerateUnprunedTuples(3);//consider all unpruned triples
        
        //go through possible triples
        for(RCTuple triple : possibleTriples) {
            double[] pairAbsE = new double[3];
            double minAbsE = Double.POSITIVE_INFINITY;
            double max2PairE = Double.NEGATIVE_INFINITY;
            
            for(int i=0; i<3; i++){
                RCTuple pair = triple.subtractMember(i);
                pairAbsE[i] = Math.abs( emat.getPairwise(pair.pos.get(0), pair.RCs.get(0), 
                        pair.pos.get(1), pair.RCs.get(1)) );
                
                minAbsE = Math.min(minAbsE,pairAbsE[i]);
            }
            
            for(int i=0; i<3; i++)
                max2PairE = Math.max( max2PairE, Math.min(pairAbsE[(i+1)%3],pairAbsE[(i+2)%3]) );
            
            
            topCliques.add( new TupE(triple,minAbsE) );
            top2Pair.add( new TupE(triple,max2PairE) );
            
            //now throw out the weakest-interacting triples if we have too many
            if(topCliques.size()>numTriplesPerType)
                topCliques.poll();//CHECK DIRECTIONS
            if(top2Pair.size()>numTriplesPerType)
                top2Pair.poll();
        }
        
        
        for(TupE tupe : topCliques)
            tupList.add(tupe.tup);
        for(TupE tupe : top2Pair)
            tupList.add(tupe.tup);
        
        return tupList;
    }

    public EnergyMatrix getEmat() {
        return emat;
    }

    public void setEmat(EnergyMatrix emat) {
        this.emat = emat;
    }
    
    
    public ArrayList<RCTuple> enumeratePLUGContactPairs(){
        //Unpruned pairs that have contacts detected by PLUG
        List<RCTuple> ans = enumerateUnprunedTuples(2).stream()
                .filter(tup->plugMat.getTupleValue(tup).hasNonBoxConstr())
                .collect(Collectors.toList());
        
        return new ArrayList(ans);
    }
    
    public ArrayList<RCTuple> enumeratePLUGCliqueTriples(){
        //Unpruned triples in which all 3 pairs have contacts detected by PLUG
        return enumeratePLUGTriples(true);
    }
    
    public ArrayList<RCTuple> enumeratePLUG2PairTriples(){
        //Unpruned triples in which >= 2 of the 3 pairs have contacts detected by PLUG
        return enumeratePLUGTriples(false);
    }
    
    public ArrayList<RCTuple> enumeratePLUGTriples(boolean clique){
        ArrayList<RCTuple> ans = new ArrayList<>();
        for(RCTuple contactPair : enumeratePLUGContactPairs()){
            for(int pos=0; pos<Math.min(contactPair.pos.get(0), contactPair.pos.get(1)); pos++){
                for(int rc=0; rc<pruneMat.getNumConfAtPos(pos); rc++){
                    RCTuple triple = contactPair.addRC(pos, rc);
                    if(!pruneMat.isPruned(triple)){
                        RCTuple pair2 = triple.subtractMember(0);
                        boolean pair2Contact = plugMat.getTupleValue(pair2).hasNonBoxConstr();
                        RCTuple pair3 = triple.subtractMember(1);
                        boolean pair3Contact = plugMat.getTupleValue(pair3).hasNonBoxConstr();
                        if(clique){
                            if(pair2Contact&&pair3Contact){
                                ans.add(triple);
                            }
                        }
                        else {
                            if(pair2Contact||pair3Contact){
                                ans.add(triple);
                            }
                        }
                    }
                }
            }
        }
        return ans;
    }
    
    public ArrayList<RCTuple> enumeratePLUGCliqueQuads(){
        //Unpruned quads in which all pairs have PLUG contact
        ArrayList<RCTuple> ans = new ArrayList<>();
        for(RCTuple triple : enumeratePLUGCliqueTriples()){
            int minPos = Collections.min(triple.pos);//avoid redundancy
            for(int pos=0; pos<minPos; pos++){
                for(int rc=0; rc<pruneMat.getNumConfAtPos(pos); rc++){
                    RCTuple quad = triple.addRC(pos, rc);
                    if(!pruneMat.isPruned(quad)){
                        boolean isClique = true;
                        //see if all pairs in the quad involving the new pos interact
                        for(int a=0; a<3; a++){
                            for(int b=0; b<2; b++){
                                RCTuple otherPair = quad.subtractMember(a).subtractMember(b);
                                if(!plugMat.getTupleValue(otherPair).hasNonBoxConstr()){
                                    isClique = false;
                                    break;
                                }
                            }
                        }
                        if(isClique)
                            ans.add(quad);
                    }
                }
            }
        }
        return ans;
    }
    
    
    
    
}
