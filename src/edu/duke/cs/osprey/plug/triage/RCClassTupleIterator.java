/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug.triage;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;

/**
 *
 * @author mhall44
 */



class MultiCounter {
    //incrementable counter with multiple indices; limits specified
    int[] counters;
    ArrayList<Integer> limits;
    boolean forceDesc;//only consider counters in descending order (for repeat-free unordered tuples)
    boolean exhausted = false;
    
    MultiCounter(boolean forceDesc, Integer... lim){//initialize to right before (0,0,...)
        this.forceDesc = forceDesc;
        limits = new ArrayList(Arrays.asList(lim));
        counters = new int[limits.size()];
        Arrays.fill(counters, 0);
        counters[limits.size()-1] = -1;
    }
    
    int get(int i){
        return counters[i];
    }
    
    private boolean counterInRange(int posInTup){
        //is given position in counter within allowed range
        if(forceDesc && posInTup>0){
            if(counters[posInTup] >= counters[posInTup-1])
                return false;
        }
        return counters[posInTup] < limits.get(posInTup);
    }
    
    void increment(){
        for(int posInTup=limits.size()-1; posInTup>=0; posInTup--){
            counters[posInTup]++;
            if(counterInRange(posInTup))//successful increment
                return;
            else//reset this counter, increment slower-moving one
                counters[posInTup] = 0;
        }
        //if we got here we have reached the limit for all dimensions
        exhausted = true;
    }
}



public class RCClassTupleIterator implements Iterator<RCClassTuple> {     
    RCClass.ClassType[] classTypes;
    int tupSize;
    ArrayList<RCClassification> classesAtRes;
    PruningLookup prooningLookup;
    
    RCClassTuple nextTup = null;
    MultiCounter nextTupPos;
    MultiCounter nextTupRCClasses;
    
    //constructor for 1-tuple
    RCClassTupleIterator(RCClass.ClassType classType, 
            ArrayList<RCClassification> classesAtRes, PruningLookup prooningLookup) {
        //intra constructor
        classTypes = new RCClass.ClassType[] {classType};
        this.classesAtRes = classesAtRes;
        this.prooningLookup = prooningLookup;
        
        tupSize = 1;
        nextTupPos = new MultiCounter(true,classesAtRes.size());
        nextTupRCClasses = new MultiCounter(false,0);//empty counter
        findNextTup();//get next tuple ready
    }
    
    RCClassTupleIterator(RCClass.ClassType classType, RCClass.ClassType classType2,
            ArrayList<RCClassification> classesAtRes, PruningLookup prooningLookup) {
        //intra constructor
        classTypes = new RCClass.ClassType[] {classType,classType2};
        this.classesAtRes = classesAtRes;
        this.prooningLookup = prooningLookup;
        
        tupSize = 2;
        nextTupPos = new MultiCounter(true,classesAtRes.size(),classesAtRes.size());
        nextTupRCClasses = new MultiCounter(false,0,0);//empty counter
        findNextTup();//get next tuple ready
    }
    
    
    private void findNextTup(){
        while(true){//repeat til find a good tuple or exhaust possibilities
            nextTupRCClasses.increment();
            while(nextTupRCClasses.exhausted){
                nextTupPos.increment();
                if(nextTupPos.exhausted){//done iterating
                    nextTup=null;
                    return;
                }
                buildRCClassCounter();
                nextTupRCClasses.increment();
            }
            buildNextTup();
            if( ! prooningLookup.checkIfPruned(nextTup) ){
                return;//nextTup is good
            }
        }
    }
    
    
    private void buildRCClassCounter(){
        //build RC class counter based on nextTupPos
        Integer[] limits = new Integer[tupSize];
        for(int i=0; i<tupSize; i++){
            limits[i] = classesAtRes.get(nextTupPos.get(i)).numClassesOfType(classTypes[i]);
        }
        nextTupRCClasses = new MultiCounter(false,limits);
    }
    
    private void buildNextTup(){
        //build nextTup from nextTupPos and nextTupRCClasses
        RCClass[] classes = new RCClass[tupSize];
        for(int i=0; i<tupSize; i++){
            classes[i] = classesAtRes.get(nextTupPos.get(i)).get(classTypes[i],nextTupRCClasses.get(i));
        }
        nextTup = new RCClassTuple(classes);
    }

    @Override
    public boolean hasNext() {
        return (nextTup!=null);
    }

    @Override
    public RCClassTuple next() {
        RCClassTuple ans = nextTup;
        findNextTup();
        return ans;
    }
    
}
