/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tupexp;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.TupleEnumerator;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class LUTECompStatus implements Serializable {
    
    String fileName;//for saving to disk
    transient SearchProblem sp;//this will be put in by the SearchProblem when loading
    
    
    //Storing best to date
    EnergyMatrix curBestMatrix = null;
    double curBestResid = Double.POSITIVE_INFINITY;
    
    //DEBUG!!  For now (PLUG testing), will go through a fixed set of expansions:
    //pairwise, 2-partner triples, 5-partner triples, PLUG clique triples, 
    //triples with 2 PLUG pairs, triples with 2 plug paris + clique quads
    //in the future may want to cut off if residual good enough
    int curStep = 0;
    final static int numSteps = 5;//6;//clique quads seem rare and not too helpful, and they make DEE complicated
    
    
    //Store expander if there's an expansion in progress, delete once it's done so we 
    //don't try to resume with the wrong one
    TupleExpander curTupleExpander = null;
    
    
    long lastSaveTime = 0;
    final static int saveInterval = 1800000;//30 min in ms.  When scoring samples, save this often
    
    
    public LUTECompStatus(SearchProblem sp, String fn){
        this.sp = sp;
        fileName = fn;
    }
    
    public void setSp(SearchProblem s){
        sp = s;
    }
    
    
    void save(TupleExpander te){
        //update the cached expander & save to disk
        curTupleExpander = te;
        saveToDisk();
    }
    
    
    boolean isTimeToSave(){
        return System.currentTimeMillis() > lastSaveTime + saveInterval;
    }
    
    private void saveToDisk(){
        ObjectIO.writeObject(this, fileName);
        lastSaveTime = System.currentTimeMillis();
    }
    
    private void deleteFromDisk(){
        //delete saved file.  Important since it might be big and also it would mess
        //up any subsequent LUTE runs
        new File(fileName).delete();
    }
    
    
    public EnergyMatrix precomputeLUTE(){
        if(curStep!=0)
            System.out.println("RESUMING LUTE AT STEP "+curStep);
        else
            System.out.println("STARTING LUTE PRECOMPUTATION AT FIRST STEP");
        
        int numStepsHere = (sp.plugMat==null) ? 3 : numSteps;//DEBUG!!
        for(; curStep<numStepsHere; curStep++){
            
            saveToDisk();//this is a good place to save so we can resume starting a new step
            double resid;
            if(curTupleExpander==null){//starting step, not resuming
                curTupleExpander = new ConfETupleExpander(sp,this);
                resid = curTupleExpander.calcExpansion(tuplesToFit());
            }
            else {//resuming step
                System.out.println("LOADING SAVED TUPLE EXPANDER AT STEP "+curStep);
                curTupleExpander.scoreUnscoredSamples();
                resid = curTupleExpander.calcExpansionFromCurSamples();
            }
            
            System.out.println("RESIDUAL FOR LUTE STEP "+curStep+" : "+resid);
            
            if(resid<curBestResid){
                System.out.println("RESIDUAL IS BEST SO FAR");
                curBestResid = resid;
                curBestMatrix = curTupleExpander.getEnergyMatrix();
            }
            
            //step is completed
            curTupleExpander = null;//no longer relevant
        }
        
        deleteFromDisk();
        
        System.out.println("FINISHED LUTE PRECOMPUTATION.  BEST RESID: "+curBestResid);
        return curBestMatrix;
    }
    
    
    ArrayList<RCTuple> tuplesToFit(){
        //always include pairs
        //based on TupExpChooser
        TupleEnumerator tupEnum = new TupleEnumerator(sp.pruneMat,sp.emat,sp.plugMat,sp.confSpace.numPos);
        
        ArrayList<RCTuple> ans = tupEnum.enumerateUnprunedTuples(2);//start with pairs
        if(ans.isEmpty())//no pairs...indicates <2 flexible positions (or everything pruned, which will be detected)
            return tupEnum.enumerateUnprunedTuples(1);//so add singles
        
        //now add higher-order stuff (except for step 0 = pairs)
        if(curStep==1 || curStep==2){
            int numPartners = (curStep==1) ? 2 : 5;
            ArrayList<ArrayList<Integer>> topPositionTriples = tupEnum.topPositionTriples(numPartners);
            ArrayList<RCTuple> topTriples = tupEnum.enumerateUnprunedTuples(topPositionTriples);
            ans.addAll(topTriples);
        }
        else if(curStep==3)
            ans.addAll(tupEnum.enumeratePLUGCliqueTriples());
        else if(curStep>3){
            ans.addAll(tupEnum.enumeratePLUG2PairTriples());
            if(curStep==5)
                ans.addAll(tupEnum.enumeratePLUGCliqueQuads());
        }
        
        System.out.println("PREPARED "+ans.size()+" tuples to fit for step "+curStep
                +": "+curStepDescription());
        return ans;
    }
    
    String curStepDescription(){
        switch(curStep){
            case 0:
                return "PAIRWISE";
            case 1:
                return "2-PARTNER TRIPLES";
            case 2:
                return "5-PARTNER TRIPLES";
            case 3:
                return "PLUG CLIQUE TRIPLES";
            case 4:
                return "PLUG 2-PAIR TRIPLES";
            case 5:
                return "PLUG 2-PAIR TRIPLES + CLIQUE QUADS";
            default:
                throw new RuntimeException("ERROR: There are only 6 steps in LUTECompStatus now");
        }
    }
    
    
}
