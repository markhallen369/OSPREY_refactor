/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug.triage;

import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.plug.RCTuplePolytope;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * 
 * This is going to be a more flexible pruning system
 * Early steps need to be very efficient
 * Needed to cut down huge, e.g. FVB, conf spaces
 * intended for huge conf spaces most of which can be pruned geometrically somehow
 * 
 * Hopefully this will cover a reasonable area
 * but it difficult, can iteratively expand (e.g.) ligand BB flex til get nice
 * space of things that fit (maybe even try LUTE fit to see how that goes)
 * 
 * Pruning is insufficient to really "triage" anything
 * and rama not even convex in dof space
 * So "lions" look infeasible...more like lion to myself about the feasibility of searching
 * big bb spaces systematically
 * see RamaLinConstr.java
 * 
 * Without Lions or substantial geometric pruning, we are heavily reliant on good energy estimation,
 * along with marginal conformational search advances + manual search, etc.
 * for progress.  
 * 
 * @author mhall44
 */
public class PLUGTriage {
    
    SearchProblem sp;
    
    PruningLookup curPruned = new PruningLookup();//various stuff (not in dense pruning matrix) currently pruned
    ArrayList<RCClassification> classesAtRes = new ArrayList<>();//all the classes available at each res
    
    
    //Pruneables:
    //RC
    //RC pair
    //Truncation/isosteric class equivalents of any of these
    
    
    //HOWEVER, we can have a pruning layer that only considers clashes at the central conf
    //or maybe other pts quickly identified!
    //(since dist opt seems to be bottleneck)
    //for full confs of course will want full mat since that's ultimately poly time
    //call this QUICKCLASHES
    //then regular mat is CLASHES, at which point we are building a PolytopeMatrix
    //hmm or at least the intra/shell parts
    
    
    //clash and Rama assumptions are fairly straightforward,
    //but there are a couple other assumptions that are more at will
    boolean demandPackedVol = false;//Some form of this is needed for biophysical feasibility
    //but will need to specify what "dry zone" we want packed
    boolean demandHBond = false;//This is just a design decision, specify what target AA need an H-bond
    //DEBUG!!!  Also should adjust clash def to exclude close H-bonds
    boolean demandFlexResContact = true;//This is largely a design decision too
    //DEBUG!!!  Will demand mutable res contact flexible res as a sort of heuristic
    //that ensures decent packing (though not all designs will want this)
    
    HashMap<RCClass,RCTuplePolytope> classQuickPolytopes = new HashMap<>();
            //for use in QUICKCLASH
    
    HashSet<Residue> mutableRes;
    
    CentralStericCache stericCache;
    
    static final boolean origBBLinks = true;//use links for altered BB classes
    
    public PLUGTriage(SearchProblem sp){
        this.sp = sp;
        for(int pos=0; pos<sp.confSpace.numPos; pos++){
            classesAtRes.add( new RCClassification(sp.confSpace.posFlex.get(pos),sp.confSpace.mutDOFs.get(pos),origBBLinks) );
        }
        stericCache = new CentralStericCache(sp.confSpace,classesAtRes,sp.shellResidues);
        pickMutableRes();
    }
    
    
    private void pickMutableRes(){
        mutableRes = new HashSet<>();
        for(int pos=0; pos<sp.confSpace.numPos; pos++){
            boolean mutable = false;
            ArrayList<RC> RCs = sp.confSpace.posFlex.get(pos).RCs;
            String baseAA = RCs.get(0).AAType;
            for(RC rc : RCs){
                if( ! rc.AAType.equalsIgnoreCase(baseAA) ){
                    mutable = true;
                }
            }
            if(mutable)
                mutableRes.add(sp.confSpace.posFlex.get(pos).res);
        }
    }
    
    
    public void doTriage(){
        
                //NEED TO BUILD A SUITABLE BUT SMALL TEST SYSTEM LIKE ONE FLOATING AA TRYING TO BIND STUFF

        
        //The number of available pruneables of each kind will strictly decrease
        //So, let's do pruning steps in order of per-pruneable time X number of TOTAL pruneables available
        //this time X frac of pruneables still standing = actual time cost (assuming negligible iterating time)
        
        //note: only worth starting matrix comp (lin constr to use at higher order)
        //when ready to prune higher order.   
        
        ArrayList<GeomPruningStep> pruningSteps = new ArrayList<>();
        //may need to sort this by est cost
        
        queueUpPruningStepLevels(PruningCriteria.STERIC, 1, pruningSteps);
        
        //ok first let's try to cut down the number of RCs some
        //start w/ basic steric pruning (based on huge overlap at orig conf, thus fast)
        /*doPruningStep(STERIC, BBCLASS, BASIC);//BASIC = only checking vox center
        doPruningStep(STERIC, BETACLASS, BASIC);
        doPruningStep(STERIC, GAMMACLASS, BASIC);
        doPruningStep(STERIC, ROT, BASIC);*/
        
        if(demandFlexResContact){
            //do the same as steric but for flex-res contact...this is per-res
            //this is the only pruning we will need for this (somewhat weak) criterion!  
            queueUpPruningStepLevels(PruningCriteria.FLEXRESCONTACT, 1, pruningSteps);
        }
        
        
        //ok now flesh out the steric pruning
        //rama only needed for bbclass
        queueUpPruningStepLevel(PruningCriteria.RAMA, RCClass.ClassType.BBCLASS, RCClass.ClassType.BBCLASS, pruningSteps);

        queueUpPruningStepLevels(PruningCriteria.STERIC, 2, pruningSteps);
        
        
        /*
        //DEBUG!!!! implement this!!!
        
        //now go on to quickclash pruning...
        queueUpPruningStepLevels(PruningCriteria.QUICKCLASH, 1, pruningSteps);
        queueUpPruningStepLevels(PruningCriteria.QUICKCLASH, 2, pruningSteps);
        //WILL NEED TO PRECOMPUTE FOR THIS
        */
        
        System.out.print("Number of RCs by pos: ");
        for(int pos=0; pos<sp.confSpace.numPos; pos++)
            System.out.print(sp.confSpace.posFlex.get(pos).RCs.size()+" ");
        System.out.println();
        
        
        for(GeomPruningStep step : pruningSteps){
            step.prune();
            pruneSinglesUsingPairs();
            //DEBUG!!  also prune classes using singles??  (to avoid having to iterate over them here)
            
            //DEBUG!!!
            System.out.println("STEP COMPLETED.  classes pruned: "+curPruned.proondClasses.size());
            int numPairsPruned = curPruned.proondPairz.size();
            System.out.println("Pairs pruned: "+numPairsPruned);
            System.out.print("RCs remaining by pos: ");
            ArrayList<Integer> remaining = new ArrayList<>();
            int pairsRemaining = 0;//assumes no pairwise pruning yet...only show if so
            for(int pos=0; pos<sp.confSpace.numPos; pos++){
                PositionConfSpace pcSpace = sp.confSpace.posFlex.get(pos);
                int remainingAtPos = pickUnprunedRCs(pcSpace.RCs, pos).size();
                System.out.print(remainingAtPos+" ");
                remaining.add(remainingAtPos);
                for(int pos2=0; pos2<pos; pos2++)
                    pairsRemaining += remainingAtPos*remaining.get(pos2);
            }
            System.out.println();
            if(numPairsPruned==0)
                System.out.println("Pairs remaining: "+pairsRemaining);
        }
        
        //full LP pruning, polytope mat construction can take place in the normal way outside triage
        //so cut out the bad RCs and make a (probably res-wise sparse) pair pruning matrix 
        removeBadRCs();
        sp.pruneMat = new PruningMatrix(sp.confSpace, Double.POSITIVE_INFINITY);//valid regardless of ival
        prunePairs();
    }
    
    private void pruneSinglesUsingPairs(){
        while(true){
            boolean prunedSomething = false;
            ArrayList<ArrayList<RC> > unprunedRCs = new ArrayList<>();
            //int pairsRemaining = 0;//DEBUG!! takes too long
            for(int pos=0; pos<sp.confSpace.numPos; pos++){
                PositionConfSpace pcSpace = sp.confSpace.posFlex.get(pos);
                unprunedRCs.add( pickUnprunedRCs(pcSpace.RCs, pos) );
            }
            
            for(int pos1=0; pos1<sp.confSpace.numPos; pos1++){
                for(RC rc1 : unprunedRCs.get(pos1)){
                    RCClass cl1 = classesAtRes.get(pos1).singletonClassLookup.get(rc1);
                    for(int pos2=0; pos2<sp.confSpace.numPos; pos2++){
                        if(pos1!=pos2){//an unpruned RC must have unpruned partners at every other res
                            boolean confPossible = false;

                            for(RC rc2 : unprunedRCs.get(pos2)){
                                RCClass cl2 = classesAtRes.get(pos2).singletonClassLookup.get(rc2);
                                if(!curPruned.checkClassPruned(cl2)){
                                    RCClassTuple pair = new RCClassTuple(cl1,cl2);
                                    if(!curPruned.checkIfPruned(pair)){
                                        confPossible = true;
                                        break;
                                        //pairsRemaining++;
                                    }
                                }
                            }
                            if(!confPossible){
                                curPruned.markAsPruned(new RCClassTuple(cl1));
                                break;
                            }
                        }
                    }
                }
            }
            
            if(!prunedSomething){
                //pairsRemaining /= 2;//double-counted
                //System.out.println("Consolidated pairs pruning into singles.  Pairs remaining: "+pairsRemaining);
                return;
            }
            //if pruned something go do another cycle, because might be able to prune something else
        }
    }
    
    private void removeBadRCs(){
        //assuming no energy calc yet in sp, remove the pruned RCs
        for(int pos=0; pos<sp.confSpace.numPos; pos++){
            PositionConfSpace pcSpace = sp.confSpace.posFlex.get(pos);
            pcSpace.RCs = pickUnprunedRCs(pcSpace.RCs, pos);
        }
    }
    
    private ArrayList<RC> pickUnprunedRCs(ArrayList<RC> rcList, int pos){
        ArrayList<RC> ans = new ArrayList<>();
        for(RC rc : rcList){
            if( ! curPruned.checkClassPruned( classesAtRes.get(pos).singletonClassLookup.get(rc) ) ){
                ans.add(rc);
            }
        }
        return ans;
    }
    
    
    private void prunePairs(){
        //Prune infeasible pairs
        for(int pos=0; pos<sp.confSpace.numPos; pos++){
            ArrayList<RC> rcList1 = sp.confSpace.posFlex.get(pos).RCs;
            for(int rc=0; rc<rcList1.size(); rc++){
                RCClass cl1 = classesAtRes.get(pos).singletonClassLookup.get(rcList1.get(rc));
                for(int pos2=0; pos2<pos; pos2++){
                    ArrayList<RC> rcList2 = sp.confSpace.posFlex.get(pos2).RCs;
                    for(int rc2=0; rc2<rcList2.size(); rc2++){
                        RCClass cl2 = classesAtRes.get(pos2).singletonClassLookup.get(rcList2.get(rc2));
                        if(curPruned.checkIfPruned(new RCClassTuple(cl1,cl2))){
                            sp.pruneMat.markAsPruned(new RCTuple(pos,rc,pos2,rc2));
                        }
                    }
                }
            }
        }
    }
    
    
    
    
    static enum PruningCriteria {
        STERIC, QUICKCLASH, FLEXRESCONTACT, RAMA
    }
    
    void queueUpPruningStepLevels(PruningCriteria criterion, int arity, ArrayList<GeomPruningStep> pruningSteps){
        for(RCClass.ClassType ctype : RCClass.ClassType.values()){
            if(arity==1)
                queueUpPruningStepLevel(criterion, ctype, pruningSteps);
            else if(arity==2){
                for(RCClass.ClassType ctype2 : RCClass.ClassType.values()){//ORDER BETTER??
                    queueUpPruningStepLevel(criterion, ctype, ctype2, pruningSteps);
                }
            }
            else
                throw new RuntimeException("ERROR: unsupported arity "+arity);
        }
    }
    
    
    void queueUpPruningStepLevel(PruningCriteria criterion, RCClass.ClassType ctype, 
            ArrayList<GeomPruningStep> pruningSteps){//queue up singles pruning
        GeomPruningStep step;
        switch(criterion){
            case STERIC:
                step = new StericPruningStep(classesAtRes, curPruned, ctype, sp.shellResidues, stericCache);
                break;
            case FLEXRESCONTACT:
                step = new FlexResContactPruningStep(classesAtRes, curPruned, ctype, sp.shellResidues, mutableRes, stericCache);
                break;
            default:
                throw new RuntimeException("ERROR CRITERION NOT YET CODED");
        }
        pruningSteps.add(step);
    }
    
    
    void queueUpPruningStepLevel(PruningCriteria criterion, RCClass.ClassType ctype, 
            RCClass.ClassType ctype2, ArrayList<GeomPruningStep> pruningSteps){//pairs
        
        GeomPruningStep step;
        switch(criterion){
            case STERIC:
                step = new StericPruningStep(classesAtRes, curPruned, ctype, ctype2, sp.shellResidues, stericCache);
                break;
            case RAMA:
                step = new RamaPruningStep(classesAtRes, curPruned, ctype, ctype2, sp.shellResidues, stericCache, sp.confSpace);
                break;
            default:
                throw new RuntimeException("ERROR CRITERION NOT YET CODED");
        }
        pruningSteps.add(step);
    }
    
        
}
