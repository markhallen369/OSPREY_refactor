/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug.triage;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.bbfree.BBFreeBlock;
import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.plug.VoxelVDWListChecker;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 *
 * @author mhall44
 */
public abstract class GeomPruningStep {
    
    ArrayList<RCClassification> classesAtRes;
    PruningLookup prooningLookup;
    RCClass.ClassType classType;
    
    int tupSize;
    RCClass.ClassType classType2;//if needed
    
    ArrayList<Residue> shellResidues;
    
    CentralStericCache stericCache;

    public GeomPruningStep(ArrayList<RCClassification> classesAtRes, PruningLookup prooningLookup, 
            RCClass.ClassType classType, ArrayList<Residue> shellResidues, CentralStericCache stericCache) {
        this.shellResidues = shellResidues;
        this.classesAtRes = classesAtRes;
        this.prooningLookup = prooningLookup;
        this.stericCache = stericCache;
        this.classType = classType;
        tupSize = 1;
    }
    
    public GeomPruningStep(ArrayList<RCClassification> classesAtRes, PruningLookup prooningLookup, 
            RCClass.ClassType classType, RCClass.ClassType classType2, 
            ArrayList<Residue> shellResidues, CentralStericCache stericCache) {
        this.shellResidues = shellResidues;
        this.classesAtRes = classesAtRes;
        this.prooningLookup = prooningLookup;
        this.stericCache = stericCache;
        this.classType = classType;
        this.classType2 = classType2;
        tupSize = 2;
    }
    
    
    
    void prune(){
        Iterator<RCClassTuple> it;
        if(tupSize==1)
            it = new RCClassTupleIterator(classType, classesAtRes, prooningLookup);
        else
            it = new RCClassTupleIterator(classType, classType2, classesAtRes, prooningLookup);
        while(it.hasNext()){
            RCClassTuple cand = it.next();
            if(canPrune(cand)){
                prooningLookup.markAsPruned(cand);
            }
        }
    }
    
    
    
    abstract boolean canPrune(RCClassTuple tup);
    
    
    static boolean applyRCCenter(RC rc, Residue res, ResidueTypeDOF mutDOF){
        //return whether successfully set interval
        /*if(rc.template != res.template){
            mutDOF.switchToTemplate(rc.template);
        }*/
        //go by names for now, many RCs don't store templates (could force storage if need more?)
        if( ! rc.AAType.equalsIgnoreCase(res.template.name) ){
            mutDOF.mutateTo(rc.AAType);
        }
        
        //DEBUG!!!  calling code hackily for now
        VoxelVDWListChecker vvlc = new VoxelVDWListChecker();
        for(int dof=0; dof<rc.DOFs.size(); dof++){
            double lb = rc.DOFmin.get(dof);
            vvlc.addDOFInterval( rc.DOFs.get(dof), lb, rc.DOFmax.get(dof)-lb );
        }
        vvlc.addVoxLinConstr(rc.linConstr);
        
        HashMap<DegreeOfFreedom,Double> cen = vvlc.calcSpecialCenter();
        //add non-special parts of center
        for(int dof=0; dof<rc.DOFs.size(); dof++){
            DegreeOfFreedom curDOF = rc.DOFs.get(dof);
            if(!cen.containsKey(curDOF)){
                cen.put( curDOF, 0.5*(rc.DOFmin.get(dof)+rc.DOFmax.get(dof)) );
            }
        }
        
        //apply BBFreeDOFs in a block together for speed
        HashSet<BBFreeBlock> blocks = new HashSet<>();
        for(DegreeOfFreedom dof : cen.keySet()){
            if(dof instanceof BBFreeDOF){
                blocks.add((BBFreeBlock)dof.getBlock());
            }
        }
        for(BBFreeBlock block : blocks){
            ArrayList<BBFreeDOF> dofList = block.getDOFs();
            DoubleMatrix1D x = DoubleFactory1D.dense.make(dofList.size());
            for(int d=0; d<dofList.size(); d++)
                x.set(d,cen.get(dofList.get(d)));
            if(!block.setDOFs(x))
                return false;
        }
        //now apply the rest
        for(DegreeOfFreedom dof : cen.keySet()){
            if( ! (dof instanceof BBFreeDOF) )
                dof.apply(cen.get(dof));
        }
        
        return true;
    }
    
    //For debugging
    void writePDBFile(Molecule m, String name){
        PDBFileWriter.writePDBFile(m, name);
    }
    
}
