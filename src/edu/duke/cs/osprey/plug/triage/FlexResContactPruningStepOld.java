/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*

package edu.duke.cs.osprey.plug.triage;

import edu.duke.cs.osprey.plug.ISAtomsIterable;
import edu.duke.cs.osprey.plug.Res2AtomsIterable;
import static edu.duke.cs.osprey.plug.VoxelVDWDistExplorer.getVDWRadius;
import static edu.duke.cs.osprey.plug.triage.GeomPruningStep.applyRCCenter;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.util.ArrayList;
import java.util.HashSet;
*/
/**
 *
 * not only does this version not use steric cache,
 * it's also doing some weird hybrid of steric & frc that doesn't make sense...
 * 
 * @author mhall44
 */


/*
public class FlexResContactPruningStepOld extends GeomPruningStep {
    
    
    static double contactBuffer = 1.5;//Mutable residue RCs are required to be
    //within this distance of non-mutable (i.e. target) residues, measured between VDW spheres, to pass this
    //pruning step.  
    //If at this distance they could probably get into contact by within-voxel minimization

    HashSet<Residue> mutableRes = null;
    
    FlexResContactPruningStep(ArrayList<RCClassification> classesAtRes, PruningLookup prooningLookup, 
            RCClass.ClassType classType, ArrayList<Residue> shellResidues,
            HashSet<Residue> mutableRes) {
        super(classesAtRes, prooningLookup, classType, shellResidues);
        this.mutableRes = mutableRes;
    }
    
    FlexResContactPruningStep(ArrayList<RCClassification> classesAtRes, PruningLookup prooningLookup, 
            RCClass.ClassType classType, RCClass.ClassType classType2, 
            ArrayList<Residue> shellResidues, HashSet<Residue> mutableRes) {
        super(classesAtRes, prooningLookup, classType, classType2, shellResidues);
        this.mutableRes = mutableRes;
    }
    
    @Override
    boolean canPrune(RCClassTuple tup) {
        //will prune if the first residue is mutable
        //but can't touch any non-mutable residues
        //yeah I know this is kind of a weird criterion
        //but it's a simple example of a packing criterion that might sometimes be useful
        
        if(tup.size()>2)
            throw new RuntimeException("ERROR: Flex-res contact pruning only makes sense for singles and pairs");
        
        Residue res1 = tup.classes.get(0).res;
        if(!mutableRes.contains(res1))
            return false;

        
        for(int cl=0; cl<tup.size(); cl++){
            RCClass c = tup.get(cl);
            if(!applyRCCenter(c.core, c.res, c.mutDOF))
                return true;
        }
        
        
        //OK now we'll look for a possible contact.  Prune if can't find one!
        //first try pivot atom
        Atom pivot = tup.get(0).pivotAtom;
        if(pivot!=null){
            double bubbleRad = tup.get(0).bubbleRad;
            for(Atom at2 : otherAtomsIterable(pivot,res1,tup)){//atoms that may contact at1
                if(!mutableRes.contains(at2.res)){
                    
                    //DEBUG!!!!!forcing contact with cxcl12
                    if(Integer.valueOf(at2.res.getPDBResNumber())<100){
                        
                        double dist = VectorAlgebra.distance(pivot.getCoords(), at2.getCoords());
                        if( dist < bubbleRad+getVDWRadius(at2)+contactBuffer )
                            return false;
                    }
                }
            }
        }
        
        for(Atom at1 : res1.atoms){
            for(Atom at2 : otherAtomsIterable(at1,res1,tup)){//atoms that may contact at1
                if(!mutableRes.contains(at2.res)){
                    
                    //DEBUG!!!!!forcing contact with cxcl12
                    if(Integer.valueOf(at2.res.getPDBResNumber())<100){
                    
                        double dist = VectorAlgebra.distance(at1.getCoords(), at2.getCoords());
                        if( dist < getVDWRadius(at1)+getVDWRadius(at2)+contactBuffer )
                            return false;
                    }
                }
            }
        }
        
        //if we get here, no contacts possible
        return true;
    }
    
    
    private Iterable<Atom> otherAtomsIterable(Atom at1, Residue res1, RCClassTuple tup){
        //other atoms that may contact atom 1
        if(tup.size()==1)
            return new ISAtomsIterable(at1,res1,shellResidues);
        else
            return new Res2AtomsIterable(at1,res1,tup.classes.get(1).res);
    }
    

    
}
*/