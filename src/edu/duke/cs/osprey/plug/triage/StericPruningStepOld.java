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
*/
/**
 *
 * @author mhall44
 */

/*
public class StericPruningStepOld extends GeomPruningStep {
    
    
    static double stericPruningOverlap = 1;//how much overlap at central conf is considered
    //to suffice for steric pruning (i.e., w/o even checking if clash can be escaped)

    StericPruningStepOld(ArrayList<RCClassification> classesAtRes, PruningLookup prooningLookup, 
            RCClass.ClassType classType, ArrayList<Residue> shellResidues) {
        super(classesAtRes, prooningLookup, classType, shellResidues);
    }
    
    StericPruningStepOld(ArrayList<RCClassification> classesAtRes, PruningLookup prooningLookup, 
            RCClass.ClassType classType, RCClass.ClassType classType2, ArrayList<Residue> shellResidues) {
        super(classesAtRes, prooningLookup, classType, classType2, shellResidues);
    }
    
    @Override
    boolean canPrune(RCClassTuple tup) {
        //will prune if the cores of the RC classes clash, indicating clash unavoidable for all 
        //RC bindings of class(es)
        
        if(tup.size()>2)
            throw new RuntimeException("ERROR: Steric pruning only makes sense for singles and pairs");
        
        for(int cl=0; cl<tup.size(); cl++){
            RCClass c = tup.get(cl);
            if(!applyRCCenter(c.core, c.res, c.mutDOF))
                return true;//invalid conf can def be pruned
        }
        
        Residue res1 = tup.classes.get(0).res;
        
        for(Atom at1 : res1.atoms){
            Iterable<Atom> otherAtoms;
            if(tup.size()==1)
                otherAtoms = new ISAtomsIterable(at1,res1,shellResidues);
            else
                otherAtoms = new Res2AtomsIterable(at1,res1,tup.classes.get(1).res);
            
            for(Atom at2 : otherAtoms){//atoms that may clash with at1
                double dist = VectorAlgebra.distance(at1.getCoords(), at2.getCoords());
                if( dist < getVDWRadius(at1)+getVDWRadius(at2)-stericPruningOverlap)
                    return true;
            }
        }
        
        //if we get here, no clashes
        return false;
    }
    
    

    
    
}
*/