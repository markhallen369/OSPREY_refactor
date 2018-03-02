/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug.triage;

import edu.duke.cs.osprey.structure.AtomNeighbors;
import edu.duke.cs.osprey.structure.AtomNeighbors.NEIGHBORTYPE;
import edu.duke.cs.osprey.structure.ProbeAtomNeighbors;
import edu.duke.cs.osprey.structure.Residue;
import java.util.BitSet;

/**
 *
 * Keeps track of the atom pairs in a pair of a residues
 * that can be close without being considered a clash
 * this is specific to particular AA types for the residues
 * 
 * @author mhall44
 */
public class ResPairClashExemptions {
 
    int firstResNumAtoms;
    BitSet exemptPairs = new BitSet();

    public ResPairClashExemptions(Residue res, Residue nres) {
        //base exemption on current residue types of res, nres
        firstResNumAtoms = res.atoms.size();
        for(int at=0; at<firstResNumAtoms; at++){
            ProbeAtomNeighbors neigh = new ProbeAtomNeighbors(res.atoms.get(at));
            for(int at2=0; at2<nres.atoms.size(); at2++){
                if(neigh.classifyAtom(nres.atoms.get(at2)) != NEIGHBORTYPE.NONBONDED)
                    exemptPairs.set(index(at,at2));
            }
        }
    }
    
    
    public boolean contains(int at1, int at2){
        //is the atom pair (at1,at2) (indices in respective residues) exempt?
        return exemptPairs.get(index(at1,at2));
    }
    
    private int index(int at1, int at2){
        //combine in single index
        return at1+at2*firstResNumAtoms;
    }
}
