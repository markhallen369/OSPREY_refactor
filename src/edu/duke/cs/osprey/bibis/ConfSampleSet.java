/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.bibis;

import edu.duke.cs.osprey.plug.PolytopeMatrix;
import java.io.Serializable;
import java.util.Iterator;

/**
 *
 * Like TESampleSet but with contDOFs too
 * 
 * Could be represented as a list of samples, but
 * some versions may have samples that are only different at 1 or 2 pos, etc.
 * 
 * @author mhall44
 */
public interface ConfSampleSet extends Serializable {
    
    public EnergiedConfSampleSet integrateEnergies(IntegrableDOF dof, EnergyModel E);
    
    public EnergiedConfSampleSet evalEnergies(EnergyModel E);//e.g. to replace AMBER w/ AMOEBA energies
    //but keep slightly augmented AMBER features
    
    public void checkPLUGViolations(PolytopeMatrix plugMat);
    
    public int getNumSamples();
    
}
