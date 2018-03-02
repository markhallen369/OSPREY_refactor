/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.bibis;

/**
 *
 * A ConfSampleSet with energies (usually partially integrated) for all the samples
 * Can be generated from a regular ConfSampleSet by calculating energies for the samples
 * 
 * @author mhall44
 */
public interface EnergiedConfSampleSet extends ConfSampleSet {
    
    void checkError(EnergyModel f);
        //Check how well f models the energies this sample set
        //print output
    
}
