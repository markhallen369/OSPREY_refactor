/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.bibis;

import java.io.Serializable;

/**
 *
 * A DOF wrt which we can evaluate a Boltzmann integral
 * 
 * @author mhall44
 */
public interface IntegrableDOF extends Serializable {
    
    public static final double RT = 1.9891/1000.0 * 298.15;//RT in kcal/mol..see PoissonBoltzmannEnergy
    
    abstract double boltzmannIntegratedEnergy(EnergyModel baseE, ConfSample samp);
    //integrate baseE over this DOF; fixed DOF values as given by samp
    
    abstract FeatureSet featureSetForInteg(FeatureSet baseFeat);
    //given a feature set for baseE, come up with a good post-integation feature set
    //note: this could take some validation and adjustemnt...if fitter breaks to terms
    //this could give us valuable information.  May only be making a candidate here, 
    //fitter may have to expand it
    
    
    abstract String description();
    
         
        
}
