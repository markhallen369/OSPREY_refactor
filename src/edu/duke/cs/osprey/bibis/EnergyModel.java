/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.bibis;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * Model of the energy as a function of some degrees of freedom
 * 
 * @author mhall44
 */
public interface EnergyModel extends Serializable {
    
    double evalEnergy(ConfSample samp);
    EnergyModel integrateModel(IntegrableDOF dof, EnergiedConfSampleSet samples);//make a version of the model integrated wrt dof
    
    Iterable<DegreeOfFreedom> getContDOFs();//what continuous DOFs does this energy model depend on?
    int getNumParams();
    
    double evalEnergy(ConfSample samp, MoleculeModifierAndScorer mms);//use the provided mms to evaluate the energy, if applicable
    MoleculeModifierAndScorer makeMMS(ConfSample samp);//build the mms.  Null if no mms used in evaluating this EnergyModel
}
