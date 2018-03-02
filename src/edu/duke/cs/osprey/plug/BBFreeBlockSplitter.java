/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;

/**
 *
 * 
 * Makes a mini-block, intended for restricted motion at some DOFs,
 * with a different center
 * intended for when the big block is too big for accurate Taylor series, but 
 * is still free of singularities (allowing its free DOFs to be used throughout a larger area)
 * 
 * 
 * @author mhall44
 */
class BBFreeBlockSplitter {

    
    void makeDOFMapping() {
        //DEBUG!!!!  not making this second block for now
    }
    //the new DOFs should hold a copy of their "parent" so they can update curFreeDOFVals
    //This suffices to keep the old block updated,
    //because each call to setDOFs will set the conf correctly to the specified free DOF vals (which may be based on curFreeDOFVals)
    //regardless of initial state
    
    //new DOFs are equivalent to old ones just intended for use differently
    //we can even first try w/o them
    //because can probably get the same effect by iteration.  

    DegreeOfFreedom mapDOF(DegreeOfFreedom dof) {
        return dof;
    }
    
}
