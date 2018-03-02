/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof;

import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import java.io.Serializable;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author mhall44
 */
public abstract class DegreeOfFreedom implements Serializable {
    //This class represents a conformational degree of freedom
    //It is used to abstract out flexibility: for purposes of all the conformational search machinery,
    //the conformation can be represented as a vector of degree-of-freedom values
    //then implementations of this class determine what that vector means conformationally
    //by applying appropriate changes to molec. 
    
    //Let's say molecule, once created, can only be changed by DegreeOfFreedom.apply!
    
    private static final long serialVersionUID = 3348198591978172994L;
    
    double curVal;
    
    public double confDOFNum = -1;//number among confDOFs in a confSpace
    
    public abstract void apply(double paramVal);//apply the given parameter value
    //(some degrees of freedom may have convenience methods to call this, e.g. mutation called by aa type)
    
    public double getCurVal() { return curVal; }
    
    
    //If this DegreeOfFreedom moves only a single residue, return that residue
    //Otherwise return null
    //Used in setting up partial energy functions (see MultiTermEnergyFunction)
    public Residue getResidue() { return null; }
    
    // enables parallel molecule manipulation without data races
    // these two methods are only implemented for perturbations that aren't part of a block
    // (DOFBlock.copyForNewMolecule handle these operations in that case)
    public DegreeOfFreedom copy() {
        throw new UnsupportedOperationException("unsupported by " + getClass().getName());
    }
    public void setMolecule(Molecule val) {
        throw new UnsupportedOperationException("unsupported by " + getClass().getName());
    }
    
    public abstract DOFBlock getBlock();//return the DOF block for a DOF (return null if none)
    
    
    public List<Residue> listAffectedResidues(){
        Residue res = getResidue();
        if(res!=null)
            return Arrays.asList(res);
        
        DOFBlock block = getBlock();
        if(block==null)
            throw new RuntimeException("ERROR: Can't list affected residues if neither single residue nor block specified");
        else
            return block.listResidues();
    }
    
    
    public abstract String getName();//make a name for this DOF
    //should suffice to uniquely identify the DOF among DOFs that appear in a given system
    //but should be the same for equivalent DOFs in different copies of a system
}
