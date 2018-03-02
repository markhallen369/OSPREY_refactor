/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug.triage;

import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.tools.RigidBodyMotion;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author mhall44
 */
class BBVoxSphereModel {
    
    
    HashMap<String,double[]> bbAtomCoords = new HashMap<>();//map name to coords
    RigidBodyMotion sidechainTransformation;

    BBVoxSphereModel(ResSphereModel firstRCModel, ResSphereModel altBBModel, ResidueTemplate template) {
        //want the backbone described in altBBModel
        //firstRCModel provided so we can see how the sidechain transforms to get to this state
        ArrayList<Atom> templAtoms = template.templateRes.atoms;
        if(firstRCModel.numAtoms!=templAtoms.size() || altBBModel.numAtoms!=templAtoms.size())
            throw new RuntimeException("ERROR: atom count mismatch!");//both models should match template restype
        
        double origNCAC[][] = new double[3][3];
        double altNCAC[][] = new double[3][3];
        
        for(int at=0; at<templAtoms.size(); at++){
            String atomName = templAtoms.get(at).name;
            if(HardCodedResidueInfo.possibleBBAtomsLookup.contains(atomName)){//backbone atom
                double coords[] = new double[3];
                System.arraycopy(altBBModel.coords, 3*at, coords, 0, 3);
                bbAtomCoords.put(atomName, coords);
            
                //sidechain is placed based on backbone N, CA, C (e.g. by idealizer)
                //CATS motions move the N-CA-C triangle as a rigid body (barring truncation error)
                //so we'll get the rigid-body transformation from N, CA and C
                if(atomName.equalsIgnoreCase("CA")){
                    System.arraycopy(coords, 0, altNCAC[0], 0, 3);
                    System.arraycopy(firstRCModel.coords, 3*at, origNCAC[0], 0, 3);
                }
                else if(atomName.equalsIgnoreCase("N")){
                    System.arraycopy(coords, 0, altNCAC[1], 0, 3);
                    System.arraycopy(firstRCModel.coords, 3*at, origNCAC[1], 0, 3);
                }
                else if(atomName.equalsIgnoreCase("C")){
                    System.arraycopy(coords, 0, altNCAC[2], 0, 3);
                    System.arraycopy(firstRCModel.coords, 3*at, origNCAC[2], 0, 3);
                }
            }
        }
        
        sidechainTransformation = new RigidBodyMotion(origNCAC, altNCAC);
    }
    
    
    ResSphereModel transform(ResSphereModel origBBModel, ResidueTemplate templ){
        //given a ResSphereModel in the orig bb conf, corresponding to given template,
        //create a transformed version based on this bb vox center
        
        double newCoords[] = origBBModel.coords.clone();
        double newRadii[] = origBBModel.radii.clone();
        int numAtoms = origBBModel.numAtoms;
        
        for(int at=0; at<numAtoms; at++){//copy correct backbone coords in
            String atName = templ.templateRes.atoms.get(at).name;

            if(bbAtomCoords.containsKey(atName)){
                System.arraycopy(bbAtomCoords.get(atName), 0, newCoords, 3*at, 3);
            }
            else
                sidechainTransformation.transform(newCoords, at);
        }
        
        return new ResSphereModel(newCoords, newRadii);
    }
}
