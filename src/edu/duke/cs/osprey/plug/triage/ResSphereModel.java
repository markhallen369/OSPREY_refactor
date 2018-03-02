/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug.triage;

import static edu.duke.cs.osprey.plug.VoxelVDWDistExplorer.getVDWRadius;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.VectorAlgebra;

/**
 *
 * @author mhall44
 */
class ResSphereModel {
    
    //represent a residue for steric-check purposes as a series of spheres (representing atoms)
    double[] coords;//coordinates of the atom centers
    double[] radii;//corresponding VDW radii
    int numAtoms;
    
    //sometimes one of the atoms can be modeled as bigger if trying to upper bound the steric extent
    //of the res
    int pivotAtom = -1;
    double bubbleRad = 0;
    
    //for quick screening, make a bounding sphere about the individual atom spheres
    double[] boundingSphereCenter;
    double boundingSphereRadius;
    
            //flexrescontact too!  and can now do quickclash, even fullclash fairly quickly
        //ways to handle linearization error if needed: error buffer; quadratic terms
        //I think this method will need to be replaced
    
    
    public ResSphereModel(double[] coords, double[] radii){
        this.coords = coords;
        this.radii = radii;
        numAtoms = radii.length;
        computeBoundingSphere();
    }

    ResSphereModel(Residue res) {
        coords = res.coords.clone();
        numAtoms = res.atoms.size();
        radii = new double[numAtoms];
        for(int at=0; at<numAtoms; at++)
            radii[at] = getVDWRadius(res.atoms.get(at));
        computeBoundingSphere();
    }
    
    ResSphereModel(Residue res, Atom pivotAtom, double bubbleRad) {
        //make a model where pivotAtom is inflated to bubbleRad (for outer-bounding purposes)
        coords = res.coords.clone();
        numAtoms = res.atoms.size();
        radii = new double[numAtoms];
        for(int at=0; at<numAtoms; at++){
            radii[at] = getVDWRadius(res.atoms.get(at));
            if(pivotAtom!=null){
                if(res.atoms.get(at).name.equalsIgnoreCase(pivotAtom.name))
                    radii[at] = Math.max(radii[at], bubbleRad);
            }
        }
        computeBoundingSphere();
    }
    
    private void computeBoundingSphere(){
        //compute bounding sphere from invidual atoms
        boundingSphereCenter = new double[3];//average atom coords for this
        for(int at=0; at<numAtoms; at++){
            for(int dim=0; dim<3; dim++){
                boundingSphereCenter[dim] += coords[3*at+dim]/numAtoms;
            }
        }
        
        boundingSphereRadius = 0;
        for(int at=0; at<numAtoms; at++){
            double dist = VectorAlgebra.distance(boundingSphereCenter, 0, coords, at);
            boundingSphereRadius = Math.max(boundingSphereRadius, dist+radii[at]);
        }
    }
    
    
    boolean clashesWith(ResSphereModel model2, double buffer, ResPairClashExemptions exemptions){
        //if spheres overlap between this and model2 by more than buffer,
        //and the overlapping spheres are not exempted, 
        //then return true
        
        double centerDist = VectorAlgebra.distance(boundingSphereCenter, model2.boundingSphereCenter);
        if( centerDist > boundingSphereRadius+model2.boundingSphereRadius-buffer ){
            //bounding spheres don't intersect (at least not within buffer),
            //so no way to have any clashes
            return false;
        }
        
        for(int at=0; at<numAtoms; at++){
            for(int at2=0; at2<model2.numAtoms; at2++){
                double dist = VectorAlgebra.distance(coords, at, model2.coords, at2);
                if( dist < radii[at]+model2.radii[at2]-buffer ){//atoms close enough to clash
                    if(exemptions!=null){//there are exemptions, see if (at,at2) is one
                        if(exemptions.contains(at,at2))
                            continue;
                    }
                    
                    return true;
                }
            }
        }
        
        return false;
    }
    
}
