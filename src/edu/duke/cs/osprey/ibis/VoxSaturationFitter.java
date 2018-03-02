/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ibis;

/**
 *
 * This fitter first first takes a bunch of voxels "saturated" with samples (requires special sample set)
 * and fits a polynomial to describe (integrated) energy within each
 * Then fits combinatorial expansion to describe these polynomials 
 * (try/compare in both coeff space, and in samp space via TupIterFitter)
 * The main purpose of this fitter is to see if residuals are mostly from intra- or inter-voxel error.
 * 
 * 
 * @author mhall44
 */
public class VoxSaturationFitter implements SLEFitter {

    @Override
    public double[] fitCoeffs(double[] initCoeffs) {
        
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
