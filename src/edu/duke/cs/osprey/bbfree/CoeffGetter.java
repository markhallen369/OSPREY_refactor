/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.bbfree;

import cern.colt.matrix.DoubleMatrix1D;

/**
 *
 * @author mhall44
 */
public class CoeffGetter {

    public static DoubleMatrix1D getCoeffs(BBFreeDOF dof){
        return dof.coeffs;
    }
}

