/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 *
 * SVD test
 * 
 * @author mhall44
 */
public class Carpenter {
    
    public static void main(String[] args){
        DoubleMatrix2D M = (DoubleMatrix2D) ObjectIO.readObject("M.dat", false);
        
        for(int i=0; i<M.rows(); i++)
            M.set(i, i, M.get(i,i)+1e-8);
        
        
        SingularValueDecomposition svd = new SingularValueDecomposition(M);
        System.out.println(svd.getS());
        System.out.println(svd.getU());
        System.out.println(svd.getV());
    }
    
}
