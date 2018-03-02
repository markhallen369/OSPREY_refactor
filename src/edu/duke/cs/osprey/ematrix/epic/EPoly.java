package edu.duke.cs.osprey.ematrix.epic;

/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.1 beta
	Copyright (C) 2001-2012 Bruce Donald Lab, Duke University

	OSPREY is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as
	published by the Free Software Foundation, either version 3 of
	the License, or (at your option) any later version.

	OSPREY is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, see:
	      <http://www.gnu.org/licenses/>.

	There are additional restrictions imposed on the use and distribution
	of this open-source code, including: (A) this header must be included
	in any modification or extension of the code; (B) you are required to
	cite our papers in any publications that use this code. The citation
	for the various different modules of our software, together with a
	complete list of requirements and restrictions are found in the
	document license.pdf enclosed with this distribution.

	Contact Info:
			Bruce Donald
			Duke University
			Department of Computer Science
			Levine Science Research Center (LSRC)
			Durham
			NC 27708-0129
			USA
			e-mail:   www.cs.duke.edu/brd/

	<signature of Bruce Donald>, Mar 1, 2012
	Bruce Donald, Professor of Computer Science
 */

///////////////////////////////////////////////////////////////////////////////////////////////
//	EPoly.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.plug.RCTuplePolytope;
import java.io.FileInputStream;
import java.io.Serializable;
import java.util.ArrayList;


public class EPoly implements Serializable {
    //energy term represented by a polynomial fit
    //possibly augmented by an additional "sparse" energy term
    
    
    //The continuously varying degrees of freedom the terms acts on
    int numDOFs;
    //bounds on the degrees of freedom
    DoubleMatrix1D DOFmax, DOFmin;
    ArrayList<DegreeOfFreedom> DOFs;//the actual degrees of freedom

    

    double coeffs[];//coefficients for series (expanded in original, centered coordinates)
    int order;//order of polynomial (can be 3 or 4)

    
    
    DoubleMatrix1D center;//the center (typically, minimum point) for relative coordinates for this term
    double minE;//energy at center
    //will be included in evaluate if includeMinE is called, else center will evalute to 0
    
    
    String fitDescription = "N/A";//description of fit used to get this EPoly, if applicable

    
    
    //Sparse VDW terms
    public SAPE sapeTerm = null;
    public double baseSAPE = 0;//value of SAPE terms at center (SAPE will be evaluated relative to this)

    RCTuplePolytope tope = null;//If not null, the EPoly is fit to be good within this polytope
    //which indicates the non-clashing region (and possibly other restrictions we desire)

    public EPoly(int numDOFs, ArrayList<DegreeOfFreedom> DOFs, DoubleMatrix1D DOFmax, DoubleMatrix1D DOFmin, 
            DoubleMatrix1D center, double minE, double[] coeffs, int order, RCTuplePolytope tope ) {
        
        this.numDOFs = numDOFs;
        this.DOFs = DOFs;
        this.DOFmax = DOFmax;
        this.DOFmin = DOFmin;
        this.center = center;
        this.minE = minE;
        this.coeffs = coeffs;
        this.order = order;
        this.tope = tope;
    }

    
    
    public double evaluate(DoubleMatrix1D x, boolean includeMinE, boolean useSharedMolec) {
        //evaluate this EPoly as a function of internal coordinates x
        //Faster if we can use a shared molecule (this could be possible if this EPoly
        //is part of a MolecEObjFunction)
        
        double serVal = evalPolynomial(x);
        
        if(includeMinE)
            serVal += minE;
                 
        if(sapeTerm!=null){//Need to include the SAPE term
            
            if(useSharedMolec)//shared molecule assumed to be in the right conformation already
                return serVal + sapeTerm.getEnergySharedMolec() - baseSAPE;
            else //Use the molecule stored in the SAPE object
                return serVal + sapeTerm.getEnergyStandalone(x) - baseSAPE;
        }
        else 
            return serVal;
    }
    
    public double evalPolynomial(DoubleMatrix1D x){
        return evalSeries(toRelCoords(x));
    }
    
    
    double evalSeries(DoubleMatrix1D z){
        //evaluate the actual series
        //(function of relative coordinates)
        return SeriesFitter.evalSeries(coeffs, z, numDOFs, false, order);
    }
    
    
    /*
     * These functions might be useful for EPIC fitting of non-pairwise energies?
     * 
     * void recenterAtMinimum(DegreeOfFreedom[] allDOFs){
        //recenter at the minimum
        
        double constraints[][] = new double[2][];
        constraints[0] = DOFmin.toArray();
        constraints[1] = DOFmax.toArray();
        
        DoubleMatrix1D DOFVals = DoubleFactory1D.dense.make(numDOFs);
        CETEnergyFunction ef = new CETEnergyFunction( DOFVals, DOFNums, new ContETerm[] {this}, true );
        CETObjFunction of = new CETObjFunction(DOFNums,constraints,allDOFs,ef,null);
        CCDMinimizer ccdMin = new CCDMinimizer(of,false);
        DoubleMatrix1D newCenter = ccdMin.minimize();
        
        recenter(newCenter);
    }
    
    
    void recenter(DoubleMatrix1D x){
        //move center to x
        //adjust coefficients (including minE=constant coeff.) accordingly
        DoubleMatrix1D dcen = x.copy().assign(center,Functions.minus);//change in center
        //if the current polynomial is f(x), the new polynomial must be f(x+dcen)
        MultivariatePoly currentPoly = new MultivariatePoly(coeffs,numDOFs,order,false,order,null);
        currentPoly.setConstant(minE);
        MultivariatePoly newPoly = currentPoly.shift(dcen);
        
        minE = newPoly.getConstant();
        center = x;
        newPoly.removeConstant();
        coeffs = newPoly.toSeriesDoubleArray(false,order,null);
    }
    
    
    void fixBaseSVE(){
        //take a polynomial in the form (series value) + (SVE value) + (minE if appropriate)
        //and convert it to the form (series value) + (SVE value) - baseSVE + (minE if appropriate)
        //so that minE is actually the energy at center
        sveOF.setDOFs(center);
        baseSVE = sapeTerm.getEnergy();
        minE += baseSVE;
    }*/
    
    public DoubleMatrix1D polynomialGradient(DoubleMatrix1D x){
        DoubleMatrix1D z = toRelCoords(x);
        
        if(this instanceof EPolyPC)
            throw new RuntimeException("ERROR: gradient for EPolyPC not currently supported");
        
        return SeriesFitter.evalSeriesGradient(coeffs, z, numDOFs, false, order, order, null);
    }
    
    public DoubleMatrix1D gradient(DoubleMatrix1D x/*, boolean useSharedMolec*/) {
        //evaluate this EPoly gradient as a function of internal coordinates x  
        if(sapeTerm!=null)
            throw new RuntimeException("ERROR: SVE gradient not currently supported");
        else 
            return polynomialGradient(x);
    }
    
    
    
    public DoubleMatrix2D polynomialHessian(DoubleMatrix1D x){
        DoubleMatrix1D z = toRelCoords(x);
        
        if(this instanceof EPolyPC)
            throw new RuntimeException("ERROR: Hessian for EPolyPC not currently supported");
        
        return SeriesFitter.evalSeriesHessian(coeffs, z, numDOFs, false, order, order, null);
    }
    
    public DoubleMatrix2D hessian(DoubleMatrix1D x) {
        if(sapeTerm!=null)
            throw new RuntimeException("ERROR: SVE Hessian not currently supported");
        else 
            return polynomialHessian(x);
    }
    
    
    
    
    
    DoubleMatrix1D toRelCoords(DoubleMatrix1D x){
        //convert coordinates to be relative to center
        DoubleMatrix1D ans = x.copy();
        ans.assign(center,Functions.minus);
        return ans;
    }
    
    
    //reverse conversion
    DoubleMatrix1D toAbsCoords(DoubleMatrix1D z){
        //convert coordinates to be relative to center
        DoubleMatrix1D ans = z.copy();
        ans.assign(center,Functions.plus);
        return ans;
    }
    
    
    public boolean valuesInRange(DoubleMatrix1D x){
        //Check if x is in range for this voxel
        for(int dof=0; dof<numDOFs; dof++){
            if(x.get(dof)>DOFmax.get(dof) || x.get(dof)<DOFmin.get(dof))
                return false;
        }
        return true;
    }

    public void setMinE(double minE) {
        this.minE = minE;
    }

    public double getMinE() {
        return minE;
    }
    
    
    
    
    public void deleteMOFStandalone(){
        //If using SAPE, delete the SAPE standalone energy evaluator to save memory
        if(sapeTerm!=null)
            sapeTerm.mofStandalone = null;
    }

    public double[] getCoeffs() {
        return coeffs;
    }

    public int getOrder() {
        return order;
    }

    public ArrayList<DegreeOfFreedom> getDOFs() {
        return DOFs;
    }

    public DoubleMatrix1D getCenter() {
        return center;
    }
    
    
    
    
    

    
    
    
    
    
}
