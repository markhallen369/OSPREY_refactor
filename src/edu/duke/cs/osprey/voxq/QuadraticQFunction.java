/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.voxq;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import org.apache.commons.math3.special.Erf;

/**
 * Q-function for Metropolis sampling in one dimension within a voxel
 * the energy surface is likely to be roughly quadratic in each dimension
 * (at least around the current point)
 * 
 * @author mhall44
 */
public class QuadraticQFunction {
    //Q = exp(ax^2+bx+c), x=DOF value will be between xLo, xHi
    double a, b, c, xLo, xHi;
    boolean useLinearPrior;//two ways to draw: constrained linear or unconstrained quadratic prior

    /*public QuadraticQFunction(MoleculeModifierAndScorer mms, int dof, double origVal) {
        //specify voxel/obj fcn, DOF being sampled (indexed in mms), and starting value of dof
        
        //start simple: convex quadratic thru origVal and edges
        //if too much autocorr or rejection will refine
        DoubleMatrix1D constr[] = mms.getConstraints();
        xLo = constr[0].get(dof);
        xHi = constr[1].get(dof);
        //constant factor, i.e. energy offset, doesn't matter
        double xLoE = mms.getValForDOF(dof, xLo);
        double xHiE = mms.getValForDOF(dof, xHi);
        
        if(xHi<xLo+1e-14)
            throw new RuntimeException("ERROR: Trying to sample a rigid DOF!");
        
        //origVal directly at the edge is not useful in building a quadratic,
        //so move our "3rd pt" 10% away from the edge in that case
        origVal = Math.max(origVal,0.9*xLo+0.1*xHi);
        origVal = Math.min(origVal,0.1*xLo+0.9*xHi);
        double origValE = mms.getValForDOF(dof, origVal);
        
        //OK let's check if origVal is below the straight line between xLo & xHi
        //if so we'll do the quadratic, else just use the straight line
        double slope = (xHiE-xLoE)/(xHi-xLo);
        
        if( origValE-xLoE < slope*(origVal-xLo) ){
            a = ( (origValE-xLoE)/(origVal-xLo) - (xHiE-xLoE)/(xHi-xLo) ) / (origVal-xHi);
            b = ( xHiE-xLoE - a*(xHi*xHi-xLo*xLo) )/ (xHi-xLo);
            c = xHiE - xHi*(b+a*xHi);
        }
        else {
            a = 0;
            b = slope;
            c = xHiE - xHi*slope;
        }
        
        //do E to -E/RT conversion
        a /= -IntraVoxelSampler.RT;
        b /= -IntraVoxelSampler.RT;
        c /= -IntraVoxelSampler.RT;
        
        //Finally, normalize the distribution
        normalizeDistr();
        
        if( a!=0 && (!erfcInvNumericsOK()) ){
            a = 0;
            b = -slope/IntraVoxelSampler.RT;
            c = -(xHiE - xHi*slope)/IntraVoxelSampler.RT;
            normalizeDistr();
        }
        
        
        //DEBUG!!
        for(double q : new double[]{a,b,c})
            if(Double.isInfinite(q) || Double.isNaN(q))
                System.out.println("What's wrong with you bluggles?");
    }*/



    public QuadraticQFunction(MoleculeModifierAndScorer mms, int dof, double origVal) {
        //specify voxel/obj fcn, DOF being sampled (indexed in mms), and starting value of dof
        
        //start simple: convex quadratic thru origVal and edges
        //if too much autocorr or rejection will refine
        DoubleMatrix1D constr[] = mms.getConstraints();
        xLo = constr[0].get(dof);
        xHi = constr[1].get(dof);
        if(xHi<xLo+1e-14)
            throw new RuntimeException("ERROR: Trying to sample a rigid DOF!");
        
        //constant factor, i.e. energy offset, doesn't matter
        double origValE = mms.getValForDOF(dof, origVal);
        double xLoE = mms.getValForDOF(dof, xLo);
        double xHiE = mms.getValForDOF(dof, xHi);
        
        //Choose x1 and x3 to roughly bracket the region from which it is useful to sample
        //(they'll be voxel edges unless need to move in to avoid clashes)
        double x1 = chooseReasonableOuterPt(origVal, origValE, xLo, xLoE, mms, dof);
        double x3 = chooseReasonableOuterPt(origVal, origValE, xHi, xHiE, mms, dof);
        
        //choose a point in the middle to build a quadratic (origVal unless too close to one of the edges,
        //in which case move 10% away from the edge)
        double x2 = Math.max(origVal,0.9*x1+0.1*x3);
        x2 = Math.min(x2,0.1*x1+0.9*x3);
        
        double E1 = getEnergyIfNeeded(x1, mms, dof, xLo, xLoE);
        double E3 = getEnergyIfNeeded(x3, mms, dof, xHi, xHiE);
        double E2 = getEnergyIfNeeded(x2, mms, dof, origVal, origValE);
        
        
        double slope = (E3-E1)/(x3-x1);
        
        if( E2-E1 < slope*(x2-x1) ){//middle point lies below line between outer points
            //build quadratic
            boolean success = setupQuadratic(x1, x2, x3, E1, E2, E3);
            if( ( a!=0 && (!erfcInvNumericsOK()) ) || (!success) ){
                //quadratic bad numerically.  First condition is for erfc numerics, second for normalization constant
                setupLinear(x1, x3, E1, E3);
            }
        }
        else
            setupLinear(x1, x3, E1, E3);
        
        for(double q : new double[]{a,b,c})
            if(Double.isInfinite(q) || Double.isNaN(q))
                throw new RuntimeException("ERROR: Infinite or NaN coefficient in sampling prior");
    }
    
    private double getEnergyIfNeeded(double x, MoleculeModifierAndScorer mms, int dof, double xDone,
            double Edone){
        //We have energy for xDone; compute x if it's different
        if(x!=xDone)
            return mms.getValForDOF(dof, x);
        else
            return Edone;
    }
    
    
    private static double chooseReasonableOuterPt(double origVal, double origValE, double outerPt, double outerPtE, 
        MoleculeModifierAndScorer mms, int dof){
        //Basically want E(outerPt) to be not too much higher than origValE
        
        double top = 20;//try to keep E(outerPt) below this
        //For a realistic well, energies this high should have negligible probability
        double bottom = 10;//if outerPt starts above top, let's keep it below bottom to make sure
        //high-energy region represented OK
       
        boolean neededAdjustment = false;//did we need to move outerPt in
        
        while(outerPtE > origValE + top){//moving inward
            neededAdjustment = true;
            outerPt = 0.5*(origVal+outerPt);
            outerPtE = mms.getValForDOF(dof, outerPt);
        }
        
        if(neededAdjustment){//move back out to get above bottom (if needed)
            double prevOuterPt = 2*outerPt - origVal;//this had energy above origValE+top
            while(outerPtE <= origValE + bottom){//OK this might overshoot top a little but no big deal
                outerPt = 0.8*outerPt + 0.2*prevOuterPt;
                outerPtE = mms.getValForDOF(dof, outerPt);
            }
        }
        
        return outerPt;
    }
    
    
    boolean setupQuadratic(double x1, double x2, double x3, double E1, double E2, double E3){
        //set up the prior by putting a quadratic through the specified points
        a = ( (E2-E1)/(x2-x1) - (E3-E1)/(x3-x1) ) / (x2-x3);
        b = ( E3-E1 - a*(x3*x3-x1*x1) )/ (x3-x1);
        c = E3 - x3*(b+a*x3);
                
        //do E to -E/RT conversion
        a /= -IntraVoxelSampler.RT;
        b /= -IntraVoxelSampler.RT;
        c /= -IntraVoxelSampler.RT;
        
        //Finally, normalize the distribution
        return normalizeDistr();
    }
    
    void setupLinear(double x1, double x3, double E1, double E3){
        a = 0;
        double slope = (E3-E1)/(x3-x1);
        b = -slope/IntraVoxelSampler.RT;
        c = -(E3 - x3*slope)/IntraVoxelSampler.RT;
        normalizeDistr();
    }
    
    
    double drawDOFValue(){
        return cumulDistrInv(Math.random());
    }
    
    double evalQ(double x){
        return Math.exp(c + x*(b+a*x));
    }
    
    
    private double cumulDistr(double x){
        //integral of exp(ax^2+bx+c) from xLo to x
        //this should be fairly resistant to under/overflow once the distribution is normalized
        if(Math.abs(a)<1e-14){//linear case
            return ( Math.exp(b*x+c) - Math.exp(b*xLo+c) )/b;
        }
        else {
            double C = 0.5*Math.exp(c-0.25*b*b/a)*Math.sqrt(-Math.PI/a);
            return C * ( Erf.erf( cdErfArg(xLo), cdErfArg(x) ) );
        }
    }
    
    private boolean normalizeDistr(){
        //want integral of exp(ax^2+bx+c) from xLo to x to equal 1
        //adjust c accordingly.  Return false if unsuccessful
        if(Math.abs(a)<1e-14){//linear case
            //integral = exp(b*xLo+c) * (exp(b*(xHi-xLo))-1) / b
            double expl = b / (Math.exp(b*(xHi-xLo))-1);
            
            if(expl==0)
                throw new RuntimeException("ERROR: UNDERFLOW");
            //treat second exp >> 1 as special case???
            
            c = Math.log(expl)-b*xLo;
        }
        else {
            double erfDiff = Erf.erf(cdErfArg(xLo), cdErfArg(xHi));
            if(erfDiff!=0){//erfDiff=0 prevents normalization; in this case will fail erfcInvNumericsTest
                double exponent = -Math.log( 0.5*erfDiff*Math.sqrt(-Math.PI/a) );
                c = exponent+0.25*b*b/a;
            }
        }
        
        double newNorm = cumulDistr(xHi);
        if(Double.isInfinite(newNorm) && Math.abs(a)>=1e-14 )
            return false;//quadratic can't be normalized, so will be replaced by linear
        else if( Math.abs(newNorm-1) > 1e-5 )
            throw new RuntimeException("ERROR: Unsuccessful normalization");
        
        return true;
    }
    
    private double cdErfArg(double x){
        double denom = 2*Math.sqrt(-a);
        double num = -b-2*a*x;
        return num/denom;
    }
    
    /*private boolean erfcInvNumericsOK(){
        //check if a quadratic distribution has OK numerics for erfcInv 
        //this fails if the minimum of the quadratic is way outside the voxel
        //in this case better use a linear representation
        for(double x : new double[] {xLo,xHi}){//if OK at both endpoints, should be OK in between
            if(!erfcInvNumericsOK(cdErfArg(x)))
                return false;
        }
        
        return true;
    }
    
    private boolean erfcInvNumericsOK(double x){
        //some points near edges of bell curve max out erfcInv
        if(x>0)
            return ! Double.isInfinite(myErfcInv(Erf.erfc(x)));
        else
            return ! Double.isInfinite(myErfcInv(Erf.erfc(-x)));
    }*/
    
    
    private boolean erfcInvNumericsOK(){
        //check if a quadratic distribution has OK numerics for erfcInv 
        //this fails if the minimum of the quadratic is way outside the voxel
        //in this case better use a linear representation
        
        double erfArgRange = Math.abs( cdErfArg(xHi) - cdErfArg(xLo) );
        
        for(double x : new double[] {xLo,xHi}){//if OK at both endpoints, should be OK in between
            if(!erfcInvNumericsOK(cdErfArg(x), erfArgRange))
                return false;
        }
        
        return true;
    }
    
    private boolean erfcInvNumericsOK(double x, double refRange){
        //some points near edges of bell curve make erfc inversion innaccurate.  Detect this.
        //refRange used to determine tolerance
        x = Math.abs(x);
        double invResult = myErfcInv(Erf.erfc(x));
        if(Double.isInfinite(invResult))//definitely going to be a problem
            return false;
        else if( Math.abs(x-invResult) > 0.01*refRange )
            return false;
        
        return true;
    }
    
    
    private double cumulDistrInv(double F){
        //inverse of cumulDistr        
        if(Math.abs(a)<1e-14){//linear case
            double ans = ( Math.log( b*F + Math.exp(b*xLo+c) ) - c ) / b;
            
            if(ans<xLo-1e-6 || ans>xHi+1e-6)//DEBUG!!!
                System.out.println("Out of range QuadraticQFunction draw...");//DEBUG!!!

            //double Fcheck = cumulDistr(ans);//DEBUG!!!
            return ans;
        }
        else {
            double C = 0.5*Math.exp(c-0.25*b*b/a)*Math.sqrt(-Math.PI/a);
            
            double erfArg1 = cdErfArg(xLo);
            double erfArg2;
            if(erfArg1>1){//these three options are mathematically equivalent,
                //but will be selected to best provide numerical stability
                double erfcVal = -F/C + Erf.erfc(erfArg1);
                erfArg2 = myErfcInv(erfcVal);
            }
            else if(erfArg1<-1){
                double oerfcVal = F/C + Erf.erfc(-erfArg1);
                erfArg2 = -myErfcInv(oerfcVal);
            }
            else {
                double erfVal = F/C + Erf.erf(erfArg1);
                erfArg2 = Erf.erfInv(erfVal);
            }
            
            double denom = 2*Math.sqrt(-a);
            double ans = (-erfArg2*denom - b)/(2*a);
                        
            //double Fcheck = cumulDistr(ans);//DEBUG!!!
            if(ans<xLo-1 || ans>xHi+1)//DEBUG!!!
                System.out.println("Out of range QuadraticQFunction draw...");//DEBUG!!!
                        
            return ans;
        }
    }
    
    
    private static double myErfcInv(double z){
        //erfcInv maxes out too early.  But I noticed that for x=4 to 27 (and probably much higher),
        // ln(erfc(x)) = -0.997021118x^2 - 0.166032837x - 1.474007920
        // with R² = 0.999999967
        //so we'll use that when regular erfcInv is infinite (or close to it)
        if(z>1)
            return -myErfcInv(2-z);
        
        if(z>1e-15)//I think Erf.erfcInv should be OK in this range
            return Erf.erfcInv(z);
        else
            return quadApproxErfcInv(z);
    }
    
    private static double quadApproxErfcInv(double z){
        double discr = 0.166032837*0.166032837 - 4 * 0.997021118 * (1.474007920+Math.log(z));
        if(discr<0)//this shouldn't happen with the small z we use here
            throw new RuntimeException("ERROR: erfc quad approx out of range!");
        
        return (Math.sqrt(discr)-0.166032837) / (2*0.997021118);
    }
    
}