/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ematrix.epic;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import cern.jet.math.Functions;
import static edu.duke.cs.osprey.ematrix.epic.EPICFitter.sampPerParam;
import edu.duke.cs.osprey.energy.derivs.EnergyDerivs;
import edu.duke.cs.osprey.bibis.BasicConfSampleSet;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.plug.RCTuplePolytope;
import edu.duke.cs.osprey.tools.ObjectIO;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import static java.util.stream.Collectors.toList;
import org.apache.commons.math3.optim.linear.LinearConstraint;

/**
 *
 * A version of EPIC fitter that uses derivatives to get coefficients
 * instead of just fitting the whole polynomial to undifferentiated energies
 * This is mainly a test case for the use of derivatives to split IBIS fits
 * Going to just deal with quadratics here
 * second derivatives are sufficient to split pairwise IBIS fits down to LUTE+EPIC size
 * 
 * Also sampling will be by PLUG
 * previous tests show this can be done effectively (even easier than usual) with rej sampling
 * but want to try LP-based sampling too as a test case for IBIS
 * 
 * @author mhall44
 */
public class DerivEPICFitter extends EPICFitter {
    
    
    //fitting is different but let's validate the same, for comparison
    //maybe even super.doFit and validate here?
    
    public DerivEPICFitter ( MoleculeModifierAndScorer mof, EPICSettings eset,
            DoubleMatrix1D cen, double me, RCTuplePolytope rtp) {
        super(mof, eset, cen, me, rtp);
    }
    
    
    
    
    
    public EPoly doFit(FitParams fp){

        int numParams = SeriesFitter.getNumParams(numDOFs, false, fp.order);
        int numSamples = sampPerParam * numParams;
        

        if(fp.PCOrder>fp.order || fp.order>2)
            throw new RuntimeException("ERROR: DerivEPICFitter not currently supporting "+fp.getDescription());
        
        //DEBUG!!! drawing these so the SAPE is generated the same as in the usual EPICFitter
        DoubleMatrix1D[] sampRel = new DoubleMatrix1D[numSamples];
        DoubleMatrix1D[] sampAbs = new DoubleMatrix1D[numSamples];
        double trueVal[] = new double[numSamples];

        generateSamples(numSamples,sampRel,sampAbs,trueVal,numSamples/2);
        
        
        boolean overCutoff[] = new boolean[numSamples];
        SAPE sapeTerm = null;
        double baseShift = 0;//SVE at center
        
        if(fp.SAPECutoff>0){         
            sapeTerm = new SAPE(objFcn,fp.SAPECutoff,sampAbs);
            
            baseShift = sapeTerm.getEnergyStandalone(center);//SAPE contribution at center
            
            for(int s=0; s<numSamples; s++){
                double shift = sapeTerm.getEnergyStandalone(sampAbs[s]);
                overCutoff[s]  = (trueVal[s]>es.EPICThresh1);
                trueVal[s] -= shift - baseShift;
            }
        }
       
        double[] seriesCoeffs = seriesCoeffsSingleDiff(numSamples, fp, overCutoff, sampRel, sampAbs, 
                sapeTerm, baseShift);
        
        EPoly ans = new EPoly(numDOFs, objFcn.getDOFs(), DOFmax, DOFmin, center,
                    minE, seriesCoeffs, fp.order, tope);
        
        //add SAPE term
        if(fp.SAPECutoff>0){
            ans.sapeTerm = sapeTerm;
            ans.baseSAPE = baseShift;
        }
        
        
        ans.fitDescription = fp.getDescription();
        
        return ans;
    }
    
    
    
    double[] seriesCoeffsCrossDeriv(int numSamples, FitParams fp, boolean[] overCutoff, DoubleMatrix1D[] sampRel,
            DoubleMatrix1D sampAbs[], SAPE sapeTerm, double baseSAPE){
        //get quadratic cross-terms by derivatives (in IBIS this would let them be separated out),
        //and the rest by fitting (in IBIS could fit differences to a single point, once cross-terms are known)
        
        double seriesCoeffs[] = new double[SeriesFitter.getNumParams(numDOFs, false, 2)];
        
        if(Double.isInfinite(fp.SAPECutoff)){//SAPE takes care of everything
                System.out.println("No fit needed: SAPE is full energy");
                //all 0's for polynomial is best, anything else is noise
                return seriesCoeffs;
        }
        else {
            //first attempt will be to take 10 pointz and all their derivs
            //if this is not good enough, draw more distinct pointz
            
            int numSampsToFit = 40;
            
            ArrayList<DoubleMatrix1D> sampsForDerivs = new ArrayList<>();
            ArrayList<DoubleMatrix1D> relSampsForDerivs = new ArrayList<>();
            
            /*for(int s=0; s<numSamples&&sampsForDerivs.size()<numSampsToFit; s++){//DEBUG!!!
                if(!overCutoff[s]){
                    sampsForDerivs.add(sampAbs[s]);
                    relSampsForDerivs.add(sampRel[s]);
                }
                if(s==numSamples-1)
                    throw new RuntimeException("ERROR: TOO FEW SUBCUTOFF SAMPLES");//we should have drawn enough...
            }*/
            //numSampsPerParam = sampsForDerivs.size();
            //DEBUG!!! Trying to get samples by LP
            for(int s=0; s<numSampsToFit; s++){
                DoubleMatrix1D x = drawLPSamp();
                sampsForDerivs.add(x);
                relSampsForDerivs.add(x.copy().assign(center,Functions.minus));
            }
            
            
            List<EnergyDerivs> eDerivs = sampsForDerivs.stream().map(x->new EnergyDerivs((MoleculeModifierAndScorer)objFcn, x)).collect(toList());
            
            List<EnergyDerivs> eDerivsSAPE = null;
            if(sapeTerm!=null){
                final MoleculeModifierAndScorer sapeMMS = sapeTerm.mofStandalone;
                eDerivsSAPE = sampsForDerivs.stream().map(x->new EnergyDerivs(sapeMMS, x)).collect(toList());
            }
            
            DoubleMatrix2D H = DoubleFactory2D.dense.make(numDOFs,numDOFs);//Hessian of our series
            //except we'll only include the cross terms
            
            
            //quadratic terms first
            //int curCoeffIndex = numDOFs;//coefficient at which quadratic terms start
            for(int i=0; i<numDOFs; i++){//order has to match the silly old SeriesFitter
                for(int j=0; j<i; j++){
                    //fit just the one parameter
                    //the best fit to y=a is a=avg(observed y).  
                    //you can weight by x if you want.  not sure if it matters tho
                    
                    double avgDeriv2 = 0;
                    double wtsum = 0;
                    for(int c=0; c<numSampsToFit; c++){
                        double deriv2 = eDerivs.get(c).getHess().get(i,j);
                        if(sapeTerm!=null)
                            deriv2 -= eDerivsSAPE.get(c).getHess().get(i,j);
                        double wt = 1;//trueVal[c]>1 ? 1./trueVal[c] : 1.;//DEBUG!!
                        avgDeriv2 += wt*deriv2;
                        wtsum += wt;
                    }
                    avgDeriv2 /= wtsum;
                    
                    H.set(i, j, avgDeriv2);
                    H.set(j, i, avgDeriv2);
                }
            }
            
            
            //fit A*seriesCoeffs = b
            DoubleMatrix2D A = DoubleFactory2D.dense.make(numSampsToFit, 2*numDOFs);
            DoubleMatrix1D b = DoubleFactory1D.dense.make(numSampsToFit);
                   
            
            for(int s=0; s<numSampsToFit; s++){
                for(int i=0; i<numDOFs; i++)
                    A.set(s, i, relSampsForDerivs.get(s).get(i));
                for(int i=0; i<numDOFs; i++){
                    A.set(s, numDOFs+i, relSampsForDerivs.get(s).get(i)*relSampsForDerivs.get(s).get(i));
                }                
                double crossTermSum = 0.5 * Algebra.DEFAULT.mult(H,relSampsForDerivs.get(s)).zDotProduct(relSampsForDerivs.get(s));
                
                double trueVal = eDerivs.get(s).getEnergy() - minE;
                if(sapeTerm!=null)
                    trueVal -= eDerivsSAPE.get(s).getEnergy() - baseSAPE;
                
                b.set(s, trueVal - crossTermSum);
            }
            
            double nonCrossSer[] = lsqSolve(A,b);
            for(int i=0; i<numDOFs; i++)
                seriesCoeffs[i] = nonCrossSer[i];
            int coeffCount = numDOFs;
            for(int i=0; i<numDOFs; i++){
                for(int j=0; j<=i; j++){
                    seriesCoeffs[coeffCount] = (i==j) ? nonCrossSer[numDOFs+i] : H.get(i,j);
                    coeffCount++;
                }
            }
        }
        
        return seriesCoeffs;
    }
    
    
    
    
    
    double[] seriesCoeffsSingleDeriv(int numSamples, FitParams fp, boolean[] overCutoff, DoubleMatrix1D[] sampRel,
            DoubleMatrix1D sampAbs[], SAPE sapeTerm, double baseSAPE){
        //fit quadratic cross-terms to cross-derivatives,
        //and the rest to 1st derivatives
        //this one does not work well (noticeably worse than cross derivs, quadratic not good enough even for regular 1CC8,
        //like byDerivs)
        
        double seriesCoeffs[] = new double[SeriesFitter.getNumParams(numDOFs, false, 2)];
        
        if(Double.isInfinite(fp.SAPECutoff)){//SAPE takes care of everything
                System.out.println("No fit needed: SAPE is full energy");
                //all 0's for polynomial is best, anything else is noise
                return seriesCoeffs;
        }
        else {
            //first attempt will be to take 10 pointz and all their derivs
            //if this is not good enough, draw more distinct pointz
            
            int numSampsToFit = 40;
            
            ArrayList<DoubleMatrix1D> sampsForDerivs = new ArrayList<>();
            ArrayList<DoubleMatrix1D> relSampsForDerivs = new ArrayList<>();
            
            /*for(int s=0; s<numSamples&&sampsForDerivs.size()<numSampsToFit; s++){//DEBUG!!!
                if(!overCutoff[s]){
                    sampsForDerivs.add(sampAbs[s]);
                    relSampsForDerivs.add(sampRel[s]);
                }
                if(s==numSamples-1)
                    throw new RuntimeException("ERROR: TOO FEW SUBCUTOFF SAMPLES");//we should have drawn enough...
            }*/
            //numSampsPerParam = sampsForDerivs.size();
            //DEBUG!!! Trying to get samples by LP
            for(int s=0; s<numSampsToFit; s++){
                DoubleMatrix1D x = drawLPSamp();
                sampsForDerivs.add(x);
                relSampsForDerivs.add(x.copy().assign(center,Functions.minus));
            }
            
            
            List<EnergyDerivs> eDerivs = sampsForDerivs.stream().map(x->new EnergyDerivs((MoleculeModifierAndScorer)objFcn, x)).collect(toList());
            
            List<EnergyDerivs> eDerivsSAPE = null;
            if(sapeTerm!=null){
                final MoleculeModifierAndScorer sapeMMS = sapeTerm.mofStandalone;
                eDerivsSAPE = sampsForDerivs.stream().map(x->new EnergyDerivs(sapeMMS, x)).collect(toList());
            }
            
            DoubleMatrix2D H = DoubleFactory2D.dense.make(numDOFs,numDOFs);//Hessian of our series
            //except we'll only include the cross terms
            
            
            //quadratic terms first
            //int curCoeffIndex = numDOFs;//coefficient at which quadratic terms start
            for(int i=0; i<numDOFs; i++){//order has to match the silly old SeriesFitter
                for(int j=0; j<i; j++){
                    //fit just the one parameter
                    //the best fit to y=a is a=avg(observed y).  
                    //you can weight by x if you want.  not sure if it matters tho
                    
                    double avgDeriv2 = 0;
                    double wtsum = 0;
                    for(int c=0; c<numSampsToFit; c++){
                        double deriv2 = eDerivs.get(c).getHess().get(i,j);
                        if(sapeTerm!=null)
                            deriv2 -= eDerivsSAPE.get(c).getHess().get(i,j);
                        double wt = 1;//trueVal[c]>1 ? 1./trueVal[c] : 1.;//DEBUG!!
                        avgDeriv2 += wt*deriv2;
                        wtsum += wt;
                    }
                    avgDeriv2 /= wtsum;
                    
                    H.set(i, j, avgDeriv2);
                    H.set(j, i, avgDeriv2);
                }
            }
            
            double selfQuadCoeffs[] = new double[numDOFs];//helpful to keep this separate from H
            
            for(int i=0; i<numDOFs; i++){//fit the two parameters for DOF i to dE/dx_i
                DoubleMatrix2D A = DoubleFactory2D.dense.make(numSampsToFit, 2);
                DoubleMatrix1D b = DoubleFactory1D.dense.make(numSampsToFit);


                for(int s=0; s<numSampsToFit; s++){
                    A.set(s, 0, 1);
                    A.set(s, 1, 2*relSampsForDerivs.get(s).get(i));
                    
                    double crossTermContrib = H.viewRow(i).zDotProduct(relSampsForDerivs.get(s));
                    double trueDeriv = eDerivs.get(s).getGrad().get(i);
                    if(sapeTerm!=null)
                        trueDeriv -= eDerivsSAPE.get(s).getGrad().get(i);

                    b.set(s, trueDeriv - crossTermContrib);
                }
                
                double nonCrossSer[] = lsqSolve(A,b);
                seriesCoeffs[i] = nonCrossSer[0];
                selfQuadCoeffs[i] = nonCrossSer[1];
            }
            
            
            //finally get the cross derivs into seriesCoeffs
            int coeffCount = numDOFs;
            for(int i=0; i<numDOFs; i++){
                for(int j=0; j<=i; j++){
                    seriesCoeffs[coeffCount] = (i==j) ? selfQuadCoeffs[i] : H.get(i,j);
                    coeffCount++;
                }
            }
        }
        
        return seriesCoeffs;
    }
    
    
    
    double[] seriesCoeffsSingleDiff(int numSamples, FitParams fp, boolean[] overCutoff, DoubleMatrix1D[] sampRel,
            DoubleMatrix1D sampAbs[], SAPE sapeTerm, double baseSAPE){
        //fit quadratic cross-terms to cross-derivatives,
        //then fit the self-terms for each DOF separately to difference involving that DOF
        //(and thus independent of all the other cross terms)
        
        double seriesCoeffs[] = new double[SeriesFitter.getNumParams(numDOFs, false, 2)];
        
        if(Double.isInfinite(fp.SAPECutoff)){//SAPE takes care of everything
                System.out.println("No fit needed: SAPE is full energy");
                //all 0's for polynomial is best, anything else is noise
                return seriesCoeffs;
        }
        else {
            //first attempt will be to take 10 pointz and all their derivs
            //if this is not good enough, draw more distinct pointz
            
            int numSampsToFit = 40;
            
            ArrayList<DoubleMatrix1D> sampsForDerivs = new ArrayList<>();
            ArrayList<DoubleMatrix1D> relSampsForDerivs = new ArrayList<>();
            
            /*for(int s=0; s<numSamples&&sampsForDerivs.size()<numSampsToFit; s++){//DEBUG!!!
                if(!overCutoff[s]){
                    sampsForDerivs.add(sampAbs[s]);
                    relSampsForDerivs.add(sampRel[s]);
                }
                if(s==numSamples-1)
                    throw new RuntimeException("ERROR: TOO FEW SUBCUTOFF SAMPLES");//we should have drawn enough...
            }*/
            //numSampsPerParam = sampsForDerivs.size();
            //DEBUG!!! Trying to get samples by LP
            for(int s=0; s<numSampsToFit; s++){
                DoubleMatrix1D x = drawLPSamp();
                sampsForDerivs.add(x);
                relSampsForDerivs.add(x.copy().assign(center,Functions.minus));
            }
            
            
            List<EnergyDerivs> eDerivs = sampsForDerivs.stream().map(x->new EnergyDerivs((MoleculeModifierAndScorer)objFcn, x)).collect(toList());
            
            List<EnergyDerivs> eDerivsSAPE = null;
            if(sapeTerm!=null){
                final MoleculeModifierAndScorer sapeMMS = sapeTerm.mofStandalone;
                eDerivsSAPE = sampsForDerivs.stream().map(x->new EnergyDerivs(sapeMMS, x)).collect(toList());
            }
            
            DoubleMatrix2D H = DoubleFactory2D.dense.make(numDOFs,numDOFs);//Hessian of our series
            //except we'll only include the cross terms
            
            
            //quadratic terms first
            //int curCoeffIndex = numDOFs;//coefficient at which quadratic terms start
            for(int i=0; i<numDOFs; i++){//order has to match the silly old SeriesFitter
                for(int j=0; j<i; j++){
                    //fit just the one parameter
                    //the best fit to y=a is a=avg(observed y).  
                    //you can weight by x if you want.  not sure if it matters tho
                    
                    double avgDeriv2 = 0;
                    double wtsum = 0;
                    for(int c=0; c<numSampsToFit; c++){
                        double deriv2 = eDerivs.get(c).getHess().get(i,j);
                        if(sapeTerm!=null)
                            deriv2 -= eDerivsSAPE.get(c).getHess().get(i,j);
                        double wt = 1;//trueVal[c]>1 ? 1./trueVal[c] : 1.;//DEBUG!!
                        avgDeriv2 += wt*deriv2;
                        wtsum += wt;
                    }
                    avgDeriv2 /= wtsum;
                    
                    H.set(i, j, avgDeriv2);
                    H.set(j, i, avgDeriv2);
                }
            }
            
            double selfQuadCoeffs[] = new double[numDOFs];//helpful to keep this separate from H
            
            for(int i=0; i<numDOFs; i++){//fit the two parameters for DOF i to dE/dx_i
                DoubleMatrix2D A = DoubleFactory2D.dense.make(numSampsToFit, 2);
                DoubleMatrix1D b = DoubleFactory1D.dense.make(numSampsToFit);


                for(int s=0; s<numSampsToFit; s++){
                    double altx = drawAltDOFVal(sampsForDerivs.get(s),i);//alt value for this degree of freedom still within polytope
                    DoubleMatrix1D altRelSamp = relSampsForDerivs.get(s).copy();
                    altRelSamp.set(i, altx-center.get(i));
                    DoubleMatrix1D altAbsSamp = sampsForDerivs.get(s).copy();
                    altAbsSamp.set(i,altx);
                    
                    A.set(s, 0, sampsForDerivs.get(s).get(i)-altx);
                    A.set(s, 1, relSampsForDerivs.get(s).get(i)*relSampsForDerivs.get(s).get(i)-altRelSamp.get(i)*altRelSamp.get(i));
                    
                    double crossTermDiff = 0.5 * Algebra.DEFAULT.mult(H,relSampsForDerivs.get(s)).zDotProduct(relSampsForDerivs.get(s))
                            - 0.5 * Algebra.DEFAULT.mult(H,altRelSamp).zDotProduct(altRelSamp);
                                        
                    double trueDiff = eDerivs.get(s).getEnergy() - objFcn.getValue(altAbsSamp);
                    if(sapeTerm!=null)
                        trueDiff -= eDerivsSAPE.get(s).getEnergy() - sapeTerm.mofStandalone.getValue(altAbsSamp);

                    b.set(s, trueDiff - crossTermDiff);
                }
                
                double nonCrossSer[] = lsqSolve(A,b);
                seriesCoeffs[i] = nonCrossSer[0];
                selfQuadCoeffs[i] = nonCrossSer[1];
            }
            
            
            //finally get the cross derivs into seriesCoeffs
            int coeffCount = numDOFs;
            for(int i=0; i<numDOFs; i++){
                for(int j=0; j<=i; j++){
                    seriesCoeffs[coeffCount] = (i==j) ? selfQuadCoeffs[i] : H.get(i,j);
                    coeffCount++;
                }
            }
        }
        
        return seriesCoeffs;
    }
    
    
    
    
    
    private DoubleMatrix1D drawLPSamp(){
        ArrayList<LinearConstraint> constr = tope.getConstr();
        ArrayList<Double> ansList = BasicConfSampleSet.sampleFeasPt(constr);
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(ansList.size());
        for(int i=0; i<ansList.size(); i++)
            ans.set(i, ansList.get(i));
        return ans;
    }
    
    private double drawAltDOFVal(DoubleMatrix1D x, int dof){
        //draw a random value for dof to be plugged into x
        //It should stay in the polytope
        double[] altx = x.copy().toArray();
        while(true){//assuming x is in tope so loop should eventually terminate
            double altVal = DOFmin.get(dof) + Math.random()*(DOFmax.get(dof)-DOFmin.get(dof));
            altx[dof] = altVal;
            if(tope.containsPoint(altx))
                return altVal;
        }
    }
    
    
    
    double[] seriesCoeffsByFit(int numSamples, FitParams fp, boolean[] overCutoff, DoubleMatrix1D[] sampRel,
            DoubleMatrix1D sampAbs[], double[] trueVal, SAPE sapeTerm){
        //get coeffs entirely by fitting, should match uje EPICFitter
        
        int numParams = SeriesFitter.getNumParams(numDOFs, false, 2);
        
        if(Double.isInfinite(fp.SAPECutoff)){//SAPE takes care of everything
                System.out.println("No fit needed: SAPE is full energy");
                //all 0's for polynomial is best, anything else is noise
                return new double[numParams];
        }
        else {
            //first attempt will be to take 10 pointz and all their derivs
            //if this is not good enough, draw more distinct pointz
            
            int numSampsToFit = 40;
            
            ArrayList<DoubleMatrix1D> sampsForDerivs = new ArrayList<>();
            ArrayList<DoubleMatrix1D> relSampsForDerivs = new ArrayList<>();
            for(int s=0; s<numSamples&&sampsForDerivs.size()<numSampsToFit; s++){//DEBUG!!!
                if(!overCutoff[s]){
                    sampsForDerivs.add(sampAbs[s]);
                    relSampsForDerivs.add(sampRel[s]);
                }
                if(s==numSamples-1)
                    throw new RuntimeException("ERROR: TOO FEW SUBCUTOFF SAMPLES");//we should have drawn enough...
            }
            //numSampsPerParam = sampsForDerivs.size();
            
            //fit M*seriesCoeffs = b
            DoubleMatrix2D A = DoubleFactory2D.dense.make(numSampsToFit, numParams);
            DoubleMatrix1D b = DoubleFactory1D.dense.make(numSampsToFit);
                   
            
            for(int s=0; s<numSampsToFit; s++){
                for(int i=0; i<numDOFs; i++)
                    A.set(s, i, relSampsForDerivs.get(s).get(i));
                int coeffCount = numDOFs;
                for(int i=0; i<numDOFs; i++){
                    for(int j=0; j<=i; j++){
                        A.set(s, coeffCount, relSampsForDerivs.get(s).get(i)*relSampsForDerivs.get(s).get(j));
                        coeffCount++;
                    }
                }
                b.set(s, trueVal[s]);
            }
            
            return lsqSolve(A,b);
        }
    }
    
    
    
    double[] lsqSolve(DoubleMatrix2D A, DoubleMatrix1D b){
        //solve Ax=b for x in the least-squares sense
        //Adapted from SeriesFitter
        DoubleMatrix2D At = Algebra.DEFAULT.transpose(A);
        DoubleMatrix2D M = Algebra.DEFAULT.mult(At, A);
        
        DoubleMatrix2D C = DoubleFactory2D.dense.make(A.columns(), 1);
        C.viewColumn(0).assign( Algebra.DEFAULT.mult(At, b) );
                
        try {
            return Algebra.DEFAULT.solve(M,C).viewColumn(0).toArray();
        }
        catch(IllegalArgumentException e){//indices singular M
            //solve in a way robust to singularities
            //basically pseudoinverse(M)*C
            SingularValueDecomposition svd = new SingularValueDecomposition(M);
            //M = U*S*V', so M^-1 = V*invS*U'
            DoubleMatrix2D invS = svd.getS().copy();
            
            //this tolerance is from SingularValueDecomposition.java (used there to compute rank)
            //here we use the special case that M is square
            double eps = Math.pow(2.0,-52.0);
            double tol = invS.rows()*invS.get(0,0)*eps;
            
            for(int i=0; i<invS.rows(); i++){
                double singVal = invS.get(i,i);
                if(singVal>tol)
                    invS.set(i, i, 1./singVal);
                else
                    invS.set(i, i, 0);
            }
            
            DoubleMatrix2D ansCol = Algebra.DEFAULT.mult( Algebra.DEFAULT.mult(
                    Algebra.DEFAULT.mult( svd.getV(), invS ), Algebra.DEFAULT.transpose(svd.getU())), C );
            
            
            //DEBUG!!!!
            //DoubleMatrix2D checkC = Algebra.DEFAULT.mult(M, ansCol);
            
            return ansCol.viewColumn(0).toArray();
        }
    }
    
    
    double[] seriesCoeffsByDerivs(int numSamples, FitParams fp, boolean[] overCutoff, DoubleMatrix1D[] sampRel,
            DoubleMatrix1D sampAbs[], double[] trueVal, SAPE sapeTerm){
        //get coeffs one param at a time by using derivatives
        
        //ok now just do the quadratic fit
        double seriesCoeffs[] = new double[SeriesFitter.getNumParams(numDOFs, false, 2)];
        
        
        if(Double.isInfinite(fp.SAPECutoff)){//SAPE takes care of everything
                System.out.println("No fit needed: SAPE is full energy");
                //all 0's for polynomial is best, anything else is noise
        }
        else {
            //first attempt will be to take 10 pointz and all their derivs
            //if this is not good enough, draw more distinct pointz
            
            int numSampsPerParam = 10;
            
            ArrayList<DoubleMatrix1D> sampsForDerivs = new ArrayList<>();
            ArrayList<DoubleMatrix1D> relSampsForDerivs = new ArrayList<>();
            for(int s=0; s<numSamples&&sampsForDerivs.size()<numSampsPerParam; s++){//DEBUG!!!
                if(!overCutoff[s]){
                    sampsForDerivs.add(sampAbs[s]);
                    relSampsForDerivs.add(sampRel[s]);
                }
                if(s==numSamples-1)
                    throw new RuntimeException("ERROR: TOO FEW SUBCUTOFF SAMPLES");//we should have drawn enough...
            }
            //numSampsPerParam = sampsForDerivs.size();
            
            
            List<EnergyDerivs> eDerivs = sampsForDerivs.stream().map(x->new EnergyDerivs((MoleculeModifierAndScorer)objFcn, x)).collect(toList());
            
            List<EnergyDerivs> eDerivsSAPE = null;
            if(sapeTerm!=null){
                final MoleculeModifierAndScorer sapeMMS = sapeTerm.mofStandalone;
                eDerivsSAPE = sampsForDerivs.stream().map(x->new EnergyDerivs(sapeMMS, x)).collect(toList());
            }
            
            DoubleMatrix2D H = DoubleFactory2D.dense.make(numDOFs,numDOFs);//Hessian of our series
            
            //quadratic terms first
            int curCoeffIndex = numDOFs;//coefficient at which quadratic terms start
            for(int i=0; i<numDOFs; i++){//order has to match the silly old SeriesFitter
                for(int j=0; j<=i; j++){
                    //fit just the one parameter
                    //the best fit to y=a is a=avg(observed y).  
                    //you can weight by x if you want.  not sure if it matters tho
                    
                    double avgDeriv2 = 0;
                    double wtsum = 0;
                    for(int c=0; c<numSampsPerParam; c++){
                        double deriv2 = eDerivs.get(c).getHess().get(i,j);
                        if(sapeTerm!=null)
                            deriv2 -= eDerivsSAPE.get(c).getHess().get(i,j);
                        double wt = trueVal[c]>1 ? 1./trueVal[c] : 1.;//DEBUG!!
                        avgDeriv2 += wt*deriv2;
                        wtsum += wt;
                    }
                    avgDeriv2 /= wtsum;
                    
                    H.set(i, j, avgDeriv2);
                    H.set(j, i, avgDeriv2);
                    
                    seriesCoeffs[curCoeffIndex] = (i==j) ? avgDeriv2/2 : avgDeriv2;
                    curCoeffIndex++;
                }
            }
            
            //now linear terms
            for(int i=0; i<numDOFs; i++){
                //grad = g + H*z.  since we know H we just take g = average(grad-Hz) where z is rel DOF vals
                
                double avgDeriv = 0;
                double wtsum = 0;
                for(int c=0; c<numSampsPerParam; c++){
                    double deriv = eDerivs.get(c).getGrad().get(i) - H.viewRow(i).zDotProduct(relSampsForDerivs.get(c));
                    if(sapeTerm!=null)
                        deriv -= eDerivsSAPE.get(c).getGrad().get(i);
                    double wt = trueVal[c]>1 ? 1./trueVal[c] : 1.;//DEBUG!!
                    avgDeriv += wt*deriv;
                    wtsum += wt;
                }
                avgDeriv /= wtsum;

                seriesCoeffs[i] = avgDeriv;
            }
            
            //For comparability to normal EPICFitter, don't adjust constant term
        }
        
        return seriesCoeffs;
    }
    
    
    
}
