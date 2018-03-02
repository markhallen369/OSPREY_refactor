/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.bibis;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.ematrix.epic.EPoly;
import edu.duke.cs.osprey.ematrix.epic.EPolyPC;
import edu.duke.cs.osprey.ematrix.epic.SAPE;
import edu.duke.cs.osprey.ematrix.epic.SeriesFitter;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import static java.util.stream.Collectors.toList;

/**
 *
 * Can be EPIC-like, LUTE-like, or sequence-based
 * 
 * @author mhall44
 */
public interface DenseFeatureSet extends Serializable {
    
    //when evaluating these per-sample functions, can assume the FeatureSet has called initSample
    //to set molecule geometry or whatever
    //and also that this DenseFeatureSet is at the right RCs to be active
    public abstract double evalEnergy(ConfSample samp, double[] coeffs, int offset, boolean freezeSAPE);
    //sum_j Aij coeffs_j (see FeatureSet, Aij is feature j evaluated at sample i)
    //offset is the first feature index (j) for this DenseFeatureSet

    public abstract double evalNoSAPEEnergy(ConfSample samp, double[] coeffs, int offset);
    
    public abstract void updateATProduct(ConfSample samp, double sampVal, double[] outVec, int offset, boolean freezeSAPE);
    //outVec_j += A_ij * sampVal
    
    public int getNumFeatures();
    public List<DegreeOfFreedom> getContDOFs();
    
    
    class Const implements DenseFeatureSet {
        //a feature set consisting of one feature, whose value is 1 (if this DenseFeatureSet is active
        //for the current RC's)
        @Override
        public double evalEnergy(ConfSample samp, double[] coeffs, int offset, boolean freezeSAPE) {
            return coeffs[offset];
        }
        @Override
        public double evalNoSAPEEnergy(ConfSample samp, double[] coeffs, int offset) {
            return coeffs[offset];
        }
        @Override
        public void updateATProduct(ConfSample samp, double sampVal, double[] outVec, int offset, boolean freezeSAPE) {
            outVec[offset] += sampVal;
        }

        @Override
        public int getNumFeatures() {
            return 1;
        }

        @Override
        public List<DegreeOfFreedom> getContDOFs() {
            return new ArrayList<>();//Const term doesn't depend on any cont dofs
        }
    }
    
    
    class EPIC implements DenseFeatureSet {
        //EPIC polynomial
        //Might have SAPE.  If so, make sure to initSample to get molecule in right conf
        //features are poly coeffs (SeriesFitter format, for dofs), followed by SAPE if not null
        SAPE sapeTerm;
        int order;
        ArrayList<DegreeOfFreedom> dofs;
        DoubleMatrix1D center;
        private int numFeatures;
        
        int[] dofIndices;//DEBUG!!!  let's see if we can speed things up this way
        
        HashMap<Integer,ResidueTemplate> moveableResTypes;
        
        
        DoubleMatrix1D dofVals;//DEBUG!!!!
        
        public EPIC(SAPE sapeTerm, int order, ArrayList<DegreeOfFreedom> dofs, ArrayList<Double> centerDOFVals,
                HashMap<Integer,ResidueTemplate> moveableResTypes, int[] dofIndices){
            this.sapeTerm = sapeTerm;
            this.order = order;
            this.dofs = dofs;
            this.moveableResTypes = moveableResTypes;
            calcNumFeatures();
            
            center = DoubleFactory1D.dense.make(centerDOFVals.size());
            for(int a=0; a<centerDOFVals.size(); a++)
                center.set(a, centerDOFVals.get(a));
            
            this.dofIndices = dofIndices;
            
            dofVals = DoubleFactory1D.dense.make(dofs.size());
        }
        
        public EPIC(EPoly p, ConfSpace cSpace, RCTuple tup){
            sapeTerm = p.sapeTerm;
            order = p.getOrder();
            if(p instanceof EPolyPC)
                throw new RuntimeException("ERROR: PC not supported yet");
            dofs = p.getDOFs();
            center = p.getCenter();
            
            calcNumFeatures();
            
            moveableResTypes = new HashMap<>();
            for(int a=0; a<tup.size(); a++){
                PositionConfSpace pcs = cSpace.posFlex.get(tup.pos.get(a));
                String resTypeName = pcs.RCs.get(tup.RCs.get(a)).AAType;
                ResidueTemplate resType = EnvironmentVars.resTemplates.getTemplateForMutation(resTypeName, pcs.res, true);
                int resIndex = pcs.res.indexInMolecule;
                moveableResTypes.put(resIndex, resType);
            }
            
            dofIndices = new int[dofs.size()];
            List<String> allDOFNames = cSpace.confDOFs.stream().map(d->d.getName()).collect(toList());
            for(int a=0; a<dofs.size(); a++)
                dofIndices[a] = allDOFNames.indexOf(dofs.get(a).getName());
            
            dofVals = DoubleFactory1D.dense.make(dofs.size());
        }
        
        private void calcNumFeatures(){
            numFeatures = SeriesFitter.getNumParams(dofs.size(), true, order);
            if(sapeTerm!=null)
                numFeatures++;
        }
        
        private DoubleMatrix1D assembleRelDOFVals(ConfSample samp){
            int numDOFs = dofs.size();
            //DoubleMatrix1D dofVals = DoubleFactory1D.dense.make(numDOFs);
            for(int d=0; d<numDOFs; d++)
                //dofVals.set(d, samp.contDOFVals.get(dofs.get(d).getName()));
                dofVals.set(d, samp.contDOFValsArr[dofIndices[d]] - center.getQuick(d));
            //dofVals.assign(center, Functions.minus);
            return dofVals;
        }
                
        
        @Override
        public double evalEnergy(ConfSample samp, double[] coeffs, int offset, boolean freezeSAPE) {
            double E = evalNoSAPEEnergy(samp,coeffs,offset);
                        
            if(sapeTerm!=null){
                double sapeCoeff = freezeSAPE ? 1 : coeffs[offset+numFeatures-1];
                E += sapeCoeff*sapeTerm.getEnergySharedMolec();//assuming initSample called...
            }

            return E;
        }
        
        @Override
        public double evalNoSAPEEnergy(ConfSample samp, double[] coeffs, int offset) {
            int numDOFs = dofs.size();
            DoubleMatrix1D relDOFVals = assembleRelDOFVals(samp);
            
            double polyCoeffs[] = new double[SeriesFitter.getNumParams(numDOFs, true, order)];
            System.arraycopy(coeffs, offset, polyCoeffs, 0, polyCoeffs.length);
            double E = SeriesFitter.evalSeries(polyCoeffs, relDOFVals, numDOFs, true, order);

            return E;
        }

        @Override
        public void updateATProduct(ConfSample samp, double sampVal, double[] outVec, int offset, boolean freezeSAPE) {
            int numDOFs = dofs.size();
            DoubleMatrix1D relDOFVals = assembleRelDOFVals(samp);
            
            int numMonomials = SeriesFitter.getNumParams(numDOFs, true, order);
            DoubleMatrix1D monomialVals = DoubleFactory1D.dense.make(numMonomials);
            SeriesFitter.calcSampParamCoeffs(monomialVals, relDOFVals, numDOFs, true, order, order, null);
            for(int m=0; m<numMonomials; m++)
                outVec[offset+m] += sampVal*monomialVals.get(m);
            
            if(sapeTerm!=null && !freezeSAPE)
                outVec[offset+numFeatures-1] += sampVal*sapeTerm.getEnergySharedMolec();//assuming initSample called...
        }

        @Override
        public int getNumFeatures() {
            return numFeatures;
        }

        @Override
        public List<DegreeOfFreedom> getContDOFs() {
            return dofs;
        }
    }
    
}





