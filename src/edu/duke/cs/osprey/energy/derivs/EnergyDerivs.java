/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.energy.derivs;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.ematrix.epic.EPICEnergyFunction;
import edu.duke.cs.osprey.ematrix.epic.EPoly;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ResPairEnergy;
import edu.duke.cs.osprey.energy.forcefield.SingleResEnergy;
import edu.duke.cs.osprey.energy.forcefield.SparseFFEnergy;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 *
 * This calculates the gradient and Hessian of the energy at a given conformation
 * in O(n^2) instead of O(n^2d) or O(n^2d^2) time for d DOFs
 * This is done via precomputing the gradient/Hessian of the atomic coordinates, 
 * which is O(nd) for grad, O(nd^2) for Hessian
 * (for EPIC terms, no conf stuff needed, and cost is O(max degree*cost of EPIC energy evaluation))
 * 
 * 
 * @author mhall44
 */
public class EnergyDerivs {
    
    MoleculeModifierAndScorer mms;//the energy as a function of internal coords
    
    HashMap<Residue,ResCoordDerivs> resCoordDerivs = new HashMap<>();
    
    DoubleMatrix1D curCoords;
    double E;
    DoubleMatrix1D Egrad = null;
    DoubleMatrix2D Ehess = null;
    
        
    
    double baseSAPE = 0;//for checking energy
    
    
    public EnergyDerivs(MoleculeModifierAndScorer mms, DoubleMatrix1D curCoords){
        this.mms = mms;
        this.curCoords = curCoords;
        precomputeCoordDerivs();
        calcEnergyAndDerivs();
    }
    
    private void precomputeCoordDerivs(){//gonna just use numerics here for now
        HashMap<Residue,ArrayList<Integer> > resDOFIndices = new HashMap<>();//indices (in mms) of the DOFs affecting each res
        
        ArrayList<DegreeOfFreedom> allDOFs = mms.getDOFs();
        
        for(int dof=0; dof<allDOFs.size(); dof++){
            List<Residue> affectedRes = allDOFs.get(dof).listAffectedResidues();
            for(Residue res : affectedRes){
                if( ! resDOFIndices.containsKey(res) )
                    resDOFIndices.put(res, new ArrayList<>());
                resDOFIndices.get(res).add(dof);
            }
        }
        
        for(Residue res : resDOFIndices.keySet()){
            resCoordDerivs.put(res, new ResCoordDerivs(res, resDOFIndices.get(res), mms));
        }
    }
    
    
    private void calcEnergyAndDerivs(){
        E = mms.getValue(curCoords);//this has two useful side effects:
        //sets the coordinates correctly, and initializes any non-initialized forcefield energies
        //sorry
        
        double Echeck = 0;//DEBUG!!!
        
        int n = curCoords.size();
        double grad[] = new double[n];
        double hess[][] = new double[n][n];
                
        baseSAPE = 0;
        
        for(ForcefieldEnergy ff: enumerateForcefields(mms.getEfunc()) ){
            double[][][] termsAndDerivs = ff.calculateEnergyTermsAndDerivs();
            Residue res1 = ff.getRes1();
            Residue res2 = ff.getRes2();
            ResCoordDerivs cd1 = getCoordDerivs(res1);
            ResCoordDerivs cd2 = getCoordDerivs(res2);
            
            ResPairEnergyDerivs ffDerivs = new ResPairEnergyDerivs(termsAndDerivs, cd1, cd2);
            Echeck += ffDerivs.E;
                        
            //add in gradient information from the term ff
            for(int dof=0; dof<ffDerivs.numDOFs; dof++){
                grad[ffDerivs.mmsDOFIndex(dof)] += ffDerivs.grad[dof];
                for(int dof2=0; dof2<ffDerivs.numDOFs; dof2++){
                    hess[ffDerivs.mmsDOFIndex(dof)][ffDerivs.mmsDOFIndex(dof2)] += ffDerivs.hess[dof][dof2];
                }
            }
        }
        
        if(mms.getEfunc() instanceof EPICEnergyFunction){
            EPICEnergyFunction eef = ((EPICEnergyFunction)mms.getEfunc());
            ArrayList<EPoly> terms = eef.getTerms();
            for(int t=0; t<terms.size(); t++){
                EPoly ep = terms.get(t);
                ArrayList<Integer> termDOFs = eef.termDOFs.get(t);

                DoubleMatrix1D coordsForTerm = DoubleFactory1D.dense.make(termDOFs.size());
                for(int dof=0; dof<termDOFs.size(); dof++)
                    coordsForTerm.set(dof, curCoords.get(termDOFs.get(dof)));
                DoubleMatrix1D termGrad = ep.polynomialGradient(coordsForTerm);
                DoubleMatrix2D termHess = ep.polynomialHessian(coordsForTerm);
                Echeck += ep.evalPolynomial(coordsForTerm);
                
                for(int dof=0; dof<termDOFs.size(); dof++){
                    grad[termDOFs.get(dof)] += termGrad.get(dof);
                    for(int dof2=0; dof2<termDOFs.size(); dof2++){
                        hess[termDOFs.get(dof)][termDOFs.get(dof2)] += termHess.get(dof,dof2);
                    }
                }
            }
        }
        
        Echeck -= baseSAPE;
        
        //DEBUG!!!  Also add in EPIC terms, others as needed
        
        Egrad = DoubleFactory1D.dense.make(grad);
        Ehess = DoubleFactory2D.dense.make(hess);
    }
    
    private ResCoordDerivs getCoordDerivs(Residue res){
        //get coord derivs, assuming they are already computed
        if(resCoordDerivs.containsKey(res))
            return resCoordDerivs.get(res);
        else//assuming coords are invariant to all DOFs
            return new ResCoordDerivs(res);
    }
    
    
    public double getEnergy(){
        return E;
    }
    
    public DoubleMatrix1D getGrad(){
        return Egrad;
    }
    
    public DoubleMatrix2D getHess(){
        return Ehess;
    }
    
    
    
    //DEBUG!!!  forcefield enumeration, should probably be made more general and include EPIC
    private ArrayList<ForcefieldEnergy> enumerateForcefields(EnergyFunction efunc){
        ArrayList<ForcefieldEnergy> ans = new ArrayList<>();
        enumerateForcefields(efunc, ans);
        return ans;
    }
    
    private void enumerateForcefields(EnergyFunction efunc, ArrayList<ForcefieldEnergy> list){
        if(efunc instanceof ResPairEnergy)
            list.add(((ResPairEnergy)efunc).getFFEnergy());
        else if(efunc instanceof SingleResEnergy)
            list.add(((SingleResEnergy)efunc).getFFEnergy());
        else if(efunc instanceof SparseFFEnergy)
            list.add(((SparseFFEnergy)efunc).getFFEnergy());
        else if(efunc instanceof EPICEnergyFunction){//take forcefields, polynomials handled separately
            if(((EPICEnergyFunction)efunc).includeMinE || !((EPICEnergyFunction)efunc).useSharedMolec)
                throw new RuntimeException("ERROR: Unsupported EPICEnergyFunction settings");
            for(EPoly ep : ((EPICEnergyFunction)efunc).getTerms() ){
                if(ep.sapeTerm!=null){
                    enumerateForcefields(ep.sapeTerm.sharedMolecEnergyFunction, list);
                    baseSAPE += ep.baseSAPE;
                }
            }
        }
        else if(efunc instanceof MultiTermEnergyFunction){
            ArrayList<Double> coeffs = ((MultiTermEnergyFunction)efunc).getCoeffs();
            ArrayList<EnergyFunction> terms = ((MultiTermEnergyFunction)efunc).getTerms();
            for(int t=0; t<terms.size(); t++){
                if(coeffs.get(t)!=1)
                    throw new RuntimeException("ERROR: not supporting nontrivial coeffs here");//but may have to...
                else
                    enumerateForcefields(terms.get(t), list);
            }
        }
    }
    
    
    
}
