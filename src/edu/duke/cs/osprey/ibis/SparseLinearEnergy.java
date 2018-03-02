/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ibis;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPoly;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tupexp.TESampleSet;
import java.util.ArrayList;

/**
 *
 * An estimate of some kind of energy as a linear function of conformational features
 * typically features will be 0 for some RC tuples; if nonzero may be a monomial or SAPE term
 * dependent on continuous conformational DOFs or just a constant
 * Can represent conformational energy or free energy (integrated over some DOFs) as a func of other DOFs
 * 
 * 
 * @author mhall44
 */
public class SparseLinearEnergy  {
    
    FeatureSet featSet;
    double[] coeffs;
    
    public SparseLinearEnergy(FeatureSet featSet, double[] coeffs){
        this.featSet = featSet;
        this.coeffs = coeffs;
    }
    
    public SparseLinearEnergy(FeatureSet featSet, ConfSampleSet trainingSet, double initCoeffs[]){
        //Fit coeffs to match the energies in trainingSet.  initCoeffs (initial guess) can be null
        this.featSet = featSet;
        SLEFitter fitter = makeFitter(featSet, trainingSet);
        coeffs = fitter.fitCoeffs(initCoeffs);
    }
    
    public SparseLinearEnergy(EPICMatrix epicMat, PruningMatrix pruneMat, PolytopeMatrix plugMat){
        //The EPIC energy is sparse and linear
        //this initialization gives us a good feature set for the unintegrated energy
        //which we can use to make good feature sets for integrated energies
        //if using a >pairwise E-func can use AMBER+EPIC to get an idea of what features are needed
        //then re-fit this energy
        featSet = new FeatureSet(epicMat.getConfSpace());//initialize the feature set
        ArrayList<Double> coeffList = new ArrayList<>();

        VoxDrawTupleExpander te = new VoxDrawTupleExpander(pruneMat,plugMat,featSet);
        TESampleSet tss = new TESampleSet(te);//for checking if tuples are feasible. DEBUG!! don't necessarily need tuple expander machinery for this
        
        if(epicMat.hasHigherOrderTerms())
            throw new RuntimeException("ERROR: Higher-order terms not supported currently in EPIC matrix");
        
        for(int pos=0; pos<epicMat.getNumPos(); pos++){
            for(int rc=0; rc<epicMat.getNumConfAtPos(pos); rc++){
                addEPICTermIfPresent(epicMat, new RCTuple(pos,rc), tss, coeffList);
                for(int pos2=0; pos2<pos; pos2++){
                    for(int rc2=0; rc2<epicMat.getNumConfAtPos(pos2); rc2++){
                        addEPICTermIfPresent(epicMat, new RCTuple(pos,rc,pos2,rc2), tss, coeffList);
                    }
                }
            }
        }
        
        coeffs = new double[coeffList.size()];
        for(int c=0; c<coeffList.size(); c++)
            coeffs[c] = coeffList.get(c);
    }
    
    
    public SparseLinearEnergy(EnergyMatrix emat, PruningMatrix pruneMat, PolytopeMatrix plugMat){
        //This is for taking a pairwise LUTE matrix and testing ability to convert to K* scores
        ConfSpace confSpace = plugMat.cSpace;
        featSet = new FeatureSet(confSpace);//initialize the feature set
        ArrayList<Double> coeffList = new ArrayList<>();

        if(emat.hasHigherOrderTerms())
            throw new RuntimeException("ERROR: Higher-order terms not supported currently in EPIC matrix");
        
        VoxDrawTupleExpander te = new VoxDrawTupleExpander(pruneMat,plugMat,featSet);
        TESampleSet tss = new TESampleSet(te);//for checking if tuples are feasible. DEBUG!! don't necessarily need tuple expander machinery for this

        
        for(int pos=0; pos<emat.getNumPos(); pos++){
            for(int rc=0; rc<emat.getNumConfAtPos(pos); rc++){
                double E = emat.getOneBody(pos,rc);
                if(Double.isFinite(E)){
                    if(tss.tupleFeasible(new RCTuple(pos,rc))){
                        //DEBUG!!!  need to figure out how to do dfs w/ plug.  Currently tuple feas checks are inconsistent...else checking finite would suffice
                        featSet.addDenseFeatureSet(new DenseFeatureSet.Const(), new RCTuple(pos,rc));
                        coeffList.add(E);
                    }
                }
                for(int pos2=0; pos2<pos; pos2++){
                    for(int rc2=0; rc2<emat.getNumConfAtPos(pos2); rc2++){
                        E = emat.getPairwise(pos, rc, pos2, rc2);
                        if(Double.isFinite(E)){
                            if(tss.tupleFeasible(new RCTuple(pos,rc,pos2,rc2))){
                                featSet.addDenseFeatureSet(new DenseFeatureSet.Const(), new RCTuple(pos,rc,pos2,rc2));
                                coeffList.add(E);
                            }
                        }
                    }
                }
            }
        }
        
        coeffs = new double[coeffList.size()];
        for(int c=0; c<coeffList.size(); c++)
            coeffs[c] = coeffList.get(c);
    }
    
    
    
    private void addEPICTermIfPresent(EPICMatrix epicMat, RCTuple tup, TESampleSet tss, ArrayList<Double> coeffList){
        //If tup can feasibly be found in good conformations, add its EPIC term to featSet and coeffList
        EPoly ep = epicMat.getTupleValue(tup);
        if(ep!=null){
            if(tss.tupleFeasible(tup)){
                coeffList.add(ep.getMinE()-ep.baseSAPE);//EPIC polys have this implied constant term
                if(ep.getDOFs().isEmpty())
                    featSet.addDenseFeatureSet(new DenseFeatureSet.Const(), tup);
                else {
                    featSet.addDenseFeatureSet(new DenseFeatureSet.EPIC(ep,epicMat.getConfSpace(),tup), tup);
                    for(double co : ep.getCoeffs())
                        coeffList.add(co);
                    if(ep.sapeTerm!=null && !featSet.freezeSAPE)
                        coeffList.add(1.);//SAPE term has coeff of 1
                }
            }
        }
    }
    
    
    public double evalEnergy(ConfSample samp){
        return featSet.evalEnergy(samp,coeffs);
    }
    
    
    public double evalEnergy(ConfSample samp, MoleculeModifierAndScorer mms) {
        return featSet.evalEnergy(samp,coeffs, mms);
    }
        
    
    private SLEFitter makeFitter(FeatureSet featSet, ConfSampleSet trainingSet){
        
        //DEBUG!!!
        if(!featSet.hasRedundantFeatures)
            return new TupIterFitter(featSet, trainingSet);
        //should be fast, but doesn't support redundant features
        //should check before fitting redundant EPIC anyway!  Discrete matrix doesn't iterate anyway
        
        if(featSet.freezeSAPE)
            return new SLEFitter.FrozenSAPECGFitter(featSet, trainingSet);
        else
            return new SLEFitter.ContBatchCGFitter(featSet, trainingSet);
    }
    
    
    public SparseLinearEnergy integrateModel(IntegrableDOF dof, ConfSampleSet samples){
        FeatureSet newFeat = dof.featureSetForInteg(featSet);
        double[] initCoeffs = initCoeffsForNewFeat(newFeat);
        SparseLinearEnergy newF = new SparseLinearEnergy(newFeat, samples, initCoeffs);
        return newF;
    }
    
    double[] initCoeffsForNewFeat(FeatureSet newFeat){
        double[] ans = new double[newFeat.numFeatures];
        for(RCTuple tup : newFeat.unprunedFeatTuples()){
            DenseFeatureSet oldFS = featSet.featMatrix.getTupleValue(tup);
            if(oldFS!=null){//if feature doesn't appear in current expansion, initialize its coeff to 0 in new one
                DenseFeatureSet newFS = newFeat.featMatrix.getTupleValue(tup);
                int oldOffset = featSet.featOffsets.getTupleValue(tup);
                int newOffset = newFeat.featOffsets.getTupleValue(tup);
                copyEquivalentCoeffs(ans, oldFS, newFS, oldOffset, newOffset);
            }
        }
        return ans;
    }
    
    private void copyEquivalentCoeffs(double[] newCoeffs, DenseFeatureSet oldFS, DenseFeatureSet newFS, int oldOffset, int newOffset){
        //Use the coeffs for oldFS to make appropriate initial-guess coeffs for newFS, put them in newCoeffs
        //DEBUG!!!  Might need to change this if make new kinds of feature sets.  Could also potentially map DOFs?  Not sure
        if(newFS.getNumFeatures()==oldFS.getNumFeatures()){//hopefully the same or very similar
            int numFeat = newFS.getNumFeatures();
            if(featSet.freezeSAPE && featSet.extractSAPE(newFS)!=null)
                numFeat--;
            System.arraycopy(coeffs, oldOffset, newCoeffs, newOffset, numFeat);
        }
        else
            newCoeffs[newOffset] = coeffs[oldOffset];//Copy the leading term; if there are other coeffs they can be refit starting from 0
    }
    
    
    public EnergyMatrix asEnergyMatrix(){
        //DEBUG!!!  This is quite similar to featSet.unprunedFeatTuples()
        //but doing it separately to get the null dense feature sets and set their energies to inf
        TupleMatrixGeneric<DenseFeatureSet> fmat = featSet.featMatrix;

        EnergyMatrix ans = new EnergyMatrix(fmat.getNumPos(), fmat.getNumConfAtPos(), Double.POSITIVE_INFINITY);
        
        if(fmat.hasHigherOrderTerms())
            throw new RuntimeException("ERROR: Higher-order terms not yet supported here");
        
        for(int pos=0; pos<fmat.getNumPos(); pos++){
            for(int rc=0; rc<fmat.getNumConfAtPos(pos); rc++){
                RCTuple single = new RCTuple(pos,rc);
                ans.setTupleValue(single, extractTupleEnergy(single));
                for(int pos2=0; pos2<pos; pos2++){
                    for(int rc2=0; rc2<fmat.getNumConfAtPos(pos2); rc2++){
                        RCTuple pair = new RCTuple(pos,rc,pos2,rc2);
                        ans.setTupleValue(pair, extractTupleEnergy(pair));
                    }
                }
            }
        }
        return ans;
    }
    
    double extractTupleEnergy(RCTuple tup){
        //convert contribution of tuple to a scalar energy for putting in an EnergyMatrix
        DenseFeatureSet fs = featSet.featMatrix.getTupleValue(tup);
        if(fs==null)//pruned
            return Double.POSITIVE_INFINITY;
        if(fs instanceof DenseFeatureSet.Const)
            return coeffs[featSet.featOffsets.getTupleValue(tup)];
        else
            throw new RuntimeException("ERROR: Can't convert non-const tuple contribution to scalar energy");
    }

    
    public Iterable<DegreeOfFreedom> getContDOFs() {
        return featSet.confSpace.confDOFs;//featSet.contDOFs;
    }  

    public int getNumParams() {
        return featSet.numFeatures;
    }

    public MoleculeModifierAndScorer makeMMS(ConfSample samp) {
        return featSet.makeMMS(samp);
    }
    
    
}
