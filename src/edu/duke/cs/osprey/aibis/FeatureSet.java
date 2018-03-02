/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.aibis;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPoly;
import edu.duke.cs.osprey.ematrix.epic.SAPE;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;


/**
 *
 * @author mhall44
 */
public class FeatureSet implements Serializable {
    
    TupleMatrixGeneric<DenseFeatureSet> featMatrix;
    TupleMatrixGeneric<Integer> featOffsets;//where in the feature vector do the features for each tuple start?
    int numFeatures;
    ConfSpace confSpace;
    
    boolean freezeSAPE = true;//Don't scale SAPE terms by coeffs
    boolean hasSAPE;
    boolean hasRedundantFeatures = false;
    
    HashMap<String,Integer> dofIndexLookup;
        
    public FeatureSet(ConfSpace confSpace){
        this.confSpace = confSpace;
        featMatrix = new TupleMatrixGeneric<>(confSpace, Double.POSITIVE_INFINITY, null);
        //if we consider our clash definitions to be solid we effectively have pruned for infinite pruning interval
        featOffsets = new TupleMatrixGeneric<>(confSpace, Double.POSITIVE_INFINITY, -1);
        numFeatures = 0;
        hasSAPE = false;
        
        dofIndexLookup = new HashMap<>();
        for(int dofIndex=0; dofIndex<confSpace.confDOFs.size(); dofIndex++)
            dofIndexLookup.put(confSpace.confDOFs.get(dofIndex).getName(), dofIndex);
    }
    
    
    final void addDenseFeatureSet(DenseFeatureSet s, RCTuple tup){
        featMatrix.setTupleValue(tup, s);
        featOffsets.setTupleValue(tup, numFeatures);//add to end of feature vector
        numFeatures += s.getNumFeatures();
        //contDOFs.addAll(s.getContDOFs());
        if(extractSAPE(s)!=null){
            hasSAPE = true;
            if(freezeSAPE)//SAPE feature shouldn't be counted
                numFeatures--;
        }
    }
    
    final void addRedundantDenseFeatureSet(DenseFeatureSet s, RCTuple tup, int offset){
        //Parameter shared w/ previous term (at specified offset), so not really a new feature
        featMatrix.setTupleValue(tup, s);
        featOffsets.setTupleValue(tup, offset);
        if(extractSAPE(s)!=null)
            hasSAPE = true;
        
        hasRedundantFeatures = true;
    }
    
    double evalEnergy(ConfSample samp, double[] coeffs){
        return evalEnergy(samp, coeffs, makeMMS(samp));
    }
        
        
    //let Aij denote feature value j for sample i
    double evalEnergy(ConfSample samp, double[] coeffs, MoleculeModifierAndScorer mms){
        //Given a set of coeffs, evaluate energy (sum_j Aij coeffs_j)
        if(hasSAPE)
            applyContDOFs(samp.contDOFValsArr, mms);
        double E = 0;
        
        for(RCTuple tup : tuplesForSamp(samp)){
            Integer offset = featOffsets.getTupleValue(tup);
            if(offset==null)//sample is pruned
                return Double.POSITIVE_INFINITY;
            else
                E += featMatrix.getTupleValue(tup).evalEnergy(samp, coeffs, offset, freezeSAPE);
        }
        return E;
    }
    
    
    
    void applyContDOFs(double[] allDOFVals, MoleculeModifierAndScorer mms){
        for(int dofNum=0; dofNum<mms.getNumDOFs(); dofNum++){
            //DEBUG!!! may benefit from settings DOFs at once
            //Look into blocks for MMS!!!
            String dofName = mms.getDOFs().get(dofNum).getName();
            int dofIndex = dofIndexLookup.get(dofName);//index in full DOF list
            mms.setDOF(dofNum, allDOFVals[dofIndex]);
            //if(contDOFVals.containsKey(dof.getName()))//we have a value for this dof
            //    mms.setDOF(dofNum, contDOFVals.get(dof.getName()));
        }
    }
    
    /*double evalNoSAPEEnergy(ConfSample samp, double[] coeffs){
        //Given a set of coeffs, evaluate energy (sum_j Aij coeffs_j)
        double E = 0;
        
        for(RCTuple tup : tuplesForSamp(samp)){
            Integer offset = featOffsets.getTupleValue(tup);
            if(offset==null)//sample is pruned
                return Double.POSITIVE_INFINITY;
            else
                E += featMatrix.getTupleValue(tup).evalNoSAPEEnergy(samp, coeffs, offset);
        }
        return E;
    }*/
    //Trying to avoid tuplesForSamp bottleneck (making RCTuples)
    double evalNoSAPEEnergy(ConfSample samp, double[] coeffs){
        //Given a set of coeffs, evaluate energy (sum_j Aij coeffs_j)
        double E = 0;
        
        for(int pos=0; pos<featMatrix.getNumPos(); pos++){
            int rc = samp.getRC(pos);//DEBUG!!!  Assumes full-conf samp
            Integer offset = featOffsets.getOneBody(pos,rc);
            if(offset==null)//sample is pruned
                return Double.POSITIVE_INFINITY;
            E += featMatrix.getOneBody(pos,rc).evalNoSAPEEnergy(samp, coeffs, offset);
            for(int pos2=0; pos2<pos; pos2++){
                int rc2 = samp.getRC(pos2);
                offset = featOffsets.getPairwise(pos,rc,pos2,rc2);
                if(offset==null)//sample is pruned
                    return Double.POSITIVE_INFINITY;
                E += featMatrix.getPairwise(pos,rc,pos2,rc2).evalNoSAPEEnergy(samp, coeffs, offset);
            }
        }
        if(featMatrix.hasHigherOrderTerms())
            throw new RuntimeException("ERROR: higher-order terms not yet handled here");
        
        return E;
    }
    
    
    double evalFrozenEnergy(ConfSample samp, MoleculeModifierAndScorer mms){
        //part of energy not dependent on coeffs
        double E = 0;

        if(hasSAPE){
            applyContDOFs(samp.contDOFValsArr, mms);
            for(SAPE sape : listSAPETerms(samp))
                E += sape.getEnergySharedMolec();
        }

        return E;
    }
        
    
    void updateATProduct(ConfSample samp, double sampVal, double[] outVec){
        MoleculeModifierAndScorer mms = freezeSAPE ? null : makeMMS(samp);//only need this if optimizing SAPE coeff
        updateATProduct(samp, sampVal, outVec, mms);
    }
    
    /*void updateATProduct(ConfSample samp, double sampVal, double[] outVec, MoleculeModifierAndScorer mms){
        //outVec is becoming A^T inVec where inVec_i = sampVal.  
        //Deal with the constribution from samp
        //by doing outVec_j += A_ij sampVal for each j
        if(hasSAPE && !freezeSAPE)//if freezing SAPE then SAPE not needed here
            samp.applyContDOFs(mms);
        for(RCTuple tup : tuplesForSamp(samp)){
            int offset = featOffsets.getTupleValue(tup);
            featMatrix.getTupleValue(tup).updateATProduct(samp, sampVal, outVec, offset, freezeSAPE);
        }
    }*/
    
    //Version without tuplesForSamp
    void updateATProduct(ConfSample samp, double sampVal, double[] outVec, MoleculeModifierAndScorer mms){
        //outVec is becoming A^T inVec where inVec_i = sampVal.  
        //Deal with the constribution from samp
        //by doing outVec_j += A_ij sampVal for each j
        if(hasSAPE && !freezeSAPE)//if freezing SAPE then SAPE not needed here
            applyContDOFs(samp.contDOFValsArr, mms);
        
        for(int pos=0; pos<featMatrix.getNumPos(); pos++){
            int rc = samp.getRC(pos);//DEBUG!!!  Assumes full-conf samp
            Integer offset = featOffsets.getOneBody(pos,rc);
            featMatrix.getOneBody(pos,rc).updateATProduct(samp, sampVal, outVec, offset, freezeSAPE);
            for(int pos2=0; pos2<pos; pos2++){
                int rc2 = samp.getRC(pos2);
                offset = featOffsets.getPairwise(pos,rc,pos2,rc2);
                featMatrix.getPairwise(pos,rc,pos2,rc2).updateATProduct(samp, sampVal, outVec, offset, freezeSAPE);
            }
        }
        if(featMatrix.hasHigherOrderTerms())
            throw new RuntimeException("ERROR: higher-order terms not yet handled here");
    }
    
    
    
    
    
    MoleculeModifierAndScorer makeMMS(ConfSample samp){
        //Build a molecule modifier and scorer for the sample's voxel, applicable to this feature's SAPE terms
        if(!hasSAPE)
            return null;
        MoleculeModifierAndScorer mms = new MoleculeModifierAndScorer(null, confSpace, samp.vox);
        //Assign SAPE terms.  Note: Calling this function on another sample may invalidate this
        for(SAPE sape : listSAPETerms(samp))
            sape.assignSharedMolecule(mms.getMolec());
        return mms;
    }
    
    
    /*void initSample(ConfSample samp){
        //set up for computations involving a sample 
        MoleculeModifierAndScorer mms = makeMMS(samp);
        samp.applyContDOFs(mms);
    }*/
    
    
    ArrayList<RCTuple> unprunedFeatTuples(){
        //Enumerate the unpruned RCTuples associated with features in this FeatureSet
        ArrayList<RCTuple> ans = new ArrayList<>();
        for(int pos=0; pos<featMatrix.getNumPos(); pos++){
            for(int rc=0; rc<featMatrix.getNumConfAtPos(pos); rc++){
                if(featMatrix.getOneBody(pos,rc)!=null)
                    ans.add(new RCTuple(pos,rc));
                
                for(int pos2=0; pos2<pos; pos2++){
                    for(int rc2=0; rc2<featMatrix.getNumConfAtPos(pos2); rc2++){
                        if(featMatrix.getPairwise(pos,rc,pos2,rc2)!=null)
                            ans.add(new RCTuple(pos,rc,pos2,rc2));
                    }
                }
            }
        }
        
        if(featMatrix.hasHigherOrderTerms())
            throw new RuntimeException("ERROR: Higher-order terms not yet supported here");
        //DEBUG!! will want to fix this and possibly also avoid storing the whole list of RCTuples
        
        return ans;
    }
    
    
    ArrayList<RCTuple> tuplesForSamp(ConfSample samp){
        //Enumerate RCTuples with features that are a part of this samp
        ArrayList<RCTuple> ans = new ArrayList<>();
        for(int pos=0; pos<featMatrix.getNumPos(); pos++){
            int rc = samp.getRC(pos);//DEBUG!!!  Assumes full-conf samp
            ans.add(new RCTuple(pos,rc));
            for(int pos2=0; pos2<pos; pos2++){
                int rc2 = samp.getRC(pos2);
                ans.add(new RCTuple(pos,rc,pos2,rc2));
            }
        }
        if(featMatrix.hasHigherOrderTerms())
            throw new RuntimeException("ERROR: higher-order terms not yet handled here");
        return ans;
    }
    
    ArrayList<SAPE> listSAPETerms(ConfSample samp){
        //Get all the SAPE terms applicable to samp
        ArrayList<SAPE> ans = new ArrayList<>();
        for(RCTuple tup : tuplesForSamp(samp)){
            SAPE sape = extractSAPE(featMatrix.getTupleValue(tup));
            if(sape!=null)
                ans.add(sape);
        }
        return ans;
    }
    
    SAPE extractSAPE(DenseFeatureSet fs){//Extract the SAPE if there is one
        if(fs instanceof DenseFeatureSet.EPIC)
            return ((DenseFeatureSet.EPIC)fs).sapeTerm;
        else//Not an EPIC feature set, so no SAPE term.  Might even be null.  
            return null;
    }

    
    
}
