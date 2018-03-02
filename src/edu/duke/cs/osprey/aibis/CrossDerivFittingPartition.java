/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.aibis;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.energy.derivs.EnergyDerivs;
import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.RealLinearOperator;

/**
 *
 * Fit a set of >=quadratic parameters (probably those associated with a particular res pair or single for energy term reasons)
 * 
 * In quadratic fit this should grab all the poly's that go into cross terms but still need SAPEs on right hand side!  And integrated crossderivs!!
 * 
 * 
 * hmm so actually the fits themselves have no benefit in being combined for different derivatives (disjoint parameter sets), at least for quadratics
 * but it's most efficient to get SAPE cross derivs together
 * so maybe get those by batch and create CrossDerivFittingPartition by batches??  
 * 
 * This doesn't apply to the other types of fitting partitions though
 * 
 * @author mhall44
 */


/*
public class CrossDerivFittingPartition extends FittingPartition {
    //DEBUG!!  Assumes only quadratic cross-terms and SAPE terms
    //have nonzero cross-2nd-derivs
    
    int dof1,dof2;
    //Only on type of continuous term, several RCTuples, so use regular ConfSamples
    ArrayList<ConfSample> samples;
    double[] sampleCrossDerivs;//cross-derivatives of the (energy-SAPE energy) corresponding to the samples
    //TO GENERATE THESE: draw samples appropriate for RCTuple with the
    //desired cross derivs, make an EnergyDerivs for the relevant energy term
    //probably only worth calculating cross-derivs that are not always within a term together. 
    
    int numRCsAtPos[];
    
    
    
    private static class RelevantPosTuple {
        //a tuple of positions for which the cross-derivative appears
        int[] pos;//the positions
        private int[] paramIndices;//Indices of RC tuples' parameters within the fit parameters
        int numRCsAtPos[];
        
        RelevantPosTuple(RCTuple tup, int numRCsAtPos[]){
            //create from a single tuple
            this.numRCsAtPos = numRCsAtPos;
            pos = new int[tup.size()];
            int paramIndicesSize = 1;
            for(int q=0; q<tup.size(); q++){
                pos[q] = tup.pos.get(q);
                paramIndicesSize *= numRCsAtPos[pos[q]];
            }
            paramIndices = new int[paramIndicesSize];
            Arrays.fill(paramIndices, -1);
            paramIndices[getIndexForRCs(tup.RCs)] = 0;
        }
        
        
        //constructor starts us off with no fit parameters
        RelevantPosTuple(int[] pos, int numRCsAtPos[]){
            this.numRCsAtPos = numRCsAtPos;
            this.pos = pos;
            int paramIndicesSize = 1;
            for(int p : pos)
                paramIndicesSize *= numRCsAtPos[p];
            paramIndices = new int[paramIndicesSize];
            Arrays.fill(paramIndices, -1);
        }
        
        int getParamIndex(int[] sample){
            //lookup in paramIndices, based on the RC assignments to the pos in sample
            int ans = 0;
            for(int a=0; a<pos.length; a++){
                if(a>0)
                    ans *= numRCsAtPos[pos[a-1]];
                ans += sample[pos[a]];
            }
            return paramIndices[ans];
        }
        
        final int getIndexForRCs(ArrayList<Integer> RCs){
            //index in paramIndices, based on the RC's for pos
            int ans = 0;
            for(int a=0; a<pos.length; a++){
                if(a>0)
                    ans *= numRCsAtPos[pos[a-1]];
                ans += RCs.get(a);
            }
            return ans;
        }
        
        RCTuple rcTupleAtIndex(int index){
            //RCTuple for an index in paramIndices
            ArrayList<Integer> RCs = new ArrayList<>();
            //we will get the RCs in reverse order
            for(int a=pos.length-1; a>=0; a--){//we actually get the RCs in reverse order
                RCs.add(0, index%pos[a]);
                index /= pos[a];
            }
            
            return new RCTuple(new ArrayList(Arrays.asList(pos)), RCs);
        }
    }
    
    
    
    
    
    
    ArrayList<RelevantPosTuple> relevantPosTuples;
    
    
    public CrossDerivFittingPartition(FeatureSet featSet, ArrayList<ConfSample> samples,
            double[] sampleCrossDerivs, int dof1, int dof2, int posWithCrossDerivs[][], 
            ArrayList<RelevantPosTuple> relevantPosTuples){
        this.featSet = featSet;
        this.samples = samples;
        this.sampleCrossDerivs = sampleCrossDerivs;
        this.relevantPosTuples = relevantPosTuples;
        numSamples = samples.size();
        this.dof1 = dof1;
        this.dof2 = dof2;
        numRCsAtPos = featSet.confSpace.getNumRCsAtPos();
        
        computeRefitParamIndices();
        numParams = refitParamIndices.size();
             
        prepareOperators();
    }
    
    
    public static ArrayList<CrossDerivFittingPartition> makeSidechainPairCrossDerivs(FeatureSet featSet,
            ConfSampleSet trainingSamples, int pos1, int pos2){
        //make the cross derivatives between a pair of sidechains's dihedrals
        //these are the only cross-derivatives we need in a rigid-backbone calculation
        int posWithCrossDerivs[][] = new int[][] {new int[] {pos1,pos2}};
        
        //hmm so in rigid-backbone case we only need single-term cross derivs.  
        //OK really need to make a choice about how to do cross derivs
        
        ArrayList<CrossDerivFittingPartition> ans = new ArrayList<>();
        int numRCsAtPos[] = featSet.confSpace.getNumRCsAtPos();
        
        for(int rc1=0; rc1<numRCsAtPos[pos1]; rc1++){
            for(int rc2=0; rc2<numRCsAtPos[pos2]; rc2++){
                RCTuple pair = new RCTuple(pos1,rc1,pos2,rc2);
                if( featSet.featMatrix.getPairwise(pos1,rc1,pos2,rc2) != null ){
                    ArrayList<ConfSample> samples = trainingSamples.getSampsForTup(pair, 40);//DEBUG!!  40 samples
                    EnergyDerivs[] ed = new EnergyDerivs[samples.size()];
                    for(int s=0; s<samples.size(); s++)
                        ed[s] = makeEnergyDerivs(samples.get(s), pair);//integrated - SAPE

                    
                    int termDOFs[] = edDOFINDICES;
                    int[][] dofPartitions = partitionDOFsForPair(pair);
                    //cross-terms will be between DOF indices in dofPartitions[0]
                    //and those in dofPartitions[1] (indexed among dofs in ed)
                    for(int dofIndex1 : dofPartitions[0]){
                        for(int dofIndex2 : dofPartitions[1]){
                            double sampleCrossDerivs[] = new double[samples.size()];
                            for(int s=0; s<samples.size(); s++)
                                sampleCrossDerivs[s] = ed[s].getHess().get(dofIndex1, dofIndex2);
                            
                            ans.add( new CrossDerivFittingPartition(featSet, samples, sampleCrossDerivs,
                                    termDOFs[dofIndex1], termDOFs[dofIndex2], posWithCrossDerivs, 
                                    Arrays.asList(new RelevantPosTuple(pair,numRCsAtPos))) );
                        }
                    }
                }
            }
        }
        
        return ans;
    }
    
    
    
    
    private void makeEnergyDerivs(ConfSample samp){
        EnergyDerivs ed = new EnergyDerivs(mms, x);
    }
    
    
    public static int[][] partitionDOFsForPair(RCTuple pair, etc){
        
    }
    
    
    
    
    //how to represent parameters?
    //there is one param per tuple, only for a small subset of tuples,
    //although fairly dense within relevant position pairs
    //For a sample, I need to look up what tuples are involved and what their param indices are
    //(local to this parition, and alsoglobal)
    
    //so have a list of pos/pairs, for each have an arr indexed by RC index (int) tuple -> param index or whatever
    
    //MAY WANT TO HAVE SPECIAL SOLVER IF OPERATOR TURNS OUT TO BE 1X1
    
    final void prepareOperators() {
        OpenMapRealMatrix directOp = new OpenMapRealMatrix(numSamples,numParams);
        for(int s=0; s<numSamples; s++){
            for(RelevantPosTuple rpt : relevantPosTuples){
                int paramIndex = rpt.getParamIndex(samples.get(s).vox);
                if(paramIndex!=-1)
                    directOp.setEntry(s, paramIndex, 1);//d^2(xy)/dxdy = 1
            }
        }
            
        directOperator = directOp;
        transposeOperator = (RealLinearOperator) directOp.transpose();
    }

    @Override
    double[] computeTargetVals(double[] coeffs) {
        return sampleCrossDerivs;
    }

    private ArrayList<Integer> computeRefitParamIndices() {
        ArrayList<Integer> ans = new ArrayList<>();
        for(RelevantPosTuple rpt : relevantPosTuples){
            for(int rcTupIndex=0; rcTupIndex<rpt.paramIndices.length; rcTupIndex++){
                if(rpt.paramIndices[rcTupIndex] != -1){//there's a fittable parameter for this RCTuple
                    RCTuple tup = rpt.rcTupleAtIndex(rcTupIndex);
                    //the parameter we want is the cross terms for dof1, dof2 in the feature for tup
                    int globalOffset = featSet.featOffsets.getTupleValue(tup);
                    DenseFeatureSet fs = featSet.featMatrix.getTupleValue(tup);
                    if(fs instanceof DenseFeatureSet.EPIC){
                        int indexInPoly = ((DenseFeatureSet.EPIC)fs).monomialIndexByDOFs(dof1,dof2);
                        ans.add( globalOffset + indexInPoly );
                    }
                }
            }
        }
        
        return ans;
    }
    
}

*/