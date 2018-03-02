/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.xminfit;

import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import java.util.ArrayList;
import java.util.Random;

/**
 * This version does not attempt to do gradients; rather it uses a bunch of samples
 * to tabulate the energy sum as a function of one degree of freedom
 * (given the others are as specified by the expansion)
 * and then updates the expansion to get this right
 * @author mhall44
 */
public class BatchCCDXMinFitter {
    
    final static int numGridpts = 20;//how many gridpts to use for each DOF
    //fine resolution seems pretty important.  1CC8 conf [5,7,7,5,0,7,4] is
    //0.2 over the CCD min with 20 gridpts, sub-0.1 when go to 40 gridpts.  
    //For [5,7,10,5,0,7,4] even that is not enough...around 0.6 at 40, 
    //drop to sub-0.1 with 300 gridpts
    
    public static void doXMinFit(ArrayList<SampleXMin> samples, XMFMatrix XMinMatrix){
        double oldESum, ESum=SampleXMin.curESum(samples,XMinMatrix);
        final double ESumConvTol = 0.01;//DEBUG!!
        
        System.out.println("Doing XMinFit.  Starting ESum: "+ESum);
        
        double gridwidths[] = new double[XMinMatrix.numCATSDOFs];
        double gridbase[] = new double[XMinMatrix.numCATSDOFs];
        buildGrid(gridwidths, gridbase, samples.get(0));
                
        do {
            oldESum = ESum;

            for(int catsDOFNum=0; catsDOFNum<XMinMatrix.numCATSDOFs; catsDOFNum++){
                
                double gridwidth = gridwidths[catsDOFNum];
                double[][] grids = sampEGrids(samples,XMinMatrix,catsDOFNum,gridbase,gridwidths);//energy at this dof's gridpts
                //Due to sparsity of feaches, lets just use the whole set
                
                int[] curGridpts = curSampleGridpts(samples,XMinMatrix,catsDOFNum,gridbase,gridwidths);
                //which gridpts each sampl is currently at
                
                int dofNumIter = 0;
                while (true){//now start updating tuples' parameters until convergence
                    //always update by one grid width   
                    boolean updatedTup = false;//tracking convergence
                    for(RCTuple tup : XMinMatrix.allTuplesForFit()){
                        double curTupESum = tupESum(tup, samples, grids, curGridpts, 0);
                        double upTupESum = tupESum(tup, samples, grids, curGridpts, 1);
                        double downTupESum = tupESum(tup, samples, grids, curGridpts, -1);
                        
                        if(curTupESum > Math.min(upTupESum,downTupESum)){
                            updatedTup = true;
                            //can improve total energy sum by adjusting this tuple's parameter
                            if(upTupESum<downTupESum){//adjust up
                                adjustMatrixValue(XMinMatrix, tup, catsDOFNum, gridwidth);
                                for(int sampNum : sampNumsForTup(tup,samples))
                                    curGridpts[sampNum]++;
                            }
                            else {//adjust down
                                adjustMatrixValue(XMinMatrix, tup, catsDOFNum, -gridwidth);
                                for(int sampNum : sampNumsForTup(tup,samples))
                                    curGridpts[sampNum]--;
                            }
                        }
                    }
                    
                    if(!updatedTup){//converged
                        System.out.println("CATS DOF "+catsDOFNum+" went through "+dofNumIter+" iterations");
                        break;
                    }
                    else
                        dofNumIter++;
                }
            }
            
            ESum = SampleXMin.curESum(samples,XMinMatrix);
            System.out.println("E sum: "+ESum);
        } while(ESum < oldESum-ESumConvTol);//convergence of overall fit
        
    }
    
    static double[][] sampEGrids(ArrayList<SampleXMin> samples, XMFMatrix XMinMatrix, int catsDOFNum, double gridbase[], double gridwidths[]){
        double[][] ans = new double[samples.size()][numGridpts];
        for(int s=0; s<samples.size(); s++){
            MoleculeModifierAndScorer mms = samples.get(s).makeObjFcn(null);
            mms.setDOFs( SampleXMin.minDOFVals(samples.get(s).makeObjFcn(XMinMatrix)).dofValues );//put rest of bb where matrix indicates, sidechain at best place given bb
            ArrayList<Integer> catsDOFIndices = SampleXMin.catsDOFIndices(mms);
            for(int g=0; g<numGridpts; g++){
                ans[s][g] = mms.getValForDOF(catsDOFIndices.get(catsDOFNum), gridbase[catsDOFNum]+g*gridwidths[catsDOFNum]);
            }
        }
        
        return ans;
    }
    
    private static int[] curSampleGridpts(ArrayList<SampleXMin> samples, XMFMatrix XMinMatrix, int catsDOFNum, double gridbase[], double gridwidths[]){
        //what gridpt is each of the sample currently at (based on XMinMatrix)?
        int ans[] = new int[samples.size()];
        for(int s=0; s<samples.size(); s++){
            DoubleMatrix1D xmin = XMinMatrix.getXMin(samples.get(s).conf);
            ans[s] = Math.round( (float) ((xmin.get(catsDOFNum)-gridbase[catsDOFNum])/gridwidths[catsDOFNum]) );
        }
        return ans;
    }
    
    
    
    private static void adjustMatrixValue(XMFMatrix XMinMatrix, RCTuple tup, 
            int catsDOFNum, double adjustment){
        DoubleMatrix1D tupVec = XMinMatrix.getTupleValue(tup);
        tupVec.set(catsDOFNum, tupVec.get(catsDOFNum)+adjustment);
    }
    
    
    private static double tupESum(RCTuple tup, ArrayList<SampleXMin> samples, 
            double[][] grids, int[] curGridpts, int offset){
        //using the grids, calculate the sum of energies for samples involved the given tuple
        double ans = 0;
        for(int s : sampNumsForTup(tup,samples)){
            int gridpt = curGridpts[s] + offset;
            //make sure gridpt is in range
            if(gridpt<0)//have a small adjustment to favor staying in range if it's a wash otherwise
                ans += grids[s][0] + 0.001;
            else if(gridpt>=grids[s].length)
                ans += grids[s][grids[s].length-1] + 0.001;
            else
                ans += grids[s][gridpt];
        }
        return ans;
    }
    
    
    private static ArrayList<Integer> sampNumsForTup(RCTuple tup, ArrayList<SampleXMin> samples){
        //DEBUG!!! This may be slow, will have to cache somehow if this is a bottleneck...
        //(possibly cache a hash set of sample indices with each rc, can check each rc in the tup against this)
        //(which takes up asymptotically the same amount of space as storign teh conf for each sample
        ArrayList<Integer> ans = new ArrayList<>();
        for(int sampNum=0; sampNum<samples.size(); sampNum++){
            boolean sampleMatchesTuple = true;
            for(int q=0; q<tup.size(); q++){
                if(samples.get(sampNum).conf[tup.pos.get(q)] != tup.RCs.get(q)){
                    sampleMatchesTuple = false;
                    break;
                }
            }
            
            if(sampleMatchesTuple)
                ans.add(sampNum);
        }
        
        return ans;
    }
    
    
    
    
    private static void buildGrid(double[] gridwidths, double[] gridbase, SampleXMin samp){
        //make a grid for each CATS degree of freedom, specified by the spacing (gridwidths) and the first DOF value (gridbase)
        MoleculeModifierAndScorer mms = samp.makeObjFcn(null);
        ArrayList<Integer> catsDOFIndices = SampleXMin.catsDOFIndices(mms);
        DoubleMatrix1D[] constr = mms.getConstraints();
        
        for(int d=0; d<catsDOFIndices.size(); d++){
            gridwidths[d] = ( constr[1].get(catsDOFIndices.get(d)) - constr[0].get(catsDOFIndices.get(d)) ) / (numGridpts-1);
            gridbase[d] = constr[0].get(catsDOFIndices.get(d));
        }
    }
    
    
    
}
