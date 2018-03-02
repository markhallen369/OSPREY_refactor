/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.xminfit;

import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.confspace.RCTuple;
import java.util.ArrayList;
import java.util.Random;

/**
 *
 * @author mhall44
 */
public class GradDescentXMinFitter {
    

    
    
    
    public static void doXMinFit(ArrayList<SampleXMin> samples, XMFMatrix XMinMatrix){
        double oldESum, ESum=SampleXMin.curESum(samples,XMinMatrix);
        final double ESumConvTol = 0.01;//DEBUG!!
        final int MINIBATCH_SIZE = samples.size();//DEBUG!!  Not really a minibatch I guess
        final double learningRate = 0.01/samples.size();//DEBUG!!
        
        Random rand = new Random();
        System.out.println("Doing XMinFit.  Starting ESum: "+ESum);
        
        do {
            oldESum = ESum;
            //maybe we can just try some SGD-type steps to get started.  
            for(int s=0; s<MINIBATCH_SIZE; s++){
                int randomSampNum = rand.nextInt(samples.size());
                SampleXMin samp = samples.get(randomSampNum);
                DoubleMatrix1D sampleCATSGrad = samp.curCATSGrad(XMinMatrix);
                //now do updates
                //dE_samp/dq = sum_samp dE_samp/dx_samp dx_samp/dq (dx_samp/dq = I(tup in samp))
                //q += learningRate * dE_samp/dq
                
                //DEBUG!!!  To avoid big gradients at clashes messing up everything
                double gradNormsq = sampleCATSGrad.zDotProduct(sampleCATSGrad);
                if(gradNormsq>10.){
                    sampleCATSGrad.assign(Functions.mult(1./Math.sqrt(gradNormsq)));
                }
                
                for(RCTuple tup : XMinMatrix.tuplesForSamp(samp.conf)){
                    DoubleMatrix1D q = XMinMatrix.getTupleValue(tup);
                    q.assign(sampleCATSGrad, Functions.plusMult(learningRate));
                }
            }
            
            ESum = SampleXMin.curESum(samples,XMinMatrix);
            System.out.println("E sum: "+ESum);
        } while(ESum < oldESum-ESumConvTol);//convergence of overall fit,
        //not of individual voxel minimizations (demands of fit could compromise those,
        //and we can test empirically if this is an issue)
        
        //if grad is not so good then can try a CCD-like approach
        //loop through all the CATS DOFs, for each go through all the samples with
        //the other DOFs frozen, build a spline or something of the energy on the voxel,
        //then minimize that :)
        
            /*
    Indeed this currently overshoots
    Doing XMinFit.  Starting ESum: 128081.52291967154
E sum: 129129.49096076888
XMinMatrix test.  
 Average signed energy difference (fullmin-matrix): -7.60514354654389
 Average unsigned energy difference (fullmin-matrix): 7.621729744282714
 Number of vox where xminfit is higher energy: 5540
 Number of vox where xminfit is lower energy: 11
 Average distance between full & matrix minima: 0.7715229651554191
    */
        
        //break this minim thing out into its own class anyway
    }
    
    
}
