/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ematrix.epic;

import java.util.ArrayList;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.minimization.PLUGEnhancedMinimizer;
import edu.duke.cs.osprey.minimization.SQPMinimizer;
import edu.duke.cs.osprey.plug.PolytopeMatrix;

/**
 *
 * A matrix of 1-body, pairwise, etc. energies as polynomials in internal coordinates (EPIC)
 * Together with the EnergyMatrix, this allows precomputation of the complete energy of a continuously flexible system.  
 * 
 * @author mhall44
 */
public class EPICMatrix extends TupleMatrixGeneric<EPoly> {
    
    ConfSpace confSpace = null;//the conformational space for which the energy polynomials are defined
    //we need this to be able to get conformational energies from the polynomials

    
    public EPICMatrix(ConfSpace confSpace, double pruningInterval) {
        super(confSpace.numPos, confSpace.getNumRCsAtPos(),pruningInterval,null);
        this.confSpace = confSpace;
    }
    
    
    public double minimizeEnergy(RCTuple RCTup, boolean useSQP, boolean includeMinE,
            PolytopeMatrix plugMat, int numConsistent){
        //Minimize EPIC energy for the tuple
        //return infinity if the tuple is pruned
        //If plugMat is not null then use the PLUGEnhancedMinimizer
        
       
        if( useSQP && plugMat!=null )
            throw new RuntimeException("ERROR: PLUG-enhanced SQP minimizer not currently available");
        if( numConsistent!=1 && plugMat==null )
            throw new RuntimeException("ERROR: numConsistent is for PLUGEnhancedMinimizer");
        
        //137DEBUG!!!  seeing if min taking too long
        long startTime = System.currentTimeMillis();
        
        EPICEnergyFunction efunc = internalEnergyFunction(RCTup,includeMinE);
        
        //137DEBUG!!
        long gotEFcnTime = System.currentTimeMillis()-startTime;
        
        if(efunc==null)//pruned!
            return Double.POSITIVE_INFINITY;
        
        MoleculeModifierAndScorer objFcn = new MoleculeModifierAndScorer(efunc,confSpace,RCTup);
        
        Minimizer minim;
        if(useSQP)
            minim = new SQPMinimizer(objFcn,null);
        else if(plugMat!=null)
            minim = new PLUGEnhancedMinimizer(objFcn,plugMat,RCTup,numConsistent);
        else
            minim = new CCDMinimizer(objFcn,false);
        
        
        //137DEBUG!!
        long preMinTime = System.currentTimeMillis()-startTime;
        
        
        DoubleMatrix1D bestDOFVals = minim.minimize().dofValues;
        double E = objFcn.getValue(bestDOFVals);
        
        //137DEBUG!!!
        long sampTime = System.currentTimeMillis() - startTime;
        //taking over 10 s is going to be an issue
        if(sampTime > 10000){
            System.out.println();
            System.out.println("Minimization took over 10 s (ms shown): "+sampTime);
            System.out.println("Time to get E Fcn: "+gotEFcnTime+" to start min: "+preMinTime);
            System.out.println("Sample: "+RCTup.stringListing());
            System.out.println("Energy: "+E);
            if(minim instanceof CCDMinimizer){
                CCDMinimizer ccdMin = (CCDMinimizer)minim;

                System.out.println("GVCounts.  Estmin: "+ccdMin.GVCountEstmin+" Edge: "+ccdMin.GVCountEdge
                        +" Bigger: "+ccdMin.GVCountBigger+" Smaller: "+ccdMin.GVCountSmaller);
            }
            
            System.out.println();
        }
                
        return E;
    }
    
    
    
    
    /*public double minContE(RCTuple RCTup, boolean useSQP){
        //compute the minimized continuous component of the energy for the specified tuple
        //(i.e., conformationally minimized sum of EPIC term values, without minE)
        //return infinity if the tuple is pruned
       
       
            //137DEBUG!!!  seeing if min taking too long
            long startTime = System.currentTimeMillis();
       
       
       
        EPICEnergyFunction efunc = internalEnergyFunction(RCTup);
        
        //137DEBUG!!
        long gotEFcnTime = System.currentTimeMillis()-startTime;
        
        
        if(efunc==null)//pruned!
            return Double.POSITIVE_INFINITY;
        
        MoleculeModifierAndScorer objFcn = new MoleculeModifierAndScorer(efunc,confSpace,RCTup);
        
        Minimizer minim;
        if(useSQP)
            minim = new SQPMinimizer(objFcn,null);
        else
            minim = new CCDMinimizer(objFcn,false);
        
        
        //137DEBUG!!
        long preMinTime = System.currentTimeMillis()-startTime;
        
        
        
        DoubleMatrix1D bestDOFVals = minim.minimize().dofValues;
        double E = objFcn.getValue(bestDOFVals);
        
        //DEBUG!!!!
        /*Minimizer ccdMin = new CCDMinimizer(objFcn,false);
        DoubleMatrix1D ccdDOFVals = ccdMin.minimize().dofValues;
        System.out.println("Energy track");
        for(int a=0; a<=20; a++){
            double frac = a/20.;
            DoubleMatrix1D x = bestDOFVals.copy().assign(Functions.mult(1-frac)).assign(ccdDOFVals,Functions.plusMult(frac));
            System.out.println(objFcn.getValue(x));
        }
        System.out.println("Energy track done");*/
        
        //DEBUG!!!!!!!!
        /*String objFcnFile = "OBJFCN"+System.currentTimeMillis()+".dat";
        ObjectIO.writeObject(objFcn, objFcnFile);
        System.out.println("OUTPUTTED OBJ FCN TO "+objFcnFile+".  E: "+E+" bestDOFVals: "+bestDOFVals);
        */
        
        //137DEBUG!!!
       /* long sampTime = System.currentTimeMillis() - startTime;
        //taking over 10 s is going to be an issue
        if(sampTime > 10000){
            System.out.println();
            System.out.println("Minimization took over 10 s (ms shown): "+sampTime);
            System.out.println("Time to get E Fcn: "+gotEFcnTime+" to start min: "+preMinTime);
            System.out.println("Sample: "+RCTup.stringListing());
            System.out.println("Energy: "+E);
            if(minim instanceof CCDMinimizer){
                CCDMinimizer ccdMin = (CCDMinimizer)minim;

                System.out.println("GVCounts.  Estmin: "+ccdMin.GVCountEstmin+" Edge: "+ccdMin.GVCountEdge
                        +" Bigger: "+ccdMin.GVCountBigger+" Smaller: "+ccdMin.GVCountSmaller);
            }
            
            System.out.println();
        }
        
        
        
        return E;
    }
   
    public double minContE(int[] conf){
        return minContE(conf,false);//CCD by default for now
    }
   
    public double minContE(int[] conf, boolean useSQP){
        ////Given a list of RCs for each position, compute the minimized cont component of the energy
        //negative RC numbers indicate undefined positions 
        //(exclude these positions: EPIC terms always positive, so this gives
        //a valid lower bound).  again, return infinity if conf pruned
        return minContE(new RCTuple(conf), useSQP);
    }*/
    
    
    public EPICEnergyFunction internalEnergyFunction(RCTuple tup, boolean includeMinE){
        //Make an energy function representing the internal energy of the tuple
        //kind of an EPIC analog of EnergyMatrix.getInternalEnergy
        //return null if the tuple is pruned
        
        int numPosInTuple = tup.pos.size();
        
        ArrayList<EPoly> terms = new ArrayList<>();
        
        for(int indexInTuple=0; indexInTuple<numPosInTuple; indexInTuple++){
            int posNum = tup.pos.get(indexInTuple);
            int RCNum = tup.RCs.get(indexInTuple);
            
            EPoly intraE = getOneBody(posNum,RCNum);
            if(intraE==null)//pruned
                return null;
            
            terms.add( intraE );
            
            for(int index2=0; index2<indexInTuple; index2++){
                int pos2 = tup.pos.get(index2);
                int rc2 = tup.RCs.get(index2);
                
                EPoly pairwiseE = getPairwise(posNum,RCNum,pos2,rc2);
                if(pairwiseE==null)//pruned
                    return null;
                
                terms.add( pairwiseE );
                
                //higher-order energies will go here
                //though currently not supporting these for EPIC
                if( getHigherOrderTerms(posNum,RCNum,pos2,rc2) != null ){
                    throw new UnsupportedOperationException("ERROR: Higher-order EPIC terms "
                            + "not yet supported in internalEnergyFunction");
                }
            }
        }
        
        EPICEnergyFunction efunc = new EPICEnergyFunction(terms, includeMinE);
        return efunc;
    }

    public ConfSpace getConfSpace() {
        return confSpace;
    }
    
    
    
    
}
