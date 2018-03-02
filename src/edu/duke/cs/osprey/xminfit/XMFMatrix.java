/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.xminfit;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import java.util.ArrayList;

/**
 *
 * Matrix of contributions to the minimized CATS coordinates from rotamer pairs
 * 
 * @author mhall44
 */
public class XMFMatrix extends TupleMatrixGeneric<DoubleMatrix1D> {
    
    int numCATSDOFs;
    
    XMFMatrix(ConfSpace cSpace){
        super(cSpace, Double.POSITIVE_INFINITY, null);
        numCATSDOFs=0;
        for(DegreeOfFreedom dof : cSpace.confDOFs){
            if(dof instanceof BBFreeDOF)
                numCATSDOFs++;
        }
    }
    
    DoubleMatrix1D getXMin(int[] conf){
        //add up the tuple contributions
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(numCATSDOFs); 
        for(RCTuple tup : tuplesForSamp(conf)){
            ans.assign(getTupleValue(tup), Functions.plus);
        }
        return ans;
    }
    
    ArrayList<RCTuple> tuplesForSamp(int[] conf){
        //which tuples in XMinMatrix are active for samp?
        ArrayList<RCTuple> ans = new ArrayList<>();
        for(int pos=0; pos<conf.length; pos++){
            if(getOneBody(pos,conf[pos])!=null)
                ans.add(new RCTuple(pos,conf[pos]));
            for(int pos2=0; pos2<pos; pos2++){
                if(getPairwise(pos,conf[pos],pos2,conf[pos2])!=null)
                    ans.add(new RCTuple(pos,conf[pos],pos2,conf[pos2]));
            }
        }
        return ans;
    }
    
    
    ArrayList<RCTuple> allTuplesForFit(){
        //which tuples in XMinMatrix are active for samp?
        ArrayList<RCTuple> ans = new ArrayList<>();
        for(int pos=0; pos<getNumPos(); pos++){
            for(int rc=0; rc<getNumConfAtPos(pos); rc++){
                if(getOneBody(pos,rc)!=null)
                    ans.add(new RCTuple(pos,rc));
                for(int pos2=0; pos2<pos; pos2++){
                    for(int rc2=0; rc2<getNumConfAtPos(pos2); rc2++){
                        if(getPairwise(pos,rc,pos2,rc2)!=null)
                            ans.add(new RCTuple(pos,rc,pos2,rc2));
                    }
                }
            }
        }
        return ans;
    }
    
    //initializing with 0
    public void initOneBody(int pos, int rc){
        setOneBody(pos,rc,DoubleFactory1D.dense.make(numCATSDOFs));
    }
    
    public void initPairwise(int pos, int rc, int pos2, int rc2){
        setPairwise(pos,rc,pos2,rc2,DoubleFactory1D.dense.make(numCATSDOFs));
    }
    
}
