/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ibis;

import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tupexp.BasicPruningTupleExpander;

/**
 *
 * @author mhall44
 */
public class VoxDrawTupleExpander extends BasicPruningTupleExpander {//for drawing discrete samples

        FeatureSet featSet;
        int numSampPerFeach = 30;//DEBUG!!!
        
        VoxDrawTupleExpander(PruningMatrix pruneMat, PolytopeMatrix plugMat, FeatureSet featSet){
            super(pruneMat,plugMat);
            this.featSet = featSet;
        }
        
        @Override
        public int numSamplesNeeded(int tup){
            //for a tuple (index in tuples), how many samples are needed?
            return numSampPerFeach*featSet.featMatrix.getTupleValue(getTuples().get(tup)).getNumFeatures();
        }
}
