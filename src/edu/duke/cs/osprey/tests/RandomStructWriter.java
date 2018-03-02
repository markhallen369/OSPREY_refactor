/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.pruning.Pruner;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tupexp.ConfETupleExpander;
import edu.duke.cs.osprey.tupexp.TupleExpander;
import java.util.ArrayList;

/**
 *
 * Given config for a GMEC design, sample some minimized conformations randomly
 * anything not sterically prunable is fair game
 * 
 * @author mhall44
 */
public class RandomStructWriter {
    
    SearchProblem searchSpace;
    int numStructs;
    double stericThresh;
    
    public RandomStructWriter(SearchProblem searchSpace, int numStructs, double stericThresh) {
        this.searchSpace = searchSpace;
        this.numStructs = numStructs;
        this.stericThresh = stericThresh;
    }
    
    
    //use args like for Main
    public static void main(String args[]){
        ConfigFileParser cfp = new ConfigFileParser(args);
        cfp.loadData();
        RandomStructWriter rsw = new RandomStructWriter(
                cfp.getSearchProblem(),
                cfp.getParams().getInt("NUMRANDSTRUCTS"),
                cfp.getParams().getDouble("STERICTHRESH")
        );
        rsw.writeStructs();
    }
    
    public void writeStructs(){
        searchSpace.loadEnergyMatrix();//The tuple enumerator wants this
        searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace, Double.POSITIVE_INFINITY);
        
        Pruner pruner = new Pruner(searchSpace, false, 0, 0, false, false);
        pruner.pruneSteric(stericThresh);
        
        //searchSpace.plugMat = new PolytopeMatrix(searchSpace,true);
        
        TupleExpander expander = new ConfETupleExpander(searchSpace, null);
        ArrayList<int[]> confs = expander.sampleRandomConfs(numStructs);
        for(int c=0; c<numStructs; c++){
            searchSpace.outputMinimizedStruct(confs.get(c), searchSpace.name + ".RANDSTRUCT"+c+".pdb");
        }
    }
    
}
