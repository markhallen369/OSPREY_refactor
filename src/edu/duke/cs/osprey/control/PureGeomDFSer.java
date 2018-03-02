/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.plug.FunnyVoxBuilder;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.plug.PureGeomDFS;
import edu.duke.cs.osprey.plug.triage.PLUGTriage;

/**
 *
 * Purely geometric depth-first searcher
 * 
 * @author mhall44
 */
public class PureGeomDFSer {
    
    SearchProblem sp;
    int bbVoxSplitNum;
    double resPoseVoxWidth;

    public PureGeomDFSer(ConfigFileParser cfp) {
        sp = cfp.getSearchProblem();
        bbVoxSplitNum = cfp.params.getInt("BBVoxSplitNum", 1);
        resPoseVoxWidth = cfp.params.getDouble("ResPoseVoxWidth", 1);
    }
    
    
    
    public void doSearch(){
        
        if(bbVoxSplitNum>1){//split the big backbone voxel
            System.out.println("Splitting into funny voxels");
            FunnyVoxBuilder fvb = new FunnyVoxBuilder(sp, bbVoxSplitNum, resPoseVoxWidth);
            fvb.buildFunnyVox();
        }
        
        //DEBUG!! This may be only necessary for fvb...
        PLUGTriage pt = new PLUGTriage(sp);
        pt.doTriage();
        
        PolytopeMatrix mat = new PolytopeMatrix(sp, false);
        //consider pruning?  Could help with DFS efficiency...
        PureGeomDFS pgd = new PureGeomDFS(mat);
        RCTuple conf = pgd.findConf();
        System.out.println("FOUND CONF BY PURE GEOM DFS: "+conf.stringListing());
    }
    
}
