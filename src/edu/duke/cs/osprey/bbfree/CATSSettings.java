/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.bbfree;

import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * These will be used to set up BB free blocks
 * 
 * @author mhall44
 */
public class CATSSettings {
    ArrayList<String[]> freeBBZones = new ArrayList<>();
    boolean fixAnchors = true;
    boolean tooBigForSeries = false;
    double voxWidth = 0;//0 if to be automatically selected

    public CATSSettings(){
        
    }
    
    public CATSSettings(ArrayList<String[]> freeBBZones, boolean fixAnchors, boolean tooBigForSeries, double voxWidth) {
        this.freeBBZones = freeBBZones;
        this.fixAnchors = fixAnchors;
        this.tooBigForSeries = tooBigForSeries;
        this.voxWidth = voxWidth;
    }
    
    public ArrayList<BBFreeBlock> buildBBFreeBlocks(ArrayList<String> flexibleRes, Molecule m){
        //create a BFB for each (start res, end res) pair.  PDB residue numbers provided.  
        ArrayList<BBFreeBlock> ans = new ArrayList<>();
        
        
        for(String[] termini : freeBBZones){
            ArrayList<Residue> curBFBRes = m.resListFromTermini(termini, flexibleRes);
            
            BBFreeBlock bfb;
            if(voxWidth==0)//automatically select voxel width
                bfb = new BBFreeBlock(curBFBRes,fixAnchors);
            else {
                int numRes = curBFBRes.size();
                int numFreeDOFs = fixAnchors ? 2*numRes-6 : 2*numRes+2;
                double voxLim[] = new double[numFreeDOFs];
                Arrays.fill(voxLim, voxWidth/2);
                bfb = new BBFreeBlock(curBFBRes,fixAnchors,voxLim,tooBigForSeries);
            }
                
            
            ans.add(bfb);
        }
        
        return ans;
    }
}
