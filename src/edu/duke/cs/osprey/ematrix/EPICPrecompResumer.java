/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import java.io.File;
import java.io.Serializable;

/**
 *
 * Saves progress in EPIC precomputation
 * 
 * @author mhall44
 */
public class EPICPrecompResumer implements Serializable {
    
    EPICMatrix partialMtx;
    int curRes = 0;
    String fileName;
    long lastSaveTime = 0;
    final static int saveInterval = 1800000;//ms
    
    
    public EPICPrecompResumer(EPICMatrix partialMtx, String fileName){
        this.partialMtx = partialMtx;
        this.fileName = fileName;
    }
    
    void saveIfTime(){
        //save to disk if saveInterval has elapsed since last save
        if(System.currentTimeMillis() > lastSaveTime + saveInterval){
            ObjectIO.writeObject(this, fileName);
            lastSaveTime = System.currentTimeMillis();
        }
    }
    
    void deleteFile(){
        if(lastSaveTime>0)//we have saved to a file, which is no longer needed and could confuse future runs
            new File(fileName).delete();
    }
}
