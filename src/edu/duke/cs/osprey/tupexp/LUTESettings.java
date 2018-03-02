/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tupexp;

import edu.duke.cs.osprey.control.ParamSet;
import java.io.Serializable;

/**
 *
 * Settings for LUTE
 * 
 * @author mhall44
 */
@SuppressWarnings("serial")
public class LUTESettings implements Serializable {
    
    boolean useLUTE=false;
    public double goalResid=0.01;
    boolean useRelWt=false;
    boolean useThreshWt=false;
    
    public boolean useSQP = false;
    public boolean minimizeWithPLUG = false;
    
    public boolean usePLUGEnhancedMinimizer = false;
    public int numConsistent = 1;//how many consistent minima the PLUG-enhanced minimizer should aim for
    //(1 means just minimize once like usual)
    
    public LUTESettings(){
        //by default, no LUTE
        useLUTE = false;
    }
    
    public LUTESettings(ParamSet params){
        //initialize from input parameter set
        useLUTE = params.getBool("USETUPEXP");
        goalResid = params.getDouble("LUTEGOALRESID");
        useRelWt = params.getBool("LUTERELWT");
        useRelWt = params.getBool("LUTETHRESHWT");
        
        //SQP has more of a tendency to find local min far from the initial point
        //and minimizing only w/i the PLUG voxel gives it a weird shape that sometimes hampers LUTE
        //so turning these off by default
        minimizeWithPLUG = params.getBool("LUTEMINIMIZEWITHPLUG",false);
        useSQP = params.getBool("LUTEUSESQP",false);
        
        //PLUG-enhanced minimizer will stick to our empirically decent CCD
        //but use PLUG to ensure we go downhill from a declashed initial pt
        usePLUGEnhancedMinimizer = params.getBool("LUTEPLUGENHANCEMIN");
        numConsistent = params.getInt("LUTENUMCONSISTENT",1);
    }
}
