/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug.triage;

import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author mhall44
 */
public class RCClassTuple {
    
    ArrayList<RCClass> classes;
    
    public RCClassTuple(RCClass... cl){
        classes = new ArrayList(Arrays.asList(cl));
    }
    
    RCClass get(int i){
        return classes.get(i);
    }
    
    
    int size(){
        return classes.size();
    }
}
