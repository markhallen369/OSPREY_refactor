/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.TermECalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.plug.RCPairVDWChecker;
import edu.duke.cs.osprey.plug.RCTuplePolytope;

/**
 *
 * Fit just one EPIC term
 * This is to improve EPIC by dealing with problem terms
 * 
 * @author mhall44
 */
public class OneEPICFit {
    
    
    
    public static void main(String[] args){
        //Specify desired system and term here
        //args = new String[] {"-c","KStar.cfg","findGMEC","System.cfg","DEE.cfg"};
        int pos1=3;
        int rc1=61;
        int pos2=12;
        int rc2=76;
        //ROT:3 1 0 12 12 10
        
        RCTuple RCs = new RCTuple(1,7,3,5);//(1,0,0,7);//(pos1,rc1,pos2,rc2);
        
        
        ConfigFileParser cfp = new ConfigFileParser(args);
        cfp.loadData();
        SearchProblem sp = cfp.getSearchProblem();
        
        //Let's see what PLUG can do here
        RCPairVDWChecker rpvc = new RCPairVDWChecker(sp.confSpace, RCs, sp.shellResidues);
        RCTuplePolytope tope = rpvc.buildPolytope();
        PolytopeMatrix plugMat = new PolytopeMatrix(sp.confSpace);//build matrix just for one entry lol
        plugMat.setTupleValue(RCs, tope);
                
        TermECalculator tec = new TermECalculator(sp.confSpace, sp.shellResidues, 
                true, false, null, new EPICSettings(), false, plugMat, 0);//pos1, pos2);
            
        tec.calcTupleEnergy(RCs);//this might crash on trying to store the fit,
        //but that's OK, just want to see fitting output
    }
    
}
