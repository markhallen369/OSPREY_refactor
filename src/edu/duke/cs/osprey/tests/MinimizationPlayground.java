/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import com.joptimizer.functions.LinearMultivariateRealFunction;
import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.bbfree.CoeffGetter;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.ematrix.epic.EPICEnergyFunction;
import edu.duke.cs.osprey.ematrix.epic.EPoly;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ResPairEnergy;
import edu.duke.cs.osprey.energy.forcefield.SingleResEnergy;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.minimization.SQPMinimizer;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import edu.duke.cs.osprey.tools.ObjectIO;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * Playing with minimization
 * 
 * @author mhall44
 */
public class MinimizationPlayground {
    
    
    private static class foo implements Serializable {
        int a;
        double b;
        transient LinearMultivariateRealFunction lmrf;
        public foo(){}
        public foo(int a, double b, double[] q, double r){
            this.a = a;
            this.b = b;
            lmrf = new LinearMultivariateRealFunction(q,r);
        }
        private void readObject(java.io.ObjectInputStream stream)
         throws IOException, ClassNotFoundException {
         stream.defaultReadObject();
         double q[] = (double[])stream.readObject();
         double r = stream.readDouble();
         lmrf = new LinearMultivariateRealFunction(q,r);
        }
        private void writeObject(java.io.ObjectOutputStream stream)
         throws IOException {
         stream.defaultWriteObject();
         stream.writeObject(lmrf.getQ().toArray());
         stream.writeDouble(lmrf.getR());
        }
    }
        
    public static void main(String args[]){
        
        /*foo f = new foo(5,3.2,new double[]{3,34,5},-1.7);
        ObjectIO.writeObject(f, "f.dat");
        foo f2 = (foo)ObjectIO.readObject("f.dat",false);
        System.exit(0);*/
        
        /*SQPMinimizer minim = (SQPMinimizer)ObjectIO.readObject("SQP_MINIMIZER_STEP90.dat", false);
        minim.minimizeFromCurPoint();
        System.exit(0);*/
        
        ConfigFileParser cfp = new ConfigFileParser(args);//args 1, 3+ are configuration files
	cfp.loadData();
        
        //listing clashes, pmat precomputed
        //DEBUG!!
        /*PolytopeMatrix pmat = (PolytopeMatrix)ObjectIO.readObject("pmat.dat", false);
        int[] conf = new int[] {5,7,12,5,0,7,4};
        RCTuple confTuple = new RCTuple(conf);
        System.out.println("Checking polytope violations for CCD opt vals: ");
        //DoubleMatrix1D optDOFVals = DoubleFactory1D.dense.make(new double[] {-57.913409, 166.684413, -0.065158, -0.229097, 0.234918, -0.202681, 0.149283, -0.069924, -0.296158, -0.455166, -183.766828, 51, -58.078234, -64.192948, -79, -56.226172, 87.764445, 9, -67.947732, -69.06033, -31, -69.087609, 168.660219});
        DoubleMatrix1D optDOFVals = DoubleFactory1D.dense.make(new double[] {-60.81941, 169.162243, -185.934169, 51, -56, -74, -79, -56.692953, 86.192216, 9, -63.004258, -63.574936, -31, -73.119956, 170.102228});
        MoleculeModifierAndScorer mms = new MoleculeModifierAndScorer(null,pmat.cSpace,confTuple);
        pmat.listClashes(conf, optDOFVals, mms.getDOFs());
        mms.setDOFs(optDOFVals);
        PDBFileWriter.writePDBFile(pmat.cSpace.m, "OPT_DOF_VALS.pdb");
        System.out.println("Checking polytope violations for SQP opt vals: ");
        //DoubleMatrix1D sqpDOFVals = DoubleFactory1D.dense.make(new double[] {-58.132046, 166.952221, -0.056139, -0.27844, 0.334519, -0.247856, -0.016941, -0.150404, -0.306353, -0.455166, -181.534193, 51, -57.582896, -64.845457, -79, -56, 87.563931, 9, -68.561964, -74, -49, -68.638399, 168.256258});
        DoubleMatrix1D sqpDOFVals = DoubleFactory1D.dense.make(new double[] {-60.81944, 169.162236, -185.927224, 51, -56, -74, -79, -56.695492, 86.192418, 9, -63.156913, -63.406383, -31, -73.11997, 170.102227});
        pmat.listClashes(conf, sqpDOFVals, mms.getDOFs());
        mms.setDOFs(sqpDOFVals);
        PDBFileWriter.writePDBFile(pmat.cSpace.m, "SQP_DOF_VALS.pdb");
        System.out.println("Done checking polytope violations");
        System.exit(0);*/
        
        
        //computing pmat
        SearchProblem sp = cfp.getSearchProblem();
        
        
        //lets us load and use a precomputed EPIC matrix
        sp.pruneMat = new PruningMatrix(sp.confSpace,0);
        
        /*
        sp.pruneMat = new PruningMatrix(sp.confSpace, Double.POSITIVE_INFINITY);//always valid
        PolytopeMatrix pmat = new PolytopeMatrix(sp, true);
        sp.plugMat = pmat;*/
        
        /*ObjectIO.writeObject(pmat, "pmat.dat");
        ObjectIO.writeObject(sp, "SEARCHPROB.dat");*/
        
        //loading pmat
        //SearchProblem sp = (SearchProblem)ObjectIO.readObject("SEARCHPROB.dat", false);
        //PolytopeMatrix pmat = sp.plugMat;
        
        sp.loadEnergyMatrix();
        //searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace, pruningInterval);
        /*PruningControl pruningControl = new PruningControl(sp, 5., false, 100, 3, true, true, false, false, false, 100);
        pruningControl.prune();*/
        
        //sp.pruneMat.setOneBody(2, 12, false);//DEBUG!!!
        //CORRESPONDINGLY ALSO IGNORING NULL TERMS IN FULL STERIC POLYTOPE!!
        
        //FOR EPIC - get right confspace
        sp.loadEPICMatrix();
        //sp.confSpace = sp.epicMat.getConfSpace();//swap in right confspace.  DEBUG!!
        //sp.plugMat = new PolytopeMatrix(sp, false);
        
        //minimization
        
        //int conf1[] = new int[15];//start with all wt rots
        //conf1[12] = 5;//TRP
        int[] conf1 = new int[] {5,7,12,5,0,7,4};
        //int[] conf1 = new int[] {5,0,10,7,0,7,6};
        /*DoubleMatrix1D xEPIC = DoubleFactory1D.dense.make(new double[]
            {-57.969192, 166.735048, -0.066281, -0.229844, 0.245356, -0.203498, 0.140798, -0.071649, -0.304984, -0.455166, -183.708158, 51, -58.074254, -64.066479, -79, -56.170042, 87.783936, 9, -68.044382, -69.201669, -31, -69.036891, 168.647869}
            );

        checkEPICParts(sp, conf1, xEPIC);
        
        
        //System.out.println("regular Energy: "+sp.minimizedEnergy(conf1));
        //System.out.println("EPIC Energy: "+sp.EPICMinimizedEnergy(conf1));
        System.exit(0);*/
        
        sp.luteSettings.minimizeWithPLUG = false;
        sp.luteSettings.useSQP = true;
        System.out.println("SQP Energy: "+sp.EPICMinimizedEnergy(conf1));
        sp.luteSettings.useSQP = false;
        System.out.println("CCD Energy: "+sp.EPICMinimizedEnergy(conf1));
        
        //System.out.println("Energy: "+sp.minimizedEnergy(conf1));
        //System.out.println("Energy: "+sp.EPICMinimizedEnergy(conf1));
        //sp.outputMinimizedStruct(conf1, "VOXMIN.pdb");
    }
    
    
    private static void checkEPICParts(SearchProblem sp, int[] conf1, DoubleMatrix1D xVals){

       /* System.out.println("SP COEFFS: ");
        for(DegreeOfFreedom dof : sp.confSpace.confDOFs){
            if(dof instanceof BBFreeDOF){
                System.out.println( CoeffGetter.getCoeffs((BBFreeDOF)dof) );
            }
        }
        System.out.println("EPIC COEFFS: ");
        for(DegreeOfFreedom dof : sp.epicMat.getConfSpace().confDOFs){
            if(dof instanceof BBFreeDOF){
                System.out.println( CoeffGetter.getCoeffs((BBFreeDOF)dof) );
            }
        }
        System.exit(0);*/
        
        
        //EnergyFunction ef, ConfSpace cSpace, RCTuple RCTup
        EPICEnergyFunction epicEFunc = sp.epicMat.internalEnergyFunction(new RCTuple(conf1),true);
        
        MoleculeModifierAndScorer mofReg = new MoleculeModifierAndScorer(sp.fullConfE, sp.confSpace, new RCTuple(conf1));
        MoleculeModifierAndScorer mofEPIC = new MoleculeModifierAndScorer(epicEFunc, sp.epicMat.getConfSpace(), new RCTuple(conf1));
        mofReg.setDOFs(xVals);
        mofEPIC.setDOFs(xVals);
        
        ArrayList<Double> epicTerms = epicEFunc.allTermValues();
        double Ereg = 0;
        double Eepic = 0;
        int termCount=0;

        ArrayList<String> flexRes = sp.confSpace.flexibleRes;
        
        for(int pos=0; pos<conf1.length; pos++){
            System.out.println("One-body for "+pos);
            double epicTerm = epicTerms.get(termCount);// + sp.emat.getOneBody(pos, conf1[pos]);
            String res1Name = flexRes.get(pos);
            double regTerm = sumEFuncTerms((MultiTermEnergyFunction)sp.fullConfE, res1Name, null, flexRes);
            System.out.println("EPIC term: "+epicTerm+" regular term: "+regTerm);
            Ereg += regTerm;
            Eepic += epicTerm;
            termCount++;
            for(int pos2=0; pos2<pos; pos2++){
                System.out.println("Pairwise for "+pos+" and "+pos2);
                epicTerm = epicTerms.get(termCount);//MINE IS NEEDED IN CASE MINIMA DIFFER//sp.emat.getPairwise(pos, conf1[pos], pos2, conf1[pos2]);
                String res2Name = flexRes.get(pos2);
                regTerm = sumEFuncTerms((MultiTermEnergyFunction)sp.fullConfE, res1Name, res2Name, flexRes);
                System.out.println("EPIC term: "+epicTerm+" regular term: "+regTerm);
                Ereg += regTerm;
                Eepic += epicTerm;
                termCount++;
            }
        }
        
        System.out.println("Total Ereg: "+Ereg);
        System.out.println("Total Eepic: "+Eepic);
        
        epicEFunc.includeMinE = false;//calculation as before
        
        System.out.println("regular Energy check: "+mofReg.getValue(xVals));
        System.out.println( "EPIC Energy check: " + (sp.emat.confE(conf1)+mofEPIC.getValue(xVals)) );
                
        
        System.out.println("CHECKING TERM 5,4");
        EPoly ep = sp.epicMat.getPairwise(5, 7, 4, 0);
        double pairBound = sp.emat.getPairwise(5, 7, 4, 0);

        int termDOFs[] = new int[] {18,19,20,2,3,4,5,6,7,8,9};
        mofEPIC.setDOFs(xVals);
        DoubleMatrix1D xTerm = DoubleFactory1D.dense.make(termDOFs.length);
        for(int d=0; d<termDOFs.length; d++)
            xTerm.set(d, xVals.get(termDOFs[d]));
        System.out.println("Base value: "+ep.evaluate(xTerm, false, true));
        mofReg.setDOFs(xVals);
        System.out.println("Base reg check: "+(sumEFuncTerms((MultiTermEnergyFunction)sp.fullConfE, flexRes.get(4), flexRes.get(5), flexRes)-pairBound));
        
        double otherDOFVals[][] = new double[][] {
            new double[] {-56, -56, -49, 0.455166, -0.107276, 0.455166, -0.292837, 0.455166, 0.089244, -0.455166, -0.003245},
            new double[] {-10.63327459165211, -9.167270316149569, 6.985844252411099, -0.24142498452620087, 0.26037350862564523, -0.05627662436052938, 0.3039165249188948, -0.5774368876815295, 0.24833056393202224, 0.3684391962210501, 0.15802742635836925},
            new double[] {-1.0732686002828693, -6.757939953877447, 4.2944963512546686, -0.7350561782243951, -0.10331189306957919, -0.791498025715407, 0.5053677398040649, -0.08497922244146539, -0.4535507546358377, 0.11221793796208426, -0.0782301225058552},
            new double[] {-12.30756222915431, -7.313838707007184, 10.519662140205526, -0.1811920284655939, 0.1922155969096836, -0.7833398167005862, 0.3349855444461576, -0.8138050517454138, 0.15060573172230207, 0.34339243937603237, 0.3805897897824038}
        };
        
        for(double otherVals[] : otherDOFVals){
            DoubleMatrix1D z = xVals.copy();
            for(int d=0; d<termDOFs.length; d++)
                z.set(termDOFs[d], otherVals[d]);
            mofEPIC.setDOFs(z);
            DoubleMatrix1D ov = DoubleFactory1D.dense.make(otherVals);
            System.out.println("EPIC term value is "+ep.evaluate(ov, false, true)+" at "+ov);
            mofReg.setDOFs(z);
            System.out.println("Reg term value: "+(sumEFuncTerms((MultiTermEnergyFunction)sp.fullConfE, flexRes.get(4), flexRes.get(5), flexRes)-pairBound));
        }
    }
    
    
    private static double sumEFuncTerms(MultiTermEnergyFunction mtef, String res1Name, String res2Name, ArrayList<String> allFlexResNames){
        //sum up terms in mtef for the desired 1-body or pairwise term
        double E = 0;
        for(int t=0; t<mtef.getTerms().size(); t++){
            EnergyFunction term = mtef.getTerms().get(t);
            if(term instanceof SingleResEnergy){
                if(res2Name==null){//only can contribute to a 1-body term
                    if( ((SingleResEnergy)term).getRes().getPDBResNumber().equalsIgnoreCase(res1Name) )
                        E += mtef.getCoeffs().get(t) * term.getEnergy();
                }
            }
            else if(term instanceof ResPairEnergy){
                boolean termApplies = false;
                String termRes1Name = ((ResPairEnergy)term).getRes1().getPDBResNumber();
                String termRes2Name = ((ResPairEnergy)term).getRes2().getPDBResNumber();
                if(res2Name==null){//want res1Name - shell interactions
                    if(res1Name.equalsIgnoreCase(termRes1Name))
                        termApplies = ! allFlexResNames.contains(termRes2Name);//ok don't worry about case...
                    else if(res1Name.equalsIgnoreCase(termRes2Name))
                        termApplies = ! allFlexResNames.contains(termRes1Name);
                }
                else {
                    if(termRes1Name.equalsIgnoreCase(res1Name))
                        termApplies = termRes2Name.equalsIgnoreCase(res2Name);
                    else if (termRes1Name.equalsIgnoreCase(res2Name))
                        termApplies = termRes2Name.equalsIgnoreCase(res1Name);
                }
                if(termApplies)
                    E += mtef.getCoeffs().get(t) * term.getEnergy();
            }
            else
                throw new RuntimeException("ERROR: Unsupport efunc type");
        }
        return E;
    }
    
}
