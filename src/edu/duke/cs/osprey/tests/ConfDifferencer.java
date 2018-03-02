/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import com.joptimizer.functions.LinearMultivariateRealFunction;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.ematrix.epic.EPICEnergyFunction;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.minimization.SQPMinimizer;
import edu.duke.cs.osprey.plug.Declasher;
import edu.duke.cs.osprey.plug.LPChecks;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import org.apache.commons.math3.optim.linear.LinearConstraint;

/**
 *
 * Assess LUTE feasibility by differencing confs
 * 
 * @author mhall44
 */
public class ConfDifferencer {
    
    static Random rand = new Random();
    
    
    static boolean declashDifferencing = false;
            
    public static void main(String args[]){
        ConfigFileParser cfp = new ConfigFileParser(args);//args 1, 3+ are configuration files
	cfp.loadData();
        
        ConfDifferencer cd = new ConfDifferencer(cfp);
        cd.testDifferencing();
        //DEBUG!!!
        //cd.testDeclashing();
    }
    
    private SearchProblem sp;
    
    public ConfDifferencer(ConfigFileParser cfp){
        sp = prepareSearchProblem(cfp);
    }
    
    public void testDeclashing(){
        int numSamps = 20;
        for(int s=0; s<numSamps; s++){
            DiffSpace ds = sampleDiffSpace(0);//no diff needed, just sample confs
            System.out.println("SAMPLE "+s);
            System.out.println("Main conf: "+ds.mainConf.stringListing());
            
            SQPMinimizer.lowerTol = true;
            sp.luteSettings.useSQP = true;
            sp.luteSettings.minimizeWithPLUG = true;
            MinimizedConf baseMC = new MinimizedConf(ds.baseConf());
            System.out.println("Minimized E (PLUG): "+baseMC.E);
            
            SQPMinimizer.lowerTol = false;
            baseMC = new MinimizedConf(ds.baseConf());
            System.out.println("Minimized E (PLUG, higher tol): "+baseMC.E);
            
            sp.luteSettings.minimizeWithPLUG = false;
            sp.luteSettings.useSQP = false;
            baseMC = new MinimizedConf(ds.baseConf());
            System.out.println("Minimized E (no PLUG): "+baseMC.E);
            
            double Edeclashed = Declasher.declashedE(sp, baseMC.conf);
            System.out.println("Declashed energy: "+Edeclashed);
        }
    }
    
    public void testDifferencing(){
        //ConfETupleExpander te = new ConfETupleExpander(sp);
        //TESampleSet tss = new TESampleSet(te);
        
        int numDiffSamps = 100;
        for(int s=0; s<numDiffSamps; s++){
            //this can be done w/ cplug, sqp or ccd
            DiffSpace ds = sampleDiffSpace(3);
            //ds.mainConf = new RCTuple(new int[]{5,7,7,5,1,7,4});//DEBUG!!!
            System.out.println("SAMPLE "+s);
            System.out.println("Main conf: "+ds.mainConf.stringListing());
            System.out.println("Alt RCs: "+ds.altRCs.stringListing());
            
            for(int bison=0; bison<2; bison++){
                if(bison==1){
                    System.out.println("FLIPPED: ");
                    ds = ds.flipped();
                }
                
                //get the confs in the diffspace
                MinimizedConf baseMC = new MinimizedConf(ds.baseConf());
                //single perturbations
                MinimizedConf mc1 = new MinimizedConf(ds.singlyPerturbedConf(0));
                MinimizedConf mc2 = new MinimizedConf(ds.singlyPerturbedConf(1));
                MinimizedConf mc3 = new MinimizedConf(ds.singlyPerturbedConf(2));
                //pairs
                MinimizedConf mc12 = new MinimizedConf(ds.pairPerturbedConf(0,1));
                MinimizedConf mc23 = new MinimizedConf(ds.pairPerturbedConf(1,2));
                MinimizedConf mc13 = new MinimizedConf(ds.pairPerturbedConf(0,2));
                //all (triple)
                MinimizedConf altMC = new MinimizedConf(ds.altConf());

                System.out.println("Base E: "+baseMC.E);
                System.out.println("Alt E: "+altMC.E);
                
                double singlewiseAltE = mc12.E + mc23.E + mc13.E - 2*baseMC.E;
                System.out.println("1-wise altE: "+singlewiseAltE);

                double pairwiseAltE = baseMC.E + mc12.E + mc23.E + mc13.E - mc1.E - mc2.E - mc3.E;
                System.out.println("Pairwise altE: "+pairwiseAltE);
                
                if(Math.abs(altMC.E-pairwiseAltE)>0.6)
                    System.out.println("Pairwiseness error > thermal");
                else if(Math.abs(altMC.E-pairwiseAltE)>0.3)
                    System.out.println("Pairwiseness error > half of thermal");
                
                
                //ok so there are actually 2 possible kinds of singles differences for min DOFs:
                //direct differencing (requires same number of DOFs),
                //replacement (requires assignment of DOFs to res)
                //thus must combine these in case of mutations with CATS
                //pairwise can come from differencing these
                System.out.println("Base min DOFs: "+baseMC.minDOFs);
                System.out.println("Alt min DOFs: "+altMC.minDOFs);
                
                for(boolean perRCDOFs : new boolean[] {true,false}){//there are multiple ways to do low-order expansions for mindofs, should try them all
                    for(boolean enforceConstr : new boolean[] {true,false}){
                        System.out.println("perRCDOFs: "+perRCDOFs+" enforceConstr: "+enforceConstr);
                
                        DoubleMatrix1D minDOFs1wise = singlewiseMinDOFs(baseMC,altMC,perRCDOFs,enforceConstr,mc1,mc2,mc3);
                        System.out.println("1-wise minDOFs: "+minDOFs1wise);
                        if(noNans(minDOFs1wise) && altMC.mms.getEfunc()!=null){//nans may arise when no constr-allowed mindofs
                            double Emd = altMC.mms.getValue(minDOFs1wise);
                            System.out.println("1-wise minDOFs E: "+Emd);
                            if(Emd<altMC.E-0.1 && enforceConstr)//beat the alt "minimum"
                                System.out.println("Min failure detected");
                        }

                        DoubleMatrix1D minDOFsPairwise = pairwiseMinDOFs(baseMC,altMC,perRCDOFs,enforceConstr,mc1,mc2,mc3,mc12,mc23,mc13);
                        System.out.println("Pairwise minDOFs: "+minDOFsPairwise);
                        if(noNans(minDOFsPairwise) && altMC.mms.getEfunc()!=null){
                            double Emd = altMC.mms.getValue(minDOFsPairwise);
                            System.out.println("Pairwise minDOFs E: "+Emd);
                            if(Emd<altMC.E-0.1 && enforceConstr)//beat the alt "minimum"
                                System.out.println("Min failure detected");
                        }
                        //if pairwiseness is failing, is pairwise/1wise better than real min in at least one direction?
                    }
                }
            }
            
            System.out.println();
            System.out.println();
        }
        
        //we find out if there may be good energies at pairwise-predictable min dofs that
        //just aren't found by minimization
        //we also find out if pairwise E is good for a lot of confs
        //and also if good 1-wise min decomp may explain good pairwise E
        /*sp.luteSettings.minimizeWithPLUG = false;
        sp.luteSettings.useSQP = true;
        System.out.println("SQP Energy: "+sp.EPICMinimizedEnergy(conf1));
        sp.luteSettings.useSQP = false;
        System.out.println("CCD Energy: "+sp.EPICMinimizedEnergy(conf1));*/
    }
    
    
    private static boolean noNans(DoubleMatrix1D x){
        for(double d : x.toArray()){
            if(Double.isNaN(d))
                return false;
        }
        return true;
    }
    
    //approximation of alt min DOFs
    //degrees of freedom in alt conf should either be shared between all confs, or appear in one of the changed rots
    //the approx is handled by pure differencing or by taking the new DOF value respectively
    private DoubleMatrix1D singlewiseMinDOFs(MinimizedConf baseMC, MinimizedConf altMC, boolean perRCDOFs, 
            boolean enforceConstr, MinimizedConf... singlePertConfs){
        //singlewise approximation of min DOFs of alt conf
        ArrayList<DegreeOfFreedom> altDOFs = altMC.mms.getDOFs();
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(altDOFs.size());
        for(int dof=0; dof<altDOFs.size(); dof++){
            DegreeOfFreedom curDOF = altDOFs.get(dof);
            double dofVal;
            if(mcHasDOF(baseMC, curDOF, altMC, perRCDOFs)){//shared DOF.  
                dofVal = -(singlePertConfs.length-1)*baseMC.dofVal(curDOF);
                for(MinimizedConf spc : singlePertConfs)
                    dofVal += spc.dofVal(curDOF);
            }
            else {
                int pertNum = whichPertHasDOF(curDOF,singlePertConfs);
                dofVal = singlePertConfs[pertNum].dofVal(curDOF);
            }
            
            if(enforceConstr){
                DoubleMatrix1D[] constr = altMC.mms.getConstraints();
                if(dofVal < constr[0].get(dof))
                    dofVal = constr[0].get(dof);
                else if(dofVal > constr[1].get(dof))
                    dofVal = constr[1].get(dof);
            }
            
            ans.set(dof, dofVal);
        }
        return ans;
    }
    
    
    private DoubleMatrix1D pairwiseMinDOFs( MinimizedConf baseMC, MinimizedConf altMC, boolean perRCDOFs, 
            boolean enforceConstr,
            MinimizedConf mc1, MinimizedConf mc2, MinimizedConf mc3, 
            MinimizedConf mc12, MinimizedConf mc23, MinimizedConf mc13 ){
        ArrayList<DegreeOfFreedom> altDOFs = altMC.mms.getDOFs();
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(altDOFs.size());
        for(int dof=0; dof<altDOFs.size(); dof++){
            DegreeOfFreedom curDOF = altDOFs.get(dof);
            double dofVal;
            if(mcHasDOF(baseMC, curDOF, altMC, perRCDOFs)){//shared DOF.  Pure differencing
                dofVal = baseMC.dofVal(curDOF)
                        + mc12.dofVal(curDOF) + mc23.dofVal(curDOF) + mc13.dofVal(curDOF)
                        - mc1.dofVal(curDOF) - mc2.dofVal(curDOF) - mc3.dofVal(curDOF);
            }
            else {
                //only one single pert should have dof
                //we can only do differencing among the three confs that have the dof
                int pertNum = whichPertHasDOF(curDOF,mc1,mc2,mc3);
                if(pertNum==0)
                    dofVal = mc12.dofVal(curDOF) + mc13.dofVal(curDOF) - mc1.dofVal(curDOF);
                else if(pertNum==1)
                    dofVal = mc12.dofVal(curDOF) + mc23.dofVal(curDOF) - mc2.dofVal(curDOF);
                else if(pertNum==2)
                    dofVal = mc23.dofVal(curDOF) + mc13.dofVal(curDOF) - mc3.dofVal(curDOF);
                else
                    throw new RuntimeException("ERROR: Some pert should have the DOF");
            }
            
            if(enforceConstr){
                DoubleMatrix1D[] constr = altMC.mms.getConstraints();
                if(dofVal < constr[0].get(dof))
                    dofVal = constr[0].get(dof);
                else if(dofVal > constr[1].get(dof))
                    dofVal = constr[1].get(dof);
            }
            
            ans.set(dof, dofVal);
        }
        return ans;
    }
    
    
    private boolean mcHasDOF(MinimizedConf baseMC, DegreeOfFreedom curDOF, MinimizedConf altMC, boolean perRCDOFs){
        //does baseMC have curDOF (one of the DOFs of altMC)?
        //This can be interpreted strictly (perRCDOFs: if curDOF is a dihedral must be same RC)
        //or loosely (match shared dihedrals even if RC differs)
        if(perRCDOFs && curDOF instanceof FreeDihedral){
            FreeDihedral dih = (FreeDihedral)curDOF;
            String PDBResNum = dih.getResidue().getPDBResNumber();
            //watch out in case molecs in mms and sp not the same...
            int pos = sp.confSpace.getDesignIndex(sp.confSpace.m.getResByPDBResNumber(PDBResNum));
            //ok now see if base and alt have the same RCs at pos
            return baseMC.conf.RCAtPos(pos) == altMC.conf.RCAtPos(pos);
        }
        
        return baseMC.hasDOF(curDOF);
    }
    
    int whichPertHasDOF(DegreeOfFreedom dof, MinimizedConf... MCs){
        //Which of the (presumably perturbed) minimized confs has the DOF?
        for(int pertNum=0; pertNum<MCs.length; pertNum++){
            if(MCs[pertNum].hasDOF(dof))
                return pertNum;
        }
        throw new RuntimeException("ERROR: didn't find the DOF");
    }
    
    
    private static SearchProblem prepareSearchProblem(ConfigFileParser cfp){
        //akch maybz the best way 2 do this is to just load a tupexp mat
        //and look at what's inf
        //proon that
        
        SearchProblem sp = cfp.getSearchProblem();
        
        //lets us load and use a precomputed EPIC matrix (and also LUTE mtx)
        sp.pruneMat = new PruningMatrix(sp.confSpace,0);
        //sp.pruneMat = new PruningMatrix(sp.confSpace, Double.POSITIVE_INFINITY);//always valid
        
        
        if(sp.luteSettings.minimizeWithPLUG){
            //PolytopeMatrix pmat = new PolytopeMatrix(sp, true);
            //sp.plugMat = pmat;
            sp.plugMat = (PolytopeMatrix)ObjectIO.readObject("pmat_070217.dat", true);
            //DEBUG!!!!
            //ObjectIO.writeObject(pmat, "pmat_070317.dat");
        }
        
        sp.loadEnergyMatrix();
        //searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace, pruningInterval);
        /*PruningControl pruningControl = new PruningControl(sp, 5., false, 100, 3, true, true, false, false, false, 100);
        pruningControl.prune();*/
        
        //sp.pruneMat.setOneBody(2, 12, false);//DEBUG!!!
        //CORRESPONDINGLY ALSO IGNORING NULL TERMS IN FULL STERIC POLYTOPE!!
        
        //FOR EPIC - get right confspace
        sp.loadEPICMatrix();
        sp.loadTupExpEMatrix();//DEBUG!!!  gonna use this as a proxy for pruning
        return sp;
    }
    
    
    private DiffSpace sampleDiffSpace(int numAlt){
        int maxNumTries = 200;
        for(int tryNum=0; tryNum<maxNumTries; tryNum++){//how many times we're willing to start from scratch
            boolean success = true;
            DiffSpace ds = new DiffSpace(numAlt);
            while(ds.hasChoicesLeft()){
                if(!ds.makeNextChoice()){
                    success = false;
                    break;
                }
            }
            if(success)
                return ds;
        }
        throw new RuntimeException("ERROR: Couldn't get a good diff space  waaaaaaaaa");
    }
    
    boolean tupleOK(RCTuple tup){//DEBUG!!!  this is sort of weirdly piggybacking on the lute mtx
        return Double.isFinite(sp.tupExpEMat.getInternalEnergy(tup));
    }
    
    
    static int pickElement(ArrayList<Integer> list){
        return list.get(rand.nextInt(list.size()));
    }
    
    int numRCsForPos(int pos){
        return sp.confSpace.posFlex.get(pos).RCs.size();
    }
    
    private class DiffSpace {
        //a little space of confs where we can do differencing bc all the confs are unpruned
        //may be only partially defined
        RCTuple mainConf = new RCTuple();
        RCTuple altRCs = new RCTuple();//not really a tuple, but a set of alternate RCss
        int numAlt;//how many positions have alternates
        
        ArrayList<Integer> unfilledPos;//positions still to fill
        int unfilledAltPos = -1;//alt position to fill (fill this before next regular position, if not -1)
                
        
        DiffSpace(){
            
        }
        
        DiffSpace(int numAlt){
            this.numAlt = numAlt;
            unfilledPos = new ArrayList<>();
            for(int pos=0; pos<sp.confSpace.numPos; pos++)
                unfilledPos.add(pos);
        }
        
        boolean hasChoicesLeft(){
            //haven't fully specified conf and alt confs, need to makeNextChoice more
            return (unfilledAltPos!=-1 || !unfilledPos.isEmpty());
        }
        
        boolean makeNextChoice(){//add another position (regular or alt) to the space
            //return false if impossible
            ArrayList<RCTuple> curTups = enumerateConfs();
            ArrayList<Integer> newRCOptions = new ArrayList();
            
            if(unfilledAltPos!=-1){
                //pick a new alternate position
                int pickPos = unfilledAltPos;
                int origRC = mainConf.RCAtPos(pickPos);
                for(int pickRC=0; pickRC<numRCsForPos(pickPos); pickRC++){
                    if(pickRC != origRC){//no point in repeating ourselves
                        for(RCTuple tup : curTups){
                            RCTuple augTup = tup.subtractPos(pickPos).addRC(pickPos, pickRC);
                            if(tupleOK(augTup))
                                newRCOptions.add(pickRC);
                        }
                    }
                }
                if(newRCOptions.isEmpty())//no options...
                    return false;
                else {
                    altRCs = altRCs.addRC(pickPos, pickElement(newRCOptions));
                    unfilledAltPos = -1;
                    return true;
                }
            }
            else if ( ! unfilledPos.isEmpty() ){
                //pick a new regular position
                int pickPos = pickElement(unfilledPos);
                for(int pickRC=0; pickRC<numRCsForPos(pickPos); pickRC++){
                    for(RCTuple tup : curTups){
                        RCTuple augTup = tup.addRC(pickPos, pickRC);
                        if(tupleOK(augTup))
                            newRCOptions.add(pickRC);
                    }
                }
                if(newRCOptions.isEmpty())//no options...
                    return false;
                else {
                    mainConf = mainConf.addRC(pickPos, pickElement(newRCOptions));
                    unfilledPos.remove(Integer.valueOf(pickPos));
                    if(altRCs.size()<numAlt){
                        //need another alternate, let's have it at this position
                        unfilledAltPos = pickPos;
                    }
                    return true;
                }
            }
            else
                throw new RuntimeException("Tried to makeNextChoice but free will is an illusion");
        }
        
        ArrayList<RCTuple> enumerateConfs(){//enumerate (possibly partial) confs
            ArrayList<RCTuple> ans = new ArrayList<>(Arrays.asList(mainConf));
            for(int altIndex=0; altIndex<altRCs.size(); altIndex++){
                ArrayList<RCTuple> groovin = new ArrayList<>();
                int curAltPos = altRCs.pos.get(altIndex);
                for(RCTuple t00p : ans){
                    groovin.add(t00p);
                    groovin.add(t00p.subtractPos(curAltPos).addRC(curAltPos, altRCs.RCs.get(altIndex)));
                }
                ans = groovin;
            }
            return ans;
        }
        
        DiffSpace flipped(){
            //flip alternates and main
            if(hasChoicesLeft())
                throw new RuntimeException("ERROR flipping is for full confs");
            DiffSpace flippy = new DiffSpace();
            flippy.unfilledPos = new ArrayList<>();
            flippy.numAlt = numAlt;
            
            flippy.mainConf = mainConf.copy();
            flippy.altRCs = altRCs.copy();
            for(int altIndex=0; altIndex<numAlt; altIndex++){
                int pos = altRCs.pos.get(altIndex);
                int mainIndex = mainConf.pos.indexOf(pos);
                flippy.altRCs.RCs.set(altIndex, mainConf.RCs.get(mainIndex));
                flippy.mainConf.RCs.set(mainIndex, altRCs.RCs.get(altIndex));
            }
            return flippy;
        }
        
        RCTuple baseConf(){
            return mainConf;
        }
        
        RCTuple altConf(){
            RCTuple ans = mainConf.copy();
            for(int index=0; index<altRCs.size(); index++){
                ans.RCs.set(mainConf.pos.indexOf(altRCs.pos.get(index)), altRCs.RCs.get(index));
            }
            return ans;
        }
        
        RCTuple singlyPerturbedConf(int index){
            //perturbed only based on the specified index in altRCs
            int pertPos = altRCs.pos.get(index);
            return mainConf.subtractPos(pertPos).addRC(pertPos, altRCs.RCs.get(index));
        }
        
        RCTuple pairPerturbedConf(int index, int index2){
            int pertPos = altRCs.pos.get(index);
            int pertPos2 = altRCs.pos.get(index2);
            return mainConf.subtractPos(pertPos).subtractPos(pertPos2).
                    addRC(pertPos,altRCs.RCs.get(index)).addRC(pertPos2,altRCs.RCs.get(index2));
        }
    }
    
    
    private class MinimizedConf {
        RCTuple conf;
        MoleculeModifierAndScorer mms;
        double E;
        DoubleMatrix1D minDOFs;
        
        MinimizedConf(RCTuple tup){
            conf = tup;
            EPICEnergyFunction efunc = sp.epicMat.internalEnergyFunction(tup,false);
            
            //efunc.includeMinE = true;//might was well just difference the cont parts...
            mms = new MoleculeModifierAndScorer(efunc,sp.epicMat.getConfSpace(),tup); 
            if(efunc==null){//indicates bad conf
                minDOFs = DoubleFactory1D.dense.make(mms.getNumDOFs(), Double.NaN);
                E = Double.POSITIVE_INFINITY;
                return;
            }
            
            Minimizer.Result minResult;
            if(declashDifferencing){
                minResult = Declasher.declash(sp, conf);
                if(minResult.dofValues == null)
                    minResult = null;
            }
            else {
                Minimizer min;

                if(sp.luteSettings.minimizeWithPLUG){
                    ArrayList<LinearConstraint> constrList = sp.plugMat.getFullStericPolytope(conf, mms.getDOFs());

                    LinearMultivariateRealFunction constr[] = new LinearMultivariateRealFunction[constrList.size()];
                    for(int c=0; c<constrList.size(); c++){
                        constr[c] = LPChecks.toLinearMultivariateRealFunction(constrList.get(c));
                    }

                    min = new SQPMinimizer(mms,constr);
                }
                else if(sp.luteSettings.useSQP)
                    min = new SQPMinimizer(mms,null);
                else
                    min = new CCDMinimizer(mms,false);

                minResult = min.minimize();
            }
            
            if(minResult==null){
                minDOFs = DoubleFactory1D.dense.make(mms.getNumDOFs(), Double.NaN);
                E = Double.POSITIVE_INFINITY;
            }
            else {
                minDOFs = minResult.dofValues;
                E = minResult.energy;
            }
        }
        
        boolean hasDOF(DegreeOfFreedom dof){
            return mms.getDOFs().contains(dof);
        }
        
        double dofVal(DegreeOfFreedom dof){
            int index = mms.getDOFs().indexOf(dof);
            if(index==-1)
                throw new RuntimeException("ERROR: Expected mms to have dof but it doesn't");
            return minDOFs.get(index);
        }
    }
    
    //Questions: How often are there significant (>~thermal) deviations from pairwiseness?
    //A: Maybe 10% for 1CC8 w/ PLUG; seems to be a big contributor to resid
    //When they happen, are they because the minimum implied by the pairwise model
    //(or rather for concreteness, the 1- or 2-wise min pert models) is really worse than full-min,
    //or because they just happen to minimize somewhere else?
    // A: Often really worse.  min pert pairwise model in general not as good as pairwise E. 
    //Sample 57: 10x thermal error, and 4-5x thermal min error where min pert pairwise E < alt E with constr on
    //But there are min failures even with PLUG.  
    //(See ConfDifferencer_plugmin_070317.txt)
    //Can 1- or 2-wise min pert be improved by not counting DOFs as shared if they're chi's from different RCs?
    //How about by handling vox maxing out better?  Can this improvement be ported to E?
    //These don't seem consistently better, sometimes I guess
    //Can PLUG help by providing better init vals for min?  Max buffer?
    //Maybe... (see Declasher)
    
    //So analysis based on pairwise (and especially 1-wise) min pert is not expected to improve PLUG.
    //The energy expansion is what it is, if resid is good enough then still can use it.  
    //Also global minimization errors are relatively common even with PLUG.  
    //How about de-clash???
}
