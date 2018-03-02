/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.xminfit;

import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.TupleEnumerator;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.plug.RCPairVDWChecker;
import edu.duke.cs.osprey.plug.RCTuplePolytope;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tupexp.BasicPruningTupleExpander;
import edu.duke.cs.osprey.tupexp.ConfETupleExpander;
import edu.duke.cs.osprey.tupexp.TESampleSet;
import edu.duke.cs.osprey.tupexp.TupExpChooser;
import edu.duke.cs.osprey.tupexp.TupleExpander;
import java.util.ArrayList;
import java.util.Random;
import org.apache.commons.math3.optim.linear.LinearConstraint;

/**
 *
 * This class is for playing with minimization of CATS DOFs not separately for each voxel,
 * but as a tuple expansion that is the best fit for a bunch of voxels
 * 
 * Going to try a pairwise expansion wrt all pairs of residues that come into contact
 * of which at least one has CATS flex -- these should be the most important terms
 * 
 * Need to compare to full minimization and check gradients at "min"
 * 
 * 
 * @author mhall44
 */
public class XMinFitPlayground {
    
    SearchProblem searchSpace;
    
    XMFMatrix XMinMatrix;
    
    static boolean oneConf = true;//testing minimizer on single conf
    int singleConf[] = new int[] {5,7,7,5,0,7,4};
    
        //use args like for Main
    public static void main(String args[]){
        ConfigFileParser cfp = new ConfigFileParser(args);
        cfp.loadData();
        XMinFitPlayground xp = new XMinFitPlayground(cfp.getSearchProblem());
        xp.testXMinFit();
    }
    
    
    public XMinFitPlayground(SearchProblem searchSpace){
        this.searchSpace = searchSpace;
        //initXMinMatrix will initialize XMinMatrix
    }
    
    
    
    
    public void testXMinFit(){
        //given a cfp-generated search problem, test LUTE matrix generation
        //with various parameters
        //only plug pruning because the problem cases are high-ival and maybe even non-pairwise
        searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace, Double.POSITIVE_INFINITY);
        
        if(oneConf){
            pruneAllButOneConf();
        }
        else {
            searchSpace.plugMat = new PolytopeMatrix(searchSpace,true);
            ObjectIO.writeObject(searchSpace.plugMat, searchSpace.name+".PLUG.dat");
            //searchSpace.plugMat = (PolytopeMatrix) ObjectIO.readObject(searchSpace.name+".PLUG.dat",false);//DEBUG!!
        }
        
        searchSpace.loadEnergyMatrix();//The tuple enumerator wants this
        searchSpace.loadEPICMatrix();
        
        if(oneConf)
            initXMinMatrixOneConf();
        else
            initXMinMatrix();//pick out tuples to use
        
        doXMinFit();//do the fit minimization
        testXMinMatrix();
    }
    
    
    
    public void initXMinMatrix(){
        //there will be a contribution for every contacting RC pair,
        //and how about for all singles too
        XMinMatrix = new XMFMatrix(searchSpace.confSpace);
        
        /*for(int pos=0; pos<XMinMatrix.getNumPos(); pos++){
            for(int rc=0; rc<XMinMatrix.getNumConfAtPos(pos); rc++){
                //singles
                XMinMatrix.initOneBody(pos,rc);
                for(int pos2=0; pos2<pos; pos2++){
                    for(int rc2=0; rc2<XMinMatrix.getNumConfAtPos(pos2); rc2++){
                        RCPairVDWChecker checker = new RCPairVDWChecker(searchSpace.confSpace,new RCTuple(pos,rc,pos2,rc2),searchSpace.shellResidues);
                        ArrayList<LinearConstraint> polytopeConstr = checker.calcPolytopeConstr(new ArrayList<>());//list of potential clashes
                        if(polytopeConstr==null)//impossible tuple...mark as null here too
                            XMinMatrix.setPairwise(pos, rc, pos2, rc2, null);
                        else if(!polytopeConstr.isEmpty()){//DEBUG!! could use a less hacky way to check for contacts
                            int numVoxConstr = 2*polytopeConstr.get(0).getCoefficients().getDimension();
                            if(polytopeConstr.size()<numVoxConstr)
                                throw new RuntimeException("ERROR: Expected more voxel constraints");
                            else if(polytopeConstr.size()>numVoxConstr)//there are contacts.  Let's use this pair
                                XMinMatrix.initPairwise(pos, rc, pos2, rc2);
                        }
                    }
                }
            }
        }*/
        //let's reuse the constraints from PLUG mtx to determine contacts
        for(int pos=0; pos<XMinMatrix.getNumPos(); pos++){
            for(int rc=0; rc<XMinMatrix.getNumConfAtPos(pos); rc++){
                //singles
                XMinMatrix.initOneBody(pos,rc);
                for(int pos2=0; pos2<pos; pos2++){
                    for(int rc2=0; rc2<XMinMatrix.getNumConfAtPos(pos2); rc2++){
                        RCTuplePolytope rtp = searchSpace.plugMat.getPairwise(pos,rc,pos2,rc2);//lists potential clashes
                        if(rtp==null)//impossible tuple...mark as null here too
                            XMinMatrix.setPairwise(pos, rc, pos2, rc2, null);
                        else if(!rtp.getConstr().isEmpty()){//DEBUG!! could use a less hacky way to check for contacts
                            ArrayList<LinearConstraint> polytopeConstr = rtp.getConstr();
                            int numVoxConstr = 2*polytopeConstr.get(0).getCoefficients().getDimension();
                            if(polytopeConstr.size()<numVoxConstr)
                                throw new RuntimeException("ERROR: Expected more voxel constraints");
                            else if(polytopeConstr.size()>numVoxConstr)//there are contacts.  Let's use this pair
                                XMinMatrix.initPairwise(pos, rc, pos2, rc2);
                        }
                    }
                }
            }
        }
        
    }
    
    public void doXMinFit(){
        ArrayList<SampleXMin> samples = drawSamples();
        /*GradDescent*/BatchCCDXMinFitter.doXMinFit(samples,XMinMatrix);//DEBUG!!!
    }
    
    
    
    public void testXMinMatrix(){
        //let's start by looking at average (signed and unsigned) differences in energy
        //between XMin and regular
        //also differences in x?
        ArrayList<SampleXMin> samples = drawSamples();
        double signedEDiffSum = 0;
        double unsignedEDiffSum = 0;
        int numXMinFitHigher = 0;
        int numXMinFitLower = 0;
        double distSum = 0;
        for(SampleXMin samp : samples){
            double[] dEandDist = samp.compareXMinToMatrix(XMinMatrix);
            double dE = dEandDist[0];
            signedEDiffSum += dE;
            unsignedEDiffSum += Math.abs(dE);
            if(dE>0)
                numXMinFitLower++;
            else
                numXMinFitHigher++;
            distSum += dEandDist[1];
        }
        
        System.out.println( "XMinMatrix test.  "
                + "\n Average signed energy difference (fullmin-matrix): "+(signedEDiffSum/samples.size())
                + "\n Average unsigned energy difference (fullmin-matrix): "+(unsignedEDiffSum/samples.size())
                + "\n Number of vox where xminfit is higher energy: "+numXMinFitHigher
                + "\n Number of vox where xminfit is lower energy: "+numXMinFitLower
                + "\n Average distance between full & matrix minima: "+(distSum/samples.size()) );
        
        testLUTEFits();
    }
    
    
    public void testLUTEFits(){
        //see how well we can do a LUTE fit
        //compare regular LUTE with XMinMatrix
        TupleExpander baseExpander = new ConfETupleExpander(searchSpace,null);
        TupleExpander xExpander = makeXTupleExpander();
        
        TupleEnumerator tupEnum = new TupleEnumerator(searchSpace.pruneMat,searchSpace.emat,
                searchSpace.plugMat,searchSpace.confSpace.numPos);
        TupExpChooser baseChooser = new TupExpChooser(baseExpander, tupEnum);//make a chooser to choose what tuples will be in the expansion
        TupExpChooser xChooser = new TupExpChooser(xExpander, tupEnum);
        
        System.out.println("BASE PAIRWISE EXPANSION: ");
        baseChooser.calcPairwiseExpansion();
        System.out.println();
        System.out.println("X PAIRWISE EXPANSION: ");
        xChooser.calcPairwiseExpansion();
        System.out.println();
        
        System.out.println("BASE 2-PARTNER EXPANSION: ");
        baseChooser.calcExpansionResTriples(2);
        System.out.println();
        System.out.println("X 2-PARTNER EXPANSION: ");
        xChooser.calcExpansionResTriples(2);
        System.out.println();
        
        System.out.println("BASE 5-PARTNER EXPANSION: ");
        baseChooser.calcExpansionResTriples(5);
        System.out.println();
        System.out.println("X 5-PARTNER EXPANSION: ");
        xChooser.calcExpansionResTriples(5);
        System.out.println();
    }
    
    
    private TupleExpander makeXTupleExpander(){
        //A tuple expander where energies are for the xminfit confs
        return new ConfETupleExpander(searchSpace,null){
            @Override
            public double scoreAssignmentList(int[] assignmentList) {
                SampleXMin samp = new SampleXMin(searchSpace,assignmentList);
                return samp.curEnergy(XMinMatrix);
            }
        };
    }
    
    
    private ArrayList<SampleXMin> drawSamples(){
        //draw samples covering the tuples currently populated in the XMinMatrix
        //use tss: Make a te, add tuples in XMinMatrix directly to tuples, then
        //make tss from that
        TupleExpander te = new BasicPruningTupleExpander(searchSpace.pruneMat,searchSpace.plugMat);//any kind of te will do, no need to eval energies
        TESampleSet tss = new TESampleSet(te);//draws actual samples
        te.setTrainingSamples(tss);
        for(RCTuple tup : XMinMatrix.allTuplesForFit()){
            te.tryAddingTuple(tup);
        }
        
        ArrayList<SampleXMin> ans = new ArrayList<>();
        for(int[] samp : tss.getSamples()){
            ans.add(new SampleXMin(searchSpace, samp));
        }
        return ans;
    }
    
    
    void pruneAllButOneConf(){
        if(singleConf.length!=searchSpace.confSpace.numPos)
            throw new RuntimeException("ERROR: singleConf is wrong length");
        int numRCs[] = searchSpace.confSpace.getNumRCsAtPos();
        for(int pos=0; pos<singleConf.length; pos++){
            for(int rc=0; rc<numRCs[pos]; rc++){//prune rc if not in singleConf
                if(rc!=singleConf[pos])
                    searchSpace.pruneMat.setOneBody(pos, rc, true);
            }
        }
    }
    
    void initXMinMatrixOneConf(){
        //just to see what can be done let's allow all pairs to optimize
        XMinMatrix = new XMFMatrix(searchSpace.confSpace);
        for(int pos=0; pos<singleConf.length; pos++){
            XMinMatrix.initOneBody(pos, singleConf[pos]);
            for(int pos2=0; pos2<pos; pos2++)
                XMinMatrix.initPairwise(pos, singleConf[pos], pos2, singleConf[pos2]);
        }
    }
    
}
