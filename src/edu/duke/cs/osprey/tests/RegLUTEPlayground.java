/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.TupleEnumerator;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICEnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tupexp.ConfETupleExpander;
import edu.duke.cs.osprey.tupexp.TupExpChooser;
import edu.duke.cs.osprey.tupexp.TupleExpander;
import java.util.ArrayList;

/**
 *
 * For playing with regularized LUTE
 * 
 * 
 * @author mhall44
 */
public class RegLUTEPlayground {
    
    SearchProblem searchSpace;
    
    //settings
    double regLambda;
    boolean sidechainMinFirst;
    boolean regCATSOnly;
    boolean includeRegEnergy;
    

    public RegLUTEPlayground(SearchProblem searchSpace, double regLambda, 
            boolean sidechainMinFirst, boolean regCATSOnly, boolean includeRegEnergy) {
        this.searchSpace = searchSpace;
        this.regLambda = regLambda;
        this.sidechainMinFirst = sidechainMinFirst;
        this.regCATSOnly = regCATSOnly;
        this.includeRegEnergy = includeRegEnergy;
                
        searchSpace.useTupExpForSearch = true;//for testing conf search
    }
    
    
    //use args like for Main
    public static void main(String args[]){
        ConfigFileParser cfp = new ConfigFileParser(args);
        cfp.loadData();
        RegLUTEPlayground rlp = new RegLUTEPlayground(
                cfp.getSearchProblem(),
                cfp.getParams().getDouble("LUTEREGLAMBDA", 0),
                cfp.getParams().getBool("LUTESCMINFIRST", false),
                cfp.getParams().getBool("LUTEREGCATSONLY", true),
                cfp.getParams().getBool("LUTEINCLREGE", false)
        );
        rlp.testRegLUTE();
    }
    
    public void testRegLUTE(){
        //given a cfp-generated search problem, test LUTE matrix generation
        //with various parameters
        //only plug pruning because the problem cases are high-ival and maybe even non-pairwise
        searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace, Double.POSITIVE_INFINITY);
        searchSpace.plugMat = new PolytopeMatrix(searchSpace,true);
        searchSpace.loadEnergyMatrix();//The tuple enumerator wants this
        searchSpace.loadEPICMatrix();
        makeAndTestLUTEMatrix();
    }
    
    
    public void makeAndTestLUTEMatrix(){
        TupleExpander expander = makeTupleExpander();
        TupleEnumerator tupEnum = new TupleEnumerator(searchSpace.pruneMat,searchSpace.emat,
                searchSpace.plugMat, searchSpace.confSpace.numPos);
        TupExpChooser chooser = new TupExpChooser(expander, tupEnum);//make a chooser to choose what tuples will be in the expansion
            
        double curResid = chooser.calcPairwiseExpansion();//start simple...
        testLUTEMatrix(expander.getEnergyMatrix());
        double goalResid = 0.01;

        if(curResid > goalResid){//go to triples if needed
            System.out.println("EXPANDING PAIRWISE EXPANSION WITH STRONGLY PAIR-INTERACTING TRIPLES (2 PARTNERS)...");
            curResid = chooser.calcExpansionResTriples(2);
            testLUTEMatrix(expander.getEnergyMatrix());
        }
        if(curResid > goalResid){//go to 5 partners if still need better resid...
            System.out.println("EXPANDING EXPANSION WITH STRONGLY PAIR-INTERACTING TRIPLES (5 PARTNERS)...");
            curResid = chooser.calcExpansionResTriples(5);
            testLUTEMatrix(expander.getEnergyMatrix());
        }
        if(curResid > goalResid){
            System.out.println("WARNING: Desired LUTE residual threshold "+
                    goalResid+" not reached; best="+curResid);
        }
        
        expander.dumpSampleSets(searchSpace.name);
    }
    
    public void testLUTEMatrix(EnergyMatrix luteMat){
        //Plug the LUTE matrix into searchSpace and see what it can do
        searchSpace.tupExpEMat = luteMat;
        System.out.println("TESTING LUTE MATRIX");
        
        ConfTree tree = ConfTree.makeFull(searchSpace);
        for(int confNum=0; confNum<5; confNum++){
            ScoredConf conf = tree.nextConf();
            if(conf!=null){
                double minE = calcMinimizedEnergy(conf.getAssignments());
                System.out.print("Conf "+confNum+" score="+conf.getScore()+" minE="+minE+": ");
                for(int rc : conf.getAssignments())
                    System.out.print(rc+" ");
                System.out.println();
            }
        }
        
        //a reference conformation for 1CC8
        try {
            int[] conf1 = new int[] {5,7,7,5,0,7,4};
            double score = luteMat.confE(conf1);
            double minE = calcMinimizedEnergy(conf1);
            System.out.println("Score for 5775074: "+score+" minE: "+minE);
        }
        catch(Exception e){
            System.out.println("LOL this is not 1CC8");
        }
    }
    
    
    private TupleExpander makeTupleExpander(){
        return new ConfETupleExpander(searchSpace,null){
            @Override
            public double scoreAssignmentList(int[] assignmentList) {
                return calcMinimizedEnergy(assignmentList);
            }
        };
    }
    
    private double calcMinimizedEnergy(int[] assignmentList){
        //Based on SearchProblem.EPICMinimizedEnergy(assignmentList)
        RCTuple tup = new RCTuple(assignmentList);
        EPICEnergyFunction epicEfunc = searchSpace.epicMat.internalEnergyFunction(tup,true);
        
        MultiTermEnergyFunction efunc = new MultiTermEnergyFunction();
        efunc.addTerm(epicEfunc);
        RegEFunc regEfunc = new RegEFunc();
        efunc.addTerm(regEfunc);
        
        MoleculeModifierAndScorer objFcn = new MoleculeModifierAndScorer(efunc,searchSpace.epicMat.getConfSpace(),tup);
        objFcn.setEfunc(epicEfunc);//init for epic efunc
        regEfunc.init(objFcn);//special init needed!  for now just doing it here, later can change NeedsInit interface
        
        CCDMinimizer ccdMin = new CCDMinimizer(objFcn,false);
        ccdMin.useInitFixableDOFs = sidechainMinFirst;
        DoubleMatrix1D bestDOFVals = ccdMin.minimize().dofValues;
        double E = objFcn.getValue(bestDOFVals);
        
        if(!includeRegEnergy)//want only the real energy in the LUTE expansion, though reg can shift the minimum
            E -= regEfunc.getEnergy();//we're still at the optimal coords now
        
        return E;
    }
    
    
    private class RegEFunc implements EnergyFunction {
        //regularization energy

        DoubleMatrix1D curDOFVals = null;
        DoubleMatrix1D dofLambdas = null;
        DoubleMatrix1D center = null;
        
        private void init(MoleculeModifierAndScorer objFcn){
            curDOFVals = objFcn.getCurDOFVals();//hold a pointer here, like EPICEnergyFunction
            
            //Penalty will be regLambda for going to the edge of a voxel in any dimension
            DoubleMatrix1D voxConstr[] = objFcn.getConstraints();
            DoubleMatrix1D voxSize = voxConstr[1].copy().assign(voxConstr[0],Functions.minus);
            center = voxConstr[0].copy().assign(voxConstr[1],Functions.plus).assign(Functions.mult(0.5));
            dofLambdas = voxSize.assign(voxSize,Functions.mult).assign(Functions.inv).assign(Functions.mult(4.*regLambda));
            
            if(regCATSOnly){
                ArrayList<DegreeOfFreedom> dofs = objFcn.getDOFs();
                for(int d=0; d<dofs.size(); d++){
                    if(!(dofs.get(d) instanceof BBFreeDOF))
                        dofLambdas.set(d,0);
                }
            }
        }
        
        @Override
        public double getEnergy() {
            DoubleMatrix1D devsq = curDOFVals.copy().assign(center,Functions.minus);
            devsq.assign(devsq,Functions.mult);
            return devsq.zDotProduct(dofLambdas);
        }
    }
    
    
    //shall we also try exporting LUTE energies for other kinds of fitting outside LUTE???
    //Might have trouble getting into A* but could use in an IBIS/K* type scheme
    //Just set a flag in TupleExpander to dump its datasets!!! (w/ or w/o its own fits)
    //(and maybe adjust sample sizes to support higher-parameter fitting too....)
}
