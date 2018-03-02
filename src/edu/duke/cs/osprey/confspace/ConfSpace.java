/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import java.io.Serializable;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import com.joptimizer.functions.LinearMultivariateRealFunction;
import edu.duke.cs.osprey.bbfree.BBFreeBlock;
import edu.duke.cs.osprey.bbfree.CATSSettings;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.dof.MoveableStrand;
import edu.duke.cs.osprey.dof.ProlinePucker;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.dof.StrandRotation;
import edu.duke.cs.osprey.dof.StrandTranslation;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.dof.deeper.perts.Perturbation;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.KSTermini;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.minimization.SQPMinimizer;
import edu.duke.cs.osprey.plug.LPChecks;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.StringParsing;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.Relationship;

/**
 *
 * @author mhall44
 */
public class ConfSpace implements Serializable {
    //This class represents the conformational search space for a design
    //used for GMEC-based design, K*, or anything else we want to do with a conformational space
    //This class just defines the conformational space itself (the molecule + all kinds of flexibility
    //and possible mutations, and how these are represented as RCs, etc.)
    //This class can be put in a SearchProblem to add annotations like what RCs are pruned,
    //what their pairwise energies are, etc.  
    
	private static final long serialVersionUID = 6414329117813457771L;    
 
	public Molecule m;
    //The molecule will be composed of residues. 
    //It will have one set of coordinates, which are stored in the residues to make mutation
    //and pairwise energy computation easy (no need for complicated partial arrays, subtracting off
    //template energies, etc., and mutation will be very quick and will only require
    //AMBER reinitialization for the affected residue)
        
    //If we need to keep a copy of, say, the original coordinates, we can have a second molecule origMolec
    //for that
    //Once loaded, the molecule can only be changed by functions overriding DegreeOfFreedom.applyValue
    //this will give us a uniform framework for applying conformational and sequence changes
    //and each DOF will store its value, so we can easily and unambiguously look up the state of the
    //molecule at any given time
    
    
    public ArrayList<DegreeOfFreedom> confDOFs = new ArrayList<>();
    //list of conformational degrees of freedom
    
    public ArrayList<ResidueTypeDOF> mutDOFs = new ArrayList<>();
    //list of sequence degrees of freedom (residue types for the mutable residues)
    
    
    public ArrayList<PositionConfSpace> posFlex = new ArrayList<>();
    //defines the flexible positions and what RCs they have
    //generally each position is a residue, but it could be more than one ("super-residue" with "super-RCs")

    public ArrayList<String> flexibleRes;
    public int numPos;//number of flexible positions
    
    public boolean useEllipses = false;
    
    /** initialize a new conformational space, defining all its flexibility
    /*   we use one residue per position here
     * 
     * @param PDBFile the structure to read from
     * @param flexibleRes list of residue numbers to be made flexible (as in PDB file)
     * @param allowedAAs list of allowed residue types at each flexible position
     * @param addWT whether to add wild-type to the allowed AA types 
     * @param contSCFlex means allow continuous sidechain flexibility
     * @param dset DEEPer Settings
     * @param moveableStrands ... ? 
     * @param freeBBZones ...? 
     * @param ellipses model ellipses
     * @param addWTRots add the wild-type 'rotamers'
     */
    public ConfSpace(String PDBFile, ArrayList<String> flexibleRes, ArrayList<ArrayList<String>> allowedAAs, 
            boolean addWT, ArrayList<String> wtRotOnlyRes, boolean contSCFlex, DEEPerSettings dset, ArrayList<String[]> moveableStrands, 
            CATSSettings catsSettings, boolean ellipses, boolean addWTRots, KSTermini termini){
    
    	useEllipses = ellipses;  	
    	this.flexibleRes = flexibleRes;
        numPos = flexibleRes.size();
        
        //read the structure and assign templates, deleting unassignable res...
        m = PDBFileReader.readPDBFile(PDBFile, termini);
        
        // before making any structure changes, capture the wt rots if needed
        List<ResidueTemplate> wtRots = new ArrayList<>(Collections.nCopies(numPos, null));
        if (addWTRots) {
        	for (int i=0; i<numPos; i++) {
        		// TODO: support alternate conformations?
        		Residue res = m.getResByPDBResNumber(flexibleRes.get(i));
        		wtRots.set(i, ResidueTemplate.makeFromResidueConfs(res));
        	}
        }
        
        //Make all the degrees of freedom
        //start with proline puckers (added to res)
        makeRingPuckers(allowedAAs, flexibleRes);
        
        //now the single-res degrees of freedom (mutations, sidechain dihedrals)
        ArrayList<ArrayList<DegreeOfFreedom>> singleResDOFs = new ArrayList<>();

        for(int pos=0; pos<numPos; pos++){
            Residue res = m.getResByPDBResNumber( flexibleRes.get(pos) );
            
            //at this point, m has all wild-type residues, so just see what res is now
            String wtName = res.template.name;
            if( ! StringParsing.containsIgnoreCase(allowedAAs.get(pos), wtName) ){//wtName not currently in allowedAAs
                if(addWT || allowedAAs.get(pos).isEmpty())//It should be (addWT or blank AA type)
                    allowedAAs.get(pos).add(wtName);
                else//it should not be...make sure wt rots aren't included
                    wtRots.set(pos, null);
            }
            
            ArrayList<DegreeOfFreedom> resDOFs = mutablePosDOFs(res,allowedAAs.get(pos));//add mutation and dihedral confDOFs
                        
            ResidueTypeDOF resMutDOF = (ResidueTypeDOF)resDOFs.remove(0);//first mutable pos DOF is the mutation-type DOF
 
            mutDOFs.add(resMutDOF);

            singleResDOFs.add(resDOFs);
        }
        
        //now rigid-body strand motions...
        ArrayList<DegreeOfFreedom> strandDOFs = strandMotionDOFs(moveableStrands,flexibleRes);
        confDOFs.addAll(strandDOFs);
        
        //...and perturbations
        //standardize conformations first since we'll record initial resBBState here
        standardizeMutatableRes(allowedAAs, flexibleRes);
        
        ArrayList<Perturbation> perts = dset.makePerturbations(m);//will make pert block here
        confDOFs.addAll(perts);
        
        ArrayList<BBFreeBlock> bfbList = getBBFreeBlocks(catsSettings,flexibleRes);
        for(BBFreeBlock bfb : bfbList)
            confDOFs.addAll( bfb.getDOFs() );
        
        //OK now make RCs using these DOFs
        for(int pos=0; pos<numPos; pos++){
            Residue res = m.getResByPDBResNumber( flexibleRes.get(pos) );
            
            ArrayList<DegreeOfFreedom> resDOFs = singleResDOFs.get(pos);
            ArrayList<DegreeOfFreedom> resStrandDOFs = strandDOFsForRes(res,strandDOFs);//and check that all res flex on moving strands!
            
            BBFreeBlock curBFB = getCurBFB(bfbList,res);
            
            boolean wtRotOnly = wtRotOnlyRes.contains(flexibleRes.get(pos));
            if(wtRotOnly && !(wtRots.get(pos)!=null&&allowedAAs.get(pos).size()==1) )
                throw new RuntimeException("ERROR: WT rot only on but residue not single AA type with wild-type rotamer");
            
            PositionConfSpace rcs = new PositionConfSpace(pos, res, resDOFs, allowedAAs.get(pos), contSCFlex,
                    resStrandDOFs, perts, dset.getPertIntervals(), dset.getPertStates(pos), curBFB, useEllipses, 
                    wtRots.get(pos), wtRotOnly);
            posFlex.add(rcs);
                        
            if (useEllipses) {
            	confDOFs.addAll(rcs.getEllipsoidalArray());
            } else {
            	confDOFs.addAll(resDOFs);
            }
            
        }
        
        
        //finally number all the confDOFs to facilitate comparisons
        for(int dofNum=0; dofNum<confDOFs.size(); dofNum++)
            confDOFs.get(dofNum).confDOFNum = dofNum;
        
        //DEBUG!!!
        /*PDBFileWriter.writePDBFile(m, "STRUCT1.pdb");
        perts.get(0).apply(5);
        PDBFileWriter.writePDBFile(m, "STRUCT2.pdb");
        perts.get(0).apply(0);
        PDBFileWriter.writePDBFile(m, "STRUCT3.pdb");*/
    }
    
    public ConfSpace(ConfSpace other) {
    	// just make a shallow copy
    	this.m = other.m;
    	this.confDOFs = new ArrayList<>(other.confDOFs);
    	this.mutDOFs = new ArrayList<>(other.mutDOFs);
		this.posFlex = new ArrayList<>(other.posFlex);
		this.flexibleRes = new ArrayList<>(other.flexibleRes);
    	this.numPos = other.numPos;
    	this.useEllipses = other.useEllipses;
    }
    
    /*private ArrayList<BBFreeBlock> getBBFreeBlocks(CATSSettings catsSettings, ArrayList<String> flexibleRes){
        //create a BFB for each (start res, end res) pair.  PDB residue numbers provided.  
        ArrayList<BBFreeBlock> ans = new ArrayList<>();
        
        for(String[] termini : catsSettings.freeBBZones){
            ArrayList<Residue> curBFBRes = resListFromTermini(termini, flexibleRes);
            BBFreeBlock bfb = new BBFreeBlock(curBFBRes,true);
            ans.add(bfb);
        }
        
        return ans;
    }*/
    private ArrayList<BBFreeBlock> getBBFreeBlocks(CATSSettings catsSettings, ArrayList<String> flexibleRes){
        return catsSettings.buildBBFreeBlocks(flexibleRes, m);
    }
    
    private BBFreeBlock getCurBFB(ArrayList<BBFreeBlock> bfbList, Residue res){
        //If res is in one of the BFB's, return that BFB; else return null;
        for(BBFreeBlock bfb : bfbList){
            if(bfb.getResidues().contains(res)){//current res is in this BFB
                return bfb;
            }
        }
        
        return null;//not in a BFB
    }
    
    
    private void makeRingPuckers(ArrayList<ArrayList<String>> allowedAAs, ArrayList<String> flexibleRes){
        //For each residue, create a Pro pucker DOF if Pro is allowed
        
        for(int pos=0; pos<numPos; pos++){
            Residue res = m.getResByPDBResNumber( flexibleRes.get(pos) );
            
            if(res.template.name.equalsIgnoreCase("PRO"))//the reisdue is a proline
                res.pucker = new ProlinePucker(res);
            else {//see if it can mutate to a proline
                for(String AAType : allowedAAs.get(pos)){
                    if(AAType.equalsIgnoreCase("PRO")){
                        res.pucker = new ProlinePucker(res);
                        break;
                    }
                }
            }
            
            if(res.pucker != null)
                confDOFs.add(res.pucker);
        }
    }
    
    
   private ArrayList<Residue> resListFromTermini(String[] termini, ArrayList<String> flexibleRes){ 
       return m.resListFromTermini(termini, flexibleRes);
   }
   
    
    private ArrayList<DegreeOfFreedom> strandMotionDOFs(ArrayList<String[]> moveableStrands,
            ArrayList<String> flexibleRes){
        //Generate all the strand rotation/translation DOFs,
        //given the termini of the moving strands
        
        ArrayList<DegreeOfFreedom> ans = new ArrayList<>();
        
        for(String[] termini : moveableStrands){
            ArrayList<Residue> curStrandRes = resListFromTermini(termini, flexibleRes);
            MoveableStrand str = new MoveableStrand(curStrandRes);
            ans.addAll(str.getDOFs());
        }
        
        return ans;
    }
    
    
    private ArrayList<DegreeOfFreedom> strandDOFsForRes(Residue res, ArrayList<DegreeOfFreedom> strandDOFs){
        //List the strandDOFs that move res
        ArrayList<DegreeOfFreedom> ans = new ArrayList<>();
        
        for(DegreeOfFreedom dof : strandDOFs){
            MoveableStrand str;
            if(dof instanceof StrandRotation)
                str = ((StrandRotation)dof).getMoveableStrand();
            else//must be translation
                str = ((StrandTranslation)dof).getMoveableStrand();
            
            if(str.getResidues().contains(res))
                ans.add(dof);
        }
        
        return ans;
    }
    
    
    private void standardizeMutatableRes(ArrayList<ArrayList<String>> allowedAAs, ArrayList<String> flexibleRes){
        //"mutate" all mutatable residues to the template version of their residue type
        //this ensures that the non-adjustable DOFs (bond angles, etc.) will be as in the template
        //(for consistency purposes)
        for(int pos=0; pos<numPos; pos++){
            Residue res = m.getResByPDBResNumber( flexibleRes.get(pos) );
            String resName = res.template.name;
            
            if(EnvironmentVars.resTemplates.getTemplateForMutation(resName, res, false) != null){
                //mutation to current residue type is possible, i.e., this position is mutatable
            //if(HardCodedResidueInfo.canMutateTo(res.template)){//this messes up for N-term mutations
                
                ResidueTypeDOF mutDOF = mutDOFs.get(pos);
                
                for(String allowedAA : allowedAAs.get(pos)){
                    if(allowedAA.equalsIgnoreCase("PRO")){
                        //mutating to and from PRO alters the sidechain a little, including idealization
                        //(because the sidechain connects in two places, it behaves a little different)
                        //Once we mutate to Pro and back we have a consistent conformation
                        mutDOF.mutateTo("PRO");
                        break;
                    }
                }                
                
                // AAO 2016: mutation assumes residue is an amino acid. throws an exception otherwise
                if(HardCodedResidueInfo.hasAminoAcidBB(res) && !res.fullName.startsWith("FOL")) {
                	mutDOF.mutateTo(resName);
                }
            }
        }
    }
    
    static ArrayList<DegreeOfFreedom> mutablePosDOFs(Residue res, ArrayList<String> allowedAAs){
        //mutation and dihedral confDOFs for the specified position
        
        ArrayList<DegreeOfFreedom> ans = new ArrayList<>();
        //assuming for now that this is a single-residue position...if multiple just add the following confDOFs for each residue
        
        //mutation DOF
        ans.add(new ResidueTypeDOF(res));
        
        int maxNumDihedrals = 0;//we need to create enough dihedral confDOFs for the allowed AA type with the most dihedrals
        
        for(String AAType : allowedAAs)
            maxNumDihedrals = Math.max( maxNumDihedrals, EnvironmentVars.resTemplates.numDihedralsForResType(AAType) );
        
        for(int dih=0; dih<maxNumDihedrals; dih++)
            ans.add( new FreeDihedral(res,dih) );
                
        if(res.pucker!=null)
            ans.add(res.pucker);
        
        return ans;
    }
    
    
    public double minimizeEnergy(int[] conf, EnergyFunction efunc, String outputPDBFile){
        //let's keep CCD default for now
        return minimizeEnergy(conf, efunc, outputPDBFile, false);
    }
    
    public double minimizeEnergy(int[] conf, EnergyFunction efunc, String outputPDBFile, boolean useSQP){
        return minimizeEnergy(conf, efunc, efunc, outputPDBFile, useSQP);
    }
    
    public double minimizeEnergy(int[] conf, EnergyFunction minEfunc, EnergyFunction evalEfunc, 
            String outputPDBFile, boolean useSQP){
        //minimize the energy of a conformation, within the DOF bounds indicated by conf (a list of RCs)
        //return the minimized energy
        //if outputPDBFile isn't null, then output the minimized conformation to that file
        //This version uses minEFunc to minimize and then evalEFunc to evaluate energy at the minimum
        
        RCTuple RCs = new RCTuple(conf);
        MoleculeModifierAndScorer energyForMin = new MoleculeModifierAndScorer(minEfunc,this,RCs);
        MoleculeModifierAndScorer energyForEval = new MoleculeModifierAndScorer(evalEfunc,this,RCs);
        
        DoubleMatrix1D optDOFVals;
        
        if(energyForMin.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
            Minimizer min;
            if(useSQP)
                min = new SQPMinimizer(energyForMin,null);//this function doesn't add additional linear constraints
            else
                min = new CCDMinimizer(energyForMin,false);
            //with the generic objective function interface we can easily include other minimizers though

            
            //DEBUG!!!!!  Timing pre-minimization without PB
            /*ArrayList<EnergyFunction> terms = ((MultiTermEnergyFunction)energy.getEfunc()).getTerms();
            ArrayList<Double> coeffs = ((MultiTermEnergyFunction)energy.getEfunc()).getCoeffs();
            if( terms.get(terms.size()-1) instanceof PoissonBoltzmannEnergy ){
                PoissonBoltzmannEnergy pbe = (PoissonBoltzmannEnergy)terms.remove(terms.size()-1);
                double pbcoeff = coeffs.remove(terms.size()-1);
                long startTime = System.currentTimeMillis();
                DoubleMatrix1D startDOFVals = min.minimize();
                long time1 = System.currentTimeMillis();
                terms.add(pbe);
                coeffs.add(pbcoeff);
                ((CCDMinimizer)min).setInitVals(startDOFVals);
                optDOFVals = min.minimize();
            }
            else//NON-DEBUG!*/
                optDOFVals = min.minimize().dofValues;
        }
        else//molecule is already in the right, rigid conformation
            optDOFVals = DoubleFactory1D.dense.make(0);
        
        double minE = energyForEval.getValue(optDOFVals);//this will put m into the minimized conformation
        
        
        
        //DEBUG!!!!!
        /*DoubleMatrix1D xEPIC = DoubleFactory1D.dense.make(new double[]
        {-57.969192, 166.735048, -0.066281, -0.229844, 0.245356, -0.203498, 0.140798, -0.071649, -0.304984, -0.455166, -183.708158, 51, -58.074254, -64.066479, -79, -56.170042, 87.783936, 9, -68.044382, -69.201669, -31, -69.036891, 168.647869}
        );
        System.out.println("E at EPIC min: "+energy.getValue(xEPIC));
        
         DoubleMatrix1D xreg = DoubleFactory1D.dense.make(new double[]
        {-57.916999, 166.69074, -0.062106, -0.232096, 0.23491, -0.201836, 0.147409, -0.074318, -0.293336, -0.455166, -183.76415, 51, -58.090563, -64.197299, -79, -56.228339, 87.764565, 9, -67.969832, -69.058391, -31, -69.092538, 168.669794}
        );
        System.out.println("E at reg min: "+energy.getValue(xreg));*/
                
        //DEBUG!!!!
        
        /*System.out.println("CCD opt DOF vals: "+optDOFVals);
        System.out.println("CCD min E: "+minE);
             
                Minimizer min2 = new SQPMinimizer(energy,null);
                Minimizer.Result sqpResult = min2.minimize();
                System.out.println("SQP opt DOF vals: "+sqpResult.dofValues);
                System.out.println("SQP min E: "+sqpResult.energy);

        System.out.println("Here are the energies along the line segment from CCD min to SQP min: ");
        int numSteps=20;
        for(int step=0; step<=20; step++){
            double stepFrac = step*1.0/numSteps;
            DoubleMatrix1D x = optDOFVals.copy().assign(Functions.mult(1-stepFrac));
            x.assign(sqpResult.dofValues, Functions.plusMult(stepFrac));
            System.out.println(energy.getValue(x));
        }
        System.out.println("Line segment energies done.");*/
                
        //System.out.println("CCD opt vals in polytope :"+pmat.isPointFeasible(conf, optDOFVals, energy.getDOFs()));
        //System.out.println("SQP opt vals in polytope :"+pmat.isPointFeasible(conf, sqpResult.dofValues, energy.getDOFs()));

        
        //DEBUG!!!!!!
        /*DoubleMatrix1D x = optDOFVals.copy();
        int numGridpts = 100;
        int graphDims[] = new int[] {9,10};
        double lb[] = new double[2];
        double step[] = new double[2];
        for(int a=0; a<2; a++){
            lb[a] = energy.getConstraints()[0].get(graphDims[a]);
            step[a] = ( energy.getConstraints()[1].get(graphDims[a]) - lb[a] ) / (numGridpts-1);
        }
        System.out.println("Graphing energy values in voxel wrt graphDims");
        System.out.println("DOF values at minimum energy: "+optDOFVals.get(graphDims[0])+" "+optDOFVals.get(graphDims[1]));
        System.out.println("Minimum energy: "+minE);
        for(int g1=0; g1<numGridpts; g1++){
            for(int g2=0; g2<numGridpts; g2++){
                x.set(graphDims[0], lb[0]+step[0]*g1);
                x.set(graphDims[1], lb[1]+step[1]*g2);
                double E = energy.getValue(x);
                System.out.println(x.get(graphDims[0])+" "+x.get(graphDims[1])+" "+E);
            }
        }
        System.out.println("Done graphing.");*/
        //DEBUG!!!
        
        
        if(outputPDBFile!=null)
            PDBFileWriter.writePDBFile(m, outputPDBFile, minE);
        
        return minE;
    }
    
    
    
    public double minimizeEnergyWithPLUG(int[] conf, EnergyFunction efunc, String outputPDBFile, PolytopeMatrix pmat){
        //minimize the energy of a conformation, within the DOF bounds indicated by conf (a list of RCs)
        //return the minimized energy
        //if outputPDBFile isn't null, then output the minimized conformation to that file
        
        //System.out.println("MINIMIZING ENERGY WITH PLUG");
        
        RCTuple RCs = new RCTuple(conf);
        MoleculeModifierAndScorer energy = new MoleculeModifierAndScorer(efunc,this,RCs);
             
        ArrayList<LinearConstraint> constrList = pmat.getFullStericPolytope(RCs, energy.getDOFs());
        
        LinearMultivariateRealFunction constr[] = new LinearMultivariateRealFunction[constrList.size()];
        for(int c=0; c<constrList.size(); c++){
            constr[c] = LPChecks.toLinearMultivariateRealFunction(constrList.get(c));
        }
        
        
                Minimizer min2 = new SQPMinimizer(energy,constr);
                Minimizer.Result sqpResult = min2.minimize();
                //System.out.println("SQP opt DOF vals: "+sqpResult.dofValues);
                //System.out.println("SQP min E: "+sqpResult.energy);

        //System.out.println("SQP opt vals in polytope :"+pmat.isPointFeasible(conf, sqpResult.dofValues, energy.getDOFs()));
                if(sqpResult==null){
                    if(outputPDBFile!=null)
                        System.out.println("Warning: No constraint-satisfying conformation found to put in "+outputPDBFile);
                    return Double.POSITIVE_INFINITY;
                }
                
                
        
        if(outputPDBFile!=null)
            PDBFileWriter.writePDBFile(m, outputPDBFile, sqpResult.energy);
        
        //System.out.println("DONE MINIMIZING ENERGY WITH PLUG");
        
        return sqpResult.energy;
    }
    
    
    
    
	public MultiTermEnergyFunction getDecomposedMinimizedEnergy(int[] conf, EnergyFunction efunc, String outputPDBFile){
		//minimize the energy of a conformation, within the DOF bounds indicated by conf (a list of RCs)
		//return the minimized energy
		//if outputPDBFile isn't null, then output the minimized conformation to that file

		RCTuple RCs = new RCTuple(conf);
		MoleculeModifierAndScorer energy = new MoleculeModifierAndScorer(efunc,this,RCs);

		DoubleMatrix1D optDOFVals;

		if(energy.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
			Minimizer min = new CCDMinimizer(energy,false);
			optDOFVals = min.minimize().dofValues;
		}
		else//molecule is already in the right, rigid conformation
			optDOFVals = DoubleFactory1D.dense.make(0);

		double minE = energy.getValue(optDOFVals);//this will put m into the minimized conformation
		
		if(outputPDBFile!=null)
			PDBFileWriter.writePDBFile(m, outputPDBFile, minE);
		
		return (MultiTermEnergyFunction)energy.getEfunc();
	}
    
    
    public int[] getNumRCsAtPos(){
        //list number of RCs at each position
        int[] numAllowed = new int[numPos];
        
        for(int pos=0; pos<numPos; pos++)
            numAllowed[pos] = posFlex.get(pos).RCs.size();
        
        return numAllowed;
    }
    
    
    public double getRCResEntropy(int pos, int rc){
        String resType = posFlex.get(pos).RCs.get(rc).AAType;
        double resEntropy = EnvironmentVars.resTemplates.getResEntropy(resType);
        return resEntropy;
    }
    
    public double getConfResEntropy(int[] conf){
        double ans = 0;
        
        for(int pos=0; pos<numPos; pos++)
            ans += getRCResEntropy(pos, conf[pos]);
        
        return ans;
    }


	public BigInteger getNumConformations() {
		BigInteger count = BigInteger.valueOf(1);
		for (int pos=0; pos<numPos; pos++) {
			count = count.multiply(BigInteger.valueOf(posFlex.get(pos).RCs.size()));
		}
		return count;
	}
        
    public ArrayList<DegreeOfFreedom> listAllDOFs(){
        ArrayList<DegreeOfFreedom> ans = new ArrayList<>();
        ans.addAll(mutDOFs);
        ans.addAll(confDOFs);
        return ans;
    }
    
    
    public int getDesignIndex(Residue res){
        for(int pos=0; pos<numPos; pos++){
            if(posFlex.get(pos).res==res)
                return pos;
        }
        return -1;
    }
    
    
    /*DoubleMatrix1D[] convertConfToDOFBounds(int[] conf){
        //Argument: RC assignments for all the flexible residues (RCs defined in resFlex)
        //return bounds (lower bounds, upper bounds) for all the degrees of freedom in the system
        RCTuple fullTuple = new RCTuple(conf);
        
        return RCTupleDOFBounds(fullTuple,confDOFs);
        //PULL FROM RC OBJECTS
    }
    
    public DoubleMatrix1D[] RCTupleDOFBounds(RCTuple tuple, ArrayList<DegreeOfFreedom> dofList){
        //Compute the bounds on the confDOFs in DOFlist imposed by the RCs in tuple
        //Create two vectors, vector 0 giving the lower bound for each DOF in dofList,
        //vector 1 giving the upper bound
        //can leave entries at 0 for AA types
        
        //WE SHOUDL TAKE MUT OUT OF DOFS ITS FUNDAMENTALLY DIFFERENT...NO DOUBLE REP, NOT CONF CHANGE
        
        //THIS IS FOR MAKING MOLEC E OBJ FUNCTION
        //MAKE IT TAKE RCs, RETURN ALL CONTINUOUS DOFS AND THEIR BOUNDS
        
        //ACTUALLY WAIT LETS MAKE THIS MOLEC E OBJ FUNCTION CONSTRUCTOR
    }*/
    
    //pairwise minimization will be similar...just only use degrees of freedom affecting the residue pair
    //and the objective function will represent the pairwise energy between the residues
    
    
    
    
}
