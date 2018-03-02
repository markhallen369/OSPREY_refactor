/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.bibis;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import java.io.Serializable;
import java.util.ArrayList;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.linear.LinearConstraint;

/**
 *
 * The beak is the tool with which the ibis eats,
 * i.e., converts the original expansion into an integrated one
 * 
 * @author mhall44
 */
public class Beak implements Serializable {
    
    ConfSampleSet trainingSamples, CVSamples;
    PruningMatrix pruneMat;
    PolytopeMatrix plugMat;

    public Beak(PruningMatrix pruneMat, PolytopeMatrix plugMat) {
        this.pruneMat = pruneMat;
        this.plugMat = plugMat;
    }//will draw samples in doIntegration
    
    public EnergyModel doIntegration(EnergyModel input, ArrayList<IntegrableDOF> dofsToIntegrate){
        trainingSamples = generateSamples(input);
        CVSamples = generateSamples(input);
        
        trainingSamples.checkPLUGViolations(plugMat);//DEBUG!!!
        
        EnergyModel curE = input;
        for(IntegrableDOF dof : dofsToIntegrate){
            System.out.println("Fitting DOF: "+dof.description()+" Starting with "+curE.getNumParams()+" parameters");
            curE = integrateWrtDOF(curE,dof);
            System.out.flush();
        }
        return curE;
    }
    
    public SparseLinearEnergy refitModel(EnergyModel input, FeatureSet featSet){
        //this is to check the ability to fit reliably when we know a good model is possible, or to check if a simpler model if possible
        trainingSamples = generateSamples(input);
        CVSamples = generateSamples(input);
        
        System.out.println("FITTING SPARSELINEARENERGY TO ENERGYMODEL DIRECTLY.  Input number of parameters: "+input.getNumParams()
                + " Number of features for fit: "+featSet.numFeatures);
        ConfSampleSet relevantTrainingSamples = selectSamples(trainingSamples,input);
        EnergiedConfSampleSet trainingEnergies = relevantTrainingSamples.evalEnergies(input);
        
        double initCoeffs[] = null;//((SparseLinearEnergy)input).coeffs;//DEBUG!!!  null would be a better test, want to see what this does though
        SparseLinearEnergy newF = new SparseLinearEnergy(featSet,trainingEnergies,initCoeffs);
        
        System.out.print("TRAINING ");
        trainingEnergies.checkError(newF);
        
        //DEBUG!!!  might find a more useful way to validate
        ConfSampleSet relevantCVSamples = selectSamples(CVSamples,input);
        EnergiedConfSampleSet CVEnergies = relevantCVSamples.evalEnergies(input);
        System.out.print("VALIDATION ");
        CVEnergies.checkError(newF);
        
        return newF;
    }
    
    public EnergyModel integrateWrtDOF(EnergyModel input, IntegrableDOF dof){
        ConfSampleSet relevantTrainingSamples = selectSamples(trainingSamples,input);
        EnergiedConfSampleSet trainingEnergies = relevantTrainingSamples.integrateEnergies(dof, input);
        EnergyModel newF = input.integrateModel(dof, trainingEnergies);
        
        System.out.print("TRAINING ");
        trainingEnergies.checkError(newF);
        
        //DEBUG!!!  might find a more useful way to validate
        ConfSampleSet relevantCVSamples = selectSamples(CVSamples,input);
        EnergiedConfSampleSet CVEnergies = relevantCVSamples.integrateEnergies(dof, input);
        System.out.print("VALIDATION ");
        CVEnergies.checkError(newF);
        
        return newF;
    }
        
    
    private ConfSampleSet selectSamples(ConfSampleSet superset, EnergyModel f){
        //DEBUG!!
        //can ultimately just use the needed number of samples per term, but for now keep it simple
        return superset;
    }
    
    private ConfSampleSet generateSamples(EnergyModel input){
        //DEBUG!!!  both sample set and breadth of inputs handled may need changes
        //although, if start w/ EPIC this is an SLE so even if use non-SLE's after integration need
        //only support SLE's here
        System.out.println("GENERATING SAMPLE SET for model with "+input.getNumParams()+" parameters");
        ConfSampleSet ans;
        if(input instanceof SparseLinearEnergy)
            //ans = new BasicConfSampleSet( ((SparseLinearEnergy)input).featSet, plugMat, pruneMat );
            ans = new ContBatchSampleSet( ((SparseLinearEnergy)input).featSet, plugMat, pruneMat );
            //DEBUG!!
        else
            throw new RuntimeException("ERROR: don't know how to draw samples for "+input.getClass());
        System.out.println("GENERATED SAMPLE SET with "+ans.getNumSamples()+" samples");
        System.out.flush();
        return ans;
    }
    
    public static void main(String args[]){
        //Initial test...just try building a LUTE matrix from an EPIC matrix,
        //see how good the fits are
        
        //luteKStarTest();
        //System.exit(0);//DEBUG!!!
        
        //DEBUG!!!!!!
        /*RealVector vec = new ArrayRealVector( new double[] {182.7184381912, 0, -65.0325540209, 113.2892091236, 0, 54.4636158988, 183.7030120587, 66, -60.2038654726, -60.2038654726, -60.2038654726, -60.2038654726, -60.2038654726, 94, -9, -56, 166} );
        int doNum = 7;
        ArrayList<LinearConstraint> polytope = (ArrayList<LinearConstraint>)ObjectIO.readObject("POLYTOPE.dat", false);
        double[] q = BasicConfSampleSet.boundSingleDOF(vec, polytope, doNum);
        System.exit(0);*/
        
        
        ConfigFileParser cfp = new ConfigFileParser(new String[] {"-c","KStar.cfg","alkluiaacuii"});
        cfp.loadData();
        
        EPICMatrix epicMat = (EPICMatrix) ObjectIO.readObject("ibis.EPICMAT.dat", true);
        PolytopeMatrix plugMat = (PolytopeMatrix) ObjectIO.readObject("ibis.PLUGMAT.dat", true);
        PruningMatrix pruneMat = (PruningMatrix) ObjectIO.readObject("ibis.PRUNEMAT.dat", true);
        
        ensureDOFsMatch(plugMat.cSpace,epicMat.getConfSpace());
        
        Beak curvedBeak = new Beak(pruneMat,plugMat);
        
        SparseLinearEnergy initSLE = new SparseLinearEnergy(epicMat,pruneMat,plugMat);//DEBUG!!  should probs put in plugMtx for pruning?
        ArrayList<IntegrableDOF> dofsToIntegrate = new ArrayList<>();//just do cont integ at first...def try K* later though
        //for(DegreeOfFreedom dof : epicMat.getConfSpace().confDOFs)
        //    dofsToIntegrate.add(new ContIntegrableDOF(dof,plugMat));
        //DEBUG!!! it's best to do chi2's before chi1's etc so doing this backwards
        for(int dofNum=epicMat.getConfSpace().confDOFs.size()-1; dofNum>=0; dofNum--){
            dofsToIntegrate.add(new ContIntegrableDOF(epicMat.getConfSpace().confDOFs.get(dofNum), plugMat, dofNum));
        }
        
        //DEBUG!!
        //curvedBeak.refitModel(initSLE, initSLE.featSet);
        
        
        EnergyModel luteModel = curvedBeak.doIntegration(initSLE, dofsToIntegrate);
        
        
        //Let's go on to K*!!
        //ObjectIO.writeObject(luteModel, "LUTE_MODEL.dat");
        dofsToIntegrate.clear();
        /*
        EnergyModel luteModel = (EnergyModel) ObjectIO.readObject("LUTE_MODEL.dat", true);
        ArrayList<IntegrableDOF> dofsToIntegrate = new ArrayList<>();*/
        
        
        for(int pos=0; pos<epicMat.getConfSpace().numPos; pos++){
            dofsToIntegrate.add(new RCIntegrableDOF(pos, epicMat.getConfSpace().posFlex.get(pos).RCs));
        }
        luteModel = curvedBeak.doIntegration(luteModel, dofsToIntegrate);
        
        
        //A LUTE model is an SLE, so we can do this:
        EnergyMatrix luteMat = ((SparseLinearEnergy)luteModel).asEnergyMatrix();
        System.out.println("The ibis has completed its work.");
    }
    
    
    public static void luteKStarTest(){
        //Integrating a LUTE matrix
        
        ConfigFileParser cfp = new ConfigFileParser(new String[] {"-c","KStar.cfg","alkluiaacuii"});
        cfp.loadData();
        
        EnergyMatrix luteMat = (EnergyMatrix) ObjectIO.readObject("ibis.LUTEMAT.dat", true);
        PolytopeMatrix plugMat = (PolytopeMatrix) ObjectIO.readObject("ibis.PLUGMAT.dat", true);
        PruningMatrix pruneMat = (PruningMatrix) ObjectIO.readObject("ibis.PRUNEMAT.dat", true);
                
        ConfSpace confSpace = plugMat.cSpace;
        Beak curvedBeak = new Beak(pruneMat,plugMat);
        SparseLinearEnergy initSLE = new SparseLinearEnergy(luteMat,pruneMat,plugMat);//DEBUG!!  should probs put in plugMtx for pruning?

        ArrayList<IntegrableDOF> dofsToIntegrate = new ArrayList<>();
        for(int pos=0; pos<confSpace.numPos; pos++){
            dofsToIntegrate.add(new RCIntegrableDOF(pos, confSpace.posFlex.get(pos).RCs));
        }
        EnergyModel integModel = curvedBeak.doIntegration(initSLE, dofsToIntegrate);
        
        EnergyMatrix integMat = ((SparseLinearEnergy)integModel).asEnergyMatrix();
        System.out.println("The ibis has completed its work.");
    }
    
    
    
    static private void ensureDOFsMatch(ConfSpace cSpace1, ConfSpace cSpace2){
        if(cSpace1.confDOFs.size() != cSpace2.confDOFs.size())
            throw new RuntimeException("ERROR: Unexpected different numbers of DOFs");
        for(int d=0; d<cSpace1.confDOFs.size(); d++){
            if(!cSpace1.confDOFs.get(d).getName().equals(cSpace2.confDOFs.get(d).getName()))
                throw new RuntimeException("ERROR: DOF mismatch");
        }
    }
}
