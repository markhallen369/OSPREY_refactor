/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug.triage;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.bbfree.BBFreeBlock;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.RamachandranChecker;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.Relationship;

/**
 *
 * Modeling Ramachandran constraints as linear constraints in CATS space
 * This will make Ramachandran pruning an easy LP thing
 * 
 * @author mhall44
 */
public class RamaLinModel {
    
    ArrayList<LinearConstraint> constr;//constraint to define the local Ramachandran-allowed region
    ArrayList<Integer> constrBres;//which bres each constr is for
    BBFreeBlock bfb;
    int numFreeDOFs;
    double cutoff = 0.02;
    int[] plotNums;//which plot to check Rama against...
    
    
    public RamaLinModel(BBFreeBlock bfb, boolean buildModel, boolean verbose){
        this.bfb = bfb;
        numFreeDOFs = bfb.getDOFs().size();
        
        //let's keep plot nums approp for current res types DEBUG!!
        //assuming gly and pro pos won't change
        List<Residue> resList = bfb.getResidues();
        plotNums = new int[resList.size()];
        for(int bres=0; bres<resList.size(); bres++){
            if(resList.get(bres).template.name.equalsIgnoreCase("GLY"))
                plotNums[bres] = 0;
            else if(resList.get(bres).template.name.equalsIgnoreCase("PRO"))
                plotNums[bres] = 1;
            else 
                plotNums[bres] = 2;
            
            //let's do pre-pro DEBUG!!!  assuming no pro right afterwards
            if(bres<resList.size()-1){
                if(resList.get(bres+1).template.name.equalsIgnoreCase("PRO"))
                    plotNums[bres] = 3;
            }
        }
        
        if(buildModel)
            buildModel(verbose);
    }
    
    
    //build from tangency
    //start with a little main here, see how well it classifies 
    public static void main(String args[]){
        
        args = new String[] {"-c","KStar.cfg","findGMEC","System.cfg","DEE.cfg"};
        ConfigFileParser cfp = new ConfigFileParser(args);
        cfp.loadData();
        
        Molecule m = PDBFileReader.readPDBFile("2K04.pep7.pdb");
        //113-117 bb dihs are defined.  don't mess with first N, last C for now, they might be messed up
        //and they don't come in a pair.  
        int mNumRes = m.residues.size();
        List<Residue> bfbRes = m.residues.subList(mNumRes-7, mNumRes);
        
        //BBFreeBlock bfb = new BBFreeBlock(bfbRes, false);
        double voxLim[] = new double[2*bfbRes.size()+2];
        Arrays.fill(voxLim, 1);
        BBFreeBlock bfb = new BBFreeBlock(bfbRes, false, voxLim, true);
        
        RamaLinModel rlm = new RamaLinModel(bfb,true,true);//false,false
        //RamaLinModel rlm = new RamaLinModel(bfb,false,false);
        System.out.println("Rama densities for res: ");
        DoubleMatrix1D allZeroes = DoubleFactory1D.dense.make( rlm.numFreeDOFs );
        for(int bres=1; bres<=3; bres++){
            System.out.println(rlm.linInterpRamaDensity(allZeroes,bres));
            System.out.println(rlm.bresConstrSatisfied(allZeroes,bres));
        }
    }
    
    
    private void buildModel(boolean verbose){
        constr = new ArrayList<>();
        constrBres = new ArrayList<>();
        int numSamp = 10;
        for(int bres=1; bres<=3; bres++){
            for(int s=0; s<numSamp; s++){
                DoubleMatrix1D samp = bfb.randomValsInVoxel();//random in voxel
                samp = scaleToBoundary(samp,bres);
                if(samp!=null){
                    constr.add(linearizedBoundary(samp, bres));
                    constrBres.add(bres);
                    checkConvexityAtBoundary(samp, bres);//DEBUG!!!!
                    //lots of values above cutoff seen here!  indicates not close to convex
                    //This explains poor results (majority of favorable pts not classified as favorable)
                    //even in small vox
                    //This is going to make lions very difficult
                    //very little pruning results from steric (<50% of RCs, seemingly less for pairs)
                    //unlikely quickclash can do much better since woudl need tons of tuples 
                    //right on boundary of good/bad
                    //plus lions basically sets up the classic local np-hard protein design case:
                    //all kinds of different pairwise energy values
                    //this can only be avoided if there's lots of pruning
                    //pruning by negative space/lack of contacts looks intractable because so many
                    //different contacts are possible for any atom...pruning requires
                    //(at contacts at2 OR at3...) AND (at' contacts at2 OR at3...)
                    //so setting up a voxel space w/ big bb range is looking intractable.  
                    //maybe manual design could help here but leans heavily on energy function.  
                }
            }
        }
        if(verbose){
            System.out.println("Picked "+constr.size()+" constraints out of "+numSamp+" samples");
            checkModel();
        }
    }
    
    
    private void checkModel(){
        System.out.println("Checking model: ");
        int numSamp = 10;
        int numCorrectFavorable = 0;
        int numCorrectUnfavorable = 0;
        int numIncorrectFavorable = 0;
        int numIncorrectUnfavorable = 0;
        for(int bres=1; bres<=3; bres++){
            for(int s=0; s<numSamp; s++){
                DoubleMatrix1D samp = bfb.randomValsInVoxel();
                boolean favorable = (linInterpRamaDensity(samp,bres)>cutoff);
                boolean linFavorable = bresConstrSatisfied(samp,bres);
                if(favorable&&linFavorable)
                    numCorrectFavorable++;
                if( (!favorable) && (!linFavorable) )
                    numCorrectUnfavorable++;
                if( (!favorable) && linFavorable )
                    numIncorrectUnfavorable++;
                else//favorable && !linFavorable
                    numIncorrectFavorable++;
            }
        }
        System.out.println("Tried "+(3*numSamp)+" samples");
        System.out.println("Favorable samples correct: "+numCorrectFavorable+" Incorrect: "+numIncorrectFavorable);
        System.out.println("Unfavorable samples correct: "+numCorrectUnfavorable+" Incorrect: "+numIncorrectUnfavorable);
    }
    
    boolean bresConstrSatisfied(DoubleMatrix1D x, int bres){
        for(int cnum=0; cnum<constr.size(); cnum++){
            if(constrBres.get(cnum)==bres){
                LinearConstraint c = constr.get(cnum);
                DoubleMatrix1D g = DoubleFactory1D.dense.make(c.getCoefficients().toArray());
                if(c.getRelationship()!=Relationship.GEQ)
                    throw new RuntimeException("ERROR: Expected GEQ constr");
                if( g.zDotProduct(x)<c.getValue() )
                    return false;
            }
        }
        return true;
    }
    
    private LinearConstraint linearizedBoundary(DoubleMatrix1D x, int bres){
        //A constraint on the free DOFs, tangent to the surface of the Ramachandran-favored region at x
        //grad(den) dot (DOFVals-x) > 0
        DoubleMatrix1D gradDen = linInterpRamaDenGrad(x, bres);
        return new LinearConstraint(gradDen.toArray(), Relationship.GEQ, gradDen.zDotProduct(x));
    }
    
    private DoubleMatrix1D linInterpRamaDenGrad(DoubleMatrix1D bbDOFVals, int bres){
        //Let's keep this simple for now
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(numFreeDOFs);
        double step = 1e-5;
        for(int b=0; b<numFreeDOFs; b++){
            DoubleMatrix1D x = bbDOFVals.copy();
            x.set(b, bbDOFVals.get(b)+step);
            double upVal = linInterpRamaDensity(x,bres);
            x.set(b, bbDOFVals.get(b)-step);
            double downVal = linInterpRamaDensity(x,bres);
            ans.set(b, (upVal-downVal)/(2*step) );
        }
        return ans;
    }
    
    private double linInterpRamaDensity(DoubleMatrix1D bbDOFVals, int bres){
        bfb.setDOFs(bbDOFVals);
        double phiPsi[] = RamachandranChecker.getPhiPsi(bfb.getResidues().get(bres));
        
        //DEBUG!!!!!  want to get rid of alpha bit since 3 bres are all in beta zone
        if(phiPsi[1]<45)
            return 0;
        
        return RamachandranChecker.getInstance().getLinInterpDensity(phiPsi[0], phiPsi[1], plotNums[bres]);
    }
    
    
    private DoubleMatrix1D scaleToBoundary(DoubleMatrix1D samp, int bres){
        //use bisection to scale samp to be at the boundary of the Rama-allowed region
        //DEBUG!!!  should really just limit to current island
        double [][] freeDOFVoxel = bfb.getFreeDOFVoxel();
        double scalingUB = Double.POSITIVE_INFINITY;//upper bound on scaling
        for(int d=0; d<numFreeDOFs; d++){
            if(samp.get(d)>0)
                scalingUB = Math.min(scalingUB, freeDOFVoxel[1][d]/samp.get(d));
            else
                scalingUB = Math.min(scalingUB, freeDOFVoxel[0][d]/samp.get(d));
        }
        
        double scalingLB = 0;
        DoubleMatrix1D x = samp.copy();
        x.assign(Functions.mult(scalingUB));
        if(linInterpRamaDensity(x,bres)>cutoff)
            return null;//already in range
        
        while(scalingUB-scalingLB>0.0001){
            double mid = (scalingUB+scalingLB)/2;
            x.assign(samp);
            x.assign(Functions.mult(mid));
            if(linInterpRamaDensity(x,bres)>cutoff)
                scalingLB = mid;
            else
                scalingUB = mid;
        }
        
        x.assign(samp);
        x.assign(Functions.mult((scalingUB+scalingLB)/2));
        return x;
    }
    
    
    private void checkConvexityAtBoundary(DoubleMatrix1D x, int bres){
        double dx = 0.01;
        DoubleMatrix1D inDir = linInterpRamaDenGrad(x,bres);
        inDir.assign(Functions.mult(1./Math.sqrt(Algebra.DEFAULT.norm2(x))));
        DoubleMatrix1D z = x.copy();
        System.out.println("Checking convexity at boundary.  cur val: "+linInterpRamaDensity(z,bres));
        z.assign(inDir,Functions.plusMult(dx));
        System.out.println("Val stepping inward: "+linInterpRamaDensity(z,bres));
        z.assign(x);
        z.assign(inDir,Functions.plusMult(-dx));
        System.out.println("Outward: "+linInterpRamaDensity(z,bres));
        for(int s=0; s<5; s++){
            z.assign(x);
            DoubleMatrix1D q = DoubleFactory1D.dense.random(x.size());
            q.assign( inDir, Functions.plusMult(-q.zDotProduct(inDir)) );
            q.assign(Functions.mult(dx/Math.sqrt(Algebra.DEFAULT.norm2(x))));
            z.assign(q, Functions.plus);
            System.out.println("Val stepping perpendicularly: "+linInterpRamaDensity(z,bres));
        }
    }
    
}
