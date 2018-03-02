/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.bbfree;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.dof.DOFBlock;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.deeper.GenChi1Calc;
import edu.duke.cs.osprey.dof.deeper.SidechainIdealizer;
import edu.duke.cs.osprey.ematrix.epic.SeriesFitter;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.RigidBodyMotion;
import edu.duke.cs.osprey.tools.RotationMatrix;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

/**
 *
 * This is a set of consecutive residues whose backbones can move freely,
 * subject to distance, angle, and omega constraints
 * If fixAnchors, then the freely moving zone is anchored fix CA's at either end
 * Else, then the free DOFs account for the BB conf from one anchor CA to the next,
 * and then anything beyond the anchors is placed onto that
 * 
 * @author mhall44
 */
public class BBFreeBlock implements Serializable, DOFBlock {
    
    List<Residue> residues;//the residues moving freely, in order 
    //including the two end ones (only carboxyl or amine moving freely)
    
    int numFreeDOFs;
    int numFullDOFs;
    
    ArrayList<BBFreeDOF> freeDOFs = new ArrayList<>();
    
    double[][] freeDOFVoxel;//voxel of allowed values for free DOFs.  We'll minimize over this.
    //Indices: 0/1 (for min/max), DOF # 
    
    double[] curFreeDOFVals;//current values, so we can update one at a time

    double[] fullDOFCenter;//values of full DOFs at center of voxel
    DoubleMatrix1D freeDOFCenter;
    
    double[][] fullDOFPolys;//Polynomials represented in double[] format (see SeriesFitter)
    
    
    PepPlaneLinModel[] pepPlanes;//linear models of BB atoms in peptide planes
    //as a function of N and CA coordinates
    
    
    DoubleMatrix2D freeDOFMatrix;//matrix whose rows are the coefficients defining free DOFs
    
    static int polyOrder = 2;
    
    
    ArrayList<ArrayList<JacDerivEntry>> jacDerivs;
    
    boolean fixAnchors = true;
    boolean tooBigForSeries = false;//voxel too big for series; will probably be split
    
    ArrayList<Double> targetConstrVals;
    
    boolean useDistSqrt = true;//DEBUG!!!  use actual distances instead of quadratic constr
    
    public BBFreeBlock(){
        
    }
    
    public BBFreeBlock(List<Residue> residues, boolean fixAnchors){
        // Calculate the free DOFs, indicate appropriate voxels for them
        // and fit polynomials to express the full set of polynomials as a function of them
        init(residues,fixAnchors);
        setVoxelBySeries();//set voxel, adjust size if needed to ensure series validity
        //to specified residual in constraints at indicated fullDOF vals (maybe 0.01 distance err?)
    }
    
    
    
    public BBFreeBlock(List<Residue> residues, boolean fixAnchors, double[] voxLim, boolean tooBigForSeries){
        // Impose custom DOF limits
        init(residues,fixAnchors);
        this.tooBigForSeries = tooBigForSeries;

        freeDOFVoxel = new double[2][numFreeDOFs];
        for(int dof=0; dof<numFreeDOFs; dof++){
            freeDOFVoxel[0][dof] = -voxLim[dof];
            freeDOFVoxel[1][dof] = voxLim[dof];
        }
        //THESE ARE RELATIVE TO CENTER
        
        VoxelSeriesChecker vsc = new VoxelSeriesChecker(residues,numFreeDOFs,numFullDOFs,
            fullDOFPolys, pepPlanes, freeDOFMatrix, freeDOFCenter, fixAnchors, useDistSqrt);
        if(tooBigForSeries){
            System.out.println("Created voxel too big for series.");
        }
        else{
            double resid = vsc.getConstraintsResid(freeDOFVoxel);
            System.out.println("Manually created BBFreeBlock voxel.  Resid: "+resid);
        }
        
        targetConstrVals = vsc.targetConstraintVals;
    }
    
    
    
    
    private void init(List<Residue> residues, boolean fixAnchors){        
        this.fixAnchors = fixAnchors;
        this.residues = residues;
        int numRes = residues.size();
        numFreeDOFs = fixAnchors ? 2*numRes-6 : 2*numRes+2;
        numFullDOFs = fixAnchors ? 6*numRes-9 : 6*numRes-3;
        
        pepPlanes = new PepPlaneLinModel[numRes-1];
        for(int planeNum=0; planeNum<numRes-1; planeNum++)
            pepPlanes[planeNum] = new PepPlaneLinModel(residues.get(planeNum),residues.get(planeNum+1));
        
        //generating free DOFs based on current geometry 
        DoubleMatrix2D constrJac = getConstrJac();//Evaluate at current CA, N coords
        DoubleMatrix2D freeDOFCoeffs = getOrthogVectors(constrJac);
        
        
        //DEBUG!!!!!!
        //freeDOFCoeffs = selectLocalizableCoeffs(constrJac,freeDOFCoeffs);
        
        freeDOFMatrix = Algebra.DEFAULT.transpose(freeDOFCoeffs);
        
        for(int freeDOF=0; freeDOF<numFreeDOFs; freeDOF++){
            freeDOFs.add( new BBFreeDOF(freeDOFCoeffs.viewColumn(freeDOF), this, freeDOFs.size()) );
        }
        
        fullDOFCenter = new double[numFullDOFs];
        
        if(fixAnchors) {
            for(int resNum=1; resNum<numRes; resNum++){
                Residue curRes = residues.get(resNum);
                //set init coords
                System.arraycopy(curRes.getCoordsByAtomName("N"), 0, fullDOFCenter, 6*(resNum-1), 3);
                if(resNum<numRes-1)//CA free too
                    System.arraycopy(curRes.getCoordsByAtomName("CA"), 0, fullDOFCenter, 6*resNum-3, 3);
            }
        }
        else {
            for(int resNum=0; resNum<numRes; resNum++){
                Residue curRes = residues.get(resNum);
                if(resNum>0)
                    System.arraycopy(curRes.getCoordsByAtomName("N"), 0, fullDOFCenter, 6*resNum-3, 3);
                System.arraycopy(curRes.getCoordsByAtomName("CA"), 0, fullDOFCenter, 6*resNum, 3);
            }
        }

        freeDOFCenter = DoubleFactory1D.dense.make(numFreeDOFs);
        for(int f=0; f<numFreeDOFs; f++)
            freeDOFCenter.set( f, freeDOFs.get(f).evalAtFullDOFs(DoubleFactory1D.dense.make(fullDOFCenter)) );
        
        DoubleMatrix2D J = DoubleFactory2D.dense.compose(
                new DoubleMatrix2D[][] { new DoubleMatrix2D[] {Algebra.DEFAULT.transpose(freeDOFCoeffs)},
                    new DoubleMatrix2D[] {constrJac}
                } );//Jacobian of free DOFs and then constr vars
        makeTaylorSeries(J);
        
        curFreeDOFVals = new double[numFreeDOFs];//freeDOFCenter.copy().toArray();//may move away from center, so copy
        //treating these as relative!  for voxel and setDOFs purposes
    }
    
    private void setVoxelBySeries(){        
        freeDOFVoxel = new double[2][numFreeDOFs];
        Arrays.fill(freeDOFVoxel[0],-1);
        Arrays.fill(freeDOFVoxel[1],1);
        //THESE ARE RELATIVE TO CENTER
        
        double sizeStepFactor = 1.3;//how much we'll try decreasing it at a time, if not good enough
        
        double targetResid = 1e-4;//max residual allowed in bond lengths or dot products
        VoxelSeriesChecker vsc = new VoxelSeriesChecker(residues,numFreeDOFs,numFullDOFs,
            fullDOFPolys, pepPlanes, freeDOFMatrix, freeDOFCenter, fixAnchors, useDistSqrt);
        
        
        double resid = vsc.getConstraintsResid(freeDOFVoxel);
        System.out.println("Init voxel resid: "+resid);
        
        while (resid>targetResid) {
            
            for(int a=0; a<2; a++){
                for(int f=0; f<numFreeDOFs; f++)
                    freeDOFVoxel[a][f] /= sizeStepFactor;
            }

            resid = vsc.getConstraintsResid(freeDOFVoxel);
            
            System.out.println("Voxel reduced resid: "+resid);
        }
        
        targetConstrVals = vsc.targetConstraintVals;
    }
    
    
       
    //Cached stuff for derivative calcs
    //i, v, w will be full-DOF ("x") indices
    //u, j, k, l, m are free-DOF indices (we ultimately differentiate wrt these)
    //If > numFreeDOFs is used as a constr-DOF index
    double[][] deriv1;//Indices: i, j.  dx_i/dy_j
    double[][][] wderiv2;//Indices: i, j, w.  d^2 x_i / dy_j dx_w
    double[][][] deriv2;//d^x x_i / dy_j dy_k
    double[][][][] wderiv3;//d^3 x_i / dy_j dx_w dy_l
    double[][][][] deriv3;//dx^4 x_i / dy_j dy_k dy_l
    double[][][][][] wderiv4;//d^4 x_i / dy_j dx_w dy_l dy_m
    double[][][][][] deriv4;//d^4 x_i / dy_j dy_k dy_l dy_m
    
    
    private void allocateDerivs(int numFullDOFs, int numFreeDOFs){
        deriv1 = new double[numFullDOFs][numFullDOFs];
        wderiv2 = new double[numFullDOFs][numFullDOFs][numFullDOFs];
        deriv2 = new double[numFullDOFs][numFullDOFs][numFreeDOFs];
        
        if(polyOrder>=3){
            wderiv3 = new double[numFullDOFs][numFullDOFs][numFullDOFs][numFreeDOFs];
            deriv3 = new double[numFullDOFs][numFullDOFs][numFreeDOFs][numFreeDOFs];
            if(polyOrder>=4){
                wderiv4 = new double[numFullDOFs][numFullDOFs][numFullDOFs][numFreeDOFs][numFreeDOFs];
                deriv4 = new double[numFullDOFs][numFullDOFs][numFreeDOFs][numFreeDOFs][numFreeDOFs];
            }
        }
    }
    
    private void calcDerivs(int numFullDOFs, int numFreeDOFs, DoubleMatrix2D Jinv){
        //Calculate and cache derivatives for use in Taylor series
        
        allocateDerivs(numFullDOFs, numFreeDOFs);
        
        //1st order
        for(int i=0; i<numFullDOFs; i++){
            for(int j=0; j<numFullDOFs; j++)
                deriv1[i][j] = Jinv.get(i, j);
        }
        
        //2nd order
        for(int i=0; i<numFullDOFs; i++){
            for(int j=0; j<numFullDOFs; j++){
                for(int w=0; w<numFullDOFs; w++){
                    wderiv2[i][j][w] = 0;
                    for(JacDerivEntry jde : jacDerivs.get(w))
                        wderiv2[i][j][w] -= jde.val * deriv1[i][jde.u]  * deriv1[jde.v][j];
                }
                
                for(int k=0; k<numFreeDOFs; k++){
                    deriv2[i][j][k] = 0;
                    for(int w=0; w<numFullDOFs; w++)
                        deriv2[i][j][k] += wderiv2[i][j][w] * deriv1[w][k];
                }
            }
        }
        
        //3rd order
        if(polyOrder>=3){
            for(int i=0; i<numFullDOFs; i++){
                for(int j=0; j<numFullDOFs; j++){
                    
                    for(int l=0; l<numFreeDOFs; l++){
                    
                        for(int w=0; w<numFullDOFs; w++){
                            wderiv3[i][j][w][l] = 0;
                            for(JacDerivEntry jde : jacDerivs.get(w)){
                                wderiv3[i][j][w][l] -= jde.val * deriv2[i][jde.u][l]  * deriv1[jde.v][j];
                                wderiv3[i][j][w][l] -= jde.val * deriv1[i][jde.u]  * deriv2[jde.v][j][l];
                            }
                        }

                        for(int k=0; k<numFreeDOFs; k++){
                            deriv3[i][j][k][l] = 0;
                            for(int w=0; w<numFullDOFs; w++){
                                deriv3[i][j][k][l] += wderiv3[i][j][w][l] * deriv1[w][k];
                                deriv3[i][j][k][l] += wderiv2[i][j][w] * deriv2[w][k][l];
                            }
                        }
                    }
                }
            }
        }
        
        //4th order
        if(polyOrder>=4){
            for(int i=0; i<numFullDOFs; i++){
                for(int j=0; j<numFullDOFs; j++){
                    
                    for(int l=0; l<numFreeDOFs; l++){
                        for(int m=0; m<numFreeDOFs; m++){
                    
                            for(int w=0; w<numFullDOFs; w++){
                                wderiv4[i][j][w][l][m] = 0;
                                for(JacDerivEntry jde : jacDerivs.get(w)){
                                    wderiv4[i][j][w][l][m] -= jde.val * deriv3[i][jde.u][l][m]  * deriv1[jde.v][j];
                                    wderiv4[i][j][w][l][m] -= jde.val * deriv2[i][jde.u][l]  * deriv2[jde.v][j][m];
                                    wderiv4[i][j][w][l][m] -= jde.val * deriv2[i][jde.u][m]  * deriv2[jde.v][j][l];
                                    wderiv4[i][j][w][l][m] -= jde.val * deriv1[i][jde.u]  * deriv3[jde.v][j][l][m];
                                }
                            }

                            for(int k=0; k<numFreeDOFs; k++){
                                deriv4[i][j][k][l][m] = 0;
                                for(int w=0; w<numFullDOFs; w++){
                                    deriv4[i][j][k][l][m] += wderiv4[i][j][w][l][m] * deriv1[w][k];
                                    deriv4[i][j][k][l][m] += wderiv3[i][j][w][l] * deriv2[w][k][m];
                                    deriv4[i][j][k][l][m] += wderiv3[i][j][w][m] * deriv2[w][k][l];
                                    deriv4[i][j][k][l][m] += wderiv2[i][j][w] * deriv3[w][k][l][m];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    
    
    private void makeTaylorSeries(DoubleMatrix2D J){
        //Build fullDOFPolys as Taylor series, centered at orig conf
        //J is d(freeDOFs,constrVars)/d(fullDOFs)
        
        
        int numParams = SeriesFitter.getNumParams(numFreeDOFs, true, polyOrder);
        
        DoubleMatrix2D Jinv = Algebra.DEFAULT.inverse(J);//d(fullDOFs)/d(freeDOFs,constrVars)
        
        calcDerivs(numFullDOFs, numFreeDOFs, Jinv);
        
        fullDOFPolys = new double[numFullDOFs][];
        
        for(int fullDOF=0; fullDOF<numFullDOFs; fullDOF++){
            
            double[] series = new double[numParams];
            series[0] = fullDOFCenter[fullDOF];
            
            for(int freeDOF=0; freeDOF<numFreeDOFs; freeDOF++){
                series[freeDOF+1] = deriv1[fullDOF][freeDOF];
            }
            
            int coeffCount = numFreeDOFs+1;
            
            if(polyOrder>=2){

                for(int dof=0; dof<numFreeDOFs; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){

                        series[coeffCount] = deriv2[fullDOF][dof][dof2] / multFacProd(dof,dof2);
                        coeffCount++;
                    }
                }
            }
            
            if(polyOrder>=3){
                for(int dof=0; dof<numFreeDOFs; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                        for(int dof3=0; dof3<=dof2; dof3++){
                            
                            series[coeffCount] = deriv3[fullDOF][dof][dof2][dof3] / multFacProd(dof,dof2,dof3);
                            coeffCount++;
                        }
                    }
                }
            }
            
            if(polyOrder>=4){
                for(int dof=0; dof<numFreeDOFs; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                        for(int dof3=0; dof3<=dof2; dof3++){
                            for(int dof4=0; dof4<=dof3; dof4++){
                            
                                series[coeffCount] = deriv4[fullDOF][dof][dof2][dof3][dof4] / multFacProd(dof,dof2,dof3,dof4);
                                coeffCount++;
                            }
                        }
                    }
                }
            }
            
                    
            fullDOFPolys[fullDOF] = series;
        }
    }
    
    
    private static int multFacProd(int... a){
        //product of factorials of multiplicites of integers
        //argument will be in nonascending order
        //this is the denominator we need for Taylor series terms (enumerated w/o repeats)
        int curFactor = 1;
        int ans = 1;
        
        for(int b=0; b<a.length-1; b++){
            if(a[b]==a[b+1])
                curFactor++;
            else
                curFactor = 1;
            
            ans *= curFactor;
        }
        
        return ans;
    }

    @Override
    public DOFBlock copyForNewMolecule(Molecule mol, LinkedHashMap<DegreeOfFreedom, DegreeOfFreedom> copiedDOFMap) {
        BBFreeBlock copiedBlock = new BBFreeBlock();
        
        copiedBlock.residues = new ArrayList<>();
        for(Residue res : residues){
            Residue newMolRes = mol.getResByPDBResNumber(res.getPDBResNumber());
            copiedBlock.residues.add(newMolRes);
        }
        
        copiedBlock.freeDOFs = new ArrayList<>();
        for(BBFreeDOF dof : freeDOFs){
            BBFreeDOF copiedDOF = new BBFreeDOF(dof.coeffs, copiedBlock, dof.indexInBlock);
            copiedDOFMap.put(dof, copiedDOF);
            copiedBlock.freeDOFs.add(copiedDOF);
        }
        
        //shallow copy is OK for read-only fields
        copiedBlock.freeDOFVoxel = freeDOFVoxel;
        copiedBlock.curFreeDOFVals = curFreeDOFVals.clone();
        copiedBlock.fullDOFCenter = fullDOFCenter;
        copiedBlock.freeDOFCenter = freeDOFCenter;
        copiedBlock.fullDOFPolys = fullDOFPolys;
        copiedBlock.pepPlanes = pepPlanes;
        copiedBlock.freeDOFMatrix = freeDOFMatrix;
        copiedBlock.jacDerivs = jacDerivs;
        
        copiedBlock.targetConstrVals = targetConstrVals;
        copiedBlock.tooBigForSeries = tooBigForSeries;
        copiedBlock.fixAnchors = fixAnchors;
        copiedBlock.numFreeDOFs = numFreeDOFs;
        copiedBlock.numFullDOFs = numFullDOFs;
        //no need to copy all the deriv1, deriv2, etc. since they're just used to construct the block
        
        return copiedBlock;
    }

    @Override
    public List<Residue> listResidues() {
        return residues;
    }
    
    
    
    private class JacDerivEntry implements Serializable {
        //The derivative of the Jacobian with respect to the full DOFs is sparse
        //Each of these entries, when placed in jacDerivs[w],
        //indicates that d^2 (free DOF u)/ d(full DOF v) d(full DOF w) = val
        int u, v;
        double val;

        public JacDerivEntry(int u, int v, double val) {
            this.u = u;
            this.v = v;
            this.val = val;
        }
    }
    
    
    private static double invMatrixDeriv(DoubleMatrix2D Minv, int i, int j, int u, int v){
        //d(M^-1)_ij / dM_uv
        return - Minv.get(i, u) * Minv.get(v, j);
    }
    
    public static DoubleMatrix2D getOrthogVectors(DoubleMatrix2D M){
        //Get a matrix whose columns are orthogonal to the row space of M
        //which expected to be nonsingular
        DoubleMatrix2D Maug = DoubleFactory2D.dense.make(M.columns(), M.columns());
        Maug.viewPart(0, 0, M.rows(), M.columns()).assign(M);
        
        SingularValueDecomposition svd = new SingularValueDecomposition(Maug);
        
        int numOrthVecs = M.columns() - M.rows();
        if(svd.rank() != M.rows()){
            throw new RuntimeException("ERROR: Singularity in constr jac.  Rank: "+svd.rank());
        }
        
        DoubleMatrix2D orthVecs = svd.getV().viewPart(0, M.rows(), M.columns(), numOrthVecs);
        
        //DEBUG!!!  Should be 0 and identity respecitvely
        /*DoubleMatrix2D check = Algebra.DEFAULT.mult(M, orthVecs);
        DoubleMatrix2D checkOrth = Algebra.DEFAULT.mult( Algebra.DEFAULT.transpose(orthVecs), orthVecs);
        */
        
        
        return orthVecs;
    }
    
    
    public DoubleMatrix2D getConstrJac(){
        //Jacobian constraints, evaluated at current (original) coordinates
        
        int numRes = residues.size();
        
        jacDerivs = new ArrayList<ArrayList<JacDerivEntry>>();
        for(int f=0; f<numFullDOFs; f++)
            jacDerivs.add(new ArrayList<JacDerivEntry>());
        
        DoubleMatrix2D ans = DoubleFactory2D.dense.make(numFullDOFs-numFreeDOFs,numFullDOFs);
        
        
        int constrCount = 0;
        
        for(int resNum=0; resNum<numRes; resNum++){
            
            Residue curRes = residues.get(resNum);
            
            double NCoord[] = curRes.getCoordsByAtomName("N");
            double CACoord[] = curRes.getCoordsByAtomName("CA");
            int NFreeIndex = fixAnchors ? 2*(resNum-1) : 2*resNum-1;//index of N among free atoms

            if(resNum>0){//these constrs act only on fixed coords for resNum=0
                
                Residue prevRes = residues.get(resNum-1);
                double prevCACoord[] = prevRes.getCoordsByAtomName("CA");
                
                //N to CA constr for resNum
                makeDistConstrJac(ans,constrCount,NCoord,CACoord,NFreeIndex,NFreeIndex+1);
                constrCount++;

                //CA to CA constr
                makeDistConstrJac(ans,constrCount,prevCACoord,CACoord,NFreeIndex-1,NFreeIndex+1);
                constrCount++;

                //last CA to N constr
                makeDistConstrJac(ans,constrCount,prevCACoord,NCoord,NFreeIndex-1,NFreeIndex);
                constrCount++;
            }
            
            //OK and finally N-CA-C' angle constr
            if(resNum==numRes-1){
                if(fixAnchors){//C' fixed: special form of constraint
                    double CCoord[] = curRes.getCoordsByAtomName("C");
                    makeLastNCACCConstrJac( ans, constrCount, CACoord, CCoord, NFreeIndex );
                    constrCount++;
                }
            }
            else if ( fixAnchors || resNum>0 ){
                Residue nextRes = residues.get(resNum+1);
                double NCoordNext[] = nextRes.getCoordsByAtomName("N");
                double CACoordNext[] = nextRes.getCoordsByAtomName("CA");

                makeNCACCConstrJac(ans, constrCount, NCoord, CACoord, 
                    NCoordNext, CACoordNext, NFreeIndex, resNum);
                constrCount++;
            }
        }
                
        return ans;
    }
    
    
    private void recordJacDeriv(int constrNum, int fullDOF1, int fullDOF2, double val){
        jacDerivs.get(fullDOF1).add( new JacDerivEntry(constrNum+numFreeDOFs,fullDOF2,val) );
    }
    
    
    private void makeDistConstrJac( DoubleMatrix2D constrJac, int constrNum,
            double coord1[], double coord2[],
            int atomNum1, int atomNum2 ) {
        //Enter the Jacobian of a distance constraint between the atoms as #constrNum
        //in constrJac.  Atom numbers (among free N, CA) and (initial) coordinates given
        //-1 for atomNum means fixed atom, so not part of gradient
        
        if(useDistSqrt){
            double dist = VectorAlgebra.distance(coord1, coord2);
            for(int dim=0; dim<3; dim++){
                if(atomNum1>=0 && 3*atomNum1<constrJac.columns()){
                    double diff = coord1[dim]-coord2[dim];
                    constrJac.set(constrNum, 3*atomNum1+dim, diff/dist );
                    
                    for(int dim2=0; dim2<3; dim2++){
                        double secondDeriv = -diff*(coord1[dim2]-coord2[dim2]) / (dist*dist*dist);
                        if(dim==dim2)
                            secondDeriv += 1/dist;
                        recordJacDeriv(constrNum, 3*atomNum1+dim, 3*atomNum1+dim2, secondDeriv);
                        if(atomNum2>=0 && 3*atomNum2<constrJac.columns())//atomNum2 free to move
                            recordJacDeriv(constrNum, 3*atomNum1+dim, 3*atomNum2+dim2, -secondDeriv);
                    }
                }

                if(atomNum2>=0 && 3*atomNum2<constrJac.columns()){
                    double diff = coord2[dim]-coord1[dim];
                    constrJac.set(constrNum, 3*atomNum2+dim, diff/dist );
                    
                    for(int dim2=0; dim2<3; dim2++){
                        double secondDeriv = -diff*(coord2[dim2]-coord1[dim2]) / (dist*dist*dist);
                        if(dim==dim2)
                            secondDeriv += 1/dist;
                        recordJacDeriv(constrNum, 3*atomNum2+dim, 3*atomNum2+dim2, secondDeriv);
                        if(atomNum1>=0 && 3*atomNum1<constrJac.columns())//atomNum1 free to move
                            recordJacDeriv(constrNum, 3*atomNum2+dim, 3*atomNum1+dim2, -secondDeriv);
                    }
                }
            }
        }
        else {
            for(int dim=0; dim<3; dim++){
                if(atomNum1>=0 && 3*atomNum1<constrJac.columns()){
                    constrJac.set(constrNum, 3*atomNum1+dim, 2*coord1[dim]-2*coord2[dim]);

                    recordJacDeriv(constrNum, 3*atomNum1+dim, 3*atomNum1+dim, 2);
                    if(atomNum2>=0 && 3*atomNum2<constrJac.columns())//atomNum2 free to move
                        recordJacDeriv(constrNum, 3*atomNum1+dim, 3*atomNum2+dim, -2);
                }

                if(atomNum2>=0 && 3*atomNum2<constrJac.columns()){
                    constrJac.set(constrNum, 3*atomNum2+dim, 2*coord2[dim]-2*coord1[dim]);

                    recordJacDeriv(constrNum, 3*atomNum2+dim, 3*atomNum2+dim, 2);
                    if(atomNum1>=0 && 3*atomNum1<constrJac.columns())//atomNum1 free to move
                        recordJacDeriv(constrNum, 3*atomNum2+dim, 3*atomNum1+dim, -2);
                }
            }
        }
    }
    
    
    private void makeNCACCConstrJac( DoubleMatrix2D constrJac, int constrNum,
            double NCoord[], double CACoord[], double NCoordNext[], double CACoordNext[],
            int atomNumN, int pepPlaneNum ){ 
        //constraint the dot product (C'-CA) dot (N-CA)
        //Uses N, CA coord variables for this residue and the next
        //(since C' coords are a linear function of this CA and next N, CA coords)
        //(For this purpose we constrain not the actual C' but its projection into the 
        //CA-N-C plane)
        
        if(useDistSqrt && atomNumN>=0){
            //non-sqrt constr is linear for the first res, so just use that
            makeNCACConstrJacSqrt(constrJac,constrNum,NCoord,CACoord,NCoordNext,
                    CACoordNext,atomNumN,pepPlaneNum);
            return;
        }
        
        //we'll need expansion coefficients for C'...
        double[] expansionCoeffs = pepPlanes[pepPlaneNum].getProjCAtomCoeffs();
        double a = expansionCoeffs[0];
        double b = expansionCoeffs[1];
        double c = expansionCoeffs[2];
        
        for(int dim=0; dim<3; dim++){
            if(atomNumN>=0){//not the first residue, which has fixed N & CA
                double NGrad = (a-1)*CACoord[dim] + b*NCoordNext[dim] + c*CACoordNext[dim];
                double CAGrad  = (a-1)*NCoord[dim] + 2*(1-a)*CACoord[dim] - b*NCoordNext[dim] - c*CACoordNext[dim];
                constrJac.set(constrNum, 3*atomNumN+dim, NGrad);
                constrJac.set(constrNum, 3*atomNumN+3+dim, CAGrad);
                
                //jac derivs for N
                recordJacDeriv(constrNum, 3*atomNumN+dim, 3*atomNumN+3+dim, a-1);
                recordJacDeriv(constrNum, 3*atomNumN+dim, 3*atomNumN+6+dim, b);
                if(3*atomNumN+9<constrJac.columns())
                    recordJacDeriv(constrNum, 3*atomNumN+dim, 3*atomNumN+9+dim, c);
                
                //and for C
                recordJacDeriv(constrNum, 3*atomNumN+3+dim, 3*atomNumN+dim, a-1);
                recordJacDeriv(constrNum, 3*atomNumN+3+dim, 3*atomNumN+3+dim, 2*(1-a));
                recordJacDeriv(constrNum, 3*atomNumN+3+dim, 3*atomNumN+6+dim, -b);
                if(3*atomNumN+9<constrJac.columns())
                    recordJacDeriv(constrNum, 3*atomNumN+3+dim, 3*atomNumN+9+dim, -c);
            }
            
            constrJac.set(constrNum, 3*atomNumN+6+dim, b*NCoord[dim]-b*CACoord[dim]);
            
            if(atomNumN>=0){
                recordJacDeriv(constrNum, 3*atomNumN+6+dim, 3*atomNumN+dim, b);
                recordJacDeriv(constrNum, 3*atomNumN+6+dim, 3*atomNumN+3+dim, -b);
            }
            
            
            if(3*atomNumN+9<constrJac.columns()){//i.e., not second-to-last res (which has fixed CANext)
                constrJac.set(constrNum, 3*atomNumN+9+dim, c*NCoord[dim]-c*CACoord[dim]);
                
                if(atomNumN>=0){
                    recordJacDeriv(constrNum, 3*atomNumN+9+dim, 3*atomNumN+dim, c);
                    recordJacDeriv(constrNum, 3*atomNumN+9+dim, 3*atomNumN+3+dim, -c);
                }
            }
        }
    }
    
    
    private void makeNCACConstrJacSqrt( DoubleMatrix2D constrJac, int constrNum,
            double NCoord[], double CACoord[], double NCoordNext[], double CACoordNext[],
            int atomNumN, int pepPlaneNum ){
        //This version just constrains the N-C' distance
         //given fixed geometry of each peptide plane,
         //this is equivalent to the other NCAC constr
         //this function assumes all four atoms except maybe last CA are free to move
         //(otherwise dot product constraint becomes linear and we can just use that)
        
        //we'll need expansion coefficients for C'...
        double[] expansionCoeffs = pepPlanes[pepPlaneNum].getProjCAtomCoeffs();
        double a = expansionCoeffs[0];
        double b = expansionCoeffs[1];
        double c = expansionCoeffs[2];
        
        //calculated project C' coordinates
        double CCoord[] = new double[3];
        for(int dim=0; dim<3; dim++)
            CCoord[dim] = a*CACoord[dim] + b*NCoordNext[dim] + c*CACoordNext[dim];
        
        
        double dist = VectorAlgebra.distance(NCoord, CCoord);
        for(int dim=0; dim<3; dim++){
            
            double diff = NCoord[dim] - CCoord[dim];
            constrJac.set(constrNum, 3*atomNumN+dim, diff/dist );

            for(int dim2=0; dim2<3; dim2++){
                double secondDeriv = -diff*(NCoord[dim2]-CCoord[dim2]) / (dist*dist*dist);
                if(dim==dim2)
                    secondDeriv += 1/dist;
                recordJacDeriv(constrNum, 3*atomNumN+dim, 3*atomNumN+dim2, secondDeriv);

                //second derivatives wrt peptide plane atoms (which determine C' coord)
                //d(dr/dx)/dy = d(dr/dx)/dC' dC'/dy
                recordJacDeriv(constrNum, 3*atomNumN+dim, 3*atomNumN+3+dim2, -a*secondDeriv);
                recordJacDeriv(constrNum, 3*atomNumN+dim, 3*atomNumN+6+dim2, -b*secondDeriv);
                if(3*atomNumN+9<constrJac.columns())//atomNum2 free to move
                    recordJacDeriv(constrNum, 3*atomNumN+dim, 3*atomNumN+9+dim2, -c*secondDeriv);
            }
            
            
            for(int pepAtomNum=0; pepAtomNum<3; pepAtomNum++){//which atom of the peptide plane containing C'
                int curAtomNum = atomNumN + pepAtomNum + 1;
                if(3*curAtomNum<constrJac.columns()){
                    diff = CCoord[dim]-NCoord[dim];
                    constrJac.set(constrNum, 3*curAtomNum+dim, expansionCoeffs[pepAtomNum]*diff/dist );

                    for(int dim2=0; dim2<3; dim2++){
                        double secondDeriv = -diff*(CCoord[dim2]-NCoord[dim2]) / (dist*dist*dist);//this is d^2r/dC'_dim dC'_dim2
                        if(dim==dim2)
                            secondDeriv += 1/dist;
                        
                        //for peptide plane atoms A and A2, d^2 r/dA_dim dA2_dim2 = d^2 r/dC'_dim dC'_dim2 dC'_dim dA_dim/dC'_dim dA2_dim2/dC'_dim2
                        for(int pepAt2=0; pepAt2<3; pepAt2++){
                            int atNum2=atomNumN+pepAt2+1;
                            if(3*atNum2<constrJac.columns()){
                                recordJacDeriv(constrNum, 3*curAtomNum+dim, 3*atNum2+dim2, expansionCoeffs[pepAtomNum]*expansionCoeffs[pepAt2]*secondDeriv);
                            }
                        }
                            
                        //also do second deriv with N coordinate
                        recordJacDeriv(constrNum, 3*curAtomNum+dim, 3*atomNumN+dim2, -expansionCoeffs[pepAtomNum]*secondDeriv);
                    }
                }
            }
        }
    }
    
    
    private void makeLastNCACCConstrJac( DoubleMatrix2D constrJac, int constrNum,
            double CACoord[], double CCoord[], int atomNumN ){ 
        //Last constraint is simpler, because CA and C' are both fixed
        //Do this even with sqrt because linear constraint can be handled exactly!
        
        for(int dim=0; dim<3; dim++)
            constrJac.set(constrNum, 3*atomNumN+dim, CCoord[dim] - CACoord[dim]);
        
        //no jac derivs because this constraint is linear
    }
    
    
    void applyNCACoords(double[] fullDOFVals){
        //place the N and CA atoms as indicated by the full DOF vals
        int numRes = residues.size();
        for(int resNum=0; resNum<numRes; resNum++){
            
            Residue curRes = residues.get(resNum);
            int NFreeIndex = fixAnchors ? 2*(resNum-1) : 2*resNum-1;//index of N among free atoms

            if(resNum>0 || !fixAnchors){//N or CA moves
                int CAIndex = curRes.getAtomIndexByName("CA");
                int NIndex = curRes.getAtomIndexByName("N");

                for(int dim=0; dim<3; dim++){
                    if( resNum<numRes-1 || !fixAnchors )//CA moves
                        curRes.coords[3*CAIndex+dim] = fullDOFVals[3*NFreeIndex+3+dim];
                    
                    if(resNum>0)//N moves (non-fixed anchor case)
                        curRes.coords[3*NIndex+dim] = fullDOFVals[3*NFreeIndex+dim];
                }
            }
        }
    }
    
    
    double[] recordNCACoords(){
        //record the full DOF vals indicated by the current conformation
        int numRes = residues.size();
        double[] fullDOFVals = new double[numFullDOFs];
        for(int resNum=0; resNum<numRes; resNum++){
            
            Residue curRes = residues.get(resNum);
            int NFreeIndex = fixAnchors ? 2*(resNum-1) : 2*resNum-1;//index of N among free atoms

            if(resNum>0 || !fixAnchors){//N or CA moves
                int CAIndex = curRes.getAtomIndexByName("CA");
                int NIndex = curRes.getAtomIndexByName("N");

                for(int dim=0; dim<3; dim++){
                    if( resNum<numRes-1 || !fixAnchors )//CA moves
                        fullDOFVals[3*NFreeIndex+3+dim] = curRes.coords[3*CAIndex+dim];
                    
                    if(resNum>0)//N moves (non-fixed anchor case)
                        fullDOFVals[3*NFreeIndex+dim] = curRes.coords[3*NIndex+dim];
                }
            }
        }
        return fullDOFVals;
    }
    
    
    
    public boolean setDOFs(DoubleMatrix1D x){
        //x: free DOFs (relative to center, so can eval polys directly)
        
        int numRes = residues.size();
        
        double fullDOFVals[] = computeFullDOFValues(x);
        
        if(fullDOFVals==null)//DEBUG!!! failed
            return false;
                    
        //record current information needed for placement of sidechain
        double genChi1[] = new double[numRes];
        //double genChi1SinCos[][] = new double[numRes][];//thought it would help, doesn't
        //turns out acos looks really expensive with sampled profiling, not with instrumented
        
        double SCTranslations[][] = new double[numRes][3];
        for(int resNum=0; resNum<numRes; resNum++){
            
            Residue curRes = residues.get(resNum);
            int NFreeIndex = fixAnchors ? 2*(resNum-1) : 2*resNum-1;//index of N among free atoms
            
            //the sidechain will be translated based on CA motion
            if( (resNum>0 && resNum<numRes-1) || !fixAnchors ){//CA moves, so translate sidechain with it
                double curCACoords[] = curRes.getCoordsByAtomName("CA");
                for(int dim=0; dim<3; dim++)
                    SCTranslations[resNum][dim] = fullDOFVals[3*NFreeIndex+3+dim] - curCACoords[dim];
            }
            
            genChi1[resNum] = GenChi1Calc.getGenChi1(curRes);
            //genChi1SinCos[resNum] = GenChi1Calc.getGenChi1SinCos(curRes);
        }
        
        //OK now place the backbone
        
        applyNCACoords(fullDOFVals);
        
        //OK now that the CA's and N's are in place, we can finish each peptide plane
        for(int pepPlaneNum=0; pepPlaneNum<numRes-1; pepPlaneNum++){
            //set each atom based on CA and N (linear relationship)
            Residue res1 = residues.get(pepPlaneNum);//residue at start of peptide plane...
            Residue res2 = residues.get(pepPlaneNum+1);//...and at end
            
            double CA1[] = res1.getCoordsByAtomName("CA");
            double N[] = res2.getCoordsByAtomName("N");
            double CA2[] = res2.getCoordsByAtomName("CA");
                        
            res1.setCoordsByAtomName("C", pepPlanes[pepPlaneNum].calcCCoords(CA1,N,CA2,false));
            res1.setCoordsByAtomName("O", pepPlanes[pepPlaneNum].calcOCoords(CA1,N,CA2));
            res2.setCoordsByAtomName("H", pepPlanes[pepPlaneNum].calcHCoords(CA1,N,CA2));
        }
            
        
        
        //OK and now that the backbone atoms are in place, we can handle the sidechains and HA's
        for(int resNum=0; resNum<numRes; resNum++){
            //first translate into place...
            RigidBodyMotion motion = new RigidBodyMotion(new double[3], RotationMatrix.identity(), SCTranslations[resNum]);
            SidechainIdealizer.moveSidechain(residues.get(resNum), motion);
            //...now idealize
            SidechainIdealizer.idealizeSidechain(residues.get(resNum));
            //and get gen chi1 to where it was before, so this BB motion commutes w/ sidechain dihedral changes
            GenChi1Calc.setGenChi1(residues.get(resNum), genChi1[resNum]);
            //GenChi1Calc.setGenChi1SinCos(residues.get(resNum), genChi1SinCos[resNum]);
        }
        
        curFreeDOFVals = x.toArray();
        return true;
    }
    
    
    public ArrayList<BBFreeDOF> getDOFs(){
        return freeDOFs;
    }

    public double[][] getFreeDOFVoxel() {
        return freeDOFVoxel;
    }

    public List<Residue> getResidues() {
        return residues;
    }
    
    
    static DoubleMatrix1D toDM1D(ArrayList<Double> b){
        double[] arr = b.stream().mapToDouble(Double::doubleValue).toArray();
        return DoubleFactory1D.dense.make(arr);
    }
    
    static DoubleMatrix1D concat(DoubleMatrix1D a, ArrayList<Double> b){
        DoubleMatrix1D[] parts = new DoubleMatrix1D[] {a, toDM1D(b)};
        return DoubleFactory1D.dense.make(parts);
    }
    
    
    //DEBUG!!!!
    ArrayList<Double> getConstrVals(double fullDOFVals[]){
        double[] backupFullDOFVals = recordNCACoords();
        applyNCACoords(fullDOFVals);
        //calculate constraints & their jac
        ArrayList<Double> constrVals = new VoxelSeriesChecker(residues,numFreeDOFs,numFullDOFs,
            fullDOFPolys, pepPlanes, freeDOFMatrix, freeDOFCenter, fixAnchors, useDistSqrt).targetConstraintVals;
        applyNCACoords(backupFullDOFVals);
        return constrVals;
    }
    
    DoubleMatrix2D getConstrJac(double fullDOFVals[]){
        //this version evaluates the Jacobian at fullDOFVals rather than current coords
        double[] backupFullDOFVals = recordNCACoords();
        applyNCACoords(fullDOFVals);
        //calculate constraints & their jac
        DoubleMatrix2D constrJac = getConstrJac();
        applyNCACoords(backupFullDOFVals);
        return constrJac;
    }
    
    ArrayList<Double> getConstrVals(double y[], DoubleMatrix1D changeDir, double coeff){
        //get ConstrVals for y + coeff*changeDir
        double fullDOFVals[] = y.clone();
        for(int a=0; a<y.length; a++)
            fullDOFVals[a] += coeff*changeDir.get(a);
        return getConstrVals(fullDOFVals);
    }
    
    
    double getConstrValsResid(double y[], DoubleMatrix1D changeDir, double coeff){
        ArrayList<Double> constrVals = getConstrVals(y,changeDir,coeff);
        double resid = 0;
        for(int d=0; d<constrVals.size(); d++){
            double dev = constrVals.get(d) - targetConstrVals.get(d);
            resid += dev*dev/numFullDOFs;//to match normalization used in computeFullDOFValues
        }
        return resid;
    }
    
    
    void numCheckConstrJac(DoubleMatrix2D constrJac, double fullDOFVals[]){
        DoubleMatrix2D jacCheck = constrJac.like();
        double step = 1e-6;
        for(int fullDOF=0; fullDOF<numFullDOFs; fullDOF++){
            DoubleMatrix1D changeDir = DoubleFactory1D.dense.make(numFullDOFs);
            changeDir.set(fullDOF, step);
            DoubleMatrix1D constrValUp = toDM1D(getConstrVals(fullDOFVals,changeDir,1));
            DoubleMatrix1D constrValDown = toDM1D(getConstrVals(fullDOFVals,changeDir,-1));
            constrValUp.assign(constrValDown,Functions.minus).assign(Functions.mult(0.5/step));
            jacCheck.viewColumn(fullDOF).assign(constrValUp);
        }
        double diff = Algebra.DEFAULT.normF( jacCheck.copy().assign(constrJac,Functions.minus) );
        int aaa = 0;//DEBUG!!
    }
    
    
    static double meanSquareDiff(ArrayList<Double> a, ArrayList<Double> b){
        double ans = 0;
        double n = a.size();
        if(n!=b.size())
            throw new RuntimeException("ERROR: List size mismatch");
        for(int j=0; j<n; j++){
            double dev = a.get(j)-b.get(j);
            ans += dev*dev;
        }
        return ans / n;
    }
    
    ///DEBUG!!!!  this is to avoid singularities
    boolean useSingularityDetector = true;
    
    private double[] computeFullDOFValues(DoubleMatrix1D freeDOFVals){
        //Compute full DOF values (N, CA coords) as a function of the free DOFs
        int numDOFsFull = fullDOFPolys.length;
        double[] fullDOFVals = new double[numDOFsFull];
        
        if(tooBigForSeries){
            //since we're not sure about series quality, we will check constraint satisfaction
            //and if needed, break up the journey from the origin (where series is initially centered) to freeDOFVals into steps
            //at each of which the series is recalculated
            
            
            //back up the initial polynomial so we can use it again later
            double[] backupFullDOFCenter = fullDOFCenter;
            double[][] backupFullDOFPolys = fullDOFPolys;
            
            
            double targetResid = 1e-4;
            //target residual in constr & free DOF values
            //between what we want and what our calculated full DOF values imply
            
            //for intermediate values x=((1-a)*freeDOFCenter + a*freeDOFVals), we demand resid<a*targetResid
            //progress is the current value of a
            double progress = 0;
            //DoubleMatrix1D x = DoubleFactory1D.dense.make(numFreeDOFs);//where we currently are in free DOF space (and where current series is centered)
            //since initial series is centered at freeDOFCenter, we effectively start there
            //we'll takes steps step*progressDir
            double step = 1;
            
            int numStepsToGo = 1;//number of steps to reach freeDOFVals from freeDOFCenter
            
            while(numStepsToGo>0){
                
                //attempt to make the step
               // DoubleMatrix1D curTarget = freeDOFVals.copy().assign(Functions.mult(progress+step));
                DoubleMatrix1D stepVector = freeDOFVals.copy().assign(Functions.mult(step));//(current target)-(current series center)
                for(int fullDOF=0; fullDOF<numDOFsFull; fullDOF++)
                    fullDOFVals[fullDOF] = SeriesFitter.evalSeries(fullDOFPolys[fullDOF], stepVector, numFreeDOFs, true, polyOrder);

                double resid = meanSquareDiff(getConstrVals(fullDOFVals),targetConstrVals);
                if( resid > (progress+step)*targetResid ){
                    //not good enough...need smaller steps
                    step /= 2;
                    numStepsToGo *= 2;
                    
                    if(numStepsToGo>10 && useSingularityDetector){//DEBUG!!! singularity likely
                        fullDOFCenter = backupFullDOFCenter;
                        fullDOFPolys = backupFullDOFPolys;
                        return null;
                    }
                }
                else {
                    //step was successful
                    progress += step;
                    //x = curTarget;
                    numStepsToGo--;
                    
                    //if not yet at freeDOFVals then prepare series centered at x to get us closer
                    if(numStepsToGo>0){
                        fullDOFCenter = fullDOFVals.clone();
                        DoubleMatrix2D J = DoubleFactory2D.dense.compose(
                            new DoubleMatrix2D[][] { new DoubleMatrix2D[] {freeDOFMatrix},//doesn't change, neither does freeDOFCenter (together they define free DOFs)
                                new DoubleMatrix2D[] {getConstrJac(fullDOFVals)}
                            } );
                        makeTaylorSeries(J);//this will create fullDOFPolys
                        
                        
                        //DEBUG!!!!
                        //System.out.println("DET J: "+Algebra.DEFAULT.det(J));
                        //detects if we cross a singularity
                        //in general, if we do, then step size tanks (reasonable, and we don't need to model such cases).  
                    }
                }
            }
            
            //revert to original polynomial for use next time
            //this gets us to a starting point wiht exactly right constraints
            //which is important for long-term stability
            fullDOFCenter = backupFullDOFCenter;
            fullDOFPolys = backupFullDOFPolys;
        }
        else {
            for(int fullDOF=0; fullDOF<numDOFsFull; fullDOF++)
                fullDOFVals[fullDOF] = SeriesFitter.evalSeries(fullDOFPolys[fullDOF], freeDOFVals, freeDOFVals.size(), true, polyOrder);
        }
        
        return fullDOFVals;
    }
    
    
    
    private double[] computeFullDOFValuesOld(DoubleMatrix1D freeDOFVals){
        //Compute full DOF values (N, CA coords) as a function of the free DOFs
        int numDOFsFull = fullDOFPolys.length;
        double[] fullDOFVals = new double[numDOFsFull];
        
        for(int fullDOF=0; fullDOF<numDOFsFull; fullDOF++)
            fullDOFVals[fullDOF] = SeriesFitter.evalSeries(fullDOFPolys[fullDOF], freeDOFVals, freeDOFVals.size(), true, polyOrder);
        
        if(tooBigForSeries){
            //Let's refine with a linear expansion
            //centered at the new x vals
            
            double targetResid = 1e-4;
            //target residual in constr & free DOF values
            //between what we want and what our calculated full DOF values imply
            
            while(true){
                //OK to evaluate constraints and their derivatives
                //we will temporarily set the N's and CA's, set them back later
                double[] backupFullDOFVals = recordNCACoords();
                applyNCACoords(fullDOFVals);
                //calculate constraints & their jac
                ArrayList<Double> constrVals = new VoxelSeriesChecker(residues,numFreeDOFs,numFullDOFs,
                    fullDOFPolys, pepPlanes, freeDOFMatrix, freeDOFCenter, fixAnchors, useDistSqrt).targetConstraintVals;
                DoubleMatrix2D constrJac = getConstrJac();      
                applyNCACoords(backupFullDOFVals);//set coords back

                numCheckConstrJac(constrJac,fullDOFVals);//DEBUG!!
                

                //calculate all the transformed vals (free and constr).  Call this x

                DoubleMatrix1D xfree = Algebra.DEFAULT.mult(freeDOFMatrix, 
                        DoubleFactory1D.dense.make(fullDOFVals));
                xfree.assign(freeDOFCenter, Functions.minus);
                DoubleMatrix1D x = concat(xfree,constrVals);

            //DEBUG!!!!!!
            //freeDOFCoeffs = selectLocalizableCoeffs(constrJac,freeDOFCoeffs);

                DoubleMatrix2D J = DoubleFactory2D.dense.compose(
                    new DoubleMatrix2D[][] { new DoubleMatrix2D[] {freeDOFMatrix},
                        new DoubleMatrix2D[] {constrJac}
                    } );


                //current full DOF vals imply a set of free coords and constraint vals
                //we can center our series there

                DoubleMatrix1D desiredTransformedVals = concat(freeDOFVals, targetConstrVals);
                //dx = desired remaining change in transformed coords
                DoubleMatrix1D dx = desiredTransformedVals.copy().assign(x, Functions.minus);
                DoubleMatrix2D Jinv = Algebra.DEFAULT.inverse(J);

                
                
                DoubleMatrix1D dy = Algebra.DEFAULT.mult(Jinv, dx);//update to full DOF vals
                
                //DEBUG!!!
                double resid0 = getConstrValsResid(fullDOFVals, dy, 0);//should match dx
                double resid1 = getConstrValsResid(fullDOFVals, dy, 1);
                double residp1 = getConstrValsResid(fullDOFVals, dy, 0.1);
                double residp01 = getConstrValsResid(fullDOFVals, dy, 0.01);
                double residp001 = getConstrValsResid(fullDOFVals, dy, 0.001);

                
                for(int fullDOF=0; fullDOF<numDOFsFull; fullDOF++)
                    fullDOFVals[fullDOF] += dy.get(fullDOF);
                
                double resid = Algebra.DEFAULT.norm2(dx) / numFullDOFs;
                if(resid <= targetResid)
                    break;
            }
        }
        
        return fullDOFVals;
    }
    
    
    DoubleMatrix2D resPoseBasis(int resNum){
        //basis for the pose of the residue (as column vectors in full-DOF space)
        //basis can be redundant.
        //CAN THIS ALSO BE DONE AT OTHER CENTERS TO GET APPROPRIATE BOUNDS?
        //RESULTING DOFS WILL AT LEAST BE ORTH TO CONSTRAINTS WHEN RESTRICTED TO OUR RES' BASIS
        //REGARDLESS OF POSE AT OTHER RES
        
        int atomNumN = fixAnchors ? 2*(resNum-1) : 2*resNum-1;
        DoubleMatrix2D ans;
        
        if(resNum>0 && resNum<residues.size()-1){//N, CA, C' all moveable
            ans = DoubleFactory2D.dense.make(numFullDOFs,9);
            double CCoeffs[] = pepPlanes[resNum].getProjCAtomCoeffs();
            
            for(int dim=0; dim<3; dim++){
                ans.set(3*atomNumN+dim, dim, 1);//N
                ans.set(3*atomNumN+3+dim, dim+3, 1);//CA
                for(int at=0; at<3; at++){//C'
                    int pAtomOffset = 3*atomNumN+3*(at+1);
                    if(pAtomOffset < numFullDOFs)
                        ans.set(pAtomOffset+dim, dim+6, CCoeffs[at]);
                }
            }
        }
        else if(resNum==0){//first residue, no N and maybe no CA
            //we do C' and O.  This scheme ensures enough DOFs to describe CB and thus sidechain motion
            int basisSize = fixAnchors ? 6 : 9;
            ans = DoubleFactory2D.dense.make(numFullDOFs,basisSize);
            double CCoeffs[] = pepPlanes[resNum].getProjCAtomCoeffs();
            double OCoeffs[] = pepPlanes[resNum].getProjOAtomCoeffs();
            for(int dim=0; dim<3; dim++){
                for(int at=0; at<3; at++){
                    int pAtomOffset = 3*atomNumN+3*(at+1);
                    if(pAtomOffset < numFullDOFs){
                        ans.set(pAtomOffset+dim, dim, CCoeffs[at]);//C'
                        ans.set(pAtomOffset+dim, dim+3, OCoeffs[at]);//O
                    }
                }
                if(!fixAnchors)
                    ans.set(3*atomNumN+3+dim, dim+6, 1);//CA
            }
        }
        else {//last residue, only N, H and maybe CA
            int basisSize = fixAnchors ? 6 : 9;
            ans = DoubleFactory2D.dense.make(numFullDOFs,basisSize);
            double HCoeffs[] = pepPlanes[resNum-1].getProjHAtomCoeffs();
            for(int dim=0; dim<3; dim++){
                ans.set(3*atomNumN+dim, dim, 1);//N
                for(int at=0; at<3; at++){
                    int pAtomOffset = 3*atomNumN+3*(at-1);
                    if(pAtomOffset < numFullDOFs){
                        ans.set(pAtomOffset+dim, dim+3, HCoeffs[at]);//H
                    }
                }
                if(!fixAnchors)
                    ans.set(3*atomNumN+3+dim, dim+6, 1);//CA
            }
        }
        
        return ans;
    }
    
    
    public DoubleMatrix2D selectResFreeDOFVectors(int resNum){
        //return a set of 6 (column) vectors in free-DOF space that determine,
        //as much as possible, the conf of the specified residue
        DoubleMatrix2D fullVec = Algebra.DEFAULT.transpose(selectResFullDOFVectors(resNum));
        return Algebra.DEFAULT.mult(freeDOFMatrix, fullVec);
    }
    
    public DoubleMatrix2D selectResFullDOFVectors(int resNum){
        //return a set of 6 (row) vectors in full-DOF space that determine,
        //as much as possible, the conf of the specified residue
        
        //collect constrJac (at central conf): needed to get info on how to localize coeffs
        DoubleMatrix2D constrJac = getConstrJac(fullDOFCenter);
        
        //So the space of motions of the N, CA, and C' (and thus CB too) of our residue has rank 6
        //Project each coord of these atoms into the space we are selecting 
        DoubleMatrix2D poseBasis = resPoseBasis(resNum);//full pose basis (as column vectors)
        
        DoubleMatrix2D curSpace = constrJac;//space to orthogonalize to
        int numConstr = curSpace.rows();
        for(int b=0; b<poseBasis.columns(); b++){//iterate through basis of space we want to represent
            DoubleMatrix2D curMz = getOrthogVectors(curSpace);
            DoubleMatrix2D newLocalizedDOF = DoubleFactory2D.dense.make(1,numFullDOFs);
            for(int freeDOF=0; freeDOF<curMz.columns(); freeDOF++){
                double elem = curMz.viewColumn(freeDOF).zDotProduct(poseBasis.viewColumn(b));
                newLocalizedDOF.viewRow(0).assign( curMz.viewColumn(freeDOF), Functions.plusMult(elem) );
            }
            
            if( Algebra.DEFAULT.normF(newLocalizedDOF) > 1e-10){//Not redundant with other localized DOFs
                newLocalizedDOF.assign(Functions.mult(1./Algebra.DEFAULT.normF(newLocalizedDOF)));
                curSpace = DoubleFactory2D.dense.appendRows(curSpace, newLocalizedDOF);
            }
        }
                
        int numLocalizedDOFs;
        if( fixAnchors && (resNum==0||resNum==residues.size()-1) )
            numLocalizedDOFs = 5;
        else
            numLocalizedDOFs = 6;
        
        if( curSpace.rows() != numConstr+numLocalizedDOFs ){
            throw new RuntimeException("ERROR: Expected "+numLocalizedDOFs+" localized DOFs for res, got "
                    +(curSpace.rows()-numConstr));
            //generally there will be 6, but if anchors fixed then first & last res will have five
        }
        
        DoubleMatrix2D resDOFs = curSpace.viewPart(numConstr, 0, numLocalizedDOFs, numFullDOFs);//given as rows here...
        
        //DEBUG!!!
        /*double checkO = Algebra.DEFAULT.normF(Algebra.DEFAULT.mult(constrJac,Algebra.DEFAULT.transpose(resDOFs)));
        DoubleMatrix2D shouldBeDiag = Algebra.DEFAULT.mult(resDOFs,Algebra.DEFAULT.transpose(resDOFs));
        for(int i=0; i<numLocalizedDOFs; i++)
            shouldBeDiag.set(i, i, shouldBeDiag.get(i,i)-1);
        double checkON = Algebra.DEFAULT.normF(shouldBeDiag);*/
        //DEBUG!!!
        
        return resDOFs;//given as rows here...
        //and in full DOF space
        /*DoubleMatrix2D newMz = curSpace.viewPart(numConstr, 0, numFreeDOFs, numFullDOFs);
        newMz = Algebra.DEFAULT.transpose(newMz);
        
        //check newMz is orthonormal and orthogonal 
        //DEBUG!!!
        double checkO = Algebra.DEFAULT.normF(Algebra.DEFAULT.mult(constrJac,newMz));
        DoubleMatrix2D shouldBeDiag = Algebra.DEFAULT.mult(Algebra.DEFAULT.transpose(newMz),newMz);
        for(int i=0; i<numFreeDOFs; i++)
            shouldBeDiag.set(i, i, 0);
        double checkON = Algebra.DEFAULT.normF(shouldBeDiag);*/
    }
    
    
    
    
    
    DoubleMatrix2D selectLocalizableCoeffs(DoubleMatrix2D constrJac, DoubleMatrix2D origMz){
        //return a set of coefficients such that the first 6 determine,
        //as much as possible, the conf of a particular residue
        
        //this particular function is more of a joke, focusing on a peptide plane rather than a res
        //and also hard-coding what peptide plane it is
        //selectResDOFs is the serious version
        int targetRes = 2;
        
        //first attempt: try just projecting the full DOFs into free DOF space
        //we want to see if we can localize to the 9 DOFs of a peptide plane
        //and get good resid for those 9
        //if so then our chosen free dofs are giving effective local control,
        //try to go to 6 etc
        int startFullDOF = 6*targetRes;//start at CA of target res,
        //will do pep plane extending into next res
        int endFullDOF = startFullDOF+9;//end of pep plane
        
        DoubleMatrix2D curSpace = constrJac;//space to orthogonalize to
        int numConstr = curSpace.rows();
        for(int targFullDOF=startFullDOF; targFullDOF<endFullDOF; targFullDOF++){
            DoubleMatrix2D curMz = getOrthogVectors(curSpace);
            DoubleMatrix2D newLocalizedDOF = DoubleFactory2D.dense.make(1,numFullDOFs);
            for(int freeDOF=0; freeDOF<curMz.columns(); freeDOF++){
                double elem = curMz.get(targFullDOF, freeDOF);
                newLocalizedDOF.viewRow(0).assign( curMz.viewColumn(freeDOF), Functions.plusMult(elem) );
            }
            
            if( Algebra.DEFAULT.normF(newLocalizedDOF) > 1e-10){//Not redundant with other localized DOFs
                newLocalizedDOF.assign(Functions.mult(1./Algebra.DEFAULT.normF(newLocalizedDOF)));
                curSpace = DoubleFactory2D.dense.appendRows(curSpace, newLocalizedDOF);
            }
        }
                
        
        //OK now the other DOFs will simply fill in the space left by the localized and other DOFs
        DoubleMatrix2D otherFreeDOFs = getOrthogVectors(curSpace);
        curSpace = DoubleFactory2D.dense.compose(
                new DoubleMatrix2D[][] { new DoubleMatrix2D[] {curSpace},
                    new DoubleMatrix2D[] {Algebra.DEFAULT.transpose(otherFreeDOFs)},
                } );
        
        
        DoubleMatrix2D newMz = curSpace.viewPart(numConstr, 0, numFreeDOFs, numFullDOFs);
        newMz = Algebra.DEFAULT.transpose(newMz);
        
        //check newMz is orthonormal and orthogonal 
        //DEBUG!!!
        double checkO = Algebra.DEFAULT.normF(Algebra.DEFAULT.mult(constrJac,newMz));
        DoubleMatrix2D shouldBeDiag = Algebra.DEFAULT.mult(Algebra.DEFAULT.transpose(newMz),newMz);
        for(int i=0; i<numFreeDOFs; i++)
            shouldBeDiag.set(i, i, 0);
        double checkON = Algebra.DEFAULT.normF(shouldBeDiag);
        
        
        return newMz;
    }
    
    
    
    DoubleMatrix2D selectLocalizableCoeffsOld(DoubleMatrix2D constrJac, DoubleMatrix2D origMz){
        //return a set of coefficients such that the first 6 determine,
        //as much as possible, the conf of a particular residue
        int targetRes = 0;
        
        //first attempt: try just projecting the full DOFs into free DOF space
        //we want to see if we can localize to the 9 DOFs of a peptide plane
        //and get good resid for those 9
        //if so then our chosen free dofs are giving effective local control,
        //try to go to 6 etc
        int startFullDOF = 6*targetRes;//start at CA of target res,
        //will do pep plane extending into next res
        int endFullDOF = startFullDOF+9;//end of pep plane
        
        DoubleMatrix2D localizedDOFs = DoubleFactory2D.dense.make(numFullDOFs,9);
        for(int targFullDOF=startFullDOF; targFullDOF<endFullDOF; targFullDOF++){
            for(int freeDOF=0; freeDOF<numFreeDOFs; freeDOF++){
                double elem = origMz.get(targFullDOF, freeDOF);
                localizedDOFs.viewColumn(targFullDOF).assign( origMz.viewColumn(freeDOF), 
                        Functions.plusMult(elem) );
            }
        }
       
        
        
        //OK now the other DOFs will simply fill in the space left by the localized and other DOFs
        DoubleMatrix2D augCJ = DoubleFactory2D.dense.compose(
                new DoubleMatrix2D[][] { new DoubleMatrix2D[] {Algebra.DEFAULT.transpose(localizedDOFs)},
                    new DoubleMatrix2D[] {constrJac}
                } );
        
        DoubleMatrix2D otherFreeDOFs = getOrthogVectors(augCJ);
        DoubleMatrix2D newMz = DoubleFactory2D.dense.compose(
                new DoubleMatrix2D[][] { new DoubleMatrix2D[] {localizedDOFs, otherFreeDOFs} } );
        
        return newMz;
        //also take a look to see how much the other full DOFs participate in other stuff
        //other possible attempt: assume full control w/i free DOF space
    }
    
    
    public DoubleMatrix1D randomValsInVoxel(){
        DoubleMatrix1D x = DoubleFactory1D.dense.make(numFreeDOFs);
        for(int d=0; d<numFreeDOFs; d++){
            double lb = freeDOFVoxel[0][d];
            double val = lb + Math.random() * (freeDOFVoxel[1][d]-lb);
            x.set(d, val);
        }
        return x;
    }

    
    
}
