/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.plug.RCTuplePolytope;
import edu.duke.cs.osprey.plug.VoxelVDWListChecker;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.awt.Image;
import java.awt.image.BufferedImage;
import static java.awt.image.BufferedImage.TYPE_INT_ARGB;
import static java.awt.image.BufferedImage.TYPE_INT_RGB;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import javax.imageio.ImageIO;

/**
 *
 * Make pictures of clashes wrt first two dihedrals
 * 
 * 
 * @author mhall44
 */
public class ClashFigMaker {
    
    PolytopeMatrix plugMat;
    static final int RESOLUTION = 1000;//per dimension
    
    public ClashFigMaker(PolytopeMatrix plugMat){
        this.plugMat = plugMat;
    }
    
    public static void main(String[] args){
        new ConfigFileParser(new String[] {"-c","KStar.cfg","arsdf"}).loadData();
        ClashFigMaker cfm = new ClashFigMaker( (PolytopeMatrix) ObjectIO.readObject("1CC8.PLUGMAT.dat", true) );
        
        int numToDraw = 6;
        
        int numDrawn = 0;
        int conf[] = new int[] {5,7,7,5,0,7,4};
        for(int pos=0; pos<conf.length; pos++){
            numDrawn += cfm.drawClashFigs(new RCTuple(pos,conf[pos]), "CLASHFIGS."+pos, numToDraw-numDrawn);
            for(int pos2=0; pos2<pos; pos2++){
                numDrawn += cfm.drawClashFigs(new RCTuple(pos,conf[pos],pos2,conf[pos2]), 
                        "CLASHFIGS."+pos+"."+pos2, numToDraw-numDrawn);
            }
        }
    }
    
    private int drawClashFigs(RCTuple tup, String nameStem, int maxToDraw){
        //draw figures for this tuple's clashes, return how many there were
        RCTuplePolytope rtp = plugMat.getTupleValue(tup);
        if(rtp.getDOFs().size()>=2){
            int numToDraw = Math.min(maxToDraw, rtp.getConstr().size()-2*rtp.getDOFs().size());
            for(int clashNum=0; clashNum<numToDraw; clashNum++)
                drawClashFig(tup, clashNum, nameStem+"."+clashNum+".png");
            return numToDraw;
        }
        else
            return 0;
    }
    
    
    private void drawClashFig(RCTuple tup, int clashNum, String fileName){
        //excepts there to be at least clashNum clash constraints, and >=2 cont DOFs
        MoleculeModifierAndScorer mms = new MoleculeModifierAndScorer(null, plugMat.cSpace, tup);
        DoubleMatrix1D constr[] = mms.getConstraints();
        Atom[] clashingAtoms = findClashingAtoms(tup, clashNum);
        
        //have DOFs #2+ stay at center
        mms.setDOFs( constr[0].copy().assign(constr[1],Functions.plus).assign(Functions.mult(0.5)) );
        boolean grid[][] = new boolean[RESOLUTION][RESOLUTION];//true for clashing pixels
        double xWidth = constr[1].get(0)-constr[0].get(0);
        double yWidth = constr[1].get(1)-constr[0].get(1);
        double minDist = VoxelVDWListChecker.getTargetDist(clashingAtoms[0], clashingAtoms[1]);

        for(int xStep=0; xStep<RESOLUTION; xStep++){
            for(int yStep=0; yStep<RESOLUTION; yStep++){
                mms.setDOF(0, constr[0].get(0)+xWidth*xStep/(RESOLUTION-1));
                mms.setDOF(1, constr[0].get(1)+yWidth*yStep/(RESOLUTION-1));
                double dist = VectorAlgebra.distance(clashingAtoms[0].getCoords(), clashingAtoms[1].getCoords());
                grid[xStep][yStep] = (dist<minDist);
            }
        }
        
        drawGrid(grid, fileName);
    }
    
    private Atom[] findClashingAtoms(RCTuple tup, int clashNum){
        //look through the clashing atom pairs mentioned in plugMat for tup,
        //return the clashNum'th
        //It is assumed there is one, will crash cryptically if there isn't
        ArrayList<String> atomPairNames = plugMat.getTupleValue(tup).getAtomPairNames();
        int constrNum=-1;
        for(int c=0; c<clashNum+1; c++){
            do {
                constrNum++;
            } while(atomPairNames.get(constrNum).startsWith("VOXEL CONSTRAINT"));
        }
        String namesSplit[] = atomPairNames.get(constrNum).split(" ; ");
        Atom ans[] = new Atom[2];
        for(int a=0; a<2; a++){
            String split2[] = namesSplit[a].split(" , ");
            ans[a] = plugMat.cSpace.m.getResByFullName(split2[0]).getAtomByName(split2[1]);
        }
        return ans;
    }
    
    private static void drawGrid(boolean[][] grid, String fileName){
        //draw a png of the grid with specified fileName; true=red, false=whites
        BufferedImage bi = new BufferedImage(RESOLUTION, RESOLUTION, TYPE_INT_ARGB);
        int red = -65536;//1 bits for A and R, 0 for G and B
        int white = -1;//all 1 bits
        for(int x=0; x<RESOLUTION; x++){
            for(int y=0; y<RESOLUTION; y++){
                bi.setRGB( x, y, grid[x][y] ? red : white );
            }
        }
        
        try {
            ImageIO.write(bi, "png", new File(fileName));
        }
        catch(IOException e){
            throw new RuntimeException(e.getMessage());
        }
    }
    
}
