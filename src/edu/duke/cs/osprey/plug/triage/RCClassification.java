/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug.triage;

import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import org.apache.commons.math3.optim.linear.LinearConstraint;

/**
 *
 * @author mhall44
 */
public class RCClassification {
    
    HashMap<RCClass.ClassType, ArrayList<RCClass>> classes = new HashMap<>();
    
    HashMap<RC,RCClass> singletonClassLookup = new HashMap<>();//look up the singleton class for an RC
    
    HashSet<String> coreAATypes = new HashSet<>();//AA types associated with the cores of all the classes
    
    
    
    //constructor will build all the RC classifications
    public RCClassification(PositionConfSpace pcSpace, ResidueTypeDOF mutDOF, boolean origBBLinks){
    
        //Each class can be defined by a set of "canonical" DOF values
        //DEBUG!!  for this version, assuming a certain DOF structure
        //normally voxelized FVB backbone classes,
        //and then sidechain dihedrals; only count these if all heavy atoms
        //also note that other kinds of classes, like based on poses of some end groups,
        //could also be beneficial, but might require working outside the tree structure
        
        //If origBBLinks then assume doing FVB and create linked instead of full classes for non-original backbones
        
        //create classes backwards
        //i.e. first create SINGLETON
        //then actually go and match intervals to rc
        //and then to each other
        //ok and then of course sillymark will want to make cores and envelopes
        
        //first make the singleton classes
        ArrayList<RCClass> singletonClasses = new ArrayList<>();
        int numBBVox = CentralStericCache.countBBVox(pcSpace.RCs);
        for(RC rc : pcSpace.RCs){
            if(origBBLinks){
                if(!isBBVoxOrig(rc.bbVoxNum,numBBVox))//skip altered backbones for now...
                    continue;
            }
            RCClass singleton = new RCClass(rc,pcSpace.res,mutDOF);
            singletonClasses.add(singleton);
            singletonClassLookup.put(rc, singleton);
        }
        classes.put(RCClass.ClassType.SINGLETON, singletonClasses);
        
        //OK put singletons in candidate form so we can merge them
        ArrayList<ClassCandidate> curClassification = singletonCandidates(singletonClasses);
        //and create the list of DOFs used to split classes.  Start with CATS + 2 sidechain dihedrals
        HashSet<DegreeOfFreedom> curSplittingDOFs = CATSandSCDOFs(pcSpace, 2);
        
        //now make bigger classes from the smaller candidates, basically removing a DOF each time 
        curClassification = createMergedClasses(curClassification, curSplittingDOFs);
        buildClasses(curClassification, RCClass.ClassType.CHI2CLASS);
        
        curSplittingDOFs = CATSandSCDOFs(pcSpace, 1);
        curClassification = createMergedClasses(curClassification, curSplittingDOFs);
        buildClasses(curClassification, RCClass.ClassType.CHI1CLASS);
        
        curSplittingDOFs = CATSandSCDOFs(pcSpace, 0);
        curClassification = createMergedClasses(curClassification, curSplittingDOFs);
        buildClasses(curClassification, RCClass.ClassType.BBCLASS);//should only move bb
        
        for(RCClass.ClassType type : RCClass.ClassType.values()){
            for(RCClass cl : classes.get(type)){
                coreAATypes.add(cl.core.AAType);
            }
        }
        
        if(origBBLinks){
            makeLinkedClasses(pcSpace);
        }
    }
    
    
    private ArrayList<ClassCandidate> singletonCandidates(ArrayList<RCClass> singletonClasses){
        //To start the merging process, make class candidates from the (completed) candidate classes
        ArrayList<ClassCandidate> ans = new ArrayList<>();
        for(RCClass cl : singletonClasses){
            ans.add(new ClassCandidate(cl));
        }
        return ans;
    }
    
    
    private HashSet<DegreeOfFreedom> CATSandSCDOFs(PositionConfSpace pcSpace, int numSCDihedrals){
        HashSet<DegreeOfFreedom> DOFs = new HashSet<>();
        for(DegreeOfFreedom dof : pcSpace.getAllConfDOFs()){
            if(dof instanceof BBFreeDOF)
                DOFs.add(dof);
            else if(dof instanceof FreeDihedral){
                if( ((FreeDihedral)dof).getDihedralNumber() < numSCDihedrals )
                    DOFs.add(dof);
            }
        }
        return DOFs;
    }
    
    
    private ArrayList<ClassCandidate> createMergedClasses(ArrayList<ClassCandidate> babies, 
            HashSet<DegreeOfFreedom> splittingDOFs){
        //create appropriate parents for the babies
        //based on their intervals for the specified classes
        ArrayList<ClassCandidate> ans = new ArrayList<>();
        for(ClassCandidate baby : babies){
            //try to find a parent for the baby
            boolean foundParent = false;
            for(ClassCandidate possibleParent : ans){
                if(possibleParent.canBeParentOf(baby)){
                    possibleParent.adopt(baby);
                    foundParent = true;
                    break;
                }
            }
            if(!foundParent){//need to create one
                ans.add(new ClassCandidate(baby, splittingDOFs));
            }
        }
        
        return ans;
    }
    
    
    
    private void buildClasses(ArrayList<ClassCandidate> cands, RCClass.ClassType type){
        //build classes from the specified candidates, put them in classes
        ArrayList<RCClass> curTypeClasses = new ArrayList<>();
        for(ClassCandidate cand : cands){
            RCClass newClass = new RCClass(cand.subClasses, type);
            cand.builtClass = newClass;
            curTypeClasses.add(newClass);
            for(RCClass subClass : cand.subClasses)
                subClass.parent = newClass;
        }
        
        classes.put(type, curTypeClasses);
    }
    
    
    private class ClassCandidate {
        //a class will be associated with a certain interval for certain DOFs
        HashMap<DegreeOfFreedom, double[]> intervals = new HashMap<>();//will be used to generate core/envelope
        //and to choose parent for candidate class
        
        ArrayList<RCClass> subClasses = null;
        RCClass builtClass = null;//when we build the class, we'll put it here
        //so class merging can refer to it
        
        //assume linConstr in the various RCs define a set of disparate voxels
        //So the candidate will have a unique set of linConstr
        ArrayList<LinearConstraint> linConstr;
        
        private ClassCandidate(RCClass cl){
            //build a ClassCandidate for cl, to help merge it into bigger classes
            builtClass = cl;//class is already built
            int numDOFs = cl.core.DOFs.size();
            for(int dof=0; dof<numDOFs; dof++){//DOFs in core define class (for singleton these are just RC DOFs and bounds)
                intervals.put( cl.core.DOFs.get(dof), new double[] {cl.core.DOFmin.get(dof),cl.core.DOFmax.get(dof)} );
            }
            linConstr = cl.core.linConstr;
        }
        
        private ClassCandidate(ClassCandidate baby, HashSet<DegreeOfFreedom> splittingDOFs){
            //build a parent for the baby that defines intervals only for the splittingDOFs
            //if any of the splitting DOFs are absent then mark null (so baby will only get siblings with those DOFs also absent)
            for(DegreeOfFreedom dof : baby.intervals.keySet()){
                if(splittingDOFs.contains(dof)){
                    intervals.put(dof, baby.intervals.get(dof));
                }
            }
            
            //add nulls as needed
            for(DegreeOfFreedom dof : splittingDOFs){
                if(!intervals.containsKey(dof)){
                    intervals.put(dof,null);
                }
            }
            
            linConstr = baby.linConstr;
           
            subClasses = new ArrayList(Arrays.asList(baby.builtClass));
            //assuming no built class yet
        }
        
        private boolean canBeParentOf(ClassCandidate baby){
            //Can we add baby as a subclass?
            
            //if( ! constrMatch(linConstr,baby.linConstr) )
            if( subClasses.get(0).core.bbVoxNum != baby.builtClass.core.bbVoxNum )
                return false;
            
            for(DegreeOfFreedom dof : intervals.keySet()){
                if(intervals.get(dof)==null){
                    if(baby.intervals.containsKey(dof)){
                        if(baby.intervals.get(dof)!=null){
                            return false;
                        }
                    }
                }
                else {
                    if( ! baby.intervals.containsKey(dof))
                        return false;
                    if( baby.intervals.get(dof)==null )
                        return false;
                    if( ! intervalsIntersect(baby.intervals.get(dof),intervals.get(dof)) )
                        return false;
                }
            }
            return true;
        }
        
        private void adopt(ClassCandidate baby){
            subClasses.add(baby.builtClass);
            for(DegreeOfFreedom dof : intervals.keySet()){
                intervals.put( dof, intervalUnion(intervals.get(dof),baby.intervals.get(dof)) );
            }
        }
    }
    
    boolean intervalsIntersect(double[] interv1, double[] interv2){
        //we basically demand nonzero intersection
        if(interv1[0]>interv2[1]-1e-14)//no intersection, and interv2 comes first
            return false;
        if(interv2[0]>interv1[1]-1e-14)//no intersection, and interv1 comes first
            return false;
        //else there must be an intersection
        return true;
    }
    
    double[] intervalUnion(double[] interv1, double[] interv2){
        if(interv1==null && interv2==null)
            return null;
        return new double[] {Math.min(interv1[0],interv2[0]), Math.max(interv1[1],interv2[1])};
    }
    
    private boolean constrMatch(ArrayList<LinearConstraint> constr1, ArrayList<LinearConstraint> constr2){
        if(constr1.size() != constr2.size())
            throw new RuntimeException("ERROR: constr sizes don't match fool!!!");
        
        for(int c=0; c<constr1.size(); c++){
            LinearConstraint c1 = constr1.get(c);
            LinearConstraint c2 = constr2.get(c);
            //demand match in form as well as substance of constraints (order, constant factors, etc.)
            if( Math.abs(c1.getValue()-c2.getValue()) > 1e-14 )
                return false;
            if( c1.getRelationship() != c2.getRelationship() )
                return false;
            if( c1.getCoefficients().getDimension() != c2.getCoefficients().getDimension() )
                return false;
            if( c1.getCoefficients().getDistance(c2.getCoefficients()) > 1e-14 )
                return false;
        }
        
        return true;
    }
    
    
    private void makeLinkedClasses(PositionConfSpace pcSpace){
        //we currently have all the orig-BB classes
        //make equivalents for other BB voxels
        //assuming RCs in FVB order (all RCs for first bb vox, for second, ...)
        
        int numOrigBBRCs = singletonClassLookup.size();
        int numRCs = pcSpace.RCs.size();
        int numBBVox = numRCs / numOrigBBRCs;
        
        HashSet<RCClass> origBBClasses = new HashSet<>();
        for(RCClass.ClassType type : RCClass.ClassType.values()){
            for(RCClass origClass : classes.get(type)){
                origBBClasses.add(origClass);
            }
        }
        
        for(int bbVoxNum=0; bbVoxNum<numBBVox; bbVoxNum++){//go over non-orig backbones
            
            if(isBBVoxOrig(bbVoxNum,numBBVox))
                continue;
            
            HashMap<RCClass,RCClass> orig2Altered = new HashMap<>();
            
            //create singleton classes
            for(int rcAtVox=0; rcAtVox<numOrigBBRCs; rcAtVox++){
                int rcNum = rcAtVox*numBBVox+bbVoxNum;
                RC rc = pcSpace.RCs.get(rcNum);
                if(rc.bbVoxNum != bbVoxNum)
                    throw new RuntimeException("ERROR: not seeing expected FVB RC numbering");
                RC origRC = pcSpace.RCs.get(rcAtVox*numBBVox+origBBVoxNum(numBBVox));
                RCClass origClass = singletonClassLookup.get(origRC);
                RCClass newClass = new RCClass(origClass, bbVoxNum);
                singletonClassLookup.put(rc, newClass);
                classes.get(RCClass.ClassType.SINGLETON).add(newClass);
                orig2Altered.put(origClass, newClass);
            }
            
            //create other classes
            for(RCClass origClass : origBBClasses){
                if(origClass.type != RCClass.ClassType.SINGLETON){
                    RCClass newClass = new RCClass(origClass, bbVoxNum);
                    classes.get(newClass.type).add(newClass);
                    orig2Altered.put(origClass, newClass);
                }
            }
            
            //now link up parents and subclasses
            //for current scheme these are within a bbvox
            for(RCClass origClass : origBBClasses){
                RCClass newClass = orig2Altered.get(origClass);
                if(origClass.parent != null)
                    newClass.parent = orig2Altered.get(origClass.parent);
                if(origClass.subclasses != null){
                    newClass.subclasses = new ArrayList<>();
                    for(RCClass origSubclass : origClass.subclasses){
                        newClass.subclasses.add(orig2Altered.get(origSubclass));
                    }
                }
            }
        }
    }
    
    
    
    int numClassesOfType(RCClass.ClassType ctype){
        return classes.getOrDefault(ctype, new ArrayList<>()).size();
    }
    
    RCClass get(RCClass.ClassType ctype, int index){
        //index'th class of given type
        return classes.get(ctype).get(index);
    }
    
    static final boolean isBBVoxOrig(int bbVoxNum, int numBBVox){
        return (bbVoxNum==-1)||(bbVoxNum==origBBVoxNum(numBBVox));
    }
    
    static int origBBVoxNum(int numBBVox){
        //number of bb vox centered at origin
        //for uniform grid, this only makes sense for odd grid widths
        switch(numBBVox){
            case 1:
                return 0;
            case 729:
                return 364;
            case 15625:
                return 7812;
            default:
                throw new RuntimeException("ERROR unsupported grid width");
                //but for numBBVox = a^6 and a odd it should be (a/2)*(1+a+...+a^5)
        }
    }
}

