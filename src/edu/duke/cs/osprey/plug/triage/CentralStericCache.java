/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug.triage;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * @author mhall44
 */
public class CentralStericCache {
    
    
    private HashMap<RCClass,ResSphereModel> origBBSphereModels = new HashMap<>();
    private HashMap<RCClass,ResSphereModel> origBBSphereBubbleModels = new HashMap<>();
    private ArrayList<ResSphereModel> shellSphereModels = new ArrayList<>();
    

    
    //exemptions don't depend on BB conf
    private HashMap<RCClass,HashMap<Residue,ResPairClashExemptions> > shellExemptions = new HashMap<>();
    private HashMap<RCClass,HashMap<Residue,HashMap<String,ResPairClashExemptions> > > pairwiseExemptions = new HashMap<>();
    
    private HashMap<Residue,ArrayList<BBVoxSphereModel> > bbVoxSphereModels = new HashMap<>();//changes to sphere model associated with a bb vox
    
    
    public CentralStericCache(ConfSpace confSpace, ArrayList<RCClassification> classesAtRes, ArrayList<Residue> shellResidues){
        //assuming FVB standard ordering here
        for(int pos=0; pos<confSpace.numPos; pos++){
            ArrayList<RC> RCs = confSpace.posFlex.get(pos).RCs;
            Residue res = confSpace.posFlex.get(pos).res;
            ResidueTypeDOF mutDOF = confSpace.mutDOFs.get(pos);
            int numBBVox = countBBVox(RCs);
            ResSphereModel firstRCModel = null;
            
            //create the orig bb models
            
            for(RCClass.ClassType type : RCClass.ClassType.values()){
                for(RCClass cl : classesAtRes.get(pos).classes.get(type)){
                    if(cl.origBBEquivalent==null){//cl is orig bb
                        if(!GeomPruningStep.applyRCCenter(cl.core, res, mutDOF))
                            throw new RuntimeException("ERROR: Couldn't apply orig-BB RC!");
                        //assuming no intra clashes for rots (else shouldn't be rots)
                        ResSphereModel model = new ResSphereModel(res);
                        origBBSphereModels.put(cl, model);
                        origBBSphereBubbleModels.put(cl, new ResSphereModel(res, cl.pivotAtom, cl.bubbleRad));
                        if(cl.core==RCs.get(RCClassification.origBBVoxNum(numBBVox)))//first orig-bb rc, ie reference for bb transformations
                            firstRCModel = model;

                        computeExemptions(cl, confSpace, classesAtRes);
                    }
                }
            }
            
            if(firstRCModel==null)
                throw new RuntimeException("ERROR: First orig-bb RC has no singleton class!");
            
            
            /*for(int rcNum=RCClassification.origBBVoxNum(numBBVox); rcNum<RCs.size(); rcNum+=numBBVox){
                RC rc = RCs.get(rcNum);
                RCClass origClass = classesAtRes.get(pos).singletonClassLookup.get(rc);
                if(RCClassification.isBBVoxOrig(rc.bbVoxNum,numBBVox)) {//orig bb
                    if(!GeomPruningStep.applyRCCenter(rc, res, mutDOF))
                        throw new RuntimeException("ERROR: Couldn't apply orig-BB RC!");
                    //assuming no intra clashes for rots (else shouldn't be rots)
                    ResSphereModel model = new ResSphereModel(res);
                    origBBSphereModels.put(origClass, model);
                    if(firstRCModel==null)
                        firstRCModel = model;
                    
                    computeExemptions(origClass, confSpace, classesAtRes);
                }
                else {
                    throw new RuntimeException("ERROR: expected orig BB here");
                }
            }
            THIS MISSES THE NON-SINGLETON CLASSES!
            */
            
            
            ArrayList<BBVoxSphereModel> resBBModels = new ArrayList<>();
            
            //create other bb vox models
            for(int bbVoxNum=0; bbVoxNum<numBBVox; bbVoxNum++){
                if(RCClassification.isBBVoxOrig(bbVoxNum,numBBVox)){
                    resBBModels.add(null);//no model needed for orig bb
                }
                else {
                    int rcNum = bbVoxNum;//first RC for the bbVoxNum
                    //assumed to correspond (except for bb motions) to first RC in first bb vox
                    if(GeomPruningStep.applyRCCenter(RCs.get(rcNum), res, mutDOF)){//RC is valid; apply it
                        ResSphereModel altBBModel = new ResSphereModel(res);
                        resBBModels.add(new BBVoxSphereModel(firstRCModel, altBBModel, res.template));
                    }
                    else//not valid; put in null
                        resBBModels.add(null);
                }
            }
            
            bbVoxSphereModels.put(res, resBBModels);
        }
        
        //now handle the shell residues
        for(Residue res : shellResidues){
            shellSphereModels.add(new ResSphereModel(res));
        }
    }
    
    
    static int countBBVox(ArrayList<RC> RCs){
        //assumes FVB ordering
        if(RCs.isEmpty())
            return 0;
        if(RCs.get(0).bbVoxNum==-1)
            return 1;
        for(int rcNum=1; rcNum<RCs.size(); rcNum++){
            if(RCs.get(rcNum).bbVoxNum==0)
                return rcNum;
        }
        throw new RuntimeException("ERROR expected FVB ordering");//shouldn't get here
    }
    
    
    
    private void computeExemptions(RCClass cl, ConfSpace confSpace, ArrayList<RCClassification> classesAtRes){
        //given a class (assumed to have its RC already applied), compute exemptions involving it
        //and add them to shellExemptions and pairwiseExemptions
        pairwiseExemptions.put(cl, new HashMap<>());
        shellExemptions.put(cl, new HashMap<>());
        
        //Exemptions only needed between residues bonded to each other
        //DEBUG!!!  this is true for proteins, in general not clear
        HashSet<Residue> resNeighbors = computeResNeighbors(cl.res);
        for(Residue nres : resNeighbors){
            int npos = confSpace.getDesignIndex(nres);
            if(npos==-1)//nres is shell
                shellExemptions.get(cl).put(nres, new ResPairClashExemptions(cl.res,nres));
            else {
                HashMap<String,ResPairClashExemptions> am = new HashMap<>();
                for( String AAType : classesAtRes.get(npos).coreAATypes ) {
                    if( ! nres.template.name.equalsIgnoreCase(AAType)){
                        confSpace.mutDOFs.get(npos).mutateTo(AAType);
                    }
                    am.put(AAType, new ResPairClashExemptions(cl.res,nres));
                }
                pairwiseExemptions.get(cl).put(nres, am);
            }
        }
    }
    
    
    private HashSet<Residue> computeResNeighbors(Residue res){
        HashSet<Residue> ans = new HashSet<>();
        for(Atom at : res.atoms){
            for(Atom bat : at.bonds){
                if(bat.res!=res){
                    ans.add(bat.res);
                }
            }
        }
        return ans;
    }
            
    
    ResSphereModel getShellSphereModel(int shellPos){
        return shellSphereModels.get(shellPos);
    }
    
    ResSphereModel getSphereModel(RCClass cl){
        return getSphereModel(cl,false);
    }
    
    ResSphereModel getSphereModel(RCClass cl, boolean includeBubble){
        if(cl.origBBEquivalent==null){//"original" BB voxel
            return getModelOrigBB(cl,includeBubble);
        }
        else {
            RCClass clOrigBB = cl.origBBEquivalent;
            ResSphereModel origBBModel = getModelOrigBB(clOrigBB,includeBubble);
            BBVoxSphereModel bbModel = bbVoxSphereModels.get(clOrigBB.res).get(cl.bbVoxNum);
            
            if(bbModel==null)//conf impossible.  (origBBModel should be possible)
                return null;
            
            return bbModel.transform(origBBModel, clOrigBB.core.template);
        }
    }
    
    private ResSphereModel getModelOrigBB(RCClass cl, boolean includeBubble){
        if(includeBubble)
            return origBBSphereBubbleModels.get(cl);
        else
            return origBBSphereModels.get(cl);
    }
    
    private RCClass origBBVersion(RCClass cl){
        if(cl.origBBEquivalent==null)
            return cl;
        else
            return cl.origBBEquivalent;
    }
    
    HashMap<Residue,ResPairClashExemptions> getShellExemptions(RCClass cl){
        return shellExemptions.get(origBBVersion(cl));
    }
    
    
    ResPairClashExemptions getRCPairExemptions(RCClass cl1, RCClass cl2){
        HashMap<Residue,HashMap<String,ResPairClashExemptions> > resExemptions = pairwiseExemptions.get(origBBVersion(cl1));
         cl2 = origBBVersion(cl2);
        if(resExemptions.containsKey(cl2.res)){//cl1, cl2 residues bonded, exemptions possible
            String resType2 = cl2.core.AAType;
            return resExemptions.get(cl2.res).get(resType2);
        }
        else
            return null;//indicates no exemptions (generally the case if cl1, cl2 not consecutive res)
    }
    
}
