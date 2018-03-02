/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug.triage;

import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.plug.RCTuplePolytope;
import static edu.duke.cs.osprey.plug.VoxelVDWDistExplorer.getVDWRadius;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.util.ArrayList;

/**
 *
 * A class of RCs
 * 
 * @author mhall44
 */
public class RCClass {
    
    
    static enum ClassType {
        //DEBUG!!!  for now just the one main hierarchy
        //ordering from less to more numerous (for efficiency in pruning
        BBCLASS, CHI1CLASS, CHI2CLASS, SINGLETON
    }
    
    Residue res = null;
    ResidueTypeDOF mutDOF = null;
    
    RC core = null;//this need not be a real RC in the confspace
    //but basically the class can be pruned by clash or Rama if the "core" can
    //represents the steric presence that all RCs in the class share
    
    //ok there also needs to be an envelope
    //but we'll just define this as core + "bubble" around "pivot" atom
    Atom pivotAtom = null;
    double bubbleRad = 0;
    
    ArrayList<RCClass> subclasses = null;//null if this is a singleton class
    
    RCClass parent = null;//parent (superclass) of this class 
    ClassType type;
    
    RCClass origBBEquivalent = null;//for linked BB class setup.  In this case core, etc. likely unnecessary
    int bbVoxNum;//origBBEquivalent and bbVoxNum together actually suffice to define an altered-BB class
    
    public RCClass(RC rc, Residue re, ResidueTypeDOF md){
        //build a singleton class from an RC
        core = rc;
        type = ClassType.SINGLETON;
        mutDOF = md;
        res = re;
        bbVoxNum = core.bbVoxNum;
        core.template = EnvironmentVars.resTemplates.getTemplateForMutation(core.AAType, res, true);
    }
    
    public RCClass(ArrayList<RCClass> subClasses, ClassType type){
        this.subclasses = subClasses;
        this.type = type;
        
        res = subClasses.get(0).res;
        mutDOF = subClasses.get(0).mutDOF;
        for(RCClass cl : subClasses){
            if(cl.res!=res){
                throw new RuntimeException("ERROR: Can't merge RC classes for different res");
            }
        }
        
        
        CoreEnvelopePooler cep = new CoreEnvelopePooler();
        for(RCClass sc : subClasses){
            ResidueTemplate templ = EnvironmentVars.resTemplates.getTemplateForMutation(sc.core.AAType, res, true);
            cep.poolTemplate(templ, sc.pivotAtom, sc.bubbleRad);
        }
        
        //DEBUG!!! using simple core and envelope forms, may refine later
        //OK actually we won't have envelope RC just do core + a radius!!!
        core = truncateRC(subClasses.get(0).core, cep.coreTemplate());
        bbVoxNum = core.bbVoxNum;
        
        pivotAtom = choosePivotAtom(core.template,type);
        bubbleRad = cep.bubbleRad;
    }
    
    
    public RCClass(RCClass origBBEquivalent, int bbVoxNum){
        //create an altered-BB equivalent to origBBEquiv
        //parents and subclasses left null now, will be put in once
        //create these equivalents for all orig-BB classes at the given bbVoxNum
        this.origBBEquivalent = origBBEquivalent;
        this.bbVoxNum = bbVoxNum;
        type = origBBEquivalent.type;
        //core, etc. assumed unneeded, left null
    }
    
    
    private Atom choosePivotAtom(ResidueTemplate templ, ClassType type){
        //choose atom as center of bubble
        if(type==ClassType.BBCLASS)
            return templ.templateRes.getAtomByName("CA");
        
        int dihNum = (type==ClassType.CHI2CLASS) ? 1 : 0;
        if(dihNum>=templ.numDihedrals)//templ doesn't have this dihedral, so no pivoting and no bubble needed
            return null;
        
        return templ.templateRes.atoms.get( templ.dihedral4Atoms[dihNum][2] );
    }
    
    private class CoreEnvelopePooler {
        //pools information on subclass cores and envelopes
        //so that we can make an overall core and envelope
        //DEBUG!!!  Some of these assumptions may need to be modified for some
        //not-natural AA res
        //At least, we need to ensure that when a res is fully defined by the splitting DOFs, we choose core/envelope accordingly
        //(should be true for natural AA)
        
        //some info about what atoms all the envelope has
        boolean allHaveCBeta = true;//CB is afaik always SP3, so sterically >=ALA
        boolean allHaveTwoCGammas = true;//sterically >= VAL
        boolean allHaveBenzRing = true;
        boolean allHaveChi1 = true;
        boolean allHaveChi2 = true;
        
        //DEBUG!!!  This could be a loose bound, but should only be especially loose
        //for small classes (where not much extends beyond pivot atom).  
        //Just bound sidechain (expect for CG's, CB's, and their H's)
        //with a sphere around pivot atom
        double bubbleRad = 0;
        Atom pivotAtom = null;
        
        
        void poolTemplate(ResidueTemplate templ, Atom subPivotAtom, double subBubbleRad){
            //process templ for inclusion in pool
            //DEBUG!!!  pretty aa-library specific
            if(templ.templateRes.getAtomByName("CB")==null)
                allHaveCBeta = false;
            if(templ.templateRes.getAtomByName("CG1")==null)
                allHaveTwoCGammas = false;
            if( ! (templ.name.equalsIgnoreCase("PHE")||templ.name.equalsIgnoreCase("TYR")) )
                allHaveBenzRing = false;
            if(templ.numDihedrals<2)
                allHaveChi2 = false;
            if(templ.numDihedrals<1)
                allHaveChi1 = false;
            
            //now figure out bubble radius
            /*int dihNum = (type==ClassType.CHI2CLASS) ? 1 : 0;
            if(type==ClassType.BBCLASS)
                pivotAtom = templ.templateRes.getAtomByName("CA");
            else
                pivotAtom = templ.templateRes.atoms.get( templ.dihedral4Atoms[dihNum][2] );
            */
            Atom curPivotAtom = choosePivotAtom(templ, type);
            if(curPivotAtom!=null){
            
                for(Atom at : templ.templateRes.atoms){
                    if( ! HardCodedResidueInfo.possibleBBAtomsLookup.contains(at.name) ){
                        if( at.name.contains("HA") )
                            continue;
                        if(type==ClassType.CHI1CLASS || type==ClassType.CHI2CLASS){
                            if( at.name.contains("CB")||at.name.contains("HB") )
                                continue;
                        }
                        if(type==ClassType.CHI2CLASS){
                            if( at.name.contains("CG")||at.name.contains("HG") )
                                continue;
                        }

                        double dist = VectorAlgebra.distance(curPivotAtom.getCoords(), at.getCoords());
                        bubbleRad = Math.max(bubbleRad, dist+getVDWRadius(at));
                    }
                }

                //finally make sure our bubble includes the subclass bubble
                if(subPivotAtom!=null){//there is a subclass bubble
                    double pivotAtDist = VectorAlgebra.distance(curPivotAtom.getCoords(), subPivotAtom.getCoords());
                    //use an upper-bound on how big we need to make the current bubble to include the subclass bubble
                    bubbleRad = Math.max(bubbleRad, pivotAtDist+subBubbleRad);
                }
            }
            
            //DEBUG!!!  ideally bubble should be on actual conf not template
        }
        
        
        ResidueTemplate coreTemplate(){
            return EnvironmentVars.resTemplates.getTemplateForMutation(coreTemplateName(), res, true);
        }
        
        String coreTemplateName(){
            //choose an appropriate core template, based on collected subclass data
            switch(type){
                case CHI2CLASS:
                    if(allHaveTwoCGammas)
                        return "ILE";
                    if(allHaveBenzRing)
                        return "PHE";
                    else if(allHaveChi2)
                        return "PLA";//truncate to a bare CD
                    //else fall through, treat as chi1 class, which it's equivalent to
                case CHI1CLASS:
                    if(allHaveTwoCGammas)
                        return "VAL";
                    else if(allHaveChi1)
                        return "ELA";//truncate to a bare CG
                case BBCLASS:
                    if(allHaveCBeta)
                        return "ALA";
                    else
                        return "GLY";
                default:
                    throw new RuntimeException("ERROR: Class type not supported: "+type);
            }
        }
    }
    
    
    private RC truncateRC(RC rc, ResidueTemplate truncTemplate){
        //cut off rc to the part that depends only on the specified class type's DOFs
        //so that restraining those DOFs will let us prune classes of that type sterically
        //just by examining the truncated part
        
        //Let's truncate with a hydrogenated (e.g. methyl) cap
        //create new template if needed
        
        //let's make a cache of truncated templates that we can use
        //then can just check equality of pointers to see if need to mutate
        //they'll all be called TRU (for truncated) but RC applier
        //DEBUG!!!  this will not interact properly with normal mutator, only can be used in screener
        //(though normal mutator will correctly replace any TRU with the normal AA it needs)
        ArrayList<DegreeOfFreedom> truncDOFs = new ArrayList<>();
        ArrayList<Double> truncDOFmax = new ArrayList<>();
        ArrayList<Double> truncDOFmin = new ArrayList<>();
        
        for(int dofNum=0; dofNum<rc.DOFs.size(); dofNum++){
            if(DOFMatchesType(rc.DOFs.get(dofNum), type)){
                truncDOFs.add(rc.DOFs.get(dofNum));
                truncDOFmin.add(rc.DOFmin.get(dofNum));
                truncDOFmax.add(rc.DOFmax.get(dofNum));
            }
        }
        
        RC newRC = new RC(truncTemplate.name, truncTemplate, -1, truncDOFs, truncDOFmin, truncDOFmax, -1);
        newRC.linConstr = RCTuplePolytope.transferConstraints(rc.DOFs, truncDOFs, rc.linConstr);//since classes don't split bb space (yet)
        newRC.bbVoxNum = rc.bbVoxNum;
        return newRC;
    }
    
    private boolean DOFMatchesType(DegreeOfFreedom dof, ClassType type){
        //Does this type of class restrict this dof?
        if(type==ClassType.SINGLETON)
            return true;
        if(dof instanceof BBFreeDOF)
            return true;
        
        if(dof instanceof FreeDihedral){
            int dihNum = ((FreeDihedral)dof).getDihedralNumber();
            switch(type){
                case BBCLASS:
                    return false;
                case CHI1CLASS:
                    return dihNum<1;
                case CHI2CLASS:
                    return dihNum<2;
                default:
                    throw new RuntimeException("ERROR: Unsupported ClassType "+type);
            }
        }
        
        return false;//other types not restricted currently
    }
    
    /*static final HashMap<ClassType,HashMap<String,ResidueTemplate>> truncTemplateCache = new HashMap<>();
    static {
        for(ClassType type : ClassType.values()){
            truncTemplateCache.put(type, new HashMap<>());
        }
    }*/
    
}
