package edu.duke.cs.osprey.kstar;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.IntStream;

import edu.duke.cs.osprey.kstar.impl.KSImplKAStar;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFTraditional;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFUB;

/*
 * TODO
 * 1) only do upper bound stability check for depth 1?
 * 2) am i re-computing unbound pfs for tight bound leaf?
 */

public class KAStarNode {

	private static KSCalc wt;
	private static HashMap<Integer, KSAllowedSeqs> strand2AllowedSeqs;
	private static KSAbstract ksObj;
	private static int numCreated = 0;
	private static int numExpanded = 0;
	private static int numPruned = 0;
	public static int numLeavesCreated = 0;
	public static int numLeavesCompleted = 0;

	public KSCalc ub;
	public KSCalc lb;
	protected double parentlbScore; //parent lower bound
	protected double lbScore;
	protected double ubScore;
	protected boolean scoreNeedsRefinement;


	public static void init( KSAbstract ksObj, 
			HashMap<Integer, KSAllowedSeqs> strand2AllowedSeqs, KSCalc wt ) {

		KAStarNode.ksObj = ksObj;
		KAStarNode.strand2AllowedSeqs = strand2AllowedSeqs;
		KAStarNode.wt = wt;
	}


	public KAStarNode( KSCalc lb, KSCalc ub, boolean scoreNeedsRefinement ) {
		this.lb = lb;
		this.ub = ub;
		parentlbScore = Double.NEGATIVE_INFINITY;
		lbScore = Double.NEGATIVE_INFINITY;
		ubScore = Double.POSITIVE_INFINITY;
		this.scoreNeedsRefinement = scoreNeedsRefinement;
	}


	public static int getNumLeavesCreated() {
		return numLeavesCreated;
	}


	public static int getNumLeavesCompleted() {
		return numLeavesCompleted;
	}


	public static int getNumCreated() {
		return numCreated;
	}


	public static int getNumExpanded() {
		return numExpanded;
	}


	public static int getNumPruned() {
		return numPruned;
	}


	// only expand if scoreNeedsRefinement
	public ArrayList<KAStarNode> expand() {

		ArrayList<KAStarNode> children = null;

		if( isFullyDefined() ) {

			if( !scoreNeedsRefinement() ) {

				if( !isFullyProcessed() ) {
					// we will use tight bounds; compute tight bound		
					if(this.lb.canContinue())
						computeScore(this);
				}

				children = new ArrayList<>();
				children.add(this);
				pruneIncompatibleSuccessors(children, null);

				return children;
			}
		}

		numExpanded++;

		// using complex, p, l convention
		ArrayList<Integer> strands = new ArrayList<>(Arrays.asList(KSTermini.LIGAND, KSTermini.PROTEIN, KSTermini.COMPLEX));

		ArrayList<Integer> nextDepths = new ArrayList<>(); 
		for( int i = 0; i < strands.size(); ++i ) nextDepths.add(0);

		ArrayList<ArrayList<String>> seqs = new ArrayList<>(); 
		for( int i = 0; i < strands.size(); ++i ) seqs.add(new ArrayList<String>());

		KSAllowedSeqs pSeqs = strand2AllowedSeqs.get(KSTermini.PROTEIN);
		KSAllowedSeqs lSeqs = strand2AllowedSeqs.get(KSTermini.LIGAND);

		// get next depths for each strand
		for( int strand : strands ) {

			if( depth() != 0 ) {
				ArrayList<String> seq = lb.getPF(strand).getSequence();
				seqs.set(strand, seq);

				if( strand != KSTermini.COMPLEX )
					nextDepths.set(strand, Math.min(seq.size()+1, strand2AllowedSeqs.get(strand).getStrandSubSeqsMaxDepth()));
			}

			else {
				KSAllowedSeqs strandSeqs = strand2AllowedSeqs.get(strand);
				// go to the first non-zero depth
				int depth = 1;

				// set next depth
				if( strand != KSTermini.COMPLEX ) {
					// prepare subsequences
					strandSeqs.getStrandSubSeqsAtDepth( depth );
					nextDepths.set(strand, depth);
				}

				else
					// prepare subsequences
					strandSeqs.getStrandSubSeqsAtDepth( depth, pSeqs, lSeqs );
			}

			if( strand == KSTermini.COMPLEX )
				nextDepths.set(strand, nextDepths.get(KSTermini.PROTEIN) + nextDepths.get(KSTermini.LIGAND));
		}

		// get all sequences at next depth
		HashSet<ArrayList<String>> nextPSeqs = new HashSet<>(strand2AllowedSeqs.get(KSTermini.PROTEIN).getStrandSubSeqsAtDepth(nextDepths.get(KSTermini.PROTEIN)));
		HashSet<ArrayList<String>> nextLSeqs = new HashSet<>(strand2AllowedSeqs.get(KSTermini.LIGAND).getStrandSubSeqsAtDepth(nextDepths.get(KSTermini.LIGAND)));
		HashSet<ArrayList<String>> nextPLSeqs = new HashSet<>(strand2AllowedSeqs.get(KSTermini.COMPLEX).getStrandSubSeqsAtDepth(nextDepths.get(KSTermini.COMPLEX)));

		// in general, we want to add a single residue to the current protein and ligand
		// so our children are valid combinations of nextProtein x nextLigand.
		// however, from the root node, all reachable nodes are valid.

		if( depth() != 0 ) {
			// reduce next L to those subsequences reachable from this
			ArrayList<String> currentLSeq = seqs.get(KSTermini.LIGAND);
			for( Iterator<ArrayList<String>> iterator = nextLSeqs.iterator(); iterator.hasNext(); ) {
				ArrayList<String> subSeq = iterator.next();	    
				if( !subSeq.subList(0, currentLSeq.size()).equals(currentLSeq) ) iterator.remove();
			}

			// reduce next P to those subsequences reachable from this
			ArrayList<String> currentPSeq = seqs.get(KSTermini.PROTEIN);
			for( Iterator<ArrayList<String>> iterator = nextPSeqs.iterator(); iterator.hasNext(); ) {
				ArrayList<String> subSeq = iterator.next();
				if( !subSeq.subList(0, currentPSeq.size()).equals(currentPSeq) ) iterator.remove();
			}
		}

		children = getChildren( nextPLSeqs, nextPSeqs, nextLSeqs );

		if(children.size() > 0) {

			// compute scores
			switch( KSImplKAStar.nodeExpansionMethod ) {

			case "serial":
				computeScoresSerial(children);
				break;

			case "parallel2":
				computeScoresComplexParallel(children);
				break;

			default:
				computeScoresSimpleParallel(children);
				break;
			}

			pruneIncompatibleSuccessors( children, nextDepths );
			//doIntermutationPruning(children);
		}

		return children;
	}


	private void pruneIncompatibleSuccessors( ArrayList<KAStarNode> children, ArrayList<Integer> nextDepths ) {
		// remove children whose upper bound unbound states are not stable wrt wt
		for( Iterator<KAStarNode> iterator = children.iterator(); iterator.hasNext(); ) {

			KAStarNode child = iterator.next();

			// true is not possible, because we are not computing the bound state. this will catch not_stable and not_possible
			if( child.ub != null && child.ub.getEpsilonStatus() != EApproxReached.FALSE ) {

				// remove sequences that contain no valid conformations
				for( int strand : Arrays.asList(KSTermini.LIGAND, KSTermini.PROTEIN) ) {
					PFAbstract pf = child.ub.getPF(strand);

					if( pf.getEpsilonStatus() == EApproxReached.NOT_POSSIBLE && (pf.getQStar().add(pf.getQPrime().add(pf.getPStar()))).compareTo(BigDecimal.ZERO) > 0 ) {
						System.out.println("ERROR: " + KSAbstract.list1D2String( pf.getSequence(), " " )  + " " + pf.getFlexibility() + " is wrong!!!");
					}

					if( pf.getEpsilonStatus() == EApproxReached.NOT_STABLE || 
							pf.getNumUnPruned().compareTo(BigInteger.ZERO) == 0 ) {

						ArrayList<String> seq = pf.getSequence();

						for( int strand2 : Arrays.asList(strand, KSTermini.COMPLEX) ) {
							if(nextDepths != null) pruneSequences(seq, strand2, nextDepths.get(strand2));
						}
					}
				}

				// epsilon is not going to be true, since we do not compute complex
				// epsilon cannot be not possible or not stable
				iterator.remove();
			}

			else if( !child.lb.doingKAStar() && !child.lb.canContinue() ) {
				iterator.remove();
			}
		}
	}


	/*
	private void doIntermutationPruning( ArrayList<KAStarNode> children ) {

		if(!KSImplKAStar.useTightBounds) 
			return;

		if(KSAbstract.interMutationConst <= 0.0) 
			return;

		KSCalc best = ksObj.getBestCalc();
		if(best == null) 
			return;

		for( Iterator<KAStarNode> iterator = children.iterator(); iterator.hasNext(); ) {

			KAStarNode child = iterator.next();

			// applies to partially or fully defined
			if(child.lb.getEpsilonStatus() != EApproxReached.TRUE || 
					child.lb.getEpsilonStatus() != EApproxReached.FALSE) continue;

			// intermutation pruning criterion
			if(ksObj.passesInterMutationPruning(child.lb))
				continue;

			else {
				if(child.lb.getEpsilonStatus() == EApproxReached.FALSE) {
					// this block only applies to leaf nodes where we are minimizing confs
					// protein and ligand strands are complete by default, so only abort complex
					PFAbstract pf = child.lb.getPF(Termini.COMPLEX);
					pf.abort(true);
					pf.cleanup();
				}

				iterator.remove();
			}
		}
	}
	 */


	private void pruneSequences( ArrayList<String> seq, int strand, int depth ) {
		HashSet<ArrayList<String>> set = null;
		int maxDepth = strand2AllowedSeqs.get(strand).getStrandSubSeqsMaxDepth();

		for( ; depth <= maxDepth; ++depth) {
			set = strand2AllowedSeqs.get(strand).getStrandSubSeqsAtDepth(depth);
			int oldSize = set.size();

			KSAllowedSeqs.deleteFromSet(seq, set);
			numPruned += (oldSize - set.size());
		}
	}


	private void computeScore(KAStarNode child) {

		if( child.scoreNeedsRefinement() ) {

			PFAbstract.suppressOutput = true;

			// compute protein and ligand upper bound pfs; prune if not stable
			// for stable pfs, compute lower bound
			for( int strand : Arrays.asList(KSTermini.LIGAND, KSTermini.PROTEIN) ) {
				if( child.ub.getEpsilonStatus() == EApproxReached.FALSE ) {
					child.ub.runPF(child.ub.getPF(strand), wt.getPF(strand), true, true);
				}
			}

			if(child.ub.getEpsilonStatus() != EApproxReached.FALSE) {
				PFAbstract.suppressOutput = false;
				return;
			}

			child.lb.run(wt, true, false);
			child.lbScore = -1.0 * child.lb.getKStarScoreLog10(true);

			checkConsistency(child);

			PFAbstract.suppressOutput = false;
		}

		else {
			KSAbstract.doCheckPoint = KSImplKAStar.useTightBounds ? true : false;

			// we process n-body leaf nodes as streams (in the complex case, anyway)
			boolean lbStabilityCheck = KSAbstract.doCheckPoint;

			if(!KSImplKAStar.useTightBounds) { // we are not doing n body minimization
				// compute protein and ligand upper bound pfs; prune if not stable
				// for stable pfs, compute lower bound
				for( int strand : Arrays.asList(KSTermini.LIGAND, KSTermini.PROTEIN) ) {
					if( child.ub.getEpsilonStatus() == EApproxReached.FALSE ) {
						child.ub.runPF(child.ub.getPF(strand), wt.getPF(strand), true, true);
					}
				}

				if(child.ub.getEpsilonStatus() != EApproxReached.FALSE) {
					KSAbstract.doCheckPoint = false;
					return;
				}
			}

			child.lb.run(wt, false, lbStabilityCheck);

			if( !child.lb.doingKAStar() && !child.lb.canContinue() ) { // epsilon is not possible or not stable
				return;
			}

			if( child.lb.getEpsilonStatus() == EApproxReached.TRUE ) {
				numLeavesCompleted = ksObj.getNumSeqsCompleted(1);
			}

			child.lbScore = -1.0 * child.lb.getKStarScoreLog10(true);

			KSAbstract.doCheckPoint = false;
		}
	}


	private class CalcParams {
		int strand; boolean complete; boolean stabilityCheck;
		public CalcParams( int strand, boolean complete, boolean stabilityCheck ) {
			this.strand = strand; this.complete = complete; this.stabilityCheck = stabilityCheck;
		}
	}


	protected void computeScoresSerial(ArrayList<KAStarNode> children) {

		if(children.size() == 0) return;
		for(KAStarNode child : children) {
			if(!child.isFullyProcessed())
				computeScore(child);
		}
	}


	protected void computeScoresSimpleParallel(ArrayList<KAStarNode> children) {

		if(children.size() == 0) return;

		boolean parallel = false;

		if(canParallelize(children))
			parallel = true;

		if(parallel) {
			children.parallelStream().forEach( child -> {
				if(!child.isFullyProcessed()) 
					computeScore(child);
			});
		}

		else
			computeScoresSerial(children);
	}


	protected void computeScoresComplexParallel( ArrayList<KAStarNode> children ) {

		if( children.size() == 0 ) return;

		// computes unique partition functions rather than aggregating as calculations
		if( children.get(0).scoreNeedsRefinement() ) {
			PFAbstract.suppressOutput = true;

			// get upper bound protein and ligands
			ArrayList<KSCalc> calcs = new ArrayList<>();
			ArrayList<CalcParams> params = new ArrayList<>();

			for( int strand : Arrays.asList(KSTermini.LIGAND, KSTermini.PROTEIN) ) {
				ArrayList<KSCalc> ans = getCalcs4Strand( children, false, strand );
				for(KSCalc calc : ans) {
					calcs.add(calc);
					params.add(new CalcParams(strand, true, true));
				}
			}
			// run all ub pfs
			IntStream.range(0, calcs.size()).parallel().forEach(i -> {
				KSCalc calc = calcs.get(i);
				CalcParams param = params.get(i);
				calc.runPF(calc.getPF(param.strand), wt.getPF(param.strand), param.complete, param.stabilityCheck);
			});

			// clear calcs, params
			calcs.clear();
			params.clear();

			ArrayList<KAStarNode> children2 = new ArrayList<>();

			// select, retain viable children for pruning
			// we don't do bound state here, so true cannot happen. false is a filter for not_possible
			for( KAStarNode child : children ) {
				if(child.ub.getEpsilonStatus() == EApproxReached.FALSE)
					children2.add(child);
			}

			// for viable children, run all pfs
			for( int strand : Arrays.asList(KSTermini.LIGAND, KSTermini.PROTEIN, KSTermini.COMPLEX) ) {
				ArrayList<KSCalc> ans = getCalcs4Strand( children2, true, strand );
				for(KSCalc calc : ans) {
					calcs.add(calc);
					params.add(new CalcParams(strand, true, false));
				}
			}
			IntStream.range(0, calcs.size()).parallel().forEach(i -> {
				KSCalc calc = calcs.get(i);
				CalcParams param = params.get(i);
				calc.runPF(calc.getPF(param.strand), wt.getPF(param.strand), param.complete, param.stabilityCheck);
			});

			// set scores for surviving children
			for( KAStarNode child : children2 ) {
				child.lbScore = -1.0 * child.lb.getKStarScoreLog10(true);

				checkConsistency(child);
			}

			PFAbstract.suppressOutput = false;
		}

		else
			computeScoresSerial(children);
	}


	private ArrayList<KSCalc> getCalcs4Strand( ArrayList<KAStarNode> children, boolean lb, int strand ) {

		// get minimal set of kscalcs that covers all desired pfs 
		HashMap<ArrayList<String>, KSCalc> seq2KSCalc = new HashMap<>();

		for( KAStarNode child : children ) {
			KSCalc calc = lb ? child.lb : child.ub;
			PFAbstract pf = calc.getPF(strand);
			seq2KSCalc.put(pf.getSequence(), calc);
		}

		ArrayList<KSCalc> calcs = new ArrayList<>(seq2KSCalc.values());
		return calcs;
	}


	private boolean canParallelize(ArrayList<KAStarNode> children) {		
		if(children.size() < 2) 
			return false;

		HashSet<String> setPL = new HashSet<>(children.size());
		HashSet<String> setP = new HashSet<>(children.size());
		HashSet<String> setL = new HashSet<>(children.size());

		for(KAStarNode child : children) {
			String plSeq = KSAbstract.list1D2String(child.lb.getPF(KSTermini.COMPLEX).getSequence(), " ");
			String pSeq = KSAbstract.list1D2String(child.lb.getPF(KSTermini.PROTEIN).getSequence(), " ");
			String lSeq = KSAbstract.list1D2String(child.lb.getPF(KSTermini.LIGAND).getSequence(), " ");

			if(setPL.contains(plSeq)) return false;
			setPL.add(plSeq);

			if(setP.contains(pSeq)) return false;
			setP.add(pSeq);

			if(setL.contains(lSeq)) return false;
			setL.add(lSeq);
		}

		return true;
	}


	private ArrayList<KAStarNode> getChildren( HashSet<ArrayList<String>> nextPLSeqs,
			HashSet<ArrayList<String>> pSeqs, HashSet<ArrayList<String>> lSeqs ) {

		ArrayList<ArrayList<String>> strandSeqs = new ArrayList<>(Arrays.asList(null, null, null));
		// lower bound pf values
		ArrayList<Boolean> lbContSCFlexVals = new ArrayList<>(Arrays.asList(false, false, true));
		ArrayList<String> lbPFImplVals = new ArrayList<>(Arrays.asList(new PFTraditional().getImpl(), new PFTraditional().getImpl(), new PFUB().getImpl()));

		// upper bound pf values
		ArrayList<Boolean> ubContSCFlexVals = new ArrayList<>(Arrays.asList(true, true, false));
		ArrayList<String> ubPFImplVals = new ArrayList<>(Arrays.asList(new PFUB().getImpl(), new PFUB().getImpl(), new PFTraditional().getImpl()));

		// minimized pf values
		ArrayList<Boolean> tightContSCFlexVals = new ArrayList<>(Arrays.asList(true, true, true));
		ArrayList<String> tightPFImplVals = new ArrayList<>(Arrays.asList(PFAbstract.getCFGImpl(), PFAbstract.getCFGImpl(), PFAbstract.getCFGImpl()));

		ArrayList<KAStarNode> ans = new ArrayList<>();

		for( ArrayList<String> pSeq : pSeqs ) {

			for( ArrayList<String> lSeq : lSeqs ) {

				ArrayList<String> putativeNextPLSeq = new ArrayList<String>();

				putativeNextPLSeq.addAll(pSeq);
				putativeNextPLSeq.addAll(lSeq);
				putativeNextPLSeq.trimToSize();

				if( !nextPLSeqs.contains(putativeNextPLSeq) ) continue;

				numCreated++;

				// create partition functions for next sequences
				strandSeqs.set(KSTermini.COMPLEX, putativeNextPLSeq);
				strandSeqs.set(KSTermini.PROTEIN, pSeq);
				strandSeqs.set(KSTermini.LIGAND, lSeq);

				if( !isFullyDefined() ) {
					// create partition functions
					ConcurrentHashMap<Integer, PFAbstract> lbPFs = ksObj.createPFs4Seqs(strandSeqs, lbContSCFlexVals, lbPFImplVals);
					ConcurrentHashMap<Integer, PFAbstract> ubPFs = ksObj.createPFs4Seqs(strandSeqs, ubContSCFlexVals, ubPFImplVals);

					// create KUStar node with lower and upper bounds
					ans.add( new KAStarNode( new KSCalc(numCreated, lbPFs), new KSCalc(numCreated, ubPFs), childScoreNeedsRefinement(lbPFs) ) );
					ans.get(ans.size()-1).parentlbScore = this.lbScore;
				}

				else if( KSImplKAStar.useTightBounds ) {
					// create a leaf node; we don't need upper or lower bounds ;we only need the minimized partition 
					// function our search problems exist, so we need only delete the lb, ub pfs from the table

					ConcurrentHashMap<Integer, PFAbstract> tightPFs = ksObj.createPFs4Seqs(strandSeqs, tightContSCFlexVals, tightPFImplVals);

					// assign sequence number from allowedSequences obj
					KSAllowedSeqs complexSeqs = ksObj.strand2AllowedSeqs.get(KSTermini.COMPLEX);
					int seqID = complexSeqs.getPosOfSeq(tightPFs.get(KSTermini.COMPLEX).getSequence());

					// create new KUStar node with tight score
					ans.add( new KAStarNode( new KSCalc(seqID, tightPFs), null, false ) );
					ans.get(ans.size()-1).parentlbScore = this.lbScore;

					// processed nodes are those whose confs will be minimized
					numLeavesCreated = ksObj.getNumSeqsCreated(1);
				}

				else
					throw new RuntimeException("ERROR: cannot expand a fully assigned node");
			}
		}

		return ans;
	}


	private int depth() {
		if( lb == null ) return 0;
		// number of assigned residues is depth
		return lb.getPF(KSTermini.COMPLEX).getSequence().size();
	}


	public boolean scoreNeedsRefinement() {
		return scoreNeedsRefinement;
	}


	public double getParentLBScore() {
		return parentlbScore;
	}


	public double getLBScore() {
		return lbScore;
	}


	public double getUBScore() {

		if(ub == null) return ubScore;

		PFAbstract.suppressOutput = true;
		//ub.run(wt, false, true);
		ub.runPF(ub.getPF(KSTermini.COMPLEX), null, true, false);
		PFAbstract.suppressOutput = false;

		return ubScore = -1.0 * ub.getKStarScoreLog10(true);
	}


	private boolean isFullyDefined() {
		int maxDepth = wt.getPF(KSTermini.COMPLEX).getSequence().size();

		if( depth() < maxDepth )
			return false;

		return true;
	}


	public void checkConsistency(KAStarNode node) {
		// score must be > lb
		double nodeLB = node.getLBScore(), parentLB = node.getParentLBScore();

		if(nodeLB == Double.NEGATIVE_INFINITY)
			return;

		if( parentLB > nodeLB ) 
			throw new RuntimeException("ERROR: parentLB: " + parentLB + " must be <= nodeLB: " + nodeLB);
	}


	private boolean childScoreNeedsRefinement(ConcurrentHashMap<Integer, PFAbstract> lbPFs) {

		PFAbstract pf = lbPFs.get(KSTermini.COMPLEX);

		if(!pf.isFullyDefined()) return true;

		if(!KSImplKAStar.useTightBounds)
			return false;

		return true;
	}


	public boolean isFullyProcessed() {
		return !scoreNeedsRefinement() && lb.getEpsilonStatus() != EApproxReached.FALSE;
	}


	//Comparator anonymous class implementation
	public static Comparator<KAStarNode> KUStarNodeComparator = new Comparator<KAStarNode>() {

		@Override
		public int compare( KAStarNode lhs, KAStarNode rhs ) {
			return lhs.lbScore >= rhs.lbScore ? 1 : -1;
		}
	};
}
