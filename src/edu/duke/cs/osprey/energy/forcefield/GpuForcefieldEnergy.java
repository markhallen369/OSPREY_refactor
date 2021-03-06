package edu.duke.cs.osprey.energy.forcefield;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.ForcefieldKernel;
import edu.duke.cs.osprey.gpu.cuda.GpuStream;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.cuda.kernels.ForcefieldKernelCuda;
import edu.duke.cs.osprey.gpu.opencl.GpuQueue;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.gpu.opencl.kernels.ForcefieldKernelOpenCL;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class GpuForcefieldEnergy implements EnergyFunction.DecomposableByDof, EnergyFunction.NeedsCleanup, EnergyFunction.ExplicitChemicalChanges {
	
	private static final long serialVersionUID = -9142317985561910731L;
	
	private class KernelBuilder {
		
		private int subsetSequenceNumber;
		private ForcefieldKernel kernel;
		
		public KernelBuilder() {
			subsetSequenceNumber = -1;
			kernel = null;
		}
		
		public ForcefieldKernel get(int expectedSequenceNumber) {
			
			// if this kernel doesn't match the sequence, wipe it
			if (kernel != null && expectedSequenceNumber != subsetSequenceNumber) {
				kernel.cleanup();
				kernel = null;
			}
			
			// make a new kernel if needed
			if (kernel == null) {
				try {
					
					if (openclQueue != null) {
						kernel = new ForcefieldKernelOpenCL(openclQueue, ffenergy);
					} else if (cudaStream != null) {
						kernel = new ForcefieldKernelCuda(cudaStream, ffenergy);
					} else {
						throw new Error("bad gpu queue/context configuration, this is a bug");
					}
					
					// update the sequence
					subsetSequenceNumber = kernel.getForcefield().getFullSubset().handleChemicalChanges();
					
				} catch (IOException ex) {
					
					// if we can't find the gpu kernel source, that's something a programmer needs to fix
					throw new Error("can't initialize gpu kernel", ex);
				}
			}
			
			return kernel;
		}
		
		public void cleanup() {
			if (kernel != null) {
				kernel.cleanup();
				kernel = null;
			}
		}
	}
	
	private BigForcefieldEnergy ffenergy;
	private BigForcefieldEnergy.Subset ffsubset;
	private GpuQueuePool openclQueuePool;
	private GpuQueue openclQueue;
	private GpuStreamPool cudaStreamPool;
	private GpuStream cudaStream;
	private KernelBuilder kernelBuilder;
	private Map<Residue,GpuForcefieldEnergy> efuncCache;
	
	public GpuForcefieldEnergy(ForcefieldParams ffparams, ForcefieldInteractions interactions, GpuQueuePool queuePool) {
		this.ffenergy = new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct);
		this.ffsubset = ffenergy.getFullSubset();
		this.openclQueuePool = queuePool;
		this.openclQueue = queuePool.checkout();
		this.cudaStreamPool = null;
		this.cudaStream = null;
		this.kernelBuilder = new KernelBuilder();
		this.efuncCache = new HashMap<>();
	}
	
	public GpuForcefieldEnergy(ForcefieldParams ffparams, ForcefieldInteractions interactions, GpuStreamPool streamPool) {
		this.ffenergy = new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct);
		this.ffsubset = ffenergy.getFullSubset();
		this.openclQueuePool = null;
		this.openclQueue = null;
		this.cudaStreamPool = streamPool;
		this.cudaStream = streamPool.checkout();
		this.kernelBuilder = new KernelBuilder();
		this.efuncCache = new HashMap<>();
	}
	
	public GpuForcefieldEnergy(GpuForcefieldEnergy parent, ForcefieldInteractions interactions) {
		this.ffenergy = parent.ffenergy;
		this.ffsubset = ffenergy.new Subset(interactions);
		if (parent.openclQueue != null) {
			this.openclQueuePool = null;
			this.openclQueue = parent.openclQueue;
			this.cudaStreamPool = null;
			this.cudaStream = null;
		} else {
			this.openclQueuePool = null;
			this.openclQueue = null;
			this.cudaStreamPool = null;
			this.cudaStream = parent.cudaStream;
		}
		this.kernelBuilder = parent.kernelBuilder;
		this.efuncCache = null;
	}
	
	public ForcefieldKernel getKernel() {
		return kernelBuilder.get(handleChemicalChanges());
	}
	
	public BigForcefieldEnergy.Subset getSubset() {
		return ffsubset;
	}
	
	@Override
	public int handleChemicalChanges() {
		return ffsubset.handleChemicalChanges();
	}
	
	@Override
	public double getEnergy() {
		
		// upload data
		ForcefieldKernel kernel = getKernel();
		kernel.setSubset(getSubset());
		kernel.uploadCoordsAsync();
		
		// compute the energies
		kernel.runAsync();
		
		// read the results
		return kernel.downloadEnergySync();
	}
	
	@Override
	public void cleanup() {
		
		kernelBuilder.cleanup();
		
		if (openclQueuePool != null) {
			openclQueuePool.release(openclQueue);
		}
		
		if (cudaStreamPool != null) {
			cudaStreamPool.release(cudaStream);
		}
		
		if (efuncCache != null) {
			for (GpuForcefieldEnergy efunc : efuncCache.values()) {
				efunc.cleanup();
			}
			efuncCache.clear();
		}
	}

	@Override
	public List<EnergyFunction> decomposeByDof(Molecule m, List<DegreeOfFreedom> dofs) {
		
		List<EnergyFunction> efuncs = new ArrayList<>();
		
		for (DegreeOfFreedom dof : dofs) {

			Residue res = dof.getResidue();
			if (res == null) {
				
				// when there's no residue at the dof, then use the whole efunc
				efuncs.add(this);
				
			} else {
				
				// otherwise, make an efunc for only that residue
				// but share efuncs between dofs in the same residue
				GpuForcefieldEnergy efunc = efuncCache.get(res);
				if (efunc == null) {
					efunc = new GpuForcefieldEnergy(this, ffenergy.getInteractions().makeSubsetByResidue(res));
					efuncCache.put(res, efunc);
				}
				efuncs.add(efunc);
			}
		}
		
		return efuncs;
	}
}
