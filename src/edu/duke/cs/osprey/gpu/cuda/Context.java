package edu.duke.cs.osprey.gpu.cuda;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

import edu.duke.cs.osprey.tools.ResourceExtractor;
import jcuda.Pointer;
import jcuda.driver.CUcontext;
import jcuda.driver.CUctx_flags;
import jcuda.driver.CUdeviceptr;
import jcuda.driver.CUfunction;
import jcuda.driver.CUmodule;
import jcuda.driver.JCudaDriver;

public class Context {
	
	private Gpu gpu;
	private CUcontext context;
	
	private Map<String,CUmodule> kernels;
	
	public Context(Gpu gpu) {
		
		this.gpu = gpu;
		
		// create the cuda context
		context = new CUcontext();
		int flags = CUctx_flags.CU_CTX_SCHED_YIELD;
		//int flags = CUctx_flags.CU_CTX_SCHED_SPIN;
		//int flags = CUctx_flags.CU_CTX_SCHED_BLOCKING_SYNC;
		JCudaDriver.cuCtxCreate(context, flags, gpu.getDevice());
		
		kernels = new HashMap<>();
	}
	
	public Gpu getGpu() {
		return gpu;
	}
	
	public synchronized CUmodule getKernel(String name)
	throws IOException {
		
		// check the cache first
		CUmodule kernel = kernels.get(name);
		if (kernel == null) {
				
			// cache miss, load the kernel
		
			// do we have a kernel binary?
			String resourcePath = String.format("kernelBinaries/%s.bin", name);
			URL url = getClass().getResource(resourcePath);
			if (url == null) {
				throw new IOException("precompiled kernel binary not found at: " + resourcePath);
			}
			
			// yup, load it
			File kernelFile = ResourceExtractor.extract(url);
			kernel = new CUmodule();
			JCudaDriver.cuModuleLoad(kernel, kernelFile.getAbsolutePath());
			
			// update the cache
			kernels.put(name, kernel);
		}
		
		return kernel;
	}
	
	public CUdeviceptr malloc(long numBytes) {
		CUdeviceptr pdBuf = new CUdeviceptr();
		JCudaDriver.cuMemAlloc(pdBuf, numBytes);
		return pdBuf;
	}
	
	public void free(CUdeviceptr pdBuf) {
		JCudaDriver.cuMemFree(pdBuf);
	}
	
	public void uploadAsync(CUdeviceptr pdBuf, Pointer phBuf, long numBytes, GpuStream stream) {
		JCudaDriver.cuMemcpyHtoDAsync(pdBuf, phBuf, numBytes, stream.getStream());
	}
	
	public void downloadAsync(Pointer phBuf, CUdeviceptr pdBuf, long numBytes, GpuStream stream) {
		JCudaDriver.cuMemcpyDtoHAsync(phBuf, pdBuf, numBytes, stream.getStream());
	}
	
	public void pinBuffer(Pointer phBuf, long numBytes) {
		JCudaDriver.cuMemHostRegister(phBuf, numBytes, 0);
	}
	
	public void unpinBuffer(Pointer phBuf) {
		JCudaDriver.cuMemHostUnregister(phBuf);
	}
	
	public void launchKernel(CUfunction func, int gridBlocks, int blockThreads, int sharedMemBytes, Pointer pArgs, GpuStream stream) {
		JCudaDriver.cuLaunchKernel(
			func,
			gridBlocks, 1, 1,
			blockThreads, 1, 1,
			sharedMemBytes,
			stream.getStream(),
			pArgs,
			null
		);
	}
	
	public void waitForGpu() {
		JCudaDriver.cuCtxSynchronize();
	}
	
	public void attachCurrentThread() {
		JCudaDriver.cuCtxSetCurrent(context);
	}
	
	public synchronized void cleanup() {
		
		for (CUmodule kernel : kernels.values()) {
			JCudaDriver.cuModuleUnload(kernel);
		}
		kernels.clear();
		
		JCudaDriver.cuCtxDestroy(context);
	}
}
