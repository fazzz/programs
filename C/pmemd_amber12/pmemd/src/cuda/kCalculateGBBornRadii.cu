/***************************************************/
/*                                                 */
/*      AMBER NVIDIA CUDA CPU IMPLEMENTATION       */
/*                 PMEMD VERSION                   */
/*                     2010                        */
/*                      by                         */
/*             Scott Le Grand (NVIDIA)             */
/*               Duncan Poole (NVIDIA)             */
/*                Ross Walker (SDSC)               */
/*                                                 */
/***************************************************/

#include <cuda.h>
#include "gpu.h"
#include "ptxmacros.h"
static __constant__ cudaSimulation cSim;
static __constant__ PMEFloat ta                 = (PMEFloat)(1.0 / 3.0);
static __constant__ PMEFloat tb                 = (PMEFloat)(2.0 / 5.0);
static __constant__ PMEFloat tc                 = (PMEFloat)(3.0 / 7.0);
static __constant__ PMEFloat td                 = (PMEFloat)(4.0 / 9.0);
static __constant__ PMEFloat tdd                = (PMEFloat)(5.0 / 11.0);

void SetkCalculateGBBornRadiiSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetkCalculateGBBornRadiiSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

__global__ void 
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_GBBORNRADII_THREADS_PER_BLOCK, SM_3X_GBBORNRADII_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_GBBORNRADII_THREADS_PER_BLOCK, SM_2X_GBBORNRADII_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_GBBORNRADII_THREADS_PER_BLOCK, SM_13_GBBORNRADII_BLOCKS_MULTIPLIER)
#endif
kCalculateGBBornRadii_kernel()
#include "kCalculateGBBornRadii.h"

#define GB_IGB78
__global__ void 
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_GBBORNRADII_THREADS_PER_BLOCK, SM_3X_GBBORNRADII_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_GBBORNRADII_THREADS_PER_BLOCK, SM_2X_GBBORNRADII_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_GBBORNRADIIIGB78_THREADS_PER_BLOCK, SM_13_GBBORNRADII_BLOCKS_MULTIPLIER)
#endif
kCalculateGBBornRadiiIGB78_kernel()
#include "kCalculateGBBornRadii.h"
#undef GB_IGB78

void kCalculateGBBornRadiiInitKernels(gpuContext gpu)
{
    cudaFuncSetSharedMemConfig(kCalculateGBBornRadii_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateGBBornRadiiIGB78_kernel, cudaSharedMemBankSizeEightByte);
    if (gpu->sm_version >= SM_3X)
    {
        cudaFuncSetCacheConfig(kCalculateGBBornRadii_kernel, cudaFuncCachePreferL1);
        cudaFuncSetCacheConfig(kCalculateGBBornRadiiIGB78_kernel, cudaFuncCachePreferL1);
    }
}

void kCalculateGBBornRadii(gpuContext gpu)
{
    if ((gpu->sim.igb == 7) || (gpu->sim.igb == 8))
        kCalculateGBBornRadiiIGB78_kernel<<<gpu->GBBornRadiiBlocks, gpu->GBBornRadiiIGB78ThreadsPerBlock>>>();
    else
        kCalculateGBBornRadii_kernel<<<gpu->GBBornRadiiBlocks, gpu->GBBornRadiiThreadsPerBlock>>>();
    LAUNCHERROR("kCalculateGBBornRadii");
}

__global__ void 
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#endif
kReduceGBBornRadii_kernel()
{
#ifndef MPI
    bool bIGB2578                                       = (cSim.igb == 2) || (cSim.igb == 5) || (cSim.igb == 7) || (cSim.igb == 8);
#endif
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x;   
    while (pos < cSim.atoms)
    {
        PMEDouble reff_i                                = (PMEDouble)0.0;
#ifndef MPI      
        PMEDouble rborn_i                               = cSim.pAtomRBorn[pos];
#endif

#ifdef use_SPFP
        reff_i                                          = (PMEDouble)cSim.pReffAccumulator[pos] * ONEOVERFORCESCALE;
#else
        unsigned int pos1                               = pos;
        while (pos1 < cSim.stride * cSim.nonbondForceBuffers)
        {
            reff_i                                     += cSim.pReffBuffer[pos1];
            pos1                                       += cSim.stride;
        }
#endif

#ifndef MPI        
        // Process final Born Radii
        PMEDouble ri                                    = rborn_i - cSim.offset;
        PMEDouble ri1i                                  = (PMEDouble)1.0 / ri;
       
        if (bIGB2578)
        {
            // apply the new Onufriev "gbalpha, gbbeta, gbgamma" correction:
            PMEDouble psi_i                             = -ri * reff_i;
            if (cSim.igb == 8)
                reff_i                                  = ri1i - tanh((cSim.pgb_alpha[pos] + cSim.pgb_gamma[pos] * psi_i * psi_i - cSim.pgb_beta[pos] * psi_i) * psi_i) / rborn_i;
            else
                reff_i                                  = ri1i - tanh((cSim.gb_alpha + cSim.gb_gamma * psi_i * psi_i - cSim.gb_beta * psi_i) * psi_i) / rborn_i;
            reff_i                                      = max(reff_i, (PMEDouble)1.0 / (PMEDouble)30.0);
            reff_i                                      = (PMEDouble)1.0 / reff_i;
            cSim.pPsi[pos]                              = psi_i;
        }
        else
        {
            // "standard" GB, including the "diagonal" term here:
            reff_i                                      = (PMEDouble)1.0 / (reff_i + ri1i); 
        }      

        cSim.pReffSP[pos]                               = reff_i;  
        cSim.pReff[pos]                                 = reff_i;      
#else
        cSim.pReffa[pos]                                = reff_i;
#endif
        pos                                            += blockDim.x * gridDim.x;
    }

}

void kReduceGBBornRadii(gpuContext gpu)
{
    kReduceGBBornRadii_kernel<<<gpu->blocks, gpu->reduceBufferThreadsPerBlock>>>();
    LAUNCHERROR("kReduceGBBornRadii");
}

#ifdef MPI
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#endif
kProcessGBBornRadii_kernel()
{
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x;
    bool bIGB2578                                       = (cSim.igb == 2) || (cSim.igb == 5) || (cSim.igb == 7) || (cSim.igb == 8);
    while (pos < cSim.atoms)
    {
        PMEDouble reff_i                                = cSim.pReffa[pos];
        PMEDouble rborn_i                               = cSim.pAtomRBorn[pos];
    
        // Process final Born Radii
        PMEDouble ri                                    = rborn_i - cSim.offset;
        PMEDouble ri1i                                  = (PMEDouble)1.0 / ri;
       
        if (bIGB2578)
        {
            // apply the new Onufriev "gbalpha, gbbeta, gbgamma" correction:
            PMEDouble psi_i                             = -ri * reff_i;
            if (cSim.igb == 8)
                reff_i                                  = ri1i - tanh((cSim.pgb_alpha[pos] + cSim.pgb_gamma[pos] * psi_i * psi_i - cSim.pgb_beta[pos] * psi_i) * psi_i) / rborn_i;
            else
                reff_i                                  = ri1i - tanh((cSim.gb_alpha + cSim.gb_gamma * psi_i * psi_i - cSim.gb_beta * psi_i) * psi_i) / rborn_i;
            reff_i                                      = max(reff_i, (PMEDouble)1.0 / (PMEDouble)30.0);
            reff_i                                      = (PMEDouble)1.0 / reff_i;
            cSim.pPsi[pos]                              = psi_i;
        }
        else
        {
            // "standard" GB, including the "diagonal" term here:
            reff_i                                      = (PMEDouble)1.0 / (reff_i + ri1i); 
        }      
        cSim.pReffSP[pos]                               = reff_i;  
        cSim.pReff[pos]                                 = reff_i;      
        pos                                            += blockDim.x * gridDim.x;
    }

}

void kProcessGBBornRadii(gpuContext gpu)
{
    kProcessGBBornRadii_kernel<<<gpu->blocks, gpu->reduceBufferThreadsPerBlock>>>();
    LAUNCHERROR("kProcessGBBornRadii");
}
#endif

__global__ void 
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#endif
kClearGBBuffers_kernel()
{
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x;   
    while (pos < cSim.atoms)
    {

#ifdef use_SPFP
        cSim.pReffAccumulator[pos]                      = (PMEAccumulator)0;
        cSim.pSumdeijdaAccumulator[pos]                 = (PMEAccumulator)0;
#else
        unsigned int pos1                               = pos;
        while (pos1 < cSim.stride * cSim.nonbondForceBuffers)
        {
            cSim.pReffBuffer[pos1]                      = 0.0;
            cSim.pSumdeijdaBuffer[pos1]                 = 0.0;
            pos1                                       += cSim.stride;
        }
#endif

        pos                                            += blockDim.x * gridDim.x;
    }
}

void kClearGBBuffers(gpuContext gpu)
{
    kClearGBBuffers_kernel<<<gpu->blocks, gpu->reduceBufferThreadsPerBlock>>>();
    LAUNCHERROR("kClearGBBuffers");
}



