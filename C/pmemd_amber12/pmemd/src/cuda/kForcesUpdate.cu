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

// kForces.cu
void SetkForcesUpdateSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetkForcesUpdateSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

#ifdef MPI
__global__ void 
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_CLEARFORCES_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CLEARFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CLEARFORCES_THREADS_PER_BLOCK, 1)
#endif
kClearNBForces_kernel()
{
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x;   
#pragma unroll 16    
    while (pos < cSim.stride3)
    {
        cSim.pNBForce[pos]                              = (PMEDouble)0.0;
        pos                                            += blockDim.x * gridDim.x;
    }
}   

void kClearNBForces(gpuContext gpu)
{

    kClearNBForces_kernel<<<gpu->blocks, gpu->NLClearForcesThreadsPerBlock>>>(); 
    LAUNCHERROR("kClearNBForces");
}
#endif

__global__ void 
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_CLEARFORCES_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CLEARFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CLEARFORCES_THREADS_PER_BLOCK, 1)
#endif
kClearForces_kernel()
{
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x; 

    // Clear GB NB kernel counters
    if (pos < 3)
        cSim.pGBBRPosition[pos]                         = cSim.GBTotalWarps[pos];
  
    if (pos < ENERGYTERMS)
        cSim.pEnergyBuffer[pos]                         = 0;

#ifdef use_SPFP  
    unsigned int count                                  = cSim.stride3 * cSim.nonbondForceBuffers;
#pragma unroll 16 
    while (pos < count)
    {
        cSim.pForceAccumulator[pos]                     = (PMEAccumulator)0;
        pos                                            += blockDim.x * gridDim.x;
    }
#else
    while (pos < cSim.stride3)
    {
        unsigned int pos1                               = pos;
#pragma unroll 16        
        while (pos1 < cSim.stride3 * cSim.maxForceBuffers)
        {
            cSim.pForceBuffer[pos1]                     = (PMEDouble)0.0;
            pos1                                       += cSim.stride3;
        }
        pos                                            += blockDim.x * gridDim.x;
    }
#endif
}   

#ifdef use_SPFP
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_CLEARFORCES_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CLEARFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CLEARFORCES_THREADS_PER_BLOCK, 1)
#endif
kNLClearForces_kernel()
{
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x; 
    unsigned int count                                  = cSim.stride3 * cSim.nonbondForceBuffers;  
    if (pos == 0)
        *(cSim.pNLPosition)                             = cSim.NLNonbondWarps;  

    if (pos < ENERGYTERMS)
        cSim.pEnergyBuffer[pos]                         = 0;

#pragma unroll 16   
    while (pos < count)
    {
        cSim.pForceAccumulator[pos]                     = (PMEAccumulator)0;
        pos                                            += blockDim.x * gridDim.x;
    }
}   


void kClearForces(gpuContext gpu)
{
    if (gpu->bNeighborList)
    {
        kNLClearForces_kernel<<<gpu->blocks, gpu->clearForcesThreadsPerBlock>>>();
    }
    else
    {
        kClearForces_kernel<<<gpu->blocks, gpu->clearForcesThreadsPerBlock>>>();  
    } 
    LAUNCHERROR("kClearForces");
}



__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#endif
kReduceForces_kernel()
{
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x;  
#pragma unroll 16
    while (pos < cSim.stride3)
    {
        PMEAccumulator inForce                          = cSim.pForceAccumulator[pos];
        PMEDouble outForce                              = (PMEDouble)inForce * ONEOVERFORCESCALE;
        cSim.pForce[pos]                                = outForce;
        pos                                            += blockDim.x * gridDim.x;
    }
}

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#endif
kNTPReduceForces_kernel()
{
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x;  
#pragma unroll 16
    while (pos < cSim.stride3)
    {
        PMEAccumulator inForce1                         = cSim.pBondedForceAccumulator[pos];
        PMEAccumulator inForce2                         = cSim.pNBForceAccumulator[pos];
        PMEDouble outForce                              = (PMEDouble)(inForce1 + inForce2) * ONEOVERFORCESCALE;
        cSim.pForce[pos]                                = outForce;
        PMEDouble NBForce                               = (PMEDouble)inForce2 * ONEOVERFORCESCALE;
        cSim.pNBForce[pos]                              = NBForce;
        pos                                            += blockDim.x * gridDim.x;
    }
}


void kReduceForces(gpuContext gpu)
{
    if (gpu->sim.ntp > 0)
    {
        kNTPReduceForces_kernel<<<gpu->blocks, gpu->reduceForcesThreadsPerBlock>>>();
    }
    else
    {
        kReduceForces_kernel<<<gpu->blocks, gpu->reduceForcesThreadsPerBlock>>>();  
    }
    LAUNCHERROR("kReduceForces");
}

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#endif
kReduceNBForces_kernel()
{
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x;  
    while (pos < cSim.stride3)
    {
        PMEAccumulator inForce                          = cSim.pNBForceAccumulator[pos];
        PMEDouble NBForce                               = (PMEDouble)inForce * ONEOVERFORCESCALE;
        cSim.pNBForce[pos]                              = NBForce;
        pos                                            += blockDim.x * gridDim.x;
    }
}


void kReduceNBForces(gpuContext gpu)
{
    kReduceNBForces_kernel<<<gpu->blocks, gpu->reduceForcesThreadsPerBlock>>>();     
    LAUNCHERROR("kReduceNBForces");
}

#else  // use_SPFP


__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_NLCLEARFORCES_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_NLCLEARFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_NLCLEARFORCES_THREADS_PER_BLOCK, 1)
#endif
kNLClearForces_kernel()
#include "kNLCF.h"

#define CLEAR_NTP
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_CLEARFORCES_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CLEARFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CLEARFORCES_THREADS_PER_BLOCK, 1)
#endif
kNLClearForcesNTP_kernel()
#include "kNLCF.h"

#define CLEAR_LARGE
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_CLEARFORCES_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CLEARFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CLEARFORCES_THREADS_PER_BLOCK, 1)
#endif
kNLClearForcesNTPLarge_kernel()
#include "kNLCF.h"
#undef CLEAR_LARGE
#undef CLEAR_NTP


#if (__CUDA_ARCH__ >= 300) 
#define CLEARFUNC(FUNCNAME, YDIVISOR, XDIVISOR) __global__ void __launch_bounds__(SM_3X_CLEARFORCES_THREADS_PER_BLOCK, 1) kNLClear##FUNCNAME##_##YDIVISOR##_##XDIVISOR##_kernel() 
#elif (__CUDA_ARCH__ >= 200) 
#define CLEARFUNC(FUNCNAME, YDIVISOR, XDIVISOR) __global__ void __launch_bounds__(SM_2X_CLEARFORCES_THREADS_PER_BLOCK, 1) kNLClear##FUNCNAME##_##YDIVISOR##_##XDIVISOR##_kernel() 
#else
#define CLEARFUNC(FUNCNAME, YDIVISOR, XDIVISOR) __global__ void __launch_bounds__(SM_13_CLEARFORCES_THREADS_PER_BLOCK, 1) kNLClear##FUNCNAME##_##YDIVISOR##_##XDIVISOR##_kernel() 
#endif 

#define CLEAR_YDIVISOR 1
#define CLEAR_XDIVISOR 1
CLEARFUNC(Forces, 1, 1)
#include "kNLCF.h"

#define CLEAR_NTP
CLEARFUNC(ForcesNTP, 1, 1)
#include "kNLCF.h"

#define CLEAR_LARGE
CLEARFUNC(ForcesNTPLarge, 1, 1)
#include "kNLCF.h"
#undef CLEAR_LARGE
#undef CLEAR_NTP
#undef CLEAR_YDIVISOR
#undef CLEAR_XDIVISOR

#define CLEAR_YDIVISOR 2
#define CLEAR_XDIVISOR 1
CLEARFUNC(Forces, 2, 1)
#include "kNLCF.h"

#define CLEAR_NTP
CLEARFUNC(ForcesNTP, 2, 1)
#include "kNLCF.h"

#define CLEAR_LARGE
CLEARFUNC(ForcesNTPLarge, 2, 1)
#include "kNLCF.h"
#undef CLEAR_LARGE
#undef CLEAR_NTP
#undef CLEAR_YDIVISOR
#undef CLEAR_XDIVISOR

#define CLEAR_YDIVISOR 3
#define CLEAR_XDIVISOR 1
CLEARFUNC(Forces, 3, 1)
#include "kNLCF.h"

#define CLEAR_NTP
CLEARFUNC(ForcesNTP, 3, 1)
#include "kNLCF.h"

#define CLEAR_LARGE
CLEARFUNC(ForcesNTPLarge, 3, 1)
#include "kNLCF.h"
#undef CLEAR_LARGE
#undef CLEAR_NTP
#undef CLEAR_YDIVISOR
#undef CLEAR_XDIVISOR

#define CLEAR_YDIVISOR 3
#define CLEAR_XDIVISOR 2
CLEARFUNC(Forces, 3, 2)
#include "kNLCF.h"

#define CLEAR_NTP
CLEARFUNC(ForcesNTP, 3, 2)
#include "kNLCF.h"

#define CLEAR_LARGE
CLEARFUNC(ForcesNTPLarge, 3, 2)
#include "kNLCF.h"
#undef CLEAR_LARGE
#undef CLEAR_NTP
#undef CLEAR_YDIVISOR
#undef CLEAR_XDIVISOR

#define CLEAR_YDIVISOR 3
#define CLEAR_XDIVISOR 5
CLEARFUNC(Forces, 3, 5)
#include "kNLCF.h"

#define CLEAR_NTP
CLEARFUNC(ForcesNTP, 3, 5)
#include "kNLCF.h"

#define CLEAR_LARGE
CLEARFUNC(ForcesNTPLarge, 3, 5)
#include "kNLCF.h"
#undef CLEAR_LARGE
#undef CLEAR_NTP
#undef CLEAR_YDIVISOR
#undef CLEAR_XDIVISOR

#define CLEAR_YDIVISOR 4
#define CLEAR_XDIVISOR 5
CLEARFUNC(Forces, 4, 5)
#include "kNLCF.h"

#define CLEAR_NTP
CLEARFUNC(ForcesNTP, 4, 5)
#include "kNLCF.h"

#define CLEAR_LARGE
CLEARFUNC(ForcesNTPLarge, 4, 5)
#include "kNLCF.h"
#undef CLEAR_LARGE
#undef CLEAR_NTP
#undef CLEAR_YDIVISOR
#undef CLEAR_XDIVISOR

#define CLEAR_YDIVISOR 4
#define CLEAR_XDIVISOR 7
CLEARFUNC(Forces, 4, 7)
#include "kNLCF.h"

#define CLEAR_NTP
CLEARFUNC(ForcesNTP, 4, 7)
#include "kNLCF.h"

#define CLEAR_LARGE
CLEARFUNC(ForcesNTPLarge, 4, 7)
#include "kNLCF.h"
#undef CLEAR_LARGE
#undef CLEAR_NTP
#undef CLEAR_YDIVISOR
#undef CLEAR_XDIVISOR

#define CLEAR_YDIVISOR 4
#define CLEAR_XDIVISOR 14
CLEARFUNC(Forces, 4, 14)
#include "kNLCF.h"

#define CLEAR_NTP
CLEARFUNC(ForcesNTP, 4, 14)
#include "kNLCF.h"

#define CLEAR_LARGE
CLEARFUNC(ForcesNTPLarge, 4, 14)
#include "kNLCF.h"
#undef CLEAR_LARGE
#undef CLEAR_NTP
#undef CLEAR_YDIVISOR
#undef CLEAR_XDIVISOR

#define CLEAR_YDIVISOR 8
#define CLEAR_XDIVISOR 14
CLEARFUNC(Forces, 8, 14)
#include "kNLCF.h"

#define CLEAR_NTP
CLEARFUNC(ForcesNTP, 8, 14)
#include "kNLCF.h"

#define CLEAR_LARGE
CLEARFUNC(ForcesNTPLarge, 8, 14)
#include "kNLCF.h"
#undef CLEAR_LARGE
#undef CLEAR_NTP
#undef CLEAR_YDIVISOR
#undef CLEAR_XDIVISOR

typedef void (*KernelPointer)();
static KernelPointer spNLReduceForcesKernel = NULL;
static KernelPointer spNLClearForcesKernel = NULL;

void SetNLClearForcesKernel(gpuContext gpu)
{
    unsigned int clearType                          = gpu->sim.NLYDivisor * 100 + gpu->sim.NLXDivisor;
    if (gpu->sim.ntp > 0)
    {
        if (gpu->sim.NLCellBuffers >= gpu->sim.maxForceBuffers)
        {
            switch (clearType)
            {
                case 101:
                    spNLClearForcesKernel                   = &kNLClearForcesNTP_1_1_kernel;
                    break;
                        
                case 201:
                    spNLClearForcesKernel                   = &kNLClearForcesNTP_2_1_kernel;
                    break;
                        
                case 301:
                    spNLClearForcesKernel                   = &kNLClearForcesNTP_3_1_kernel;
                    break;
                        
                case 302:
                    spNLClearForcesKernel                   = &kNLClearForcesNTP_3_2_kernel;
                    break;
                case 305:
                    spNLClearForcesKernel                   = &kNLClearForcesNTP_3_5_kernel;
                    break;
                        
                case 405:
                    spNLClearForcesKernel                   = &kNLClearForcesNTP_4_5_kernel;
                    break;
                        
                case 407:
                    spNLClearForcesKernel                   = &kNLClearForcesNTP_4_7_kernel;
                    break;
                        
                case 414:
                    spNLClearForcesKernel                   = &kNLClearForcesNTP_4_14_kernel;
                    break;  
                                      
                case 814:
                    spNLClearForcesKernel                   = &kNLClearForcesNTP_8_14_kernel;
                    break;                                
                }                       
        }
        else
        {           
                switch (clearType)
                {
                case 101:
                    spNLClearForcesKernel                   = &kNLClearForcesNTPLarge_1_1_kernel;
                    break;
                        
                case 201:
                    spNLClearForcesKernel                   = &kNLClearForcesNTPLarge_2_1_kernel;
                    break;
                        
                case 301:
                    spNLClearForcesKernel                   = &kNLClearForcesNTPLarge_3_1_kernel;
                    break;
                        
                case 302:
                    spNLClearForcesKernel                   = &kNLClearForcesNTPLarge_3_2_kernel;
                        break;
                case 305:
                    spNLClearForcesKernel                   = &kNLClearForcesNTPLarge_3_5_kernel;
                    break;
                        
                case 405:
                    spNLClearForcesKernel                   = &kNLClearForcesNTPLarge_4_5_kernel;
                        break;
                        
                case 407:
                    spNLClearForcesKernel                   = &kNLClearForcesNTPLarge_4_7_kernel;
                    break;
                        
                case 414:
                    spNLClearForcesKernel                   = &kNLClearForcesNTPLarge_4_14_kernel;
                    break;  
                                     
                case 814:
                    spNLClearForcesKernel                   = &kNLClearForcesNTPLarge_8_14_kernel;
                    break;                                
            }                               
        }
    
    }
    else
    {
        switch (clearType)
        {
            case 101:
                spNLClearForcesKernel                       = &kNLClearForces_1_1_kernel;
                break;
                
            case 201:
                spNLClearForcesKernel                       = &kNLClearForces_2_1_kernel;
                break;
                
            case 301:
                spNLClearForcesKernel                       = &kNLClearForces_3_1_kernel;
                break;
                
            case 302:
                spNLClearForcesKernel                       = &kNLClearForces_3_2_kernel;
                break;
            case 305:
                spNLClearForcesKernel                       = &kNLClearForces_3_5_kernel;
                break;
                
            case 405:
                spNLClearForcesKernel                       = &kNLClearForces_4_5_kernel;
                break;
                
            case 407:
                spNLClearForcesKernel                       = &kNLClearForces_4_7_kernel;
                break;
                
            case 414:
                spNLClearForcesKernel                       = &kNLClearForces_4_14_kernel;
                break;  
                              
            case 814:
                spNLClearForcesKernel                       = &kNLClearForces_8_14_kernel;
                break;                                
        }
        
    }
    
}



void kClearForces(gpuContext gpu)
{
    if (gpu->bNeighborList)
    {
        spNLClearForcesKernel<<<gpu->blocks, gpu->NLClearForcesThreadsPerBlock>>>();
    }
    else
    {
        kClearForces_kernel<<<gpu->blocks, gpu->clearForcesThreadsPerBlock>>>();  
    }  
    LAUNCHERROR("kClearForces");
}

#if (__CUDA_ARCH__ >= 300) 
#define REDUCEFUNC(FUNCNAME, YDIVISOR, XDIVISOR) __global__ void __launch_bounds__(SM_3X_REDUCEFORCES_THREADS_PER_BLOCK, 1) kNLReduce##FUNCNAME##_##YDIVISOR##_##XDIVISOR##_kernel() 
#elif (__CUDA_ARCH__ >= 200) 
#define REDUCEFUNC(FUNCNAME, YDIVISOR, XDIVISOR) __global__ void __launch_bounds__(SM_2X_REDUCEFORCES_THREADS_PER_BLOCK, 1) kNLReduce##FUNCNAME##_##YDIVISOR##_##XDIVISOR##_kernel() 
#else
#define REDUCEFUNC(FUNCNAME, YDIVISOR, XDIVISOR) __global__ void __launch_bounds__(SM_13_REDUCEFORCES_THREADS_PER_BLOCK, 1) kNLReduce##FUNCNAME##_##YDIVISOR##_##XDIVISOR##_kernel() 
#endif 

#define REDUCE_YDIVISOR 1
#define REDUCE_XDIVISOR 1
REDUCEFUNC(Forces, 1, 1)
#include "kNLRF.h"

#define REDUCE_NTP
REDUCEFUNC(ForcesNTP, 1, 1)
#include "kNLRF.h"

#define REDUCE_LARGE
REDUCEFUNC(ForcesNTPLarge, 1, 1)
#include "kNLRF.h"

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPLargeNode0, 1, 1)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_LARGE

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPNode0, 1, 1)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_NTP
#undef REDUCE_YDIVISOR
#undef REDUCE_XDIVISOR




#define REDUCE_YDIVISOR 2
#define REDUCE_XDIVISOR 1
REDUCEFUNC(Forces, 2, 1)
#include "kNLRF.h"

#define REDUCE_NTP
REDUCEFUNC(ForcesNTP, 2, 1)
#include "kNLRF.h"

#define REDUCE_LARGE
REDUCEFUNC(ForcesNTPLarge, 2, 1)
#include "kNLRF.h"

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPLargeNode0, 2, 1)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_LARGE

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPNode0, 2, 1)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_NTP
#undef REDUCE_YDIVISOR
#undef REDUCE_XDIVISOR




#define REDUCE_YDIVISOR 3
#define REDUCE_XDIVISOR 1
REDUCEFUNC(Forces, 3, 1)
#include "kNLRF.h"

#define REDUCE_NTP
REDUCEFUNC(ForcesNTP, 3, 1)
#include "kNLRF.h"

#define REDUCE_LARGE
REDUCEFUNC(ForcesNTPLarge, 3, 1)
#include "kNLRF.h"

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPLargeNode0, 3, 1)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_LARGE

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPNode0, 3, 1)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_NTP
#undef REDUCE_YDIVISOR
#undef REDUCE_XDIVISOR



#define REDUCE_YDIVISOR 3
#define REDUCE_XDIVISOR 2
REDUCEFUNC(Forces, 3, 2)
#include "kNLRF.h"

#define REDUCE_NTP
REDUCEFUNC(ForcesNTP, 3, 2)
#include "kNLRF.h"

#define REDUCE_LARGE
REDUCEFUNC(ForcesNTPLarge, 3, 2)
#include "kNLRF.h"

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPLargeNode0, 3, 2)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_LARGE

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPNode0, 3, 2)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_NTP
#undef REDUCE_YDIVISOR
#undef REDUCE_XDIVISOR




#define REDUCE_YDIVISOR 3
#define REDUCE_XDIVISOR 5
REDUCEFUNC(Forces, 3, 5)
#include "kNLRF.h"

#define REDUCE_NTP
REDUCEFUNC(ForcesNTP, 3, 5)
#include "kNLRF.h"

#define REDUCE_LARGE
REDUCEFUNC(ForcesNTPLarge, 3, 5)
#include "kNLRF.h"

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPLargeNode0, 3, 5)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_LARGE

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPNode0, 3, 5)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_NTP
#undef REDUCE_YDIVISOR
#undef REDUCE_XDIVISOR


#define REDUCE_YDIVISOR 4
#define REDUCE_XDIVISOR 5
REDUCEFUNC(Forces, 4, 5)
#include "kNLRF.h"

#define REDUCE_NTP
REDUCEFUNC(ForcesNTP, 4, 5)
#include "kNLRF.h"

#define REDUCE_LARGE
REDUCEFUNC(ForcesNTPLarge, 4, 5)
#include "kNLRF.h"

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPLargeNode0, 4, 5)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_LARGE

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPNode0, 4, 5)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_NTP
#undef REDUCE_YDIVISOR
#undef REDUCE_XDIVISOR


#define REDUCE_YDIVISOR 4
#define REDUCE_XDIVISOR 7
REDUCEFUNC(Forces, 4, 7)
#include "kNLRF.h"

#define REDUCE_NTP
REDUCEFUNC(ForcesNTP, 4, 7)
#include "kNLRF.h"

#define REDUCE_LARGE
REDUCEFUNC(ForcesNTPLarge, 4, 7)
#include "kNLRF.h"

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPLargeNode0, 4, 7)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_LARGE

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPNode0, 4, 7)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_NTP
#undef REDUCE_YDIVISOR
#undef REDUCE_XDIVISOR


#define REDUCE_YDIVISOR 4
#define REDUCE_XDIVISOR 14
REDUCEFUNC(Forces, 4, 14)
#include "kNLRF.h"

#define REDUCE_NTP
REDUCEFUNC(ForcesNTP, 4, 14)
#include "kNLRF.h"

#define REDUCE_LARGE
REDUCEFUNC(ForcesNTPLarge, 4, 14)
#include "kNLRF.h"

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPLargeNode0, 4, 14)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_LARGE

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPNode0, 4, 14)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_NTP
#undef REDUCE_YDIVISOR
#undef REDUCE_XDIVISOR

#define REDUCE_YDIVISOR 8
#define REDUCE_XDIVISOR 14
REDUCEFUNC(Forces, 8, 14)
#include "kNLRF.h"

#define REDUCE_NTP
REDUCEFUNC(ForcesNTP, 8, 14)
#include "kNLRF.h"

#define REDUCE_LARGE
REDUCEFUNC(ForcesNTPLarge, 8, 14)
#include "kNLRF.h"

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPLargeNode0, 8, 14)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_LARGE

#ifdef MPI
#define MPI_REDUCE_NODE0
REDUCEFUNC(ForcesNTPNode0, 8, 14)
#include "kNLRF.h"
#undef MPI_REDUCE_NODE0
#endif
#undef REDUCE_NTP
#undef REDUCE_YDIVISOR
#undef REDUCE_XDIVISOR

__global__ void kReduceForces_kernel()
{
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x;      
    while (pos < cSim.stride3)
    {
        unsigned int pos1                               = pos;
        PMEDouble force                                 = (PMEDouble)0.0;
        while (pos1 < cSim.stride3 * cSim.maxForceBuffers)
        {
            force                                      += cSim.pForceBuffer[pos1];
            pos1                                       += cSim.stride3;
        }
        cSim.pForce[pos]                                = force;
        pos                                            += blockDim.x * gridDim.x;
    }
}

void SetNLReduceForcesKernel(gpuContext gpu)
{
    unsigned int reductionType                          = gpu->sim.NLYDivisor * 100 + gpu->sim.NLXDivisor;
    if (gpu->sim.ntp > 0)
    {
        if (gpu->sim.NLCellBuffers >= gpu->sim.maxForceBuffers)
        {
#ifdef MPI
            if ((gpu->gpuID == 0) & !gpu->sim.bIPSActive)
            {
                switch (reductionType)
                {
                    case 101:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPNode0_1_1_kernel;
                        break;
                        
                    case 201:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPNode0_2_1_kernel;
                        break;
                        
                    case 301:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPNode0_3_1_kernel;
                        break;
                        
                    case 302:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPNode0_3_2_kernel;
                        break;
                    case 305:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPNode0_3_5_kernel;
                        break;
                        
                    case 405:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPNode0_4_5_kernel;
                        break;
                        
                    case 407:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPNode0_4_7_kernel;
                        break;
                        
                    case 414:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPNode0_4_14_kernel;
                        break;  
                                      
                    case 814:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPNode0_8_14_kernel;
                        break;                                
                }
            }
            else
#endif            
            {
                switch (reductionType)
                {
                    case 101:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTP_1_1_kernel;
                        break;
                        
                    case 201:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTP_2_1_kernel;
                        break;
                        
                    case 301:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTP_3_1_kernel;
                        break;
                        
                    case 302:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTP_3_2_kernel;
                        break;
                    case 305:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTP_3_5_kernel;
                        break;
                        
                    case 405:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTP_4_5_kernel;
                        break;
                        
                    case 407:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTP_4_7_kernel;
                        break;
                        
                    case 414:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTP_4_14_kernel;
                        break;  
                                      
                    case 814:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTP_8_14_kernel;
                        break;                                
                }                       
            }
        }
        else
        {
#ifdef MPI
            if ((gpu->gpuID == 0) & !gpu->sim.bIPSActive)
            {
                switch (reductionType)
                {
                    case 101:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLargeNode0_1_1_kernel;
                        break;
                        
                    case 201:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLargeNode0_2_1_kernel;
                        break;
                        
                    case 301:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLargeNode0_3_1_kernel;
                        break;
                        
                    case 302:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLargeNode0_3_2_kernel;
                        break;
                    case 305:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLargeNode0_3_5_kernel;
                        break;
                        
                    case 405:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLargeNode0_4_5_kernel;
                        break;
                        
                    case 407:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLargeNode0_4_7_kernel;
                        break;
                        
                    case 414:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLargeNode0_4_14_kernel;
                        break;  
                                      
                    case 814:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLargeNode0_8_14_kernel;
                        break;                                
                }
            }
            else
#endif            
            {
                switch (reductionType)
                {
                    case 101:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLarge_1_1_kernel;
                        break;
                        
                    case 201:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLarge_2_1_kernel;
                        break;
                        
                    case 301:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLarge_3_1_kernel;
                        break;
                        
                    case 302:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLarge_3_2_kernel;
                        break;
                    case 305:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLarge_3_5_kernel;
                        break;
                        
                    case 405:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLarge_4_5_kernel;
                        break;
                        
                    case 407:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLarge_4_7_kernel;
                        break;
                        
                    case 414:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLarge_4_14_kernel;
                        break;  
                                      
                    case 814:
                        spNLReduceForcesKernel              = &kNLReduceForcesNTPLarge_8_14_kernel;
                        break;                                
                }                       
            }        
        
        }
    
    }
    else
    {
        switch (reductionType)
        {
            case 101:
                spNLReduceForcesKernel                      = &kNLReduceForces_1_1_kernel;
                break;
                
            case 201:
                spNLReduceForcesKernel                      = &kNLReduceForces_2_1_kernel;
                break;
                
            case 301:
                spNLReduceForcesKernel                      = &kNLReduceForces_3_1_kernel;
                break;
                
            case 302:
                spNLReduceForcesKernel                      = &kNLReduceForces_3_2_kernel;
                break;
            case 305:
                spNLReduceForcesKernel                      = &kNLReduceForces_3_5_kernel;
                break;
                
            case 405:
                spNLReduceForcesKernel                      = &kNLReduceForces_4_5_kernel;
                break;
                
            case 407:
                spNLReduceForcesKernel                      = &kNLReduceForces_4_7_kernel;
                break;
                
            case 414:
                spNLReduceForcesKernel                      = &kNLReduceForces_4_14_kernel;
                break;  
                              
            case 814:
                spNLReduceForcesKernel                      = &kNLReduceForces_8_14_kernel;
                break;                                
        }
        
    }
    
}


void kReduceForces(gpuContext gpu)
{
    

    if (gpu->bNeighborList)
    {
        spNLReduceForcesKernel<<<gpu->blocks, gpu->NLReduceForcesThreadsPerBlock>>>();
    }
    else
    {
        kReduceForces_kernel<<<gpu->blocks, gpu->reduceForcesThreadsPerBlock>>>();
    }
    LAUNCHERROR("kReduceForces"); 
}
#endif // use_SPFP

#ifdef MPI
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#endif
__global__ void kTransposeForces_kernel()
{
#if (__CUDA_ARCH__ >= 300)
    __shared__ volatile PMEDouble sForce[3 * SM_3X_REDUCEFORCES_THREADS_PER_BLOCK];
#elif (__CUDA_ARCH__ >= 200)
    __shared__ volatile PMEDouble sForce[3 * SM_2X_REDUCEFORCES_THREADS_PER_BLOCK];
#else
    __shared__ volatile PMEDouble sForce[3 * SM_13_REDUCEFORCES_THREADS_PER_BLOCK];
#endif 
    int opos                                            = blockIdx.x * blockDim.x + threadIdx.x;
    while (opos < cSim.reducedAtoms)
    {
        int rpos                                        = opos + cSim.minReducedAtom;
        rpos                                            = (rpos < cSim.paddedNumberOfAtoms) ? rpos : rpos - cSim.paddedNumberOfAtoms;
        PMEDouble forceX                                = cSim.pForceX[rpos];
        PMEDouble forceY                                = cSim.pForceY[rpos];
        PMEDouble forceZ                                = cSim.pForceZ[rpos];
        volatile PMEDouble* psForce                     = &sForce[3 * threadIdx.x];
        *psForce++                                      = forceX;
        *psForce++                                      = forceY;
        *psForce                                        = forceZ;        
        int tgx                                         = threadIdx.x & GRIDBITSMASK;
        int tbx                                         = threadIdx.x - tgx;
        psForce                                         = &sForce[3 * tbx + tgx];     
        PMEDouble* pOutForce                            = &cSim.pOutForce[3 * (opos - tgx) + tgx];      
        *pOutForce                                      = *psForce;
        pOutForce                                      += GRID;
        psForce                                        += GRID;
        *pOutForce                                      = *psForce;
        pOutForce                                      += GRID;
        psForce                                        += GRID;
        *pOutForce                                      = *psForce;
        opos                                           += blockDim.x * gridDim.x;
    }
}

void kTransposeForces(gpuContext gpu)
{
    //printf("%06d %06d %06d %06d\n", gpu->gpuID, gpu->sim.minReducedAtom, gpu->sim.maxReducedAtom, gpu->sim.reducedAtoms);
    kTransposeForces_kernel<<<gpu->blocks, gpu->reduceForcesThreadsPerBlock>>>();
    LAUNCHERROR("kTransposeForces");  
}
#endif


// KUpdate.cu
static __constant__ PMEDouble boltz2                    = 0.00831441 * 0.5 / 4.184;

struct COM {
    PMEFloat xmin;
    PMEFloat ymin;
    PMEFloat zmin;
    PMEFloat xmax;
    PMEFloat ymax;
    PMEFloat zmax;
};

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_UPDATE_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kRecenter_Molecule1_kernel()
{
#if (__CUDA_ARCH__ >= 300)
    __shared__ COM sA[SM_3X_UPDATE_THREADS_PER_BLOCK];
#elif (__CUDA_ARCH__ >= 200)
    __shared__ COM sA[SM_2X_UPDATE_THREADS_PER_BLOCK];
#else
    __shared__ COM sA[SM_13_UPDATE_THREADS_PER_BLOCK];
#endif 
    PMEFloat xmin                                       = (PMEFloat)999999999999.0;
    PMEFloat ymin                                       = (PMEFloat)999999999999.0;
    PMEFloat zmin                                       = (PMEFloat)999999999999.0;
    PMEFloat xmax                                       = (PMEFloat)-999999999999.0;
    PMEFloat ymax                                       = (PMEFloat)-999999999999.0;
    PMEFloat zmax                                       = (PMEFloat)-999999999999.0;
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x;
    
    
    // Perform individual sums
    while (pos < cSim.atoms)
    {
        PMEFloat2 xy                                    = cSim.pAtomXYSP[pos];
        PMEFloat z                                      = cSim.pAtomZSP[pos]; 
        xmax                                            = max(xy.x, xmax);
        xmin                                            = min(xy.x, xmin);
        ymax                                            = max(xy.y, ymax);
        ymin                                            = min(xy.y, ymin);
        zmax                                            = max(z,    zmax);
        zmin                                            = min(z,    zmin);
        pos                                            += blockDim.x * gridDim.x;
    }
    
    // Perform local reduction to thread 0
    sA[threadIdx.x].xmin                                = xmin;
    sA[threadIdx.x].ymin                                = ymin;
    sA[threadIdx.x].zmin                                = zmin;
    sA[threadIdx.x].xmax                                = xmax;
    sA[threadIdx.x].ymax                                = ymax;
    sA[threadIdx.x].zmax                                = zmax;
    __syncthreads();
    unsigned int m                                      = 1;
    while (m < blockDim.x)
    {
        int p                                           = threadIdx.x + m;
        PMEFloat xmin                                   = ((p < blockDim.x) ? sA[p].xmin : (PMEFloat)9999999999.0);
        PMEFloat ymin                                   = ((p < blockDim.x) ? sA[p].ymin : (PMEFloat)9999999999.0);
        PMEFloat zmin                                   = ((p < blockDim.x) ? sA[p].zmin : (PMEFloat)9999999999.0);
        PMEFloat xmax                                   = ((p < blockDim.x) ? sA[p].xmax : (PMEFloat)-9999999999.0);
        PMEFloat ymax                                   = ((p < blockDim.x) ? sA[p].ymax : (PMEFloat)-9999999999.0);
        PMEFloat zmax                                   = ((p < blockDim.x) ? sA[p].zmax : (PMEFloat)-9999999999.0);
        __syncthreads();
        sA[threadIdx.x].xmin                            = min(sA[threadIdx.x].xmin, xmin);
        sA[threadIdx.x].ymin                            = min(sA[threadIdx.x].ymin, ymin);
        sA[threadIdx.x].zmin                            = min(sA[threadIdx.x].zmin, zmin);
        sA[threadIdx.x].xmax                            = max(sA[threadIdx.x].xmax, xmax);
        sA[threadIdx.x].ymax                            = max(sA[threadIdx.x].ymax, ymax);
        sA[threadIdx.x].zmax                            = max(sA[threadIdx.x].zmax, zmax);
        __syncthreads();
        m                                              *= 2;
    }
    
    // Output sum if thread 0
    if (threadIdx.x == 0)
    {
        cSim.pXMin[blockIdx.x]                          = sA[threadIdx.x].xmin;
        cSim.pYMin[blockIdx.x]                          = sA[threadIdx.x].ymin;
        cSim.pZMin[blockIdx.x]                          = sA[threadIdx.x].zmin; 
        cSim.pXMax[blockIdx.x]                          = sA[threadIdx.x].xmax;
        cSim.pYMax[blockIdx.x]                          = sA[threadIdx.x].ymax;
        cSim.pZMax[blockIdx.x]                          = sA[threadIdx.x].zmax; 
    }
    
 
}

__global__ void 
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_UPDATE_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kRecenter_Molecule2_kernel()
{
#if (__CUDA_ARCH__ >= 300)
    __shared__ COM sA[SM_3X_UPDATE_THREADS_PER_BLOCK];
#elif (__CUDA_ARCH__ >= 200)
    __shared__ COM sA[SM_2X_UPDATE_THREADS_PER_BLOCK];
#else
    __shared__ COM sA[SM_13_UPDATE_THREADS_PER_BLOCK];
#endif 
    // Read in local offsets
    unsigned int pos                                    = threadIdx.x;
    while (pos < gridDim.x)
    {
        sA[pos].xmin                                    = cSim.pXMin[pos]; 
        sA[pos].ymin                                    = cSim.pYMin[pos]; 
        sA[pos].zmin                                    = cSim.pZMin[pos]; 
        sA[pos].xmax                                    = cSim.pXMax[pos]; 
        sA[pos].ymax                                    = cSim.pYMax[pos]; 
        sA[pos].zmax                                    = cSim.pZMax[pos]; 
        pos                                            += blockDim.x;
    }
    __syncthreads();
    
    // Perform local reduction to thread 0
    unsigned int m                                      = 1;
    while (m < gridDim.x)
    {
        int p                                           = threadIdx.x + m;
        PMEFloat xmin                                   = ((p < gridDim.x) ? sA[p].xmin : (PMEFloat)9999999999.0);
        PMEFloat ymin                                   = ((p < gridDim.x) ? sA[p].ymin : (PMEFloat)9999999999.0);
        PMEFloat zmin                                   = ((p < gridDim.x) ? sA[p].zmin : (PMEFloat)9999999999.0);
        PMEFloat xmax                                   = ((p < gridDim.x) ? sA[p].xmax : (PMEFloat)-9999999999.0);
        PMEFloat ymax                                   = ((p < gridDim.x) ? sA[p].ymax : (PMEFloat)-9999999999.0);
        PMEFloat zmax                                   = ((p < gridDim.x) ? sA[p].zmax : (PMEFloat)-9999999999.0);
        __syncthreads();
        sA[threadIdx.x].xmin                            = min(sA[threadIdx.x].xmin, xmin);
        sA[threadIdx.x].ymin                            = min(sA[threadIdx.x].ymin, ymin);
        sA[threadIdx.x].zmin                            = min(sA[threadIdx.x].zmin, zmin);
        sA[threadIdx.x].xmax                            = max(sA[threadIdx.x].xmax, xmax);
        sA[threadIdx.x].ymax                            = max(sA[threadIdx.x].ymax, ymax);
        sA[threadIdx.x].zmax                            = max(sA[threadIdx.x].zmax, zmax);
        __syncthreads();
        m                                              *= 2;
    }
    PMEDouble xcenter                                   = (PMEFloat)-0.5 * (sA[0].xmin + sA[0].xmax);
    PMEDouble ycenter                                   = (PMEFloat)-0.5 * (sA[0].ymin + sA[0].ymax);
    PMEDouble zcenter                                   = (PMEFloat)-0.5 * (sA[0].zmin + sA[0].zmax);
      
    // Perform individual sums
    pos                                                 = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.atoms)
    {
        PMEDouble x                                     = cSim.pAtomX[pos];
        PMEDouble y                                     = cSim.pAtomY[pos];
        PMEDouble z                                     = cSim.pAtomZ[pos];
        x                                              += xcenter;
        y                                              += ycenter; 
        z                                              += zcenter; 
        
        PMEFloat2 xy                                    = {x, y};
        cSim.pAtomX[pos]                                = x;
        cSim.pAtomY[pos]                                = y;
        cSim.pAtomZ[pos]                                = z;
        cSim.pAtomXYSP[pos]                             = xy;
        cSim.pAtomZSP[pos]                              = z;
        pos                                            += blockDim.x * gridDim.x;
    }  
    
    // Fix restraints
    pos                                                 = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.constraints)
    {
        PMEDouble2 constraint1                          = cSim.pConstraint1[pos];
        PMEDouble2 constraint2                          = cSim.pConstraint2[pos];
        constraint1.y                                  += xcenter;
        constraint2.x                                  += ycenter;
        constraint2.y                                  += zcenter;
        cSim.pConstraint1[pos]                          = constraint1;
        cSim.pConstraint2[pos]                          = constraint2;
        pos                                            += blockDim.x * gridDim.x;
    }
    
}


extern "C" void kRecenter_Molecule(gpuContext gpu)
{
    kRecenter_Molecule1_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>();
    LAUNCHERROR("kRecenter_Molecule");
    kRecenter_Molecule2_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>();
    LAUNCHERROR("kRecenter_Molecule");
}

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_UPDATE_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kUpdate_kernel(PMEDouble dt)
#include "kU.h"

#define LANGEVIN
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_UPDATE_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kLangevinUpdate_kernel(PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln)
#include "kU.h"
#undef LANGEVIN

#define PME
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_UPDATE_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMEUpdate_kernel(PMEDouble dt)
#include "kU.h"

#define LANGEVIN
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_UPDATE_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMELangevinUpdate_kernel(PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln)
#include "kU.h"
#undef LANGEVIN
#undef PME



void kUpdate(gpuContext gpu, PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln)
{
    // Choose Langevin update if necessary
    if (gpu->ntt == 3)
    {
        // Update random numbers if necessary
        if (gpu->randomCounter >= gpu->sim.randomSteps)
        {
#ifdef CPU_RANDOMS
            cpu_kRandom(gpu);
#else
            kRandom(gpu);
#endif
            gpu->randomCounter = 0;
        }
        
        if (gpu->bNeighborList)
            kPMELangevinUpdate_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(dt, temp0, gamma_ln);
        else
            kLangevinUpdate_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(dt, temp0, gamma_ln);
            
        gpu->randomCounter++; 
    }
    else
    {
        if (gpu->bNeighborList)
            kPMEUpdate_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(dt);
        else
            kUpdate_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(dt);
    }   
    LAUNCHERROR("kUpdate"); 
}


__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_UPDATE_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kResetVelocities_kernel(PMEDouble temp, PMEDouble half_dtx)
{
    unsigned int pos                        = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int increment                  = gridDim.x * blockDim.x;
    unsigned int rpos                       = cSim.pRandomPos[blockIdx.x];
    PMEDouble boltz                         = 8.31441e-3 * temp / 4.184;
    while (pos < cSim.atoms)
    {  
        PMEDouble invMass                   = cSim.pAtomInvMass[pos];    
        PMEDouble forceX                    = cSim.pForceX[pos];
        PMEDouble forceY                    = cSim.pForceY[pos];
        PMEDouble forceZ                    = cSim.pForceZ[pos];
        PMEDouble velX, velY, velZ;
        
        // Zero velocities if it's really cold
        if (temp < 1.0e-6)
        {
            velX                            = (PMEDouble)0.0;
            velY                            = (PMEDouble)0.0;
            velZ                            = (PMEDouble)0.0;
        }
        else
        {
            PMEDouble gaussX                = cSim.pRandomX[rpos + pos];
            PMEDouble gaussY                = cSim.pRandomY[rpos + pos];
            PMEDouble gaussZ                = cSim.pRandomZ[rpos + pos]; 
            PMEDouble sd                    = sqrt(boltz * invMass);        
            velX                            = sd * gaussX;
            velY                            = sd * gaussY;
            velZ                            = sd * gaussZ;        
        }
             
        // Back velocities up a half-step
        PMEDouble wfac                      = invMass * half_dtx;
        velX                               -= forceX * wfac;
        velY                               -= forceY * wfac;
        velZ                               -= forceZ * wfac;
       
        // Write final velocities
        cSim.pVelX[pos]                     = velX;
        cSim.pVelY[pos]                     = velY;
        cSim.pVelZ[pos]                     = velZ;
        pos                                += increment;
    }
    
    // Update RNG position
    __syncthreads();
    if (threadIdx.x == 0)
        cSim.pRandomPos[blockIdx.x]         = rpos + cSim.paddedNumberOfAtoms;    
}

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_UPDATE_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kNLResetVelocities_kernel(PMEDouble temp, PMEDouble half_dtx)
{
    unsigned int pos                        = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int increment                  = gridDim.x * blockDim.x;
    unsigned int rpos                       = cSim.pRandomPos[blockIdx.x];
    PMEDouble boltz                         = 8.31441e-3 * temp / 4.184;
    while (pos < cSim.atoms)
    {  
        PMEDouble invMass                   = cSim.pImageInvMass[pos];    
        PMEDouble forceX                    = cSim.pForceX[pos];
        PMEDouble forceY                    = cSim.pForceY[pos];
        PMEDouble forceZ                    = cSim.pForceZ[pos];
        PMEDouble velX, velY, velZ;
        
        // Zero velocities if it's really cold
        if (temp < 1.0e-6)
        {
            velX                            = (PMEDouble)0.0;
            velY                            = (PMEDouble)0.0;
            velZ                            = (PMEDouble)0.0;
        }
        else
        {
            PMEDouble gaussX                = cSim.pRandomX[rpos + pos];
            PMEDouble gaussY                = cSim.pRandomY[rpos + pos];
            PMEDouble gaussZ                = cSim.pRandomZ[rpos + pos]; 
            PMEDouble sd                    = sqrt(boltz * invMass);        
            velX                            = sd * gaussX;
            velY                            = sd * gaussY;
            velZ                            = sd * gaussZ;        
        }
        
        
        // Back velocities up a half-step
        PMEDouble wfac                      = invMass * half_dtx;
        velX                               -= forceX * wfac;
        velY                               -= forceY * wfac;
        velZ                               -= forceZ * wfac;
        
        // Write final velocities
        cSim.pImageVelX[pos]                = velX;
        cSim.pImageVelY[pos]                = velY;
        cSim.pImageVelZ[pos]                = velZ;
        pos                                += increment;
    }
    
    // Update RNG position
    __syncthreads();
    if (threadIdx.x == 0)
        cSim.pRandomPos[blockIdx.x]         = rpos + cSim.paddedNumberOfAtoms;    
}


void kResetVelocities(gpuContext gpu, double temp, double half_dtx)
{
  
    // Update random numbers if necessary
    if (gpu->randomCounter >= gpu->sim.randomSteps)
    {
#ifdef CPU_RANDOMS
        cpu_kRandom(gpu);
#else
        kRandom(gpu);
#endif
        gpu->randomCounter = 0;
    }
        
    if (gpu->bNeighborList)
        kNLResetVelocities_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(temp, half_dtx);
    else
        kResetVelocities_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(temp, half_dtx);
            
    gpu->randomCounter++; 
    LAUNCHERROR("kResetVelocities"); 
}

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_UPDATE_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kRecalculateVelocities_kernel(PMEDouble dtx_inv)
{
    unsigned int pos                        = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int increment                  = gridDim.x * blockDim.x;

    while (pos < cSim.atoms)
    {  
        PMEDouble oldAtomX                  = cSim.pForceX[pos];
        PMEDouble atomX                     = cSim.pAtomX[pos];
        PMEDouble oldAtomY                  = cSim.pForceY[pos];
        PMEDouble atomY                     = cSim.pAtomY[pos];
        PMEDouble oldAtomZ                  = cSim.pForceZ[pos];
        PMEDouble atomZ                     = cSim.pAtomZ[pos];   
        PMEDouble velX                      = (atomX - oldAtomX) * dtx_inv;
        PMEDouble velY                      = (atomY - oldAtomY) * dtx_inv;
        PMEDouble velZ                      = (atomZ - oldAtomZ) * dtx_inv;
        cSim.pVelX[pos]                     = velX;
        cSim.pVelY[pos]                     = velY;
        cSim.pVelZ[pos]                     = velZ;
        pos                                += increment;
    }
}

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_UPDATE_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMERecalculateVelocities_kernel(PMEDouble dtx_inv)
{
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int increment                              = gridDim.x * blockDim.x;

    while (pos < cSim.atoms)
    {  
        PMEDouble oldAtomX                  = cSim.pForceX[pos];
        PMEDouble atomX                     = cSim.pImageX[pos];
        PMEDouble oldAtomY                  = cSim.pForceY[pos];
        PMEDouble atomY                     = cSim.pImageY[pos];
        PMEDouble oldAtomZ                  = cSim.pForceZ[pos];
        PMEDouble atomZ                     = cSim.pImageZ[pos];   
        PMEDouble velX                      = (atomX - oldAtomX) * dtx_inv;
        PMEDouble velY                      = (atomY - oldAtomY) * dtx_inv;
        PMEDouble velZ                      = (atomZ - oldAtomZ) * dtx_inv;
        cSim.pImageVelX[pos]                = velX;
        cSim.pImageVelY[pos]                = velY;
        cSim.pImageVelZ[pos]                = velZ;
        pos                                += increment;
    }
}


void kRecalculateVelocities(gpuContext gpu, PMEDouble dtx_inv)
{
    if (gpu->bNeighborList)
    {
        kPMERecalculateVelocities_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(dtx_inv);
    }
    else
    {
        kRecalculateVelocities_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(dtx_inv);
    }
    LAUNCHERROR("kRecalculateVelocities");  
}

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kCalculateKineticEnergy_kernel(PMEFloat c_ave)
{
extern __shared__ KineticEnergy sE[];

    PMEFloat eke                            = (PMEFloat)0.0;
    PMEFloat ekph                           = (PMEFloat)0.0;
    PMEFloat ekpbs                          = (PMEFloat)0.0;
    unsigned int pos                        = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Sum up kinetic energies
    while (pos < cSim.atoms)
    {
        PMEFloat mass                       = cSim.pAtomMass[pos];
        PMEFloat vx                         = cSim.pVelX[pos];
        PMEFloat lvx                        = cSim.pLVelX[pos];
        PMEFloat vy                         = cSim.pVelY[pos];
        PMEFloat lvy                        = cSim.pLVelY[pos];
        PMEFloat vz                         = cSim.pVelZ[pos];       
        PMEFloat lvz                        = cSim.pLVelZ[pos];
        PMEFloat svx                        = vx + lvx;
        PMEFloat svy                        = vy + lvy;
        PMEFloat svz                        = vz + lvz;
        eke                                += mass * (svx * svx + svy * svy + svz * svz);
        ekpbs                              += mass * (vx * lvx + vy * lvy + vz * lvz);
        ekph                               += mass * (vx * vx + vy * vy + vz * vz);
        pos                                += blockDim.x * gridDim.x;
    }       
    eke                                    *= (PMEFloat)0.125 * c_ave;
    ekph                                   *= (PMEFloat)0.5;
    ekpbs                                  *= (PMEFloat)0.5;
    sE[threadIdx.x].KE.EKE                  = eke;
    sE[threadIdx.x].KE.EKPH                 = ekph;
    sE[threadIdx.x].KE.EKPBS                = ekpbs;
        

    // Reduce per-thread kinetic energies
    __syncthreads();
    unsigned int m                          = 1;
    while (m < blockDim.x)
    {
        int p                               = threadIdx.x + m;
        eke                                 = ((p < blockDim.x) ? sE[p].KE.EKE : (PMEFloat)0.0f);
        ekph                                = ((p < blockDim.x) ? sE[p].KE.EKPH : (PMEFloat)0.0f);
        ekpbs                               = ((p < blockDim.x) ? sE[p].KE.EKPBS : (PMEFloat)0.0f);
        __syncthreads();
        sE[threadIdx.x].KE.EKE                += eke;
        sE[threadIdx.x].KE.EKPH               += ekph;
        sE[threadIdx.x].KE.EKPBS              += ekpbs;
        __syncthreads();
        m                                  *= 2;
    }

    // Save result
    if (threadIdx.x < 3)
    {
        cSim.pKineticEnergy[blockIdx.x].array[threadIdx.x]
                                            = sE[0].array[threadIdx.x];
    }
}

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kPMECalculateKineticEnergy_kernel(PMEFloat c_ave)
{
extern __shared__ KineticEnergy sE[];

    PMEFloat eke                            = (PMEFloat)0.0;
    PMEFloat ekph                           = (PMEFloat)0.0;
    PMEFloat ekpbs                          = (PMEFloat)0.0;
    unsigned int pos                        = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Sum up kinetic energies
    while (pos < cSim.atoms)
    {
        PMEFloat mass                       = cSim.pImageMass[pos];
        PMEFloat vx                         = cSim.pImageVelX[pos];
        PMEFloat lvx                        = cSim.pImageLVelX[pos];
        PMEFloat vy                         = cSim.pImageVelY[pos];
        PMEFloat lvy                        = cSim.pImageLVelY[pos];
        PMEFloat vz                         = cSim.pImageVelZ[pos];       
        PMEFloat lvz                        = cSim.pImageLVelZ[pos];
        PMEFloat svx                        = vx + lvx;
        PMEFloat svy                        = vy + lvy;
        PMEFloat svz                        = vz + lvz;
        eke                                += mass * (svx * svx + svy * svy + svz * svz);
        ekpbs                              += mass * (vx * lvx + vy * lvy + vz * lvz);
        ekph                               += mass * (vx * vx + vy * vy + vz * vz);
        pos                                += blockDim.x * gridDim.x;
    }       
    eke                                    *= (PMEFloat)0.125 * c_ave;
    ekph                                   *= (PMEFloat)0.5;
    ekpbs                                  *= (PMEFloat)0.5;
    sE[threadIdx.x].KE.EKE                  = eke;
    sE[threadIdx.x].KE.EKPH                 = ekph;
    sE[threadIdx.x].KE.EKPBS                = ekpbs;   

    // Reduce per-thread kinetic energies
    __syncthreads();
    unsigned int m                          = 1;
    while (m < blockDim.x)
    {
        int p                               = threadIdx.x + m;
        eke                                 = ((p < blockDim.x) ? sE[p].KE.EKE : (PMEFloat)0.0);
        ekph                                = ((p < blockDim.x) ? sE[p].KE.EKPH : (PMEFloat)0.0);
        ekpbs                               = ((p < blockDim.x) ? sE[p].KE.EKPBS : (PMEFloat)0.0);
        __syncthreads();
        sE[threadIdx.x].KE.EKE             += eke;
        sE[threadIdx.x].KE.EKPH            += ekph;
        sE[threadIdx.x].KE.EKPBS           += ekpbs;
        __syncthreads();
        m                                  *= 2;
    }       
     
    // Save result
    if (threadIdx.x < 3)
    {
        cSim.pKineticEnergy[blockIdx.x].array[threadIdx.x]
                                            = sE[0].array[threadIdx.x];
    }
}

void kCalculateKineticEnergy(gpuContext gpu, PMEFloat c_ave)
{
    if (gpu->bNeighborList)
    {
        kPMECalculateKineticEnergy_kernel<<<gpu->blocks, gpu->threadsPerBlock, gpu->threadsPerBlock * sizeof(KineticEnergy)>>>(c_ave);
    }
    else
    {
        kCalculateKineticEnergy_kernel<<<gpu->blocks, gpu->threadsPerBlock, gpu->threadsPerBlock * sizeof(KineticEnergy)>>>(c_ave);
    }
    LAUNCHERROR("kCalculateKineticEnergy");
}

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kScaleVelocities_kernel(PMEDouble scale)
{
    unsigned int pos                        = blockIdx.x * blockDim.x + threadIdx.x;
    
    while (pos < cSim.atoms)
    {
        double vx                           = cSim.pVelX[pos];
        double vy                           = cSim.pVelY[pos];
        double vz                           = cSim.pVelZ[pos];
        vx                                 *= scale;
        vy                                 *= scale;
        vz                                 *= scale;
        cSim.pVelX[pos]                     = vx;
        cSim.pVelY[pos]                     = vy;
        cSim.pVelZ[pos]                     = vz;
        pos                                += blockDim.x * gridDim.x;       
    }
}

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kPMEScaleVelocities_kernel(PMEDouble scale)
{
    unsigned int pos                        = blockIdx.x * blockDim.x + threadIdx.x;
    
    while (pos < cSim.atoms)
    {
        double vx                           = cSim.pImageVelX[pos];
        double vy                           = cSim.pImageVelY[pos];
        double vz                           = cSim.pImageVelZ[pos];
        vx                                 *= scale;
        vy                                 *= scale;
        vz                                 *= scale;
        cSim.pImageVelX[pos]                = vx;
        cSim.pImageVelY[pos]                = vy;
        cSim.pImageVelZ[pos]                = vz;
        pos                                += blockDim.x * gridDim.x;       
    }
}

void kScaleVelocities(gpuContext gpu, PMEDouble scale)
{
    if (gpu->bNeighborList)
    {
        kPMEScaleVelocities_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>(scale);
    }
    else
    {
        kScaleVelocities_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>(scale);
    }
    LAUNCHERROR("kScaleVelocities");
}

#define KPMECALCULATECOM_KERNEL kPMECalculateCOM_kernel
#define KPMECALCULATESOLUTECOM_KERNEL kPMECalculateSoluteCOM_kernel
#define KPMECALCULATECOMKINETICENERGY_KERNEL kPMECalculateCOMKineticEnergy_kernel
#define KCALCULATEMOLECULARVIRIAL_KERNEL kCalculateMolecularVirial_kernel
#define KPRESSURESCALECOORDINATES_KERNEL kPressureScaleCoordinates_kernel
#include "kNTPKernels.h"
#undef KPMECALCULATECOM_KERNEL
#undef KPMECALCULATESOLUTECOM_KERNEL
#undef KPMECALCULATECOMKINETICENERGY_KERNEL
#undef KCALCULATEMOLECULARVIRIAL_KERNEL
#undef KPRESSURESCALECOORDINATES_KERNEL

#define NTP_LOTSOFMOLECULES
#define KPMECALCULATECOM_KERNEL kPMECalculateCOMLarge_kernel
#define KPMECALCULATESOLUTECOM_KERNEL kPMECalculateSoluteCOMLarge_kernel
#define KPMECALCULATECOMKINETICENERGY_KERNEL kPMECalculateCOMKineticEnergyLarge_kernel
#define KCALCULATEMOLECULARVIRIAL_KERNEL kCalculateMolecularVirialLarge_kernel
#define KPRESSURESCALECOORDINATES_KERNEL kPressureScaleCoordinatesLarge_kernel
#include "kNTPKernels.h"
#undef KPMECALCULATECOM_KERNEL
#undef KPMECALCULATESOLUTECOM_KERNEL
#undef KPMECALCULATECOMKINETICENERGY_KERNEL
#undef KCALCULATEMOLECULARVIRIAL_KERNEL
#undef KPRESSURESCALECOORDINATES_KERNEL
#undef NTP_LOTSOFMOLECULES

#include "kNTPCalls.h"
#include "kRandom.h"

#define EP_NEIGHBORLIST
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kNLOrientForces_kernel()
#include "kOrientForcesKernel.h"

#define EP_VIRIAL
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kNLOrientForcesVirial_kernel()
#include "kOrientForcesKernel.h"
#undef EP_VIRIAL
#undef EP_NEIGHBORLIST

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kOrientForces_kernel()
#include "kOrientForcesKernel.h"

void kOrientForces(gpuContext gpu)
{
    if (gpu->bNeighborList)
    {
        if (gpu->sim.ntp > 0)
            kNLOrientForcesVirial_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
        else
            kNLOrientForces_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    }
    else
        kOrientForces_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kOrientForces");
}


#define EP_NEIGHBORLIST
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kNLLocalToGlobal_kernel()
#include "kLocalToGlobalKernel.h"
#undef EP_NEIGHBORLIST

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kLocalToGlobal_kernel()
#include "kLocalToGlobalKernel.h"

void kLocalToGlobal(gpuContext gpu)
{
    if (gpu->bNeighborList)
    {
        kNLLocalToGlobal_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    }
    else
        kLocalToGlobal_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kLocalToGlobal");
}

