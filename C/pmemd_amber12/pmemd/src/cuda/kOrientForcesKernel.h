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



{
struct Virial
{
    PMEDouble v11;
    PMEDouble v22;
    PMEDouble v33;
};
#ifdef EP_VIRIAL
#if (__CUDA_ARCH__ >= 300)
__shared__ Virial sV[SM_3X_THREADS_PER_BLOCK];
#elif (__CUDA_ARCH__ >= 200)
__shared__ Virial sV[SM_2X_THREADS_PER_BLOCK];
#else
__shared__ Virial sV[SM_13_THREADS_PER_BLOCK];
#endif
    PMEDouble v11                           = 0.0;
    PMEDouble v22                           = 0.0;
    PMEDouble v33                           = 0.0;
#endif

    unsigned int pos                        = blockIdx.x * blockDim.x + threadIdx.x;
#define EP_ONEPOINT
    while (pos < cSim.EP11Offset)
    {    
#define EP_TYPE1
        if (pos < cSim.EP11s)
        {
#include "kOrientForcesLoop.h"
        }
#undef EP_TYPE1
        pos                                += blockDim.x * gridDim.x;
    }
    
    while (pos < cSim.EP12Offset)
    {
        pos                                -= cSim.EP11Offset;
#define EP_TYPE2
        if (pos < cSim.EP12s)
        {    
#include "kOrientForcesLoop.h"
        }
#undef EP_TYPE2
        pos                                += cSim.EP11Offset + blockDim.x * gridDim.x;
    }
#undef EP_ONEPOINT

#define EP_TWOPOINTS
    while (pos < cSim.EP21Offset)
    {  
        pos                                -= cSim.EP12Offset;      
#define EP_TYPE1
        if (pos < cSim.EP21s)
        {
#include "kOrientForcesLoop.h"
        }
#undef EP_TYPE1
        pos                                += cSim.EP12Offset + blockDim.x * gridDim.x;
    }
    
    while (pos < cSim.EP22Offset)
    {
        pos                                -= cSim.EP21Offset;
#define EP_TYPE2
        if (pos < cSim.EP22s)
        {    
#include "kOrientForcesLoop.h"
        }
#undef EP_TYPE2
        pos                                += cSim.EP21Offset + blockDim.x * gridDim.x;
    }
#undef EP_TWOPOINTS   

#ifdef EP_VIRIAL
    sV[threadIdx.x].v11                     = v11;
    sV[threadIdx.x].v22                     = v22;
    sV[threadIdx.x].v33                     = v33;
    __syncthreads();
    unsigned int m                          = 1;    
    while (m < blockDim.x)
    {
        int p                               = threadIdx.x + m;
        PMEDouble v11                       = ((p < blockDim.x) ? sV[p].v11 : (PMEDouble)0.0);
        PMEDouble v22                       = ((p < blockDim.x) ? sV[p].v22 : (PMEDouble)0.0);
        PMEDouble v33                       = ((p < blockDim.x) ? sV[p].v33 : (PMEDouble)0.0);
        __syncthreads();
        sV[threadIdx.x].v11                += v11;
        sV[threadIdx.x].v11                += v22;
        sV[threadIdx.x].v33                += v33;        
        __syncthreads();
        m                                  *= 2;
    } 
    if (threadIdx.x == 0)
    {   
        unsigned long long int val11            = llitoulli(lliroundd(sV[threadIdx.x].v11 * FORCESCALE));
        unsigned long long int val22            = llitoulli(lliroundd(sV[threadIdx.x].v22 * FORCESCALE));
        unsigned long long int val33            = llitoulli(lliroundd(sV[threadIdx.x].v33 * FORCESCALE));    
        atomicAdd(cSim.pVirial_11, val11);
        atomicAdd(cSim.pVirial_22, val22);  
        atomicAdd(cSim.pVirial_33, val33);
    }
#endif

}

