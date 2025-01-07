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
static __constant__ cudaSimulation cSim;

void SetkCalculateGBNonbondEnergy1Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetkCalculateGBNonBondEnergy1Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_GBNONBONDENERGY1_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_GBNONBONDENERGY1_THREADS_PER_BLOCK, 1)
#endif
kCalculateGBNonbondForces1_kernel()
#include "kCalculateGBNonbondEnergy1.h"

extern "C" void kCalculateGBNonbondForces1(gpuContext gpu)
{
    // Refresh texture if necessary
    
    kCalculateGBNonbondForces1_kernel<<<gpu->blocks, gpu->GBNonbondEnergy1ThreadsPerBlock>>>();   
    LAUNCHERROR("kCalculateGBNonbondForces1");
}

#define GB_ENERGY
__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_GBNONBONDENERGY1_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_GBNONBONDENERGY1_THREADS_PER_BLOCK, 1)
#endif
void kCalculateGBNonbondEnergy1_kernel()
#include "kCalculateGBNonbondEnergy1.h"
#undef GB_ENERGY


extern "C" void kCalculateGBNonbondEnergy1(gpuContext gpu)
{
    // Refresh texture if necessary
  
    kCalculateGBNonbondEnergy1_kernel<<<gpu->blocks, gpu->GBNonbondEnergy1ThreadsPerBlock>>>();   
    LAUNCHERROR("kCalculateGBNonbondEnergy1");
}

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#endif
kReduceGBTemp7_kernel()
{
#ifndef MPI
    bool bIGB2578                                       = (cSim.igb == 2) || (cSim.igb == 5) || (cSim.igb == 7) || (cSim.igb == 8);    
#endif
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x;   
    while (pos < cSim.atoms)
    {
        unsigned int pos1                               = pos;
#ifndef MPI        
        PMEFloat reff_i                                 = cSim.pReffSP[pos];
        PMEFloat psi_i                                  = cSim.pPsi[pos];
        PMEFloat rborn_i                                = cSim.pAtomRBorn[pos];
        PMEFloat qi                                     = cSim.pAtomChargeSP[pos];
#endif        
        PMEDouble sumdeijda_i                           = (PMEDouble)0.0;
        while (pos1 < cSim.stride * cSim.nonbondForceBuffers)
        {
            sumdeijda_i                                += cSim.pSumdeijdaBuffer[pos1];
            pos1                                       += cSim.stride;
        }       
#ifndef MPI        
        // Process Temp7 component       
        PMEFloat expmkf                                = exp(-cSim.gb_kappa * reff_i) * cSim.extdiel_inv;
        PMEFloat dl                                    = cSim.intdiel_inv - expmkf;
        PMEFloat qi2h                                  = (PMEFloat)0.50 * qi * qi;
        PMEFloat qid2h                                 = qi2h * dl;
        sumdeijda_i                                    = -sumdeijda_i + qid2h - cSim.gb_kappa * qi2h * expmkf * reff_i;
        if (cSim.alpb == 0)
        {
            // egb                                        -= qid2h / reff_i;
           
        }
        else
        {
            // egb                                        -= qid2h * (1.0 / reff_i + cSim.one_arad_beta);
            sumdeijda_i                                *= ((PMEFloat)1.0 + cSim.one_arad_beta * reff_i);
        }
         
        if (bIGB2578)
        {
            
            // new onufriev: we have to later scale values by a
            //               alpha,beta,gamma -dependent factor:           
            PMEFloat thi, thi2;
            if (cSim.igb == 8)
            {
                PMEFloat alpha                          = cSim.pgb_alpha[pos];
                PMEFloat gamma                          = cSim.pgb_gamma[pos];
                PMEFloat beta                           = cSim.pgb_beta[pos];
                thi                                     = tanh((alpha + gamma * psi_i * psi_i - beta * psi_i) * psi_i);
                thi2                                    = (alpha + (PMEFloat)3.0 * gamma * psi_i * psi_i - (PMEFloat)2.0 * beta * psi_i) * ((PMEFloat)1.0 - thi * thi) * (rborn_i - cSim.offset) / rborn_i;
            }
            else
            {
                thi                                     = tanh((cSim.gb_alpha + cSim.gb_gamma * psi_i * psi_i - cSim.gb_beta * psi_i) * psi_i);
                thi2                                    = (cSim.gb_alpha + (PMEFloat)3.0 * cSim.gb_gamma * psi_i * psi_i - (PMEFloat)2.0 * cSim.gb_beta * psi_i) * ((PMEFloat)1.0 - thi * thi) * (rborn_i - cSim.offset) / rborn_i;
            }
            sumdeijda_i                                *= thi2;
        }
        cSim.pTemp7[pos]                                = sumdeijda_i;
#else
        cSim.pTemp7a[pos]                               = sumdeijda_i;
#endif
        pos                                            += blockDim.x * gridDim.x;
    }
}

extern "C" void kReduceGBTemp7(gpuContext gpu)
{
    kReduceGBTemp7_kernel<<<gpu->blocks, gpu->reduceBufferThreadsPerBlock>>>();
    LAUNCHERROR("kReduceGBTemp7");
}

#ifdef MPI
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#endif
kProcessGBTemp7_kernel()
{
    bool bIGB2578                                       = (cSim.igb == 2) || (cSim.igb == 5) || (cSim.igb == 7) || (cSim.igb == 8);    
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x;   
    while (pos < cSim.atoms)
    {
        PMEFloat reff_i                                 = cSim.pReffSP[pos];
        PMEFloat psi_i                                  = cSim.pPsi[pos];
        PMEFloat rborn_i                                = cSim.pAtomRBorn[pos];
        PMEFloat qi                                     = cSim.pAtomChargeSP[pos];
        PMEDouble sumdeijda_i                           = cSim.pTemp7a[pos];

        // Process Temp7 component       
        PMEFloat expmkf                                = exp(-cSim.gb_kappa * reff_i) * cSim.extdiel_inv;
        PMEFloat dl                                    = cSim.intdiel_inv - expmkf;
        PMEFloat qi2h                                  = (PMEFloat)0.50 * qi * qi;
        PMEFloat qid2h                                 = qi2h * dl;
        sumdeijda_i                                    = -sumdeijda_i + qid2h - cSim.gb_kappa * qi2h * expmkf * reff_i;
        if (cSim.alpb == 0)
        {
            // egb                                        -= qid2h / reff_i;
           
        }
        else
        {
            // egb                                        -= qid2h * (1.0 / reff_i + cSim.one_arad_beta);
            sumdeijda_i                                *= ((PMEFloat)1.0 + cSim.one_arad_beta * reff_i);
        }
         
        if (bIGB2578)
        {
            
            // new onufriev: we have to later scale values by a
            //               alpha,beta,gamma -dependent factor:
            PMEFloat thi, thi2;
            if (cSim.igb == 8)
            {
                PMEFloat alpha                          = cSim.pgb_alpha[pos];
                PMEFloat gamma                          = cSim.pgb_gamma[pos];
                PMEFloat beta                           = cSim.pgb_beta[pos];
                thi                                     = tanh((alpha + gamma * psi_i * psi_i - beta * psi_i) * psi_i);
                thi2                                    = (alpha + (PMEFloat)3.0 * gamma * psi_i * psi_i - (PMEFloat)2.0 * beta * psi_i) * ((PMEFloat)1.0 - thi * thi) * (rborn_i - cSim.offset) / rborn_i;
            }
            else
            {
                thi                                     = tanh((cSim.gb_alpha + cSim.gb_gamma * psi_i * psi_i - cSim.gb_beta * psi_i) * psi_i);
                thi2                                    = (cSim.gb_alpha + (PMEFloat)3.0 * cSim.gb_gamma * psi_i * psi_i - (PMEFloat)2.0 * cSim.gb_beta * psi_i) * ((PMEFloat)1.0 - thi * thi) * (rborn_i - cSim.offset) / rborn_i;
            }
            sumdeijda_i                                *= thi2;
        }
      
        cSim.pTemp7[pos]                                = sumdeijda_i;
        pos                                            += blockDim.x * gridDim.x;
    }
}

extern "C" void kProcessGBTemp7(gpuContext gpu)
{
    kProcessGBTemp7_kernel<<<gpu->blocks, gpu->reduceBufferThreadsPerBlock>>>();
    LAUNCHERROR("kProcessGBTemp7");
}
#endif


__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#endif
kReduceGBTemp7Energy_kernel()
{
#ifndef MPI
    volatile __shared__ PMEDouble sE[1024];
    bool bIGB2578                                       = (cSim.igb == 2) || (cSim.igb == 5) || (cSim.igb == 7) || (cSim.igb == 8);    
#endif
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x; 
#ifndef MPI
    PMEDouble egb                                       = (PMEDouble)0.0;  
#endif

    while (pos < cSim.atoms)
    {
        unsigned int pos1                               = pos;
#ifndef MPI
        PMEFloat reff_i                                 = cSim.pReffSP[pos];
        PMEFloat psi_i                                  = cSim.pPsi[pos];
        PMEFloat rborn_i                                = cSim.pAtomRBorn[pos];
        PMEFloat qi                                     = cSim.pAtomChargeSP[pos];
#endif        
        PMEDouble sumdeijda_i                           = (PMEDouble)0.0;
        while (pos1 < cSim.stride * cSim.nonbondForceBuffers)
        {
            sumdeijda_i                                += cSim.pSumdeijdaBuffer[pos1];
            pos1                                       += cSim.stride;
        }       
#ifndef MPI        
        // Process Temp7 component       
        PMEFloat expmkf                                = exp(-cSim.gb_kappa * reff_i) * cSim.extdiel_inv;
        PMEFloat dl                                    = cSim.intdiel_inv - expmkf;
        PMEFloat qi2h                                  = (PMEFloat)0.50 * qi * qi;
        PMEFloat qid2h                                 = qi2h * dl;
        sumdeijda_i                                    = -sumdeijda_i + qid2h - cSim.gb_kappa * qi2h * expmkf * reff_i;
        if (cSim.alpb == 0)
        {
           egb                                        -= qid2h / reff_i;
        }
        else
        {
            egb                                        -= qid2h * ((PMEFloat)1.0 / reff_i + cSim.one_arad_beta);
            sumdeijda_i                                *= ((PMEFloat)1.0 + cSim.one_arad_beta * reff_i);
        }
         
        if (bIGB2578)
        {
            
            // new onufriev: we have to later scale values by a
            //               alpha,beta,gamma-dependent factor:
            PMEFloat thi, thi2;
            if (cSim.igb == 8)
            {
                PMEFloat alpha                          = cSim.pgb_alpha[pos];
                PMEFloat gamma                          = cSim.pgb_gamma[pos];
                PMEFloat beta                           = cSim.pgb_beta[pos];
                thi                                     = tanh((alpha + gamma * psi_i * psi_i - beta * psi_i) * psi_i);
                thi2                                    = (alpha + (PMEFloat)3.0 * gamma * psi_i * psi_i - (PMEFloat)2.0 * beta * psi_i) * ((PMEFloat)1.0 - thi * thi) * (rborn_i - cSim.offset) / rborn_i;
            }
            else
            {
                thi                                     = tanh((cSim.gb_alpha + cSim.gb_gamma * psi_i * psi_i - cSim.gb_beta * psi_i) * psi_i);
                thi2                                    = (cSim.gb_alpha + (PMEFloat)3.0 * cSim.gb_gamma * psi_i * psi_i - (PMEFloat)2.0 * cSim.gb_beta * psi_i) * ((PMEFloat)1.0 - thi * thi) * (rborn_i - cSim.offset) / rborn_i;
            }
            sumdeijda_i                                *= thi2;
        }
        cSim.pTemp7[pos]                                = sumdeijda_i;
#else     
        cSim.pTemp7a[pos]                               = sumdeijda_i;
#endif
        pos                                            += blockDim.x * gridDim.x;
    }
    
#ifndef MPI    
    // Reduce Generalized Born energy
    sE[threadIdx.x]                                     = egb;
    __syncthreads();
    unsigned int m                                      = 1;
    while (m < blockDim.x)
    {
        int p                                           = threadIdx.x + m;
        PMEDouble d                                     = ((p < blockDim.x) ? sE[p] : (PMEDouble)0.0);
        __syncthreads();
        sE[threadIdx.x]                                += d;
        __syncthreads();
        m                                              *= 2;
    }
    egb                                                 = sE[threadIdx.x];
    unsigned long long int val                          = (unsigned long long int)(fabs(egb) * ENERGYSCALE + (PMEDouble)0.5);
    if (egb < (PMEDouble)0.0)
        val                                             = 0ull - val;
    if (threadIdx.x == 0)
        atomicAdd(cSim.pEGB, val);  
#endif
}


extern "C" void kReduceGBTemp7Energy(gpuContext gpu)
{
    kReduceGBTemp7Energy_kernel<<<gpu->blocks, gpu->reduceBufferThreadsPerBlock>>>();
    LAUNCHERROR("kReduceGBTemp7Energy");
}

#ifdef MPI
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_REDUCEBUFFER_THREADS_PER_BLOCK, 1)
#endif
kProcessGBTemp7Energy_kernel()
{
    volatile __shared__ PMEDouble sE[1024];
    bool bIGB2578                                       = (cSim.igb == 2) || (cSim.igb == 5) || (cSim.igb == 7) || (cSim.igb == 8);    
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x; 
    PMEDouble egb                                       = (PMEDouble)0.0;  
    while (pos < cSim.atoms)
    {
        PMEFloat reff_i                                 = cSim.pReffSP[pos];
        PMEFloat psi_i                                  = cSim.pPsi[pos];
        PMEFloat rborn_i                                = cSim.pAtomRBorn[pos];
        PMEFloat qi                                     = cSim.pAtomChargeSP[pos];        
        PMEDouble sumdeijda_i                           = cSim.pTemp7a[pos];
      
        // Process Temp7 component       
        PMEFloat expmkf                                = exp(-cSim.gb_kappa * reff_i) * cSim.extdiel_inv;
        PMEFloat dl                                    = cSim.intdiel_inv - expmkf;
        PMEFloat qi2h                                  = (PMEFloat)0.50 * qi * qi;
        PMEFloat qid2h                                 = qi2h * dl;
        sumdeijda_i                                    = -sumdeijda_i + qid2h - cSim.gb_kappa * qi2h * expmkf * reff_i;
        if (cSim.alpb == 0)
        {
           egb                                        -= qid2h / reff_i;
           
        }
        else
        {
            egb                                        -= qid2h * ((PMEFloat)1.0 / reff_i + cSim.one_arad_beta);
            sumdeijda_i                                *= ((PMEFloat)1.0 + cSim.one_arad_beta * reff_i);
        }
         
        if (bIGB2578)
        {
            
            // new onufriev: we have to later scale values by a
            //               alpha,beta,gamma -dependent factor:          
            PMEFloat thi, thi2;
            if (cSim.igb == 8)
            {
                PMEFloat alpha                          = cSim.pgb_alpha[pos];
                PMEFloat gamma                          = cSim.pgb_gamma[pos];
                PMEFloat beta                           = cSim.pgb_beta[pos];
                thi                                     = tanh((alpha + gamma * psi_i * psi_i - beta * psi_i) * psi_i);
                thi2                                    = (alpha + (PMEFloat)3.0 * gamma * psi_i * psi_i - (PMEFloat)2.0 * beta * psi_i) * ((PMEFloat)1.0 - thi * thi) * (rborn_i - cSim.offset) / rborn_i;
            }
            else
            {
                thi                                     = tanh((cSim.gb_alpha + cSim.gb_gamma * psi_i * psi_i - cSim.gb_beta * psi_i) * psi_i);
                thi2                                    = (cSim.gb_alpha + (PMEFloat)3.0 * cSim.gb_gamma * psi_i * psi_i - (PMEFloat)2.0 * cSim.gb_beta * psi_i) * ((PMEFloat)1.0 - thi * thi) * (rborn_i - cSim.offset) / rborn_i;
            }

            sumdeijda_i                                *= thi2;
        }      
        cSim.pTemp7[pos]                                = sumdeijda_i;
        pos                                            += blockDim.x * gridDim.x;
    }
   
    // Reduce Generalized Born energy
    sE[threadIdx.x]                                     = egb;
    __syncthreads();
    unsigned int m                                      = 1;
    while (m < blockDim.x)
    {
        int p                                           = threadIdx.x + m;
        PMEDouble d                                     = ((p < blockDim.x) ? sE[p] : (PMEDouble)0.0);
        __syncthreads();
        sE[threadIdx.x]                                += d;
        __syncthreads();
        m                                              *= 2;
    }
    egb                                                 = sE[threadIdx.x];
    unsigned long long int val                          = (unsigned long long int)(fabs(egb) * ENERGYSCALE + (PMEDouble)0.5);
    if (egb < (PMEDouble)0.0)
        val                                             = 0ull - val;
    if (threadIdx.x == 0)
        atomicAdd(cSim.pEGB, val);  
}


extern "C" void kProcessGBTemp7Energy(gpuContext gpu)
{
    kProcessGBTemp7Energy_kernel<<<gpu->blocks, gpu->reduceBufferThreadsPerBlock>>>();
    LAUNCHERROR("kProcessGBTemp7Energy");
}
#endif

