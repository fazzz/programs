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
static __constant__ PMEFloat te                 = (PMEFloat)(4.0 / 3.0);
static __constant__ PMEFloat tf                 = (PMEFloat)(12.0 / 5.0);
static __constant__ PMEFloat tg                 = (PMEFloat)(24.0 / 7.0);
static __constant__ PMEFloat th                 = (PMEFloat)(40.0 / 9.0);
static __constant__ PMEFloat thh                = (PMEFloat)(60.0 / 11.0);

void SetkCalculateGBNonbondEnergy2Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetkCalculateGBNonBondEnergy2Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

__global__ void 
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_GBNONBONDENERGY2_THREADS_PER_BLOCK, SM_3X_GBNONBONDENERGY2_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_GBNONBONDENERGY2_THREADS_PER_BLOCK, SM_2X_GBNONBONDENERGY2_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_GBNONBONDENERGY2_THREADS_PER_BLOCK, SM_13_GBNONBONDENERGY2_BLOCKS_MULTIPLIER)
#endif
kCalculateGBNonbondEnergy2_kernel()
#include "kCalculateGBNonbondEnergy2.h"

#define GB_IGB78
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_GBNONBONDENERGY2_THREADS_PER_BLOCK, SM_3X_GBNONBONDENERGY2_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_GBNONBONDENERGY2_THREADS_PER_BLOCK, SM_2X_GBNONBONDENERGY2_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_GBNONBONDENERGY2_THREADS_PER_BLOCK, SM_13_GBNONBONDENERGY2_BLOCKS_MULTIPLIER)
#endif
kCalculateGBNonbondEnergy2IGB78_kernel()
#include "kCalculateGBNonbondEnergy2.h"
#undef IGB78

void kCalculateGBNonbondEnergy2InitKernels(gpuContext gpu)
{
    cudaFuncSetSharedMemConfig(kCalculateGBNonbondEnergy2_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateGBNonbondEnergy2IGB78_kernel, cudaSharedMemBankSizeEightByte);
    if (gpu->sm_version >= SM_3X)
    {
        cudaFuncSetCacheConfig(kCalculateGBNonbondEnergy2_kernel, cudaFuncCachePreferEqual);
    }
}

void kCalculateGBNonbondEnergy2(gpuContext gpu)
{
    if ((gpu->sim.igb == 7) || (gpu->sim.igb == 8))
        kCalculateGBNonbondEnergy2IGB78_kernel<<<gpu->GBNonbondEnergy2Blocks, gpu->GBNonbondEnergy2IGB78ThreadsPerBlock>>>();   
    else
        kCalculateGBNonbondEnergy2_kernel<<<gpu->GBNonbondEnergy2Blocks, gpu->GBNonbondEnergy2ThreadsPerBlock>>>();
    LAUNCHERROR("kCalculateGBNonbondEnergy2");
}
