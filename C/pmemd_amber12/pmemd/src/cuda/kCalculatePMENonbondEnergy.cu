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

//#define PME_VIRIAL
//#define PME_ENERGY

static __constant__ cudaSimulation cSim;

#ifndef use_DPDP
texture<float2, 1, cudaReadModeElementType> xytexref;
texture<float, 1, cudaReadModeElementType> ztexref;
texture<float, 1, cudaReadModeElementType> qtexref;
texture<float2, 1, cudaReadModeElementType> sigepstexref;
#endif

void SetkCalculatePMENonbondEnergySim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetkCalculatePMENonBondEnergySim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}


#ifndef use_DPDP
static __forceinline__ __device__ float __internal_fmad(float a, float b, float c)
{
#if __CUDA_ARCH__ >= 200
  return __fmaf_rn (a, b, c);
#else /* __CUDA_ARCH__ >= 200 */
  return a * b + c;
#endif /* __CUDA_ARCH__ >= 200 */
}

// Faster ERFC approximation courtesy of Norbert Juffa. NVIDIA Corporation
static __forceinline__ __device__ PMEFloat fasterfc(PMEFloat a) 
{
  /* approximate log(erfc(a)) with rel. error < 7e-9 */
  PMEFloat t, x = a;
  t =                       (PMEFloat)-1.6488499458192755E-006;
  t = __internal_fmad(t, x, (PMEFloat)2.9524665006554534E-005);
  t = __internal_fmad(t, x, (PMEFloat)-2.3341951153749626E-004);
  t = __internal_fmad(t, x, (PMEFloat)1.0424943374047289E-003);
  t = __internal_fmad(t, x, (PMEFloat)-2.5501426008983853E-003);
  t = __internal_fmad(t, x, (PMEFloat)3.1979939710877236E-004);
  t = __internal_fmad(t, x, (PMEFloat)2.7605379075746249E-002);
  t = __internal_fmad(t, x, (PMEFloat)-1.4827402067461906E-001);
  t = __internal_fmad(t, x, (PMEFloat)-9.1844764013203406E-001);
  t = __internal_fmad(t, x, (PMEFloat)-1.6279070384382459E+000);
  t = t * x;
  return exp2f(t);
}
#endif




// Nonbond kernels

#define PME_ATOMS_PER_WARP (32)
#define PME_VIRIAL
__global__ void 
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMENonbondForcesVirial32_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMENonbondEnergyVirial32_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY

#define PME_IS_ORTHOGONAL
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMEOrthogonalNonbondForcesVirial32_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMEOrthogonalNonbondEnergyVirial32_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMENonbondForces32_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMENonbondEnergy32_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY

#define PME_IS_ORTHOGONAL
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMEOrthogonalNonbondForces32_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMEOrthogonalNonbondEnergy32_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP



#define PME_ATOMS_PER_WARP (16)
#define PME_VIRIAL

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMENonbondForcesVirial16_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMENonbondEnergyVirial16_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY

#define PME_IS_ORTHOGONAL
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMEOrthogonalNonbondForcesVirial16_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMEOrthogonalNonbondEnergyVirial16_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMENonbondForces16_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMENonbondEnergy16_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY

#define PME_IS_ORTHOGONAL
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMEOrthogonalNonbondForces16_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculatePMEOrthogonalNonbondEnergy16_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP





extern "C" void kCalculatePMENonbondForces(gpuContext gpu)
{   
    
#ifndef use_DPDP
    // Bind textures
    xytexref.normalized             = 0;
    xytexref.filterMode             = cudaFilterModePoint;
    xytexref.addressMode[0]         = cudaAddressModeClamp;
    xytexref.channelDesc.x          = 32;       
    xytexref.channelDesc.y          = 32;     
    xytexref.channelDesc.z          = 0;
    xytexref.channelDesc.w          = 0;
    cudaBindTexture(NULL, xytexref, gpu->sim.pAtomXYSP, gpu->sim.stride * sizeof(float2));        
    ztexref.normalized              = 0;
    ztexref.filterMode              = cudaFilterModePoint;
    ztexref.addressMode[0]          = cudaAddressModeClamp;
    ztexref.channelDesc.x           = 32;       
    ztexref.channelDesc.y           = 0;     
    ztexref.channelDesc.z           = 0;
    ztexref.channelDesc.w           = 0;
    cudaBindTexture(NULL, ztexref, gpu->sim.pAtomZSP, gpu->sim.stride * sizeof(float));   
    qtexref.normalized              = 0;
    qtexref.filterMode              = cudaFilterModePoint;
    qtexref.addressMode[0]          = cudaAddressModeClamp;
    qtexref.channelDesc.x           = 32;       
    qtexref.channelDesc.y           = 0;     
    qtexref.channelDesc.z           = 0;
    qtexref.channelDesc.w           = 0;
    cudaBindTexture(NULL, qtexref, gpu->sim.pAtomChargeSP, gpu->sim.stride * sizeof(float));     
    sigepstexref.normalized         = 0;
    sigepstexref.filterMode         = cudaFilterModePoint;
    sigepstexref.addressMode[0]     = cudaAddressModeClamp;
    sigepstexref.channelDesc.x      = 32;       
    sigepstexref.channelDesc.y      = 32;     
    sigepstexref.channelDesc.z      = 0;
    sigepstexref.channelDesc.w      = 0;
    cudaBindTexture(NULL, sigepstexref, gpu->sim.pImageSigEps, gpu->sim.stride * sizeof(float2)); 
#endif

    if (gpu->sim.NLAtomsPerWarp == 32)
    {
        if (gpu->sim.ntp > 0)
        {
             if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondForcesVirial32_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculatePMENonbondForcesVirial32_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();      
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondForces32_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculatePMENonbondForces32_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
        }    
    }
    else
    {
        if (gpu->sim.ntp > 0)
        {
             if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondForcesVirial16_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculatePMENonbondForcesVirial16_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();      
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondForces16_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculatePMENonbondForces16_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
        }    
    }
    LAUNCHERROR("kCalculatePMENonbondForces"); 

#ifndef use_DPDP
    // Unbind textures
    cudaUnbindTexture(xytexref);
    cudaUnbindTexture(ztexref);
    cudaUnbindTexture(qtexref);
    cudaUnbindTexture(sigepstexref);
#endif
}


extern "C" void kCalculatePMENonbondEnergy(gpuContext gpu)
{
#ifndef use_DPDP
    // Bind textures
    xytexref.normalized             = 0;
    xytexref.filterMode             = cudaFilterModePoint;
    xytexref.addressMode[0]         = cudaAddressModeClamp;
    xytexref.channelDesc.x          = 32;       
    xytexref.channelDesc.y          = 32;     
    xytexref.channelDesc.z          = 0;
    xytexref.channelDesc.w          = 0;
    cudaBindTexture(NULL, xytexref, gpu->sim.pAtomXYSP, gpu->sim.stride * sizeof(float2));        
    ztexref.normalized              = 0;
    ztexref.filterMode              = cudaFilterModePoint;
    ztexref.addressMode[0]          = cudaAddressModeClamp;
    ztexref.channelDesc.x           = 32;       
    ztexref.channelDesc.y           = 0;     
    ztexref.channelDesc.z           = 0;
    ztexref.channelDesc.w           = 0;
    cudaBindTexture(NULL, ztexref, gpu->sim.pAtomZSP, gpu->sim.stride * sizeof(float));   
    qtexref.normalized              = 0;
    qtexref.filterMode              = cudaFilterModePoint;
    qtexref.addressMode[0]          = cudaAddressModeClamp;
    qtexref.channelDesc.x           = 32;       
    qtexref.channelDesc.y           = 0;     
    qtexref.channelDesc.z           = 0;
    qtexref.channelDesc.w           = 0;
    cudaBindTexture(NULL, qtexref, gpu->sim.pAtomChargeSP, gpu->sim.stride * sizeof(float));   
    sigepstexref.normalized         = 0;
    sigepstexref.filterMode         = cudaFilterModePoint;
    sigepstexref.addressMode[0]     = cudaAddressModeClamp;
    sigepstexref.channelDesc.x      = 32;       
    sigepstexref.channelDesc.y      = 32;     
    sigepstexref.channelDesc.z      = 0;
    sigepstexref.channelDesc.w      = 0;
    cudaBindTexture(NULL, sigepstexref, gpu->sim.pImageSigEps, gpu->sim.stride * sizeof(float2)); 
#endif

    if (gpu->sim.NLAtomsPerWarp == 32)
    {
        if (gpu->sim.ntp > 0)
        {
            if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondEnergyVirial32_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculatePMENonbondEnergyVirial32_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();    
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondEnergy32_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculatePMENonbondEnergy32_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>(); 
        }    
    }
    else
    {
        if (gpu->sim.ntp > 0)
        {
            if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondEnergyVirial16_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculatePMENonbondEnergyVirial16_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();    
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondEnergy16_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculatePMENonbondEnergy16_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondEnergyThreadsPerBlock>>>(); 
        }
    }
    LAUNCHERROR("kCalculatePMENonbondEnergy");

#ifndef use_DPDP
    // Unbind textures
    cudaUnbindTexture(xytexref);
    cudaUnbindTexture(ztexref);
    cudaUnbindTexture(qtexref);
    cudaUnbindTexture(sigepstexref);
#endif
}


#define IPS_ATOMS_PER_WARP (32)
#define IPS_VIRIAL
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSNonbondForcesVirial32_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSNonbondEnergyVirial32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

#define IPS_IS_ORTHOGONAL
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSOrthogonalNonbondForcesVirial32_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSOrthogonalNonbondEnergyVirial32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_VIRIAL

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSNonbondForces32_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSNonbondEnergy32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

#define IPS_IS_ORTHOGONAL
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSOrthogonalNonbondForces32_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSOrthogonalNonbondEnergy32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_ATOMS_PER_WARP



#define IPS_ATOMS_PER_WARP (16)
#define IPS_VIRIAL

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSNonbondForcesVirial16_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSNonbondEnergyVirial16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

#define IPS_IS_ORTHOGONAL
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSOrthogonalNonbondForcesVirial16_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSOrthogonalNonbondEnergyVirial16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_VIRIAL

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSNonbondForces16_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSNonbondEnergy16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

#define IPS_IS_ORTHOGONAL
__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSOrthogonalNonbondForces16_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
kCalculateIPSOrthogonalNonbondEnergy16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_ATOMS_PER_WARP

extern "C" void kCalculatePMENonbondEnergyInitKernels(gpuContext gpu)
{
    cudaFuncSetSharedMemConfig(kCalculatePMENonbondForcesVirial32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculatePMENonbondEnergyVirial32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculatePMEOrthogonalNonbondForcesVirial32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculatePMEOrthogonalNonbondEnergyVirial32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculatePMENonbondForces32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculatePMENonbondEnergy32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculatePMEOrthogonalNonbondForces32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculatePMEOrthogonalNonbondEnergy32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculatePMENonbondForcesVirial16_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculatePMENonbondEnergyVirial16_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculatePMEOrthogonalNonbondForcesVirial16_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculatePMEOrthogonalNonbondEnergyVirial16_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculatePMENonbondForces16_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculatePMENonbondEnergy16_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculatePMEOrthogonalNonbondForces16_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculatePMEOrthogonalNonbondEnergy16_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSNonbondForcesVirial32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSNonbondEnergyVirial32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSOrthogonalNonbondForcesVirial32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSOrthogonalNonbondEnergyVirial32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSNonbondForces32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSNonbondEnergy32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSOrthogonalNonbondForces32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSOrthogonalNonbondEnergy32_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSNonbondForcesVirial16_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSNonbondEnergyVirial16_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSOrthogonalNonbondForcesVirial16_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSOrthogonalNonbondEnergyVirial16_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSNonbondForces16_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSNonbondEnergy16_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSOrthogonalNonbondForces16_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kCalculateIPSOrthogonalNonbondEnergy16_kernel, cudaSharedMemBankSizeEightByte);
}

extern "C" void kCalculateIPSNonbondForces(gpuContext gpu)
{  
#ifndef use_DPDP
    // Bind textures
    xytexref.normalized             = 0;
    xytexref.filterMode             = cudaFilterModePoint;
    xytexref.addressMode[0]         = cudaAddressModeClamp;
    xytexref.channelDesc.x          = 32;       
    xytexref.channelDesc.y          = 32;     
    xytexref.channelDesc.z          = 0;
    xytexref.channelDesc.w          = 0;
    cudaBindTexture(NULL, xytexref, gpu->sim.pAtomXYSP, gpu->sim.stride * sizeof(float2));        
    ztexref.normalized              = 0;
    ztexref.filterMode              = cudaFilterModePoint;
    ztexref.addressMode[0]          = cudaAddressModeClamp;
    ztexref.channelDesc.x           = 32;       
    ztexref.channelDesc.y           = 0;     
    ztexref.channelDesc.z           = 0;
    ztexref.channelDesc.w           = 0;
    cudaBindTexture(NULL, ztexref, gpu->sim.pAtomZSP, gpu->sim.stride * sizeof(float));   
    qtexref.normalized              = 0;
    qtexref.filterMode              = cudaFilterModePoint;
    qtexref.addressMode[0]          = cudaAddressModeClamp;
    qtexref.channelDesc.x           = 32;       
    qtexref.channelDesc.y           = 0;     
    qtexref.channelDesc.z           = 0;
    qtexref.channelDesc.w           = 0;
    cudaBindTexture(NULL, qtexref, gpu->sim.pAtomChargeSP, gpu->sim.stride * sizeof(float));     
    sigepstexref.normalized         = 0;
    sigepstexref.filterMode         = cudaFilterModePoint;
    sigepstexref.addressMode[0]     = cudaAddressModeClamp;
    sigepstexref.channelDesc.x      = 32;       
    sigepstexref.channelDesc.y      = 32;     
    sigepstexref.channelDesc.z      = 0;
    sigepstexref.channelDesc.w      = 0;
    cudaBindTexture(NULL, sigepstexref, gpu->sim.pImageSigEps, gpu->sim.stride * sizeof(float2)); 
#endif 

    if (gpu->sim.NLAtomsPerWarp == 32)
    {
        if (gpu->sim.ntp > 0)
        {
             if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondForcesVirial32_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculateIPSNonbondForcesVirial32_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();      
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondForces32_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculateIPSNonbondForces32_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
        }    
    }
    else
    {
        if (gpu->sim.ntp > 0)
        {
             if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondForcesVirial16_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculateIPSNonbondForcesVirial16_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();      
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondForces16_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculateIPSNonbondForces16_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
        }    
    }    
    LAUNCHERROR("kCalculateIPSNonbondForces");

#ifndef use_DPDP
    // Unbind textures
    cudaUnbindTexture(xytexref);
    cudaUnbindTexture(ztexref);
    cudaUnbindTexture(qtexref);
    cudaUnbindTexture(sigepstexref);
#endif  
}


extern "C" void kCalculateIPSNonbondEnergy(gpuContext gpu)
{
#ifndef use_DPDP
    // Bind textures
    xytexref.normalized             = 0;
    xytexref.filterMode             = cudaFilterModePoint;
    xytexref.addressMode[0]         = cudaAddressModeClamp;
    xytexref.channelDesc.x          = 32;       
    xytexref.channelDesc.y          = 32;     
    xytexref.channelDesc.z          = 0;
    xytexref.channelDesc.w          = 0;
    cudaBindTexture(NULL, xytexref, gpu->sim.pAtomXYSP, gpu->sim.stride * sizeof(float2));        
    ztexref.normalized              = 0;
    ztexref.filterMode              = cudaFilterModePoint;
    ztexref.addressMode[0]          = cudaAddressModeClamp;
    ztexref.channelDesc.x           = 32;       
    ztexref.channelDesc.y           = 0;     
    ztexref.channelDesc.z           = 0;
    ztexref.channelDesc.w           = 0;
    cudaBindTexture(NULL, ztexref, gpu->sim.pAtomZSP, gpu->sim.stride * sizeof(float));   
    qtexref.normalized              = 0;
    qtexref.filterMode              = cudaFilterModePoint;
    qtexref.addressMode[0]          = cudaAddressModeClamp;
    qtexref.channelDesc.x           = 32;       
    qtexref.channelDesc.y           = 0;     
    qtexref.channelDesc.z           = 0;
    qtexref.channelDesc.w           = 0;
    cudaBindTexture(NULL, qtexref, gpu->sim.pAtomChargeSP, gpu->sim.stride * sizeof(float));     
    sigepstexref.normalized         = 0;
    sigepstexref.filterMode         = cudaFilterModePoint;
    sigepstexref.addressMode[0]     = cudaAddressModeClamp;
    sigepstexref.channelDesc.x      = 32;       
    sigepstexref.channelDesc.y      = 32;     
    sigepstexref.channelDesc.z      = 0;
    sigepstexref.channelDesc.w      = 0;
    cudaBindTexture(NULL, sigepstexref, gpu->sim.pImageSigEps, gpu->sim.stride * sizeof(float2)); 
#endif

    if (gpu->sim.NLAtomsPerWarp == 32)
    {
        if (gpu->sim.ntp > 0)
        {
            if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondEnergyVirial32_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculateIPSNonbondEnergyVirial32_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();    
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondEnergy32_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculateIPSNonbondEnergy32_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>(); 
        }    
    }
    else
    {
        if (gpu->sim.ntp > 0)
        {
            if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondEnergyVirial16_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculateIPSNonbondEnergyVirial16_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();    
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondEnergy16_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculateIPSNonbondEnergy16_kernel<<<gpu->IPSNonbondBlocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>(); 
        }
    }
    LAUNCHERROR("kCalculateIPSNonbondEnergy");
  
#ifndef use_DPDP
    // Unbind textures
    cudaUnbindTexture(xytexref);
    cudaUnbindTexture(ztexref);
    cudaUnbindTexture(qtexref);
    cudaUnbindTexture(sigepstexref);
#endif
}






