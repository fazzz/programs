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

//#define PME_VIRIAL
//#define PME_ENERGY


static __constant__ cudaSimulation cSim;

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
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMENonbondForcesVirial32_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMENonbondEnergyVirial32_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY

#define PME_IS_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMEOrthogonalNonbondForcesVirial32_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMEOrthogonalNonbondEnergyVirial32_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMENonbondForces32_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMENonbondEnergy32_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY

#define PME_IS_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMEOrthogonalNonbondForces32_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMEOrthogonalNonbondEnergy32_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP



#define PME_ATOMS_PER_WARP (16)
#define PME_VIRIAL

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMENonbondForcesVirial16_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMENonbondEnergyVirial16_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY

#define PME_IS_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMEOrthogonalNonbondForcesVirial16_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMEOrthogonalNonbondEnergyVirial16_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMENonbondForces16_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMENonbondEnergy16_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY

#define PME_IS_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMEOrthogonalNonbondForces16_kernel()
#include "kNLCPNE.h"

#define PME_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMEOrthogonalNonbondEnergy16_kernel()
#include "kNLCPNE.h"

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP





extern "C" void kCalculatePMENonbondForces(gpuContext gpu)
{   
    if (gpu->sim.NLAtomsPerWarp == 32)
    {
        if (gpu->sim.ntp > 0)
        {
             if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondForcesVirial32_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculatePMENonbondForcesVirial32_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();      
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondForces32_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculatePMENonbondForces32_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
        }    
    }
    else
    {
        if (gpu->sim.ntp > 0)
        {
             if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondForcesVirial16_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculatePMENonbondForcesVirial16_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();      
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondForces16_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculatePMENonbondForces16_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
        }    
    }


    
    
    
#ifdef MPI 
    LAUNCHERROR_NONBLOCKING("kCalculatePMENonbondForces");
#else
    LAUNCHERROR("kCalculatePMENonbondForces");
#endif  
}


extern "C" void kCalculatePMENonbondEnergy(gpuContext gpu)
{
    if (gpu->sim.NLAtomsPerWarp == 32)
    {
        if (gpu->sim.ntp > 0)
        {
            if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondEnergyVirial32_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculatePMENonbondEnergyVirial32_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();    
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondEnergy32_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculatePMENonbondEnergy32_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>(); 
        }    
    }
    else
    {
        if (gpu->sim.ntp > 0)
        {
            if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondEnergyVirial16_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculatePMENonbondEnergyVirial16_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();    
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculatePMEOrthogonalNonbondEnergy16_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculatePMENonbondEnergy16_kernel<<<gpu->blocks, gpu->PMENonbondEnergyThreadsPerBlock>>>(); 
        }
    }
#ifdef MPI 
    LAUNCHERROR_NONBLOCKING("kCalculatePMENonbondEnergy");
#else
    LAUNCHERROR("kCalculatePMENonbondEnergy");
#endif  
}


#define IPS_ATOMS_PER_WARP (32)
#define IPS_VIRIAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSNonbondForcesVirial32_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSNonbondEnergyVirial32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

#define IPS_IS_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSOrthogonalNonbondForcesVirial32_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSOrthogonalNonbondEnergyVirial32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_VIRIAL

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSNonbondForces32_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSNonbondEnergy32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

#define IPS_IS_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSOrthogonalNonbondForces32_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSOrthogonalNonbondEnergy32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_ATOMS_PER_WARP



#define IPS_ATOMS_PER_WARP (16)
#define IPS_VIRIAL

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSNonbondForcesVirial16_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSNonbondEnergyVirial16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

#define IPS_IS_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSOrthogonalNonbondForcesVirial16_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSOrthogonalNonbondEnergyVirial16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_VIRIAL

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSNonbondForces16_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSNonbondEnergy16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

#define IPS_IS_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSOrthogonalNonbondForces16_kernel()
#include "kNLCINE.h"

#define IPS_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK, 1)
#endif
kCalculateIPSOrthogonalNonbondEnergy16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_ATOMS_PER_WARP





extern "C" void kCalculateIPSNonbondForces(gpuContext gpu)
{   
    if (gpu->sim.NLAtomsPerWarp == 32)
    {
        if (gpu->sim.ntp > 0)
        {
             if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondForcesVirial32_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculateIPSNonbondForcesVirial32_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();      
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondForces32_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculateIPSNonbondForces32_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
        }    
    }
    else
    {
        if (gpu->sim.ntp > 0)
        {
             if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondForcesVirial16_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculateIPSNonbondForcesVirial16_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();      
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondForces16_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else    
                kCalculateIPSNonbondForces16_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
        }    
    }


    
    
    
#ifdef MPI 
    LAUNCHERROR_NONBLOCKING("kCalculateIPSNonbondForces");
#else
    LAUNCHERROR("kCalculateIPSNonbondForces");
#endif  
}


extern "C" void kCalculateIPSNonbondEnergy(gpuContext gpu)
{
    if (gpu->sim.NLAtomsPerWarp == 32)
    {
        if (gpu->sim.ntp > 0)
        {
            if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondEnergyVirial32_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculateIPSNonbondEnergyVirial32_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();    
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondEnergy32_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculateIPSNonbondEnergy32_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>(); 
        }    
    }
    else
    {
        if (gpu->sim.ntp > 0)
        {
            if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondEnergyVirial16_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculateIPSNonbondEnergyVirial16_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();    
        }
        else
        {
            if (gpu->sim.is_orthog)
                kCalculateIPSOrthogonalNonbondEnergy16_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>();   
            else
                kCalculateIPSNonbondEnergy16_kernel<<<gpu->blocks, gpu->IPSNonbondEnergyThreadsPerBlock>>>(); 
        }
    }
#ifdef MPI 
    LAUNCHERROR_NONBLOCKING("kCalculateIPSNonbondEnergy");
#else
    LAUNCHERROR("kCalculateIPSNonbondEnergy");
#endif  
}






