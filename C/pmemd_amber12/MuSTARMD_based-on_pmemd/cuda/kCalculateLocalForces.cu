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

//extern __shared__ Vectors sV[];
static __constant__ cudaSimulation cSim;
static __constant__ PMEDouble pt999             = (PMEDouble)(0.9990);
static __constant__ PMEDouble tm06              = (PMEDouble)(1.0e-06);
static __constant__ PMEDouble tenm3             = (PMEDouble)(1.0e-03);
static __constant__ PMEDouble tm24              = (PMEDouble)(1.0e-18);
static __constant__ PMEDouble one               = (PMEDouble)(1.0);
static __constant__ PMEDouble zero              = (PMEDouble)(0.0);
static __constant__ PMEFloat rad_to_deg_coeff   = (PMEDouble)180.0 / ((PMEDouble)CMAPSTEPSIZE * (PMEDouble)PI_VAL);

// Texture reference for PMEDouble-precision coordinates (disguised as int2 to work around HW limitations)
#ifndef use_SPSP
texture<int2, 1, cudaReadModeElementType> texref;
#else
texture<float, 1, cudaReadModeElementType> texref;
#endif

void SetkCalculateLocalForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetkCalculateLocalForcesSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

/*
 * This special version of sincos is designed for |a| < 6*PI. On a GTX 285
 * it is about 25% faster than sincos from the CUDA math library. Also uses
 * 8 fewer registers than the CUDA math library's sincos. Maximum observed
 * error is 2 ulps across range stated above. Infinities and negative zero
 * are not handled according to C99 specifications. NaNs are handled fine.
 */
__device__ void faster_sincos(double a, double *sptr, double *cptr) {
  double t, u, s, c, j, a2;
  int i;

  i = __double2int_rn (a * 6.3661977236758138e-1);
  j = (double)i;
  a = __fma_rn (-j, 1.57079632679489660e+000, a); /* PIO2_HI */
  a = __fma_rn (-j, 6.12323399573676600e-017, a); /* PIO2_LO */
  a2 = a * a;
  u =                  -1.136788825395985E-011;   
  u = __fma_rn (u, a2,  2.087588480545065E-009);
  u = __fma_rn (u, a2, -2.755731555403950E-007);
  u = __fma_rn (u, a2,  2.480158729365970E-005);
  u = __fma_rn (u, a2, -1.388888888888074E-003);
  u = __fma_rn (u, a2,  4.166666666666664E-002);
  u = __fma_rn (u, a2, -5.000000000000000E-001);
  u = __fma_rn (u, a2,  1.000000000000000E+000);
  t =                   1.5896230157221844E-010;
  t = __fma_rn (t, a2, -2.5050747762850355E-008);
  t = __fma_rn (t, a2,  2.7557313621385676E-006);
  t = __fma_rn (t, a2, -1.9841269829589539E-004);
  t = __fma_rn (t, a2,  8.3333333333221182E-003);
  t = __fma_rn (t, a2, -1.6666666666666630E-001);
  t = t * a2;
  t = __fma_rn(t, a, a);
  if (i & 1) {
    s = u;
    c = t;
  } else {
    s = t;
    c = u;
  }
  if (i & 2) {
    s = -s;
  }
  i++;
  if (i & 2) {
    c = -c;
  }
  *sptr = s;
  *cptr = c;
}

struct Energy {
    PMEDouble bond;
    PMEDouble angle;
    PMEDouble dihedral;
    PMEDouble el14;
    PMEDouble nb14;
    PMEDouble restraint;
};

struct Virial {
    PMEDouble vir_11;
    PMEDouble vir_22;
    PMEDouble vir_33;
};

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculateLocalForces_kernel()
#include "kCLF.h"

#define LOCAL_ENERGY
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculateLocalEnergy_kernel()
#include "kCLF.h"
#undef LOCAL_ENERGY

#define LOCAL_NEIGHBORLIST
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMELocalForces_kernel()
#include "kCLF.h"

#define LOCAL_ENERGY
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMELocalEnergy_kernel()
#include "kCLF.h"
#undef LOCAL_ENERGY

#define LOCAL_VIRIAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMELocalForcesVirial_kernel()
#include "kCLF.h"

#define LOCAL_ENERGY
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMELocalEnergyVirial_kernel()
#include "kCLF.h"
#undef LOCAL_ENERGY
#undef LOCAL_VIRIAL
#undef LOCAL_NEIGHBORLIST

// Consumer Fermi kernels
#define NODPTEXTURE
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculateLocalForcesFermi_kernel()
#include "kCLF.h"

#define LOCAL_ENERGY
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculateLocalEnergyFermi_kernel()
#include "kCLF.h"
#undef LOCAL_ENERGY

#define LOCAL_NEIGHBORLIST
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMELocalForcesFermi_kernel()
#include "kCLF.h"

#define LOCAL_ENERGY
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMELocalEnergyFermi_kernel()
#include "kCLF.h"
#undef LOCAL_ENERGY

#define LOCAL_VIRIAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMELocalForcesVirialFermi_kernel()
#include "kCLF.h"

#define LOCAL_ENERGY
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMELocalEnergyVirialFermi_kernel()
#include "kCLF.h"
#undef LOCAL_ENERGY
#undef LOCAL_VIRIAL
#undef LOCAL_NEIGHBORLIST
#undef NODPTEXTURE


extern "C" void kCalculateLocalForces(gpuContext gpu)
{
    if (gpu->bLocalInteractions)
    {
        cudaFuncSetCacheConfig(kCalculateLocalForces_kernel, cudaFuncCachePreferL1);
        cudaFuncSetCacheConfig(kCalculatePMELocalForces_kernel, cudaFuncCachePreferL1);
        cudaFuncSetCacheConfig(kCalculatePMELocalForcesVirial_kernel, cudaFuncCachePreferL1);

#ifndef use_SPSP        
        // Pathetic kludge to workaround GTX4xx/GTX5xx HW bug
        if (!gpu->bECCSupport && (gpu->sm_version == SM_2X))
        {
            if (gpu->bNeighborList)
            {
                if (gpu->sim.ntp > 0)
                    kCalculatePMELocalForcesVirialFermi_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
                else    
                    kCalculatePMELocalForcesFermi_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
            }
            else
                kCalculateLocalForcesFermi_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
#ifdef MPI 
            LAUNCHERROR_NONBLOCKING("kCalculateLocalForcesFermi");
#else
            LAUNCHERROR("kCalculateLocalForcesFermi");
#endif          
        }
        else
#endif
        {      
            texref.normalized       = 0;
            texref.filterMode       = cudaFilterModePoint;
            texref.addressMode[0]   = cudaAddressModeClamp;
            texref.channelDesc.x    = 32;
#ifndef use_SPSP        
            texref.channelDesc.y    = 32;
#else
            texref.channelDesc.y    = 0;
#endif        
            texref.channelDesc.z    = 0;
            texref.channelDesc.w    = 0;
#ifndef use_SPSP       
            int2* pX;
            if (gpu->bNeighborList)
                pX                  = (int2*)gpu->sim.pImageX;
            else
                pX                  = (int2*)gpu->sim.pAtomX;
            cudaBindTexture(NULL, texref, pX, gpu->sim.stride3 * sizeof(int2));
#else
            float* pX;
            if (gpu->bNeighborList)
                pX                  = (float*)gpu->sim.pImageX;
            else
                pX                  = (float*)gpu->sim.pAtomX;
            cudaBindTexture(NULL, texref, pX, gpu->sim.stride3 * sizeof(float));
#endif            
            if (gpu->bNeighborList)
            {
                if (gpu->sim.ntp > 0)
                    kCalculatePMELocalForcesVirial_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
                else    
                    kCalculatePMELocalForces_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
            }
            else
                kCalculateLocalForces_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
#ifdef MPI 
            LAUNCHERROR_NONBLOCKING("kCalculateLocalForces");
#else
            LAUNCHERROR("kCalculateLocalForces");
#endif  
            cudaUnbindTexture(texref);
        }
    }
}

extern "C" void kCalculateLocalEnergy(gpuContext gpu)
{
    if (gpu->bLocalInteractions)
    {
#ifndef use_SPSP    
        if (!gpu->bECCSupport && (gpu->sm_version == SM_2X))
        {
            if (gpu->bNeighborList)
            {
                if (gpu->sim.ntp > 0)
                    kCalculatePMELocalEnergyVirialFermi_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
                else    
                    kCalculatePMELocalEnergyFermi_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
            }
            else
                kCalculateLocalEnergyFermi_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
#ifdef MPI 
            LAUNCHERROR_NONBLOCKING("kCalculateLocalEnergyFermi");
#else
            LAUNCHERROR("kCalculateLocalEnergyFermi");
#endif          
        }
        else
#endif
        {
            texref.normalized       = 0;
            texref.filterMode       = cudaFilterModePoint;
            texref.addressMode[0]   = cudaAddressModeClamp;
            texref.channelDesc.x    = 32;
#ifndef use_SPSP        
            texref.channelDesc.y    = 32;
#else
            texref.channelDesc.y    = 0;
#endif        
            texref.channelDesc.z    = 0;
            texref.channelDesc.w    = 0;
        
#ifndef use_SPSP       
            int2* pX;
            if (gpu->bNeighborList)
                pX                  = (int2*)gpu->sim.pImageX;
            else
                pX                  = (int2*)gpu->sim.pAtomX;
            cudaBindTexture(NULL, texref, pX, gpu->sim.stride3 * sizeof(int2));
#else
            float* pX;
            if (gpu->bNeighborList)
                pX                  = (float*)gpu->sim.pImageX;
            else
                pX                  = (float*)gpu->sim.pAtomX;
            cudaBindTexture(NULL, texref, pX, gpu->sim.stride3 * sizeof(float));
#endif        
            if (gpu->bNeighborList)
            {
                if (gpu->sim.ntp > 0)
                    kCalculatePMELocalEnergyVirial_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
                else    
                    kCalculatePMELocalEnergy_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
            }
            else
                kCalculateLocalEnergy_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
#ifdef MPI 
            LAUNCHERROR_NONBLOCKING("kCalculateLocalEnergy");
#else
            LAUNCHERROR("kCalculateLocalEnergy");
#endif 
            cudaUnbindTexture(texref);
        }
    }
}



#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#endif
__global__ void kCalculateCHARMMForces_kernel()
#include "kCCF.h"

#define CHARMM_ENERGY
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#endif
__global__ void kCalculateCHARMMEnergy_kernel()
#include "kCCF.h"
#undef CHARMM_ENERGY

#define CHARMM_NEIGHBORLIST
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#endif
__global__ void kCalculatePMECHARMMForces_kernel()
#include "kCCF.h"

#define CHARMM_ENERGY
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#endif
__global__ void kCalculatePMECHARMMEnergy_kernel()
#include "kCCF.h"
#undef CHARMM_ENERGY

#define CHARMM_VIRIAL
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#endif
__global__ void kCalculatePMECHARMMForcesVirial_kernel()
#include "kCCF.h"

#define CHARMM_ENERGY
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#endif
__global__ void kCalculatePMECHARMMEnergyVirial_kernel()
#include "kCCF.h"
#undef CHARMM_ENERGY
#undef CHARMM_VIRIAL
#undef CHARMM_NEIGHBORLIST


#define NODPTEXTURE
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#endif
__global__ void kCalculateCHARMMForcesFermi_kernel()
#include "kCCF.h"

#define CHARMM_ENERGY
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#endif
__global__ void kCalculateCHARMMEnergyFermi_kernel()
#include "kCCF.h"
#undef CHARMM_ENERGY

#define CHARMM_NEIGHBORLIST
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#endif
__global__ void kCalculatePMECHARMMForcesFermi_kernel()
#include "kCCF.h"

#define CHARMM_ENERGY
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#endif
__global__ void kCalculatePMECHARMMEnergyFermi_kernel()
#include "kCCF.h"
#undef CHARMM_ENERGY

#define CHARMM_VIRIAL
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#endif
__global__ void kCalculatePMECHARMMForcesVirialFermi_kernel()
#include "kCCF.h"

#define CHARMM_ENERGY
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_CHARMMFORCES_THREADS_PER_BLOCK, 1)
#endif
__global__ void kCalculatePMECHARMMEnergyVirialFermi_kernel()
#include "kCCF.h"
#undef CHARMM_ENERGY
#undef CHARMM_NEIGHBORLIST
#undef NODPTEXTURE



extern "C" void kCalculateCHARMMForces(gpuContext gpu)
{
    if (gpu->bCharmmInteractions)
    {
#ifndef use_SPSP
        if (!gpu->bECCSupport && (gpu->sm_version == SM_2X))
        {
            if (gpu->bNeighborList)
            {
                if (gpu->sim.ntp > 0)
                    kCalculatePMECHARMMForcesVirialFermi_kernel<<<gpu->blocks, gpu->CHARMMForcesThreadsPerBlock>>>();
                else    
                    kCalculatePMECHARMMForcesFermi_kernel<<<gpu->blocks, gpu->CHARMMForcesThreadsPerBlock>>>();
            }
            else
                kCalculateCHARMMForcesFermi_kernel<<<gpu->blocks, gpu->CHARMMForcesThreadsPerBlock>>>();
#ifdef MPI 
            LAUNCHERROR_NONBLOCKING("kCalculateCHARMMForcesFermi");
#else
            LAUNCHERROR("kCalculateCHARMMForcesFermi");
#endif         
        }
        else
#endif
        {
            texref.normalized       = 0;
            texref.filterMode       = cudaFilterModePoint;
            texref.addressMode[0]   = cudaAddressModeClamp;
            texref.channelDesc.x    = 32;
#ifndef use_SPSP        
            texref.channelDesc.y    = 32;
#else
            texref.channelDesc.y    = 0;
#endif        
            texref.channelDesc.z    = 0;
            texref.channelDesc.w    = 0;
#ifndef use_SPSP       
            int2* pX;
            if (gpu->bNeighborList)
                pX                  = (int2*)gpu->sim.pImageX;
            else
                pX                  = (int2*)gpu->sim.pAtomX;
            cudaBindTexture(NULL, texref, pX, gpu->sim.stride3 * sizeof(int2));
#else
            float* pX;
            if (gpu->bNeighborList)
                pX                  = (float*)gpu->sim.pImageX;
            else
                pX                  = (float*)gpu->sim.pAtomX;
            cudaBindTexture(NULL, texref, pX, gpu->sim.stride3 * sizeof(float));
#endif            
            if (gpu->bNeighborList)
            {
                if (gpu->sim.ntp > 0)
                    kCalculatePMECHARMMForcesVirial_kernel<<<gpu->blocks, gpu->CHARMMForcesThreadsPerBlock>>>();
                else
                    kCalculatePMECHARMMForces_kernel<<<gpu->blocks, gpu->CHARMMForcesThreadsPerBlock>>>();
            }
            else
                kCalculateCHARMMForces_kernel<<<gpu->blocks, gpu->CHARMMForcesThreadsPerBlock>>>();
#ifdef MPI 
            LAUNCHERROR_NONBLOCKING("kCalculateCHARMMForces");
#else
            LAUNCHERROR("kCalculateCHARMMForces");
#endif 
            cudaUnbindTexture(texref);
        }
    }
}

extern "C" void kCalculateCHARMMEnergy(gpuContext gpu)
{    
    if (gpu->bCharmmInteractions)
    {
#ifndef use_SPSP   
        if (!gpu->bECCSupport && (gpu->sm_version == SM_2X))
        {
            if (gpu->bNeighborList)
            {
                if (gpu->sim.ntp > 0)
                    kCalculatePMECHARMMEnergyVirialFermi_kernel<<<gpu->blocks, gpu->CHARMMForcesThreadsPerBlock>>>();
                else
                    kCalculatePMECHARMMEnergyFermi_kernel<<<gpu->blocks, gpu->CHARMMForcesThreadsPerBlock>>>();
            }
            else
                kCalculateCHARMMEnergyFermi_kernel<<<gpu->blocks, gpu->CHARMMForcesThreadsPerBlock>>>();
#ifdef MPI 
            LAUNCHERROR_NONBLOCKING("kCalculateCHARMMEnergyFermi");
#else
            LAUNCHERROR("kCalculateCHARMMEnergyFermi");
#endif         
        }
        else
#endif
        {    
            texref.normalized       = 0;
            texref.filterMode       = cudaFilterModePoint;
            texref.addressMode[0]   = cudaAddressModeClamp;
            texref.channelDesc.x    = 32;
#ifndef use_SPSP        
            texref.channelDesc.y    = 32;
#else
            texref.channelDesc.y    = 0;
#endif        
            texref.channelDesc.z    = 0;
            texref.channelDesc.w    = 0;
        
#ifndef use_SPSP       
            int2* pX;
            if (gpu->bNeighborList)
                pX                  = (int2*)gpu->sim.pImageX;
            else
                pX                  = (int2*)gpu->sim.pAtomX;
            cudaBindTexture(NULL, texref, pX, gpu->sim.stride3 * sizeof(int2));
#else
            float* pX;
            if (gpu->bNeighborList)
                pX                  = (float*)gpu->sim.pImageX;
            else
                pX                  = (float*)gpu->sim.pAtomX;
            cudaBindTexture(NULL, texref, pX, gpu->sim.stride3 * sizeof(float));
#endif        
            if (gpu->bNeighborList)
            {
                if (gpu->sim.ntp > 0)
                    kCalculatePMECHARMMEnergyVirial_kernel<<<gpu->blocks, gpu->CHARMMForcesThreadsPerBlock>>>();
                else
                    kCalculatePMECHARMMEnergy_kernel<<<gpu->blocks, gpu->CHARMMForcesThreadsPerBlock>>>();
            }
            else
                kCalculateCHARMMEnergy_kernel<<<gpu->blocks, gpu->CHARMMForcesThreadsPerBlock>>>();
#ifdef MPI 
            LAUNCHERROR_NONBLOCKING("kCalculateCHARMMEnergy");
#else
            LAUNCHERROR("kCalculateCHARMMEnergy");
#endif 
            cudaUnbindTexture(texref);
        }
    }
}

