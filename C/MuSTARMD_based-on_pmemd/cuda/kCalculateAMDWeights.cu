/***************************************************/
/*                                                 */
/*      AMBER NVIDIA CUDA CPU IMPLEMENTATION       */
/*                 PMEMD VERSION                   */
/*                     2011                        */
/*                      by                         */
/*                Romelia Salomon (SDSC)           */
/*                                                 */
/***************************************************/

#include <cuda.h>
#include "gpu.h"

static __constant__ cudaSimulation cSim;
static __constant__ PMEDouble tenm3             = (PMEDouble)(1.0e-03);
static __constant__ PMEDouble tm24              = (PMEDouble)(1.0e-18);
static __constant__ PMEDouble one               = (PMEDouble)(1.0);
static __constant__ PMEDouble zero              = (PMEDouble)(0.0);

// Texture reference for PMEDouble-precision coordinates (disguised as int2 to work around HW limitations)
#ifndef use_SPSP
texture<int2, 1, cudaReadModeElementType> texref;
#else
texture<float, 1, cudaReadModeElementType> texref;
#endif

void SetkCalculateAMDWeightsSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetkCalculateAMDWeightssSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}


__device__ void faster_sincos2(double a, double *sptr, double *cptr) 
#include "kFastCosSin.h"


__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kAMDCalcWeightAndScaleForces_kernel(PMEDouble pot_ene_tot, PMEDouble dih_ene_tot, PMEDouble fwgt)
{

  //calculate AMD weight, seting dihedral boost (tboost) to zero for now

    unsigned int pos                        = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int increment                              = gridDim.x * blockDim.x;    
    while (pos < cSim.atoms)
    {
        PMEDouble forceX                            = cSim.pForceX[pos];
        PMEDouble forceY                            = cSim.pForceY[pos];
        PMEDouble forceZ                            = cSim.pForceZ[pos];

        forceX                                 *= fwgt;
        forceY                                 *= fwgt;
        forceZ                                 *= fwgt;

        cSim.pForceX[pos]                = forceX;
        cSim.pForceY[pos]                = forceY;
        cSim.pForceZ[pos]                = forceZ;
        pos                                += increment;       
    }
}

void kAMDCalcWeightAndScaleForces(gpuContext gpu, PMEDouble pot_ene_tot, PMEDouble dih_ene_tot, PMEDouble fwgt)
{

      kAMDCalcWeightAndScaleForces_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(pot_ene_tot,dih_ene_tot,fwgt);

    LAUNCHERROR("kAMDScaleForces");

}

#define LOCAL_ENERGY
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculateAmdDihedralEnergy_kernel()
#include "kCLFdih.h"
#undef LOCAL_ENERGY

#define LOCAL_NEIGHBORLIST
#define LOCAL_ENERGY
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMEAmdDihedralEnergy_kernel()
#include "kCLFdih.h"
#undef LOCAL_ENERGY
#undef LOCAL_NEIGHBORLIST

// Consumer Fermi kernels
#define NODPTEXTURE
#define LOCAL_ENERGY
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculateAmdDihedralEnergyFermi_kernel()
#include "kCLFdih.h"
#undef LOCAL_ENERGY

#define LOCAL_NEIGHBORLIST
#define LOCAL_ENERGY
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
kCalculatePMEAmdDihedralEnergyFermi_kernel()
#include "kCLFdih.h"
#undef LOCAL_ENERGY
#undef LOCAL_NEIGHBORLIST
#undef NODPTEXTURE


extern "C" void kCalculateAmdDihedralEnergy(gpuContext gpu)
{

//Calculate dihedral energy
    if (gpu->bLocalInteractions)
    {
#ifndef use_SPSP    
        if (!gpu->bECCSupport && (gpu->sm_version == SM_2X))
        {
            if (gpu->bNeighborList)
            {
                    kCalculatePMEAmdDihedralEnergyFermi_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
            }
            else
                kCalculateAmdDihedralEnergyFermi_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
#ifdef MPI 
            LAUNCHERROR_NONBLOCKING("kCalculateAmdDihedralEnergyFermi");
#else
            LAUNCHERROR("kCalculateAmdDihedralEnergyFermi");
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
                kCalculatePMEAmdDihedralEnergy_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
            }
            else
                kCalculateAmdDihedralEnergy_kernel<<<gpu->blocks, gpu->localForcesThreadsPerBlock>>>();
#ifdef MPI 
            LAUNCHERROR_NONBLOCKING("kCalculateAmdDihedralEnergy");
#else
            LAUNCHERROR("kCalculateAmdDihedralEnergy");
#endif  
            cudaUnbindTexture(texref);
        }
    }

}

