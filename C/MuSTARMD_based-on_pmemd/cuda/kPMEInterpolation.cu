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

void SetkPMEInterpolationSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetkPMEInterpolationSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMEGetGridWeights_kernel()
#include "kPGGW.h"

#define PME_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMEGetGridWeightsOrthogonal_kernel()
#include "kPGGW.h"
#undef PME_ORTHOGONAL

#define PME_NTP
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMEGetGridWeightsNTP_kernel()
#include "kPGGW.h"

#define PME_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMEGetGridWeightsOrthogonalNTP_kernel()
#include "kPGGW.h"
#undef PME_ORTHOGONAL
#undef PME_NTP

#define PME_SMALLBOX
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMEGetGridWeightsSmall_kernel()
#include "kPGGW.h"

#define PME_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMEGetGridWeightsSmallOrthogonal_kernel()
#include "kPGGW.h"
#undef PME_ORTHOGONAL

#define PME_NTP
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMEGetGridWeightsSmallNTP_kernel()
#include "kPGGW.h"

#define PME_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMEGetGridWeightsSmallOrthogonalNTP_kernel()
#include "kPGGW.h"
#undef PME_ORTHOGONAL
#undef PME_NTP
#undef PME_SMALLBOX


extern "C" void kPMEGetGridWeights(gpuContext gpu)
{

    if (gpu->bSmallBox)
    {
        if (gpu->sim.ntp > 0)
        {
            if (gpu->sim.is_orthog)    
                kPMEGetGridWeightsSmallOrthogonalNTP_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>();  
            else
                kPMEGetGridWeightsSmallNTP_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>();         
        }
        else
        {
            if (gpu->sim.is_orthog)    
                kPMEGetGridWeightsSmallOrthogonal_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>();  
            else
                kPMEGetGridWeightsSmall_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(); 
        }     
    }
    else 
    {
        if (gpu->sim.ntp > 0)
        {
            if (gpu->sim.is_orthog)    
                kPMEGetGridWeightsOrthogonalNTP_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>();  
            else
                kPMEGetGridWeightsNTP_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>();         
        }
        else
        {
            if (gpu->sim.is_orthog)    
                kPMEGetGridWeightsOrthogonal_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>();  
            else
                kPMEGetGridWeights_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(); 
        } 
    }
#ifdef MPI 
    LAUNCHERROR_NONBLOCKING("kPMEGetGridWeights");
#else
    LAUNCHERROR("kPMEGetGridWeights");
#endif  

#if 0
    cudaThreadSynchronize();
    gpu->pbAtomXYSP->Download();
    gpu->pbAtomZSP->Download();
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        printf("%06d %20.10f %20.10f %20.10f\n", i, gpu->pbAtomXYSP->_pSysData[i].x, gpu->pbAtomXYSP->_pSysData[i].y, gpu->pbAtomZSP->_pSysData[i]);
    }
   // exit(-1);
#endif
}

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMEClearChargeGridBuffer27_kernel()
{
extern __shared__ int sOddBufferOverlapFlag[];
    int* psOddXBufferOverlapFlag                = &sOddBufferOverlapFlag[0];
    int* psOddYBufferOverlapFlag                = &sOddBufferOverlapFlag[cSim.nfft1];
    int* psOddZBufferOverlapFlag                = &sOddBufferOverlapFlag[cSim.nfft1 + cSim.nfft2];


    // Read axis dependent buffer flags
    unsigned int pos                            = threadIdx.x;
    unsigned int end                            = cSim.nfft1 + cSim.nfft2 + cSim.nfft3;
    while (pos < end)
    {
        sOddBufferOverlapFlag[pos]              = cSim.pNLOddBufferOverlapFlag[pos];
        pos                                    += blockDim.x;
    }
    __syncthreads();

    pos                                         = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int limit                          = cSim.nfft1 * cSim.nfft2 * cSim.nfft3;

    while (pos < limit)
    {
        // Calculate buffer count
        int z                                   = pos / cSim.nfft1xnfft2;
        int y                                   = (pos - z * cSim.nfft1xnfft2) / cSim.nfft1;
        int x                                   = pos - z * cSim.nfft1xnfft2 - y * cSim.nfft1;
        int oddBufferCount                      = psOddZBufferOverlapFlag[z] + psOddYBufferOverlapFlag[y] + psOddXBufferOverlapFlag[x];
        int extraBuffers                        = cSim.extraChargeGridBuffers[oddBufferCount]; 
    
        // Clear first 8
        PMEFloat* pFloat                        = &cSim.pXYZ_q[pos];
        *pFloat                                 = (PMEFloat)0.0;
        pFloat                                 += cSim.XYZStride;
        *pFloat                                 = (PMEFloat)0.0;
        pFloat                                 += cSim.XYZStride;        
         *pFloat                                 = (PMEFloat)0.0;
        pFloat                                 += cSim.XYZStride;
        *pFloat                                 = (PMEFloat)0.0;
        pFloat                                 += cSim.XYZStride;
        *pFloat                                 = (PMEFloat)0.0;
        pFloat                                 += cSim.XYZStride;
        *pFloat                                 = (PMEFloat)0.0;
        pFloat                                 += cSim.XYZStride;
        *pFloat                                 = (PMEFloat)0.0;
        pFloat                                 += cSim.XYZStride;
        *pFloat                                 = (PMEFloat)0.0;
        pFloat                                 += cSim.XYZStride;        
        // Clear extra buffers
        while (extraBuffers > 3)
        
       
        {
            *pFloat                             = (PMEFloat)0.0;
            pFloat                             += cSim.XYZStride;
            *pFloat                             = (PMEFloat)0.0;
            pFloat                             += cSim.XYZStride;        
            *pFloat                             = (PMEFloat)0.0;
            pFloat                             += cSim.XYZStride;        
            *pFloat                             = (PMEFloat)0.0;
            pFloat                             += cSim.XYZStride;        
            extraBuffers                       -= 4;
        }       
        
        while (extraBuffers > 0)
        {
            *pFloat                             = (PMEFloat)0.0;
            pFloat                             += cSim.XYZStride;
            extraBuffers--;
        }
 
        pos                                    += blockDim.x * gridDim.x;
    }
}

extern "C" void kPMEClearChargeGridBuffer(gpuContext gpu)
{    
    if (gpu->bOddNLCells)
        kPMEClearChargeGridBuffer27_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock, (gpu->sim.nfft1 + gpu->sim.nfft2 + gpu->sim.nfft3) * sizeof(int)>>>(); 
    else
        cudaMemset(gpu->sim.pXYZ_q, 0, 8 * gpu->sim.XYZStride * sizeof(PMEFloat));
#ifdef MPI
    LAUNCHERROR_NONBLOCKING("kPMEClearChargeGridBuffer");
#else        
    LAUNCHERROR("kPMEClearChargeGridBuffer");
#endif
}

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#endif
kPMEReduceChargeGridBuffer8_kernel()
{
    unsigned int pos                            = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int limit                          = cSim.nfft1 * cSim.nfft2 * cSim.nfft3;

    while (pos < limit)
    {
        PMEFloat* pFloat                        = &cSim.pXYZ_q[pos];
        PMEFloat value1                         = *pFloat;
        pFloat                                 += cSim.XYZStride;
        PMEFloat value2                         = *pFloat;
        pFloat                                 += cSim.XYZStride;
        PMEFloat value3                         = *pFloat;
        pFloat                                 += cSim.XYZStride;
        PMEFloat value4                         = *pFloat;
        pFloat                                 += cSim.XYZStride;
        PMEFloat value5                         = *pFloat;
        pFloat                                 += cSim.XYZStride;
        PMEFloat value6                         = *pFloat;
        pFloat                                 += cSim.XYZStride;
        PMEFloat value7                         = *pFloat;
        pFloat                                 += cSim.XYZStride;
        PMEFloat value8                         = *pFloat;
        cSim.pXYZ_q[pos]                        = value1 + value2 + value3 + value4 + value5 + value6 + value7 + value8;       
        pos                                    += blockDim.x * gridDim.x;
    }
}

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_REDUCEFORCES_THREADS_PER_BLOCK, 1)
#endif
kPMEReduceChargeGridBuffer27_kernel()
{
extern __shared__ int sOddBufferOverlapFlag[];
    int* psOddXBufferOverlapFlag                = &sOddBufferOverlapFlag[0];
    int* psOddYBufferOverlapFlag                = &sOddBufferOverlapFlag[cSim.nfft1];
    int* psOddZBufferOverlapFlag                = &sOddBufferOverlapFlag[cSim.nfft1 + cSim.nfft2];


    // Read axis dependent buffer flags
    unsigned int pos                            = threadIdx.x;
    unsigned int end                            = cSim.nfft1 + cSim.nfft2 + cSim.nfft3;
    while (pos < end)
    {
        sOddBufferOverlapFlag[pos]              = cSim.pNLOddBufferOverlapFlag[pos];
        pos                                    += blockDim.x;
    }
    __syncthreads();

    pos                                         = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int limit                          = cSim.nfft1 * cSim.nfft2 * cSim.nfft3;
    
    while (pos < limit)
    {
        // Calculate buffer count
        int z                                   = pos / cSim.nfft1xnfft2;
        int y                                   = (pos - z * cSim.nfft1xnfft2) / cSim.nfft1;
        int x                                   = pos - z * cSim.nfft1xnfft2 - y * cSim.nfft1;
        int oddBufferCount                      = psOddZBufferOverlapFlag[z] + psOddYBufferOverlapFlag[y] + psOddXBufferOverlapFlag[x];
        int extraBuffers                        = cSim.extraChargeGridBuffers[oddBufferCount]; 
    
        PMEFloat sum                            = (PMEFloat)0.0;
        PMEFloat* pFloat                        = &cSim.pXYZ_q[pos];
        PMEFloat value1                         = *pFloat;
        pFloat                                 += cSim.XYZStride;
        PMEFloat value2                         = *pFloat;
        pFloat                                 += cSim.XYZStride;
        PMEFloat value3                         = *pFloat;
        pFloat                                 += cSim.XYZStride;
        PMEFloat value4                         = *pFloat;
        pFloat                                 += cSim.XYZStride;
        PMEFloat value5                         = *pFloat;
        pFloat                                 += cSim.XYZStride;
        PMEFloat value6                         = *pFloat;
        pFloat                                 += cSim.XYZStride;
        PMEFloat value7                         = *pFloat;
        pFloat                                 += cSim.XYZStride;
        PMEFloat value8                         = *pFloat;
        pFloat                                 += cSim.XYZStride;
        
        while (extraBuffers > 3)
        {
            PMEFloat val1                       = *pFloat;
            pFloat                             += cSim.XYZStride;
            PMEFloat val2                       = *pFloat;
            pFloat                             += cSim.XYZStride;
            PMEFloat val3                       = *pFloat;
            pFloat                             += cSim.XYZStride;
            PMEFloat val4                       = *pFloat;
            pFloat                             += cSim.XYZStride;
            sum                                += val1 + val2 + val3 + val4; 
            extraBuffers                       -= 4;
        }
        
        while (extraBuffers > 0)
        {
            PMEFloat val                        = *pFloat;
            pFloat                             += cSim.XYZStride;
            sum                                += val;
            extraBuffers--;
        }        
        cSim.pXYZ_q[pos]                        = value1 + value2 + value3 + value4 + value5 + value6 + value7 + value8 + sum;
   
        pos                                    += blockDim.x * gridDim.x;
    }
}

extern "C" void kPMEReduceChargeGridBuffer(gpuContext gpu)
{



#if 0
    int* pOddXBufferOverlapFlag                 = &(gpu->pbNLOddBufferOverlapFlag->_pSysData[0]);
    int* pOddYBufferOverlapFlag                 = &(gpu->pbNLOddBufferOverlapFlag->_pSysData[gpu->sim.nfft1]);
    int* pOddZBufferOverlapFlag                 = &(gpu->pbNLOddBufferOverlapFlag->_pSysData[gpu->sim.nfft1 + gpu->sim.nfft2]);
#if 0   
    for (int i = 0; i < gpu->sim.nfft1; i++)
    {
        printf("%d %d %d %d\n", i, pOddXBufferOverlapFlag[i], pOddYBufferOverlapFlag[i], pOddZBufferOverlapFlag[i]);
    }
#endif    
    gpu->pbXYZ_q->Download();
    float* pXYZ_q = gpu->pbXYZ_q->_pSysData;
    for (int i = 0; i < gpu->sim.nfft1 * gpu->sim.nfft2 * gpu->sim.nfft3; i++)
    {
        int z                                   = i / gpu->sim.nfft1xnfft2;
        int y                                   = (i - z * gpu->sim.nfft1xnfft2) / gpu->sim.nfft1;
        int x                                   = i - z * gpu->sim.nfft1xnfft2 - y * gpu->sim.nfft1; 
        int oddBufferCount                      = pOddZBufferOverlapFlag[z] + pOddYBufferOverlapFlag[y] + pOddXBufferOverlapFlag[x];
        int buffers                             = gpu->sim.extraChargeGridBuffers[oddBufferCount] + 8; 
        
        
        int maxBuffer = 0;
        for (int j = 0; j < 27; j++)
            if (pXYZ_q[i + j * gpu->sim.XYZStride] > 0.0f)
                maxBuffer = j;
        if (maxBuffer > buffers)
            printf("%3d %3d %3d %3d %3d %3d %3d %3d\n", x, y, z, pOddXBufferOverlapFlag[x], pOddYBufferOverlapFlag[y], pOddZBufferOverlapFlag[z], maxBuffer, buffers);
        
    }
   // exit(-1);

#endif

    if (gpu->bOddNLCells)
        kPMEReduceChargeGridBuffer27_kernel<<<gpu->blocks, gpu->reduceForcesThreadsPerBlock, (gpu->sim.nfft1 + gpu->sim.nfft2 + gpu->sim.nfft3) * sizeof(int)>>>();
    else
        kPMEReduceChargeGridBuffer8_kernel<<<gpu->blocks, gpu->reduceForcesThreadsPerBlock>>>();   
    LAUNCHERROR("kPMEReduceChargeGridBuffer");


#if 0
    gpu->pbXYZ_qc->Download();
    cufftComplex* xyz_qc = gpu->pbXYZ_qc->_pSysData;
  
    float sum = 0.0;
    for (int i = 0; i < gpu->sim.nfft1; i++)
    {
        for (int j = 0; j < gpu->sim.nfft2; j++)
        {
            for (int k = 0; k < gpu->sim.nfft3; k++)
            {
                printf("%3d %3d %3d %32.15f\n", i, j, k, xyz_qc[(k * gpu->sim.nfft2 + j) * gpu->sim.nfft1 + i].x); 
                sum += xyz_qc[(k * gpu->sim.nfft2 + j) * gpu->sim.nfft1 + i].x;
            }
        }
    }
    printf("%f\n", sum);
    printf("%d %6.1f\n", gpu->sim.atoms, sum / gpu->sim.atoms);
    exit(-1);
#endif

#if 0
    printf("%d %d %d %d\n", gpu->sim.xcells, gpu->sim.ycells, gpu->sim.zcells, gpu->sim.xycells);
    for (int pos = 0; pos < gpu->sim.cells; pos++)
    {
        int zCell                               = pos / gpu->sim.xycells;
        int yCell                               = (pos - zCell * gpu->sim.xycells) / gpu->sim.xcells;
        int xCell                               = pos - zCell * gpu->sim.xycells - yCell * gpu->sim.xcells; 
        printf("%d %d %d %d %d\n", pos, xCell, yCell, zCell, gpu->pbNLChargeGridBufferOffset->_pSysData[pos] / gpu->sim.XYZStride);
    }
    
    printf("%d\n", gpu->sim.XYZStride);
    printf("%d\n", gpu->sim.nfft1xnfft2);
    exit(-1);
#endif
}

struct FillChargeGridAtomData
{
    int ix;
    int iy;
    int iz;
    PMEFloat tx[4];
    PMEFloat ty[4];
    PMEFloat tz[4];
};
static const int LOADSIZE = 32;
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(64, 8)
#else
__launch_bounds__(64, 8)
#endif
kPMEFillChargeGridBuffer8_kernel()
{
__shared__ volatile FillChargeGridAtomData sAtom[LOADSIZE];

    // Determine grid offsets   
    const int tOffsetX                          = threadIdx.x & 0x03; 
    const int tOffsetY                          = (threadIdx.x & 0x0f) >> 2;
    const int tOffsetZ                          = threadIdx.x >> 4;
    const int iOffsetX                          = tOffsetX;
    const int iOffsetY                          = tOffsetY;
    const int iOffsetZ                          = tOffsetZ;
        
        
    // Load cell data
    if (threadIdx.x == 0)
    {
        uint2 cell                              = cSim.pNLNonbondCellStartEnd[blockIdx.x];
        sAtom[0].ix                             = cell.x;
        sAtom[0].iy                             = cell.y;
    }
    __syncthreads();
    uint2 cellStartEnd;
    cellStartEnd.x                              = sAtom[0].ix;
    cellStartEnd.y                              = sAtom[0].iy;
    __syncthreads();
        
    // Calculate buffer offset
    int zCell                                   = blockIdx.x / cSim.xycells;
    int yCell                                   = (blockIdx.x - zCell * cSim.xycells) / cSim.xcells;
    int xCell                                   = blockIdx.x - zCell * cSim.xycells - yCell * cSim.xcells; 
    int bufferOffset                            = (4 * (zCell & 0x1) + 2 * (yCell & 0x1) + (xCell & 0x1)) * cSim.XYZStride;
    PMEFloat* pXYZ_q                            = cSim.pXYZ_q + bufferOffset;

    // Iterate through cell
    unsigned int pos                            = cellStartEnd.x;
    while (pos  < cellStartEnd.y)
    {
        
        // Read Atom Data
        unsigned int maxatom                    = min(pos + LOADSIZE, cellStartEnd.y);
        unsigned int pos1                       = pos + threadIdx.x;
        if (pos1 < maxatom)
        {
            PMEFloat charge                     = cSim.pAtomChargeSP[pos1];
            int ix                              = cSim.pIFractX[pos1];
            int iy                              = cSim.pIFractY[pos1];
            int iz                              = cSim.pIFractZ[pos1];
            PMEFloat4 tx                        = cSim.pThetaX[pos1];
            PMEFloat4 ty                        = cSim.pThetaY[pos1];
            PMEFloat4 tz                        = cSim.pThetaZ[pos1];
            sAtom[threadIdx.x].ix               = ix;
            sAtom[threadIdx.x].iy               = iy;
            sAtom[threadIdx.x].iz               = iz;
            sAtom[threadIdx.x].tx[0]            = tx.x * charge;
            sAtom[threadIdx.x].tx[1]            = tx.y * charge;
            sAtom[threadIdx.x].tx[2]            = tx.z * charge;
            sAtom[threadIdx.x].tx[3]            = tx.w * charge;
            sAtom[threadIdx.x].ty[0]            = ty.x;
            sAtom[threadIdx.x].ty[1]            = ty.y;
            sAtom[threadIdx.x].ty[2]            = ty.z;
            sAtom[threadIdx.x].ty[3]            = ty.w;
            sAtom[threadIdx.x].tz[0]            = tz.x;
            sAtom[threadIdx.x].tz[1]            = tz.y;
            sAtom[threadIdx.x].tz[2]            = tz.z;
            sAtom[threadIdx.x].tz[3]            = tz.w;
        }
        __syncthreads();
      
        // Interpolate onto grid
        pos1                                    = 0;
        unsigned int lastAtom                   = min(LOADSIZE, cellStartEnd.y - pos);
        while (pos1 < lastAtom)
        {
            // Calculate values
            int ix                              = sAtom[pos1].ix + iOffsetX;
            int iy                              = sAtom[pos1].iy + iOffsetY;
            int iz                              = sAtom[pos1].iz + iOffsetZ;
            
            // Insure coordinates stay in bounds
            if (ix >= cSim.nfft1)
                ix                             -= cSim.nfft1;
            if (iy >= cSim.nfft2)
                iy                             -= cSim.nfft2;
            if (iz >= cSim.nfft3)
                iz                             -= cSim.nfft3;  
                
            // Calculate interpolation values and destinations    
            int gpos                            = (iz * cSim.nfft2 + iy) * cSim.nfft1 + ix; 
            PMEFloat qvalue                     = pXYZ_q[gpos];        
            PMEFloat value                      = sAtom[pos1].tx[tOffsetX] * sAtom[pos1].ty[tOffsetY] * sAtom[pos1].tz[tOffsetZ];
                    
            // Write memory and sync all threads          
            pXYZ_q[gpos]                        = qvalue + value;
            __threadfence_block();
            __syncthreads();
 
            pos1++;
         }
         pos                                   += LOADSIZE;
    }      
}

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(64, 8)
#else
__launch_bounds__(64, 8)
#endif
kPMEFillChargeGridBuffer27_kernel()
{
__shared__ volatile FillChargeGridAtomData sAtom[LOADSIZE];

    // Determine grid offsets   
    const int tOffsetX                          = threadIdx.x & 0x03; 
    const int tOffsetY                          = (threadIdx.x & 0x0f) >> 2;
    const int tOffsetZ                          = threadIdx.x >> 4;
    const int iOffsetX                          = tOffsetX;
    const int iOffsetY                          = tOffsetY;
    const int iOffsetZ                          = tOffsetZ;
        
        
    // Load cell data
    if (threadIdx.x == 0)
    {
        uint2 cell                              = cSim.pNLNonbondCellStartEnd[blockIdx.x];
        int bufferOffset                        = cSim.pNLChargeGridBufferOffset[blockIdx.x];
        sAtom[0].ix                             = cell.x;
        sAtom[0].iy                             = cell.y;
        sAtom[0].iz                             = bufferOffset;
    }
    __syncthreads();
    PMEFloat* pXYZ_q                            = cSim.pXYZ_q + sAtom[0].iz;
    uint2 cellStartEnd;
    cellStartEnd.x                              = sAtom[0].ix;
    cellStartEnd.y                              = sAtom[0].iy;        
    __syncthreads();
        
        
    // Iterate through cell
    unsigned int pos                            = cellStartEnd.x;
    while (pos  < cellStartEnd.y)
    {
        
        // Read Atom Data
        unsigned int maxatom                    = min(pos + LOADSIZE, cellStartEnd.y);
        unsigned int pos1                       = pos + threadIdx.x;
        if (pos1 < maxatom)
        {
            PMEFloat charge                     = cSim.pAtomChargeSP[pos1];
            int ix                              = cSim.pIFractX[pos1];
            int iy                              = cSim.pIFractY[pos1];
            int iz                              = cSim.pIFractZ[pos1];
            PMEFloat4 tx                        = cSim.pThetaX[pos1];
            PMEFloat4 ty                        = cSim.pThetaY[pos1];
            PMEFloat4 tz                        = cSim.pThetaZ[pos1];
            sAtom[threadIdx.x].ix               = ix;
            sAtom[threadIdx.x].iy               = iy;
            sAtom[threadIdx.x].iz               = iz;
            sAtom[threadIdx.x].tx[0]            = tx.x * charge;
            sAtom[threadIdx.x].tx[1]            = tx.y * charge;
            sAtom[threadIdx.x].tx[2]            = tx.z * charge;
            sAtom[threadIdx.x].tx[3]            = tx.w * charge;
            sAtom[threadIdx.x].ty[0]            = ty.x;
            sAtom[threadIdx.x].ty[1]            = ty.y;
            sAtom[threadIdx.x].ty[2]            = ty.z;
            sAtom[threadIdx.x].ty[3]            = ty.w;
            sAtom[threadIdx.x].tz[0]            = tz.x;
            sAtom[threadIdx.x].tz[1]            = tz.y;
            sAtom[threadIdx.x].tz[2]            = tz.z;
            sAtom[threadIdx.x].tz[3]            = tz.w;
        }
        __syncthreads();
      
        // Interpolate onto grid
        pos1                                    = 0;
        unsigned int lastAtom                   = min(LOADSIZE, cellStartEnd.y - pos);
        while (pos1 < lastAtom)
        {
            // Calculate values
            int ix                              = sAtom[pos1].ix + iOffsetX;
            int iy                              = sAtom[pos1].iy + iOffsetY;
            int iz                              = sAtom[pos1].iz + iOffsetZ;
            
            // Insure coordinates stay in bounds
            if (ix >= cSim.nfft1)
                ix                             -= cSim.nfft1;
            if (iy >= cSim.nfft2)
                iy                             -= cSim.nfft2;
            if (iz >= cSim.nfft3)
                iz                             -= cSim.nfft3;  
                
            // Calculate interpolation values and destinations    
            int gpos                            = (iz * cSim.nfft2 + iy) * cSim.nfft1 + ix; 
            PMEFloat qvalue                     = pXYZ_q[gpos];        
            PMEFloat value                      = sAtom[pos1].tx[tOffsetX] * sAtom[pos1].ty[tOffsetY] * sAtom[pos1].tz[tOffsetZ];
                    
            // Write memory and sync all threads          
            pXYZ_q[gpos]                        = qvalue + value;
            __threadfence_block();
            __syncthreads();
 
            pos1++;
         }
         pos                                   += LOADSIZE;
    }
}

extern "C" void PMEInitKernels(gpuContext gpu)
{
	if (gpu->sm_version >= SM_2X)
	{
        cudaFuncSetCacheConfig(kPMEFillChargeGridBuffer27_kernel, cudaFuncCachePreferL1);
        cudaFuncSetCacheConfig(kPMEFillChargeGridBuffer8_kernel, cudaFuncCachePreferL1);       
    }
}


extern "C" void kPMEFillChargeGridBuffer(gpuContext gpu)
{
	if (gpu->sm_version >= SM_2X)
	{ 
	    if (gpu->bOddNLCells)
    	    kPMEFillChargeGridBuffer27_kernel<<<gpu->sim.cells, 64>>>(); 
   	 	else
        	kPMEFillChargeGridBuffer8_kernel<<<gpu->sim.cells, 64>>>();
    }
	else {
 		if (gpu->bOddNLCells)
        	kPMEFillChargeGridBuffer27_kernel<<<gpu->sim.cells, 64>>>(); 
    	else
        	kPMEFillChargeGridBuffer8_kernel<<<gpu->sim.cells, 64>>>();  
	}
#ifdef MPI
    LAUNCHERROR_NONBLOCKING("kPMEFillChargeGridBuffer");
#else	
    LAUNCHERROR("kPMEFillChargeGridBuffer");
#endif
#if 0  
    gpu->pbNLNonbondCellStartEnd->Download();
    gpu->pbIFract->Download();
    uint2* pNLCellStartEnd = gpu->pbNLNonbondCellStartEnd->_pSysData;
    int* pIFractX = gpu->pbIFract->_pSysData;
    int* pIFractY = pIFractX + gpu->sim.stride;
    int* pIFractZ = pIFractX + gpu->sim.stride2;
    
    for (int i = 0; i < gpu->sim.cells; i++)
    {
        int start = pNLCellStartEnd[i].x;
        int end = pNLCellStartEnd[i].y;
        int zCell                               = i / gpu->sim.xycells;
        int yCell                               = (i - zCell * gpu->sim.xycells) / gpu->sim.xcells;
        int xCell                               = i - zCell * gpu->sim.xycells - yCell * gpu->sim.xcells; 
        //int bufferOffset                        = (4 * (zCell & 0x1) + 2 * (yCell & 0x1) + (xCell & 0x1));
        int xmin                                = xCell       * gpu->sim.nfft1 / gpu->sim.xcells - 4;
        if (xmin < 0)
            xmin += gpu->sim.nfft1;
        int xmax                                = (xCell + 1) * gpu->sim.nfft1 / gpu->sim.xcells;
        int ymin                                = yCell       * gpu->sim.nfft2 / gpu->sim.ycells - 4;
        if (ymin < 0)
            ymin += gpu->sim.nfft2;
        int ymax                                = (yCell + 1) * gpu->sim.nfft2 / gpu->sim.ycells;
        int zmin                                = zCell       * gpu->sim.nfft3 / gpu->sim.zcells - 4;
        if (zmin < 0)
            zmin += gpu->sim.nfft3;
        int zmax                                = (zCell + 1) * gpu->sim.nfft3 / gpu->sim.zcells;
        printf("C: %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d\n", i, xCell, yCell, zCell, xmin, xmax, ymin, ymax, zmin, zmax);
        
        for (int j = start; j < end; j++)
        {
            int ix = pIFractX[j];
            int iy = pIFractY[j];
            int iz = pIFractZ[j];
            bool valid = true;
            if (xmin < xmax)
            {
                if ((ix < xmin) | (ix > xmax))
                    valid = false;
            }
            else
            {
                if ((ix > xmax) && (ix < xmin))
                    valid = false; 
            }
            
            if (ymin < ymax)
            {
                if ((iy < ymin) | (iy > ymax))
                    valid = false;
            }
            else
            {
                if ((iy > ymax) && (iy < ymin))
                    valid = false; 
            }
            
            if (zmin < zmax)
            {
                if ((iz < zmin) | (iz > zmax))
                    valid = false;
            }
            else
            {
                if ((iz > zmax) && (iz < zmin))
                    valid = false; 
            }
            if (!valid)
            {
                printf("F: %6d %3d %3d %3d\n", j , ix, iy, iz);
            }
           // printf("A: %6d %3d %3d %3d\n", j , ix, iy, iz);
        }
        
        
        
    }
    exit(-1);
#endif    
    
    
#if 0    
    for (int i = 0; i < gpu->sim.xcells; i++)
    {
        int ixstart = i * gpu->sim.nfft1 / gpu->sim.xcells - 4;
        int ixend = (i + 1) * gpu->sim.nfft1 / gpu->sim.xcells;
        printf("%d %d %d\n", i, ixstart, ixend);
    }
    exit(-1);
#endif
}
 
#define PME_ENERGY
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMEScalarSumRCEnergy_kernel(PMEDouble ewaldcof, PMEDouble vol)
#include "kPSSE.h"
#undef PME_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMEScalarSumRC_kernel(PMEDouble ewaldcof, PMEDouble vol)
#include "kPSSE.h"

#define PME_VIRIAL
#define PME_ENERGY
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMEScalarSumRCEnergyVirial_kernel(PMEDouble ewaldcof, PMEDouble vol)
#include "kPSSE.h"
#undef PME_ENERGY

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kPMEScalarSumRCVirial_kernel(PMEDouble ewaldcof, PMEDouble vol)
#include "kPSSE.h"
#undef PME_VIRIAL

extern "C" void kPMEScalarSumRC(gpuContext gpu, PMEDouble ewaldcof, PMEDouble vol)
{
    if (gpu->sim.ntp > 0)
        kPMEScalarSumRCVirial_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(ewaldcof, vol);     
    else
        kPMEScalarSumRC_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(ewaldcof, vol);   
#ifdef MPI
    LAUNCHERROR_NONBLOCKING("kPMEScalarSumRC");
#else          
    LAUNCHERROR("kPMEScalarSumRC");
#endif
}

extern "C" void kPMEScalarSumRCEnergy(gpuContext gpu, PMEDouble ewaldcof, PMEDouble vol)
{
    if (gpu->sim.ntp > 0)
        kPMEScalarSumRCEnergyVirial_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(ewaldcof, vol);
    else
        kPMEScalarSumRCEnergy_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(ewaldcof, vol);
#ifdef MPI
    LAUNCHERROR_NONBLOCKING("kPMEScalarSumRCEnergy");
#else        
    LAUNCHERROR("kPMEScalarSumRCEnergy");
#endif
}

#ifdef use_DPDP
texture<int2, 1, cudaReadModeElementType> texref;
#else
texture<PMEFloat, 1, cudaReadModeElementType> texref;
#endif
static const int GRADSUMTHREADS = 64;

__global__ void 
#if (__CUDA_ARCH__ >= 200)
#ifdef use_DPDP
__launch_bounds__(GRADSUMTHREADS, 8)
#else
__launch_bounds__(GRADSUMTHREADS, 8)
#endif
#else
#ifdef use_DPDP
__launch_bounds__(GRADSUMTHREADS, 3)
#else
__launch_bounds__(GRADSUMTHREADS, 4)
#endif
#endif
kPMEGradSum64_kernel()
#include "kPGS.h"

#define PME_VIRIAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
#ifdef use_DPDP
__launch_bounds__(GRADSUMTHREADS, 8)
#else
__launch_bounds__(GRADSUMTHREADS, 8)
#endif
#else
#ifdef use_DPDP
__launch_bounds__(GRADSUMTHREADS, 3)
#else
__launch_bounds__(GRADSUMTHREADS, 4)
#endif
#endif
kPMEGradSum64Virial_kernel()
#include "kPGS.h"
#undef PME_VIRIAL


extern "C" void kPMEGradSum(gpuContext gpu)
{
    texref.normalized = 0;
    texref.filterMode = cudaFilterModePoint;
    texref.addressMode[0] = cudaAddressModeClamp;
    texref.channelDesc.x = 32;
#ifdef use_DPDP    
    texref.channelDesc.y = 32;
#else
    texref.channelDesc.y = 0;
#endif
    texref.channelDesc.z = 0;
    texref.channelDesc.w = 0;
#ifdef use_DPDP
    cudaBindTexture(NULL, texref, (int2*)(gpu->sim.pXYZ_q), gpu->sim.nfft1 * gpu->sim.nfft2 * gpu->sim.nfft3 * sizeof(int2));
#else
    cudaBindTexture(NULL, texref, (PMEFloat*)(gpu->sim.pXYZ_q), gpu->sim.nfft1 * gpu->sim.nfft2 * gpu->sim.nfft3 * sizeof(PMEFloat));
#endif

    int blocks;
    if (gpu->sm_version >= SM_2X)
        blocks = (gpu->sim.atoms + 32 - 1) / 32;
    else
        blocks = (gpu->sim.atoms + 16 - 1) / 16; 
    if (gpu->sim.ntp > 0)
        kPMEGradSum64Virial_kernel<<<blocks, GRADSUMTHREADS>>>();   
    else
        kPMEGradSum64_kernel<<<blocks, GRADSUMTHREADS>>>();   

    LAUNCHERROR("kPMEGradSum");
    cudaUnbindTexture(texref);
}


