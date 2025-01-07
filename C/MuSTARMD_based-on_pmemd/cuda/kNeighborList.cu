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

struct Atom {
    PMEFloat xmin;
    PMEFloat xmax;
    PMEFloat ymin;
    PMEFloat ymax;
    PMEFloat zmin;
    PMEFloat zmax;
};

void SetkNeighborListSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetkNeighborListSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

struct stuff
{
    double fx;
    double fy;
    double fz;
    unsigned int ix;
    unsigned int iy;
    unsigned int iz;
    unsigned int ox;
    unsigned int oy;
    unsigned int oz;
    unsigned int cellHash;
    unsigned int hash;
};

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kNLGenerateSpatialHash_kernel()
{
    __shared__ unsigned int sCellHash[CELLHASHCELLS];
    __shared__ PMEFloat sRecipf[9];
    unsigned int pos                            = blockIdx.x * blockDim.x + threadIdx.x; 
    unsigned int increment                      = gridDim.x * blockDim.x;
    
    // Clear atom list/exclusion mask space counter
    if (pos == 0)
    {
        *(cSim.pNLTotalOffset)                  = 0;
        *(cSim.pNLPosition)                     = cSim.NLBuildWarps;
    }
    
    // Read cell hash
    if (threadIdx.x < CELLHASHCELLS)
    {
        sCellHash[threadIdx.x]                  = cSim.pNLCellHash[threadIdx.x];
    }

    if (cSim.ntp > 0)
    {
        if (threadIdx.x < 9)
        {
            sRecipf[threadIdx.x]                = cSim.pNTPData->recipf[threadIdx.x];
        }    
         __syncthreads();
    
        while (pos < cSim.atoms)
        {
            PMEFloat x                          = cSim.pImageX[pos];
            PMEFloat y                          = cSim.pImageY[pos];
            PMEFloat z                          = cSim.pImageZ[pos];
            
            // Orthogonal/nonorthogonal handled in the same code (3 single precision multiplies and adds? Who cares and why?)
            PMEFloat fx                         = sRecipf[0] * x + sRecipf[3] * y + sRecipf[6] * z;
            PMEFloat fy                         =                  sRecipf[4] * y + sRecipf[7] * z;
            PMEFloat fz                         =                                   sRecipf[8] * z;

            // Account for minimum image convention  
            fx                                  = fx - round(fx) + (PMEFloat)0.5;
            fy                                  = fy - round(fy) + (PMEFloat)0.5;
            fz                                  = fz - round(fz) + (PMEFloat)0.5;
            fx                                  = (fx < (PMEFloat)1.0 ? fx : (PMEFloat)0.0);
            fy                                  = (fy < (PMEFloat)1.0 ? fy : (PMEFloat)0.0);
            fz                                  = (fz < (PMEFloat)1.0 ? fz : (PMEFloat)0.0);
                
            // Generate box coordinates
            cSim.pImageIndex[pos]               = pos;
            unsigned int ix                     = fx * cSim.xcells;
            unsigned int iy                     = fy * cSim.ycells;
            unsigned int iz                     = fz * cSim.zcells;
            cSim.pImageCellID[pos]              = ix + (iy << CELLIDYSHIFT) + (iz << CELLIDZSHIFT);
            unsigned int ox                     = min(CELLHASHX - 1, (unsigned int)((PMEFloat)CELLHASHX * ((fx - ix * cSim.oneOverXcellsf) * cSim.xcells))); 
            unsigned int oy                     = min(CELLHASHY - 1, (unsigned int)((PMEFloat)CELLHASHY * ((fy - iy * cSim.oneOverYcellsf) * cSim.ycells))) * CELLHASHX; 
            unsigned int oz                     = min(CELLHASHZ - 1, (unsigned int)((PMEFloat)CELLHASHZ * ((fz - iz * cSim.oneOverZcellsf) * cSim.zcells))) * CELLHASHXY; 
            unsigned int cellHash               = sCellHash[ox + oy + oz];
            unsigned int hash                   = (((iz * cSim.ycells + iy) * cSim.xcells + ix) << CELLHASHBITS) | cellHash;
            cSim.pImageHash[pos]                = hash;
            pos                                += increment;
        }   
    }
    else
    {
        __syncthreads();
        while (pos < cSim.atoms)
        {
            PMEFloat x                          = cSim.pImageX[pos];
            PMEFloat y                          = cSim.pImageY[pos];
            PMEFloat z                          = cSim.pImageZ[pos];
            
            // Orthogonal/nonorthogonal handled in the same code (3 single precision multiplies and adds? Who cares and why?)
            PMEFloat fx                         = cSim.recipf[0][0] * x + cSim.recipf[1][0] * y + cSim.recipf[2][0] * z;
            PMEFloat fy                         =                         cSim.recipf[1][1] * y + cSim.recipf[2][1] * z;
            PMEFloat fz                         =                                                 cSim.recipf[2][2] * z;

            // Account for minimum image convention  
            fx                                  = fx - round(fx) + (PMEFloat)0.5;
            fy                                  = fy - round(fy) + (PMEFloat)0.5;
            fz                                  = fz - round(fz) + (PMEFloat)0.5;
            fx                                  = (fx < (PMEFloat)1.0 ? fx : (PMEFloat)0.0);
            fy                                  = (fy < (PMEFloat)1.0 ? fy : (PMEFloat)0.0);
            fz                                  = (fz < (PMEFloat)1.0 ? fz : (PMEFloat)0.0);
                
            // Generate box coordinates
            cSim.pImageIndex[pos]               = pos;
            unsigned int ix                     = fx * cSim.xcells;
            unsigned int iy                     = fy * cSim.ycells;
            unsigned int iz                     = fz * cSim.zcells;
            cSim.pImageCellID[pos]              = ix + (iy << CELLIDYSHIFT) + (iz << CELLIDZSHIFT);
            unsigned int ox                     = min(CELLHASHX - 1, (unsigned int)((PMEFloat)CELLHASHX * ((fx - ix * cSim.oneOverXcellsf) * cSim.xcells))); 
            unsigned int oy                     = min(CELLHASHY - 1, (unsigned int)((PMEFloat)CELLHASHY * ((fy - iy * cSim.oneOverYcellsf) * cSim.ycells))) * CELLHASHX; 
            unsigned int oz                     = min(CELLHASHZ - 1, (unsigned int)((PMEFloat)CELLHASHZ * ((fz - iz * cSim.oneOverZcellsf) * cSim.zcells))) * CELLHASHXY; 
            unsigned int cellHash               = sCellHash[ox + oy + oz];
            unsigned int hash                   = (((iz * cSim.ycells + iy) * cSim.xcells + ix) << CELLHASHBITS) | cellHash;
            cSim.pImageHash[pos]                = hash;
            pos                                += increment;
        }   
    }
    
    
    
}

extern "C" void kNLGenerateSpatialHash(gpuContext gpu)
{
    kNLGenerateSpatialHash_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();   
    LAUNCHERROR("kNLGenerateSpatialHash");
    
#if 0
    gpu->pbFraction->Download();
    gpu->pbImageIndex->Download();
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        PMEFloat fx = gpu->pbFraction->_pSysData[i];
        PMEFloat fy = gpu->pbFraction->_pSysData[i + gpu->sim.stride];
        PMEFloat fz = gpu->pbFraction->_pSysData[i + gpu->sim.stride2];
        unsigned int ix                         = fx * gpu->sim.xcells;
        unsigned int iy                         = fy * gpu->sim.ycells;
        unsigned int iz                         = fz * gpu->sim.zcells;
        unsigned int ox                         = (unsigned int)((PMEFloat)CELLHASHX * ((fx - ix * gpu->sim.oneOverXcells) * gpu->sim.xcells)); 
        unsigned int oy                         = (unsigned int)((PMEFloat)CELLHASHY * ((fy - iy * gpu->sim.oneOverYcells) * gpu->sim.ycells)) * CELLHASHX; 
        unsigned int oz                         = (unsigned int)((PMEFloat)CELLHASHZ * ((fz - iz * gpu->sim.oneOverZcells) * gpu->sim.zcells)) * CELLHASHXY; 
    
    
        printf("%6d 0x%08x %10.7f %10.7f %10.7f |  %6d %6d %6d | %6d %6d %6d\n", i, gpu->pbImageIndex->_pSysData[i + gpu->sim.stride3], fx, fy, fz, ix, iy, iz, ox, oy, oz); 
    }
    exit(-1);
#endif
    
    
}



__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kNLRemapImage_kernel(unsigned int* pImageIndex)
{
    unsigned int pos                            = blockIdx.x * blockDim.x + threadIdx.x; 
    unsigned int increment                      = gridDim.x * blockDim.x;
    unsigned int index                          = 0;
    unsigned int newindex;
    
    if (pos < cSim.atoms)
        index                                   = pImageIndex[pos];

    while (pos < cSim.atoms)
    {
        unsigned int newpos                     = pos + increment;
        if (newpos < cSim.atoms)
            newindex                            = pImageIndex[newpos];
      
        // Read new data
        unsigned int atom                       = cSim.pImageAtom[index];
        PMEDouble x                             = cSim.pImageX[index];
        PMEDouble y                             = cSim.pImageY[index];
        PMEDouble z                             = cSim.pImageZ[index];
        PMEDouble vx                            = cSim.pImageVelX[index];
        PMEDouble vy                            = cSim.pImageVelY[index];
        PMEDouble vz                            = cSim.pImageVelZ[index];
        PMEDouble lvx                           = cSim.pImageLVelX[index];
        PMEDouble lvy                           = cSim.pImageLVelY[index];
        PMEDouble lvz                           = cSim.pImageLVelZ[index];
        PMEDouble q                             = cSim.pImageCharge[index];
        PMEDouble m                             = cSim.pImageMass[index];
        PMEFloat2 sigeps                        = cSim.pImageSigEps[index];
        unsigned int outputBuffers              = cSim.pImageOutputBuffers[index];
        unsigned int cellID                     = cSim.pImageCellID[index];
        cSim.pImageX2[pos]                      = x;
        cSim.pImageY2[pos]                      = y;
        cSim.pImageZ2[pos]                      = z;
        cSim.pImageAtom2[pos]                   = atom;
        cSim.pImageAtomLookup[atom]             = pos;
        PMEFloat2 xy                            = {x, y};
        cSim.pAtomXYSaveSP[pos]                 = xy;
        cSim.pAtomZSaveSP[pos]                  = z;
        cSim.pImageVelX2[pos]                   = vx;
        cSim.pImageVelY2[pos]                   = vy;
        cSim.pImageVelZ2[pos]                   = vz;
        cSim.pImageLVelX2[pos]                  = lvx;
        cSim.pImageLVelY2[pos]                  = lvy;
        cSim.pImageLVelZ2[pos]                  = lvz;
        cSim.pImageCharge2[pos]                 = q;
        cSim.pAtomChargeSP[pos]                 = q;
        cSim.pImageMass2[pos]                   = m;
        cSim.pImageInvMass2[pos]                = (m != (PMEDouble)0.0 ? (PMEDouble)1.0 / m : (PMEDouble)0.0);
        cSim.pImageSigEps2[pos]                 = sigeps;
        cSim.pImageOutputBuffers2[pos]          = outputBuffers;
        cSim.pImageCellID2[pos]                 = cellID;
        
        // Advance to next atom
        index                                   = newindex;    
        pos                                     = newpos;    
    }
}


extern "C" void kNLRemapImage(gpuContext gpu)
{
    kNLRemapImage_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>(gpu->sim.pImageIndex);
    LAUNCHERROR("kNLRemapImage");

    unsigned int *pUint;
    PMEDouble* pPMEDouble;
    PMEFloat2* pPMEFloat2;
    
    // Remap constant memory pointers
    pUint                           = gpu->sim.pImageAtom;
    gpu->sim.pImageAtom             = gpu->sim.pImageAtom2;
    gpu->sim.pImageAtom2            = pUint;
    pPMEDouble                      = gpu->sim.pImageX;
    gpu->sim.pImageX                = gpu->sim.pImageX2;
    gpu->sim.pImageX2               = pPMEDouble;
    pPMEDouble                      = gpu->sim.pImageY;
    gpu->sim.pImageY                = gpu->sim.pImageY2;
    gpu->sim.pImageY2               = pPMEDouble;
    pPMEDouble                      = gpu->sim.pImageZ;
    gpu->sim.pImageZ                = gpu->sim.pImageZ2;
    gpu->sim.pImageZ2               = pPMEDouble; 
    pPMEDouble                      = gpu->sim.pImageVelX;
    gpu->sim.pImageVelX             = gpu->sim.pImageVelX2;
    gpu->sim.pImageVelX2            = pPMEDouble;
    pPMEDouble                      = gpu->sim.pImageVelY;
    gpu->sim.pImageVelY             = gpu->sim.pImageVelY2;
    gpu->sim.pImageVelY2            = pPMEDouble;
    pPMEDouble                      = gpu->sim.pImageVelZ;
    gpu->sim.pImageVelZ             = gpu->sim.pImageVelZ2;
    gpu->sim.pImageVelZ2            = pPMEDouble;
    pPMEDouble                      = gpu->sim.pImageLVelX;
    gpu->sim.pImageLVelX            = gpu->sim.pImageLVelX2;
    gpu->sim.pImageLVelX2           = pPMEDouble;
    pPMEDouble                      = gpu->sim.pImageLVelY;
    gpu->sim.pImageLVelY            = gpu->sim.pImageLVelY2;
    gpu->sim.pImageLVelY2           = pPMEDouble;
    pPMEDouble                      = gpu->sim.pImageLVelZ;
    gpu->sim.pImageLVelZ            = gpu->sim.pImageLVelZ2;
    gpu->sim.pImageLVelZ2           = pPMEDouble;    
    pPMEDouble                      = gpu->sim.pImageCharge;
    gpu->sim.pImageCharge           = gpu->sim.pImageCharge2;
    gpu->sim.pImageCharge2          = pPMEDouble;
    pPMEDouble                      = gpu->sim.pImageMass;
    gpu->sim.pImageMass             = gpu->sim.pImageMass2;
    gpu->sim.pImageMass2            = pPMEDouble;
    pPMEDouble                      = gpu->sim.pImageInvMass;
    gpu->sim.pImageInvMass          = gpu->sim.pImageInvMass2;
    gpu->sim.pImageInvMass2         = pPMEDouble;
    pPMEFloat2                      = gpu->sim.pImageSigEps;
    gpu->sim.pImageSigEps           = gpu->sim.pImageSigEps2;
    gpu->sim.pImageSigEps2          = pPMEFloat2;
    pUint                           = gpu->sim.pImageOutputBuffers;
    gpu->sim.pImageOutputBuffers    = gpu->sim.pImageOutputBuffers2;
    gpu->sim.pImageOutputBuffers2   = pUint;
    pUint                           = gpu->sim.pImageCellID;
    gpu->sim.pImageCellID           = gpu->sim.pImageCellID2;
    gpu->sim.pImageCellID2          = pUint;
}

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kNLRemapLocalInteractions_kernel()
{
    unsigned int pos                = blockIdx.x * blockDim.x + threadIdx.x;    
    
    // Remap bond forces
    while (pos < cSim.bondOffset)
    {
        if (pos < cSim.bonds)
        {
            int4 atom                       = cSim.pBondID[pos];
            atom.x                          = cSim.pImageAtomLookup[atom.x];
            atom.y                          = cSim.pImageAtomLookup[atom.y];
            atom.z                         += atom.x;
            atom.w                         += atom.y;         
            cSim.pImageBondID[pos]          = atom;
        }
        pos                                += blockDim.x * gridDim.x;
    }

    
    // Remap bond angle forces
    while (pos < cSim.bondAngleOffset)
    {
        pos                                -= cSim.bondOffset;
        if (pos < cSim.bondAngles)
        {
            int4 atom1                      = cSim.pBondAngleID1[pos];
            int2 atom2                      = cSim.pBondAngleID2[pos];
            atom1.x                         = cSim.pImageAtomLookup[atom1.x];
            atom1.y                         = cSim.pImageAtomLookup[atom1.y];
            atom1.z                         = cSim.pImageAtomLookup[atom1.z];
            atom1.w                        += atom1.x;
            atom2.x                        += atom1.y;
            atom2.y                        += atom1.z;
            cSim.pImageBondAngleID1[pos]    = atom1;
            cSim.pImageBondAngleID2[pos]    = atom2;
        }
        pos                                += cSim.bondOffset + blockDim.x * gridDim.x;
    }
    
    // Remap dihedral forces
    while (pos < cSim.dihedralOffset)
    {
        pos                                -= cSim.bondAngleOffset;
        if (pos < cSim.dihedrals)
        {
            int4 atom1                      = cSim.pDihedralID1[pos];
            int4 atom2                      = cSim.pDihedralID2[pos];
            atom1.x                         = cSim.pImageAtomLookup[atom1.x];
            atom1.y                         = cSim.pImageAtomLookup[atom1.y];
            atom1.z                         = cSim.pImageAtomLookup[atom1.z];
            atom1.w                         = cSim.pImageAtomLookup[atom1.w];
            atom2.x                        += atom1.x;
            atom2.y                        += atom1.y;
            atom2.z                        += atom1.z;
            atom2.w                        += atom1.w;
            cSim.pImageDihedralID1[pos]     = atom1;
            cSim.pImageDihedralID2[pos]     = atom2;
        }
        pos                                += cSim.bondAngleOffset + blockDim.x * gridDim.x;
    }
   
   
    // Remap 1-4 forces
    while (pos < cSim.nb14Offset)
    {
        pos                                -= cSim.dihedralOffset;
        if (pos < cSim.nb14s)
        {
            int4 atom                       = cSim.pNb14ID[pos];
            atom.x                          = cSim.pImageAtomLookup[atom.x];
            atom.y                          = cSim.pImageAtomLookup[atom.y];
            atom.z                         += atom.x;
            atom.w                         += atom.y;
            cSim.pImageNb14ID[pos]          = atom;
        }
        pos                                += cSim.dihedralOffset + blockDim.x * gridDim.x;
    }    
    
    // Remap Constraint forces
    while (pos < cSim.constraintOffset)
    {
        pos                                -= cSim.nb14Offset;
        if (pos < cSim.constraints)
        {
            int2 atom                       = cSim.pConstraintID[pos];
            atom.x                          = cSim.pImageAtomLookup[atom.x];
            atom.y                         += atom.x;
            cSim.pImageConstraintID[pos]    = atom;
        }
        pos                                += cSim.nb14Offset + blockDim.x * gridDim.x;
    }
    pos                                    -= cSim.constraintOffset;
    
    // Remap Charmm interactions
    while (pos < cSim.UBAngleOffset)
    {
        if (pos < cSim.UBAngles)
        {
            int4 atom                       = cSim.pUBAngleID[pos];
            atom.x                          = cSim.pImageAtomLookup[atom.x];
            atom.y                          = cSim.pImageAtomLookup[atom.y];
            atom.z                         += atom.x;
            atom.w                         += atom.y;         
            cSim.pImageUBAngleID[pos]       = atom;
        }
        pos                                += blockDim.x * gridDim.x;    
    }
    
    while (pos < cSim.impDihedralOffset)
    {
        pos                                -= cSim.UBAngleOffset;
        if (pos < cSim.impDihedrals)
        {
            int4 atom1                      = cSim.pImpDihedralID1[pos];
            int4 atom2                      = cSim.pImpDihedralID2[pos];
            atom1.x                         = cSim.pImageAtomLookup[atom1.x];
            atom1.y                         = cSim.pImageAtomLookup[atom1.y];
            atom1.z                         = cSim.pImageAtomLookup[atom1.z];
            atom1.w                         = cSim.pImageAtomLookup[atom1.w];
            atom2.x                        += atom1.x;
            atom2.y                        += atom1.y;
            atom2.z                        += atom1.z;
            atom2.w                        += atom1.w;
            cSim.pImageImpDihedralID1[pos]  = atom1;
            cSim.pImageImpDihedralID2[pos]  = atom2;
        }
        pos                                += cSim.UBAngleOffset + blockDim.x * gridDim.x;    
    }
    
    while (pos < cSim.cmapOffset)
    {
        pos                                -= cSim.impDihedralOffset;
        if (pos < cSim.cmaps)
        {       
            int4 atom1                      = cSim.pCmapID1[pos];
            int4 atom2                      = cSim.pCmapID2[pos];
            int2 atom3                      = cSim.pCmapID3[pos];
            atom1.x                         = cSim.pImageAtomLookup[atom1.x];
            atom1.y                         = cSim.pImageAtomLookup[atom1.y];
            atom1.z                         = cSim.pImageAtomLookup[atom1.z];
            atom1.w                         = cSim.pImageAtomLookup[atom1.w];
            atom2.x                         = cSim.pImageAtomLookup[atom2.x];
            atom2.y                        += atom1.x;
            atom2.z                        += atom1.y;
            atom2.w                        += atom1.z;
            atom3.x                        += atom1.w;
            atom3.y                        += atom2.x;
            cSim.pImageCmapID1[pos]         = atom1;
            cSim.pImageCmapID2[pos]         = atom2;
            cSim.pImageCmapID3[pos]         = atom3;
        }
        pos                                += cSim.impDihedralOffset + blockDim.x * gridDim.x;    
    }    
    pos                                    -= cSim.cmapOffset;

    // Remap Shake constraints
    while (pos < cSim.shakeOffset)
    {
        if (pos < cSim.shakeConstraints)
        {
            int4 atom                       = cSim.pShakeID[pos];
            atom.x                          = cSim.pImageAtomLookup[atom.x];
            atom.y                          = cSim.pImageAtomLookup[atom.y];
            if (atom.z != -1)
                atom.z                      = cSim.pImageAtomLookup[atom.z];
            if (atom.w != -1)
                atom.w                      = cSim.pImageAtomLookup[atom.w];
            cSim.pImageShakeID[pos]         = atom;
        }
        pos                                += blockDim.x * gridDim.x;
    }
    pos                                    -= cSim.shakeOffset;

    while (pos < cSim.fastShakeOffset)
    {
        if (pos < cSim.fastShakeConstraints)
        {
            int4 atom                       = cSim.pFastShakeID[pos];
            atom.x                          = cSim.pImageAtomLookup[atom.x];
            atom.y                          = cSim.pImageAtomLookup[atom.y];
            atom.z                          = cSim.pImageAtomLookup[atom.z];
            cSim.pImageFastShakeID[pos]     = atom;
        }
        pos                                += blockDim.x * gridDim.x;
    }
    pos                                    -= cSim.fastShakeOffset;

    while (pos < cSim.slowShakeOffset)
    {
        if (pos < cSim.slowShakeConstraints)
        {
            int atom1                       = cSim.pSlowShakeID1[pos];
            int4 atom2                      = cSim.pSlowShakeID2[pos];
            atom1                           = cSim.pImageAtomLookup[atom1];
            atom2.x                         = cSim.pImageAtomLookup[atom2.x];
            atom2.y                         = cSim.pImageAtomLookup[atom2.y];
            atom2.z                         = cSim.pImageAtomLookup[atom2.z];
            atom2.w                         = cSim.pImageAtomLookup[atom2.w];
            cSim.pImageSlowShakeID1[pos]    = atom1;
            cSim.pImageSlowShakeID2[pos]    = atom2;
        }
        pos                                += blockDim.x * gridDim.x;
    }
    pos                                    -= cSim.slowShakeOffset;
    
    while (pos < cSim.soluteAtoms)          // Solute atoms already padded to warp width
    {
        int atom                            = cSim.pSoluteAtomID[pos];
        if (atom != -1)
            atom                            = cSim.pImageAtomLookup[atom];
        cSim.pImageSoluteAtomID[pos]        = atom;         
        pos                                += blockDim.x * gridDim.x;
    }
    pos                                    -= cSim.soluteAtoms;

    while (pos < cSim.solventMoleculeStride)
    {
        if (pos < cSim.solventMolecules)
        {
            int4 atom                           = cSim.pSolventAtomID[pos];
            atom.x                              = cSim.pImageAtomLookup[atom.x];
            if (atom.y != -1)
                atom.y                          = cSim.pImageAtomLookup[atom.y];
            if (atom.z != -1)
                atom.z                          = cSim.pImageAtomLookup[atom.z];
            if (atom.w != -1)
                atom.w                          = cSim.pImageAtomLookup[atom.w];
            cSim.pImageSolventAtomID[pos]       = atom;
        }
        pos                                    += blockDim.x * gridDim.x;
    }
    pos                                        -= cSim.solventMoleculeStride;

    while (pos < cSim.EP11Offset)
    {  
        if (pos < cSim.EP11s)
        {
            int4 frame                          = cSim.pExtraPoint11Frame[pos];
            int index                           = cSim.pExtraPoint11Index[pos];   
            frame.x                             = cSim.pImageAtomLookup[frame.x];
            frame.y                             = cSim.pImageAtomLookup[frame.y];
            frame.z                             = cSim.pImageAtomLookup[frame.z];
            frame.w                             = cSim.pImageAtomLookup[frame.w];
            index                               = cSim.pImageAtomLookup[index];
            cSim.pImageExtraPoint11Frame[pos]   = frame;
            cSim.pImageExtraPoint11Index[pos]   = index;
        }
        pos                                    += blockDim.x * gridDim.x;
    }
    while (pos < cSim.EP12Offset)
    {
        pos                                    -= cSim.EP11Offset;
        if (pos < cSim.EP12s)
        {
            int4 frame                          = cSim.pExtraPoint12Frame[pos];
            int index                           = cSim.pExtraPoint12Index[pos];   
            frame.x                             = cSim.pImageAtomLookup[frame.x];
            frame.y                             = cSim.pImageAtomLookup[frame.y];
            frame.z                             = cSim.pImageAtomLookup[frame.z];
            frame.w                             = cSim.pImageAtomLookup[frame.w];
            index                               = cSim.pImageAtomLookup[index];
            cSim.pImageExtraPoint12Frame[pos]   = frame;
            cSim.pImageExtraPoint12Index[pos]   = index;
        }
        pos                                    += cSim.EP11Offset + blockDim.x * gridDim.x;
    }
    while (pos < cSim.EP21Offset)
    {
        pos                                    -= cSim.EP12Offset;        
        if (pos < cSim.EP21s)
        {
            int4 frame                          = cSim.pExtraPoint21Frame[pos];
            int2 index                          = cSim.pExtraPoint21Index[pos];   
            frame.x                             = cSim.pImageAtomLookup[frame.x];
            frame.y                             = cSim.pImageAtomLookup[frame.y];
            frame.z                             = cSim.pImageAtomLookup[frame.z];
            frame.w                             = cSim.pImageAtomLookup[frame.w];
            index.x                             = cSim.pImageAtomLookup[index.x];
            index.y                             = cSim.pImageAtomLookup[index.y];
            cSim.pImageExtraPoint21Frame[pos]   = frame;
            cSim.pImageExtraPoint21Index[pos]   = index;
        }
        pos                                    += cSim.EP12Offset + blockDim.x * gridDim.x;
    }   
    while (pos < cSim.EP22Offset)
    {
        pos                                    -= cSim.EP21Offset;
        if (pos < cSim.EP22s)
        {
            int4 frame                          = cSim.pExtraPoint22Frame[pos];
            int2 index                          = cSim.pExtraPoint22Index[pos];   
            frame.x                             = cSim.pImageAtomLookup[frame.x];
            frame.y                             = cSim.pImageAtomLookup[frame.y];
            frame.z                             = cSim.pImageAtomLookup[frame.z];
            frame.w                             = cSim.pImageAtomLookup[frame.w];
            index.x                             = cSim.pImageAtomLookup[index.x];
            index.y                             = cSim.pImageAtomLookup[index.y];
            cSim.pImageExtraPoint22Frame[pos]   = frame;
            cSim.pImageExtraPoint22Index[pos]   = index;
        }
        pos                                    += cSim.EP21Offset + blockDim.x * gridDim.x;
    }
}

extern "C" void kNLRemapLocalInteractions(gpuContext gpu)
{
    kNLRemapLocalInteractions_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kNLRemapLocalInteractions");   
}

// Clear all cell boundaries in case some are empty
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kNLClearCellBoundaries_kernel()
{
    unsigned int pos                        = blockIdx.x * blockDim.x + threadIdx.x;
    uint2 nulldata                          = {0, 0};
    while (pos < cSim.cells)
    {
        cSim.pNLNonbondCellStartEnd[pos]    = nulldata;
        pos                                += blockDim.x * gridDim.x;
    }
}

extern "C" void kNLClearCellBoundaries(gpuContext gpu)
{
    kNLClearCellBoundaries_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kNLClearCellBoundaries");   
}


__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kNLCalculateCellBoundaries_kernel(unsigned int* pImageHash)
{
#if (__CUDA_ARCH__ >= 200)
const int cSpan = 12000;
#else
const int cSpan = 4000;
#endif
__shared__ unsigned int sHash[cSpan];
    int pos                                                 = ((cSim.atoms + 1) * blockIdx.x) / gridDim.x - 1;
    int end                                                 = ((cSim.atoms + 1) * (blockIdx.x + 1)) / gridDim.x;
    while (pos < end)
    {
        // Read span to check for transitions
        int pos1                                            = pos + threadIdx.x;
        int spos                                            = threadIdx.x;
        int span                                            = min(end, pos + cSpan);
        int lSpan                                           = min(end, pos + cSpan) - pos;
        while (pos1 < span)
        {
            // Read hash data or 0s on either end to force transitions
            if ((pos1 >= 0) && (pos1 < cSim.atoms)) 
                sHash[spos]                                 = pImageHash[pos1];
            else
                sHash[spos]                                 = 0;
            pos1                                           += blockDim.x;
            spos                                           += blockDim.x;
        }
        __syncthreads();
        spos                                                = threadIdx.x + 1;
        while (spos < lSpan)
        {
            int oldHash                                     = sHash[spos - 1] >> CELLHASHBITS;
            int newHash                                     = sHash[spos] >> CELLHASHBITS;
            if (oldHash != newHash)
            {
                if (pos + spos != cSim.atoms)
                    cSim.pNLNonbondCellStartEnd[newHash].x  = pos + spos;
                if (pos + spos != 0)
                    cSim.pNLNonbondCellStartEnd[oldHash].y  = pos + spos;   
            }
            spos                                           += blockDim.x;
        }      
        __syncthreads();  
        pos                                                 = min(end, pos + cSpan - 1);        
    }
}

extern "C" void kNLCalculateCellBoundaries(gpuContext gpu)
{
    kNLCalculateCellBoundaries_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>(gpu->sim.pImageHash);   
    LAUNCHERROR("kNLCalculateCellBoundaries");
}



__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kNLClearExclusionMasks_kernel()
{

}

extern "C" void kNLClearExclusionMasks(gpuContext gpu)
{
    kNLClearExclusionMasks_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();   
    LAUNCHERROR("kNLClearExclusionMasks");
}

struct GEData {
    unsigned int workUnit;
    unsigned int xCellStart;
    unsigned int xCellEnd;
    unsigned int yCellStart;
    unsigned int yCellEnd;
    unsigned int exclusionMap;
    unsigned int imageAtom;
    unsigned int atom;
};

#if 0
__global__ void
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_NLCALCULATE_OFFSETS_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_NLCALCULATE_OFFSETS_THREADS_PER_BLOCK, 1)
#endif
kNLCalculateOffsets_kernel()
{
#if (__CUDA_ARCH__ >= 200)
__shared__ volatile SNLRecord sNLRecord[SM_2X_NLCALCULATE_OFFSETS_THREADS_PER_BLOCK];
#else
__shared__ volatile SNLRecord sNLRecord[SM_13_NLCALCULATE_OFFSETS_THREADS_PER_BLOCK];
#endif
    unsigned int pos                            = blockIdx.x * blockDim.x;
    unsigned int warps                          = blockDim.x >> 4;

    while (pos < cSim.NLSize)
    {
        // Read NLRecord entries
        unsigned int rpos                       = threadIdx.x >> 4;
        unsigned int tpos                       = threadIdx.x - (rpos << 4);
        unsigned int reads                      = min(blockDim.x, cSim.NLSize - pos);
        
        while (rpos < reads)
        {
            sNLRecord[rpos].array[tpos]         = cSim.pNLRecord[pos + rpos].array[tpos];
            rpos                               += warps;
        }
        __syncthreads();
        
        // Read y cell data
        if (pos + threadIdx.x < cSim.NLSize)
        {
            uint2 homeCell                      = cSim.pNLNonbondCellStartEnd[sNLRecord[threadIdx.x].NL.homeCell >> NLCELLSHIFT];
            int startOffset                     = sNLRecord[threadIdx.x].NL.homeCell & NLBUFFERMASK;
            int ysize                           = max(0, (int)(((homeCell.y - homeCell.x + (GRID - 1)) >> GRIDBITS)) - startOffset) / cSim.NLYDivisor;
            
            // Calculate maximum required space
            int xsize                           = 0;
            for (int i = 0; i < sNLRecord[i].NL.neighborCells; i++)
            {
                uint2 cell                      = cSim.pNLNonbondCellStartEnd[sNLRecord[threadIdx.x].NL.neighborCell[i] >> NLCELLSHIFT];
                xsize                          += cell.y - cell.x;
            }
            xsize                               = ((xsize + GRID - 1) >> GRIDBITS) << GRIDBITS;
            cSim.pNLOffset[pos + threadIdx.x]   =  xsize * ysize;
        }
    
        pos                                    += blockDim.x * gridDim.x;
    }
}

extern "C" void kNLCalculateOffsets(gpuContext gpu)
{
    kNLCalculateOffsets_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();   
    LAUNCHERROR("kNLCalculateOffsets");
    
    cudaThreadSynchronize();
    gpu->pbNLOffset->Download();
    for (int i = 0; i < gpu->sim.NLSize; i++)
        printf("%4d %6d\n", i, gpu->pbNLOffset->_pSysData[i]);
#if 0        
    gpu->pbNLNonbondCellStartEnd->Download();  
    for (int i = 0; i < gpu->sim.cells; i++)
        printf("%4d %6d\n", i, gpu->pbNLNonbondCellStartEnd->_pSysData[i].y - gpu->pbNLNonbondCellStartEnd->_pSysData[i].x);
#endif
        
    exit(-1);
    
}
#endif

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kNLCalculateCellCoordinates_kernel()
#include "kCCC.h"

#define PME_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kNLCalculateCellCoordinatesOrthogonal_kernel()
#include "kCCC.h"
#undef PME_ORTHOGONAL

#define PME_NTP
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kNLCalculateCellCoordinatesNTP_kernel()
#include "kCCC.h"

#define PME_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kNLCalculateCellCoordinatesOrthogonalNTP_kernel()
#include "kCCC.h"
#undef PME_ORTHOGONAL
#undef PME_NTP

extern "C" void kNLCalculateCellCoordinates(gpuContext gpu)
{
    if (gpu->sim.ntp > 0)
    {
        if (gpu->sim.is_orthog)    
            kNLCalculateCellCoordinatesOrthogonalNTP_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>();  
        else
            kNLCalculateCellCoordinatesNTP_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>();         
    }
    else
    {
        if (gpu->sim.is_orthog)    
            kNLCalculateCellCoordinatesOrthogonal_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>();  
        else
            kNLCalculateCellCoordinates_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>(); 
    } 
    cudaThreadSynchronize();

#if 0
    gpu->pbAtomXYSP->Download();
    gpu->pbAtomZSP->Download();
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        printf("%6d: %16.8f %16.8f %16.8f\n", i, gpu->pbAtomXYSP->_pSysData[i].x, gpu->pbAtomXYSP->_pSysData[i].y, gpu->pbAtomZSP->_pSysData[i]);
    }
#endif    
    
    LAUNCHERROR_NONBLOCKING("kNLCalculateCellCoordinates");
}


#define PME_ATOMS_PER_WARP (32)
#define PME_IS_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK, 1)
#endif
kNLBuildNeighborListOrthogonal32_kernel()
#include "kBNL.h"

#define PME_VIRIAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK, 1)
#endif
kNLBuildNeighborListOrthogonalNTP32_kernel()
#include "kBNL.h"
#undef PME_VIRIAL
#undef PME_IS_ORTHOGONAL

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK, 1)
#endif
kNLBuildNeighborList32_kernel()
#include "kBNL.h"

#define PME_VIRIAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK, 1)
#endif
kNLBuildNeighborListNTP32_kernel()
#include "kBNL.h"
#undef PME_VIRIAL
#undef PME_ATOMS_PER_WARP


#define PME_ATOMS_PER_WARP (16)
#define PME_IS_ORTHOGONAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK + 64, 1)
#else
__launch_bounds__(SM_13_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK, 1)
#endif
kNLBuildNeighborListOrthogonal16_kernel()
#include "kBNL.h"

#define PME_VIRIAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK + 64, 1)
#else
__launch_bounds__(SM_13_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK, 1)
#endif
kNLBuildNeighborListOrthogonalNTP16_kernel()
#include "kBNL.h"
#undef PME_VIRIAL
#undef PME_IS_ORTHOGONAL

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK + 64, 1)
#else
__launch_bounds__(SM_13_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK, 1)
#endif
kNLBuildNeighborList16_kernel()
#include "kBNL.h"

#define PME_VIRIAL
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK + 64, 1)
#else
__launch_bounds__(SM_13_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK, 1)
#endif
kNLBuildNeighborListNTP16_kernel()
#include "kBNL.h"
#undef PME_VIRIAL
#undef PME_ATOMS_PER_WARP

extern "C" void kNLBuildNeighborList(gpuContext gpu)
{
//    cudaThreadSetLimit(cudaLimitPrintfFifoSize, 200000000);   
//    printf("%d %d %d\n", gpu->sim.NLAtomsPerWarp, gpu->sim.ntp, gpu->sim.is_orthog);
//    printf("%06d %06d\n", gpu->NLBuildNeighborList16ThreadsPerBlock, sizeof(uint) * gpu->sim.NLMaxExclusionsPerWarp * gpu->NLBuildNeighborList16ThreadsPerBlock / GRID);
//    printf("%06d %06d\n", gpu->NLBuildNeighborList32ThreadsPerBlock, sizeof(uint) * gpu->sim.NLMaxExclusionsPerWarp * gpu->NLBuildNeighborList32ThreadsPerBlock / GRID);
    if (gpu->sim.NLAtomsPerWarp == 32)
    {
        if (gpu->sim.ntp > 0)
        {
            if (gpu->sim.is_orthog)    
                kNLBuildNeighborListOrthogonalNTP32_kernel<<<gpu->blocks, gpu->NLBuildNeighborList32ThreadsPerBlock, sizeof(uint) * gpu->sim.NLMaxExclusionsPerWarp * gpu->NLBuildNeighborList32ThreadsPerBlock / GRID>>>();  
            else
                kNLBuildNeighborListNTP32_kernel<<<gpu->blocks, gpu->NLBuildNeighborList32ThreadsPerBlock, sizeof(uint) * gpu->sim.NLMaxExclusionsPerWarp * gpu->NLBuildNeighborList32ThreadsPerBlock / GRID>>>();         
        }
        else
        {
            if (gpu->sim.is_orthog)    
                kNLBuildNeighborListOrthogonal32_kernel<<<gpu->blocks, gpu->NLBuildNeighborList32ThreadsPerBlock, sizeof(uint) * gpu->sim.NLMaxExclusionsPerWarp * gpu->NLBuildNeighborList32ThreadsPerBlock / GRID>>>();  
            else
                kNLBuildNeighborList32_kernel<<<gpu->blocks, gpu->NLBuildNeighborList32ThreadsPerBlock, sizeof(uint) * gpu->sim.NLMaxExclusionsPerWarp * gpu->NLBuildNeighborList32ThreadsPerBlock / GRID>>>(); 
        }     
    }
    else
    {
        if (gpu->sim.ntp > 0)
        {
            if (gpu->sim.is_orthog)    
                kNLBuildNeighborListOrthogonalNTP16_kernel<<<gpu->blocks, gpu->NLBuildNeighborList16ThreadsPerBlock, sizeof(uint) * gpu->sim.NLMaxExclusionsPerWarp * gpu->NLBuildNeighborList16ThreadsPerBlock / GRID>>>();  
            else
                kNLBuildNeighborListNTP16_kernel<<<gpu->blocks, gpu->NLBuildNeighborList16ThreadsPerBlock, sizeof(uint) * gpu->sim.NLMaxExclusionsPerWarp * gpu->NLBuildNeighborList16ThreadsPerBlock / GRID>>>();         
        }
        else
        {
            if (gpu->sim.is_orthog)    
                kNLBuildNeighborListOrthogonal16_kernel<<<gpu->blocks, gpu->NLBuildNeighborList16ThreadsPerBlock, sizeof(uint) * gpu->sim.NLMaxExclusionsPerWarp * gpu->NLBuildNeighborList16ThreadsPerBlock / GRID>>>();  
            else
                kNLBuildNeighborList16_kernel<<<gpu->blocks, gpu->NLBuildNeighborList16ThreadsPerBlock, sizeof(uint) * gpu->sim.NLMaxExclusionsPerWarp * gpu->NLBuildNeighborList16ThreadsPerBlock / GRID>>>(); 
        }     
    }

    LAUNCHERROR("kNLBuildNeighborList");    

#if 0
    cudaThreadSynchronize();  
    gpu->pbNLTotalOffset->Download();
    fprintf(stdout, "Total: %d, maximum %d\n", gpu->pbNLTotalOffset->_pSysData[0], gpu->sim.NLMaxTotalOffset);
#endif
    
#if 0   
    cudaThreadSynchronize();
    gpu->pbNLOffset->Download();  
    for (int i = 0; i < gpu->sim.NLSize; i++)
        printf("%4d %12u\n", i, gpu->pbNLOffset->_pSysData[i]);
#endif

#if 0        
    gpu->pbNLNonbondCellStartEnd->Download();  
    for (int i = 0; i < gpu->sim.cells; i++)
        printf("%4d %6d\n", i, gpu->pbNLNonbondCellStartEnd->_pSysData[i].y - gpu->pbNLNonbondCellStartEnd->_pSysData[i].x);
#endif   
}

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_UPDATE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_UPDATE_THREADS_PER_BLOCK, 1)
#endif
kNLSkinTest_kernel()
{
    __shared__ volatile bool sbFail;
    __shared__ PMEFloat sOne_half_nonbond_skin_squared;
    unsigned int pos                        = blockIdx.x * blockDim.x + threadIdx.x;
    sbFail                                  = false; 
    if ((cSim.ntp > 0) && (threadIdx.x == 0))
        sOne_half_nonbond_skin_squared      = cSim.pNTPData->one_half_nonbond_skin_squared;
    __syncthreads();
    
    if (cSim.ntp > 0)
    {
        while (pos < cSim.atoms)
        {
            PMEFloat x                      = cSim.pImageX[pos];
            PMEFloat y                      = cSim.pImageY[pos];
            PMEFloat2 oldxy                 = cSim.pAtomXYSaveSP[pos];
            PMEFloat z                      = cSim.pImageZ[pos];
            PMEFloat oldz                   = cSim.pAtomZSaveSP[pos];
            PMEFloat dx                     = x - oldxy.x;
            PMEFloat dy                     = y - oldxy.y;
            PMEFloat dz                     = z - oldz;
            
            PMEFloat r2                     = dx * dx + dy * dy + dz * dz;
            if (r2 >= sOne_half_nonbond_skin_squared)
                sbFail                      = true;                  

            pos                            += blockDim.x * gridDim.x;    
        }    
    }
    else
    {
        while (pos < cSim.atoms)
        {
            PMEFloat x                      = cSim.pImageX[pos];
            PMEFloat y                      = cSim.pImageY[pos];
            PMEFloat2 oldxy                 = cSim.pAtomXYSaveSP[pos];
            PMEFloat z                      = cSim.pImageZ[pos];
            PMEFloat oldz                   = cSim.pAtomZSaveSP[pos];
            
            PMEFloat dx                     = x - oldxy.x;
            PMEFloat dy                     = y - oldxy.y;
            PMEFloat dz                     = z - oldz;
            
            PMEFloat r2                     = dx * dx + dy * dy + dz * dz;
            if (r2 >= cSim.one_half_nonbond_skin_squared)
                sbFail                      = true;                  

            pos                            += blockDim.x * gridDim.x;    
        }
    }
    
    __syncthreads();
    if ((threadIdx.x == 0) && sbFail)
        *cSim.pNLbSkinTestFail              = true;
}

extern "C" void kNLSkinTest(gpuContext gpu)
{
    static int counter = 0;
 
    *(gpu->pbNLbSkinTestFail->_pSysData) = false;    
    kNLSkinTest_kernel<<<gpu->blocks, gpu->updateThreadsPerBlock>>>();   
    LAUNCHERROR_NONBLOCKING("kNLSkinTest");
    cudaThreadSynchronize();
    
    if (*(gpu->pbNLbSkinTestFail->_pSysData))
    {
     //   printf("%d %f\n", i, gpu->sim.one_half_nonbond_skin_squared);
        gpu->bNeedNewNeighborList = true;
     //   fprintf(stdout, "%06d yes build\n", counter);   
    }
    // else    
    //    fprintf(stdout, "%06d no build\n", counter);

    counter++;
}

