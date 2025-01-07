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
#if (__CUDA_ARCH__ >= 200)
const int GRADSUMLOADSIZE = 32;          // Must be warpsize or smaller to support MPI output
#else
const int GRADSUMLOADSIZE = 16;          // Must be warpsize or smaller to support MPI output
#endif
const int GRADSUMGRID = 8;
const int GRADSUMGRIDBITS = 3; 
struct GradSumAtomData
{
    PMEFloat charge;
    int ix;
    int iy;
    int iz;
    PMEFloat tx[4];
    PMEFloat ty[4];
    PMEFloat tz0;
    PMEFloat tz1;
    PMEFloat tz2;
    PMEFloat tz3;
    PMEFloat dtx[4];
    PMEFloat dty[4];
    PMEFloat dtz0;
    PMEFloat dtz1;
    PMEFloat dtz2;
    PMEFloat dtz3;
    PMEFloat fx;
    PMEFloat fy;
    PMEFloat fz;
};

struct Factor {
    PMEFloat dfx;
    PMEFloat dfy;
    PMEFloat dfz;
};

__shared__ volatile GradSumAtomData sAtom[GRADSUMLOADSIZE];
__shared__ volatile Factor sFactor[GRADSUMTHREADS];
#ifdef PME_VIRIAL
__shared__ PMEFloat sRecipf[9];
    if (threadIdx.x < 9)
        sRecipf[threadIdx.x]                = cSim.pNTPData->recipf[threadIdx.x];
    __syncthreads();
#endif
#ifdef MPI
__shared__ volatile PMEFloat sPMEForce[3 * GRADSUMLOADSIZE];
#endif

    // Determine warp constants
    unsigned int tgx                        = threadIdx.x & (GRADSUMGRID - 1);
    volatile Factor* psFactor               = &sFactor[(threadIdx.x >> GRADSUMGRIDBITS) << GRADSUMGRIDBITS];

    // Determine grid offsets   
    const int tOffsetX                      = tgx & 0x03; 
    const int tOffsetY                      = (tgx & 0x07) >> 2;
    const int iOffsetX                      = tOffsetX;
    const int iOffsetY                      = tOffsetY;

    unsigned int pos                        = blockIdx.x * GRADSUMLOADSIZE; 
    
    // Read batch of atoms
    unsigned int maxatom                    = min(pos + GRADSUMLOADSIZE, cSim.atoms);
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
        PMEFloat4 dtx                       = cSim.pDThetaX[pos1];
        PMEFloat4 dty                       = cSim.pDThetaY[pos1];
        PMEFloat4 dtz                       = cSim.pDThetaZ[pos1];
        sAtom[threadIdx.x].charge           = charge;
        sAtom[threadIdx.x].ix               = ix;
        sAtom[threadIdx.x].iy               = iy;
        sAtom[threadIdx.x].iz               = iz;
        sAtom[threadIdx.x].tx[0]            = tx.x;
        sAtom[threadIdx.x].tx[1]            = tx.y;
        sAtom[threadIdx.x].tx[2]            = tx.z;
        sAtom[threadIdx.x].tx[3]            = tx.w;
        sAtom[threadIdx.x].ty[0]            = ty.x;
        sAtom[threadIdx.x].ty[1]            = ty.y;
        sAtom[threadIdx.x].ty[2]            = ty.z;
        sAtom[threadIdx.x].ty[3]            = ty.w;
        sAtom[threadIdx.x].tz0              = tz.x;
        sAtom[threadIdx.x].tz1              = tz.y;
        sAtom[threadIdx.x].tz2              = tz.z;
        sAtom[threadIdx.x].tz3              = tz.w;
        sAtom[threadIdx.x].dtx[0]           = dtx.x;
        sAtom[threadIdx.x].dtx[1]           = dtx.y;
        sAtom[threadIdx.x].dtx[2]           = dtx.z;
        sAtom[threadIdx.x].dtx[3]           = dtx.w;
        sAtom[threadIdx.x].dty[0]           = dty.x;
        sAtom[threadIdx.x].dty[1]           = dty.y;
        sAtom[threadIdx.x].dty[2]           = dty.z;
        sAtom[threadIdx.x].dty[3]           = dty.w;
        sAtom[threadIdx.x].dtz0             = dtz.x;
        sAtom[threadIdx.x].dtz1             = dtz.y;
        sAtom[threadIdx.x].dtz2             = dtz.z;
        sAtom[threadIdx.x].dtz3             = dtz.w;        
    }
    __syncthreads();
        
    // Process batch 
    pos1                                    = threadIdx.x / GRADSUMGRID;
    unsigned int lastAtom                   = min(GRADSUMLOADSIZE, cSim.atoms - pos);
    while (pos1 < lastAtom)
    {
        // Calculate values
        int ix                              = sAtom[pos1].ix + iOffsetX;
        int iy0                             = sAtom[pos1].iy + iOffsetY;
        int iy1                             = iy0 + 2;
        int iz0                             = sAtom[pos1].iz;
        int iz1                             = iz0 + 1;
        int iz2                             = iz0 + 2;
        int iz3                             = iz0 + 3;
            
        // Insure coordinates stay in bounds
        if (ix >= cSim.nfft1)
            ix                             -= cSim.nfft1;
        if (iy0 >= cSim.nfft2)
            iy0                            -= cSim.nfft2;
        if (iy1 >= cSim.nfft2)
            iy1                            -= cSim.nfft2;
        if (iz0 >= cSim.nfft3)
            iz0                            -= cSim.nfft3;  
        if (iz1 >= cSim.nfft3)
            iz1                            -= cSim.nfft3;
        if (iz2 >= cSim.nfft3)
            iz2                            -= cSim.nfft3;  
        if (iz3 >= cSim.nfft3)
            iz3                            -= cSim.nfft3;  
                
            // Calculate interpolation values and destinations            
#ifdef use_DPDP     
        int2 i2term0                        = tex1Dfetch(texref, ((iz0 * cSim.nfft2 + iy0) * cSim.nfft1 + ix)); 
        int2 i2term1                        = tex1Dfetch(texref, ((iz1 * cSim.nfft2 + iy0) * cSim.nfft1 + ix)); 
        int2 i2term2                        = tex1Dfetch(texref, ((iz2 * cSim.nfft2 + iy0) * cSim.nfft1 + ix)); 
        int2 i2term3                        = tex1Dfetch(texref, ((iz3 * cSim.nfft2 + iy0) * cSim.nfft1 + ix));
        int2 i2term4                        = tex1Dfetch(texref, ((iz0 * cSim.nfft2 + iy1) * cSim.nfft1 + ix)); 
        int2 i2term5                        = tex1Dfetch(texref, ((iz1 * cSim.nfft2 + iy1) * cSim.nfft1 + ix)); 
        int2 i2term6                        = tex1Dfetch(texref, ((iz2 * cSim.nfft2 + iy1) * cSim.nfft1 + ix)); 
        int2 i2term7                        = tex1Dfetch(texref, ((iz3 * cSim.nfft2 + iy1) * cSim.nfft1 + ix)); 
        PMEFloat qterm0                     = __hiloint2double(i2term0.y, i2term0.x);
        PMEFloat qterm1                     = __hiloint2double(i2term1.y, i2term1.x);
        PMEFloat qterm2                     = __hiloint2double(i2term2.y, i2term2.x);
        PMEFloat qterm3                     = __hiloint2double(i2term3.y, i2term3.x);
        PMEFloat qterm4                     = __hiloint2double(i2term4.y, i2term4.x);
        PMEFloat qterm5                     = __hiloint2double(i2term5.y, i2term5.x); 
        PMEFloat qterm6                     = __hiloint2double(i2term6.y, i2term6.x);
        PMEFloat qterm7                     = __hiloint2double(i2term7.y, i2term7.x);     
#else                           
        PMEFloat qterm0                     = tex1Dfetch(texref, (iz0 * cSim.nfft2 + iy0) * cSim.nfft1 + ix); 
        PMEFloat qterm1                     = tex1Dfetch(texref, (iz1 * cSim.nfft2 + iy0) * cSim.nfft1 + ix); 
        PMEFloat qterm2                     = tex1Dfetch(texref, (iz2 * cSim.nfft2 + iy0) * cSim.nfft1 + ix); 
        PMEFloat qterm3                     = tex1Dfetch(texref, (iz3 * cSim.nfft2 + iy0) * cSim.nfft1 + ix);
        PMEFloat qterm4                     = tex1Dfetch(texref, (iz0 * cSim.nfft2 + iy1) * cSim.nfft1 + ix); 
        PMEFloat qterm5                     = tex1Dfetch(texref, (iz1 * cSim.nfft2 + iy1) * cSim.nfft1 + ix); 
        PMEFloat qterm6                     = tex1Dfetch(texref, (iz2 * cSim.nfft2 + iy1) * cSim.nfft1 + ix); 
        PMEFloat qterm7                     = tex1Dfetch(texref, (iz3 * cSim.nfft2 + iy1) * cSim.nfft1 + ix);
#endif
        PMEFloat wx                         = sAtom[pos1].tx[tOffsetX];
        PMEFloat wy0                        = sAtom[pos1].ty[tOffsetY];
        PMEFloat wy1                        = sAtom[pos1].ty[tOffsetY + 2];
        PMEFloat dwx                        = sAtom[pos1].dtx[tOffsetX];
        PMEFloat dwy0                       = sAtom[pos1].dty[tOffsetY];
        PMEFloat dwy1                       = sAtom[pos1].dty[tOffsetY + 2];
        PMEFloat tz0                        = sAtom[pos1].tz0;
        PMEFloat dtz0                       = sAtom[pos1].dtz0;
        PMEFloat tz1                        = sAtom[pos1].tz1;
        PMEFloat dtz1                       = sAtom[pos1].dtz1;
        PMEFloat tz2                        = sAtom[pos1].tz2;
        PMEFloat dtz2                       = sAtom[pos1].dtz2;
        PMEFloat tz3                        = sAtom[pos1].tz3;
        PMEFloat dtz3                       = sAtom[pos1].dtz3;
        PMEFloat dwxwy0                     = dwx * wy0;
        PMEFloat wxdwy0                     = wx  * dwy0;
        PMEFloat wxwy0                      = wx  * wy0;
        PMEFloat qterm0a                    = qterm0 * tz0;
        PMEFloat qterm1a                    = qterm1 * tz1;
        PMEFloat qterm2a                    = qterm2 * tz2;
        PMEFloat qterm3a                    = qterm3 * tz3;
        PMEFloat f1                         = -qterm0a * dwxwy0;
        PMEFloat f2                         = -qterm0a * wxdwy0;
        PMEFloat f3                         = -qterm0  * wxwy0  * dtz0;  
        f1                                 -=  qterm1a * dwxwy0;
        f2                                 -=  qterm1a * wxdwy0;
        f3                                 -=  qterm1  * wxwy0  * dtz1;
        f1                                 -=  qterm2a * dwxwy0;
        f2                                 -=  qterm2a * wxdwy0;
        f3                                 -=  qterm2  * wxwy0  * dtz2;   
        f1                                 -=  qterm3a * dwxwy0;
        f2                                 -=  qterm3a * wxdwy0;
        f3                                 -=  qterm3  * wxwy0  * dtz3;
        PMEFloat dwxwy1                     = dwx * wy1;
        PMEFloat wxdwy1                     = wx  * dwy1;
        PMEFloat wxwy1                      = wx  * wy1;
        PMEFloat qterm4a                    = qterm4 * tz0;
        PMEFloat qterm5a                    = qterm5 * tz1;
        PMEFloat qterm6a                    = qterm6 * tz2;
        PMEFloat qterm7a                    = qterm7 * tz3;            
        f1                                 -=  qterm4a * dwxwy1;
        f2                                 -=  qterm4a * wxdwy1;
        f3                                 -=  qterm4  * wxwy1  * dtz0;  
        f1                                 -=  qterm5a * dwxwy1;
        f2                                 -=  qterm5a * wxdwy1;
        f3                                 -=  qterm5  * wxwy1  * dtz1;
        f1                                 -=  qterm6a * dwxwy1;
        f2                                 -=  qterm6a * wxdwy1;
        f3                                 -=  qterm6  * wxwy1  * dtz2;   
        f1                                 -=  qterm7a * dwxwy1;
        f2                                 -=  qterm7a * wxdwy1;
        f3                                 -=  qterm7  * wxwy1  * dtz3;
        psFactor[tgx].dfx                   = f1;
        psFactor[tgx].dfy                   = f2;
        psFactor[tgx].dfz                   = f3;
                                 
 
      
        // Reduce within warps and save result
        if (tgx < 4)
        {
            psFactor[tgx].dfx              += psFactor[tgx + 4].dfx;
            psFactor[tgx].dfy              += psFactor[tgx + 4].dfy;
            psFactor[tgx].dfz              += psFactor[tgx + 4].dfz;
        }
            
        if (tgx < 2)
        {
            psFactor[tgx].dfx              += psFactor[tgx + 2].dfx;
            psFactor[tgx].dfy              += psFactor[tgx + 2].dfy;
            psFactor[tgx].dfz              += psFactor[tgx + 2].dfz;                               
        }
          
        if (tgx == 0)
        {
            sAtom[pos1].fx                  = psFactor->dfx + psFactor[1].dfx;
            sAtom[pos1].fy                  = psFactor->dfy + psFactor[1].dfy;
            sAtom[pos1].fz                  = psFactor->dfz + psFactor[1].dfz;
        }      
        
        pos1                               += GRADSUMTHREADS / GRADSUMGRID;
    }
    __syncthreads();
        
    // Set up one store for results
    if (threadIdx.x < lastAtom)
    {
        PMEFloat charge                     = sAtom[threadIdx.x].charge;
        PMEFloat fx                         = cSim.nfft1 * charge * sAtom[threadIdx.x].fx;
        PMEFloat fy                         = cSim.nfft2 * charge * sAtom[threadIdx.x].fy;
        PMEFloat fz                         = cSim.nfft3 * charge * sAtom[threadIdx.x].fz;

#ifdef MPI 
#ifdef PME_VIRIAL            
        fz                                  = sRecipf[6] * fx + sRecipf[7] * fy + sRecipf[8] * fz;
        fy                                  = sRecipf[3] * fx + sRecipf[4] * fy;
        fx                                  = sRecipf[0] * fx;
#else
        fz                                  = cSim.recipf[2][0] * fx + cSim.recipf[2][1] * fy + cSim.recipf[2][2] * fz;
        fy                                  = cSim.recipf[1][0] * fx + cSim.recipf[1][1] * fy;
        fx                                  = cSim.recipf[0][0] * fx;
#endif
        sPMEForce[threadIdx.x * 3]          = fx;
        sPMEForce[threadIdx.x * 3 + 1]      = fy;
        sPMEForce[threadIdx.x * 3 + 2]      = fz;
#ifdef PME_VIRIAL            
        unsigned int pos1                   = pos + threadIdx.x;
        cSim.pNBForceX[pos1]                = fx;
        cSim.pNBForceY[pos1]                = fy;
        cSim.pNBForceZ[pos1]                = fz;
#endif            
    }
    __syncthreads();
        
        
    
#if (__CUDA_ARCH__ >= 200)
    pos1                                    = pos * 3 + threadIdx.x;
    cSim.pPMEForce[pos1]                    = sPMEForce[threadIdx.x];
    if (threadIdx.x < 32)
    {
        pos1                               += 64;
        cSim.pPMEForce[pos1]                = sPMEForce[threadIdx.x + 64];
    }
#else
    if (threadIdx.x < 48)
    { 
        pos1                                = pos * 3 + threadIdx.x;
        cSim.pPMEForce[pos1]                = sPMEForce[threadIdx.x];
    }
#endif
#else           
        unsigned int pos1                   = pos + threadIdx.x; 
#ifdef PME_VIRIAL            
        cSim.pForceXBuffer[pos1]            = sRecipf[0] * fx;
        cSim.pForceYBuffer[pos1]            = sRecipf[3] * fx + sRecipf[4] * fy;
        cSim.pForceZBuffer[pos1]            = sRecipf[6] * fx + sRecipf[7] * fy + sRecipf[8] * fz;
#else
        cSim.pForceXBuffer[pos1]            = cSim.recipf[0][0] * fx;
        cSim.pForceYBuffer[pos1]            = cSim.recipf[1][0] * fx + cSim.recipf[1][1] * fy;
        cSim.pForceZBuffer[pos1]            = cSim.recipf[2][0] * fx + cSim.recipf[2][1] * fy + cSim.recipf[2][2] * fz;
#endif
    }
#endif

}
