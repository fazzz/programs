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
#ifdef PME_NTP
    __shared__ volatile PMEDouble sRecip[9];
#ifdef PME_ORTHOGONAL
    __shared__ volatile PMEDouble sUcell[9];
    if (threadIdx.x < 9)
        sUcell[threadIdx.x]                     = cSim.pNTPData->ucell[threadIdx.x];    
#endif    
    if (threadIdx.x < 9)
        sRecip[threadIdx.x]                     = cSim.pNTPData->recip[threadIdx.x];
    __syncthreads();
#endif
    unsigned int pos                            = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int increment                      = gridDim.x * blockDim.x;   
    while (pos < cSim.atoms)
    {
        double x                                = cSim.pImageX[pos];
        double y                                = cSim.pImageY[pos];
        double z                                = cSim.pImageZ[pos];
        unsigned int cellID                     = cSim.pImageCellID[pos];
#if defined(PME_NTP)
#if defined(PME_ORTHOGONAL)
        double dfx                              = sRecip[0] * x;
        double dfy                              = sRecip[4] * y;
        double dfz                              = sRecip[8] * z;
#else
        double dfx                              = sRecip[0] * x + sRecip[3] * y + sRecip[6] * z;
        double dfy                              =                 sRecip[4] * y + sRecip[7] * z;
        double dfz                              =                                 sRecip[8] * z; 
#endif
#else
#if defined(PME_ORTHOGONAL)
        double dfx                              = cSim.recip[0][0] * x;
        double dfy                              = cSim.recip[1][1] * y;
        double dfz                              = cSim.recip[2][2] * z;
#else
        double dfx                              = cSim.recip[0][0] * x + cSim.recip[1][0] * y + cSim.recip[2][0] * z;
        double dfy                              =                        cSim.recip[1][1] * y + cSim.recip[2][1] * z;
        double dfz                              =                                               cSim.recip[2][2] * z;            
#endif
#endif

        // Account for minimum image convention  
        dfx                                     = (dfx - round(dfx) + 0.5);
        dfy                                     = (dfy - round(dfy) + 0.5);
        dfz                                     = (dfz - round(dfz) + 0.5);
        dfx                                     = (dfx < 1.0 ? dfx : 0.0);
        dfy                                     = (dfy < 1.0 ? dfy : 0.0);
        dfz                                     = (dfz < 1.0 ? dfz : 0.0);
        
        // Can measure relative to cell edges for 3 or more cells on each axis, otherwise measure relative to cell center
        double xCell                            = (double)(cellID & CELLIDMASK) * cSim.oneOverXcells;
        double yCell                            = (double)((cellID >> CELLIDYSHIFT) & CELLIDMASK) * cSim.oneOverYcells;
        double zCell                            = (double)((cellID >> CELLIDZSHIFT) & CELLIDMASK) * cSim.oneOverZcells;

        x                                       = dfx - xCell;
        y                                       = dfy - yCell;
        z                                       = dfz - zCell;  
        if (x > 0.5)
            x                                  -= 1.0;
        if (x <= -0.5)
            x                                  += 1.0;
        if (y > 0.5)
            y                                  -= 1.0;
        if (y <= -0.5)
            y                                  += 1.0;            
        if (z > 0.5)
            z                                  -= 1.0;
        if (z <= -0.5)
            z                                  += 1.0;   
#ifdef PME_ORTHOGONAL
#ifdef PME_NTP
        PMEFloat2 xy                            = {(PMEFloat)(sUcell[0] * x), (PMEFloat)(sUcell[4] * y)};
        cSim.pAtomXYSP[pos]                     = xy;
        cSim.pAtomZSP[pos]                      = (PMEFloat)(sUcell[8] * z);
#else                 
        PMEFloat2 xy                            = {(PMEFloat)(cSim.a * x), (PMEFloat)(cSim.b * y)};
        cSim.pAtomXYSP[pos]                     = xy;
        cSim.pAtomZSP[pos]                      = (PMEFloat)(cSim.c * z);
#endif
#else
        PMEFloat2 xy                            = {x, y};
        cSim.pAtomXYSP[pos]                     = xy;
        cSim.pAtomZSP[pos]                      = z;
#endif

       pos                                    += increment;              
    }
}
