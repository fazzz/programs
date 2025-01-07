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
const PMEFloat ONEHALF                          = 0.5;
const PMEFloat ONETHIRD                         = (1.0 / 3.0);
const PMEFloat ONE                              = 1.0; 
const PMEFloat TWO                              = 2.0; 
const PMEFloat THREE                            = 3.0; 

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
#ifdef PME_SMALLBOX
        if (x > cSim.maxCellX)
            x                                  -= 1.0;
        if (x < cSim.minCellX)
            x                                  += 1.0;
        if (y > cSim.maxCellY)
            y                                  -= 1.0;
        if (y < cSim.minCellY)
            y                                  += 1.0;            
        if (z > cSim.maxCellZ)
            z                                  -= 1.0;
        if (z <= cSim.minCellZ)
            z                                  += 1.0;   
#else        
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
#endif
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
        PMEFloat fx                             = dfx;
        PMEFloat fy                             = dfy;
        PMEFloat fz                             = dfz;
        fx                                     *= cSim.nfft1;
        fy                                     *= cSim.nfft2;
        fz                                     *= cSim.nfft3;
        int ifx                                 = int(fx);
        int ify                                 = int(fy);
        int ifz                                 = int(fz);
        fx                                     -= ifx;
        fy                                     -= ify;
        fz                                     -= ifz;
        ifx                                    -= cSim.order - 1;
        ify                                    -= cSim.order - 1;
        ifz                                    -= cSim.order - 1;
        if (ifx < 0)
            ifx                                += cSim.nfft1;
        if (ify < 0)
            ify                                += cSim.nfft2;
        if (ifz < 0)
            ifz                                += cSim.nfft3;    
        cSim.pIFractX[pos]                      = ifx;
        cSim.pIFractY[pos]                      = ify;
        cSim.pIFractZ[pos]                      = ifz;
            
        // Calculate weights
        PMEFloat4 w, d;
            
        // Order 2
        w.x                                     = ONE - fx;
        w.y                                     = fx;
            
        // Order 3
        w.z                                     = ONEHALF * fx * w.y;
        w.y                                     = ONEHALF * ((fx + ONE) * w.x + (TWO - fx) * w.y);
        w.x                                     = ONEHALF * (ONE - fx)  * w.x;
            
        // Differences
        d.x                                     = -w.x;
        d.y                                     =  w.x - w.y;
        d.z                                     =  w.y - w.z;
        d.w                                     =  w.z;
        cSim.pDThetaX[pos]                      = d;         
            
        // Order 4
        w.w                                     = ONETHIRD * fx * w.z;
        w.z                                     = ONETHIRD * ((fx + ONE) * w.y + (THREE - fx) * w.z);
        w.y                                     = ONETHIRD * ((fx + TWO) * w.x + (TWO - fx) * w.y);
        w.x                                     = ONETHIRD *  (ONE - fx) * w.x;
           
        cSim.pThetaX[pos]                       = w;
          
                
        // Order 2
        w.x                                     = ONE - fy;
        w.y                                     = fy;

        // Order 3
        w.z                                     = ONEHALF * fy * w.y;
        w.y                                     = ONEHALF * ((fy + ONE) * w.x + (TWO - fy) * w.y);
        w.x                                     = ONEHALF * (ONE - fy)  * w.x;
         
        // Differences
        d.x                                     = -w.x;
        d.y                                     =  w.x - w.y;
        d.z                                     =  w.y - w.z;
        d.w                                     =  w.z;  

        cSim.pDThetaY[pos]                      = d;        

        // Order 4
        w.w                                     = ONETHIRD * fy * w.z;
        w.z                                     = ONETHIRD * ((fy + ONE) * w.y + (THREE - fy) * w.z);
        w.y                                     = ONETHIRD * ((fy + TWO) * w.x + (TWO - fy) * w.y);
        w.x                                     = ONETHIRD *  (ONE - fy) * w.x;
            
        cSim.pThetaY[pos]                       = w;
         
        // Order 2
        w.x                                     = ONE - fz;
        w.y                                     = fz;
            
        // Order 3
        w.z                                     = ONEHALF * fz * w.y;
        w.y                                     = ONEHALF * ((fz + ONE) * w.x + (TWO - fz) * w.y);
        w.x                                     = ONEHALF * (ONE - fz)  * w.x;
            
        // Differences
        d.x                                     = -w.x;
        d.y                                     =  w.x - w.y;
        d.z                                     =  w.y - w.z;
        d.w                                     =  w.z;  

        cSim.pDThetaZ[pos]                      = d;          
            
        // Order 4
        w.w                                     = ONETHIRD * fz * w.z;
        w.z                                     = ONETHIRD * ((fz + ONE) * w.y + (THREE - fz) * w.z);
        w.y                                     = ONETHIRD * ((fz + TWO) * w.x + (TWO - fz) * w.y);
        w.x                                     = ONETHIRD *  (ONE - fz) * w.x;
            
        cSim.pThetaZ[pos]                       = w;
            
        pos                                    += increment;              
    }
}
