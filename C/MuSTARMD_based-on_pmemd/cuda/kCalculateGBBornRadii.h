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
struct Atom {
    PMEFloat x;
    PMEFloat y;
    PMEFloat z;
    PMEFloat r;
    PMEFloat r1i;
    PMEFloat s;
    PMEFloat s2;
};
#if (__CUDA_ARCH__ >= 200)
volatile __shared__ Atom sA[SM_2X_GBBORNRADII_THREADS_PER_BLOCK];
volatile __shared__ PMEDouble sReff[SM_2X_GBBORNRADII_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_2X_GBBORNRADII_THREADS_PER_BLOCK / GRID];
#else
#ifdef GB_IGB78
volatile __shared__ Atom sA[SM_13_GBBORNRADIIIGB78_THREADS_PER_BLOCK];
volatile __shared__ PMEDouble sReff[SM_13_GBBORNRADIIIGB78_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_13_GBBORNRADIIIGB78_THREADS_PER_BLOCK / GRID];
#else
volatile __shared__ Atom sA[SM_13_GBBORNRADII_THREADS_PER_BLOCK];
volatile __shared__ PMEDouble sReff[SM_13_GBBORNRADII_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_13_GBBORNRADII_THREADS_PER_BLOCK / GRID];
#endif
#endif
#ifdef GB_IGB78
__shared__ PMEFloat2 sNeckMaxValPos[21 * 21];
#endif
volatile __shared__ unsigned int sNext[GRID];

#ifdef GB_IGB78
    // Read neck lookup table
    unsigned int pos                        = threadIdx.x;
    while (pos < 21 * 21)
    {
        sNeckMaxValPos[pos]                 = cSim.pNeckMaxValPos[pos];
        pos                                += blockDim.x;
    }
#endif

    if (threadIdx.x < GRID)
    {
        sNext[threadIdx.x]                  = (threadIdx.x + 1) & (GRID - 1);
    }
    __syncthreads();

    // Initialize queue position
    volatile unsigned int* psPos            = &sPos[threadIdx.x >> GRIDBITS];
    *psPos                                  = (blockIdx.x * blockDim.x + threadIdx.x) >> GRIDBITS;
   
    while (*psPos < cSim.workUnits)
    {  
        // Extract cell coordinates from appropriate work unit
        unsigned int x                      = cSim.pWorkUnit[*psPos];
        unsigned int y                      = ((x >> 2) & 0x7fff) << GRIDBITS;
        x                                   = (x >> 17) << GRIDBITS;
        unsigned int tgx                    = threadIdx.x & (GRID - 1);
        unsigned int i                      = x + tgx;
        PMEFloat2 xyi                       = cSim.pAtomXYSP[i];
        PMEFloat zi                         = cSim.pAtomZSP[i];
        PMEFloat ri                         = cSim.pAtomRBorn[i];
        PMEFloat si                         = cSim.pAtomS[i];
        PMEFloat si2                        = si * si;
        PMEDouble reff_i                    = (PMEDouble)0.0;
        unsigned int tbx                    = threadIdx.x - tgx;
        volatile Atom* psA                  = &sA[tbx]; 
       
        if (x == y) // Handle diagonals uniquely at 50% efficiency, skipping i == j interactions
        { 
            PMEFloat xi                     = xyi.x;
            PMEFloat yi                     = xyi.y;
            sA[threadIdx.x].x               = xi;
            sA[threadIdx.x].y               = yi;
            sA[threadIdx.x].z               = zi;
            ri                             -= cSim.offset;
            PMEFloat ri1i                   = (PMEFloat)1.0 / ri;
            sA[threadIdx.x].r               = ri;
            sA[threadIdx.x].s               = si;
            sA[threadIdx.x].s2              = si2;
            sA[threadIdx.x].r1i             = ri1i;

            for (unsigned int j = sNext[tgx]; j != tgx; j = sNext[j])
            {
                PMEFloat xij                = xi - psA[j].x; 
                PMEFloat yij                = yi - psA[j].y; 
                PMEFloat zij                = zi - psA[j].z;
                PMEFloat r2                 = xij * xij + yij * yij + zij * zij;
                PMEFloat dij                = sqrt(r2);  
                       
                PMEFloat sj                 = psA[j].s;       
                if (dij < cSim.rgbmax + sj)
                {
                    PMEFloat dij1i          = (PMEFloat)1.0 / dij;
                    PMEFloat dij2i          = dij1i * dij1i;
                    PMEFloat dr;
               
                    if (dij > cSim.rgbmax - sj)
                    {
                        PMEFloat uij        = (PMEFloat)1.0 / (dij - sj);
                        dr                  = (PMEFloat)0.125 * dij1i * ((PMEFloat)1.0 + (PMEFloat)2.0 * dij * uij + 
                                              cSim.rgbmax2i * (r2 - (PMEFloat)4.0 * cSim.rgbmax * dij - psA[j].s2) + 
                                              (PMEFloat)2.0 * log((dij - sj) * cSim.rgbmax1i));
                    }
                    else if (dij > (PMEFloat)4.0 * sj)
                    {            
                        PMEFloat tmpsd      = psA[j].s2 * dij2i;
                        PMEFloat dumbo      = ta + tmpsd *  (tb + tmpsd * (tc + tmpsd * (td + tmpsd * tdd)));
                        dr                  = tmpsd * sj * dij2i * dumbo;
                    }
                    else
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (dij + sj);                 
                        if (dij > ri + sj)
                        {
                            PMEFloat v4     = log(v2 * (dij - sj));
                            dr              = (PMEFloat)0.5 * (sj / (r2 - psA[j].s2) + (PMEFloat)0.5 * dij1i * v4);
                        }
                        else if (dij > fabs(ri - sj))
                        {
                            PMEFloat v4     = log(v2 * ri);
                            PMEFloat theta  = (PMEFloat)0.5 * ri1i * dij1i * (r2 + ri * ri - psA[j].s2);
                            dr              = (PMEFloat)0.25 * (ri1i * ((PMEFloat)2.0 - theta) - 
                                              v2 + dij1i * v4);
                     
                        }        
                        else if (ri < sj)
                        {
                            PMEFloat v4     = log(v2 * (sj - dij));
                            dr              = (PMEFloat)0.5 * (sj / (r2 - psA[j].s2) + (PMEFloat)2.0 * ri1i + 
                                              (PMEFloat)0.5 * dij1i * v4);
                        }
                    }                    
#ifdef GB_IGB78
                    if (dij < cSim.gb_neckcut + ri + psA[j].r)
                    {
                        
                        unsigned int ii     = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        unsigned int jj     = round((psA[j].r - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        PMEFloat2 neckValPos = sNeckMaxValPos[ii * 21 + jj];                      
                        PMEFloat mdist      = dij - neckValPos.y;
                        PMEFloat mdist2     = mdist * mdist;
                        PMEFloat mdist6     = mdist2 * mdist;
                        mdist6              = mdist6 * mdist6;
                        PMEFloat neck       = neckValPos.x / ((PMEFloat)1.0 + mdist2 + (PMEFloat)0.3 * mdist6);
                        dr                 += cSim.gb_neckscale * neck;
                    }                    
#endif                              
                    reff_i                 -= (PMEDouble)dr;
                }            
            }

            int offset                      = x + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pReffBuffer[offset]        = reff_i;
        }
        else
        {
            int j                           = y + tgx;            
            PMEFloat2 xyj                   = cSim.pAtomXYSP[j];
            sA[threadIdx.x].s               = cSim.pAtomS[j];  
            sA[threadIdx.x].r               = cSim.pAtomRBorn[j];          
            sA[threadIdx.x].z               = cSim.pAtomZSP[j];
            sReff[threadIdx.x]              = (PMEFloat)0.0;
            PMEFloat xi                     = xyi.x;
            PMEFloat yi                     = xyi.y;           
            ri                             -= cSim.offset;
            PMEFloat ri1i                   = (PMEFloat)1.0 / ri;
            sA[threadIdx.x].x               = xyj.x;
            sA[threadIdx.x].y               = xyj.y;
            sA[threadIdx.x].s2              = sA[threadIdx.x].s * sA[threadIdx.x].s;
            sA[threadIdx.x].r              -= cSim.offset;
            sA[threadIdx.x].r1i             = (PMEFloat)1.0 / sA[threadIdx.x].r;
            volatile PMEDouble* psReff      = &sReff[tbx];
            j                               = tgx;
            do
            {
                PMEFloat xij                = xi - psA[j].x; 
                PMEFloat yij                = yi - psA[j].y; 
                PMEFloat zij                = zi - psA[j].z;
                PMEFloat r2                 = xij * xij + yij * yij + zij * zij;
                PMEFloat dij                = sqrt(r2);
                PMEFloat dij1i              = (PMEFloat)1.0 / dij;
                PMEFloat dij2i              = dij1i * dij1i;


                PMEFloat sj                 = psA[j].s;
                if (dij < cSim.rgbmax + sj)
                {
                    PMEFloat dr;
               
                    if (dij > cSim.rgbmax - sj)
                    {
                        PMEFloat uij        = (PMEFloat)1.0 / (dij - sj);
                        dr                  = (PMEFloat)0.125 * dij1i * ((PMEFloat)1.0 + (PMEFloat)2.0 * dij * uij + 
                                              cSim.rgbmax2i * (r2 - (PMEFloat)4.0 * cSim.rgbmax * dij - psA[j].s2) + 
                                              (PMEFloat)2.0 * log((dij - sj) * cSim.rgbmax1i));
                    }
                    else if (dij > (PMEFloat)4.0 * sj)
                    {            
                        PMEFloat tmpsd      = psA[j].s2 * dij2i;
                        PMEFloat dumbo      = ta + tmpsd *  (tb + tmpsd * (tc + tmpsd * (td + tmpsd * tdd)));
                        dr                  = tmpsd * sj * dij2i * dumbo;
                    }
                    else
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (dij + sj);                 
                        if (dij > ri + sj)
                        {
                            PMEFloat v4     = log(v2 * (dij - sj));
                            dr              = (PMEFloat)0.5 * (sj / (r2 - psA[j].s2) + (PMEFloat)0.5 * dij1i * v4);
                        }
                        else if (dij > fabs(ri - sj))
                        {
                            PMEFloat v4     = log(v2 * ri);
                            PMEFloat theta  = (PMEFloat)0.5 * ri1i * dij1i * (r2 + ri * ri - psA[j].s2);
                            dr              = (PMEFloat)0.25 * (ri1i * ((PMEFloat)2.0 - theta) - 
                                              v2 + dij1i * v4);
                     
                        }        
                        else if (ri < sj)
                        {
                            PMEFloat v4     = log(v2 * (sj - dij));
                            dr              = (PMEFloat)0.5 * (sj / (r2 - psA[j].s2) + (PMEFloat)2.0 * ri1i + 
                                              (PMEFloat)0.5 * dij1i * v4);
                        }
                    }
#ifdef GB_IGB78
                    if (dij < cSim.gb_neckcut + ri + psA[j].r)
                    {
                        
                        unsigned int ii     = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        unsigned int jj     = round((psA[j].r - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        PMEFloat2 neckValPos = sNeckMaxValPos[ii * 21 + jj];
                        PMEFloat mdist      = dij - neckValPos.y;
                        PMEFloat mdist2     = mdist * mdist;
                        PMEFloat mdist6     = mdist2 * mdist;
                        mdist6              = mdist6 * mdist6;
                        PMEFloat neck       = neckValPos.x / ((PMEFloat)1.0 + mdist2 + (PMEFloat)0.3 * mdist6);
                        dr                 += cSim.gb_neckscale * neck;
                    }                    
#endif                           
                    reff_i                 -= (PMEDouble)dr;                    
                }

                if (dij < cSim.rgbmax + si)
                {
                    PMEFloat dr;
               
                    if (dij > cSim.rgbmax - si)
                    {
                        PMEFloat uij        = (PMEFloat)1.0 / (dij - si);
                        dr                  = (PMEFloat)0.125 * dij1i * ((PMEFloat)1.0 + (PMEFloat)2.0 * dij * uij + 
                                              cSim.rgbmax2i * (r2 - (PMEFloat)4.0 * cSim.rgbmax * dij - si2) + 
                                              (PMEFloat)2.0 * log((dij - si) * cSim.rgbmax1i));
                    }
                    else if (dij > (PMEFloat)4.0 * si)
                    {            
                        PMEFloat tmpsd      = si2 * dij2i;
                        PMEFloat dumbo      = ta + tmpsd *  (tb + tmpsd * (tc + tmpsd * (td + tmpsd * tdd)));
                        dr                  = tmpsd * si * dij2i * dumbo;
                    }
                    else
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (dij + si);                 
                        if (dij > psA[j].r + si)
                        {
                            PMEFloat v4     = log(v2 * (dij - si));
                            dr              = (PMEFloat)0.5 * (si / (r2 - si2) + (PMEFloat)0.5 * dij1i * v4);
                        }
                        else if (dij > fabs(psA[j].r - si))
                        {
                            PMEFloat v4     = log(v2 * psA[j].r);
                            PMEFloat theta  = (PMEFloat)0.5 * psA[j].r1i * dij1i * (r2 + psA[j].r * psA[j].r - si2);
                            dr              = (PMEFloat)0.25 * (psA[j].r1i * ((PMEFloat)2.0 - theta) - 
                                              v2 + dij1i * v4);
                     
                        }        
                        else if (psA[j].r < si)
                        {
                            PMEFloat v4     = log(v2 * (si - dij));
                            dr              = (PMEFloat)0.5 * (si / (r2 - si2) + (PMEFloat)2.0 * psA[j].r1i + 
                                              (PMEFloat)0.5 * dij1i * v4);
                        }
                    }
#ifdef GB_IGB78
                    if (dij < cSim.gb_neckcut + ri + psA[j].r)
                    {
                        
                        unsigned int ii     = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        unsigned int jj     = round((psA[j].r - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        PMEFloat2 neckValPos = sNeckMaxValPos[jj * 21 + ii];
                        PMEFloat mdist      = dij - neckValPos.y;
                        PMEFloat mdist2     = mdist * mdist;
                        PMEFloat mdist6     = mdist2 * mdist;
                        mdist6              = mdist6 * mdist6;
                        PMEFloat neck       = neckValPos.x / ((PMEFloat)1.0 + mdist2 + (PMEFloat)0.3 * mdist6);
                        dr                 += cSim.gb_neckscale * neck;
                    }                    
#endif                           
                    psReff[j]              -= (PMEDouble)dr;
                }                
                j                           = sNext[j];            
            }
            while (j != tgx);
       
            int offset                      = x + tgx + (y >> GRIDBITS) * cSim.stride;
            cSim.pReffBuffer[offset]        = reff_i;
            offset                          = y + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pReffBuffer[offset]        = sReff[threadIdx.x];
        }
        if (tgx == 0)
            *psPos                          = atomicAdd(cSim.pGBBRPosition, 1);
    }
}

