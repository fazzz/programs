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

#if (__CUDA_ARCH__ >= 300)
#define PSATOMX(i)   shAtom.x
#define PSATOMY(i)   shAtom.y
#define PSATOMZ(i)   shAtom.z
#define PSATOMR(i)   shAtom.r
#define PSATOMR1I(i) shAtom.r1i
#define PSATOMS(i)   shAtom.s
#define PSATOMS2(i)  shAtom.s2
#else
#define PSATOMX(i)   psA[i].x
#define PSATOMY(i)   psA[i].y
#define PSATOMZ(i)   psA[i].z
#define PSATOMR(i)   psA[i].r
#define PSATOMR1I(i) psA[i].r1i
#define PSATOMS(i)   psA[i].s
#define PSATOMS2(i)  psA[i].s2
#endif

#if (__CUDA_ARCH__ >= 300)
Atom shAtom;
volatile __shared__ PMEForce sReff[SM_3X_GBBORNRADII_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_3X_GBBORNRADII_THREADS_PER_BLOCK / GRID];
#elif (__CUDA_ARCH__ >= 200)
volatile __shared__ Atom sA[SM_2X_GBBORNRADII_THREADS_PER_BLOCK];
volatile __shared__ PMEForce sReff[SM_2X_GBBORNRADII_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_2X_GBBORNRADII_THREADS_PER_BLOCK / GRID];
#else
#ifdef GB_IGB78
volatile __shared__ Atom sA[SM_13_GBBORNRADIIIGB78_THREADS_PER_BLOCK];
volatile __shared__ PMEForce sReff[SM_13_GBBORNRADIIIGB78_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_13_GBBORNRADIIIGB78_THREADS_PER_BLOCK / GRID];
#else
volatile __shared__ Atom sA[SM_13_GBBORNRADII_THREADS_PER_BLOCK];
volatile __shared__ PMEForce sReff[SM_13_GBBORNRADII_THREADS_PER_BLOCK];
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
        PMEForce reff_i                     = (PMEForce)0;
        unsigned int tbx                    = threadIdx.x - tgx;
#if (__CUDA_ARCH__ >= 300)
        unsigned int shIdx                  = sNext[tgx];
#else
        volatile Atom* psA                  = &sA[tbx];
#endif
        volatile PMEForce* psReff          = &sReff[tbx];
       
        if (x == y) // Handle diagonals uniquely at 50% efficiency, skipping i == j interactions
        { 
            PMEFloat xi                     = xyi.x;
            PMEFloat yi                     = xyi.y;
            PSATOMX(tgx)                    = xi;
            PSATOMY(tgx)                    = yi;
            PSATOMZ(tgx)                    = zi;
            ri                             -= cSim.offset;
            PMEFloat ri1i                   = (PMEFloat)1.0 / ri;
            PSATOMR(tgx)                    = ri;
            PSATOMS(tgx)                    = si;
            PSATOMS2(tgx)                   = si2;
            PSATOMR1I(tgx)                  = ri1i;
#if (__CUDA_ARCH__ >= 300)
            shAtom.x                        = __shfl(shAtom.x, shIdx);
            shAtom.y                        = __shfl(shAtom.y, shIdx);
            shAtom.z                        = __shfl(shAtom.z, shIdx);
            shAtom.r                        = __shfl(shAtom.r, shIdx);
            shAtom.s                        = __shfl(shAtom.s, shIdx);
            shAtom.s2                       = __shfl(shAtom.s2, shIdx);
            shAtom.r1i                      = __shfl(shAtom.r1i, shIdx);
#endif
            for (unsigned int j = sNext[tgx]; j != tgx; j = sNext[j])
            {
                PMEFloat xij                = xi - PSATOMX(j);
                PMEFloat yij                = yi - PSATOMY(j);
                PMEFloat zij                = zi - PSATOMZ(j);
                PMEFloat r2                 = xij * xij + yij * yij + zij * zij;
                PMEFloat dij                = sqrt(r2);  
                       
                PMEFloat sj                 = PSATOMS(j);     
                if (dij < cSim.rgbmax + sj)
                {
                    PMEFloat dij1i          = (PMEFloat)1.0 / dij;
                    PMEFloat dij2i          = dij1i * dij1i;
                    PMEFloat dr;
               
                    if (dij > cSim.rgbmax - sj)
                    {
                        PMEFloat uij        = (PMEFloat)1.0 / (dij - sj);
                        dr                  = (PMEFloat)0.125 * dij1i * ((PMEFloat)1.0 + (PMEFloat)2.0 * dij * uij + 
                                              cSim.rgbmax2i * (r2 - (PMEFloat)4.0 * cSim.rgbmax * dij - PSATOMS2(j)) + 
                                              (PMEFloat)2.0 * log((dij - sj) * cSim.rgbmax1i));
                    }
                    else if (dij > (PMEFloat)4.0 * sj)
                    {            
                        PMEFloat tmpsd      = PSATOMS2(j) * dij2i;
                        PMEFloat dumbo      = ta + tmpsd *  (tb + tmpsd * (tc + tmpsd * (td + tmpsd * tdd)));
                        dr                  = tmpsd * sj * dij2i * dumbo;
                    }
                    else
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (dij + sj);                 
                        if (dij > ri + sj)
                        {
                            PMEFloat v4     = log(v2 * (dij - sj));
                            dr              = (PMEFloat)0.5 * (sj / (r2 - PSATOMS2(j)) + (PMEFloat)0.5 * dij1i * v4);
                        }
                        else if (dij > fabs(ri - sj))
                        {
                            PMEFloat v4     = log(v2 * ri);
                            PMEFloat theta  = (PMEFloat)0.5 * ri1i * dij1i * (r2 + ri * ri - PSATOMS2(j));
                            dr              = (PMEFloat)0.25 * (ri1i * ((PMEFloat)2.0 - theta) - 
                                              v2 + dij1i * v4);
                     
                        }        
                        else if (ri < sj)
                        {
                            PMEFloat v4     = log(v2 * (sj - dij));
                            dr              = (PMEFloat)0.5 * (sj / (r2 - PSATOMS2(j)) + (PMEFloat)2.0 * ri1i + 
                                              (PMEFloat)0.5 * dij1i * v4);
                        }
                    }                    
#ifdef GB_IGB78
                    if (dij < cSim.gb_neckcut + ri + PSATOMR(j))
                    {
                        
                        unsigned int ii     = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        unsigned int jj     = round((PSATOMR(j) - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        PMEFloat2 neckValPos = sNeckMaxValPos[ii * 21 + jj];                      
                        PMEFloat mdist      = dij - neckValPos.y;
                        PMEFloat mdist2     = mdist * mdist;
                        PMEFloat mdist6     = mdist2 * mdist;
                        mdist6              = mdist6 * mdist6;
                        PMEFloat neck       = neckValPos.x / ((PMEFloat)1.0 + mdist2 + (PMEFloat)0.3 * mdist6);
                        dr                 += cSim.gb_neckscale * neck;
                    }                    
#endif                          
#ifdef use_SPFP
                    reff_i                 -= lliroundf(FORCESCALEF * dr);
#else    
                    reff_i                 -= (PMEForce)dr;
#endif
                }

#if (__CUDA_ARCH__ >= 300)
                shAtom.x                    = __shfl(shAtom.x, shIdx);
                shAtom.y                    = __shfl(shAtom.y, shIdx);
                shAtom.z                    = __shfl(shAtom.z, shIdx);
                shAtom.r                    = __shfl(shAtom.r, shIdx);
                shAtom.s                    = __shfl(shAtom.s, shIdx);
                shAtom.s2                   = __shfl(shAtom.s2, shIdx);
                shAtom.r1i                  = __shfl(shAtom.r1i, shIdx);
#endif  
            }
#ifdef use_SPFP
            int offset                      = x + tgx;
            atomicAdd((unsigned long long int*)&cSim.pReffAccumulator[offset], llitoulli(reff_i));
#else
            int offset                      = x + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pReffBuffer[offset]        = reff_i;
#endif
        }
        else
        {
            int j                           = y + tgx;            
            PMEFloat2 xyj                   = cSim.pAtomXYSP[j];
            PSATOMS(tgx)                    = cSim.pAtomS[j];  
            PSATOMR(tgx)                    = cSim.pAtomRBorn[j];          
            PSATOMZ(tgx)                    = cSim.pAtomZSP[j];
            psReff[tgx]                     = (PMEForce)0;
            PMEFloat xi                     = xyi.x;
            PMEFloat yi                     = xyi.y;           
            ri                             -= cSim.offset;
            PMEFloat ri1i                   = (PMEFloat)1.0 / ri;
            PSATOMX(tgx)                    = xyj.x;
            PSATOMY(tgx)                    = xyj.y;
            PSATOMS2(tgx)                   = PSATOMS(tgx) * PSATOMS(tgx);
            PSATOMR(tgx)                   -= cSim.offset;
            PSATOMR1I(tgx)                  = (PMEFloat)1.0 / PSATOMR(tgx);
            j                               = tgx;
            do
            {
                PMEFloat xij                = xi - PSATOMX(j);
                PMEFloat yij                = yi - PSATOMY(j);
                PMEFloat zij                = zi - PSATOMZ(j);
                PMEFloat r2                 = xij * xij + yij * yij + zij * zij;
                PMEFloat dij                = sqrt(r2);
                PMEFloat dij1i              = (PMEFloat)1.0 / dij;
                PMEFloat dij2i              = dij1i * dij1i;


                PMEFloat sj                 = PSATOMS(j);
                if (dij < cSim.rgbmax + sj)
                {
                    PMEFloat dr;
               
                    if (dij > cSim.rgbmax - sj)
                    {
                        PMEFloat uij        = (PMEFloat)1.0 / (dij - sj);
                        dr                  = (PMEFloat)0.125 * dij1i * ((PMEFloat)1.0 + (PMEFloat)2.0 * dij * uij + 
                                              cSim.rgbmax2i * (r2 - (PMEFloat)4.0 * cSim.rgbmax * dij - PSATOMS2(j)) + 
                                              (PMEFloat)2.0 * log((dij - sj) * cSim.rgbmax1i));
                    }
                    else if (dij > (PMEFloat)4.0 * sj)
                    {            
                        PMEFloat tmpsd      = PSATOMS2(j) * dij2i;
                        PMEFloat dumbo      = ta + tmpsd *  (tb + tmpsd * (tc + tmpsd * (td + tmpsd * tdd)));
                        dr                  = tmpsd * sj * dij2i * dumbo;
                    }
                    else
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (dij + sj);                 
                        if (dij > ri + sj)
                        {
                            PMEFloat v4     = log(v2 * (dij - sj));
                            dr              = (PMEFloat)0.5 * (sj / (r2 - PSATOMS2(j)) + (PMEFloat)0.5 * dij1i * v4);
                        }
                        else if (dij > fabs(ri - sj))
                        {
                            PMEFloat v4     = log(v2 * ri);
                            PMEFloat theta  = (PMEFloat)0.5 * ri1i * dij1i * (r2 + ri * ri - PSATOMS2(j));
                            dr              = (PMEFloat)0.25 * (ri1i * ((PMEFloat)2.0 - theta) - 
                                              v2 + dij1i * v4);
                     
                        }        
                        else if (ri < sj)
                        {
                            PMEFloat v4     = log(v2 * (sj - dij));
                            dr              = (PMEFloat)0.5 * (sj / (r2 - PSATOMS2(j)) + (PMEFloat)2.0 * ri1i + 
                                              (PMEFloat)0.5 * dij1i * v4);
                        }
                    }
#ifdef GB_IGB78
                    if (dij < cSim.gb_neckcut + ri + PSATOMR(j))
                    {
                        
                        unsigned int ii     = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        unsigned int jj     = round((PSATOMR(j) - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        PMEFloat2 neckValPos = sNeckMaxValPos[ii * 21 + jj];
                        PMEFloat mdist      = dij - neckValPos.y;
                        PMEFloat mdist2     = mdist * mdist;
                        PMEFloat mdist6     = mdist2 * mdist;
                        mdist6              = mdist6 * mdist6;
                        PMEFloat neck       = neckValPos.x / ((PMEFloat)1.0 + mdist2 + (PMEFloat)0.3 * mdist6);
                        dr                 += cSim.gb_neckscale * neck;
                    }                    
#endif 
#ifdef use_SPFP
                    reff_i                 -= lliroundf(FORCESCALEF * dr);
#else                          
                    reff_i                 -= (PMEForce)dr; 
#endif                   
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
                        if (dij > PSATOMR(j) + si)
                        {
                            PMEFloat v4     = log(v2 * (dij - si));
                            dr              = (PMEFloat)0.5 * (si / (r2 - si2) + (PMEFloat)0.5 * dij1i * v4);
                        }
                        else if (dij > fabs(PSATOMR(j) - si))
                        {
                            PMEFloat v4     = log(v2 * PSATOMR(j));
                            PMEFloat theta  = (PMEFloat)0.5 * PSATOMR1I(j) * dij1i * (r2 + PSATOMR(j) * PSATOMR(j) - si2);
                            dr              = (PMEFloat)0.25 * (PSATOMR1I(j) * ((PMEFloat)2.0 - theta) - 
                                              v2 + dij1i * v4);
                     
                        }        
                        else if (PSATOMR(j) < si)
                        {
                            PMEFloat v4     = log(v2 * (si - dij));
                            dr              = (PMEFloat)0.5 * (si / (r2 - si2) + (PMEFloat)2.0 * PSATOMR1I(j) + 
                                              (PMEFloat)0.5 * dij1i * v4);
                        }
                    }
#ifdef GB_IGB78
                    if (dij < cSim.gb_neckcut + ri + PSATOMR(j))
                    {
                        
                        unsigned int ii     = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        unsigned int jj     = round((PSATOMR(j) - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        PMEFloat2 neckValPos = sNeckMaxValPos[jj * 21 + ii];
                        PMEFloat mdist      = dij - neckValPos.y;
                        PMEFloat mdist2     = mdist * mdist;
                        PMEFloat mdist6     = mdist2 * mdist;
                        mdist6              = mdist6 * mdist6;
                        PMEFloat neck       = neckValPos.x / ((PMEFloat)1.0 + mdist2 + (PMEFloat)0.3 * mdist6);
                        dr                 += cSim.gb_neckscale * neck;
                    }                    
#endif                           
#ifdef use_SPFP
                    psReff[j]              -= lliroundf(FORCESCALEF * dr);
#else
                    psReff[j]              -= (PMEForce)dr;
#endif
                }     

#if (__CUDA_ARCH__ >= 300)
                shAtom.x                    = __shfl(shAtom.x, shIdx);
                shAtom.y                    = __shfl(shAtom.y, shIdx);
                shAtom.z                    = __shfl(shAtom.z, shIdx);
                shAtom.r                    = __shfl(shAtom.r, shIdx);
                shAtom.s                    = __shfl(shAtom.s, shIdx);
                shAtom.s2                   = __shfl(shAtom.s2, shIdx);
                shAtom.r1i                  = __shfl(shAtom.r1i, shIdx);
#endif            
                j                           = sNext[j];            
            }
            while (j != tgx);
#ifdef use_SPFP
            int offset                      = x + tgx;
            atomicAdd((unsigned long long int*)&cSim.pReffAccumulator[offset], llitoulli(reff_i));
            offset                          = y + tgx;
            atomicAdd((unsigned long long int*)&cSim.pReffAccumulator[offset], llitoulli(psReff[tgx]));
#else       
            int offset                      = x + tgx + (y >> GRIDBITS) * cSim.stride;
            cSim.pReffBuffer[offset]        = reff_i;
            offset                          = y + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pReffBuffer[offset]        = psReff[tgx];
#endif
        }
        if (tgx == 0)
            *psPos                          = atomicAdd(cSim.pGBBRPosition, 1);
    }
}

