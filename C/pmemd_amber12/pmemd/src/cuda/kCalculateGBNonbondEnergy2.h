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

struct NLForce {
    PMEForce x;
    PMEForce y;
    PMEForce z;
};

struct Atom {
    PMEFloat x;
    PMEFloat y;
    PMEFloat z;
    PMEFloat r;
    PMEFloat r1i;
    PMEFloat s;
    PMEFloat temp7;
};


#if (__CUDA_ARCH__ >= 300)
#define PSATOMX(i)     shAtom.x
#define PSATOMY(i)     shAtom.y
#define PSATOMZ(i)     shAtom.z
#define PSATOMR(i)     shAtom.r 
#define PSATOMR1I(i)   shAtom.r1i
#define PSATOMS(i)     shAtom.s
#define PSATOMTEMP7(i) shAtom.temp7
#define PSFX(i)        psF[i].x
#define PSFY(i)        psF[i].y
#define PSFZ(i)        psF[i].z
#else
#define PSATOMX(i)     psA[i].x
#define PSATOMY(i)     psA[i].y
#define PSATOMZ(i)     psA[i].z
#define PSATOMR(i)     psA[i].r
#define PSATOMR1I(i)   psA[i].r1i
#define PSATOMS(i)     psA[i].s
#define PSATOMTEMP7(i) psA[i].temp7
#define PSFX(i)        fx_j
#define PSFY(i)        fy_j
#define PSFZ(i)        fz_j
#endif

#if (__CUDA_ARCH__ >= 300)
Atom shAtom;
volatile __shared__ NLForce sForce[SM_3X_GBNONBONDENERGY2_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_3X_GBNONBONDENERGY2_THREADS_PER_BLOCK / GRID];
#elif (__CUDA_ARCH__ >= 200)
volatile __shared__ Atom sA[SM_2X_GBNONBONDENERGY2_THREADS_PER_BLOCK];
volatile __shared__ PMEForce sE[SM_2X_GBNONBONDENERGY2_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_2X_GBNONBONDENERGY2_THREADS_PER_BLOCK / GRID];
#else
#ifdef GB_IGB78
volatile __shared__ Atom sA[SM_13_GBNONBONDENERGY2IGB78_THREADS_PER_BLOCK];
volatile __shared__ PMEForce sE[SM_13_GBNONBONDENERGY2IGB78_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_13_GBNONBONDENERGY2IGB78_THREADS_PER_BLOCK / GRID];
#else
volatile __shared__ Atom sA[SM_13_GBNONBONDENERGY2_THREADS_PER_BLOCK];
volatile __shared__ PMEForce sE[SM_13_GBNONBONDENERGY2_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_13_GBNONBONDENERGY2_THREADS_PER_BLOCK / GRID];
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
        PMEFloat temp7_i                    = cSim.pTemp7[i];
        PMEForce fx_i                       = (PMEForce)0;
        PMEForce fy_i                       = (PMEForce)0;
        PMEForce fz_i                       = (PMEForce)0;
        unsigned int tbx                    = threadIdx.x - tgx;
#if (__CUDA_ARCH__ >= 300)
        volatile NLForce* psF               = &sForce[tbx];
        unsigned int shIdx                  = sNext[tgx];
#else
        volatile Atom* psA                  = &sA[tbx];
        unsigned int next                   = tbx + sNext[tgx];
#endif
 
        if (x == y) // Handle diagonals uniquely at 50% efficiency, skipping i == j interactions
        {
#if (__CUDA_ARCH__ >= 300)
            PSFX(tgx)                       = (PMEForce)0;
            PSFY(tgx)                       = (PMEForce)0;
            PSFZ(tgx)                       = (PMEForce)0;
#else
            PMEForce fx_j                   = (PMEForce)0;
            PMEForce fy_j                   = (PMEForce)0;
            PMEForce fz_j                   = (PMEForce)0;
#endif
            PMEFloat xi                     = xyi.x;
            PMEFloat yi                     = xyi.y;
            PSATOMX(tgx)                    = xi;
            PSATOMY(tgx)                    = yi;
            PSATOMZ(tgx)                    = zi;
            ri                             -= cSim.offset;
            PMEFloat ri1i                   = (PMEFloat)1.0 / ri;
            PSATOMR(tgx)                    = ri;
            PSATOMR1I(tgx)                  = ri1i;
            PSATOMS(tgx)                    = si;
            
#if (__CUDA_ARCH__ >= 300)
            shAtom.x                            = __shfl(shAtom.x, shIdx);
            shAtom.y                            = __shfl(shAtom.y, shIdx);
            shAtom.z                            = __shfl(shAtom.z, shIdx);
            shAtom.r                            = __shfl(shAtom.r, shIdx);
            shAtom.r1i                          = __shfl(shAtom.r1i, shIdx);
            shAtom.s                            = __shfl(shAtom.s, shIdx);
#endif

            for (unsigned int j = sNext[tgx]; j != tgx; j = sNext[j])
            {
                PMEFloat xij                = xi - PSATOMX(j);
                PMEFloat yij                = yi - PSATOMY(j);
                PMEFloat zij                = zi - PSATOMZ(j);
                PMEFloat r2                 = xij * xij + yij * yij + zij * zij;
                PMEFloat dij                = sqrt(r2);
                
                PMEFloat sj                 = PSATOMS(j);
                if (dij <= cSim.rgbmax + sj)
                {
                    PMEFloat dij1i          = (PMEFloat)1.0 / dij;
                    PMEFloat dij2i          = dij1i * dij1i;
                    PMEFloat dij3i          = dij2i * dij1i;
                    PMEFloat datmp          = (PMEFloat)0.0;
                    PMEFloat v3             = (PMEFloat)1.0 / (dij + sj);
                    
                    if (dij > cSim.rgbmax - sj)
                    {         
                        PMEFloat temp1      = (PMEFloat)1.0 / (dij - sj);
                        datmp               = (PMEFloat)0.125 * dij3i * ((r2 + sj * sj) *
                                              (temp1 * temp1 - cSim.rgbmax2i) - (PMEFloat)2.0 * log(cSim.rgbmax * temp1));
                    }
                    else if (dij > (PMEFloat)4.0 * sj)
                    {               
                        PMEFloat tmpsd      = sj * sj * dij2i;
                        PMEFloat dumbo      = te + tmpsd * (tf + tmpsd * (tg + tmpsd * (th + tmpsd * thh)));
                        datmp               = tmpsd * sj * dij2i * dij2i * dumbo;
                    }
                    else if (dij > ri + sj)
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (r2 - sj * sj);
                        PMEFloat v4         = log(v3 * (dij - sj));
                        datmp               = v2 * sj * ((PMEFloat)-0.5 * dij2i + v2) + 
                                              (PMEFloat)0.25 * dij3i * v4;
                    }            
                    else if (dij > abs(ri - sj))
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (dij + sj);
                        PMEFloat v4         = log(v3 * ri);
                        datmp               = (PMEFloat)-0.25 * ((PMEFloat)-0.5 * (r2 - ri * ri + sj * sj) *
                                              dij3i * ri1i * ri1i + dij1i * v2 * 
                                              (v2 - dij1i) - dij3i * v4);
                    }    
                    else if (ri < sj)
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (r2 - sj * sj);
                        PMEFloat v4         = log(v3 * (sj - dij));
                        datmp               = (PMEFloat)-0.5 * (sj * dij2i * v2 -
                                               (PMEFloat)2.0 * sj * v2 * v2 - 
                                               (PMEFloat)0.5 * dij3i * v4);
                    }
                    
#ifdef GB_IGB78                    
                    if (dij < ri + PSATOMR(j) + cSim.gb_neckcut)
                    {
                        unsigned int ii     = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        unsigned int jj     = round((PSATOMR(j) - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        PMEFloat2 neckValPos = sNeckMaxValPos[ii * 21 + jj];
                        PMEFloat mdist      = dij - neckValPos.y;
                        PMEFloat mdist2     = mdist * mdist;
                        PMEFloat mdist3     = mdist2 * mdist;
                        PMEFloat mdist5     = mdist2 * mdist3;
                        PMEFloat mdist6     = mdist3 * mdist3;
                        PMEFloat temp1      = (PMEFloat)1.0 + mdist2 + (PMEFloat)0.3 * mdist6;
                        temp1              *= temp1 * dij;
                        datmp              += (((PMEFloat)2.0 * mdist + (PMEFloat)(9.0/5.0) * mdist5) * neckValPos.x * cSim.gb_neckscale) / temp1;
                    }
#endif              
                    
                    datmp                  *= -temp7_i;
#ifdef use_SPFP
                    long long int f_x       = lliroundf(FORCESCALEF * xij * datmp);
#else
                    PMEForce f_x            = (PMEForce)(xij * datmp);
#endif
                    PSFX(j)                += f_x;
                    fx_i                   -= f_x;
#ifdef use_SPFP
                    long long int f_y       = lliroundf(FORCESCALEF * yij * datmp);
#else
                    PMEForce f_y            = (PMEForce)(yij * datmp);
#endif
                    PSFY(j)                += f_y;
                    fy_i                   -= f_y;
#ifdef use_SPFP
                    long long int f_z       = lliroundf(FORCESCALEF * zij * datmp);
#else
                    PMEForce f_z            = (PMEForce)(zij * datmp);
#endif
                    PSFZ(j)                += f_z;
                    fz_i                   -= f_z;        
                }


#if (__CUDA_ARCH__ >= 300)
                shAtom.x                    = __shfl(shAtom.x, shIdx);
                shAtom.y                    = __shfl(shAtom.y, shIdx);
                shAtom.z                    = __shfl(shAtom.z, shIdx);
                shAtom.r                    = __shfl(shAtom.r, shIdx);
                shAtom.r1i                  = __shfl(shAtom.r1i, shIdx);
                shAtom.s                    = __shfl(shAtom.s, shIdx);
#else                
                // Shuffle forces to next thread
                sE[threadIdx.x]             = fx_j;
                fx_j                        = sE[next];    
                sE[threadIdx.x]             = fy_j;
                fy_j                        = sE[next];
                sE[threadIdx.x]             = fz_j;
                fz_j                        = sE[next];                 
#endif
            }
            
            fx_i                           += PSFX(tgx);
            fy_i                           += PSFY(tgx);
            fz_i                           += PSFZ(tgx);
#ifdef use_SPFP
            int offset                      = x + tgx;
            atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset], llitoulli(fx_i));
            atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset], llitoulli(fy_i));
            atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset], llitoulli(fz_i));
#else
            int offset                      = x + tgx + (cSim.nonbondForceBuffers + (x >> GRIDBITS)) * cSim.stride3;
            cSim.pForceXBuffer[offset]      = fx_i;
            cSim.pForceYBuffer[offset]      = fy_i;
            cSim.pForceZBuffer[offset]      = fz_i;
#endif
        }
        else
        {
            unsigned int j                  = y + tgx;
            PMEFloat2 xyj                   = cSim.pAtomXYSP[j];
            PSATOMZ(tgx)                    = cSim.pAtomZSP[j];
            PSATOMR(tgx)                    = cSim.pAtomRBorn[j];
            PSATOMS(tgx)                    = cSim.pAtomS[j];
            PSATOMTEMP7(tgx)                = cSim.pTemp7[j];
#if (__CUDA_ARCH__ >= 300)
            PSFX(tgx)                       = (PMEForce)0;
            PSFY(tgx)                       = (PMEForce)0;
            PSFZ(tgx)                       = (PMEForce)0;
#else
            PMEForce fx_j                   = (PMEForce)0;
            PMEForce fy_j                   = (PMEForce)0;
            PMEForce fz_j                   = (PMEForce)0;
#endif
            PMEFloat xi                     = xyi.x;
            PMEFloat yi                     = xyi.y;
            ri                             -= cSim.offset;
            PMEFloat ri1i                   = (PMEFloat)1.0 / ri;
            PMEFloat si2                    = si * si;     
            PSATOMX(tgx)                    = xyj.x;
            PSATOMY(tgx)                    = xyj.y;
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
                PMEFloat dij3i              = dij2i * dij1i;
                
                // Atom i forces
                PMEFloat sj                 = PSATOMS(j);
                if (dij <= cSim.rgbmax + sj)
                {
                    PMEFloat datmp          = (PMEFloat)0.0;
                    PMEFloat v3             = (PMEFloat)1.0 / (dij + sj);
                    
                    if (dij > cSim.rgbmax - sj)
                    {         
                        PMEFloat temp1      = (PMEFloat)1.0 / (dij - sj);
                        datmp               = (PMEFloat)0.125 * dij3i * ((r2 + sj * sj) *
                                              (temp1 * temp1 - cSim.rgbmax2i) - (PMEFloat)2.0 * log(cSim.rgbmax * temp1));
                    }
                    else if (dij > (PMEFloat)4.0 * sj)
                    {               
                        PMEFloat tmpsd      = sj * sj * dij2i;
                        PMEFloat dumbo      = te + tmpsd * (tf + tmpsd * (tg + tmpsd * (th + tmpsd * thh)));
                        datmp               = tmpsd * sj * dij2i * dij2i * dumbo;
                    }
                    else if (dij > ri + sj)
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (r2 - sj * sj);
                        PMEFloat v4         = log(v3 * (dij - sj));
                        datmp               = v2 * sj * ((PMEFloat)-0.5 * dij2i + v2) + 
                                              (PMEFloat)0.25 * dij3i * v4;
                    }            
                    else if (dij > abs(ri - sj))
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (dij + sj);
                        PMEFloat v4         = log(v3 * ri);
                        datmp               = (PMEFloat)-0.25 * ((PMEFloat)-0.5 * (r2 - ri * ri + sj * sj) *
                                              dij3i * ri1i * ri1i + dij1i * v2 * 
                                              (v2 - dij1i) - dij3i * v4);
                    }    
                    else if (ri < sj)
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (r2 - sj * sj);
                        PMEFloat v4         = log(v3 * (sj - dij));
                        datmp               = (PMEFloat)-0.5 * (sj * dij2i * v2 -
                                               (PMEFloat)2.0 * sj * v2 * v2 - 
                                               (PMEFloat)0.5 * dij3i * v4);
                    }

#ifdef GB_IGB78                    
                    if (dij < ri + PSATOMR(j) + cSim.gb_neckcut)
                    {
                        unsigned int ii     = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        unsigned int jj     = round((PSATOMR(j) - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        PMEFloat2 neckValPos = sNeckMaxValPos[ii * 21 + jj];
                        PMEFloat mdist      = dij - neckValPos.y;
                        PMEFloat mdist2     = mdist * mdist;
                        PMEFloat mdist3     = mdist2 * mdist;
                        PMEFloat mdist5     = mdist2 * mdist3;
                        PMEFloat mdist6     = mdist3 * mdist3;
                        PMEFloat temp1      = (PMEFloat)1.0 + mdist2 + (PMEFloat)0.3 * mdist6;
                        temp1              *= temp1 * dij;
                        datmp              += (((PMEFloat)2.0 * mdist + (PMEFloat)(9.0/5.0) * mdist5) * neckValPos.x * cSim.gb_neckscale) / temp1;
                    }
#endif              
                    
                    datmp                  *= -temp7_i;
#ifdef use_SPFP
                    long long int f_x       = lliroundf(FORCESCALEF * xij * datmp);
#else
                    PMEForce f_x            = (PMEForce)(xij * datmp);
#endif
                    PSFX(j)                += f_x;
                    fx_i                   -= f_x;
#ifdef use_SPFP
                    long long int f_y       = lliroundf(FORCESCALEF * yij * datmp);
#else
                    PMEForce f_y            = (PMEForce)(yij * datmp);
#endif
                    PSFY(j)                += f_y;
                    fy_i                   -= f_y;
#ifdef use_SPFP
                    long long int f_z       = lliroundf(FORCESCALEF * zij * datmp);
#else
                    PMEForce f_z            = (PMEForce)(zij * datmp);
#endif
                    PSFZ(j)                += f_z;            
                    fz_i                   -= f_z;        
                }
                
                // Atom j forces
                if (dij <= cSim.rgbmax + si)
                {
                    PMEFloat datmp          = (PMEFloat)0.0;
                    PMEFloat v3             = (PMEFloat)1.0 / (dij + si);
                    
                    if (dij > cSim.rgbmax - si)
                    {         
                        PMEFloat temp1      = (PMEFloat)1.0 / (dij - si);
                        datmp               = (PMEFloat)0.125 * dij3i * ((r2 + si2) *
                                              (temp1 * temp1 - cSim.rgbmax2i) - (PMEFloat)2.0 * log(cSim.rgbmax * temp1));
                    }
                    else if (dij > (PMEFloat)4.0 * si)
                    {               
                        PMEFloat tmpsd      = si2 * dij2i;
                        PMEFloat dumbo      = te + tmpsd * (tf + tmpsd * (tg + tmpsd * (th + tmpsd * thh)));
                        datmp               = tmpsd * si * dij2i * dij2i * dumbo;
                    }
                    else if (dij > PSATOMR(j) + si)
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (r2 - si2);
                        PMEFloat v4         = log(v3 * (dij - si));
                        datmp               = v2 * si * ((PMEFloat)-0.5 * dij2i + v2) + 
                                              (PMEFloat)0.25 * dij3i * v4;
                    }            
                    else if (dij > abs(PSATOMR(j) - si))
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (dij + si);
                        PMEFloat v4         = log(v3 * PSATOMR(j));
                        datmp               = (PMEFloat)-0.25 * ((PMEFloat)-0.5 * (r2 - PSATOMR(j) * PSATOMR(j) + si2) *
                                              dij3i * PSATOMR1I(j) * PSATOMR1I(j) + dij1i * v2 * 
                                              (v2 - dij1i) - dij3i * v4);
                    }    
                    else if (PSATOMR(j) < si)
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (r2 - si2);
                        PMEFloat v4         = log(v3 * (si - dij));
                        datmp               = (PMEFloat)-0.5 * (si * dij2i * v2 -
                                               (PMEFloat)2.0 * si * v2 * v2 - 
                                               (PMEFloat)0.5 * dij3i * v4);
                    }
 
#ifdef GB_IGB78                    
                    if (dij < ri + PSATOMR(j) + cSim.gb_neckcut)
                    {
                        unsigned int ii     = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        unsigned int jj     = round((PSATOMR(j) - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        PMEFloat2 neckValPos = sNeckMaxValPos[jj * 21 + ii];
                        PMEFloat mdist      = dij - neckValPos.y;
                        PMEFloat mdist2     = mdist * mdist;
                        PMEFloat mdist3     = mdist2 * mdist;
                        PMEFloat mdist5     = mdist2 * mdist3;
                        PMEFloat mdist6     = mdist3 * mdist3;
                        PMEFloat temp1      = (PMEFloat)1.0 + mdist2 + (PMEFloat)0.3 * mdist6;
                        temp1              *= temp1 * dij;
                        datmp              += (((PMEFloat)2.0 * mdist + (PMEFloat)(9.0/5.0) * mdist5) * neckValPos.x * cSim.gb_neckscale) / temp1;
                    }
#endif              
 
                    datmp                  *= PSATOMTEMP7(j);
#ifdef use_SPFP
                    long long int f_x       = lliroundf(FORCESCALEF * xij * datmp);
#else
                    PMEForce f_x           = (PMEForce)(xij * datmp);
#endif
                    fx_i                   += f_x;
                    PSFX(j)                -= f_x;
#ifdef use_SPFP
                    long long int f_y       = lliroundf(FORCESCALEF * yij * datmp);
#else
                    PMEForce f_y           = (PMEForce)(yij * datmp);
#endif
                    fy_i                   += f_y;           
                    PSFY(j)                -= f_y;
#ifdef use_SPFP
                    long long int f_z       = lliroundf(FORCESCALEF * zij * datmp);
#else        
                    PMEForce f_z            = (PMEForce)(zij * datmp);
#endif
                    fz_i                   += f_z;
                    PSFZ(j)                -= f_z;        
                }

#if (__CUDA_ARCH__ >= 300)
                shAtom.x                    = __shfl(shAtom.x, shIdx);
                shAtom.y                    = __shfl(shAtom.y, shIdx);
                shAtom.z                    = __shfl(shAtom.z, shIdx);
                shAtom.r                    = __shfl(shAtom.r, shIdx);
                shAtom.r1i                  = __shfl(shAtom.r1i, shIdx);
                shAtom.s                    = __shfl(shAtom.s, shIdx);
                shAtom.temp7                = __shfl(shAtom.temp7, shIdx);
#else                       
                // Shuffle forces to next thread
                sE[threadIdx.x]             = fx_j;
                fx_j                        = sE[next];    
                sE[threadIdx.x]             = fy_j;
                fy_j                        = sE[next];
                sE[threadIdx.x]             = fz_j;
                fz_j                        = sE[next];
#endif               
                j                           = sNext[j];
            }
            while (j != tgx);
            
            // Write forces
#ifdef use_SPFP
            int offset                      = x + tgx;
            atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset], llitoulli(fx_i));
            atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset], llitoulli(fy_i));
            atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset], llitoulli(fz_i));
            offset                          = y + tgx;
            atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset], llitoulli(PSFX(j)));
            atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset], llitoulli(PSFY(j)));
            atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset], llitoulli(PSFZ(j)));    
#else
            int offset                      = x + tgx + (cSim.nonbondForceBuffers + (y >> GRIDBITS)) * cSim.stride3;
            cSim.pForceXBuffer[offset]      = fx_i;
            cSim.pForceYBuffer[offset]      = fy_i;
            cSim.pForceZBuffer[offset]      = fz_i;
            offset                          = y + tgx + (cSim.nonbondForceBuffers + (x >> GRIDBITS)) * cSim.stride3;
            cSim.pForceXBuffer[offset]      = PSFX(j);
            cSim.pForceYBuffer[offset]      = PSFY(j);
            cSim.pForceZBuffer[offset]      = PSFZ(j);    
#endif
        }           
        if (tgx == 0)
            *psPos                          = atomicAdd(cSim.pGBNB2Position, 1);
    }
}

