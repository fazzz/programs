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
    PMEFloat q;
    PMEFloat sig;
    PMEFloat eps;
    PMEFloat r;
};

#if (__CUDA_ARCH__ >= 200)
volatile __shared__ Atom sA[SM_2X_GBNONBONDENERGY1_THREADS_PER_BLOCK];
volatile __shared__ PMEDouble sE[SM_2X_GBNONBONDENERGY1_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_2X_GBNONBONDENERGY1_THREADS_PER_BLOCK / GRID];
#if defined(GB_ENERGY) && defined(use_SPSP)
volatile __shared__ double sEnergy[SM_2X_GBNONBONDENERGY1_THREADS_PER_BLOCK];
#endif
#else
volatile __shared__ Atom sA[SM_13_GBNONBONDENERGY1_THREADS_PER_BLOCK];
volatile __shared__ PMEDouble sE[SM_13_GBNONBONDENERGY1_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_13_GBNONBONDENERGY1_THREADS_PER_BLOCK / GRID];
#if defined(GB_ENERGY) && defined(use_SPSP)
volatile __shared__ double sEnergy[SM_13_GBNONBONDENERGY1_THREADS_PER_BLOCK];
#endif
#endif
volatile __shared__ unsigned int sNext[GRID];
    if (threadIdx.x < GRID)
    {
        sNext[threadIdx.x]                      = (threadIdx.x + 1) & (GRID - 1);
    }
    __syncthreads();

#ifdef GB_ENERGY
    double egb                                  = 0.0;
    double eelt                                 = 0.0;
    double evdw                                 = 0.0;
#endif

    // Initialize queue position
    volatile unsigned int* psPos                = &sPos[threadIdx.x >> GRIDBITS];
    *psPos                                      = (blockIdx.x * blockDim.x + threadIdx.x) >> GRIDBITS;         

    while (*psPos < cSim.workUnits)
    {  
        // Extract cell coordinates from appropriate work unit
        unsigned int x                          = cSim.pWorkUnit[*psPos];
        unsigned int y                          = ((x >> 2) & 0x7fff) << GRIDBITS;
        x                                       = (x >> 17) << GRIDBITS;
        unsigned int tgx                        = threadIdx.x & (GRID - 1);
        unsigned int i                          = x + tgx;
        
        PMEFloat2 xyi                           = cSim.pAtomXYSP[i];
        PMEFloat2 sigepsi                       = cSim.pAtomSigEps[i];  
        PMEFloat zi                             = cSim.pAtomZSP[i];      
        PMEFloat qi                             = cSim.pAtomChargeSP[i];
        PMEFloat ri                             = cSim.pReffSP[i];
        unsigned int excl                       = 0xffffffff;
        if (*psPos < cSim.excludedWorkUnits)
            excl                                = cSim.pExclusion[*psPos * GRID + tgx];
        PMEDouble fx_i                          = (PMEDouble)0.0;
        PMEDouble fy_i                          = (PMEDouble)0.0;
        PMEDouble fz_i                          = (PMEDouble)0.0;
        PMEDouble sumdeijda_i                   = (PMEDouble)0.0;
        unsigned int tbx                        = threadIdx.x - tgx;
        volatile Atom* psA                      = &sA[tbx];   
        unsigned int next                       = tbx + sNext[tgx];
        if (x == y)
        {
            PMEFloat xi                         = xyi.x;
            PMEFloat yi                         = xyi.y;
            sA[threadIdx.x].x                   = xi;
            sA[threadIdx.x].y                   = yi;
            PMEFloat sigi                       = sigepsi.x;
            PMEFloat epsi                       = sigepsi.y;
            sA[threadIdx.x].z                   = zi;
            sA[threadIdx.x].q                   = qi;
            sA[threadIdx.x].sig                 = sigi;
            sA[threadIdx.x].eps                 = epsi;
            sA[threadIdx.x].r                   = ri;
            
            for (unsigned int j = sNext[tgx]; j != tgx; j = sNext[j])
            {
                PMEFloat xij                    = xi - psA[j].x; 
                PMEFloat yij                    = yi - psA[j].y; 
                PMEFloat zij                    = zi - psA[j].z;
                PMEFloat r2                     = xij * xij + yij * yij + zij * zij;
                PMEFloat qiqj                   = qi * psA[j].q;
               
                PMEFloat rj                     = psA[j].r;
                PMEFloat v1                     = exp(-r2 / ((PMEFloat)4.0 * ri * rj));
                PMEFloat v3                     = r2 + rj * ri * v1;
                PMEFloat v2                     = rsqrt(v3);
                PMEFloat v5                     = rsqrt(r2);
                PMEFloat expmkf                 = cSim.extdiel_inv;
                PMEFloat fgbk                   = (PMEFloat)0.0;
                PMEFloat fgbi                   = v2;
#ifdef GB_ENERGY
                PMEFloat mul                    = fgbi;
#endif
                            
                if (cSim.gb_kappa != (PMEFloat)0.0)
                {
                    v3                          = -cSim.gb_kappa / v2;
                    PMEFloat v4                 = exp(v3);
                    expmkf                     *= v4;
                    fgbk                        = v3 * expmkf;
                    if (cSim.alpb == 1)
                    {
                        fgbk                   += fgbk * cSim.one_arad_beta * (-v3 * cSim.gb_kappa_inv);
#ifdef GB_ENERGY
                        mul                    += cSim.one_arad_beta;
#endif
                    }
                }              

                // vectmp1 = exp(-rij^2/[4*ai*aj])
                // vectmp2 = 1/fij
                // vectmp3 = -kappa*fij - if kappa .ne. 0.0f, otherwise .eq. fij
                // vectmp4 = exp(-kappa*fij)
                // vectmp5 = 1/rij


              
                PMEFloat dl                     = cSim.intdiel_inv - expmkf; 
#ifdef GB_ENERGY
                PMEFloat e                      = -qiqj * dl * mul;
                egb                            += (double)((PMEFloat)0.5 * e);
#endif
                                        
                PMEFloat temp4                  = fgbi * fgbi * fgbi;           // 1.0 / fij^3
                
                // [here, and in the gas-phase part, "de" contains -(1/r)(dE/dr)]
                
                PMEFloat temp6                  = -qiqj * temp4 * (dl + fgbk);
                
                // -qiqj/fij^3*[1/Ein - e(-Kfij)/Eout) -kappa*fij*
                // exp(-kappa*fij)(1 + fij*a*b/A ) /Eout]
                
                PMEFloat temp1                  = v1;
                
                PMEFloat de                     = temp6 * ((PMEFloat)1.0 - (PMEFloat)0.25 * temp1);
                
                PMEFloat temp5                  = (PMEFloat)0.50 * temp1 * temp6 * (ri * rj + (PMEFloat)0.25 * r2);
                
                sumdeijda_i                    += (PMEDouble)(ri * temp5);
              
                PMEFloat rinv                   = v5; // 1.0 / rij
                PMEFloat r2inv                  = rinv * rinv;
                PMEFloat eel                    = cSim.intdiel_inv * qiqj * rinv;
                PMEFloat r6inv                  = r2inv * r2inv * r2inv;
                PMEFloat eps                    = epsi * psA[j].eps;
                PMEFloat sig                    = sigi + psA[j].sig;
                sig                            *= sig * sig;
                sig                            *= sig * r6inv;
                PMEFloat f6                     = eps * sig;
                PMEFloat f12                    = f6 * sig;              
                if (excl & 0x1)
                {
                    de                         += (((PMEFloat)12.0 * f12 - (PMEFloat)6.0 * f6) + eel) * r2inv;
#ifdef GB_ENERGY
                    eelt                       += (double)((PMEFloat)0.5 * eel);    // Necessary to workaround compiler scheduling, ack
                    evdw                       += (double)((PMEFloat)0.5 * (f12 - f6));
#endif
                }
                PMEFloat dedx                   = de * xij;
                PMEFloat dedy                   = de * yij;
                PMEFloat dedz                   = de * zij;
                fx_i                           += (PMEDouble)dedx;
                fy_i                           += (PMEDouble)dedy;
                fz_i                           += (PMEDouble)dedz;
                excl                          >>= 1;
            }
            int offset                          = x + tgx + (x >> GRIDBITS) * cSim.stride3;
            cSim.pForceXBuffer[offset]          = fx_i;
            cSim.pForceYBuffer[offset]          = fy_i;
            cSim.pForceZBuffer[offset]          = fz_i;
            offset                              = x + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pSumdeijdaBuffer[offset]       = sumdeijda_i;
        }       
        else
        {
            int j                               = y + tgx;
            PMEFloat2 xyj                       = cSim.pAtomXYSP[j];
            PMEFloat2 sigepsj                   = cSim.pAtomSigEps[j];
            sA[threadIdx.x].z                   = cSim.pAtomZSP[j];
            sA[threadIdx.x].q                   = cSim.pAtomChargeSP[j];
            sA[threadIdx.x].r                   = cSim.pReffSP[j];            
            PMEDouble fx_j                      = (PMEDouble)0.0;
            PMEDouble fy_j                      = (PMEDouble)0.0;
            PMEDouble fz_j                      = (PMEDouble)0.0;
            PMEDouble sumdeijda_j               = (PMEDouble)0.0;   
            PMEFloat xi                         = xyi.x;
            PMEFloat yi                         = xyi.y;
            PMEFloat sigi                       = sigepsi.x;
            PMEFloat epsi                       = sigepsi.y;
            sA[threadIdx.x].x                   = xyj.x;
            sA[threadIdx.x].y                   = xyj.y;
            sA[threadIdx.x].sig                 = sigepsj.x;
            sA[threadIdx.x].eps                 = sigepsj.y;
            j                                   = tgx;
#if defined(GB_ENERGY) && defined(MPI)
            double oldeelt                      = eelt;
            double oldevdw                      = evdw;
            double oldegb                       = egb;
#endif            
            do
            {
                PMEFloat xij                    = xi - psA[j].x; 
                PMEFloat yij                    = yi - psA[j].y; 
                PMEFloat zij                    = zi - psA[j].z;
                PMEFloat r2                     = xij * xij + yij * yij + zij * zij;
                PMEFloat qiqj                   = qi * psA[j].q;
               
                PMEFloat rj                     = psA[j].r;
                PMEFloat v1                     = exp(-r2 / ((PMEFloat)4.0 * ri * rj));
                PMEFloat v3                     = r2 + rj * ri * v1;
                PMEFloat v2                     = rsqrt(v3);
                PMEFloat v5                     = rsqrt(r2);
                PMEFloat expmkf                 = cSim.extdiel_inv;
                PMEFloat fgbk                   = (PMEFloat)0.0;
                PMEFloat fgbi                   = v2;
#ifdef GB_ENERGY
                PMEFloat mul                    = fgbi;
#endif                              
                if (cSim.gb_kappa != (PMEFloat)0.0)
                {
                    v3                          = -cSim.gb_kappa / v2;
                    PMEFloat v4                 = exp(v3);
                    expmkf                     *= v4;
                    fgbk                        = v3 * expmkf;
                    if (cSim.alpb == 1)
                    {
                        fgbk                   += fgbk * cSim.one_arad_beta * (-v3 * cSim.gb_kappa_inv);
#ifdef GB_ENERGY
                        mul                    += cSim.one_arad_beta;
#endif
                    }
                }       
                // vectmp1 = exp(-rij^2/[4*ai*aj])
                // vectmp2 = 1/fij
                // vectmp3 = -kappa*fij - if kappa .ne. 0.0f, otherwise .eq. fij
                // vectmp4 = exp(-kappa*fij)
                // vectmp5 = 1/rij


              
                PMEFloat dl                     = cSim.intdiel_inv - expmkf;
#ifdef GB_ENERGY         
                PMEFloat e                      = -qiqj * dl * mul;
                egb                            += (double)e;
#endif
                PMEFloat temp4                  = fgbi * fgbi * fgbi;           // 1.0 / fij^3
                
                // [here, and in the gas-phase part, "de" contains -(1/r)(dE/dr)]
                
                PMEFloat temp6                  = -qiqj * temp4 * (dl + fgbk);
                
                // -qiqj/fij^3*[1/Ein - e(-Kfij)/Eout) -kappa*fij*
                // exp(-kappa*fij)(1 + fij*a*b/A ) /Eout]
                
                PMEFloat temp1                  = v1;
                
                PMEFloat de                     = temp6 * ((PMEFloat)1.0 - (PMEFloat)0.25 * temp1);
                
                PMEFloat temp5                  = (PMEFloat)0.50 * temp1 * temp6 * (ri * rj + (PMEFloat)0.25 * r2);
                
                sumdeijda_i                    += (PMEDouble)(ri * temp5);
                sumdeijda_j                    += (PMEDouble)(rj * temp5);

                PMEFloat rinv                   = v5; // 1.0 / rij
                PMEFloat r2inv                  = rinv * rinv;
                PMEFloat eel                    = cSim.intdiel_inv * qiqj * rinv;
                PMEFloat r6inv                  = r2inv * r2inv * r2inv;
                PMEFloat eps                    = epsi * psA[j].eps;
                PMEFloat sig                    = sigi + psA[j].sig;
                sig                            *= sig * sig;
                sig                            *= sig * r6inv;
                PMEFloat f6                     = eps * sig;
                PMEFloat f12                    = f6 * sig;
                if (excl & 0x1)
                {
#ifdef GB_ENERGY    
                    eelt                       += (double)eel;
                    evdw                       += (double)(f12 - f6);
#endif
                    de                         += (((PMEFloat)12.0 * f12 - (PMEFloat)6.0 * f6) + eel) * r2inv;
                }
                PMEDouble dedx                  = (double)(de * xij);
                PMEDouble dedy                  = (double)(de * yij);
                PMEDouble dedz                  = (double)(de * zij);
                fx_i                           += dedx;
                fy_i                           += dedy;
                fz_i                           += dedz;
                fx_j                           -= dedx;
                fy_j                           -= dedy;
                fz_j                           -= dedz;

                excl                          >>= 1;
                       
                // Shuffle forces to next thread
                sE[threadIdx.x]                 = fx_j;
                fx_j                            = sE[next];    
                sE[threadIdx.x]                 = fy_j;
                fy_j                            = sE[next];
                sE[threadIdx.x]                 = fz_j;
                fz_j                            = sE[next];
                sE[threadIdx.x]                 = sumdeijda_j;
                sumdeijda_j                     = sE[next];    
                j                               = sNext[j];         
            }
            while (j != tgx);
#if defined(GB_ENERGY) && defined(MPI)
            // Insanely stupid workaround for wacky compiler re-ordering if one tries to not calculate
            // these terms for tiles outside the GPU's local region (because they're not needed).  Stupid
            // but necessary solution: throw them out.  Heaven forbid the damned thing just kept the same ordering.
            if ((x < cSim.minLocalAtom) || (x >= cSim.maxLocalAtom))
            {
                eelt                            = oldeelt;
                evdw                            = oldevdw;
                egb                             = oldegb;
            }
#endif            
            
            
            // Write forces
            int offset                          = x + tgx + (y >> GRIDBITS) * cSim.stride3;
            cSim.pForceXBuffer[offset]          = fx_i;
            cSim.pForceYBuffer[offset]          = fy_i;
            cSim.pForceZBuffer[offset]          = fz_i;
            offset                              = x + tgx + (y >> GRIDBITS) * cSim.stride;
            cSim.pSumdeijdaBuffer[offset]       = sumdeijda_i;
            offset                              = y + tgx + (x >> GRIDBITS) * cSim.stride3;
            cSim.pForceXBuffer[offset]          = fx_j;
            cSim.pForceYBuffer[offset]          = fy_j;
            cSim.pForceZBuffer[offset]          = fz_j;
            offset                              = y + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pSumdeijdaBuffer[offset]       = sumdeijda_j;
        }
        if (tgx == 0)
            *psPos                              = atomicAdd(cSim.pGBNB1Position, 1);
    }

#ifdef GB_ENERGY
#ifdef use_SPSP
#define sENERGY(x) sEnergy[x]
#else
#define sENERGY(x) sE[x]
#endif
    // Reduce and write energies
    sENERGY(threadIdx.x)                        = egb;
    __syncthreads();
    unsigned int m                              = 1;
    while (m < blockDim.x)
    {
        int p                                   = threadIdx.x + m;
        double d                                = ((p < blockDim.x) ? sENERGY(p) : (double)0.0);
        __syncthreads();
        sENERGY(threadIdx.x)                   += d;
        __syncthreads();
        m                                      *= 2;
    }
    egb                                         = sENERGY(threadIdx.x);
    unsigned long long int val1                 = (unsigned long long int)(fabs(egb) * ENERGYSCALE + (double)0.5);
    if (egb < (double)0.0)
        val1                                    = 0ull - val1;
   
    sENERGY(threadIdx.x)                        = eelt;
    __syncthreads();
    m                                           = 1;
    while (m < blockDim.x)
    {
        int p                                   = threadIdx.x + m;
        double d                                = ((p < blockDim.x) ? sENERGY(p) : (double)0.0);
        __syncthreads();
        sENERGY(threadIdx.x)                   += d;
        __syncthreads();
        m                                      *= 2;
    }
    eelt                                        = sENERGY(threadIdx.x);
    unsigned long long int val2                 = (unsigned long long int)(fabs(eelt) * ENERGYSCALE + (double)0.5);
    if (eelt < (double)0.0)
        val2                                    = 0ull - val2;
        
    sENERGY(threadIdx.x)                        = evdw;
    __syncthreads();
    m                                           = 1;
    while (m < blockDim.x)
    {
        int p                                   = threadIdx.x + m;
        double d                                = ((p < blockDim.x) ? sENERGY(p) : (double)0.0);
        __syncthreads();
        sENERGY(threadIdx.x)                   += d;
        __syncthreads();
        m                                      *= 2;
    }
    evdw                                        = sENERGY(threadIdx.x);
    unsigned long long int val3                 = (unsigned long long int)(fabs(evdw) * ENERGYSCALE + (double)0.5);
    if (evdw < (double)0.0)
        val3                                    = 0ull - val3;
    
    // Write out energies
    if (threadIdx.x == 0)
    {
       atomicAdd(cSim.pEGB, val1);              
       atomicAdd(cSim.pEELT, val2);              
       atomicAdd(cSim.pEVDW, val3);              
    }
#ifdef use_SPSP
#undef sENERGY
#endif
#endif
}

