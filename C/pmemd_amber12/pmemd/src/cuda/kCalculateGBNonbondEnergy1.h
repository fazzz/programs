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

struct NLForce
{
    PMEForce x;
    PMEForce y;
    PMEForce z;
};

#if (__CUDA_ARCH__ >= 300)
#define PSATOMX(i)      shAtom.x
#define PSATOMY(i)      shAtom.y
#define PSATOMZ(i)      shAtom.z
#define PSATOMQ(i)      shAtom.q
#define PSATOMSIG(i)    shAtom.sig
#define PSATOMEPS(i)    shAtom.eps
#define PSATOMR(i)      shAtom.r
#define PSFX(i)         psF[i].x
#define PSFY(i)         psF[i].y
#define PSFZ(i)         psF[i].z
#define PSSUMDEIJDA(i)  psS[i]
#else
#define PSATOMX(i)      psA[i].x
#define PSATOMY(i)      psA[i].y
#define PSATOMZ(i)      psA[i].z
#define PSATOMQ(i)      psA[i].q
#define PSATOMSIG(i)    psA[i].sig
#define PSATOMEPS(i)    psA[i].eps
#define PSATOMR(i)      psA[i].r
#define PSFX(i)         fx_j
#define PSFY(i)         fy_j
#define PSFZ(i)         fz_j
#define PSSUMDEIJDA(i)  sumdeijda_j
#endif


#if (__CUDA_ARCH__ >= 300)
Atom shAtom;
volatile __shared__ NLForce sForce[SM_3X_GBNONBONDENERGY1_THREADS_PER_BLOCK];
volatile __shared__ PMEForce sSumdeijda[SM_3X_GBNONBONDENERGY1_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_3X_GBNONBONDENERGY1_THREADS_PER_BLOCK / GRID];
#elif (__CUDA_ARCH__ >= 200)
volatile __shared__ Atom sA[SM_2X_GBNONBONDENERGY1_THREADS_PER_BLOCK];
volatile __shared__ PMEForce sE[SM_2X_GBNONBONDENERGY1_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_2X_GBNONBONDENERGY1_THREADS_PER_BLOCK / GRID];
#else
volatile __shared__ Atom sA[SM_13_GBNONBONDENERGY1_THREADS_PER_BLOCK];
volatile __shared__ PMEForce sE[SM_13_GBNONBONDENERGY1_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_13_GBNONBONDENERGY1_THREADS_PER_BLOCK / GRID];
#endif
volatile __shared__ unsigned int sNext[GRID];
    if (threadIdx.x < GRID)
    {
        sNext[threadIdx.x]                      = (threadIdx.x + 1) & (GRID - 1);
    }
    __syncthreads();

#ifdef GB_ENERGY
    PMEForce egb                                = (PMEForce)0;
    PMEForce eelt                               = (PMEForce)0;
    PMEForce evdw                               = (PMEForce)0;
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
        PMEForce fx_i                           = (PMEForce)0;
        PMEForce fy_i                           = (PMEForce)0;
        PMEForce fz_i                           = (PMEForce)0;
        PMEForce sumdeijda_i                    = (PMEForce)0;
        unsigned int tbx                        = threadIdx.x - tgx;
#if (__CUDA_ARCH__ >= 300)
        volatile NLForce* psF                   = &sForce[tbx];
        volatile PMEForce* psS                  = &sSumdeijda[tbx];
        unsigned int shIdx                      = sNext[tgx]; 
#else
        volatile Atom* psA                      = &sA[tbx];   
        unsigned int next                       = tbx + sNext[tgx];
#endif
        if (x == y)
        {
            PMEFloat xi                         = xyi.x;
            PMEFloat yi                         = xyi.y;
            PSATOMX(tgx)                        = xi;
            PSATOMY(tgx)                        = yi;
            PMEFloat sigi                       = sigepsi.x;
            PMEFloat epsi                       = sigepsi.y;
            PSATOMZ(tgx)                        = zi;
            PSATOMQ(tgx)                        = qi;
            PSATOMSIG(tgx)                      = sigi;
            PSATOMEPS(tgx)                      = epsi;
            PSATOMR(tgx)                        = ri;

#if (__CUDA_ARCH__ >= 300)
            shAtom.x                            = __shfl(shAtom.x, shIdx);
            shAtom.y                            = __shfl(shAtom.y, shIdx);
            shAtom.z                            = __shfl(shAtom.z, shIdx);
            shAtom.q                            = __shfl(shAtom.q, shIdx);
            shAtom.r                            = __shfl(shAtom.r, shIdx);
            shAtom.sig                          = __shfl(shAtom.sig, shIdx);
            shAtom.eps                          = __shfl(shAtom.eps, shIdx);
#endif
            
            for (unsigned int j = sNext[tgx]; j != tgx; j = sNext[j])
            {
                PMEFloat xij                    = xi - PSATOMX(j);
                PMEFloat yij                    = yi - PSATOMY(j);
                PMEFloat zij                    = zi - PSATOMZ(j);
                PMEFloat r2                     = xij * xij + yij * yij + zij * zij;
                PMEFloat qiqj                   = qi * PSATOMQ(j);
               
                PMEFloat rj                     = PSATOMR(j);
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
#ifdef use_SPFP
                egb                            += lliroundf((PMEFloat)0.5 * ENERGYSCALEF * e);
#else
                egb                            += (PMEForce)((PMEFloat)0.5 * e);
#endif
#endif
                                        
                PMEFloat temp4                  = fgbi * fgbi * fgbi;           // 1.0 / fij^3
                
                // [here, and in the gas-phase part, "de" contains -(1/r)(dE/dr)]
                
                PMEFloat temp6                  = -qiqj * temp4 * (dl + fgbk);
                
                // -qiqj/fij^3*[1/Ein - e(-Kfij)/Eout) -kappa*fij*
                // exp(-kappa*fij)(1 + fij*a*b/A ) /Eout]
                
                PMEFloat temp1                  = v1;
                
                
#ifdef use_SPFP
                PMEFloat de                     = temp6 * (FORCESCALEF - (PMEFloat)0.25 * FORCESCALEF * temp1);
                PMEFloat temp5                  = (PMEFloat)0.50 * FORCESCALEF * temp1 * temp6 * (ri * rj + (PMEFloat)0.25 * r2);
                sumdeijda_i                    += lliroundf(ri * temp5);
#else
                PMEFloat de                     = temp6 * ((PMEFloat)1.0 - (PMEFloat)0.25 * temp1);
                PMEFloat temp5                  = (PMEFloat)0.50 * temp1 * temp6 * (ri * rj + (PMEFloat)0.25 * r2);
                sumdeijda_i                    += (PMEDouble)(ri * temp5);
#endif
              
                PMEFloat rinv                   = v5; // 1.0 / rij
                PMEFloat r2inv                  = rinv * rinv;
                PMEFloat eel                    = cSim.intdiel_inv * qiqj * rinv;
                PMEFloat r6inv                  = r2inv * r2inv * r2inv;
                PMEFloat eps                    = epsi * PSATOMEPS(j);
                PMEFloat sig                    = sigi + PSATOMSIG(j);
                sig                            *= sig * sig;
                sig                            *= sig * r6inv;
                PMEFloat f6                     = eps * sig;
                PMEFloat f12                    = f6 * sig;              
                if (excl & 0x1)
                {
#ifdef use_SPFP
                    de                         += FORCESCALEF * (((PMEFloat)12.0 * f12 - (PMEFloat)6.0 * f6) + eel) * r2inv;
#ifdef GB_ENERGY
                    eelt                       += lliroundf((PMEFloat)0.5 * ENERGYSCALEF * eel);    // Necessary to workaround compiler scheduling, ack
                    evdw                       += lliroundf((PMEFloat)0.5 * ENERGYSCALEF * (f12 - f6));
#endif
                }
                PMEFloat dedx                   = de * xij;
                PMEFloat dedy                   = de * yij;
                PMEFloat dedz                   = de * zij;
                fx_i                           += lliroundf(dedx);
                fy_i                           += lliroundf(dedy);
                fz_i                           += lliroundf(dedz);
#else
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
#endif
                excl                          >>= 1;

#if (__CUDA_ARCH__ >= 300)
                shAtom.x                        = __shfl(shAtom.x, shIdx);
                shAtom.y                        = __shfl(shAtom.y, shIdx);
                shAtom.z                        = __shfl(shAtom.z, shIdx);
                shAtom.q                        = __shfl(shAtom.q, shIdx);
                shAtom.r                        = __shfl(shAtom.r, shIdx);
                shAtom.sig                      = __shfl(shAtom.sig, shIdx);
                shAtom.eps                      = __shfl(shAtom.eps, shIdx);
#endif
            }

#ifdef use_SPFP
            int offset                          = x + tgx;
            atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset], llitoulli(fx_i));
            atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset], llitoulli(fy_i));
            atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset], llitoulli(fz_i));
            atomicAdd((unsigned long long int*)&cSim.pSumdeijdaAccumulator[offset], llitoulli(sumdeijda_i));
#else
            int offset                          = x + tgx + (x >> GRIDBITS) * cSim.stride3;
            cSim.pForceXBuffer[offset]          = fx_i;
            cSim.pForceYBuffer[offset]          = fy_i;
            cSim.pForceZBuffer[offset]          = fz_i;
            offset                              = x + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pSumdeijdaBuffer[offset]       = sumdeijda_i;
#endif
        }       
        else
        {
            int j                               = y + tgx;
            PMEFloat2 xyj                       = cSim.pAtomXYSP[j];
            PMEFloat2 sigepsj                   = cSim.pAtomSigEps[j];
            PSATOMZ(tgx)                        = cSim.pAtomZSP[j];
            PSATOMQ(tgx)                        = cSim.pAtomChargeSP[j];
            PSATOMR(tgx)                        = cSim.pReffSP[j];
#ifdef use_SPFP
#if (__CUDA_ARCH__ >= 300)
            PSFX(tgx)                           = 0;
            PSFY(tgx)                           = 0;
            PSFZ(tgx)                           = 0;
            PSSUMDEIJDA(tgx)                    = 0;   
#else           
            long long int fx_j                  = 0;
            long long int fy_j                  = 0;
            long long int fz_j                  = 0;
            long long int sumdeijda_j           = 0;   
#endif
#else
#if (__CUDA_ARCH__ >= 300)
            PSFX(tgx)                           = (PMEForce)0;
            PSFY(tgx)                           = (PMEForce)0;
            PSFZ(tgx)                           = (PMEForce)0;
            PSSUMDEIJDA(tgx)                    = (PMEForce)0;   
#else           
            PMEForce fx_j                       = (PMEForce)0;
            PMEForce fy_j                       = (PMEForce)0;
            PMEForce fz_j                       = (PMEForce)0;
            PMEForce sumdeijda_j                = (PMEForce)0;   
#endif
#endif
            PMEFloat xi                         = xyi.x;
            PMEFloat yi                         = xyi.y;
            PMEFloat sigi                       = sigepsi.x;
            PMEFloat epsi                       = sigepsi.y;
            PSATOMX(tgx)                        = xyj.x;
            PSATOMY(tgx)                        = xyj.y;
            PSATOMSIG(tgx)                      = sigepsj.x;
            PSATOMEPS(tgx)                      = sigepsj.y;
            j                                   = tgx;
#if defined(GB_ENERGY) && defined(MPI)
            PMEForce oldeelt                    = eelt;
            PMEForce oldevdw                    = evdw;
            PMEForce oldegb                     = egb;
#endif            
            do
            {
                PMEFloat xij                    = xi - PSATOMX(j);
                PMEFloat yij                    = yi - PSATOMY(j);
                PMEFloat zij                    = zi - PSATOMZ(j);
                PMEFloat r2                     = xij * xij + yij * yij + zij * zij;
                PMEFloat qiqj                   = qi * PSATOMQ(j);
               
                PMEFloat rj                     = PSATOMR(j);
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
#ifdef use_SPFP
                egb                            += lliroundf(ENERGYSCALEF * e);
#else
                egb                            += (PMEForce)e;
#endif
#endif
                PMEFloat temp4                  = fgbi * fgbi * fgbi;           // 1.0 / fij^3
                
                // [here, and in the gas-phase part, "de" contains -(1/r)(dE/dr)]
                
                PMEFloat temp6                  = -qiqj * temp4 * (dl + fgbk);
                
                // -qiqj/fij^3*[1/Ein - e(-Kfij)/Eout) -kappa*fij*
                // exp(-kappa*fij)(1 + fij*a*b/A ) /Eout]
                
                PMEFloat temp1                  = v1;

#ifdef use_SPFP
                PMEFloat de                     = temp6 * (FORCESCALEF - (PMEFloat)0.25 * FORCESCALEF * temp1);                
                PMEFloat temp5                  = (PMEFloat)0.50 * FORCESCALEF * temp1 * temp6 * (ri * rj + (PMEFloat)0.25 * r2);           
                sumdeijda_i                    += lliroundf(ri * temp5);
                PSSUMDEIJDA(j)                 += lliroundf(rj * temp5);
#else
                PMEFloat de                     = temp6 * ((PMEFloat)1.0 - (PMEFloat)0.25 * temp1);                
                PMEFloat temp5                  = (PMEFloat)0.50 * temp1 * temp6 * (ri * rj + (PMEFloat)0.25 * r2);           
                sumdeijda_i                    += (PMEForce)(ri * temp5);
                PSSUMDEIJDA(j)                 += (PMEForce)(rj * temp5);
#endif

                PMEFloat rinv                   = v5; // 1.0 / rij
                PMEFloat r2inv                  = rinv * rinv;
                PMEFloat eel                    = cSim.intdiel_inv * qiqj * rinv;
                PMEFloat r6inv                  = r2inv * r2inv * r2inv;
                PMEFloat eps                    = epsi * PSATOMEPS(j);
                PMEFloat sig                    = sigi + PSATOMSIG(j);
                sig                            *= sig * sig;
                sig                            *= sig * r6inv;
                PMEFloat f6                     = eps * sig;
                PMEFloat f12                    = f6 * sig;
                if (excl & 0x1)
                {

#ifdef use_SPFP    
#ifdef GB_ENERGY
                    eelt                       += lliroundf(ENERGYSCALEF * eel);
                    evdw                       += lliroundf(ENERGYSCALEF * (f12 - f6));
#endif
                    de                         += FORCESCALEF * (((PMEFloat)12.0 * f12 - (PMEFloat)6.0 * f6) + eel) * r2inv;
                }
                long long int dedx              = lliroundf(de * xij);
                long long int dedy              = lliroundf(de * yij);
                long long int dedz              = lliroundf(de * zij);
#else
#ifdef GB_ENERGY
                    eelt                       += (PMEForce)eel;
                    evdw                       += (PMEForce)(f12 - f6);
#endif
                    de                         += (((PMEFloat)12.0 * f12 - (PMEFloat)6.0 * f6) + eel) * r2inv;
                }
                PMEForce dedx                   = (PMEForce)(de * xij);
                PMEForce dedy                   = (PMEForce)(de * yij);
                PMEForce dedz                   = (PMEForce)(de * zij);
#endif
                fx_i                           += dedx;
                fy_i                           += dedy;
                fz_i                           += dedz;
                PSFX(j)                        -= dedx;
                PSFY(j)                        -= dedy;
                PSFZ(j)                        -= dedz;

                excl                          >>= 1;

#if (__CUDA_ARCH__ >= 300)
                shAtom.x                        = __shfl(shAtom.x, shIdx);
                shAtom.y                        = __shfl(shAtom.y, shIdx);
                shAtom.z                        = __shfl(shAtom.z, shIdx);
                shAtom.q                        = __shfl(shAtom.q, shIdx);
                shAtom.r                        = __shfl(shAtom.r, shIdx);
                shAtom.sig                      = __shfl(shAtom.sig, shIdx);
                shAtom.eps                      = __shfl(shAtom.eps, shIdx);
#else                       
                // Shuffle forces to next thread
                sE[threadIdx.x]                 = fx_j;
                fx_j                            = sE[next];    
                sE[threadIdx.x]                 = fy_j;
                fy_j                            = sE[next];
                sE[threadIdx.x]                 = fz_j;
                fz_j                            = sE[next];
                sE[threadIdx.x]                 = sumdeijda_j;
                sumdeijda_j                     = sE[next];
#endif
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

#ifdef use_SPFP
            int offset                          = x + tgx;
            atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset], llitoulli(fx_i));
            atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset], llitoulli(fy_i));
            atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset], llitoulli(fz_i));
            atomicAdd((unsigned long long int*)&cSim.pSumdeijdaAccumulator[offset], llitoulli(sumdeijda_i));
            offset                              = y + tgx;
            atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset], llitoulli(PSFX(tgx)));
            atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset], llitoulli(PSFY(tgx)));
            atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset], llitoulli(PSFZ(tgx)));
            atomicAdd((unsigned long long int*)&cSim.pSumdeijdaAccumulator[offset], llitoulli(PSSUMDEIJDA(tgx)));
#else
            int offset                          = x + tgx + (y >> GRIDBITS) * cSim.stride3;
            cSim.pForceXBuffer[offset]          = fx_i;
            cSim.pForceYBuffer[offset]          = fy_i;
            cSim.pForceZBuffer[offset]          = fz_i;
            offset                              = x + tgx + (y >> GRIDBITS) * cSim.stride;
            cSim.pSumdeijdaBuffer[offset]       = sumdeijda_i;
            offset                              = y + tgx + (x >> GRIDBITS) * cSim.stride3;
            cSim.pForceXBuffer[offset]          = PSFX(tgx);
            cSim.pForceYBuffer[offset]          = PSFY(tgx);
            cSim.pForceZBuffer[offset]          = PSFZ(tgx);
            offset                              = y + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pSumdeijdaBuffer[offset]       = PSSUMDEIJDA(tgx);
#endif
        }
        if (tgx == 0)
            *psPos                              = atomicAdd(cSim.pGBNB1Position, 1);
    }

#ifdef GB_ENERGY
#if (__CUDA_ARCH__ >= 300)
#define sENERGY(x) sSumdeijda[x]
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
        PMEForce d                              = ((p < blockDim.x) ? sENERGY(p) : (PMEForce)0);
        __syncthreads();
        sENERGY(threadIdx.x)                   += d;
        __syncthreads();
        m                                      *= 2;
    }
    egb                                         = sENERGY(threadIdx.x);
#ifdef use_SPFP
    unsigned long long int val1                 = llitoulli(egb);
#else
    unsigned long long int val1                 = llitoulli(lliroundd(egb * ENERGYSCALE));
#endif
   
    sENERGY(threadIdx.x)                        = eelt;
    __syncthreads();
    m                                           = 1;
    while (m < blockDim.x)
    {
        int p                                   = threadIdx.x + m;
        PMEForce d                              = ((p < blockDim.x) ? sENERGY(p) : (PMEForce)0);
        __syncthreads();
        sENERGY(threadIdx.x)                   += d;
        __syncthreads();
        m                                      *= 2;
    }
    eelt                                        = sENERGY(threadIdx.x);
#ifdef use_SPFP
    unsigned long long int val2                 = llitoulli(eelt);
#else
    unsigned long long int val2                 = llitoulli(lliroundd(eelt * ENERGYSCALE));
#endif
        
    sENERGY(threadIdx.x)                        = evdw;
    __syncthreads();
    m                                           = 1;
    while (m < blockDim.x)
    {
        int p                                   = threadIdx.x + m;
        PMEForce d                              = ((p < blockDim.x) ? sENERGY(p) : (PMEForce)0);
        __syncthreads();
        sENERGY(threadIdx.x)                   += d;
        __syncthreads();
        m                                      *= 2;
    }
    evdw                                        = sENERGY(threadIdx.x);
#ifdef use_SPFP
    unsigned long long int val3                 = llitoulli(evdw);
#else
    unsigned long long int val3                 = llitoulli(lliroundd(evdw * ENERGYSCALE));
#endif
    
    // Write out energies
    if (threadIdx.x == 0)
    {
       atomicAdd(cSim.pEGB, val1);              
       atomicAdd(cSim.pEELT, val2);              
       atomicAdd(cSim.pEVDW, val3);              
    }
#undef sENERGY
#endif
}

