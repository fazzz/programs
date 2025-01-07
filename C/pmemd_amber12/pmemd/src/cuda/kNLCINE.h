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
// #defines: IPS_ENERGY, IPS_VIRIAL, IPS_IS_ORTHOGONAL, IPS_ATOMS_PER_WARP

struct NLAtom
{
    PMEFloat x;
    PMEFloat y;
    PMEFloat z;
    PMEFloat q;
    PMEFloat sig;
    PMEFloat eps;
    unsigned int ID;
};

struct NLForce
{
    PMEForce x;
    PMEForce y;
    PMEForce z;
};

struct NLVirial
{
    long long int vir_11;
    long long int vir_22;
    long long int vir_33;
};

struct NLWarp
{
    NLEntry nlEntry;
    uint pos;
    uint xmax;
    uint ypos;
    uint ymax;
    uint yend;
    uint nlpos;
    uint offset;
    uint bufferOffset;
    uint homeCellBuffer;
};

#if (__CUDA_ARCH__ >= 300)
#define PSATOMX(i) shAtom.x
#define PSATOMY(i) shAtom.y
#define PSATOMZ(i) shAtom.z
#define PSATOMQ(i) shAtom.q
#define PSATOMSIG(i) shAtom.sig
#define PSATOMEPS(i) shAtom.eps
#define PSATOMID(i) shAtom.ID
#else
#define PSATOMX(i) psA[i].x
#define PSATOMY(i) psA[i].y
#define PSATOMZ(i) psA[i].z
#define PSATOMQ(i) psA[i].q
#define PSATOMSIG(i) psA[i].sig
#define PSATOMEPS(i) psA[i].eps
#define PSATOMID(i) psA[i].ID
#endif

#if defined(IPS_VIRIAL)
__shared__ PMEFloat sUcellf[9];
#endif
__shared__ unsigned int sNext[GRID];
#if (IPS_ATOMS_PER_WARP == 16)
__shared__ unsigned int sStart[GRID];
__shared__ unsigned int sEnd[GRID];
#endif
#if (__CUDA_ARCH__ >= 300)
__shared__ volatile NLWarp sWarp[SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK / GRID];
NLAtom shAtom;
__shared__ volatile NLForce sForce[SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK];
#ifdef IPS_VIRIAL
__shared__ volatile NLVirial sWarpVirial[SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK / GRID];
#endif
#elif (__CUDA_ARCH__ >= 200)
__shared__ volatile NLWarp sWarp[SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK / GRID];
__shared__ volatile NLAtom sAtom[SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK];
__shared__ volatile NLForce sForce[SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK];
#ifdef IPS_VIRIAL
__shared__ volatile NLVirial sWarpVirial[SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK / GRID];
#endif
#else
__shared__ volatile NLWarp sWarp[SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK / GRID];
__shared__ volatile NLAtom sAtom[SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK];
__shared__ volatile NLForce sForce[SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK];
#ifdef IPS_VIRIAL
__shared__ volatile NLVirial sWarpVirial[SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK / GRID];
#endif
#endif

    // Read static data
    if (threadIdx.x < GRID)
    {
        unsigned int offset                             = cSim.NLAtomsPerWarp * (threadIdx.x >> cSim.NLAtomsPerWarpBits);
        sNext[threadIdx.x]                              = ((threadIdx.x + 1) & cSim.NLAtomsPerWarpBitsMask) + offset;

#if (IPS_ATOMS_PER_WARP == 16)
        sStart[threadIdx.x]                             = (threadIdx.x + 1 + 8 * (threadIdx.x >> cSim.NLAtomsPerWarpBits)) & cSim.NLAtomsPerWarpBitsMask; // SIZE DEPENDENT 2, 8
        if (threadIdx.x < GRID - cSim.NLAtomsPerWarp)           // SIZE DEPENDENT (3, 2, 2) (9, 8, 8)
        {
            sEnd[threadIdx.x]                           = (((threadIdx.x + 9) & cSim.NLAtomsPerWarpBitsMask) +  
                                                          8 * (threadIdx.x >> cSim.NLAtomsPerWarpBits)) & cSim.NLAtomsPerWarpBitsMask;
        }
        else
        {
            sEnd[threadIdx.x]                           = (((threadIdx.x + 8) & cSim.NLAtomsPerWarpBitsMask) +  
                                                          8 * (threadIdx.x >> cSim.NLAtomsPerWarpBits)) & cSim.NLAtomsPerWarpBitsMask;
        }
#endif    
    }
     
#ifdef IPS_VIRIAL
    if (threadIdx.x < 9)
        sUcellf[threadIdx.x]                            = cSim.pNTPData->ucellf[threadIdx.x];
#if (__CUDA_ARCH__ >= 300)        
    if (threadIdx.x < SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK / GRID)
#elif (__CUDA_ARCH__ >= 200)        
    if (threadIdx.x < SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK / GRID)
#else
    if (threadIdx.x < SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK / GRID)
#endif
    {
        sWarpVirial[threadIdx.x].vir_11                 = 0;
        sWarpVirial[threadIdx.x].vir_22                 = 0;
        sWarpVirial[threadIdx.x].vir_33                 = 0;                
    }
#endif    

#ifdef IPS_ENERGY    
    PMEForce eed                                        = (PMEForce)0;
    PMEForce evdw                                       = (PMEForce)0;
#endif    
    volatile NLWarp* psWarp                             = &sWarp[threadIdx.x / GRID];
    psWarp->pos                                         = (blockIdx.x * blockDim.x + threadIdx.x) / GRID;
    __syncthreads();
 
    while (psWarp->pos < cSim.NLSize)
    {
#ifdef IPS_VIRIAL
        PMEForce vir_11                                 = (PMEForce)0;
        PMEForce vir_22                                 = (PMEForce)0;
        PMEForce vir_33                                 = (PMEForce)0;                
#endif

        // Read Neighbor List entry
        unsigned int tgx                                = threadIdx.x & (GRID - 1);
        unsigned int tbx                                = threadIdx.x - tgx;
#if (__CUDA_ARCH__ < 300)
        volatile NLAtom* psA                            = &sAtom[tbx];
#endif
        volatile NLForce* psF                           = &sForce[tbx];
        psWarp->nlEntry.array[tgx]                      = cSim.pNLEntry[psWarp->pos].array[tgx];
        if (tgx == 0)
            psWarp->offset                              = cSim.pNLOffset[psWarp->pos];
            
        psWarp->ypos                                    = psWarp->nlEntry.NL.yAtom;
        psWarp->yend                                    = psWarp->nlEntry.NL.yEnd >> NLCELLSHIFT;
        psWarp->nlpos                                   = 0;
        psWarp->homeCellBuffer                          = ((psWarp->nlEntry.NL.xyBufferOffset >> NLENTRYXBUFFEROFFSETSHIFT) & NLENTRYXBUFFEROFFSETMASK) + cSim.NLHomeCellBuffer;
        
        // Calculate output buffer offset based on position
        psWarp->bufferOffset                            = NEIGHBORCELLS * (psWarp->nlEntry.NL.xyBufferOffset & NLENTRYYBUFFEROFFSETMASK) * cSim.stride3;    
   //     if (tgx == 0)
   //         printf("%06d 0x%08x %06d %06d\n", psWarp->pos, psWarp->nlEntry.NL.xyBufferOffset, (psWarp->nlEntry.NL.xyBufferOffset & NLENTRYYBUFFEROFFSETMASK), ((psWarp->nlEntry.NL.xyBufferOffset >> NLENTRYXBUFFEROFFSETSHIFT) & NLENTRYXBUFFEROFFSETMASK));

        while (psWarp->ypos < psWarp->yend)
        {                
            // Read y atoms into registers
            psWarp->ymax                                = min(psWarp->ypos + cSim.NLAtomsPerWarp, psWarp->yend);
            PMEFloat xi;
            PMEFloat yi;
            PMEFloat zi;
            PMEFloat qi;
            PMEFloat sigi;
            PMEFloat epsi;
            unsigned int exclusion;
            PMEForce fx_i                               = (PMEForce)0;
            PMEForce fy_i                               = (PMEForce)0;
            PMEForce fz_i                               = (PMEForce)0;
            unsigned int index                          = psWarp->ypos + (tgx & cSim.NLAtomsPerWarpBitsMask);            
            if (index < psWarp->ymax)
            {
                PMEFloat2 xy                            = cSim.pAtomXYSP[index];  
                PMEFloat2 sigeps                        = cSim.pImageSigEps[index];                       
                zi                                      = cSim.pAtomZSP[index];   
                qi                                      = cSim.pAtomChargeSP[index];                   
                xi                                      = xy.x;
                yi                                      = xy.y;
                sigi                                    = sigeps.x;
                epsi                                    = sigeps.y;
            }
            else
            {
                xi                                      = (PMEFloat)10000.0 * index;
                yi                                      = (PMEFloat)10000.0 * index;
                zi                                      = (PMEFloat)10000.0 * index;
                qi                                      = (PMEFloat)0.0;
                sigi                                    = (PMEFloat)0.0;
                epsi                                    = (PMEFloat)0.0;
            }

#ifndef IPS_IS_ORTHOGONAL   
            // Transform into cartesian space     
#ifdef IPS_VIRIAL
            xi                                          = sUcellf[0] * xi + sUcellf[1] * yi + sUcellf[2] * zi;
            yi                                          =                   sUcellf[4] * yi + sUcellf[5] * zi;
            zi                                          =                                     sUcellf[8] * zi;              
#else                
            xi                                          = cSim.ucell[0][0] * xi + cSim.ucell[0][1] * yi + cSim.ucell[0][2] * zi;
            yi                                          =                         cSim.ucell[1][1] * yi + cSim.ucell[1][2] * zi;
            zi                                          =                                                 cSim.ucell[2][2] * zi;
#endif
#endif 
            // Special-case first tile
            // Copy register data into shared memory
            if (psWarp->nlEntry.NL.xyBufferOffset & NLENTRYHOMECELLMASK)
            {
                exclusion                               = cSim.pNLAtomList[psWarp->offset + (tgx & cSim.NLAtomsPerWarpBitsMask)];    
                psWarp->offset                         += cSim.NLAtomsPerWarp;
#if (IPS_ATOMS_PER_WARP == 16)
                exclusion                             >>= 8 * (tgx >> cSim.NLAtomsPerWarpBits); // SIZE DEPENDENT 2, 8
#endif                        
                PSATOMX(tgx)                            = xi;
                PSATOMY(tgx)                            = yi;
                PSATOMZ(tgx)                            = zi;
                PSATOMQ(tgx)                            = qi;
                PSATOMSIG(tgx)                          = sigi;
                PSATOMEPS(tgx)                          = epsi;
                
                // Set up iteration counts
#if (IPS_ATOMS_PER_WARP == 32)
                unsigned int j                          = sNext[tgx];
#if (__CUDA_ARCH__ >= 300)
                unsigned int shIdx                      = j;
                shAtom.x                                = __shfl(shAtom.x, j);
                shAtom.y                                = __shfl(shAtom.y, j);
                shAtom.z                                = __shfl(shAtom.z, j);
                shAtom.q                                = __shfl(shAtom.q, j);
                shAtom.sig                              = __shfl(shAtom.sig, j);
                shAtom.eps                              = __shfl(shAtom.eps, j);
#endif
                while (j != tgx)
#else
                unsigned int j                          = sStart[tgx];
                unsigned int end                        = sEnd[tgx];
#if (__CUDA_ARCH__ >= 300)
                unsigned int shIdx                      = sNext[tgx];
                shAtom.x                                = __shfl(shAtom.x, j);
                shAtom.y                                = __shfl(shAtom.y, j);
                shAtom.z                                = __shfl(shAtom.z, j);
                shAtom.q                                = __shfl(shAtom.q, j);
                shAtom.sig                              = __shfl(shAtom.sig, j);
                shAtom.eps                              = __shfl(shAtom.eps, j);            
#endif
                while (j != end)
#endif                
                {
                    PMEFloat xij                        = xi - PSATOMX(j);
                    PMEFloat yij                        = yi - PSATOMY(j);
                    PMEFloat zij                        = zi - PSATOMZ(j);          
                    PMEFloat r2                         = xij * xij + yij * yij + zij * zij;                            
                                
                    if (r2 < cSim.cut2)
                    {   
                        PMEFloat rinv                   = rsqrt(r2);
                        PMEFloat sig                    = sigi + PSATOMSIG(j);
                        PMEFloat r2inv                  = rinv * rinv;
                        PMEFloat eps                    = epsi * PSATOMEPS(j);
                        sig                            *= sig * sig;
                        sig                            *= sig;
                        PMEFloat f6                     = eps * sig;
                        PMEFloat f12                    = f6 * sig;
                        PMEFloat b0                     = qi * PSATOMQ(j);
                        PMEFloat r4                     = r2 * r2;
                        PMEFloat r6inv                  = r2inv * r2inv * r2inv;
                        PMEFloat r12inv                 = r6inv * r6inv;       
                        PMEFloat dpipse                 = -r2 * (cSim.bipse1  + r2 * (cSim.bipse2  + r2 * cSim.bipse3)); 
                        PMEFloat dvcu                   = -r2 * (cSim.bipsvc1 + r2 * (cSim.bipsvc2 + r2 * cSim.bipsvc3));            
                        PMEFloat dvau                   = -r4 * (cSim.bipsva1 + r4 * (cSim.bipsva2 + r4 * cSim.bipsva3));       
#ifdef IPS_ENERGY              
                        PMEFloat pipse                  = cSim.aipse0  + r2 * (cSim.aipse1  + r2 * (cSim.aipse2  + r2 * cSim.aipse3))  - cSim.pipsec;  
                        PMEFloat pvc                    = cSim.aipsvc0 + r2 * (cSim.aipsvc1 + r2 * (cSim.aipsvc2 + r2 * cSim.aipsvc3)) - cSim.pipsvcc;                
                        PMEFloat pva                    = cSim.aipsva0 + r4 * (cSim.aipsva1 + r4 * (cSim.aipsva2 + r4 * cSim.aipsva3)) - cSim.pipsvac;                  
#endif                            
                        if (!(exclusion & 0x1))
                        {                                        
                            dpipse                     += rinv; 
                            dvcu                       += (PMEFloat)6.0  * r6inv;            
                            dvau                       += (PMEFloat)12.0 * r12inv;       
#ifdef IPS_ENERGY              
                            pipse                      += rinv;  
                            pvc                        += r6inv;                
                            pva                        += r12inv;  
#endif                       
                        }
#ifdef IPS_ENERGY
#ifdef use_SPFP
                        evdw                           += lliroundf((PMEFloat)0.5 * ENERGYSCALEF * (f12 * pva - f6 * pvc));                       
                        eed                            += lliroundf((PMEFloat)0.5 * ENERGYSCALEF * b0 * pipse); 
#else
                        evdw                           += (double)((PMEFloat)0.5 * (f12 * pva - f6 * pvc));                       
                        eed                            += (double)((PMEFloat)0.5 * b0 * pipse);                 
#endif
#endif 
                        PMEFloat df                     = (f12 * dvau - f6 * dvcu + b0 * dpipse) * r2inv;
#ifdef use_SPFP
                        df                             *= FORCESCALEF;
#endif
                        PMEFloat dfdx                   = df * xij;
                        PMEFloat dfdy                   = df * yij;
                        PMEFloat dfdz                   = df * zij;

		                // Accumulate into registers only
#ifdef use_SPFP
                        fx_i                           += lliroundf(dfdx);
                        fy_i                           += lliroundf(dfdy);
                        fz_i                           += lliroundf(dfdz);
#ifdef IPS_VIRIAL                   
                        vir_11                         -= lliroundf((PMEFloat)0.5 * xij * dfdx);
                        vir_22                         -= lliroundf((PMEFloat)0.5 * yij * dfdy);
                        vir_33                         -= lliroundf((PMEFloat)0.5 * zij * dfdz);
#endif
#else
                        fx_i                           += (PMEDouble)dfdx;
                        fy_i                           += (PMEDouble)dfdy;
                        fz_i                           += (PMEDouble)dfdz;
#ifdef IPS_VIRIAL                   
                        vir_11                         -= (PMEDouble)((PMEFloat)0.5 * xij * dfdx);
                        vir_22                         -= (PMEDouble)((PMEFloat)0.5 * yij * dfdy);
                        vir_33                         -= (PMEDouble)((PMEFloat)0.5 * zij * dfdz);
#endif
#endif                  
                    }                                  
                    exclusion                         >>= 1;
#if (__CUDA_ARCH__ >= 300)
                    shAtom.x                            = __shfl(shAtom.x, shIdx);
                    shAtom.y                            = __shfl(shAtom.y, shIdx);
                    shAtom.z                            = __shfl(shAtom.z, shIdx);
                    shAtom.q                            = __shfl(shAtom.q, shIdx);
                    shAtom.sig                          = __shfl(shAtom.sig, shIdx);
                    shAtom.eps                          = __shfl(shAtom.eps, shIdx);
#endif
                    j                                   = sNext[j];
                }
            }

            // Handle remainder of line
            int tx                                      = 0;
            while (tx < psWarp->nlEntry.NL.xAtoms[psWarp->nlpos])
            {

                // Read atom ID and exclusion data
                PSATOMID(tgx)                           = cSim.pNLAtomList[psWarp->offset + tgx];
                psWarp->offset                         += GRID;
                exclusion                               = cSim.pNLAtomList[psWarp->offset + (tgx & cSim.NLAtomsPerWarpBitsMask)];
                psWarp->offset                         += cSim.NLAtomsPerWarp;  
#if (IPS_ATOMS_PER_WARP == 16)                       
                exclusion                             >>= cSim.NLAtomsPerWarp * (tgx >> cSim.NLAtomsPerWarpBits);
#endif

                // Clear shared memory forces
#ifdef use_SPFP
                psF[tgx].x                              = 0;
                psF[tgx].y                              = 0;
                psF[tgx].z                              = 0;
#else
                psF[tgx].x                              = (PMEDouble)0.0;
                psF[tgx].y                              = (PMEDouble)0.0;
                psF[tgx].z                              = (PMEDouble)0.0;
#endif
                // Read shared memory data
                if (tx + tgx < psWarp->nlEntry.NL.xAtoms[psWarp->nlpos])
                {      
                    unsigned int atom                   = PSATOMID(tgx) >> NLCELLSHIFT;
#ifndef use_DPDP
                    PMEFloat2 xy                        = tex1Dfetch(xytexref, atom);
                    PMEFloat2 sigeps                    = tex1Dfetch(sigepstexref, atom);
                    PSATOMZ(tgx)                        = tex1Dfetch(ztexref, atom);                 
                    PSATOMQ(tgx)                        = tex1Dfetch(qtexref, atom);
#else
                    PMEFloat2 xy                        = cSim.pAtomXYSP[atom];  
                    PMEFloat2 sigeps                    = cSim.pImageSigEps[atom];
                    PSATOMZ(tgx)                        = cSim.pAtomZSP[atom];                 
                    PSATOMQ(tgx)                        = cSim.pAtomChargeSP[atom]; 
#endif    
                    PSATOMX(tgx)                        = xy.x;  
                    PSATOMY(tgx)                        = xy.y;
                    PSATOMSIG(tgx)                      = sigeps.x;
                    PSATOMEPS(tgx)                      = sigeps.y;
                } 
                else
                {
                    PSATOMX(tgx)                        = (PMEFloat)-10000.0 * tgx;
                    PSATOMY(tgx)                        = (PMEFloat)-10000.0 * tgx;
                    PSATOMZ(tgx)                        = (PMEFloat)-10000.0 * tgx;
                    PSATOMQ(tgx)                        = (PMEFloat)0.0;
                    PSATOMSIG(tgx)                      = (PMEFloat)0.0;
                    PSATOMEPS(tgx)                      = (PMEFloat)0.0;                    
                }   
                                 
                // Translate all atoms into a local coordinate system within one unit
                // cell of the first atom read to avoid PBC handling within inner loops  
                int cell                                = PSATOMID(tgx) & NLCELLTYPEMASK;                 
#if defined(IPS_VIRIAL) && defined(IPS_IS_ORTHOGONAL)
                PSATOMX(tgx)                           += sUcellf[0] * cSim.cellOffset[cell][0];
                PSATOMY(tgx)                           += sUcellf[4] * cSim.cellOffset[cell][1];
                PSATOMZ(tgx)                           += sUcellf[8] * cSim.cellOffset[cell][2];
#else
                PSATOMX(tgx)                           += cSim.cellOffset[cell][0];
                PSATOMY(tgx)                           += cSim.cellOffset[cell][1];
                PSATOMZ(tgx)                           += cSim.cellOffset[cell][2];
#endif

#ifndef IPS_IS_ORTHOGONAL
#ifdef IPS_VIRIAL
                PSATOMX(tgx)                            = sUcellf[0] * PSATOMX(tgx) + sUcellf[1] * PSATOMY(tgx) + sUcellf[2] * PSATOMZ(tgx);
                PSATOMY(tgx)                            =                             sUcellf[4] * PSATOMY(tgx) + sUcellf[5] * PSATOMZ(tgx);
                PSATOMZ(tgx)                            =                                                         sUcellf[8] * PSATOMZ(tgx);
#else
                PSATOMX(tgx)                            = cSim.ucell[0][0] * PSATOMX(tgx) + cSim.ucell[0][1] * PSATOMY(tgx) + cSim.ucell[0][2] * PSATOMZ(tgx);
                PSATOMY(tgx)                            =                                   cSim.ucell[1][1] * PSATOMY(tgx) + cSim.ucell[1][2] * PSATOMZ(tgx);
                PSATOMZ(tgx)                            =                                                                     cSim.ucell[2][2] * PSATOMZ(tgx);
#endif
#endif                               
                int j                                   = tgx;
#if (__CUDA_ARCH__ >= 300)
                int shIdx                               = sNext[tgx];
#endif
                if (__any(exclusion))
                {
                    do 
                    {
                        PMEFloat xij                    = xi - PSATOMX(j); 
                        PMEFloat yij                    = yi - PSATOMY(j); 
                        PMEFloat zij                    = zi - PSATOMZ(j);                        
                        PMEFloat r2                     = xij * xij + yij * yij + zij * zij;                            
          
                        if (r2 < cSim.cut2)
                        {                                       
                            PMEFloat rinv               = rsqrt(r2);
                            PMEFloat sig                = sigi + PSATOMSIG(j);
                            PMEFloat r2inv              = rinv * rinv;
                            PMEFloat eps                = epsi * PSATOMEPS(j);
                            sig                        *= sig * sig;
                            sig                        *= sig;
                            PMEFloat f6                 = eps * sig;
                            PMEFloat f12                = f6 * sig;
                            PMEFloat b0                 = qi * PSATOMQ(j);
                            PMEFloat r4                 = r2 * r2;
                            PMEFloat r6inv              = r2inv * r2inv * r2inv;
                            PMEFloat r12inv             = r6inv * r6inv;       
                            PMEFloat dpipse             = -r2 * (cSim.bipse1  + r2 * (cSim.bipse2  + r2 * cSim.bipse3)); 
                            PMEFloat dvcu               = -r2 * (cSim.bipsvc1 + r2 * (cSim.bipsvc2 + r2 * cSim.bipsvc3));            
                            PMEFloat dvau               = -r4 * (cSim.bipsva1 + r4 * (cSim.bipsva2 + r4 * cSim.bipsva3));       
#ifdef IPS_ENERGY              
                            PMEFloat pipse              = cSim.aipse0  + r2 * (cSim.aipse1  + r2 * (cSim.aipse2  + r2 * cSim.aipse3)) - cSim.pipsec;  
                            PMEFloat pvc                = cSim.aipsvc0 + r2 * (cSim.aipsvc1 + r2 * (cSim.aipsvc2 + r2 * cSim.aipsvc3)) - cSim.pipsvcc;                
                            PMEFloat pva                = cSim.aipsva0 + r4 * (cSim.aipsva1 + r4 * (cSim.aipsva2 + r4 * cSim.aipsva3)) - cSim.pipsvac;                  
#endif                            
                            if (!(exclusion & 0x1))
                            {                                        
                                dpipse                 += rinv; 
                                dvcu                   += (PMEFloat)6.0  * r6inv;            
                                dvau                   += (PMEFloat)12.0 * r12inv;       
#ifdef IPS_ENERGY              
                                pipse                  += rinv;  
                                pvc                    += r6inv;                
                                pva                    += r12inv;  
#endif                       
                            }
#ifdef IPS_ENERGY
#ifdef use_SPFP
                            evdw                       += lliroundf(ENERGYSCALEF * (f12 * pva - f6 * pvc));                       
                            eed                        += lliroundf(ENERGYSCALEF * b0 * pipse); 
#else
                            evdw                       += (PMEForce)(f12 * pva - f6 * pvc);                       
                            eed                        += (PMEForce)(b0 * pipse);      
#endif
#endif                                                                                                                                                               
                            PMEFloat df                 = (f12 * dvau - f6 * dvcu + b0 * dpipse) * r2inv;
#ifdef use_SPFP
                            df                         *= FORCESCALEF;
#endif
                            PMEFloat dfdx               = df * xij;
                            PMEFloat dfdy               = df * yij;
                            PMEFloat dfdz               = df * zij;
#ifdef use_SPFP
                            long long int dfdx1         = lliroundf(dfdx);
                            long long int dfdy1         = lliroundf(dfdy);
                            long long int dfdz1         = lliroundf(dfdz);
                            fx_i                       += dfdx1;
                            fy_i                       += dfdy1;
                            fz_i                       += dfdz1;
                            psF[j].x                   -= dfdx1;
                            psF[j].y                   -= dfdy1;
                            psF[j].z                   -= dfdz1;
#ifdef IPS_VIRIAL                   
                            vir_11                     -= lliroundf(xij * dfdx);
                            vir_22                     -= lliroundf(yij * dfdy);
                            vir_33                     -= lliroundf(zij * dfdz);
#endif
#else
                            PMEForce dfdx1              = (PMEForce)dfdx;
                            PMEForce dfdy1              = (PMEForce)dfdy;
                            PMEForce dfdz1              = (PMEForce)dfdz;
                            fx_i                       += dfdx1;
                            fy_i                       += dfdy1;
                            fz_i                       += dfdz1;
                            psF[j].x                   -= dfdx1;
                            psF[j].y                   -= dfdy1;
                            psF[j].z                   -= dfdz1;
#ifdef IPS_VIRIAL                   
                            vir_11                     -= (PMEForce)(xij * dfdx);
                            vir_22                     -= (PMEForce)(yij * dfdy);
                            vir_33                     -= (PMEForce)(zij * dfdz);
#endif
#endif                                       
                        }
                        exclusion                     >>= 1;
#if (__CUDA_ARCH__ >= 300)
                        shAtom.x                        = __shfl(shAtom.x, shIdx);
                        shAtom.y                        = __shfl(shAtom.y, shIdx);
                        shAtom.z                        = __shfl(shAtom.z, shIdx);
                        shAtom.q                        = __shfl(shAtom.q, shIdx);
                        shAtom.sig                      = __shfl(shAtom.sig, shIdx);
                        shAtom.eps                      = __shfl(shAtom.eps, shIdx);
#endif
                        j                               = sNext[j];   
                    }
                    while (j != tgx);
                }
                else
                {
                    do 
                    {
                        PMEFloat xij                    = xi - PSATOMX(j); 
                        PMEFloat yij                    = yi - PSATOMY(j); 
                        PMEFloat zij                    = zi - PSATOMZ(j);                      
                        PMEFloat r2                     = xij * xij + yij * yij + zij * zij;                            
          
                        if (r2 < cSim.cut2)
                        {                                       
                            PMEFloat rinv               = rsqrt(r2);
                            PMEFloat sig                = sigi + PSATOMSIG(j);
                            PMEFloat r2inv              = rinv * rinv;
                            PMEFloat eps                = epsi * PSATOMEPS(j);
                            sig                        *= sig * sig;
                            sig                        *= sig;
                            PMEFloat f6                 = eps * sig;
                            PMEFloat f12                = f6 * sig;
                            PMEFloat b0                 = qi * PSATOMQ(j);
                            PMEFloat r4                 = r2 * r2;
                            PMEFloat r6inv              = r2inv * r2inv * r2inv;
                            PMEFloat r12inv             = r6inv * r6inv;       
                            PMEFloat dpipse             = rinv                    - r2 * (cSim.bipse1  + r2 * (cSim.bipse2  + r2 * cSim.bipse3)); 
                            PMEFloat dvcu               = (PMEFloat)6.0  * r6inv  - r2 * (cSim.bipsvc1 + r2 * (cSim.bipsvc2 + r2 * cSim.bipsvc3));            
                            PMEFloat dvau               = (PMEFloat)12.0 * r12inv - r4 * (cSim.bipsva1 + r4 * (cSim.bipsva2 + r4 * cSim.bipsva3));  
                            PMEFloat df                 = (f12 * dvau - f6 * dvcu + b0 * dpipse) * r2inv;                                 
#ifdef IPS_ENERGY              
                            PMEFloat pipse              = rinv   + cSim.aipse0  + r2 * (cSim.aipse1  + r2 * (cSim.aipse2  + r2 * cSim.aipse3))  - cSim.pipsec;  
                            PMEFloat pvc                = r6inv  + cSim.aipsvc0 + r2 * (cSim.aipsvc1 + r2 * (cSim.aipsvc2 + r2 * cSim.aipsvc3)) - cSim.pipsvcc;                
                            PMEFloat pva                = r12inv + cSim.aipsva0 + r4 * (cSim.aipsva1 + r4 * (cSim.aipsva2 + r4 * cSim.aipsva3)) - cSim.pipsvac;  
#ifdef use_SPFP
                            evdw                       += lliroundf(ENERGYSCALEF * (f12 * pva - f6 * pvc));                       
                            eed                        += lliroundf(ENERGYSCALEF * b0 * pipse); 
#else    
                            evdw                       += (PMEForce)(f12 * pva - f6 * pvc);                       
                            eed                        += (PMEForce)(b0 * pipse); 
#endif                                               
#endif
#ifdef use_SPFP
                            df                         *= FORCESCALEF;
#endif                                                       
                            PMEFloat dfdx               = df * xij;
                            PMEFloat dfdy               = df * yij;
                            PMEFloat dfdz               = df * zij;
#ifdef use_SPFP
                            long long int dfdx1         = lliroundf(dfdx);
                            long long int dfdy1         = lliroundf(dfdy);
                            long long int dfdz1         = lliroundf(dfdz);
                            fx_i                       += dfdx1;
                            fy_i                       += dfdy1;
                            fz_i                       += dfdz1;
                            psF[j].x                   -= dfdx1;
                            psF[j].y                   -= dfdy1;
                            psF[j].z                   -= dfdz1;
#ifdef IPS_VIRIAL                   
                            vir_11                     -= lliroundf(xij * dfdx);
                            vir_22                     -= lliroundf(yij * dfdy);
                            vir_33                     -= lliroundf(zij * dfdz);
#endif
#else
                            PMEForce dfdx1              = (PMEForce)dfdx;
                            PMEForce dfdy1              = (PMEForce)dfdy;
                            PMEForce dfdz1              = (PMEForce)dfdz;
                            fx_i                       += dfdx1;
                            fy_i                       += dfdy1;
                            fz_i                       += dfdz1;
                            psF[j].x                   -= dfdx1;
                            psF[j].y                   -= dfdy1;
                            psF[j].z                   -= dfdz1;
#ifdef IPS_VIRIAL                   
                            vir_11                     -= (PMEForce)(xij * dfdx);
                            vir_22                     -= (PMEForce)(yij * dfdy);
                            vir_33                     -= (PMEForce)(zij * dfdz);
#endif
#endif  
                        }
#if (__CUDA_ARCH__ >= 300)
                        shAtom.x                        = __shfl(shAtom.x, shIdx);
                        shAtom.y                        = __shfl(shAtom.y, shIdx);
                        shAtom.z                        = __shfl(shAtom.z, shIdx);
                        shAtom.q                        = __shfl(shAtom.q, shIdx);
                        shAtom.sig                      = __shfl(shAtom.sig, shIdx);
                        shAtom.eps                      = __shfl(shAtom.eps, shIdx);
#endif
                        j                               = sNext[j];
                    }
                    while (j != tgx);                
                }
                      
                // Dump shared memory forces
                if (tx + tgx < psWarp->nlEntry.NL.xAtoms[psWarp->nlpos])
                {
#ifdef use_SPFP
                    int offset                          = (PSATOMID(tgx) >> NLCELLSHIFT);    
                    atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset], llitoulli(psF[tgx].x));
                    atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset], llitoulli(psF[tgx].y));
                    atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset], llitoulli(psF[tgx].z));
#else
                    int offset                          = (PSATOMID(tgx) >> NLCELLSHIFT) + (PSATOMID(tgx) & NLCELLTYPEMASK) * cSim.stride3 + psWarp->bufferOffset;    
                    cSim.pForceXBuffer[offset]         += (PMEDouble)psF[tgx].x;
                    cSim.pForceYBuffer[offset]         += (PMEDouble)psF[tgx].y;
                    cSim.pForceZBuffer[offset]         += (PMEDouble)psF[tgx].z;
#endif
                }   
                                
                // Advance to next x tile
                tx                                     += GRID;
            }
            
            // Reduce register forces if necessary
            if (cSim.NLAtomsPerWarp == 16)
            {
                psF[tgx].x                              = fx_i;
                psF[tgx].y                              = fy_i;
                psF[tgx].z                              = fz_i;
                if (tgx < 16)
                {
                    fx_i                               += psF[tgx + 16].x;
                    fy_i                               += psF[tgx + 16].y;
                    fz_i                               += psF[tgx + 16].z;
                }
            }
            
            // Dump register forces
            if (psWarp->ypos + tgx < psWarp->ymax)
            {
#ifdef use_SPFP 
                int offset                              = psWarp->ypos + tgx;  
                atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset], llitoulli(fx_i));
                atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset], llitoulli(fy_i));
                atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset], llitoulli(fz_i));
#else
                int offset                              = psWarp->ypos + tgx + psWarp->homeCellBuffer * cSim.stride3;  
                cSim.pForceXBuffer[offset]             += fx_i;
                cSim.pForceYBuffer[offset]             += fy_i;
                cSim.pForceZBuffer[offset]             += fz_i;
#endif
            }                
         
             // Advance to next set of register atoms
            psWarp->ypos                               += cSim.NLYStride;
            psWarp->nlpos++;                       
        }        
        
#ifdef IPS_VIRIAL
        // Reduce virial per warp and convert to fixed point if necessary
        volatile NLVirial* psV                          = &sWarpVirial[threadIdx.x / GRID];
        psF[tgx].x                                      = vir_11;
        psF[tgx].y                                      = vir_22;
        psF[tgx].z                                      = vir_33;
        if (tgx < 16)
        {
            psF[tgx].x                                 += psF[tgx + 16].x; 
            psF[tgx].y                                 += psF[tgx + 16].y; 
            psF[tgx].z                                 += psF[tgx + 16].z; 
        }
        if (tgx < 8)
        {
            psF[tgx].x                                 += psF[tgx + 8].x; 
            psF[tgx].y                                 += psF[tgx + 8].y; 
            psF[tgx].z                                 += psF[tgx + 8].z; 
        }
        if (tgx < 4)
        {
            psF[tgx].x                                 += psF[tgx + 4].x; 
            psF[tgx].y                                 += psF[tgx + 4].y; 
            psF[tgx].z                                 += psF[tgx + 4].z; 
        }
        if (tgx < 2)
        {
            psF[tgx].x                                 += psF[tgx + 2].x; 
            psF[tgx].y                                 += psF[tgx + 2].y; 
            psF[tgx].z                                 += psF[tgx + 2].z; 
        }
        if (tgx == 0)
        {
#ifdef use_SPFP
            psV->vir_11                                += psF[tgx].x + psF[tgx + 1].x; 
            psV->vir_22                                += psF[tgx].y + psF[tgx + 1].y; 
            psV->vir_33                                += psF[tgx].z + psF[tgx + 1].z;
#else
            psV->vir_11                                += lliroundf((psF[tgx].x + psF[tgx + 1].x) * FORCESCALE); 
            psV->vir_22                                += lliroundf((psF[tgx].y + psF[tgx + 1].y) * FORCESCALE); 
            psV->vir_33                                += lliroundf((psF[tgx].z + psF[tgx + 1].z) * FORCESCALE);
#endif
        }
#endif

        // Get next Neighbor List entry
        if (tgx == 0)
            psWarp->pos                                 = atomicAdd(cSim.pNLPosition, 1);
    }
              
    // Do virial and energy reductions if necessary
#ifdef IPS_VIRIAL
    __syncthreads();
#if (__CUDA_ARCH__ >= 300) && defined(use_DPDP) // Reduce from 1024 to 512 threads
    if (threadIdx.x < 16)
    {
        sWarpVirial[threadIdx.x].vir_11                += sWarpVirial[threadIdx.x + 16].vir_11;
        sWarpVirial[threadIdx.x].vir_22                += sWarpVirial[threadIdx.x + 16].vir_22;
        sWarpVirial[threadIdx.x].vir_33                += sWarpVirial[threadIdx.x + 16].vir_33;          
    }
#elif (__CUDA_ARCH__ >= 300) // Reduce from 640 to 512 threads
    if (threadIdx.x < 4)
    {
        sWarpVirial[threadIdx.x].vir_11                += sWarpVirial[threadIdx.x + 16].vir_11;
        sWarpVirial[threadIdx.x].vir_22                += sWarpVirial[threadIdx.x + 16].vir_22;
        sWarpVirial[threadIdx.x].vir_33                += sWarpVirial[threadIdx.x + 16].vir_33;  
    }
#elif (__CUDA_ARCH__ >= 200) && !defined(use_DPDP) // Reduce from 768 to 512 threads
    if (threadIdx.x < 8)
    {
        sWarpVirial[threadIdx.x].vir_11                += sWarpVirial[threadIdx.x + 16].vir_11;
        sWarpVirial[threadIdx.x].vir_22                += sWarpVirial[threadIdx.x + 16].vir_22;
        sWarpVirial[threadIdx.x].vir_33                += sWarpVirial[threadIdx.x + 16].vir_33;   
    }
#endif

#if (__CUDA_ARCH__ >= 200) // Reduce from 512 to 256 threads
    if (threadIdx.x < 8)
    {
        sWarpVirial[threadIdx.x].vir_11                += sWarpVirial[threadIdx.x + 8].vir_11;
        sWarpVirial[threadIdx.x].vir_22                += sWarpVirial[threadIdx.x + 8].vir_22;
        sWarpVirial[threadIdx.x].vir_33                += sWarpVirial[threadIdx.x + 8].vir_33;  
    }
#endif
#if (__CUDA_ARCH__ >= 200) || !defined(use_DPDP) // Reduce from 256 to 128 threads
    if (threadIdx.x < 4)
    {
        sWarpVirial[threadIdx.x].vir_11                += sWarpVirial[threadIdx.x + 4].vir_11;
        sWarpVirial[threadIdx.x].vir_22                += sWarpVirial[threadIdx.x + 4].vir_22;
        sWarpVirial[threadIdx.x].vir_33                += sWarpVirial[threadIdx.x + 4].vir_33;   
    }
#else
    if (threadIdx.x == 0)   // Reduce from 160 to 128 thrads
    {
        sWarpVirial[threadIdx.x].vir_11                += sWarpVirial[threadIdx.x + 4].vir_11;
        sWarpVirial[threadIdx.x].vir_22                += sWarpVirial[threadIdx.x + 4].vir_22;
        sWarpVirial[threadIdx.x].vir_33                += sWarpVirial[threadIdx.x + 4].vir_33;   
    }    
#endif
    if (threadIdx.x < 2)
    {
        sWarpVirial[threadIdx.x].vir_11                += sWarpVirial[threadIdx.x + 2].vir_11;
        sWarpVirial[threadIdx.x].vir_22                += sWarpVirial[threadIdx.x + 2].vir_22;
        sWarpVirial[threadIdx.x].vir_33                += sWarpVirial[threadIdx.x + 2].vir_33;   
    }
    if (threadIdx.x < 1)
    {
        sWarpVirial[threadIdx.x].vir_11                += sWarpVirial[threadIdx.x + 1].vir_11;
        sWarpVirial[threadIdx.x].vir_22                += sWarpVirial[threadIdx.x + 1].vir_22;
        sWarpVirial[threadIdx.x].vir_33                += sWarpVirial[threadIdx.x + 1].vir_33;   
    } 
#endif   

#ifdef IPS_ENERGY
#define sEED(i)  sForce[i].x
#define sEVDW(i) sForce[i].y
    sEED(threadIdx.x)                                   = eed;
    sEVDW(threadIdx.x)                                  = evdw;
    __syncthreads();
    unsigned int m                                      = 1;
    while (m < blockDim.x)
    {
        int p                                           = threadIdx.x + m;
#ifdef use_SPFP
        long long int dx                                = ((p < blockDim.x) ? sEED(p)  : 0);
        long long int dy                                = ((p < blockDim.x) ? sEVDW(p) : 0);
#else
        double dx                                       = ((p < blockDim.x) ? sEED(p)  : 0.0);
        double dy                                       = ((p < blockDim.x) ? sEVDW(p) : 0.0);
#endif
        __syncthreads();
        sEED(threadIdx.x)                              += dx;
        sEVDW(threadIdx.x)                             += dy;
        __syncthreads();
        m                                              *= 2;
    }
#endif

#ifdef IPS_VIRIAL
    if (threadIdx.x == 0)
    {
        unsigned long long int val1                     = llitoulli(sWarpVirial[0].vir_11);
        unsigned long long int val2                     = llitoulli(sWarpVirial[0].vir_22);
        unsigned long long int val3                     = llitoulli(sWarpVirial[0].vir_33);
        atomicAdd(cSim.pVirial_11, val1);
        atomicAdd(cSim.pVirial_22, val2);  
        atomicAdd(cSim.pVirial_33, val3);
    }
#endif
#ifdef IPS_ENERGY
    if (threadIdx.x == 0)
    {
        eed                                             = sEED(0);
        evdw                                            = sEVDW(0);
#ifdef use_SPFP
        unsigned long long int val4                     = llitoulli(eed);
        unsigned long long int val5                     = llitoulli(evdw);
#else
        unsigned long long int val4                     = llitoulli(lliroundd(eed * ENERGYSCALE));
        unsigned long long int val5                     = llitoulli(lliroundd(evdw * ENERGYSCALE));        
#endif
        atomicAdd(cSim.pEED, val4);
        atomicAdd(cSim.pEVDW, val5);        
    }
#undef sEED
#undef sEVDW    
#endif       
}
