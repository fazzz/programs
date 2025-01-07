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
// #defines: PME_ENERGY, PME_VIRIAL, PME_IS_ORTHOGONAL, PME_ATOMS_PER_WARP 
#define uint unsigned int

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
    PMEDouble x;
    PMEDouble y;
    PMEDouble z;
};

#if defined(PME_ENERGY) && defined(use_SPSP)
struct NLEnergy
{
    double x;
    double y;
};
#endif

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


#if defined(PME_VIRIAL)
__shared__ volatile PMEFloat sUcellf[9];
#endif
__shared__ unsigned int sNext[GRID];
#if (PME_ATOMS_PER_WARP == 16)
__shared__ unsigned int sStart[GRID];
__shared__ unsigned int sEnd[GRID];
#endif
#if (__CUDA_ARCH__ >= 200)
__shared__ volatile NLWarp sWarp[SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK / GRID];
__shared__ volatile NLAtom sAtom[SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK];
__shared__ volatile NLForce sForce[SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK];
#ifdef PME_VIRIAL
__shared__ volatile NLVirial sWarpVirial[SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK / GRID];
#endif
#if defined(PME_ENERGY) && defined(use_SPSP)
__shared__ volatile NLEnergy sEnergy[SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK];
#endif
#else
__shared__ volatile NLWarp sWarp[SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK / GRID];
__shared__ volatile NLAtom sAtom[SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK];
__shared__ volatile NLForce sForce[SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK];
#ifdef PME_VIRIAL
__shared__ volatile NLVirial sWarpVirial[SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK / GRID];
#endif
#if defined(PME_ENERGY) && defined(use_SPSP)
__shared__ volatile NLEnergy sEnergy[SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK];
#endif
#endif

    // Read static data
    if (threadIdx.x < GRID)
    {
        unsigned int offset                             = cSim.NLAtomsPerWarp * (threadIdx.x >> cSim.NLAtomsPerWarpBits);
        sNext[threadIdx.x]                              = ((threadIdx.x + 1) & cSim.NLAtomsPerWarpBitsMask) + offset;

#if (PME_ATOMS_PER_WARP == 16)
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
     
#ifdef PME_VIRIAL
    if (threadIdx.x < 9)
        sUcellf[threadIdx.x]                            = cSim.pNTPData->ucellf[threadIdx.x];
#if (__CUDA_ARCH__ >= 200)        
    if (threadIdx.x < SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK / GRID)
#else
    if (threadIdx.x < SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK / GRID)
#endif
    {
        sWarpVirial[threadIdx.x].vir_11                 = 0;
        sWarpVirial[threadIdx.x].vir_22                 = 0;
        sWarpVirial[threadIdx.x].vir_33                 = 0;                
    }
#endif    

#ifdef PME_ENERGY    
    double eed                                          = 0.0;
    double evdw                                         = 0.0; 
#endif    
    volatile NLWarp* psWarp                             = &sWarp[threadIdx.x / GRID];
    psWarp->pos                                         = (blockIdx.x * blockDim.x + threadIdx.x) / GRID;
    __syncthreads();
 
    while (psWarp->pos < cSim.NLSize)
    { 
#ifdef PME_VIRIAL
        PMEDouble vir_11                                = (PMEDouble)0.0;
        PMEDouble vir_22                                = (PMEDouble)0.0;
        PMEDouble vir_33                                = (PMEDouble)0.0;                
#endif

        // Read Neighbor List entry
        unsigned int tgx                                = threadIdx.x & (GRID - 1);
        unsigned int tbx                                = threadIdx.x - tgx;
        volatile NLAtom* psA                            = &sAtom[tbx];
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
            PMEDouble fx_i                              = (PMEDouble)0.0;
            PMEDouble fy_i                              = (PMEDouble)0.0;
            PMEDouble fz_i                              = (PMEDouble)0.0;
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

#ifndef PME_IS_ORTHOGONAL   
            // Transform into cartesian space     
#ifdef PME_VIRIAL
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
#if (PME_ATOMS_PER_WARP == 16)
                exclusion                             >>= 8 * (tgx >> cSim.NLAtomsPerWarpBits); // SIZE DEPENDENT 2, 8
#endif            
            
                psA[tgx].x                              = xi;
                psA[tgx].y                              = yi;
                psA[tgx].z                              = zi;
                psA[tgx].q                              = qi;
                psA[tgx].sig                            = sigi;
                psA[tgx].eps                            = epsi; 
                
                // Set up iteration counts
#if (PME_ATOMS_PER_WARP == 32)
                unsigned int j                          = sNext[tgx];
                while (j != tgx)
#else
                unsigned int j                          = sStart[tgx];
                unsigned int end                        = sEnd[tgx];
                while (j != end)
#endif                
                {
                    PMEFloat xij                        = xi - psA[j].x; 
                    PMEFloat yij                        = yi - psA[j].y; 
                    PMEFloat zij                        = zi - psA[j].z;              
                    PMEFloat r2                         = xij * xij + yij * yij + zij * zij;                            
                                
                    if (r2 < cSim.cut2)
                    {    
                        PMEFloat rinv                   = rsqrt(r2);
                        PMEFloat eps                    = epsi * psA[j].eps;
                        PMEFloat r                      = r2 * rinv;   
                        PMEFloat r2inv                  = rinv * rinv;
                        PMEFloat r6inv                  = r2inv * r2inv * r2inv;
                        PMEFloat sig                    = sigi + psA[j].sig;
                        sig                            *= sig * sig;
                        sig                            *= sig * r6inv;
                        PMEFloat f6                     = eps * sig;
                        PMEFloat f12                    = f6 * sig;
                        PMEFloat df                     = (PMEFloat)0.0;
                        PMEFloat qiqj                   = qi * psA[j].q;         
#ifdef use_DPDP
                        PMEFloat swtch                  = erfc(cSim.ew_coeff * r);
#else                                     
                        PMEFloat swtch                  = fasterfc(cSim.ew_coeff * r);
#endif
                        PMEFloat d_swtch_dx             = cSim.negTwoEw_coeffRsqrtPI * exp(-cSim.ew_coeff2 * r2);                
                                        
                        if (!(exclusion & 0x1))
                        {
                            df                         += ((PMEFloat)12.0 * f12 - (PMEFloat)6.0 * f6) * r2inv;
#ifdef PME_ENERGY
                            evdw                       += (double)((PMEFloat)0.5 * (f12 - f6));
#endif
                        }
                        else
                            swtch                      -= (PMEFloat)1.0;           
#ifdef PME_ENERGY
                        PMEFloat b0                     = qiqj * swtch * rinv;
                        PMEFloat b1                     = b0 - qiqj * d_swtch_dx;
                        df                             += b1 * r2inv;
                        eed                            += (double)(b0 * (PMEFloat)0.5);    
#else
                        df                             += qiqj * (swtch * rinv - d_swtch_dx) * r2inv;                                                  
#endif                         
                        PMEFloat dfdx                   = df * xij;
                        PMEFloat dfdy                   = df * yij;
                        PMEFloat dfdz                   = df * zij;

		                // Accumulate into registers only      
                        fx_i                           += dfdx;
                        fy_i                           += dfdy;
                        fz_i                           += dfdz;
#ifdef PME_VIRIAL                   
                        vir_11                         -= (PMEFloat)0.5 * xij * dfdx;
                        vir_22                         -= (PMEFloat)0.5 * yij * dfdy;
                        vir_33                         -= (PMEFloat)0.5 * zij * dfdz;
#endif                  
                    }                                  
                    exclusion                         >>= 1;
                    j                                   = sNext[j];
                }
            }

            // Handle remainder of line
            int tx                                      = 0;
            while (tx < psWarp->nlEntry.NL.xAtoms[psWarp->nlpos])
            {

                // Read atom ID and exclusion data
                psA[tgx].ID                             = cSim.pNLAtomList[psWarp->offset + tgx];
                psWarp->offset                         += GRID;
                exclusion                               = cSim.pNLAtomList[psWarp->offset + (tgx & cSim.NLAtomsPerWarpBitsMask)];
                psWarp->offset                         += cSim.NLAtomsPerWarp;  
#if (PME_ATOMS_PER_WARP == 16)                       
                exclusion                             >>= cSim.NLAtomsPerWarp * (tgx >> cSim.NLAtomsPerWarpBits);
#endif

                // Clear shared memory forces
                psF[tgx].x                              = (PMEDouble)0.0;
                psF[tgx].y                              = (PMEDouble)0.0;
                psF[tgx].z                              = (PMEDouble)0.0;

                // Read shared memory data
                if (tx + tgx < psWarp->nlEntry.NL.xAtoms[psWarp->nlpos])
                {         
                    unsigned int atom                   = psA[tgx].ID >> NLCELLSHIFT;
                    PMEFloat2 xy                        = cSim.pAtomXYSP[atom];  
                    PMEFloat2 sigeps                    = cSim.pImageSigEps[atom];  
                    psA[tgx].z                          = cSim.pAtomZSP[atom];                 
                    psA[tgx].q                          = cSim.pAtomChargeSP[atom];         
                    psA[tgx].x                          = xy.x;  
                    psA[tgx].y                          = xy.y;
                    psA[tgx].sig                        = sigeps.x;
                    psA[tgx].eps                        = sigeps.y;
                } 
                else
                {
                    psA[tgx].x                          = (PMEFloat)-10000.0 * tgx;
                    psA[tgx].y                          = (PMEFloat)-10000.0 * tgx;
                    psA[tgx].z                          = (PMEFloat)-10000.0 * tgx;
                    psA[tgx].q                          = (PMEFloat)0.0;
                    psA[tgx].sig                        = (PMEFloat)0.0;
                    psA[tgx].eps                        = (PMEFloat)0.0;                    
                }   
                                 
                // Translate all atoms into a local coordinate system within one unit
                // cell of the first atom read to avoid PBC handling within inner loops  
                int cell                                = psA[tgx].ID & NLCELLTYPEMASK;                 
#if defined(PME_VIRIAL) && defined(PME_IS_ORTHOGONAL)
                psA[tgx].x                             += sUcellf[0] * cSim.cellOffset[cell][0];
                psA[tgx].y                             += sUcellf[4] * cSim.cellOffset[cell][1];
                psA[tgx].z                             += sUcellf[8] * cSim.cellOffset[cell][2];
#else
                psA[tgx].x                             += cSim.cellOffset[cell][0];
                psA[tgx].y                             += cSim.cellOffset[cell][1];
                psA[tgx].z                             += cSim.cellOffset[cell][2];
#endif

#ifndef PME_IS_ORTHOGONAL
#ifdef PME_VIRIAL
                psA[tgx].x                              = sUcellf[0] * psA[tgx].x + sUcellf[1] * psA[tgx].y + sUcellf[2] * psA[tgx].z;
                psA[tgx].y                              =                           sUcellf[4] * psA[tgx].y + sUcellf[5] * psA[tgx].z;
                psA[tgx].z                              =                                                     sUcellf[8] * psA[tgx].z;
#else
                psA[tgx].x                              = cSim.ucell[0][0] * psA[tgx].x + cSim.ucell[0][1] * psA[tgx].y + cSim.ucell[0][2] * psA[tgx].z;
                psA[tgx].y                              =                                 cSim.ucell[1][1] * psA[tgx].y + cSim.ucell[1][2] * psA[tgx].z;
                psA[tgx].z                              =                                                                 cSim.ucell[2][2] * psA[tgx].z;
#endif
#endif                                        
                int j                                   = tgx;
                if (__any(exclusion))
                {
                    do 
                    {
                        PMEFloat xij                    = xi - psA[j].x; 
                        PMEFloat yij                    = yi - psA[j].y; 
                        PMEFloat zij                    = zi - psA[j].z;                        
                        PMEFloat r2                     = xij * xij + yij * yij + zij * zij;                            
          
                        if (r2 < cSim.cut2)
                        {                                       
                            PMEFloat rinv               = rsqrt(r2);
                            PMEFloat sig                = sigi + psA[j].sig;
                            PMEFloat r                  = r2 * rinv;      
                            PMEFloat r2inv              = rinv * rinv;
                            PMEFloat r6inv              = r2inv * r2inv * r2inv;
                            PMEFloat eps                = epsi * psA[j].eps;
                            sig                        *= sig * sig;
                            sig                        *= sig * r6inv;
                            PMEFloat f6                 = eps * sig;
                            PMEFloat f12                = f6 * sig;
                            PMEFloat df                 = (PMEFloat)0.0;
                            PMEFloat qiqj               = qi * psA[j].q;                                                
#ifdef use_DPDP
                            PMEFloat swtch              = erfc(cSim.ew_coeff * r);
#else                                     
                            PMEFloat swtch              = fasterfc(cSim.ew_coeff * r);
#endif          
                            PMEFloat d_swtch_dx         = cSim.negTwoEw_coeffRsqrtPI * exp(-cSim.ew_coeff2 * r2);
                            if (!(exclusion & 0x1))
                            {
                                df                     += ((PMEFloat)12.0 * f12 - (PMEFloat)6.0 * f6) * r2inv;
#ifdef PME_ENERGY                            
                                evdw                   += (double)(f12 - f6);
#endif                          
                            }
                            else
                                swtch                  -= (PMEFloat)1.0;
#ifdef PME_ENERGY                                       
                            PMEFloat b0                 = qiqj * swtch * rinv;
                            PMEFloat b1                 = b0 - qiqj * d_swtch_dx;
                            df                         += b1 * r2inv;                                                  
                            eed                        += (double)b0;     
#else
                            df                         += qiqj * (swtch * rinv - d_swtch_dx) * r2inv;                                                                  
#endif                                                       
                            PMEFloat dfdx               = df * xij;
                            PMEFloat dfdy               = df * yij;
                            PMEFloat dfdz               = df * zij;
                            fx_i                       += (PMEDouble)dfdx;
                            fy_i                       += (PMEDouble)dfdy;
                            fz_i                       += (PMEDouble)dfdz;
                            psF[j].x                   -= (PMEDouble)dfdx;
                            psF[j].y                   -= (PMEDouble)dfdy;
                            psF[j].z                   -= (PMEDouble)dfdz;
#ifdef PME_VIRIAL                   
                            vir_11                     -= xij * dfdx;
                            vir_22                     -= yij * dfdy;
                            vir_33                     -= zij * dfdz;
#endif                                         
                        }
                        exclusion                     >>= 1;
                        j                               = sNext[j];   
                    }
                    while (j != tgx);
                }
                else
                {
                    do 
                    {
                        PMEFloat xij                    = xi - psA[j].x; 
                        PMEFloat yij                    = yi - psA[j].y; 
                        PMEFloat zij                    = zi - psA[j].z;                        
                        PMEFloat r2                     = xij * xij + yij * yij + zij * zij;                            
          
                        if (r2 < cSim.cut2)
                        {                                       
                            PMEFloat rinv               = rsqrt(r2);
                            PMEFloat sig                = sigi + psA[j].sig;
                            PMEFloat r                  = r2 * rinv;      
                            PMEFloat r2inv              = rinv * rinv;
                            PMEFloat r6inv              = r2inv * r2inv * r2inv;
                            PMEFloat eps                = epsi * psA[j].eps;
                            sig                        *= sig * sig;
                            sig                        *= sig * r6inv;
                            PMEFloat f6                 = eps * sig;
                            PMEFloat f12                = f6 * sig;
                            PMEFloat qiqj               = qi * psA[j].q;                                                
#ifdef use_DPDP
                            PMEFloat swtch              = erfc(cSim.ew_coeff * r);
#else                                     
                            PMEFloat swtch              = fasterfc(cSim.ew_coeff * r);
#endif         
                            PMEFloat d_swtch_dx         = cSim.negTwoEw_coeffRsqrtPI * exp(-cSim.ew_coeff2 * r2);
                            PMEFloat df                 = ((PMEFloat)12.0 * f12 - (PMEFloat)6.0 * f6) * r2inv;
#ifdef PME_ENERGY                            
                            evdw                       += (double)(f12 - f6);
                            PMEFloat b0                 = qiqj * swtch * rinv;
                            PMEFloat b1                 = b0 - qiqj * d_swtch_dx;
                            df                         += b1 * r2inv;                                                  
                            eed                        += (double)b0; 
#else
                            df                         += qiqj * (swtch * rinv - d_swtch_dx) * r2inv;                                                                         
#endif                                                       
                            PMEFloat dfdx               = df * xij;
                            PMEFloat dfdy               = df * yij;
                            PMEFloat dfdz               = df * zij;
                            fx_i                       += (PMEDouble)dfdx;
                            fy_i                       += (PMEDouble)dfdy;
                            fz_i                       += (PMEDouble)dfdz;
                            psF[j].x                   -= (PMEDouble)dfdx;
                            psF[j].y                   -= (PMEDouble)dfdy;
                            psF[j].z                   -= (PMEDouble)dfdz;
#ifdef PME_VIRIAL                   
                            vir_11                     -= xij * dfdx;
                            vir_22                     -= yij * dfdy;
                            vir_33                     -= zij * dfdz;
#endif                                         
                        }
                        j                               = sNext[j];   
                    }
                    while (j != tgx);                
                }
                      
                // Dump shared memory forces
                if (tx + tgx < psWarp->nlEntry.NL.xAtoms[psWarp->nlpos])
                {
                    int offset                          = (psA[tgx].ID >> NLCELLSHIFT) + (psA[tgx].ID & NLCELLTYPEMASK) * cSim.stride3 + psWarp->bufferOffset;       
                    cSim.pForceXBuffer[offset]         += psF[tgx].x;
                    cSim.pForceYBuffer[offset]         += psF[tgx].y;
                    cSim.pForceZBuffer[offset]         += psF[tgx].z;
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
                int offset                              = psWarp->ypos + tgx + psWarp->homeCellBuffer * cSim.stride3;   
                cSim.pForceXBuffer[offset]             += fx_i;
                cSim.pForceYBuffer[offset]             += fy_i;
                cSim.pForceZBuffer[offset]             += fz_i;
            }                
         
             // Advance to next set of register atoms
            psWarp->ypos                               += cSim.NLYStride;
            psWarp->nlpos++;                       
        }        
        
#ifdef PME_VIRIAL
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
            psV->vir_11                                += (long long int)((psF[tgx].x + psF[tgx + 1].x) * ENERGYSCALE + (PMEDouble)0.5); 
            psV->vir_22                                += (long long int)((psF[tgx].y + psF[tgx + 1].y) * ENERGYSCALE + (PMEDouble)0.5); 
            psV->vir_33                                += (long long int)((psF[tgx].z + psF[tgx + 1].z) * ENERGYSCALE + (PMEDouble)0.5); 
        }
#endif

        // Get next Neighbor List entry
        if (tgx == 0)
            psWarp->pos                                 = atomicAdd(cSim.pNLPosition, 1);
    }
              
    // Do virial and energy reductions if necessary
#ifdef PME_VIRIAL
__syncthreads();
#if (__CUDA_ARCH__ >= 200)
    if (threadIdx.x < 8)
    {
#ifdef use_DPDP    
        sWarpVirial[threadIdx.x].vir_11                += sWarpVirial[threadIdx.x + 8].vir_11;
        sWarpVirial[threadIdx.x].vir_22                += sWarpVirial[threadIdx.x + 8].vir_22;
        sWarpVirial[threadIdx.x].vir_33                += sWarpVirial[threadIdx.x + 8].vir_33;   
#else
        sWarpVirial[threadIdx.x].vir_11                += sWarpVirial[threadIdx.x + 8].vir_11 + sWarpVirial[threadIdx.x + 16].vir_11;
        sWarpVirial[threadIdx.x].vir_22                += sWarpVirial[threadIdx.x + 8].vir_22 + sWarpVirial[threadIdx.x + 16].vir_22;
        sWarpVirial[threadIdx.x].vir_33                += sWarpVirial[threadIdx.x + 8].vir_33 + sWarpVirial[threadIdx.x + 16].vir_33;
#endif
    }   
#endif
#if (__CUDA_ARCH__ >= 200) || !defined(use_DPDP)
    if (threadIdx.x < 4)
    {
        sWarpVirial[threadIdx.x].vir_11                += sWarpVirial[threadIdx.x + 4].vir_11;
        sWarpVirial[threadIdx.x].vir_22                += sWarpVirial[threadIdx.x + 4].vir_22;
        sWarpVirial[threadIdx.x].vir_33                += sWarpVirial[threadIdx.x + 4].vir_33;   
    }
#else    
    if (threadIdx.x == 0)
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

#ifdef PME_ENERGY
#ifdef use_SPSP
#define sEED(i)  sEnergy[i].x
#define sEVDW(i) sEnergy[i].y
#else
#define sEED(i)  sForce[i].x
#define sEVDW(i) sForce[i].y
#endif
    sEED(threadIdx.x)                                   = eed;
    sEVDW(threadIdx.x)                                  = evdw;
    __syncthreads();
    unsigned int m                                      = 1;
    while (m < blockDim.x)
    {
        int p                                           = threadIdx.x + m;
        double dx                                       = ((p < blockDim.x) ? sEED(p)  : 0.0);
        double dy                                       = ((p < blockDim.x) ? sEVDW(p) : 0.0);
        __syncthreads();
        sEED(threadIdx.x)                              += dx;
        sEVDW(threadIdx.x)                             += dy;
        __syncthreads();
        m                                              *= 2;
    }
#endif

#ifdef PME_VIRIAL
    if (threadIdx.x == 0)
    {
        unsigned long long int val1                     = (unsigned long long int)abs(sWarpVirial[0].vir_11);
        unsigned long long int val2                     = (unsigned long long int)abs(sWarpVirial[0].vir_22);
        unsigned long long int val3                     = (unsigned long long int)abs(sWarpVirial[0].vir_33);
        if (sWarpVirial[0].vir_11 < 0)
            val1                                        = 0ull - val1;
        atomicAdd(cSim.pVirial_11, val1);
        if (sWarpVirial[0].vir_22 < 0)
            val2                                        = 0ull - val2;
        atomicAdd(cSim.pVirial_22, val2);  
        if (sWarpVirial[0].vir_33 < 0)
            val3                                        = 0ull - val3;
        atomicAdd(cSim.pVirial_33, val3);
    }
#endif
#ifdef PME_ENERGY
    if (threadIdx.x == 0)
    {
        eed                                             = sEED(0);
        evdw                                            = sEVDW(0);
        unsigned long long int val4                     = (unsigned long long int)(fabs(eed) * ENERGYSCALE + (PMEDouble)0.5);
        unsigned long long int val5                     = (unsigned long long int)(fabs(evdw) * ENERGYSCALE + (PMEDouble)0.5);
        if (eed < (PMEDouble)0.0)
            val4                                        = 0ull - val4;
        if (evdw < (PMEDouble)0.0)
            val5                                        = 0ull - val5;        
        atomicAdd(cSim.pEED, val4);
        atomicAdd(cSim.pEVDW, val5);        
    }
#undef sEED
#undef sEVDW    
#endif       
}
