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
    PMEFloat temp7;
};
#if (__CUDA_ARCH__ >= 200)
volatile __shared__ Atom sA[SM_2X_GBNONBONDENERGY2_THREADS_PER_BLOCK];
volatile __shared__ PMEDouble sE[SM_2X_GBNONBONDENERGY2_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_2X_GBNONBONDENERGY2_THREADS_PER_BLOCK / GRID];
#else
#ifdef GB_IGB78
volatile __shared__ Atom sA[SM_13_GBNONBONDENERGY2IGB78_THREADS_PER_BLOCK];
volatile __shared__ PMEDouble sE[SM_13_GBNONBONDENERGY2IGB78_THREADS_PER_BLOCK];
volatile __shared__ unsigned int sPos[SM_13_GBNONBONDENERGY2IGB78_THREADS_PER_BLOCK / GRID];
#else
volatile __shared__ Atom sA[SM_13_GBNONBONDENERGY2_THREADS_PER_BLOCK];
volatile __shared__ PMEDouble sE[SM_13_GBNONBONDENERGY2_THREADS_PER_BLOCK];
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
        PMEDouble fx_i                      = (PMEDouble)0.0;
        PMEDouble fy_i                      = (PMEDouble)0.0;
        PMEDouble fz_i                      = (PMEDouble)0.0;
        unsigned int tbx                    = threadIdx.x - tgx;
        volatile Atom* psA                  = &sA[tbx];
        unsigned int next                   = tbx + sNext[tgx];   
        if (x == y) // Handle diagonals uniquely at 50% efficiency, skipping i == j interactions
        {
            PMEDouble fx_j                  = (PMEDouble)0.0;
            PMEDouble fy_j                  = (PMEDouble)0.0;
            PMEDouble fz_j                  = (PMEDouble)0.0;
            PMEFloat xi                     = xyi.x;
            PMEFloat yi                     = xyi.y;
            sA[threadIdx.x].x               = xi;
            sA[threadIdx.x].y               = yi;
            sA[threadIdx.x].z               = zi;
            ri                             -= cSim.offset;
            PMEFloat ri1i                   = (PMEFloat)1.0 / ri;
            sA[threadIdx.x].r               = ri;
            sA[threadIdx.x].r1i             = ri1i;
            sA[threadIdx.x].s               = si;
            
            for (unsigned int j = sNext[tgx]; j != tgx; j = sNext[j])
            {
                PMEFloat xij                = xi - psA[j].x; 
                PMEFloat yij                = yi - psA[j].y; 
                PMEFloat zij                = zi - psA[j].z;
                PMEFloat r2                 = xij * xij + yij * yij + zij * zij;
                PMEFloat dij                = sqrt(r2);
                
                PMEFloat sj                 = psA[j].s;
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
                    if (dij < ri + psA[j].r + cSim.gb_neckcut)
                    {
                        unsigned int ii     = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        unsigned int jj     = round((psA[j].r - cSim.gb_neckoffset) * (PMEFloat)20.0);
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
                    PMEDouble f_x           = (PMEDouble)(xij * datmp);
                    fx_j                   += f_x;
                    fx_i                   -= f_x;
                    PMEDouble f_y           = (PMEDouble)(yij * datmp);
                    fy_j                   += f_y;
                    fy_i                   -= f_y;
                    PMEDouble f_z           = (PMEDouble)(zij * datmp);
                    fz_j                   += f_z;
                    fz_i                   -= f_z;        
                }
                
                // Shuffle forces to next thread
                sE[threadIdx.x]             = fx_j;
                fx_j                        = sE[next];    
                sE[threadIdx.x]             = fy_j;
                fy_j                        = sE[next];
                sE[threadIdx.x]             = fz_j;
                fz_j                        = sE[next];                 
            }
            
            fx_i                           += fx_j;
            fy_i                           += fy_j;
            fz_i                           += fz_j;
            int offset                      = x + tgx + (cSim.nonbondForceBuffers + (x >> GRIDBITS)) * cSim.stride3;
            cSim.pForceXBuffer[offset]      = fx_i;
            cSim.pForceYBuffer[offset]      = fy_i;
            cSim.pForceZBuffer[offset]      = fz_i;
        }
        else
        {
            unsigned int j                  = y + tgx;
            PMEFloat2 xyj                   = cSim.pAtomXYSP[j];
            sA[threadIdx.x].z               = cSim.pAtomZSP[j];
            sA[threadIdx.x].r               = cSim.pAtomRBorn[j];
            sA[threadIdx.x].s               = cSim.pAtomS[j];
            sA[threadIdx.x].temp7           = cSim.pTemp7[j];
            PMEDouble fx_j                  = (PMEDouble)0.0;
            PMEDouble fy_j                  = (PMEDouble)0.0;
            PMEDouble fz_j                  = (PMEDouble)0.0;
            PMEFloat xi                     = xyi.x;
            PMEFloat yi                     = xyi.y;
            ri                             -= cSim.offset;
            PMEFloat ri1i                   = (PMEFloat)1.0 / ri;
            PMEFloat si2                    = si * si;     
            sA[threadIdx.x].x               = xyj.x;
            sA[threadIdx.x].y               = xyj.y;
            sA[threadIdx.x].r              -= cSim.offset;
            sA[threadIdx.x].r1i             = (PMEFloat)1.0 / sA[threadIdx.x].r;
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
                PMEFloat dij3i              = dij2i * dij1i;
                
                // Atom i forces
                PMEFloat sj                 = psA[j].s;
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
                    if (dij < ri + psA[j].r + cSim.gb_neckcut)
                    {
                        unsigned int ii     = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        unsigned int jj     = round((psA[j].r - cSim.gb_neckoffset) * (PMEFloat)20.0);
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
                    PMEDouble f_x           = (PMEDouble)(xij * datmp);
                    fx_j                   += f_x;
                    fx_i                   -= f_x;
                    PMEDouble f_y           = (PMEDouble)(yij * datmp);
                    fy_j                   += f_y;
                    fy_i                   -= f_y;
                    PMEDouble f_z           = (PMEDouble)(zij * datmp);
                    fz_j                   += f_z;            
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
                    else if (dij > psA[j].r + si)
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (r2 - si2);
                        PMEFloat v4         = log(v3 * (dij - si));
                        datmp               = v2 * si * ((PMEFloat)-0.5 * dij2i + v2) + 
                                              (PMEFloat)0.25 * dij3i * v4;
                    }            
                    else if (dij > abs(psA[j].r - si))
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (dij + si);
                        PMEFloat v4         = log(v3 * psA[j].r);
                        datmp               = (PMEFloat)-0.25 * ((PMEFloat)-0.5 * (r2 - psA[j].r * psA[j].r + si2) *
                                              dij3i * psA[j].r1i * psA[j].r1i + dij1i * v2 * 
                                              (v2 - dij1i) - dij3i * v4);
                    }    
                    else if (psA[j].r < si)
                    {
                        PMEFloat v2         = (PMEFloat)1.0 / (r2 - si2);
                        PMEFloat v4         = log(v3 * (si - dij));
                        datmp               = (PMEFloat)-0.5 * (si * dij2i * v2 -
                                               (PMEFloat)2.0 * si * v2 * v2 - 
                                               (PMEFloat)0.5 * dij3i * v4);
                    }
 
#ifdef GB_IGB78                    
                    if (dij < ri + psA[j].r + cSim.gb_neckcut)
                    {
                        unsigned int ii     = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
                        unsigned int jj     = round((psA[j].r - cSim.gb_neckoffset) * (PMEFloat)20.0);
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
 
                    datmp                  *= psA[j].temp7;
                    PMEDouble f_x           = (PMEDouble)(xij * datmp);
                    fx_i                   += f_x;
                    fx_j                   -= f_x;
                    PMEDouble f_y           = (PMEDouble)(yij * datmp);
                    fy_i                   += f_y;           
                    fy_j                   -= f_y;         
                    PMEDouble f_z           = (PMEDouble)(zij * datmp);
                    fz_i                   += f_z;
                    fz_j                   -= f_z;        
                }
                             
                // Shuffle forces to next thread
                sE[threadIdx.x]             = fx_j;
                fx_j                        = sE[next];    
                sE[threadIdx.x]             = fy_j;
                fy_j                        = sE[next];
                sE[threadIdx.x]             = fz_j;
                fz_j                        = sE[next];               
                j                           = sNext[j];
            }
            while (j != tgx);
            
            // Write forces
            int offset                      = x + tgx + (cSim.nonbondForceBuffers + (y >> GRIDBITS)) * cSim.stride3;
            cSim.pForceXBuffer[offset]      = fx_i;
            cSim.pForceYBuffer[offset]      = fy_i;
            cSim.pForceZBuffer[offset]      = fz_i;
            offset                          = y + tgx + (cSim.nonbondForceBuffers + (x >> GRIDBITS)) * cSim.stride3;
            cSim.pForceXBuffer[offset]      = fx_j;
            cSim.pForceYBuffer[offset]      = fy_j;
            cSim.pForceZBuffer[offset]      = fz_j;       
        }           
        if (tgx == 0)
            *psPos                          = atomicAdd(cSim.pGBNB2Position, 1);
    }
}

