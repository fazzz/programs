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
// CLEAR_NTP, CLEAR_LARGE, CLEAR_YDIVISOR. CLEAR_XDIVISOR, MPI
{
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x; 
    if (pos == 0)
        *(cSim.pNLPosition)                             = cSim.NLNonbondWarps;    
    if (pos < ENERGYTERMS)
        cSim.pEnergyBuffer[pos]                         = 0;      
#ifdef MPI
    while (pos < cSim.reducedAtoms3)
    {
        unsigned int bpos                               = pos;
        unsigned int offset                             = 0;
        if (bpos >= cSim.reducedAtoms)
        {
           bpos                                        -= cSim.reducedAtoms;
           offset                                      += cSim.stride;
        }     
        if (bpos >= cSim.reducedAtoms)
        {
            bpos                                       -= cSim.reducedAtoms;   
            offset                                     += cSim.stride;
        }
        bpos                                           += cSim.minReducedAtom;
        if ((bpos >= cSim.atoms) && (cSim.processedAtoms != cSim.atoms))
            bpos                                       -= cSim.atoms;
        unsigned int pos1                               = bpos + offset;        
#else
    while (pos < cSim.stride3)
    {   
        unsigned int pos1                               = pos;
#if defined(CLEAR_LARGE) || !defined(CLEAR_NTP)        
        unsigned int bpos                               = pos;    
        if (bpos >= cSim.stride)
            bpos                                       -= cSim.stride;
        if (bpos >= cSim.stride)
            bpos                                       -= cSim.stride;   
#endif     
#endif     

#if defined(CLEAR_LARGE) || !defined(CLEAR_NTP)
        unsigned int maxBuffers                         = cSim.pImageOutputBuffers[bpos]; 
#endif

#ifdef CLEAR_NTP     
        // Clear local force accumulator
        cSim.pBondedForceAccumulator[pos1]              = 0ull;
#endif
        
        // Unroll initial 14 + 1
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;  
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;          
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        

#if (CLEAR_YDIVISOR >= 2)
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;  
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;          
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
#endif

#if (CLEAR_YDIVISOR >= 3)
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;  
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;          
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
#endif

#if (CLEAR_YDIVISOR >= 4)
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;  
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;          
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
#endif

#if (CLEAR_YDIVISOR >= 8)
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;  
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;          
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;  
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;          
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;  
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;          
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;  
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;          
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;   
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3; 
#endif

#if (CLEAR_XDIVISOR >= 2)
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;  
#endif

#if (CLEAR_XDIVISOR >= 5)
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3; 
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;  
#endif

#if (CLEAR_XDIVISOR >= 7)
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3; 
#endif

#if (CLEAR_XDIVISOR >= 14)
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3; 
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;  
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3; 
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;  
        cSim.pForceBuffer[pos1]                         = (PMEDouble)0.0;
        pos1                                           += cSim.stride3;                  
#endif

#if defined(CLEAR_LARGE) || !defined(CLEAR_NTP)
        while (pos1 < maxBuffers)
        {
            cSim.pForceBuffer[pos1]                     = (PMEDouble)0.0;
            pos1                                       += cSim.stride3;
        }       
#endif   
     
        pos                                            += blockDim.x * gridDim.x;
    }
}   

