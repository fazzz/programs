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
// REDUCE_NTP, REDUCE_LARGE, REDUCE_XDIVISOR, REDUCE_YDIVISOR, MPI, MPI_REDUCE_NODE0
{
    unsigned int pos                                    = blockIdx.x * blockDim.x + threadIdx.x; 
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
        if (bpos >= cSim.atoms)
            bpos                                       -= cSim.atoms;
        unsigned int pos1                               = bpos + offset;
        unsigned int apos                               = pos1;          
#else
    while (pos < cSim.stride3)
    {
 
        unsigned int pos1                               = pos;
        unsigned int apos                               = pos;

#if defined(REDUCE_LARGE) || !defined(REDUCE_NTP)
        unsigned int bpos                               = pos;   
        if (bpos >= cSim.stride)
            bpos                                       -= cSim.stride;
        if (bpos >= cSim.stride)
            bpos                                       -= cSim.stride;   
#endif
#endif 
#if defined(REDUCE_LARGE) || !defined(REDUCE_NTP)
        unsigned int maxBuffers                         = cSim.pImageOutputBuffers[bpos];
#endif 
        PMEDouble force                                 = (PMEDouble)0.0;
        
        // Unroll initial 14 + 1
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;          
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;    
        
#if (REDUCE_YDIVISOR >= 2)
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;          
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;                
#endif

#if (REDUCE_YDIVISOR >= 3)
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;          
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;                
#endif

#if (REDUCE_YDIVISOR >= 4)
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;          
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;                
#endif

#if (REDUCE_YDIVISOR >= 8)
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;          
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;                
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;          
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;                
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;          
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;                
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;          
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;        
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;                
#endif


#if (REDUCE_XDIVISOR >= 2)
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;                
#endif

#if (REDUCE_XDIVISOR >= 5)
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;     
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;                      
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;   
#endif

#if (REDUCE_XDIVISOR >= 7)
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;                
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;                        
#endif

#if (REDUCE_XDIVISOR == 14)
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;                
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;              
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;                
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;              
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;                
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;              
        force                                          += cSim.pForceBuffer[pos1];
        pos1                                           += cSim.stride3;                     
#endif

#ifdef REDUCE_NTP        
#ifdef MPI_REDUCE_NODE0
        cSim.pNBForce[apos]                            += force;
#else
        cSim.pNBForce[apos]                             = force;
#endif
         
        // Read local forces
        unsigned long long val                          = cSim.pUllForce[apos];
        if (val >= 0x8000000000000000ull)
        {
            force                                      -= (PMEDouble)(val ^ 0xffffffffffffffffull) / FORCESCALE;
        }
        else
        {
            force                                      += (PMEDouble)val / FORCESCALE;
        }
#endif

#if defined(REDUCE_LARGE) || !defined(REDUCE_NTP)
        // Reduce remainder         
        while (pos1 < maxBuffers)
        {
            force                                      += cSim.pForceBuffer[pos1];
            pos1                                       += cSim.stride3;
        }        
#endif

        cSim.pForce[apos]                               = force;        
        pos                                            += blockDim.x * gridDim.x;
    }
}
