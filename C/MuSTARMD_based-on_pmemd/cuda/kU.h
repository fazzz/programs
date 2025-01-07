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
#ifdef MPI
#if (__CUDA_ARCH__ >= 200)
__shared__ volatile PMEDouble sForce[3 * SM_2X_UPDATE_THREADS_PER_BLOCK];
#else
__shared__ volatile PMEDouble sForce[3 * SM_13_UPDATE_THREADS_PER_BLOCK];
#endif
    unsigned int tgx                                = threadIdx.x & GRIDBITSMASK;;
    unsigned int tbx                                = threadIdx.x - tgx;
    volatile PMEDouble* psForce                     = &sForce[3 * tbx];
#endif

    unsigned int increment                          = gridDim.x * blockDim.x;
    PMEDouble dtx                                   = dt * (PMEDouble)20.455;
#ifdef LANGEVIN   
    PMEDouble gammai                                = gamma_ln / (PMEDouble)20.455;
    PMEDouble half_dtx                              = dtx * (PMEDouble)0.5;
    PMEDouble c_implic                              = (PMEDouble)1.0 / ((PMEDouble)1.0 + gammai * half_dtx);
    PMEDouble c_explic                              = (PMEDouble)1.0 - gammai * half_dtx;
    PMEDouble sdfac                                 = (PMEDouble)4.0 * gammai * boltz2 * temp0 / dtx;
    unsigned int rpos                               = cSim.pRandomPos[blockIdx.x];    
#endif
	unsigned int pos                                = blockIdx.x * blockDim.x + threadIdx.x;
#ifdef MPI
    while (pos < cSim.paddedNumberOfAtoms)
#else
    while (pos < cSim.atoms)
#endif
    {
#ifdef PME
        PMEDouble atomX                             = cSim.pImageX[pos];
        PMEDouble atomY                             = cSim.pImageY[pos];
        PMEDouble atomZ                             = cSim.pImageZ[pos];
#else
        PMEDouble atomX                             = cSim.pAtomX[pos];
        PMEDouble atomY                             = cSim.pAtomY[pos];
        PMEDouble atomZ                             = cSim.pAtomZ[pos];
#endif            
#ifdef PME
        PMEDouble invMass                           = cSim.pImageInvMass[pos];
#ifdef LANGEVIN
        PMEDouble aamass                            = cSim.pImageMass[pos];
#endif
#ifdef MPI      
        PMEDouble* pInForce                         = &cSim.pInForce[3 * (pos - tgx) + tgx];
        psForce[tgx]                                = *pInForce;
        pInForce                                   += GRID;
        psForce[tgx + GRID]                         = *pInForce;
        pInForce                                   += GRID;
        psForce[tgx + GRID * 2]                     = *pInForce;
        
        volatile PMEDouble* psForce1                = &psForce[3 * tgx];
        PMEDouble forceX                            = *psForce1++;
        PMEDouble forceY                            = *psForce1++;
        PMEDouble forceZ                            = *psForce1;
#else
        PMEDouble forceX                            = cSim.pForceX[pos];
        PMEDouble forceY                            = cSim.pForceY[pos];
        PMEDouble forceZ                            = cSim.pForceZ[pos];
#endif
        PMEDouble velX                              = cSim.pImageVelX[pos];
        PMEDouble velY                              = cSim.pImageVelY[pos];
        PMEDouble velZ                              = cSim.pImageVelZ[pos];
#else
        PMEDouble invMass                           = cSim.pAtomInvMass[pos];
#ifdef LANGEVIN
        PMEDouble aamass                            = cSim.pAtomMass[pos];
#endif
#ifdef MPI
        PMEDouble* pInForce                         = &cSim.pInForce[3 * (pos - tgx) + tgx];
        psForce[tgx]                                = *pInForce;
        pInForce                                   += GRID;
        psForce[tgx + GRID]                         = *pInForce;
        pInForce                                   += GRID;
        psForce[tgx + GRID * 2]                     = *pInForce;
        
        volatile PMEDouble* psForce1                = &psForce[3 * tgx];
        PMEDouble forceX                            = *psForce1++;
        PMEDouble forceY                            = *psForce1++;
        PMEDouble forceZ                            = *psForce1;
#else
        PMEDouble forceX                            = cSim.pForceX[pos];
        PMEDouble forceY                            = cSim.pForceY[pos];
        PMEDouble forceZ                            = cSim.pForceZ[pos]; 
#endif      
        PMEDouble velX                              = cSim.pVelX[pos];
        PMEDouble velY                              = cSim.pVelY[pos];
        PMEDouble velZ                              = cSim.pVelZ[pos];
#endif                
#ifdef LANGEVIN
        PMEDouble gaussX                            = cSim.pRandomX[rpos + pos];
        PMEDouble gaussY                            = cSim.pRandomY[rpos + pos];
        PMEDouble gaussZ                            = cSim.pRandomZ[rpos + pos]; 
        PMEDouble rsd                               = sqrt(sdfac * aamass);
#endif
        PMEDouble wfac                              = invMass * dtx;                               
            
        // Save previous velocities
#ifdef PME                
        cSim.pImageLVelX[pos]                       = velX;
        cSim.pImageLVelY[pos]                       = velY;
        cSim.pImageLVelZ[pos]                       = velZ;
#else     
        cSim.pLVelX[pos]                            = velX;
        cSim.pLVelY[pos]                            = velY;
        cSim.pLVelZ[pos]                            = velZ;
#endif  
            
        // Update velocities          
#ifdef LANGEVIN
        velX                                        = (velX * c_explic + (forceX + rsd * gaussX) * wfac) * c_implic;
        velY                                        = (velY * c_explic + (forceY + rsd * gaussY) * wfac) * c_implic;
        velZ                                        = (velZ * c_explic + (forceZ + rsd * gaussZ) * wfac) * c_implic;
#else                
        velX                                       += forceX * wfac;
        velY                                       += forceY * wfac;
        velZ                                       += forceZ * wfac;
#endif            
        // Test for vlimit if active
        // Vlimit not supported on GPU 
            
        // Save new velocity
#ifdef PME
        cSim.pImageVelX[pos]                        = velX;
        cSim.pImageVelY[pos]                        = velY;
        cSim.pImageVelZ[pos]                        = velZ;
#else                
        cSim.pVelX[pos]                             = velX;
        cSim.pVelY[pos]                             = velY;
        cSim.pVelZ[pos]                             = velZ;
#endif
            
        // Update positions and save old position for SHAKE and kinetic energy kernel      
        cSim.pForceX[pos]                           = atomX;
        cSim.pForceY[pos]                           = atomY;
        cSim.pForceZ[pos]                           = atomZ;                
        PMEDouble newAtomX                          = atomX + velX * dtx;
        PMEDouble newAtomY                          = atomY + velY * dtx;
        PMEDouble newAtomZ                          = atomZ + velZ * dtx;
#ifdef PME
        cSim.pImageX[pos]                           = newAtomX;
        cSim.pImageY[pos]                           = newAtomY;
        cSim.pImageZ[pos]                           = newAtomZ;
#else                
        PMEFloat2 xy;
        xy.x                                        = newAtomX;
        xy.y                                        = newAtomY;
        cSim.pAtomX[pos]                            = newAtomX;
        cSim.pAtomY[pos]                            = newAtomY;
        cSim.pAtomZ[pos]                            = newAtomZ;
        cSim.pAtomXYSP[pos]                         = xy;
        cSim.pAtomZSP[pos]                          = newAtomZ;
#endif     
        pos                                        += increment;
    }

#ifdef LANGEVIN    
    // Update random number position
    __syncthreads();
    if (threadIdx.x == 0)
        cSim.pRandomPos[blockIdx.x]                 = rpos + cSim.paddedNumberOfAtoms;
#endif            
}
