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
   
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
KPMECALCULATECOM_KERNEL
()
{
struct COM1 
{
    double COMX;
    double COMY;
    double COMZ;
};

#if (__CUDA_ARCH__ >= 200)
struct UllCOM1 
{
    PMEUllInt COMX;
    PMEUllInt COMY;
    PMEUllInt COMZ;
};
__shared__ volatile COM1 sC[SM_2X_THREADS_PER_BLOCK];
__shared__ volatile UllCOM1 sCOM[SM_2X_MAXMOLECULES];
#ifdef NTP_LOTSOFMOLECULES
const int maxMolecules                       = SM_2X_MAXMOLECULES;
#endif
#else
__shared__ volatile COM1 sC[SM_13_THREADS_PER_BLOCK];
#endif
    unsigned int pos;
#if (__CUDA_ARCH__ >= 200)     
    // Clear shared memory COM
    pos                                             = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES 
    while (pos < maxMolecules)
#else       
    while (pos < cSim.soluteMolecules)
#endif
    {
        sCOM[pos].COMX                              = (PMEUllInt)0;
        sCOM[pos].COMY                              = (PMEUllInt)0;
        sCOM[pos].COMZ                              = (PMEUllInt)0;
        pos                                        += blockDim.x;
    }
    __syncthreads();
#endif

    // Calculate solute COM
    pos                                             = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int tgx                                = threadIdx.x & (GRID - 1);
    volatile COM1* psC                              = &sC[threadIdx.x - tgx];
    int oldMoleculeID                               = -2;
    int atomID, moleculeID;
    PMEDouble CX                                    = (PMEDouble)0.0;
    PMEDouble CY                                    = (PMEDouble)0.0;
    PMEDouble CZ                                    = (PMEDouble)0.0;

    while (pos < cSim.soluteAtoms)
    {
        atomID                                      = cSim.pImageSoluteAtomID[pos];
        moleculeID                                  = cSim.pSoluteAtomMoleculeID[pos];
        PMEDouble mass                              = (PMEDouble)0.0;
        PMEDouble x, y, z;
        if (atomID != -1)
        {
            mass                                    = cSim.pSoluteAtomMass[pos];
            x                                       = cSim.pImageX[atomID];
            y                                       = cSim.pImageY[atomID];
            z                                       = cSim.pImageZ[atomID];
        }

        
        // Output COM upon changed status
        if (moleculeID != oldMoleculeID)
        {
            psC[tgx].COMX                           = CX;
            psC[tgx].COMY                           = CY;
            psC[tgx].COMZ                           = CZ;   
            if (oldMoleculeID >= 0)      
            {
                if (tgx < 16)
                {
                    psC[tgx].COMX                      += psC[tgx + 16].COMX;
                    psC[tgx].COMY                      += psC[tgx + 16].COMY;
                    psC[tgx].COMZ                      += psC[tgx + 16].COMZ;
                }
                if (tgx < 8)
                {
                    psC[tgx].COMX                      += psC[tgx + 8].COMX;
                    psC[tgx].COMY                      += psC[tgx + 8].COMY;
                    psC[tgx].COMZ                      += psC[tgx + 8].COMZ;
                }
                if (tgx < 4)
                {
                    psC[tgx].COMX                      += psC[tgx + 4].COMX;
                    psC[tgx].COMY                      += psC[tgx + 4].COMY;
                    psC[tgx].COMZ                      += psC[tgx + 4].COMZ;
                }
                if (tgx < 2)
                {
                    psC[tgx].COMX                      += psC[tgx + 2].COMX;
                    psC[tgx].COMY                      += psC[tgx + 2].COMY;
                    psC[tgx].COMZ                      += psC[tgx + 2].COMZ;
                }
                if (tgx == 0)
                {           
                    psC->COMX                          += psC[1].COMX;
                    psC->COMY                          += psC[1].COMY;
                    psC->COMZ                          += psC[1].COMZ;
                    unsigned long long int val1         = (unsigned long long int)(fabs(psC->COMX) * ENERGYSCALE + (PMEDouble)0.5);
                    unsigned long long int val2         = (unsigned long long int)(fabs(psC->COMY) * ENERGYSCALE + (PMEDouble)0.5);
                    unsigned long long int val3         = (unsigned long long int)(fabs(psC->COMZ) * ENERGYSCALE + (PMEDouble)0.5);
                    if (psC->COMX < (PMEDouble)0.0)
                        val1                            = 0ull - val1;
                    if (psC->COMY < (PMEDouble)0.0)
                        val2                            = 0ull - val2;
                    if (psC->COMZ < (PMEDouble)0.0)
                        val3                            = 0ull - val3;
#if (__CUDA_ARCH__ < 200)                    
                    atomicAdd(&cSim.pSoluteUllCOMX[oldMoleculeID], val1);
                    atomicAdd(&cSim.pSoluteUllCOMY[oldMoleculeID], val2);
                    atomicAdd(&cSim.pSoluteUllCOMZ[oldMoleculeID], val3);
#else
#ifdef NTP_LOTSOFMOLECULES
                    if (oldMoleculeID < maxMolecules)
                    {
                        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
                        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
                        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
                    }
                    else
                    {
                        atomicAdd(&cSim.pSoluteUllCOMX[oldMoleculeID], val1);
                        atomicAdd(&cSim.pSoluteUllCOMY[oldMoleculeID], val2);
                        atomicAdd(&cSim.pSoluteUllCOMZ[oldMoleculeID], val3);                   
                    }
#else
                    atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
                    atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
                    atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
#endif
#endif
                }
            }
                
            CX                                      = (PMEDouble)0.0;
            CY                                      = (PMEDouble)0.0;
            CZ                                      = (PMEDouble)0.0;
        }
        
        
        oldMoleculeID                               = moleculeID;     
        if (mass != 0.0)
        {
            CX                                     += mass * x;
            CY                                     += mass * y;
            CZ                                     += mass * z;
        }
        
        pos                                        += blockDim.x * gridDim.x;               
    } 
   
    // Dump last batch of solute data to shared memory
    psC[tgx].COMX                               = CX;
    psC[tgx].COMY                               = CY;
    psC[tgx].COMZ                               = CZ;
    if (oldMoleculeID >= 0)
    {
        if (tgx < 16)
        {
            psC[tgx].COMX                          += psC[tgx + 16].COMX;
            psC[tgx].COMY                          += psC[tgx + 16].COMY;
            psC[tgx].COMZ                          += psC[tgx + 16].COMZ;
        }
        if (tgx < 8)
        {
            psC[tgx].COMX                          += psC[tgx + 8].COMX;
            psC[tgx].COMY                          += psC[tgx + 8].COMY;
            psC[tgx].COMZ                          += psC[tgx + 8].COMZ;
        }
        if (tgx < 4)
        {
            psC[tgx].COMX                          += psC[tgx + 4].COMX;
            psC[tgx].COMY                          += psC[tgx + 4].COMY;
            psC[tgx].COMZ                          += psC[tgx + 4].COMZ;
        }
        if (tgx < 2)
        {
            psC[tgx].COMX                          += psC[tgx + 2].COMX;
            psC[tgx].COMY                          += psC[tgx + 2].COMY;
            psC[tgx].COMZ                          += psC[tgx + 2].COMZ;
        }
        if (tgx == 0)
        {           
            psC->COMX                              += psC[1].COMX;
            psC->COMY                              += psC[1].COMY;
            psC->COMZ                              += psC[1].COMZ;
            unsigned long long int val1             = (unsigned long long int)(fabs(psC->COMX) * ENERGYSCALE + (PMEDouble)0.5);
            unsigned long long int val2             = (unsigned long long int)(fabs(psC->COMY) * ENERGYSCALE + (PMEDouble)0.5);
            unsigned long long int val3             = (unsigned long long int)(fabs(psC->COMZ) * ENERGYSCALE + (PMEDouble)0.5);
            if (psC->COMX < (PMEDouble)0.0)
                val1                                = 0ull - val1;
            if (psC->COMY < (PMEDouble)0.0)
                val2                                = 0ull - val2;
            if (psC->COMZ < (PMEDouble)0.0)
                val3                                = 0ull - val3;
#if (__CUDA_ARCH__ < 200)                    
            atomicAdd(&cSim.pSoluteUllCOMX[oldMoleculeID], val1);
            atomicAdd(&cSim.pSoluteUllCOMY[oldMoleculeID], val2);
            atomicAdd(&cSim.pSoluteUllCOMZ[oldMoleculeID], val3);
#else
#ifdef NTP_LOTSOFMOLECULES
            if (oldMoleculeID < maxMolecules)
            {
                atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
                atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
                atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
            }
            else
            {
                atomicAdd(&cSim.pSoluteUllCOMX[oldMoleculeID], val1);
                atomicAdd(&cSim.pSoluteUllCOMY[oldMoleculeID], val2);
                atomicAdd(&cSim.pSoluteUllCOMZ[oldMoleculeID], val3);                   
            }
#else
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
#endif
#endif
        }

    }
#if (__CUDA_ARCH__ >= 200)
    // Dump solute data to main memory 
    __syncthreads();
    pos                                             = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
    while (pos < maxMolecules)
#else    
    while (pos < cSim.soluteMolecules)
#endif
    {   
        if (sCOM[pos].COMX != 0)
            atomicAdd(&cSim.pSoluteUllCOMX[pos], sCOM[pos].COMX);
        if (sCOM[pos].COMY != 0)
            atomicAdd(&cSim.pSoluteUllCOMY[pos], sCOM[pos].COMY);
        if (sCOM[pos].COMZ != 0)
            atomicAdd(&cSim.pSoluteUllCOMZ[pos], sCOM[pos].COMZ);
        pos                                        += blockDim.x;
    }    
#endif

    // Solvent atoms
    pos                                             = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.solventMolecules)
    {
        int4 atomID                                 = cSim.pImageSolventAtomID[pos];          
        PMEDouble invMass                           = cSim.pSolventInvMass[pos];
        PMEDouble mass                              = cSim.pSolventAtomMass1[pos];
        PMEDouble x                                 = cSim.pImageX[atomID.x];
        PMEDouble y                                 = cSim.pImageY[atomID.x];
        PMEDouble z                                 = cSim.pImageZ[atomID.x];
        PMEDouble CX                                = mass * x;
        PMEDouble CY                                = mass * y;
        PMEDouble CZ                                = mass * z;
        if (atomID.y != -1)
        {
            PMEDouble mass                          = cSim.pSolventAtomMass2[pos];
            PMEDouble x                             = cSim.pImageX[atomID.y];
            PMEDouble y                             = cSim.pImageY[atomID.y];
            PMEDouble z                             = cSim.pImageZ[atomID.y];
            CX                                     += mass * x;
            CY                                     += mass * y;
            CZ                                     += mass * z; 
        }
        if (atomID.z != -1)
        {
            PMEDouble mass                          = cSim.pSolventAtomMass3[pos];
            PMEDouble x                             = cSim.pImageX[atomID.z];
            PMEDouble y                             = cSim.pImageY[atomID.z];
            PMEDouble z                             = cSim.pImageZ[atomID.z];
            CX                                     += mass * x;
            CY                                     += mass * y;
            CZ                                     += mass * z; 
        }
        if (atomID.w != -1)
        {
            PMEDouble mass                          = cSim.pSolventAtomMass4[pos];
            PMEDouble x                             = cSim.pImageX[atomID.w];
            PMEDouble y                             = cSim.pImageY[atomID.w];
            PMEDouble z                             = cSim.pImageZ[atomID.w];
            CX                                     += mass * x;
            CY                                     += mass * y;
            CZ                                     += mass * z; 
        }        
            
        // Sum up results
        cSim.pSolventCOMX[pos]                      = invMass * CX;
        cSim.pSolventCOMY[pos]                      = invMass * CY;
        cSim.pSolventCOMZ[pos]                      = invMass * CZ;
        pos                                        += blockDim.x * gridDim.x;
    }
}

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
KPMECALCULATESOLUTECOM_KERNEL
()
{
struct COM1 
{
    double COMX;
    double COMY;
    double COMZ;
};

#if (__CUDA_ARCH__ >= 200)
struct UllCOM1 
{
    PMEUllInt COMX;
    PMEUllInt COMY;
    PMEUllInt COMZ;
};
__shared__ volatile COM1 sC[SM_2X_THREADS_PER_BLOCK];
__shared__ volatile UllCOM1 sCOM[SM_2X_MAXMOLECULES];
#ifdef NTP_LOTSOFMOLECULES
const int maxMolecules                       = SM_2X_MAXMOLECULES;
#endif
#else
__shared__ volatile COM1 sC[SM_13_THREADS_PER_BLOCK];
#endif
    unsigned int pos;
#if (__CUDA_ARCH__ >= 200)     
    // Clear shared memory COM
    pos                                             = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
    while (pos < maxMolecules)
#else    
    while (pos < cSim.soluteMolecules)
#endif
    {
        sCOM[pos].COMX                              = (PMEUllInt)0;
        sCOM[pos].COMY                              = (PMEUllInt)0;
        sCOM[pos].COMZ                              = (PMEUllInt)0;
        pos                                        += blockDim.x;
    }
    __syncthreads();
#endif

    // Calculate solute COM
    pos                                             = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int tgx                                = threadIdx.x & (GRID - 1);
    volatile COM1* psC                              = &sC[threadIdx.x - tgx];
    int oldMoleculeID                               = -2;
    int atomID, moleculeID;
    PMEDouble mass, x, y, z;
    PMEDouble CX                                    = (PMEDouble)0.0;
    PMEDouble CY                                    = (PMEDouble)0.0;
    PMEDouble CZ                                    = (PMEDouble)0.0;

    while (pos < cSim.soluteAtoms)
    {
        atomID                                      = cSim.pImageSoluteAtomID[pos];
        moleculeID                                  = cSim.pSoluteAtomMoleculeID[pos];
        mass                                        = (PMEDouble)0.0;
        if (atomID != -1)
        {
            mass                                    = cSim.pSoluteAtomMass[pos];
            x                                       = cSim.pImageX[atomID];
            y                                       = cSim.pImageY[atomID];
            z                                       = cSim.pImageZ[atomID];
        }

        
        // Output COM upon changed status
        if (moleculeID != oldMoleculeID)
        {
            psC[tgx].COMX                           = CX;
            psC[tgx].COMY                           = CY;
            psC[tgx].COMZ                           = CZ;    
            if (oldMoleculeID >= 0) 
            {    
                if (tgx < 16)
                {
                    psC[tgx].COMX                      += psC[tgx + 16].COMX;
                    psC[tgx].COMY                      += psC[tgx + 16].COMY;
                    psC[tgx].COMZ                      += psC[tgx + 16].COMZ;
                }
                if (tgx < 8)
                {
                    psC[tgx].COMX                      += psC[tgx + 8].COMX;
                    psC[tgx].COMY                      += psC[tgx + 8].COMY;
                    psC[tgx].COMZ                      += psC[tgx + 8].COMZ;
                }
                if (tgx < 4)
                {
                    psC[tgx].COMX                      += psC[tgx + 4].COMX;
                    psC[tgx].COMY                      += psC[tgx + 4].COMY;
                    psC[tgx].COMZ                      += psC[tgx + 4].COMZ;
                }
                if (tgx < 2)
                {
                    psC[tgx].COMX                      += psC[tgx + 2].COMX;
                    psC[tgx].COMY                      += psC[tgx + 2].COMY;
                    psC[tgx].COMZ                      += psC[tgx + 2].COMZ;
                }
                if (tgx == 0)
                {           
                    psC->COMX                          += psC[1].COMX;
                    psC->COMY                          += psC[1].COMY;
                    psC->COMZ                          += psC[1].COMZ;
                    unsigned long long int val1         = (unsigned long long int)(fabs(psC->COMX) * ENERGYSCALE + (PMEDouble)0.5);
                    unsigned long long int val2         = (unsigned long long int)(fabs(psC->COMY) * ENERGYSCALE + (PMEDouble)0.5);
                    unsigned long long int val3         = (unsigned long long int)(fabs(psC->COMZ) * ENERGYSCALE + (PMEDouble)0.5);
                    if (psC->COMX < (PMEDouble)0.0)
                        val1                            = 0ull - val1;
                    if (psC->COMY < (PMEDouble)0.0)
                        val2                            = 0ull - val2;
                    if (psC->COMZ < (PMEDouble)0.0)
                        val3                            = 0ull - val3;
#if (__CUDA_ARCH__ < 200)                    
                    atomicAdd(&cSim.pSoluteUllCOMX[oldMoleculeID], val1);
                    atomicAdd(&cSim.pSoluteUllCOMY[oldMoleculeID], val2);
                    atomicAdd(&cSim.pSoluteUllCOMZ[oldMoleculeID], val3);
#else
#ifdef NTP_LOTSOFMOLECULES
                    if (oldMoleculeID < maxMolecules)
                    {
                        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
                        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
                        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
                    }
                    else
                    {
                        atomicAdd(&cSim.pSoluteUllCOMX[oldMoleculeID], val1);
                        atomicAdd(&cSim.pSoluteUllCOMY[oldMoleculeID], val2);
                        atomicAdd(&cSim.pSoluteUllCOMZ[oldMoleculeID], val3);                   
                    }
#else
                    atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
                    atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
                    atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
#endif
#endif
                }
            }   
            CX                                      = (PMEDouble)0.0;
            CY                                      = (PMEDouble)0.0;
            CZ                                      = (PMEDouble)0.0;
        }
        
        
        oldMoleculeID                               = moleculeID;     
        if (atomID != -1)
        {
            CX                                     += mass * x;
            CY                                     += mass * y;
            CZ                                     += mass * z;
        }
     
        pos                                        += blockDim.x * gridDim.x;               
    } 
   
    // Dump last batch of solute data to shared memory
    psC[tgx].COMX                                   = CX;
    psC[tgx].COMY                                   = CY;
    psC[tgx].COMZ                                   = CZ;
    if (oldMoleculeID >= 0)
    {
        if (tgx < 16)
        {
            psC[tgx].COMX                          += psC[tgx + 16].COMX;
            psC[tgx].COMY                          += psC[tgx + 16].COMY;
            psC[tgx].COMZ                          += psC[tgx + 16].COMZ;
        }
        if (tgx < 8)
        {
            psC[tgx].COMX                          += psC[tgx + 8].COMX;
            psC[tgx].COMY                          += psC[tgx + 8].COMY;
            psC[tgx].COMZ                          += psC[tgx + 8].COMZ;
        }
        if (tgx < 4)
        {
            psC[tgx].COMX                          += psC[tgx + 4].COMX;
            psC[tgx].COMY                          += psC[tgx + 4].COMY;
            psC[tgx].COMZ                          += psC[tgx + 4].COMZ;
        }
        if (tgx < 2)
        {
            psC[tgx].COMX                          += psC[tgx + 2].COMX;
            psC[tgx].COMY                          += psC[tgx + 2].COMY;
            psC[tgx].COMZ                          += psC[tgx + 2].COMZ;
        }
        if (tgx == 0)
        {           
            psC->COMX                              += psC[1].COMX;
            psC->COMY                              += psC[1].COMY;
            psC->COMZ                              += psC[1].COMZ;
            unsigned long long int val1             = (unsigned long long int)(fabs(psC->COMX) * ENERGYSCALE + (PMEDouble)0.5);
            unsigned long long int val2             = (unsigned long long int)(fabs(psC->COMY) * ENERGYSCALE + (PMEDouble)0.5);
            unsigned long long int val3             = (unsigned long long int)(fabs(psC->COMZ) * ENERGYSCALE + (PMEDouble)0.5);
            if (psC->COMX < (PMEDouble)0.0)
                val1                                = 0ull - val1;
            if (psC->COMY < (PMEDouble)0.0)
                val2                                = 0ull - val2;
            if (psC->COMZ < (PMEDouble)0.0)
                val3                                = 0ull - val3;
#if (__CUDA_ARCH__ < 200)                    
            atomicAdd(&cSim.pSoluteUllCOMX[oldMoleculeID], val1);
            atomicAdd(&cSim.pSoluteUllCOMY[oldMoleculeID], val2);
            atomicAdd(&cSim.pSoluteUllCOMZ[oldMoleculeID], val3);
#else
#ifdef NTP_LOTSOFMOLECULES
            if (oldMoleculeID < maxMolecules)
            {
                atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
                atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
                atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
            }
            else
            {
                atomicAdd(&cSim.pSoluteUllCOMX[oldMoleculeID], val1);
                atomicAdd(&cSim.pSoluteUllCOMY[oldMoleculeID], val2);
                atomicAdd(&cSim.pSoluteUllCOMZ[oldMoleculeID], val3);                   
            }
#else
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
#endif
#endif
        }

    }
#if (__CUDA_ARCH__ >= 200)
    // Dump solute data to main memory 
    __syncthreads();
    pos                                             = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
    while (pos < maxMolecules)
#else    
    while (pos < cSim.soluteMolecules)
#endif
    {   
        if (sCOM[pos].COMX != 0)
            atomicAdd(&cSim.pSoluteUllCOMX[pos], sCOM[pos].COMX);
        if (sCOM[pos].COMY != 0)
            atomicAdd(&cSim.pSoluteUllCOMY[pos], sCOM[pos].COMY);
        if (sCOM[pos].COMZ != 0)
            atomicAdd(&cSim.pSoluteUllCOMZ[pos], sCOM[pos].COMZ);
        pos                                        += blockDim.x;
    }    
#endif
}
   
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
KPMECALCULATECOMKINETICENERGY_KERNEL
()
{
struct COMKineticEnergy 
{
    double EKCOMX;
    double EKCOMY;
    double EKCOMZ;
};

#if (__CUDA_ARCH__ >= 200)
struct UllCOMKineticEnergy 
{
    PMEUllInt EKCOMX;
    PMEUllInt EKCOMY;
    PMEUllInt EKCOMZ;
};
__shared__ volatile COMKineticEnergy sE[SM_2X_THREADS_PER_BLOCK];
__shared__ volatile UllCOMKineticEnergy sEKCOM[SM_2X_MAXMOLECULES];
#ifdef NTP_LOTSOFMOLECULES
const int maxMolecules                       = SM_2X_MAXMOLECULES;
#endif
#else
__shared__ volatile COMKineticEnergy sE[SM_13_THREADS_PER_BLOCK];
#endif
    unsigned int pos;
#if (__CUDA_ARCH__ >= 200)
    // Clear EKCOM
    pos                                             = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
    while (pos < maxMolecules)
#else    
    while (pos < cSim.soluteMolecules)
#endif
    {
        sEKCOM[pos].EKCOMX                          = (PMEUllInt)0;
        sEKCOM[pos].EKCOMY                          = (PMEUllInt)0;
        sEKCOM[pos].EKCOMZ                          = (PMEUllInt)0;
        pos                                        += blockDim.x;
    }
    __syncthreads();
#endif

    // Calculate solute kinetic energy
    pos                                             = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int tgx                                = threadIdx.x & (GRID - 1);
    volatile COMKineticEnergy* psE                  = &sE[threadIdx.x - tgx];
    int oldMoleculeID                               = -2;
    PMEDouble EX                                    = (PMEDouble)0.0;
    PMEDouble EY                                    = (PMEDouble)0.0;
    PMEDouble EZ                                    = (PMEDouble)0.0;

    while (pos < cSim.soluteAtoms)
    {
        int atomID                                  = cSim.pImageSoluteAtomID[pos];
        int moleculeID                              = cSim.pSoluteAtomMoleculeID[pos];
        PMEDouble mass                              = (PMEDouble)0.0;
        PMEDouble vx, vy, vz;
        
        if (atomID != -1)
        {
            mass                                    = cSim.pSoluteAtomMass[pos];
            vx                                      = cSim.pImageVelX[atomID];
            vy                                      = cSim.pImageVelY[atomID];
            vz                                      = cSim.pImageVelZ[atomID];
        }

        
        // Output EKCOM upon changed status
        if ((moleculeID != oldMoleculeID) && (oldMoleculeID != -2))
        {
            psE[tgx].EKCOMX                         = EX;
            psE[tgx].EKCOMY                         = EY;
            psE[tgx].EKCOMZ                         = EZ;
            if (tgx < 16)
            {
                psE[tgx].EKCOMX                    += psE[tgx + 16].EKCOMX;
                psE[tgx].EKCOMY                    += psE[tgx + 16].EKCOMY;
                psE[tgx].EKCOMZ                    += psE[tgx + 16].EKCOMZ;
            }
            if (tgx < 8)
            {
                psE[tgx].EKCOMX                    += psE[tgx + 8].EKCOMX;
                psE[tgx].EKCOMY                    += psE[tgx + 8].EKCOMY;
                psE[tgx].EKCOMZ                    += psE[tgx + 8].EKCOMZ;
            }
            if (tgx < 4)
            {
                psE[tgx].EKCOMX                    += psE[tgx + 4].EKCOMX;
                psE[tgx].EKCOMY                    += psE[tgx + 4].EKCOMY;
                psE[tgx].EKCOMZ                    += psE[tgx + 4].EKCOMZ;
            }
            if (tgx < 2)
            {
                psE[tgx].EKCOMX                    += psE[tgx + 2].EKCOMX;
                psE[tgx].EKCOMY                    += psE[tgx + 2].EKCOMY;
                psE[tgx].EKCOMZ                    += psE[tgx + 2].EKCOMZ;
            }
            if (tgx == 0)
            {           
                psE->EKCOMX                        += psE[1].EKCOMX;
                psE->EKCOMY                        += psE[1].EKCOMY;
                psE->EKCOMZ                        += psE[1].EKCOMZ;
                unsigned long long int val1         = (unsigned long long int)(fabs(psE->EKCOMX) * ENERGYSCALE + (PMEDouble)0.5);
                unsigned long long int val2         = (unsigned long long int)(fabs(psE->EKCOMY) * ENERGYSCALE + (PMEDouble)0.5);
                unsigned long long int val3         = (unsigned long long int)(fabs(psE->EKCOMZ) * ENERGYSCALE + (PMEDouble)0.5);
                if (psE->EKCOMX < (PMEDouble)0.0)
                    val1                            = 0ull - val1;
                if (psE->EKCOMY < (PMEDouble)0.0)
                    val2                            = 0ull - val2;
                if (psE->EKCOMZ < (PMEDouble)0.0)
                    val3                            = 0ull - val3;
#if (__CUDA_ARCH__ < 200)                    
                atomicAdd(&cSim.pSoluteUllEKCOMX[oldMoleculeID], val1);
                atomicAdd(&cSim.pSoluteUllEKCOMY[oldMoleculeID], val2);
                atomicAdd(&cSim.pSoluteUllEKCOMZ[oldMoleculeID], val3);
#else
#ifdef NTP_LOTSOFMOLECULES
                if (oldMoleculeID < maxMolecules)
                {
                    atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMX, val1);
                    atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMY, val2);
                    atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMZ, val3);
                }
                else
                {
                    atomicAdd(&cSim.pSoluteUllEKCOMX[oldMoleculeID], val1);
                    atomicAdd(&cSim.pSoluteUllEKCOMY[oldMoleculeID], val2);
                    atomicAdd(&cSim.pSoluteUllEKCOMZ[oldMoleculeID], val3);
                }
#else    
                atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMX, val1);
                atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMY, val2);
                atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMZ, val3);
#endif
#endif                
            }
                
            EX                                      = (PMEDouble)0.0;
            EY                                      = (PMEDouble)0.0;
            EZ                                      = (PMEDouble)0.0;
        }
        
        
        oldMoleculeID                               = moleculeID;  
        if (atomID != -1)
        {   
            EX                                     += mass * vx;
            EY                                     += mass * vy;
            EZ                                     += mass * vz;
        }
        pos                                        += blockDim.x * gridDim.x;               
    }  
   
    // Dump last batch of solute data to shared memory
    psE[tgx].EKCOMX                                 = EX;
    psE[tgx].EKCOMY                                 = EY;
    psE[tgx].EKCOMZ                                 = EZ;       
    if (tgx < 16)
    {
        psE[tgx].EKCOMX                            += psE[tgx + 16].EKCOMX;
        psE[tgx].EKCOMY                            += psE[tgx + 16].EKCOMY;
        psE[tgx].EKCOMZ                            += psE[tgx + 16].EKCOMZ;
    }
    if (tgx < 8)
    {
        psE[tgx].EKCOMX                            += psE[tgx + 8].EKCOMX;
        psE[tgx].EKCOMY                            += psE[tgx + 8].EKCOMY;
        psE[tgx].EKCOMZ                            += psE[tgx + 8].EKCOMZ;
    }
    if (tgx < 4)
    {
        psE[tgx].EKCOMX                            += psE[tgx + 4].EKCOMX;
        psE[tgx].EKCOMY                            += psE[tgx + 4].EKCOMY;
        psE[tgx].EKCOMZ                            += psE[tgx + 4].EKCOMZ;
    }
    if (tgx < 2)
    {
        psE[tgx].EKCOMX                            += psE[tgx + 2].EKCOMX;
        psE[tgx].EKCOMY                            += psE[tgx + 2].EKCOMY;
        psE[tgx].EKCOMZ                            += psE[tgx + 2].EKCOMZ;
    }
    if ((tgx == 0) && (oldMoleculeID != -2))
    {           
        psE->EKCOMX                                += psE[1].EKCOMX;
        psE->EKCOMY                                += psE[1].EKCOMY;
        psE->EKCOMZ                                += psE[1].EKCOMZ;
        unsigned long long int val1                 = (unsigned long long int)(fabs(psE->EKCOMX) * ENERGYSCALE + (PMEDouble)0.5);
        unsigned long long int val2                 = (unsigned long long int)(fabs(psE->EKCOMY) * ENERGYSCALE + (PMEDouble)0.5);
        unsigned long long int val3                 = (unsigned long long int)(fabs(psE->EKCOMZ) * ENERGYSCALE + (PMEDouble)0.5);
        if (psE->EKCOMX < (PMEDouble)0.0)
            val1                                    = 0ull - val1;
        if (psE->EKCOMY < (PMEDouble)0.0)
            val2                                    = 0ull - val2;
        if (psE->EKCOMZ < (PMEDouble)0.0)
            val3                                    = 0ull - val3; 
#if (__CUDA_ARCH__ < 200)                            
        atomicAdd(&cSim.pSoluteUllEKCOMX[oldMoleculeID], val1);
        atomicAdd(&cSim.pSoluteUllEKCOMY[oldMoleculeID], val2);
        atomicAdd(&cSim.pSoluteUllEKCOMZ[oldMoleculeID], val3);
#else
#ifdef NTP_LOTSOFMOLECULES
        if (oldMoleculeID < maxMolecules)
        {
            atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMX, val1);
            atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMY, val2);
            atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMZ, val3);
        }
        else
        {
            atomicAdd(&cSim.pSoluteUllEKCOMX[oldMoleculeID], val1);
            atomicAdd(&cSim.pSoluteUllEKCOMY[oldMoleculeID], val2);
            atomicAdd(&cSim.pSoluteUllEKCOMZ[oldMoleculeID], val3);
        }
#else    
        atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMX, val1);
        atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMY, val2);
        atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMZ, val3);
#endif
#endif
    }

#if (__CUDA_ARCH__ >= 200)
    // Dump solute atoms to memory
    __syncthreads();
    pos                                             = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
    while (pos < maxMolecules)
#else    
    while (pos < cSim.soluteMolecules)
#endif
    {   
        if (sEKCOM[pos].EKCOMX != 0)
            atomicAdd(&cSim.pSoluteUllEKCOMX[pos], sEKCOM[pos].EKCOMX);
        if (sEKCOM[pos].EKCOMY != 0)
            atomicAdd(&cSim.pSoluteUllEKCOMY[pos], sEKCOM[pos].EKCOMY);
        if (sEKCOM[pos].EKCOMZ != 0)
            atomicAdd(&cSim.pSoluteUllEKCOMZ[pos], sEKCOM[pos].EKCOMZ);
        pos                                        += blockDim.x;
    }    
#endif

    // Solvent atoms
    pos                                             = blockIdx.x * blockDim.x + threadIdx.x;
    PMEDouble EKCOMX                                = (PMEDouble)0.0;
    PMEDouble EKCOMY                                = (PMEDouble)0.0;
    PMEDouble EKCOMZ                                = (PMEDouble)0.0;
    while (pos < cSim.solventMolecules)
    {
        int4 atomID                                 = cSim.pImageSolventAtomID[pos];          
        PMEDouble invMass                           = cSim.pSolventInvMass[pos];
        PMEDouble mass                              = cSim.pSolventAtomMass1[pos];
        PMEDouble vx                                = cSim.pImageVelX[atomID.x];
        PMEDouble vy                                = cSim.pImageVelY[atomID.x];
        PMEDouble vz                                = cSim.pImageVelZ[atomID.x];
        PMEDouble EX                                = mass * vx;
        PMEDouble EY                                = mass * vy;
        PMEDouble EZ                                = mass * vz;
        if (atomID.y != -1)
        {
            PMEDouble mass                          = cSim.pSolventAtomMass2[pos];
            PMEDouble vx                            = cSim.pImageVelX[atomID.y];
            PMEDouble vy                            = cSim.pImageVelY[atomID.y];
            PMEDouble vz                            = cSim.pImageVelZ[atomID.y];
            EX                                     += mass * vx;
            EY                                     += mass * vy;
            EZ                                     += mass * vz; 
        }
        if (atomID.z != -1)
        {
            PMEDouble mass                          = cSim.pSolventAtomMass3[pos];
            PMEDouble vx                            = cSim.pImageVelX[atomID.z];
            PMEDouble vy                            = cSim.pImageVelY[atomID.z];
            PMEDouble vz                            = cSim.pImageVelZ[atomID.z];
            EX                                     += mass * vx;
            EY                                     += mass * vy;
            EZ                                     += mass * vz; 
        }
        if (atomID.w != -1)
        {
            PMEDouble mass                          = cSim.pSolventAtomMass4[pos];
            PMEDouble vx                            = cSim.pImageVelX[atomID.w];
            PMEDouble vy                            = cSim.pImageVelY[atomID.w];
            PMEDouble vz                            = cSim.pImageVelZ[atomID.w];
            EX                                     += mass * vx;
            EY                                     += mass * vy;
            EZ                                     += mass * vz; 
        }        
            
        // Sum up results
        EKCOMX                                     += invMass * EX * EX;
        EKCOMY                                     += invMass * EY * EY;
        EKCOMZ                                     += invMass * EZ * EZ;
        pos                                        += blockDim.x * gridDim.x;
    }
    
    sE[threadIdx.x].EKCOMX                          = EKCOMX;
    sE[threadIdx.x].EKCOMY                          = EKCOMY;
    sE[threadIdx.x].EKCOMZ                          = EKCOMZ;
    __syncthreads();
    
    unsigned int m                                  = 1;
    while (m < blockDim.x)
    {
        int p                                       = threadIdx.x + m;    
        PMEDouble EX                                = ((p < blockDim.x) ? sE[p].EKCOMX : (PMEDouble)0.0);
        PMEDouble EY                                = ((p < blockDim.x) ? sE[p].EKCOMY : (PMEDouble)0.0);
        PMEDouble EZ                                = ((p < blockDim.x) ? sE[p].EKCOMZ : (PMEDouble)0.0);
        __syncthreads();
        sE[threadIdx.x].EKCOMX                     += EX;
        sE[threadIdx.x].EKCOMY                     += EY;
        sE[threadIdx.x].EKCOMZ                     += EZ;
        __syncthreads();
        m                                          *= 2;
    }    
 
    if (threadIdx.x == 0)
    {
        unsigned long long int val1                         = (unsigned long long int)(fabs(sE[0].EKCOMX) * ENERGYSCALE + (PMEDouble)0.5);
        unsigned long long int val2                         = (unsigned long long int)(fabs(sE[0].EKCOMY) * ENERGYSCALE + (PMEDouble)0.5);
        unsigned long long int val3                         = (unsigned long long int)(fabs(sE[0].EKCOMZ) * ENERGYSCALE + (PMEDouble)0.5);
        if (sE[0].EKCOMX < (PMEDouble)0.0)
            val1                                            = 0ull - val1;
        atomicAdd(cSim.pEKCOMX, val1);
        if (sE[0].EKCOMY < (PMEDouble)0.0)
            val2                                            = 0ull - val2;
        atomicAdd(cSim.pEKCOMY, val2);
        if (sE[0].EKCOMZ < (PMEDouble)0.0)
            val3                                            = 0ull - val3;
        atomicAdd(cSim.pEKCOMZ, val3);
   }
}



__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
KCALCULATEMOLECULARVIRIAL_KERNEL
()
{
struct COM1 
{
    double COMX;
    double COMY;
    double COMZ;
};
struct Virial {
    PMEDouble vir_11;
    PMEDouble vir_22;
    PMEDouble vir_33;
};

#if (__CUDA_ARCH__ >= 200)
__shared__ volatile COM1 sCOM[SM_2X_MAXMOLECULES];
__shared__ volatile Virial sV[SM_2X_THREADS_PER_BLOCK];
#ifdef NTP_LOTSOFMOLECULES
const int maxMolecules                       = SM_2X_MAXMOLECULES;
#endif
#else
__shared__ volatile COM1 sCOM[SM_13_MAXMOLECULES];
__shared__ volatile Virial sV[SM_13_THREADS_PER_BLOCK];
#ifdef NTP_LOTSOFMOLECULES
const int maxMolecules                       = SM_13_MAXMOLECULES;
#endif
#endif
    
    // Read solute COMs into shared memory
    unsigned int pos                                = threadIdx.x;
    while (pos < cSim.soluteMolecules)
    {
#ifdef NTP_LOTSOFMOLECULES
        if (pos < maxMolecules)
        {
            sCOM[pos].COMX                          = cSim.pSoluteCOMX[pos];
            sCOM[pos].COMY                          = cSim.pSoluteCOMY[pos];
            sCOM[pos].COMZ                          = cSim.pSoluteCOMZ[pos];
        }
#else
        sCOM[pos].COMX                              = cSim.pSoluteCOMX[pos];
        sCOM[pos].COMY                              = cSim.pSoluteCOMY[pos];
        sCOM[pos].COMZ                              = cSim.pSoluteCOMZ[pos];
#endif
        if (blockIdx.x == 0)
        {
            cSim.pSoluteUllCOMX[pos]                = 0;
            cSim.pSoluteUllCOMY[pos]                = 0;
            cSim.pSoluteUllCOMZ[pos]                = 0;
        }
        pos                                        += blockDim.x;
    } 
    __syncthreads();
    
    // Calculate per-thread virial
    PMEDouble vir_11                                = 0.0;
    PMEDouble vir_22                                = 0.0;
    PMEDouble vir_33                                = 0.0;
    
    pos                                             = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.soluteAtoms)
    {
        int atomID                                  = cSim.pImageSoluteAtomID[pos];     
        int moleculeID                              = cSim.pSoluteAtomMoleculeID[pos]; 
        if (atomID != -1)
        {
            PMEDouble fx                            = cSim.pNBForceX[atomID];
            PMEDouble x                             = cSim.pImageX[atomID];
            PMEDouble fy                            = cSim.pNBForceY[atomID];
            PMEDouble y                             = cSim.pImageY[atomID];
            PMEDouble fz                            = cSim.pNBForceZ[atomID];
            PMEDouble z                             = cSim.pImageZ[atomID];
#ifdef NTP_LOTSOFMOLECULES
            double cx, cy, cz;
            if (moleculeID < maxMolecules)
            {
                cx                                  = sCOM[moleculeID].COMX;
                cy                                  = sCOM[moleculeID].COMY;
                cz                                  = sCOM[moleculeID].COMZ;
            }
            else
            {
                cx                                  = cSim.pSoluteCOMX[moleculeID];
                cy                                  = cSim.pSoluteCOMY[moleculeID];
                cz                                  = cSim.pSoluteCOMZ[moleculeID];
            }
            vir_11                                 += fx * (x - cx);
            vir_22                                 += fy * (y - cy);
            vir_33                                 += fz * (z - cz);
#else            
            vir_11                                 += fx * (x - sCOM[moleculeID].COMX);
            vir_22                                 += fy * (y - sCOM[moleculeID].COMY);
            vir_33                                 += fz * (z - sCOM[moleculeID].COMZ);
#endif            
        }
        pos                                        += blockDim.x * gridDim.x;
    }
    
    pos                                             = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.solventMolecules)
    {
        int4 atomID                                 = cSim.pImageSolventAtomID[pos];    
        PMEDouble COMX                              = cSim.pSolventCOMX[pos];
        PMEDouble COMY                              = cSim.pSolventCOMY[pos];
        PMEDouble COMZ                              = cSim.pSolventCOMZ[pos];
        PMEDouble fx                                = cSim.pNBForceX[atomID.x];
        PMEDouble x                                 = cSim.pImageX[atomID.x];
        PMEDouble fy                                = cSim.pNBForceY[atomID.x];
        PMEDouble y                                 = cSim.pImageY[atomID.x];
        PMEDouble fz                                = cSim.pNBForceZ[atomID.x];
        PMEDouble z                                 = cSim.pImageZ[atomID.x];
        vir_11                                     += fx * (x - COMX);
        vir_22                                     += fy * (y - COMY);
        vir_33                                     += fz * (z - COMZ);
        if (atomID.y != -1)
        {
            fx                                      = cSim.pNBForceX[atomID.y];
            x                                       = cSim.pImageX[atomID.y];
            fy                                      = cSim.pNBForceY[atomID.y];
            y                                       = cSim.pImageY[atomID.y];
            fz                                      = cSim.pNBForceZ[atomID.y];
            z                                       = cSim.pImageZ[atomID.y];
            vir_11                                 += fx * (x - COMX);
            vir_22                                 += fy * (y - COMY);
            vir_33                                 += fz * (z - COMZ);
        }
        if (atomID.z != -1)
        {
            fx                                      = cSim.pNBForceX[atomID.z];
            x                                       = cSim.pImageX[atomID.z];
            fy                                      = cSim.pNBForceY[atomID.z];
            y                                       = cSim.pImageY[atomID.z];
            fz                                      = cSim.pNBForceZ[atomID.z];
            z                                       = cSim.pImageZ[atomID.z];
            vir_11                                 += fx * (x - COMX);
            vir_22                                 += fy * (y - COMY);
            vir_33                                 += fz * (z - COMZ);
        }
        if (atomID.w != -1)
        {
            fx                                      = cSim.pNBForceX[atomID.w];
            x                                       = cSim.pImageX[atomID.w];
            fy                                      = cSim.pNBForceY[atomID.w];
            y                                       = cSim.pImageY[atomID.w];
            fz                                      = cSim.pNBForceZ[atomID.w];
            z                                       = cSim.pImageZ[atomID.w];
            vir_11                                 += fx * (x - COMX);
            vir_22                                 += fy * (y - COMY);
            vir_33                                 += fz * (z - COMZ);
        }    
        pos                                        += blockDim.x * gridDim.x;
    }
    
    // Reduce virial
    sV[threadIdx.x].vir_11                          = vir_11;
    sV[threadIdx.x].vir_22                          = vir_22;
    sV[threadIdx.x].vir_33                          = vir_33;
    __syncthreads();
    unsigned int m                                  = 1;
    while (m < blockDim.x)
    {
        int p                                       = threadIdx.x + m;    
        PMEDouble vir_11                            = ((p < blockDim.x) ? sV[p].vir_11 : (PMEDouble)0.0);
        PMEDouble vir_22                            = ((p < blockDim.x) ? sV[p].vir_22 : (PMEDouble)0.0);
        PMEDouble vir_33                            = ((p < blockDim.x) ? sV[p].vir_33 : (PMEDouble)0.0);
        __syncthreads();
        sV[threadIdx.x].vir_11                     += vir_11;
        sV[threadIdx.x].vir_22                     += vir_22;
        sV[threadIdx.x].vir_33                     += vir_33;
        __syncthreads();
        m                                          *= 2;
    }    
    if (threadIdx.x == 0)
    {
        unsigned long long int val1                 = (unsigned long long int)(fabs(sV[0].vir_11) * ENERGYSCALE + (PMEDouble)0.5);
        unsigned long long int val2                 = (unsigned long long int)(fabs(sV[0].vir_22) * ENERGYSCALE + (PMEDouble)0.5);
        unsigned long long int val3                 = (unsigned long long int)(fabs(sV[0].vir_33) * ENERGYSCALE + (PMEDouble)0.5);
        if (sV[0].vir_11 < (PMEDouble)0.0)
            val1                                    = 0ull - val1;
        if (val1 != 0)
            atomicAdd(cSim.pVirial_11, val1);
        if (sV[0].vir_22 < (PMEDouble)0.0)
            val2                                    = 0ull - val2;
        if (val2 != 0)
            atomicAdd(cSim.pVirial_22, val2);
        if (sV[0].vir_33 < (PMEDouble)0.0)
            val3                                    = 0ull - val3;
        if (val3 != 0)
            atomicAdd(cSim.pVirial_33, val3);
   }  
}

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
KPRESSURESCALECOORDINATES_KERNEL
()
{
struct COM1 
{
    double COMX;
    double COMY;
    double COMZ;
};

#if (__CUDA_ARCH__ >= 200)
__shared__ volatile COM1 sCOM[SM_2X_MAXPSMOLECULES];
#ifdef NTP_LOTSOFMOLECULES
const int maxMolecules                       = SM_2X_MAXPSMOLECULES;
#endif
#else
__shared__ volatile COM1 sCOM[SM_13_MAXPSMOLECULES];
#ifdef NTP_LOTSOFMOLECULES
const int maxMolecules                       = SM_13_MAXPSMOLECULES;
#endif
#endif
__shared__ volatile PMEDouble sLast_recip[9];
__shared__ volatile PMEDouble sUcell[9];
    // Read transformation matrices
    if (threadIdx.x < 9)
    {
        sLast_recip[threadIdx.x]                    = cSim.pNTPData->last_recip[threadIdx.x];
        sUcell[threadIdx.x]                         = cSim.pNTPData->ucell[threadIdx.x];
    }
    __syncthreads();

    
    // Read solute COMs into shared memory
    unsigned int pos                                = threadIdx.x;
    while (pos < cSim.soluteMolecules)
    {
        // Read solute COM
        PMEDouble invMass                           = cSim.pSoluteInvMass[pos];        
        PMEUllInt ullCOMX                           = cSim.pSoluteUllCOMX[pos];
        PMEUllInt ullCOMY                           = cSim.pSoluteUllCOMY[pos];
        PMEUllInt ullCOMZ                           = cSim.pSoluteUllCOMZ[pos];
        invMass                                    *= ONEOVERENERGYSCALE;
        PMEDouble ox, oy, oz;
        if (ullCOMX >= 0x8000000000000000ull)
            ox                                      = -(PMEDouble)(ullCOMX ^ 0xffffffffffffffffull);
        else   
            ox                                      =  (PMEDouble)ullCOMX;
        if (ullCOMY >= 0x8000000000000000ull)
            oy                                      = -(PMEDouble)(ullCOMY ^ 0xffffffffffffffffull);
        else
            oy                                      =  (PMEDouble)ullCOMY;
        if (ullCOMZ >= 0x8000000000000000ull)
            oz                                      = -(PMEDouble)(ullCOMZ ^ 0xffffffffffffffffull);
        else
            oz                                      =  (PMEDouble)ullCOMZ;
        ox                                         *= invMass;
        oy                                         *= invMass;
        oz                                         *= invMass;
    
        // Calculate COM displacement
        PMEDouble cx                                = sLast_recip[0] * ox + sLast_recip[3] * oy + sLast_recip[6] * oz;
        PMEDouble cy                                =                       sLast_recip[4] * oy + sLast_recip[7] * oz;
        PMEDouble cz                                =                                             sLast_recip[8] * oz;  
        cx                                          = sUcell[0] * cx + sUcell[1] * cy + sUcell[2] * cz;
        cy                                          =                  sUcell[4] * cy + sUcell[5] * cz;
        cz                                          =                                   sUcell[8] * cz;
        if (blockIdx.x == 0)
        {
            cSim.pSoluteCOMX[pos]                   = cx;
            cSim.pSoluteCOMY[pos]                   = cy;
            cSim.pSoluteCOMZ[pos]                   = cz;
        }
        cx                                         -= ox;
        cy                                         -= oy;
        cz                                         -= oz;
#ifdef NTP_LOTSOFMOLECULES
        if (pos < maxMolecules)
        {
            sCOM[pos].COMX                          = cx;
            sCOM[pos].COMY                          = cy;
            sCOM[pos].COMZ                          = cz;        
        }
        else
        {
            cSim.pSoluteDeltaCOMX[pos]              = cx;
            cSim.pSoluteDeltaCOMY[pos]              = cy;
            cSim.pSoluteDeltaCOMZ[pos]              = cz;        
        }
#else        
        sCOM[pos].COMX                              = cx;
        sCOM[pos].COMY                              = cy;
        sCOM[pos].COMZ                              = cz;
#endif
        pos                                        += blockDim.x;
    } 
    __syncthreads();
    
    
    pos                                             = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.soluteAtoms)
    {
        int atomID                                  = cSim.pImageSoluteAtomID[pos];     
        int moleculeID                              = cSim.pSoluteAtomMoleculeID[pos]; 
        if (atomID != -1)
        {
            PMEDouble x                             = cSim.pImageX[atomID];
            PMEDouble y                             = cSim.pImageY[atomID];
            PMEDouble z                             = cSim.pImageZ[atomID];
#ifdef NTP_LOTSOFMOLECULES
            if (moleculeID < maxMolecules)
            {
                x                                  += sCOM[moleculeID].COMX;
                y                                  += sCOM[moleculeID].COMY;
                z                                  += sCOM[moleculeID].COMZ;
            }
            else
            {
                x                                  += cSim.pSoluteDeltaCOMX[moleculeID];
                y                                  += cSim.pSoluteDeltaCOMY[moleculeID];
                z                                  += cSim.pSoluteDeltaCOMZ[moleculeID];            
            }
#else            
            x                                      += sCOM[moleculeID].COMX;
            y                                      += sCOM[moleculeID].COMY;
            z                                      += sCOM[moleculeID].COMZ;
#endif
            cSim.pImageX[atomID]                    = x;
            cSim.pImageY[atomID]                    = y;
            cSim.pImageZ[atomID]                    = z;
        }
        pos                                        += blockDim.x * gridDim.x;
    }
    
    pos                                             = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.solventMolecules)
    {
        int4 atomID                                 = cSim.pImageSolventAtomID[pos];    
        PMEDouble m1                                = cSim.pSolventAtomMass1[pos];
        PMEDouble invMass                           = cSim.pSolventInvMass[pos];
        PMEDouble x1                                = cSim.pImageX[atomID.x];
        PMEDouble y1                                = cSim.pImageY[atomID.x];
        PMEDouble z1                                = cSim.pImageZ[atomID.x];
        PMEDouble ox                                = m1 * x1;
        PMEDouble oy                                = m1 * y1;
        PMEDouble oz                                = m1 * z1;
        PMEDouble x2, y2, z2;
        if (atomID.y != -1)
        {
            PMEDouble m2                            = cSim.pSolventAtomMass2[pos];
            x2                                      = cSim.pImageX[atomID.y];
            y2                                      = cSim.pImageY[atomID.y];
            z2                                      = cSim.pImageZ[atomID.y];
            ox                                     += m2 * x2;
            oy                                     += m2 * y2;
            oz                                     += m2 * z2;
        }
        double x3, y3, z3;
        if (atomID.z != -1)
        {   
            PMEDouble m3                            = cSim.pSolventAtomMass3[pos];
            x3                                      = cSim.pImageX[atomID.z];
            y3                                      = cSim.pImageY[atomID.z];
            z3                                      = cSim.pImageZ[atomID.z];
            ox                                     += m3 * x3;
            oy                                     += m3 * y3;
            oz                                     += m3 * z3;
        }
        double x4, y4, z4;
        if (atomID.w != -1)
        {
            PMEDouble m4                            = cSim.pSolventAtomMass4[pos];
            x4                                      = cSim.pImageX[atomID.w];
            y4                                      = cSim.pImageY[atomID.w];
            z4                                      = cSim.pImageZ[atomID.w];
            ox                                     += m4 * x4;
            oy                                     += m4 * y4;
            oz                                     += m4 * z4;
        }    
        
        // Calculate change in COM
        ox                                         *= invMass;
        oy                                         *= invMass;
        oz                                         *= invMass;
        PMEDouble cx                                = sLast_recip[0] * ox + sLast_recip[3] * oy + sLast_recip[6] * oz;
        PMEDouble cy                                =                       sLast_recip[4] * oy + sLast_recip[7] * oz;
        PMEDouble cz                                =                                             sLast_recip[8] * oz;  
        cx                                          = sUcell[0] * cx + sUcell[1] * cy + sUcell[2] * cz;
        cy                                          =                  sUcell[4] * cy + sUcell[5] * cz;
        cz                                          =                                   sUcell[8] * cz;
        cSim.pSolventCOMX[pos]                      = cx;
        cSim.pSolventCOMY[pos]                      = cy;
        cSim.pSolventCOMZ[pos]                      = cz;
        cx                                         -= ox;
        cy                                         -= oy;
        cz                                         -= oz;
        
        // Shift atomic coordinates
        cSim.pImageX[atomID.x]                      = x1 + cx;
        cSim.pImageY[atomID.x]                      = y1 + cy;
        cSim.pImageZ[atomID.x]                      = z1 + cz;
        if (atomID.y != -1)
        {
            cSim.pImageX[atomID.y]                  = x2 + cx;
            cSim.pImageY[atomID.y]                  = y2 + cy;
            cSim.pImageZ[atomID.y]                  = z2 + cz;
        }
        if (atomID.z != -1)
        {
            cSim.pImageX[atomID.z]                  = x3 + cx;
            cSim.pImageY[atomID.z]                  = y3 + cy;
            cSim.pImageZ[atomID.z]                  = z3 + cz;
        }
        if (atomID.w != -1)
        {
            cSim.pImageX[atomID.w]                  = x4 + cx;
            cSim.pImageY[atomID.w]                  = y4 + cy;
            cSim.pImageZ[atomID.w]                  = z4 + cz;
        }
        
        pos                                        += blockDim.x * gridDim.x;
    }
}

// Constraints kernels
__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
KPMECALCULATECONSTRAINTSCOM_KERNEL
()
{
struct COM1 
{
    double COMX;
    double COMY;
    double COMZ;
};

#if (__CUDA_ARCH__ >= 200)
struct UllCOM1 
{
    PMEUllInt COMX;
    PMEUllInt COMY;
    PMEUllInt COMZ;
};
__shared__ volatile COM1 sC[SM_2X_THREADS_PER_BLOCK];
__shared__ volatile UllCOM1 sCOM[SM_2X_MAXMOLECULES];
#ifdef NTP_LOTSOFMOLECULES
const int maxMolecules                       = SM_2X_MAXMOLECULES;
#endif
#else
__shared__ volatile COM1 sC[SM_13_THREADS_PER_BLOCK];
#endif
    unsigned int pos;
#if (__CUDA_ARCH__ >= 200)     
    // Clear shared memory COM
    pos                                             = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
    while (pos < maxMolecules)
#else    
    while (pos < cSim.constraintSoluteMolecules)
#endif
    {
        sCOM[pos].COMX                              = (PMEUllInt)0;
        sCOM[pos].COMY                              = (PMEUllInt)0;
        sCOM[pos].COMZ                              = (PMEUllInt)0;
        pos                                        += blockDim.x;
    }
    __syncthreads();
#endif

    // Calculate solute COM
    pos                                             = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int tgx                                = threadIdx.x & (GRID - 1);
    volatile COM1* psC                              = &sC[threadIdx.x - tgx];
    int oldMoleculeID                               = -2;
    int moleculeID;
    PMEDouble mass, x, y, z;
    PMEDouble CX                                    = (PMEDouble)0.0;
    PMEDouble CY                                    = (PMEDouble)0.0;
    PMEDouble CZ                                    = (PMEDouble)0.0;

    while (pos < cSim.constraintSoluteAtoms)
    {
        moleculeID                                  = cSim.pConstraintSoluteAtomMoleculeID[pos];
        mass                                        = cSim.pConstraintSoluteAtomMass[pos];
        x                                           = cSim.pConstraintSoluteAtomX[pos];
        y                                           = cSim.pConstraintSoluteAtomY[pos];
        z                                           = cSim.pConstraintSoluteAtomZ[pos];

        // Output COM upon changed status
        if (moleculeID != oldMoleculeID)
        {
            psC[tgx].COMX                           = CX;
            psC[tgx].COMY                           = CY;
            psC[tgx].COMZ                           = CZ;    
      		if (oldMoleculeID >= 0)
            {
                if (tgx < 16)
                {
                    psC[tgx].COMX                      += psC[tgx + 16].COMX;
                    psC[tgx].COMY                      += psC[tgx + 16].COMY;
                    psC[tgx].COMZ                      += psC[tgx + 16].COMZ;
                }
                if (tgx < 8)
                {
                    psC[tgx].COMX                      += psC[tgx + 8].COMX;
                    psC[tgx].COMY                      += psC[tgx + 8].COMY;
                    psC[tgx].COMZ                      += psC[tgx + 8].COMZ;
                }
                if (tgx < 4)
                {
                    psC[tgx].COMX                      += psC[tgx + 4].COMX;
                    psC[tgx].COMY                      += psC[tgx + 4].COMY;
                    psC[tgx].COMZ                      += psC[tgx + 4].COMZ;
                }
                if (tgx < 2)
                {
                    psC[tgx].COMX                      += psC[tgx + 2].COMX;
                    psC[tgx].COMY                      += psC[tgx + 2].COMY;
                    psC[tgx].COMZ                      += psC[tgx + 2].COMZ;
                }
                if (tgx == 0)
                {           
                    psC->COMX                          += psC[1].COMX;
                    psC->COMY                          += psC[1].COMY;
                    psC->COMZ                          += psC[1].COMZ;
                    unsigned long long int val1         = (unsigned long long int)(fabs(psC->COMX) * ENERGYSCALE + (PMEDouble)0.5);
                    unsigned long long int val2         = (unsigned long long int)(fabs(psC->COMY) * ENERGYSCALE + (PMEDouble)0.5);
                    unsigned long long int val3         = (unsigned long long int)(fabs(psC->COMZ) * ENERGYSCALE + (PMEDouble)0.5);
                    if (psC->COMX < (PMEDouble)0.0)
                        val1                            = 0ull - val1;
                    if (psC->COMY < (PMEDouble)0.0)
                        val2                            = 0ull - val2;
                    if (psC->COMZ < (PMEDouble)0.0)
                        val3                            = 0ull - val3;
#if (__CUDA_ARCH__ < 200)                    
                    atomicAdd(&cSim.pConstraintSoluteUllCOMX[oldMoleculeID], val1);
                    atomicAdd(&cSim.pConstraintSoluteUllCOMY[oldMoleculeID], val2);
                    atomicAdd(&cSim.pConstraintSoluteUllCOMZ[oldMoleculeID], val3);
#else
#ifdef NTP_LOTSOFMOLECULES
                    if (oldMoleculeID < maxMolecules)
                    {
                        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
                        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
                        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
                    }
                    else
                    {
                        atomicAdd(&cSim.pConstraintSoluteUllCOMX[oldMoleculeID], val1);
                        atomicAdd(&cSim.pConstraintSoluteUllCOMY[oldMoleculeID], val2);
                        atomicAdd(&cSim.pConstraintSoluteUllCOMZ[oldMoleculeID], val3);                   
                    }
#else
                    atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
                    atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
                    atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
#endif
#endif
                }
            }
                
            CX                                      = (PMEDouble)0.0;
            CY                                      = (PMEDouble)0.0;
            CZ                                      = (PMEDouble)0.0;
        }
        
        
        oldMoleculeID                               = moleculeID;     
        CX                                         += mass * x;
        CY                                         += mass * y;
        CZ                                         += mass * z;
        
        pos                                        += blockDim.x * gridDim.x;               
    } 
   
    // Dump last batch of solute data to shared memory
    psC[tgx].COMX                               = CX;
    psC[tgx].COMY                               = CY;
    psC[tgx].COMZ                               = CZ;
    if (oldMoleculeID >= 0)
    {
        if (tgx < 16)
        {
            psC[tgx].COMX                          += psC[tgx + 16].COMX;
            psC[tgx].COMY                          += psC[tgx + 16].COMY;
            psC[tgx].COMZ                          += psC[tgx + 16].COMZ;
        }
        if (tgx < 8)
        {
            psC[tgx].COMX                          += psC[tgx + 8].COMX;
            psC[tgx].COMY                          += psC[tgx + 8].COMY;
            psC[tgx].COMZ                          += psC[tgx + 8].COMZ;
        }
        if (tgx < 4)
        {
            psC[tgx].COMX                          += psC[tgx + 4].COMX;
            psC[tgx].COMY                          += psC[tgx + 4].COMY;
            psC[tgx].COMZ                          += psC[tgx + 4].COMZ;
        }
        if (tgx < 2)
        {
            psC[tgx].COMX                          += psC[tgx + 2].COMX;
            psC[tgx].COMY                          += psC[tgx + 2].COMY;
            psC[tgx].COMZ                          += psC[tgx + 2].COMZ;
        }
        if (tgx == 0)
        {           
            psC->COMX                              += psC[1].COMX;
            psC->COMY                              += psC[1].COMY;
            psC->COMZ                              += psC[1].COMZ;
            unsigned long long int val1             = (unsigned long long int)(fabs(psC->COMX) * ENERGYSCALE + (PMEDouble)0.5);
            unsigned long long int val2             = (unsigned long long int)(fabs(psC->COMY) * ENERGYSCALE + (PMEDouble)0.5);
            unsigned long long int val3             = (unsigned long long int)(fabs(psC->COMZ) * ENERGYSCALE + (PMEDouble)0.5);
            if (psC->COMX < (PMEDouble)0.0)
                val1                                = 0ull - val1;
            if (psC->COMY < (PMEDouble)0.0)
                val2                                = 0ull - val2;
            if (psC->COMZ < (PMEDouble)0.0)
                val3                                = 0ull - val3;
#if (__CUDA_ARCH__ < 200)                    
            atomicAdd(&cSim.pConstraintSoluteUllCOMX[oldMoleculeID], val1);
            atomicAdd(&cSim.pConstraintSoluteUllCOMY[oldMoleculeID], val2);
            atomicAdd(&cSim.pConstraintSoluteUllCOMZ[oldMoleculeID], val3);
#else
#ifdef NTP_LOTSOFMOLECULES
            if (oldMoleculeID < maxMolecules)
            {
                atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
                atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
                atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
            }
            else
            {
                atomicAdd(&cSim.pConstraintSoluteUllCOMX[oldMoleculeID], val1);
                atomicAdd(&cSim.pConstraintSoluteUllCOMY[oldMoleculeID], val2);
                atomicAdd(&cSim.pConstraintSoluteUllCOMZ[oldMoleculeID], val3);                   
            }
#else
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
#endif
#endif
        }

    }
#if (__CUDA_ARCH__ >= 200)
    // Dump solute data to main memory 
    __syncthreads();
    pos                                             = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
    while (pos < maxMolecules)
#else    
    while (pos < cSim.constraintSoluteMolecules)
#endif
    {   
        if (sCOM[pos].COMX != 0)
            atomicAdd(&cSim.pConstraintSoluteUllCOMX[pos], sCOM[pos].COMX);
        if (sCOM[pos].COMY != 0)
            atomicAdd(&cSim.pConstraintSoluteUllCOMY[pos], sCOM[pos].COMY);
        if (sCOM[pos].COMZ != 0)
            atomicAdd(&cSim.pConstraintSoluteUllCOMZ[pos], sCOM[pos].COMZ);
        pos                                        += blockDim.x;
    }    
#endif

    // Solvent atoms
    pos                                             = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.constraintSolventMolecules)
    {       
        int atoms                                   = cSim.pConstraintSolventAtoms[pos];
        PMEDouble invMass                           = cSim.pConstraintSolventInvMass[pos];
        PMEDouble mass                              = cSim.pConstraintSolventAtomMass1[pos];
        PMEDouble x                                 = cSim.pConstraintSolventAtomX1[pos];
        PMEDouble y                                 = cSim.pConstraintSolventAtomY1[pos];
        PMEDouble z                                 = cSim.pConstraintSolventAtomZ1[pos];
        PMEDouble CX                                = mass * x;
        PMEDouble CY                                = mass * y;
        PMEDouble CZ                                = mass * z;
        if (atoms >= 2)
        {
            PMEDouble mass                          = cSim.pConstraintSolventAtomMass2[pos];
            PMEDouble x                             = cSim.pConstraintSolventAtomX2[pos];
            PMEDouble y                             = cSim.pConstraintSolventAtomY2[pos];
            PMEDouble z                             = cSim.pConstraintSolventAtomZ2[pos];
            CX                                     += mass * x;
            CY                                     += mass * y;
            CZ                                     += mass * z; 
        }
        if (atoms >= 3)
        {
            PMEDouble mass                          = cSim.pConstraintSolventAtomMass3[pos];
            PMEDouble x                             = cSim.pConstraintSolventAtomX3[pos];
            PMEDouble y                             = cSim.pConstraintSolventAtomY3[pos];
            PMEDouble z                             = cSim.pConstraintSolventAtomZ3[pos];
            CX                                     += mass * x;
            CY                                     += mass * y;
            CZ                                     += mass * z; 
        }
        if (atoms == 4)
        {
            PMEDouble mass                          = cSim.pConstraintSolventAtomMass4[pos];
            PMEDouble x                             = cSim.pConstraintSolventAtomX4[pos];
            PMEDouble y                             = cSim.pConstraintSolventAtomY4[pos];
            PMEDouble z                             = cSim.pConstraintSolventAtomZ4[pos];
            CX                                     += mass * x;
            CY                                     += mass * y;
            CZ                                     += mass * z; 
        }        
            
        // Sum up results
        cSim.pConstraintSolventCOMX[pos]            = invMass * CX;
        cSim.pConstraintSolventCOMY[pos]            = invMass * CY;
        cSim.pConstraintSolventCOMZ[pos]            = invMass * CZ;
        pos                                        += blockDim.x * gridDim.x;
    }
}

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
KPMECALCULATESOLUTECONSTRAINTSCOM_KERNEL
()
{
struct COM1 
{
    double COMX;
    double COMY;
    double COMZ;
};

#if (__CUDA_ARCH__ >= 200)
struct UllCOM1 
{
    PMEUllInt COMX;
    PMEUllInt COMY;
    PMEUllInt COMZ;
};
__shared__ volatile COM1 sC[SM_2X_THREADS_PER_BLOCK];
__shared__ volatile UllCOM1 sCOM[SM_2X_MAXMOLECULES];
#ifdef NTP_LOTSOFMOLECULES
const int maxMolecules                       = SM_2X_MAXMOLECULES;
#endif
#else
__shared__ volatile COM1 sC[SM_13_THREADS_PER_BLOCK];
#endif
    unsigned int pos;
#if (__CUDA_ARCH__ >= 200)     
    // Clear shared memory COM
    pos                                             = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
    while (pos < maxMolecules)
#else    
    while (pos < cSim.constraintSoluteMolecules)
#endif
    {
        sCOM[pos].COMX                              = (PMEUllInt)0;
        sCOM[pos].COMY                              = (PMEUllInt)0;
        sCOM[pos].COMZ                              = (PMEUllInt)0;
        pos                                        += blockDim.x;
    }
    __syncthreads();
#endif

    // Calculate solute COM
    pos                                             = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int tgx                                = threadIdx.x & (GRID - 1);
    volatile COM1* psC                              = &sC[threadIdx.x - tgx];
    int oldMoleculeID                               = -2;
    int moleculeID;
    PMEDouble mass, x, y, z;
    PMEDouble CX                                    = (PMEDouble)0.0;
    PMEDouble CY                                    = (PMEDouble)0.0;
    PMEDouble CZ                                    = (PMEDouble)0.0;

    while (pos < cSim.constraintSoluteAtoms)
    {
        moleculeID                                  = cSim.pConstraintSoluteAtomMoleculeID[pos];
        mass                                        = cSim.pConstraintSoluteAtomMass[pos];
        x                                           = cSim.pConstraintSoluteAtomX[pos];
        y                                           = cSim.pConstraintSoluteAtomY[pos];
        z                                           = cSim.pConstraintSoluteAtomZ[pos];
        
        // Output COM upon changed status
        if (moleculeID != oldMoleculeID)
        {     
            psC[tgx].COMX                           = CX;
            psC[tgx].COMY                           = CY;
            psC[tgx].COMZ                           = CZ;            
            if (oldMoleculeID >= 0)
            {
                if (tgx < 16)
                {
                    psC[tgx].COMX                  += psC[tgx + 16].COMX;
                    psC[tgx].COMY                  += psC[tgx + 16].COMY;
                    psC[tgx].COMZ                  += psC[tgx + 16].COMZ;
                }
                if (tgx < 8)
                {
                    psC[tgx].COMX                  += psC[tgx + 8].COMX;
                    psC[tgx].COMY                  += psC[tgx + 8].COMY;
                    psC[tgx].COMZ                  += psC[tgx + 8].COMZ;
                }
                if (tgx < 4)
                {
                    psC[tgx].COMX                  += psC[tgx + 4].COMX;
                    psC[tgx].COMY                  += psC[tgx + 4].COMY;
                    psC[tgx].COMZ                  += psC[tgx + 4].COMZ;
                }
                if (tgx < 2)
                {
                    psC[tgx].COMX                  += psC[tgx + 2].COMX;
                    psC[tgx].COMY                  += psC[tgx + 2].COMY;
                    psC[tgx].COMZ                  += psC[tgx + 2].COMZ;
                }
                if (tgx == 0)
                {           
                    psC->COMX                      += psC[1].COMX;
                    psC->COMY                      += psC[1].COMY;
                    psC->COMZ                      += psC[1].COMZ;
                    unsigned long long int val1     = (unsigned long long int)(fabs(psC->COMX) * ENERGYSCALE + (PMEDouble)0.5);
                    unsigned long long int val2     = (unsigned long long int)(fabs(psC->COMY) * ENERGYSCALE + (PMEDouble)0.5);
                    unsigned long long int val3     = (unsigned long long int)(fabs(psC->COMZ) * ENERGYSCALE + (PMEDouble)0.5);
                    if (psC->COMX < (PMEDouble)0.0)
                        val1                        = 0ull - val1;
                    if (psC->COMY < (PMEDouble)0.0)
                        val2                        = 0ull - val2;
                    if (psC->COMZ < (PMEDouble)0.0)
                        val3                        = 0ull - val3;
#if (__CUDA_ARCH__ < 200)           
         
                    atomicAdd(&cSim.pConstraintSoluteUllCOMX[oldMoleculeID], val1);
                    atomicAdd(&cSim.pConstraintSoluteUllCOMY[oldMoleculeID], val2);
                    atomicAdd(&cSim.pConstraintSoluteUllCOMZ[oldMoleculeID], val3);
#else
#ifdef NTP_LOTSOFMOLECULES
                    if (moleculeID < maxMolecules)
                    {
                        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
                        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
                        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
                    }
                    else
                    {
                        atomicAdd(&cSim.pConstraintSoluteUllCOMX[oldMoleculeID], val1);
                        atomicAdd(&cSim.pConstraintSoluteUllCOMY[oldMoleculeID], val2);
                        atomicAdd(&cSim.pConstraintSoluteUllCOMZ[oldMoleculeID], val3);
                    }
#else
                    atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
                    atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
                    atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
#endif
#endif
                }
            }
                
            CX                                      = (PMEDouble)0.0;
            CY                                      = (PMEDouble)0.0;
            CZ                                      = (PMEDouble)0.0;
        }
        
        
        oldMoleculeID                               = moleculeID;     
        CX                                         += mass * x;
        CY                                         += mass * y;
        CZ                                         += mass * z;
     
        pos                                        += blockDim.x * gridDim.x;               
    } 
   
    // Dump last batch of solute data to shared memory
    psC[tgx].COMX                                   = CX;
    psC[tgx].COMY                                   = CY;
    psC[tgx].COMZ                                   = CZ;
    if (oldMoleculeID >= 0)
    {
        if (tgx < 16)
        {
            psC[tgx].COMX                          += psC[tgx + 16].COMX;
            psC[tgx].COMY                          += psC[tgx + 16].COMY;
            psC[tgx].COMZ                          += psC[tgx + 16].COMZ;
        }
        if (tgx < 8)
        {
            psC[tgx].COMX                          += psC[tgx + 8].COMX;
            psC[tgx].COMY                          += psC[tgx + 8].COMY;
            psC[tgx].COMZ                          += psC[tgx + 8].COMZ;
        }
        if (tgx < 4)
        {
            psC[tgx].COMX                          += psC[tgx + 4].COMX;
            psC[tgx].COMY                          += psC[tgx + 4].COMY;
            psC[tgx].COMZ                          += psC[tgx + 4].COMZ;
        }
        if (tgx < 2)
        {
            psC[tgx].COMX                          += psC[tgx + 2].COMX;
            psC[tgx].COMY                          += psC[tgx + 2].COMY;
            psC[tgx].COMZ                          += psC[tgx + 2].COMZ;
        }
        if (tgx == 0)
        {           
            psC->COMX                              += psC[1].COMX;
            psC->COMY                              += psC[1].COMY;
            psC->COMZ                              += psC[1].COMZ;
            unsigned long long int val1             = (unsigned long long int)(fabs(psC->COMX) * ENERGYSCALE + (PMEDouble)0.5);
            unsigned long long int val2             = (unsigned long long int)(fabs(psC->COMY) * ENERGYSCALE + (PMEDouble)0.5);
            unsigned long long int val3             = (unsigned long long int)(fabs(psC->COMZ) * ENERGYSCALE + (PMEDouble)0.5);
            if (psC->COMX < (PMEDouble)0.0)
                val1                                = 0ull - val1;
            if (psC->COMY < (PMEDouble)0.0)
                val2                                = 0ull - val2;
            if (psC->COMZ < (PMEDouble)0.0)
                val3                                = 0ull - val3;
#if (__CUDA_ARCH__ < 200)                    
            atomicAdd(&cSim.pConstraintSoluteUllCOMX[oldMoleculeID], val1);
            atomicAdd(&cSim.pConstraintSoluteUllCOMY[oldMoleculeID], val2);
            atomicAdd(&cSim.pConstraintSoluteUllCOMZ[oldMoleculeID], val3);
#else
#ifdef NTP_LOTSOFMOLECULES
            if (moleculeID < maxMolecules)
            {
                atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
                atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
                atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
            }
            else
            {
                atomicAdd(&cSim.pConstraintSoluteUllCOMX[oldMoleculeID], val1);
                atomicAdd(&cSim.pConstraintSoluteUllCOMY[oldMoleculeID], val2);
                atomicAdd(&cSim.pConstraintSoluteUllCOMZ[oldMoleculeID], val3);
            }
#else
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
#endif
#endif
        }

    }
#if (__CUDA_ARCH__ >= 200)
    // Dump solute data to main memory 
    __syncthreads();
    pos                                             = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
    while (pos < maxMolecules)
#else    
    while (pos < cSim.constraintSoluteMolecules)
#endif
    {   
        if (sCOM[pos].COMX != 0)
            atomicAdd(&cSim.pConstraintSoluteUllCOMX[pos], sCOM[pos].COMX);
        if (sCOM[pos].COMY != 0)
            atomicAdd(&cSim.pConstraintSoluteUllCOMY[pos], sCOM[pos].COMY);
        if (sCOM[pos].COMZ != 0)
            atomicAdd(&cSim.pConstraintSoluteUllCOMZ[pos], sCOM[pos].COMZ);
        pos                                        += blockDim.x;
    }    
#endif
}

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
KPRESSURESCALECONSTRAINTCOORDINATES_KERNEL
()
{
struct COM1 
{
    double COMX;
    double COMY;
    double COMZ;
};
#if (__CUDA_ARCH__ >= 200)
__shared__ volatile COM1 sCOM[SM_2X_MAXPSMOLECULES];
#ifdef NTP_LOTSOFMOLECULES
const int maxMolecules                       = SM_2X_MAXPSMOLECULES;
#endif
#else
__shared__ volatile COM1 sCOM[SM_13_MAXPSMOLECULES];
#ifdef NTP_LOTSOFMOLECULES
const int maxMolecules                       = SM_13_MAXPSMOLECULES;
#endif
#endif
__shared__ volatile PMEDouble sLast_recip[9];
__shared__ volatile PMEDouble sUcell[9];
    // Read transformation matrices
    if (threadIdx.x < 9)
    {
        sLast_recip[threadIdx.x]                    = cSim.pNTPData->last_recip[threadIdx.x];
        sUcell[threadIdx.x]                         = cSim.pNTPData->ucell[threadIdx.x];
    }
    __syncthreads();

    
    // Read solute COMs into shared memory
    unsigned int pos                                = threadIdx.x;
    while (pos < cSim.constraintSoluteMolecules)
    {
        // Read solute COM    
        PMEDouble ox                                = cSim.pConstraintSoluteCOMX[pos];
        PMEDouble oy                                = cSim.pConstraintSoluteCOMY[pos];
        PMEDouble oz                                = cSim.pConstraintSoluteCOMZ[pos];

        // Calculate COM displacement
        PMEDouble cx                                = sLast_recip[0] * ox + sLast_recip[3] * oy + sLast_recip[6] * oz;
        PMEDouble cy                                =                       sLast_recip[4] * oy + sLast_recip[7] * oz;
        PMEDouble cz                                =                                             sLast_recip[8] * oz;  
        cx                                          = sUcell[0] * cx + sUcell[1] * cy + sUcell[2] * cz;
        cy                                          =                  sUcell[4] * cy + sUcell[5] * cz;
        cz                                          =                                   sUcell[8] * cz;
        cx                                         -= ox;
        cy                                         -= oy;
        cz                                         -= oz;
#ifdef NTP_LOTSOFMOLECULES
        if (pos < maxMolecules)
        {
            sCOM[pos].COMX                          = cx;
            sCOM[pos].COMY                          = cy;
            sCOM[pos].COMZ                          = cz;        
        }
        else
        {
            cSim.pConstraintSoluteDeltaCOMX[pos]    = cx;
            cSim.pConstraintSoluteDeltaCOMY[pos]    = cy;
            cSim.pConstraintSoluteDeltaCOMZ[pos]    = cz;        
        }
#else        
        sCOM[pos].COMX                              = cx;
        sCOM[pos].COMY                              = cy;
        sCOM[pos].COMZ                              = cz;
#endif
        pos                                        += blockDim.x;
    } 
    __syncthreads();
    
    pos                                             = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.constraintSoluteAtoms)
    {  
        int moleculeID                              = cSim.pConstraintSoluteAtomMoleculeID[pos];
        PMEDouble x                                 = cSim.pConstraintSoluteAtomX[pos];
        PMEDouble y                                 = cSim.pConstraintSoluteAtomY[pos];
        PMEDouble z                                 = cSim.pConstraintSoluteAtomZ[pos];
        int constraint                              = cSim.pConstraintSoluteAtomConstraintID[pos];
#ifdef NTP_LOTSOFMOLECULES
        if (moleculeID < maxMolecules)
        {
            x                                      += sCOM[moleculeID].COMX;
            y                                      += sCOM[moleculeID].COMY;
            z                                      += sCOM[moleculeID].COMZ;        
        }
        else
        {
            x                                      += cSim.pConstraintSoluteDeltaCOMX[moleculeID];
            y                                      += cSim.pConstraintSoluteDeltaCOMY[moleculeID];
            z                                      += cSim.pConstraintSoluteDeltaCOMZ[moleculeID];                
        }
#else        
        x                                          += sCOM[moleculeID].COMX;
        y                                          += sCOM[moleculeID].COMY;
        z                                          += sCOM[moleculeID].COMZ;
#endif
        cSim.pConstraintSoluteAtomX[pos]            = x;
        cSim.pConstraintSoluteAtomY[pos]            = y;
        cSim.pConstraintSoluteAtomZ[pos]            = z;
        if (constraint != -1)
        {
            PMEDouble2 constraint1                  = cSim.pConstraint1[constraint];
            PMEDouble2 constraint2                  = {y, z};
            constraint1.y                           = x;
            cSim.pConstraint2[constraint]           = constraint2;
            cSim.pConstraint1[constraint]           = constraint1;                
        }
        pos                                        += blockDim.x * gridDim.x;
    }
   
    pos                                             = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.constraintSolventMolecules)
    {
        int atoms                                   = cSim.pConstraintSolventAtoms[pos];
        PMEDouble invMass                           = cSim.pConstraintSolventInvMass[pos];
        PMEDouble m1                                = cSim.pConstraintSolventAtomMass1[pos];
        PMEDouble x1                                = cSim.pConstraintSolventAtomX1[pos];
        PMEDouble y1                                = cSim.pConstraintSolventAtomY1[pos];
        PMEDouble z1                                = cSim.pConstraintSolventAtomZ1[pos];
        int4 constraint                             = cSim.pConstraintSolventConstraint[pos];          
        PMEDouble ox                                = m1 * x1;
        PMEDouble oy                                = m1 * y1;
        PMEDouble oz                                = m1 * z1;
        PMEDouble x2, y2, z2;
        if (atoms >= 2)
        {
            PMEDouble m2                            = cSim.pConstraintSolventAtomMass2[pos];
            x2                                      = cSim.pConstraintSolventAtomX2[pos];
            y2                                      = cSim.pConstraintSolventAtomY2[pos];
            z2                                      = cSim.pConstraintSolventAtomZ2[pos];
            ox                                     += m2 * x2;
            oy                                     += m2 * y2;
            oz                                     += m2 * z2;
        }
        double x3, y3, z3;
        if (atoms >= 3)
        {   
            PMEDouble m3                            = cSim.pConstraintSolventAtomMass3[pos];
            x3                                      = cSim.pConstraintSolventAtomX3[pos];
            y3                                      = cSim.pConstraintSolventAtomY3[pos];
            z3                                      = cSim.pConstraintSolventAtomZ3[pos];
            ox                                     += m3 * x3;
            oy                                     += m3 * y3;
            oz                                     += m3 * z3;
        }
        double x4, y4, z4;
        if (atoms == 4)
        {
            PMEDouble m4                            = cSim.pConstraintSolventAtomMass4[pos];
            x4                                      = cSim.pConstraintSolventAtomX4[pos];
            y4                                      = cSim.pConstraintSolventAtomY4[pos];
            z4                                      = cSim.pConstraintSolventAtomZ4[pos];
            ox                                     += m4 * x4;
            oy                                     += m4 * y4;
            oz                                     += m4 * z4;
        }    
        
        // Calculate change in COM
        ox                                         *= invMass;
        oy                                         *= invMass;
        oz                                         *= invMass;
        PMEDouble cx                                = sLast_recip[0] * ox + sLast_recip[3] * oy + sLast_recip[6] * oz;
        PMEDouble cy                                =                       sLast_recip[4] * oy + sLast_recip[7] * oz;
        PMEDouble cz                                =                                             sLast_recip[8] * oz;  
        cx                                          = sUcell[0] * cx + sUcell[1] * cy + sUcell[2] * cz;
        cy                                          =                  sUcell[4] * cy + sUcell[5] * cz;
        cz                                          =                                   sUcell[8] * cz;
        cSim.pConstraintSolventCOMX[pos]            = cx;
        cSim.pConstraintSolventCOMY[pos]            = cy;
        cSim.pConstraintSolventCOMZ[pos]            = cz;
        cx                                         -= ox;
        cy                                         -= oy;
        cz                                         -= oz;
        
        // Shift atomic coordinates
        cSim.pConstraintSolventAtomX1[pos]          = x1 + cx;
        cSim.pConstraintSolventAtomY1[pos]          = y1 + cy;
        cSim.pConstraintSolventAtomZ1[pos]          = z1 + cz;
        if (constraint.x != -1)
        {
            PMEDouble2 constraint1                  = cSim.pConstraint1[constraint.x];
            PMEDouble2 constraint2                  = {y1 + cy, z1 + cz};
            constraint1.y                           = x1 + cx;
            cSim.pConstraint2[constraint.x]         = constraint2;
            cSim.pConstraint1[constraint.x]         = constraint1;                
        }

        if (atoms >= 2)
        {
            cSim.pConstraintSolventAtomX2[pos]      = x2 + cx;
            cSim.pConstraintSolventAtomY2[pos]      = y2 + cy;
            cSim.pConstraintSolventAtomZ2[pos]      = z2 + cz;
            if (constraint.y != -1)
            {
                PMEDouble2 constraint1              = cSim.pConstraint1[constraint.y];
                PMEDouble2 constraint2              = {y2 + cy, z2 + cz};
                constraint1.y                       = x2 + cx;
                cSim.pConstraint2[constraint.y]     = constraint2;
                cSim.pConstraint1[constraint.y]     = constraint1;                
            }
        }
        if (atoms >= 3)
        {
            cSim.pConstraintSolventAtomX3[pos]      = x3 + cx;
            cSim.pConstraintSolventAtomY3[pos]      = y3 + cy;
            cSim.pConstraintSolventAtomZ3[pos]      = z3 + cz;
            if (constraint.z != -1)
            {
                PMEDouble2 constraint1              = cSim.pConstraint1[constraint.z];
                PMEDouble2 constraint2              = {y3 + cy, z3 + cz};
                constraint1.y                       = x3 + cx;
                cSim.pConstraint2[constraint.z]     = constraint2;
                cSim.pConstraint1[constraint.z]     = constraint1;                
            }
        }
        if (atoms == 4)
        {
            cSim.pConstraintSolventAtomX4[pos]      = x4 + cx;
            cSim.pConstraintSolventAtomY4[pos]      = y4 + cy;
            cSim.pConstraintSolventAtomZ4[pos]      = z4 + cz;
            if (constraint.w != -1)
            {
                PMEDouble2 constraint1              = cSim.pConstraint1[constraint.w];
                PMEDouble2 constraint2              = {y4 + cy, z4 + cz};
                constraint1.y                       = x4 + cx;
                cSim.pConstraint2[constraint.w]     = constraint2;
                cSim.pConstraint1[constraint.w]     = constraint1;                
            }
        }
        
        pos                                        += blockDim.x * gridDim.x;
    }
}

