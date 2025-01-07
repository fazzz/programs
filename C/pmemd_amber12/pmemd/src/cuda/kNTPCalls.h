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

void kCalculateCOM(gpuContext gpu)
{
    if (gpu->sim.soluteMolecules <= gpu->maxSoluteMolecules)
        kPMECalculateCOM_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    else
        kPMECalculateCOMLarge_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kCalculateCOM");
}

__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kReduceSoluteCOM_kernel()
{
    unsigned int pos                                = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.soluteMolecules)
    {
        PMEDouble invMass                           = cSim.pSoluteInvMass[pos];        
        PMEUllInt ullCOMX                           = cSim.pSoluteUllCOMX[pos];
        PMEUllInt ullCOMY                           = cSim.pSoluteUllCOMY[pos];
        PMEUllInt ullCOMZ                           = cSim.pSoluteUllCOMZ[pos];
        cSim.pSoluteUllCOMX[pos]                    = 0;
        cSim.pSoluteUllCOMY[pos]                    = 0;
        cSim.pSoluteUllCOMZ[pos]                    = 0;
        invMass                                    *= ONEOVERENERGYSCALE;
        PMEDouble CX, CY, CZ;
        if (ullCOMX >= 0x8000000000000000ull)
            CX                                      = -(PMEDouble)(ullCOMX ^ 0xffffffffffffffffull);
        else   
            CX                                      =  (PMEDouble)ullCOMX;
        cSim.pSoluteCOMX[pos]                       = invMass * CX;
        if (ullCOMY >= 0x8000000000000000ull)
            CY                                      = -(PMEDouble)(ullCOMY ^ 0xffffffffffffffffull);
        else
            CY                                      =  (PMEDouble)ullCOMY;
        cSim.pSoluteCOMY[pos]                       = invMass * CY;
        if (ullCOMZ >= 0x8000000000000000ull)
            CZ                                      = -(PMEDouble)(ullCOMZ ^ 0xffffffffffffffffull);
        else
            CZ                                      =  (PMEDouble)ullCOMZ;
        cSim.pSoluteCOMZ[pos]                       = invMass * CZ;
        pos                                        += blockDim.x * gridDim.x;
    }
}

void kReduceSoluteCOM(gpuContext gpu)
{
    kReduceSoluteCOM_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kReduceSoluteCOM");   
}




void kCalculateSoluteCOM(gpuContext gpu)
{
    if (gpu->sim.soluteMolecules <= gpu->maxSoluteMolecules)
        kPMECalculateSoluteCOM_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    else
        kPMECalculateSoluteCOMLarge_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kCalculateSoluteCOM");
}

void kCalculateCOMKineticEnergy(gpuContext gpu)
{
    if (gpu->sim.soluteMolecules <= gpu->maxSoluteMolecules)
        kPMECalculateCOMKineticEnergy_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    else
        kPMECalculateCOMKineticEnergyLarge_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kCalculateCOMKineticEnergy"); 
}


__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kReduceCOMKineticEnergy_kernel()
{
struct COMKineticEnergy 
{
    double EKCOMX;
    double EKCOMY;
    double EKCOMZ;
};
#if (__CUDA_ARCH__ >= 300)
__shared__ volatile COMKineticEnergy sE[SM_3X_THREADS_PER_BLOCK];
#elif (__CUDA_ARCH__ >= 200)
__shared__ volatile COMKineticEnergy sE[SM_2X_THREADS_PER_BLOCK];
#else
__shared__ volatile COMKineticEnergy sE[SM_13_THREADS_PER_BLOCK];
#endif
    unsigned int pos                                = blockIdx.x * blockDim.x + threadIdx.x;
    PMEDouble EKCOMX                                = (PMEDouble)0.0;
    PMEDouble EKCOMY                                = (PMEDouble)0.0;
    PMEDouble EKCOMZ                                = (PMEDouble)0.0;
    while (pos < cSim.soluteMolecules)
    {
        PMEDouble invMass                           = cSim.pSoluteInvMass[pos];        
        PMEUllInt ullEKCOMX                         = cSim.pSoluteUllEKCOMX[pos];
        PMEUllInt ullEKCOMY                         = cSim.pSoluteUllEKCOMY[pos];
        PMEUllInt ullEKCOMZ                         = cSim.pSoluteUllEKCOMZ[pos];
        cSim.pSoluteUllEKCOMX[pos]                  = 0;
        cSim.pSoluteUllEKCOMY[pos]                  = 0;
        cSim.pSoluteUllEKCOMZ[pos]                  = 0;
        invMass                                    *= ONEOVERFORCESCALESQUARED;
        PMEDouble EX, EY, EZ;
        if (ullEKCOMX >= 0x8000000000000000ull)
            EX                                      = -(PMEDouble)(ullEKCOMX ^ 0xffffffffffffffffull);
        else   
            EX                                      =  (PMEDouble)ullEKCOMX;
        EKCOMX                                     += invMass * EX * EX;
        if (ullEKCOMY >= 0x8000000000000000ull)
            EY                                      = -(PMEDouble)(ullEKCOMY ^ 0xffffffffffffffffull);
        else
            EY                                      =  (PMEDouble)ullEKCOMY;
        EKCOMY                                     += invMass * EY * EY;
        if (ullEKCOMZ >= 0x8000000000000000ull)
            EZ                                      = -(PMEDouble)(ullEKCOMZ ^ 0xffffffffffffffffull);
        else
            EZ                                      =  (PMEDouble)ullEKCOMZ;
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
        unsigned long long int val1                 = llitoulli(lliroundd(sE[0].EKCOMX * FORCESCALE));
        unsigned long long int val2                 = llitoulli(lliroundd(sE[0].EKCOMY * FORCESCALE));
        unsigned long long int val3                 = llitoulli(lliroundd(sE[0].EKCOMZ * FORCESCALE));
        if (val1 != 0)
            atomicAdd(cSim.pEKCOMX, val1);
        if (val2 != 0)
            atomicAdd(cSim.pEKCOMY, val2);
        if (val3 != 0)
            atomicAdd(cSim.pEKCOMZ, val3);
   }
}

void kReduceCOMKineticEnergy(gpuContext gpu)
{  
    kReduceCOMKineticEnergy_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kReduceCOMKineticEnergy");   
}

void kCalculateMolecularVirial(gpuContext gpu)
{
    if (gpu->sim.soluteMolecules <= gpu->maxSoluteMolecules)
        kCalculateMolecularVirial_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    else
        kCalculateMolecularVirialLarge_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kCalculateMolecularVirial");
}

void kPressureScaleCoordinates(gpuContext gpu)
{
    if (gpu->sim.soluteMolecules <= gpu->maxPSSoluteMolecules)
        kPressureScaleCoordinates_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    else
        kPressureScaleCoordinatesLarge_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kPressureScaleCoordinates");
}


__global__ void
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kPressureScaleConstraintCoordinates_kernel()
{
__shared__ PMEDouble sUcell[9];
    // Read transformation matrices
    if (threadIdx.x < 9)
    {
        sUcell[threadIdx.x]                         = cSim.pNTPData->ucell[threadIdx.x];
    }
    __syncthreads();


    // Iterate over constraints
    unsigned int pos                                = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.constraints)
    {
        // Read COM (pre-scaled into fractional coordinates so we're always starting with fresh data
        // in order to avoid heinous FPRE over long simulations in SPFP and SPDP mode.  Also always scale 
        // original input restraint coordinates that have been pre-centered so that the system does not 
        // gradually shrink towards the origin and subtly increase pressure over 1M iterations or so.  
        PMEDouble fx                                = cSim.pConstraintCOMX[pos];
        PMEDouble fy                                = cSim.pConstraintCOMY[pos];
        PMEDouble fz                                = cSim.pConstraintCOMZ[pos];
        PMEDouble x                                 = cSim.pConstraintAtomX[pos];
        PMEDouble y                                 = cSim.pConstraintAtomY[pos];
        PMEDouble z                                 = cSim.pConstraintAtomZ[pos];
        PMEDouble2 constraint1                      = cSim.pConstraint1[pos];

        // Calculate COM displacement
        PMEDouble cx                                = sUcell[0] * fx + sUcell[1] * fy + sUcell[2] * fz;
        PMEDouble cy                                =                  sUcell[4] * fy + sUcell[5] * fz;
        PMEDouble cz                                =                                   sUcell[8] * fz;
        x                                          += cx;
        y                                          += cy;
        z                                          += cz;
        PMEDouble2 constraint2                      = {y, z};
        constraint1.y                               = x;
        cSim.pConstraint2[pos]                      = constraint2;
        cSim.pConstraint1[pos]                      = constraint1;    
        pos                                        += blockDim.x * gridDim.x;
    }
}


void kPressureScaleConstraintCoordinates(gpuContext gpu)
{
    kPressureScaleConstraintCoordinates_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kPressureScaleConstraintCoordinates");    
}
