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
#if (__CUDA_ARCH__ >= 200)
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
#if (__CUDA_ARCH__ >= 200)
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
#if (__CUDA_ARCH__ >= 200)
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
        invMass                                    *= ONEOVERENERGYSCALESQUARED;
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
        unsigned long long int val1                         = (unsigned long long int)(fabs(sE[0].EKCOMX) * ENERGYSCALE + (PMEDouble)0.5);
        unsigned long long int val2                         = (unsigned long long int)(fabs(sE[0].EKCOMY) * ENERGYSCALE + (PMEDouble)0.5);
        unsigned long long int val3                         = (unsigned long long int)(fabs(sE[0].EKCOMZ) * ENERGYSCALE + (PMEDouble)0.5);
        if (sE[0].EKCOMX < (PMEDouble)0.0)
            val1                                            = 0ull - val1;
        if (val1 != 0)
            atomicAdd(cSim.pEKCOMX, val1);
        if (sE[0].EKCOMY < (PMEDouble)0.0)
            val2                                            = 0ull - val2;
        if (val2 != 0)
            atomicAdd(cSim.pEKCOMY, val2);
        if (sE[0].EKCOMZ < (PMEDouble)0.0)
            val3                                            = 0ull - val3;
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

void kCalculateConstraintsCOM(gpuContext gpu)
{
    if (gpu->sim.constraintSoluteMolecules <= gpu->maxSoluteMolecules)
        kPMECalculateConstraintsCOM_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    else
        kPMECalculateConstraintsCOMLarge_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kCalculateConstraintsCOM");    
}

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_THREADS_PER_BLOCK, 1)
#endif
kReduceSoluteConstraintsCOM_kernel()
{
    unsigned int pos                                = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.constraintSoluteMolecules)
    {
        PMEDouble invMass                           = cSim.pConstraintSoluteInvMass[pos];        
        PMEUllInt ullCOMX                           = cSim.pConstraintSoluteUllCOMX[pos];
        PMEUllInt ullCOMY                           = cSim.pConstraintSoluteUllCOMY[pos];
        PMEUllInt ullCOMZ                           = cSim.pConstraintSoluteUllCOMZ[pos];
        cSim.pConstraintSoluteUllCOMX[pos]          = 0;
        cSim.pConstraintSoluteUllCOMY[pos]          = 0;
        cSim.pConstraintSoluteUllCOMZ[pos]          = 0;
        invMass                                    *= ONEOVERENERGYSCALE;
        PMEDouble CX, CY, CZ;
        if (ullCOMX >= 0x8000000000000000ull)
            CX                                      = -(PMEDouble)(ullCOMX ^ 0xffffffffffffffffull);
        else   
            CX                                      =  (PMEDouble)ullCOMX;
        cSim.pConstraintSoluteCOMX[pos]             = invMass * CX;
        if (ullCOMY >= 0x8000000000000000ull)
            CY                                      = -(PMEDouble)(ullCOMY ^ 0xffffffffffffffffull);
        else
            CY                                      =  (PMEDouble)ullCOMY;
        cSim.pConstraintSoluteCOMY[pos]             = invMass * CY;
        if (ullCOMZ >= 0x8000000000000000ull)
            CZ                                      = -(PMEDouble)(ullCOMZ ^ 0xffffffffffffffffull);
        else
            CZ                                      =  (PMEDouble)ullCOMZ;
        cSim.pConstraintSoluteCOMZ[pos]             = invMass * CZ;
        pos                                        += blockDim.x * gridDim.x;
    }
}

void kReduceSoluteConstraintsCOM(gpuContext gpu)
{
    kReduceSoluteConstraintsCOM_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kReduceSoluteConstraintsCOM");    
}

void kCalculateSoluteConstraintsCOM(gpuContext gpu)
{
    if (gpu->sim.constraintSoluteMolecules <= gpu->maxSoluteMolecules)
        kPMECalculateSoluteConstraintsCOM_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    else
        kPMECalculateSoluteConstraintsCOMLarge_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kCalculateSoluteConstraintsCOM");    
}

void kPressureScaleConstraintCoordinates(gpuContext gpu)
{
    if (gpu->sim.constraintSoluteMolecules <= gpu->maxPSSoluteMolecules)
        kPressureScaleConstraintCoordinates_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    else
        kPressureScaleConstraintCoordinatesLarge_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    LAUNCHERROR("kPressureScaleConstraintCoordinates");    
}
