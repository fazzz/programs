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
struct ScalarSumAtom {
    PMEDouble energy;
#ifdef PME_VIRIAL    
    PMEDouble vir_11;
    PMEDouble vir_22;
    PMEDouble vir_33;
#endif
};
#if defined(PME_ENERGY) || defined(PME_VIRIAL)
#if (__CUDA_ARCH__ >= 200)
__shared__ ScalarSumAtom sE[SM_2X_UPDATE_THREADS_PER_BLOCK];
#else
__shared__ ScalarSumAtom sE[SM_13_UPDATE_THREADS_PER_BLOCK];
#endif
#endif

#if (__CUDA_ARCH__ >= 200) && defined(use_DPDP)
#define PREFACSIZE 2048
#elif (__CUDA_ARCH__ >= 200)
#define PREFACSIZE 4096
#elif defined(PME_VIRIAL) && defined(use_DPDP)
#define PREFACSIZE 1010
#elif defined(PME_VIRIAL)
#define PREFACSIZE 2030
#elif defined(use_DPDP)
#define PREFACSIZE 1536
#else
#define PREFACSIZE 3072
#endif
__shared__ PMEFloat sPrefac[PREFACSIZE];
#ifdef PME_VIRIAL
__shared__ PMEFloat sRecipf[9];
#endif
    int increment                               = blockDim.x * gridDim.x;
    int zIncrement                              = increment / (cSim.fft_y_dim * cSim.fft_x_dim);
    increment                                  -= zIncrement * cSim.fft_y_dim * cSim.fft_x_dim;
    int yIncrement                              = increment / cSim.fft_x_dim;
    int xIncrement                              = increment - yIncrement * cSim.fft_x_dim;
    xIncrement--;
    yIncrement--;
    zIncrement--;
    
    // Insure q[0][0][0] == 0
    if ((blockIdx.x == 0) && (threadIdx.x == 0))
    {
        PMEComplex cmplx                        = {(PMEFloat)0.0, (PMEFloat)0.0};   
        cSim.pXYZ_qt[0]                         = cmplx;  
    }

#ifdef PME_ENERGY
    // Initialize energy parameters
    sE[threadIdx.x].energy                      = (PMEDouble)0.0;
#endif

#ifdef PME_VIRIAL  
    // Initialize energy parameters  
    sE[threadIdx.x].vir_11                      = 0.0;
    sE[threadIdx.x].vir_22                      = 0.0;
    sE[threadIdx.x].vir_33                      = 0.0;
    if (threadIdx.x < 9)
        sRecipf[threadIdx.x]                    = cSim.pNTPData->recipf[threadIdx.x];
    __syncthreads();
#endif    
    
    if (cSim.nSum <= PREFACSIZE)
    {    
        // Load prefac data
        unsigned int pos = threadIdx.x;
        while (pos < cSim.nSum)
        {
            sPrefac[pos]                            = cSim.pPrefac1[pos];
            pos                                    += blockDim.x;
        }
        __syncthreads();
        PMEFloat* psPrefac2                         = &sPrefac[cSim.n2Offset];
        PMEFloat* psPrefac3                         = &sPrefac[cSim.n3Offset];

        // Calculate initial position
        int offset                                  = blockIdx.x * blockDim.x + threadIdx.x + 1;
        int k3                                      = offset / (cSim.fft_y_dim * cSim.fft_x_dim);
        offset                                     -= k3 * cSim.fft_y_dim * cSim.fft_x_dim;
        int k2                                      = offset / cSim.fft_x_dim;
        int k1                                      = offset - k2 * cSim.fft_x_dim;
        int m1, m2, m3;
     
        while (k3 < cSim.fft_z_dim)
        {
        
            // Read data       
            unsigned int pos                        = (k3 * cSim.fft_y_dim + k2) * cSim.fft_x_dim + k1;
            PMEComplex q                            = cSim.pXYZ_qt[pos];        
            k1++;
            k2++;
            k3++;
        
            // Generate internal coordinates   
            m1                                      = k1 - 1;
            m2                                      = k2 - 1;
            m3                                      = k3 - 1;
            if (k1 > cSim.nf1)
                m1                                 -= cSim.nfft1;
            if (k2 > cSim.nf2)
                m2                                 -= cSim.nfft2;
            if (k3 > cSim.nf3)
                m3                                 -= cSim.nfft3;
#ifdef PME_VIRIAL
            PMEFloat mhat1                          = sRecipf[0] * m1;
            PMEFloat mhat2                          = sRecipf[3] * m1 + sRecipf[4] * m2;
            PMEFloat mhat3                          = sRecipf[6] * m1 + sRecipf[7] * m2 + sRecipf[8] * m3;
#else
            PMEFloat mhat1                          = cSim.recipf[0][0] * m1;
            PMEFloat mhat2                          = cSim.recipf[1][0] * m1 + cSim.recipf[1][1] * m2;
            PMEFloat mhat3                          = cSim.recipf[2][0] * m1 + cSim.recipf[2][1] * m2 + cSim.recipf[2][2] * m3;
#endif

            PMEFloat msq                            = mhat1 * mhat1 + mhat2 * mhat2 + mhat3 * mhat3;

            PMEFloat msq_inv                        = (PMEFloat)1.0 / msq;

            // Getting the exponential via table lookup is currently not done
            // for nonorthogonal unit cells.

            PMEFloat eterm                          = exp(-cSim.fac * msq) * sPrefac[k1] * psPrefac2[k2] * psPrefac3[k3] * cSim.pi_vol_inv * msq_inv;
#ifdef PME_VIRIAL        
            PMEFloat vterm                          = cSim.fac2 + 2.0 * msq_inv;
#endif

#if defined(PME_ENERGY) || defined (PME_VIRIAL)
            PMEFloat eterm_struc2                   = eterm * (q.x * q.x + q.y * q.y);
#endif
#ifdef PME_ENERGY
            sE[threadIdx.x].energy                 += eterm_struc2;
#endif
#ifdef PME_VIRIAL
            sE[threadIdx.x].vir_11                 += eterm_struc2 * (vterm * mhat1 * mhat1 - (PMEDouble)1.0);
            sE[threadIdx.x].vir_22                 += eterm_struc2 * (vterm * mhat2 * mhat2 - (PMEDouble)1.0);
            sE[threadIdx.x].vir_33                 += eterm_struc2 * (vterm * mhat3 * mhat3 - (PMEDouble)1.0);
#endif
#if defined(PME_ENERGY) || defined(PME_VIRIAL)
            if ((k1 > 1) && (k1 <= cSim.nfft1))
            {
                int k1s, k2s, k3s;
                int m1s, m2s, m3s;
                
                k1s                                 = cSim.nfft1 - k1 + 2;
                k2s                                 = ((cSim.nfft2 - k2 + 1) % cSim.nfft2) + 1;
                k3s                                 = ((cSim.nfft3 - k3 + 1) % cSim.nfft3) + 1;

               
                m1s                                 = k1s - 1;
                m2s                                 = k2s - 1;
                m3s                                 = k3s - 1;
                if (k1s > cSim.nf1)
                    m1s                            -= cSim.nfft1; 
                if (k2s > cSim.nf2)
                    m2s                            -= cSim.nfft2; 
                if (k3s > cSim.nf3)
                    m3s                            -= cSim.nfft3;

#ifdef PME_VIRIAL
                PMEFloat mhat1s                     = sRecipf[0] * m1s; 

                PMEFloat mhat2s                     = sRecipf[3] * m1s +
                                                      sRecipf[4] * m2s;

                PMEFloat mhat3s                     = sRecipf[6] * m1s +
                                                      sRecipf[7] * m2s +
                                                      sRecipf[8] * m3s;
#else
                PMEFloat mhat1s                     = cSim.recipf[0][0] * m1s; 

                PMEFloat mhat2s                     = cSim.recipf[1][0] * m1s +
                                                      cSim.recipf[1][1] * m2s;

                PMEFloat mhat3s                     = cSim.recipf[2][0] * m1s +
                                                      cSim.recipf[2][1] * m2s +
                                                      cSim.recipf[2][2] * m3s;
#endif

                PMEFloat msqs                       = mhat1s * mhat1s + mhat2s * mhat2s + mhat3s * mhat3s;

                PMEFloat msqs_inv                   = (PMEFloat)1.0 / msqs;

                // Getting the exponential via table lookup is currently not done
                // for nonorthogonal unit cells.

                PMEFloat eterms                     = exp(-cSim.fac * msqs) * sPrefac[k1s] * psPrefac2[k2s] * psPrefac3[k3s] * cSim.pi_vol_inv * msqs_inv;

#ifdef PME_VIRIAL
                PMEFloat vterms                     = cSim.fac2 + (PMEFloat)2.0 * msqs_inv;
#endif
                PMEFloat eterms_struc2s             = eterms * (q.x * q.x + q.y * q.y);
#ifdef PME_ENERGY
                sE[threadIdx.x].energy             += eterms_struc2s;
#endif
#ifdef PME_VIRIAL
                sE[threadIdx.x].vir_11             += eterms_struc2s * (vterms * mhat1s * mhat1s - (PMEDouble)1.0);
                sE[threadIdx.x].vir_22             += eterms_struc2s * (vterms * mhat2s * mhat2s - (PMEDouble)1.0);
                sE[threadIdx.x].vir_33             += eterms_struc2s * (vterms * mhat3s * mhat3s - (PMEDouble)1.0);
#endif            
            }
#endif // defined(PME_VIRIAL) || defined(PME_ENERGY)
            q.x                                    *= eterm;
            q.y                                    *= eterm;  
            cSim.pXYZ_qt[pos]                       = q;    
        
            // Increment position
            k1                                     += xIncrement;
            if (k1 >= cSim.fft_x_dim)
            {
                k1                                 -= cSim.fft_x_dim;
                k2++;
            }
            k2                                     += yIncrement;
            if (k2 >= cSim.fft_y_dim)
            {
                k2                                 -= cSim.fft_y_dim;
                k3++;
            }
            k3                                     += zIncrement;
        } 
    }
    else
    {
        // Calculate initial position
        int offset                                  = blockIdx.x * blockDim.x + threadIdx.x + 1;
        int k3                                      = offset / (cSim.fft_y_dim * cSim.fft_x_dim);
        offset                                     -= k3 * cSim.fft_y_dim * cSim.fft_x_dim;
        int k2                                      = offset / cSim.fft_x_dim;
        int k1                                      = offset - k2 * cSim.fft_x_dim;
        int m1, m2, m3;
     
        while (k3 < cSim.fft_z_dim)
        {
        
            // Read data     
            unsigned int pos                        = (k3 * cSim.fft_y_dim + k2) * cSim.fft_x_dim + k1;
            PMEComplex q                            = cSim.pXYZ_qt[pos];
       
            k1++;
            k2++;
            k3++;
        
            // Generate internal coordinates   
            m1                                      = k1 - 1;
            m2                                      = k2 - 1;
            m3                                      = k3 - 1;
            if (k1 > cSim.nf1)
                m1                                 -= cSim.nfft1;
            if (k2 > cSim.nf2)
                m2                                 -= cSim.nfft2;
            if (k3 > cSim.nf3)
                m3                                 -= cSim.nfft3;
      
#ifdef PME_VIRIAL
            PMEFloat mhat1                          = sRecipf[0] * m1;
            PMEFloat mhat2                          = sRecipf[3] * m1 + sRecipf[4] * m2;
            PMEFloat mhat3                          = sRecipf[6] * m1 + sRecipf[7] * m2 + sRecipf[8] * m3;
#else
            PMEFloat mhat1                          = cSim.recipf[0][0] * m1;
            PMEFloat mhat2                          = cSim.recipf[1][0] * m1 + cSim.recipf[1][1] * m2;
            PMEFloat mhat3                          = cSim.recipf[2][0] * m1 + cSim.recipf[2][1] * m2 + cSim.recipf[2][2] * m3;
#endif

            PMEFloat msq                            = mhat1 * mhat1 + mhat2 * mhat2 + mhat3 * mhat3;

            PMEFloat msq_inv                        = (PMEFloat)1.0 / msq;

            // Getting the exponential via table lookup is currently not done
            // for nonorthogonal unit cells.

            PMEFloat eterm                          = exp(-cSim.fac * msq) * cSim.pPrefac1[k1] * cSim.pPrefac2[k2] * cSim.pPrefac3[k3] * cSim.pi_vol_inv * msq_inv;

#ifdef PME_VIRIAL        
            PMEFloat vterm                          = cSim.fac2 + (PMEFloat)2.0 * msq_inv;
#endif
#if defined(PME_VIRIAL) || defined(PME_ENERGY)
            PMEFloat eterm_struc2                   = eterm * (q.x * q.x + q.y * q.y);
#endif
#ifdef PME_ENERGY
            sE[threadIdx.x].energy                 += eterm_struc2;
#endif
#ifdef PME_VIRIAL
            sE[threadIdx.x].vir_11                 += eterm_struc2 * (vterm * mhat1 * mhat1 - (PMEDouble)1.0);
            sE[threadIdx.x].vir_22                 += eterm_struc2 * (vterm * mhat2 * mhat2 - (PMEDouble)1.0);
            sE[threadIdx.x].vir_33                 += eterm_struc2 * (vterm * mhat3 * mhat3 - (PMEDouble)1.0);
#endif

#if defined(PME_VIRIAL) || defined(PME_ENERGY)
            if ((k1 > 1) && (k1 <= cSim.nfft1))
            {
                int k1s, k2s, k3s;
                int m1s, m2s, m3s;
                
                k1s                                 = cSim.nfft1 - k1 + 2;
                k2s                                 = ((cSim.nfft2 - k2 + 1) % cSim.nfft2) + 1;
                k3s                                 = ((cSim.nfft3 - k3 + 1) % cSim.nfft3) + 1;

               
                m1s                                 = k1s - 1;
                m2s                                 = k2s - 1;
                m3s                                 = k3s - 1;
                if (k1s > cSim.nf1)
                    m1s                            -= cSim.nfft1; 
                if (k2s > cSim.nf2)
                    m2s                            -= cSim.nfft2; 
                if (k3s > cSim.nf3)
                    m3s                            -= cSim.nfft3;

#ifdef PME_VIRIAL
                PMEFloat mhat1s                     = sRecipf[0] * m1s;

                PMEFloat mhat2s                     = sRecipf[3] * m1s +
                                                      sRecipf[4] * m2s;

                PMEFloat mhat3s                     = sRecipf[6] * m1s +
                                                      sRecipf[7] * m2s +
                                                      sRecipf[8] * m3s;
#else
                PMEFloat mhat1s                     = cSim.recipf[0][0] * m1s;

                PMEFloat mhat2s                     = cSim.recipf[1][0] * m1s +
                                                      cSim.recipf[1][1] * m2s;

                PMEFloat mhat3s                     = cSim.recipf[2][0] * m1s +
                                                      cSim.recipf[2][1] * m2s +
                                                      cSim.recipf[2][2] * m3s;
#endif

                PMEFloat msqs                       = mhat1s * mhat1s + mhat2s * mhat2s + mhat3s * mhat3s;

                PMEFloat msqs_inv                   = (PMEFloat)1.0 / msqs;

                // Getting the exponential via table lookup is currently not done
                // for nonorthogonal unit cells.

                PMEFloat eterms                     = exp(-cSim.fac * msqs) * cSim.pPrefac1[k1s] * cSim.pPrefac2[k2s] * cSim.pPrefac3[k3s] * cSim.pi_vol_inv * msqs_inv;

#ifdef PME_VIRIAL
                PMEFloat vterms                     = cSim.fac2 + (PMEFloat)2.0 * msqs_inv;
#endif
                PMEFloat eterms_struc2s             = eterms * (q.x * q.x + q.y * q.y);
#ifdef PME_ENERGY
                sE[threadIdx.x].energy             += eterms_struc2s;
#endif
#ifdef PME_VIRIAL
                sE[threadIdx.x].vir_11             += eterms_struc2s * (vterms * mhat1s * mhat1s - (PMEDouble)1.0);
                sE[threadIdx.x].vir_22             += eterms_struc2s * (vterms * mhat2s * mhat2s - (PMEDouble)1.0);
                sE[threadIdx.x].vir_33             += eterms_struc2s * (vterms * mhat3s * mhat3s - (PMEDouble)1.0);
#endif            
            }
#endif // defined(PME_VIRIAL) || defined(PME_ENERGY)
            q.x                                    *= eterm;
            q.y                                    *= eterm;    
            cSim.pXYZ_qt[pos]                       = q;    
        
            // Increment position
            k1                                     += xIncrement;
            if (k1 >= cSim.fft_x_dim)
            {
                k1                                 -= cSim.fft_x_dim;
                k2++;
            }
            k2                                     += yIncrement;
            if (k2 >= cSim.fft_y_dim)
            {
                k2                                 -= cSim.fft_y_dim;
                k3++;
            }
            k3                                     += zIncrement;
        }    
    }
        
#if defined(PME_ENERGY) || defined(PME_VIRIAL)        
    // Reduce virial and energy
    __syncthreads();
    int m                                           = 1;
    while (m < blockDim.x)
    {
        int p                                       = threadIdx.x + m;
#ifdef PME_ENERGY  
        PMEDouble energy                            = ((p < blockDim.x) ? sE[p].energy : (PMEDouble)0.0);
#endif
#ifdef PME_VIRIAL      
        PMEDouble vir_11                            = ((p < blockDim.x) ? sE[p].vir_11 : (PMEDouble)0.0);
        PMEDouble vir_22                            = ((p < blockDim.x) ? sE[p].vir_22 : (PMEDouble)0.0);
        PMEDouble vir_33                            = ((p < blockDim.x) ? sE[p].vir_33 : (PMEDouble)0.0);
#endif        
        __syncthreads();
#ifdef PME_ENERGY
        sE[threadIdx.x].energy                     += energy;
#endif
#ifdef PME_VIRIAL
        sE[threadIdx.x].vir_11                     += vir_11;
        sE[threadIdx.x].vir_22                     += vir_22;
        sE[threadIdx.x].vir_33                     += vir_33;
#endif        
        __syncthreads();
        m                                      *= 2;
    }
#ifdef PME_ENERGY
    unsigned long long int eer                  = (unsigned long long int)(fabs((PMEDouble)0.5 * sE[threadIdx.x].energy) * ENERGYSCALE + (PMEDouble)0.5);
    if (sE[threadIdx.x].energy < (PMEDouble)0.0)
        eer                                     = 0ull - eer;
#endif
#ifdef PME_VIRIAL        
    unsigned long long int vir_11               = (unsigned long long int)(fabs((PMEDouble)0.5 * sE[threadIdx.x].vir_11) * ENERGYSCALE + (PMEDouble)0.5);
    if (sE[threadIdx.x].vir_11 < (PMEDouble)0.0)
        vir_11                                  = 0ull - vir_11; 
    unsigned long long int vir_22               = (unsigned long long int)(fabs((PMEDouble)0.5 * sE[threadIdx.x].vir_22) * ENERGYSCALE + (PMEDouble)0.5);
    if (sE[threadIdx.x].vir_22 < (PMEDouble)0.0)
        vir_22                                  = 0ull - vir_22; 
    unsigned long long int vir_33               = (unsigned long long int)(fabs((PMEDouble)0.5 * sE[threadIdx.x].vir_33) * ENERGYSCALE + (PMEDouble)0.5);
    if (sE[threadIdx.x].vir_33 < (PMEDouble)0.0)
        vir_33                                  = 0ull - vir_33;                        
#endif

    // Write out energies
    if (threadIdx.x == 0)
    {
#ifdef PME_ENERGY
       atomicAdd(cSim.pEER, eer);
#endif
#ifdef PME_VIRIAL             
       atomicAdd(cSim.pVirial_11, vir_11);
       atomicAdd(cSim.pVirial_22, vir_22);
       atomicAdd(cSim.pVirial_33, vir_33);     
#endif                        
    }
#endif // defined(PME_ENERGY) || defined(PME_VIRIAL)      
#undef PREFACSIZE
}
