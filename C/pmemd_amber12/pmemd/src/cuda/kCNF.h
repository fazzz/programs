/***************************************************/
/*                                                  */
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
#define NMR_TIMED
#ifdef NMR_ENERGY
struct NMREnergy {
    long long int distance;
    long long int angle;
    long long int torsion;
};

#if (__CUDA_ARCH__ >= 300)
    __shared__ NMREnergy sE[SM_3X_LOCALFORCES_THREADS_PER_BLOCK];
#elif (__CUDA_ARCH__ >= 200)
    __shared__ NMREnergy sE[SM_2X_LOCALFORCES_THREADS_PER_BLOCK];
#else
    __shared__ NMREnergy sE[SM_13_LOCALFORCES_THREADS_PER_BLOCK];
#endif

    sE[threadIdx.x].distance                = (PMEDouble)0.0;
    sE[threadIdx.x].angle                   = (PMEDouble)0.0;
    sE[threadIdx.x].torsion                 = (PMEDouble)0.0;
#endif
    int pos                                 = blockIdx.x * blockDim.x + threadIdx.x; 

    // Calculate distance restraints
    if (pos < cSim.NMRDistanceOffset)
    {
        if (pos < cSim.NMRDistances)
        {
#ifdef NMR_NEIGHBORLIST
            int4 atom                       = cSim.pImageNMRDistanceID[pos];
#else
            int4 atom                       = cSim.pNMRDistanceID[pos];
#endif       

            PMEDouble2 R1R2;
            PMEDouble2 R3R4;
            PMEDouble2 K2K3;
#ifdef NMR_TIMED
            int2 Step                       = cSim.pNMRDistanceStep[pos];
            int Inc                         = cSim.pNMRDistanceInc[pos];
            // Skip restraint if not active
            if ((step < Step.x) || ((step > Step.y) && (Step.y > 0)))
                goto exit;

            // Read timed data
            PMEDouble2 R1R2Slp              = cSim.pNMRDistanceR1R2Slp[pos];
            PMEDouble2 R1R2Int              = cSim.pNMRDistanceR1R2Int[pos];
            PMEDouble2 R3R4Slp              = cSim.pNMRDistanceR3R4Slp[pos];
            PMEDouble2 R3R4Int              = cSim.pNMRDistanceR3R4Int[pos];
            PMEDouble2 K2K3Slp              = cSim.pNMRDistanceK2K3Slp[pos];
            PMEDouble2 K2K3Int              = cSim.pNMRDistanceK2K3Int[pos];
            
            // Calculate increment
            double dstep                    = step - (double)((step - Step.x) % abs(Inc));

            // Calculate restraint values
            R1R2.x                          = R1R2Slp.x * dstep + R1R2Int.x;
            R1R2.y                          = R1R2Slp.y * dstep + R1R2Int.y;
            R3R4.x                          = R3R4Slp.x * dstep + R3R4Int.x;
            R3R4.y                          = R3R4Slp.y * dstep + R3R4Int.y;
            if (Inc > 0)
            {
                K2K3.x                      = K2K3Slp.x * dstep + K2K3Int.x;
                K2K3.y                      = K2K3Slp.y * dstep + K2K3Int.y;
            }
            else
            {
                int nstepu                  = (step - Step.x) / abs(Inc);
                K2K3.x                      = K2K3Int.x * pow(K2K3Slp.x, nstepu);
                K2K3.y                      = K2K3Int.y * pow(K2K3Slp.y, nstepu);
            }
#else
            R1R2                            = cSim.pNMRDistanceR1R2[pos];
            R3R4                            = cSim.pNMRDistanceR3R4[pos];
            K2K3                            = cSim.pNMRDistanceK2K3[pos];
#endif

            //printf("%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n", R1R2.x, R1R2.y, R3R4.x, R3R4.y, K2K3.x, K2K3.y); 

#if defined(NODPTEXTURE)
#ifdef NMR_NEIGHBORLIST
            PMEDouble atomIX                = cSim.pImageX[atom.x];
            PMEDouble atomJX                = cSim.pImageX[atom.y];
            PMEDouble atomIY                = cSim.pImageY[atom.x];
            PMEDouble atomJY                = cSim.pImageY[atom.y];
            PMEDouble atomIZ                = cSim.pImageZ[atom.x];
            PMEDouble atomJZ                = cSim.pImageZ[atom.y];
#else
            PMEDouble atomIX                = cSim.pAtomX[atom.x];
            PMEDouble atomJX                = cSim.pAtomX[atom.y];
            PMEDouble atomIY                = cSim.pAtomY[atom.x];
            PMEDouble atomJY                = cSim.pAtomY[atom.y];
            PMEDouble atomIZ                = cSim.pAtomZ[atom.x];
            PMEDouble atomJZ                = cSim.pAtomZ[atom.y];
#endif
            PMEDouble xij                   = atomIX - atomJX;
            PMEDouble yij                   = atomIY - atomJY;    
#else        
            int2 iatomIX                    = tex1Dfetch(texref, atom.x);
            int2 iatomJX                    = tex1Dfetch(texref, atom.y);
            int2 iatomIY                    = tex1Dfetch(texref, atom.x + cSim.stride);
            int2 iatomJY                    = tex1Dfetch(texref, atom.y + cSim.stride);
            int2 iatomIZ                    = tex1Dfetch(texref, atom.x + cSim.stride2);
            int2 iatomJZ                    = tex1Dfetch(texref, atom.y + cSim.stride2);  
            PMEDouble atomIX                = __hiloint2double(iatomIX.y, iatomIX.x);
            PMEDouble atomJX                = __hiloint2double(iatomJX.y, iatomJX.x);
            PMEDouble xij                   = atomIX - atomJX;
            PMEDouble atomIY                = __hiloint2double(iatomIY.y, iatomIY.x);
            PMEDouble atomJY                = __hiloint2double(iatomJY.y, iatomJY.x);
            PMEDouble yij                   = atomIY - atomJY;            
            PMEDouble atomIZ                = __hiloint2double(iatomIZ.y, iatomIZ.x);
            PMEDouble atomJZ                = __hiloint2double(iatomJZ.y, iatomJZ.x);    
#endif                
            PMEDouble zij                   = atomIZ - atomJZ;
            PMEDouble rij                   = sqrt(xij * xij + yij * yij + zij * zij);
            PMEDouble df;
#ifdef NMR_ENERGY
            PMEDouble e;
#endif
            if (rij < R1R2.x)
            {
                PMEDouble dif               = R1R2.x - R1R2.y;
                df                          = 2.0 * K2K3.x * dif;
#ifdef NMR_ENERGY
                e                           = df * (rij - R1R2.x) + K2K3.x * dif * dif;
#endif
            }
            else if (rij < R1R2.y)
            {
                PMEDouble dif               = rij - R1R2.y;
                df                          = 2.0 * K2K3.x * dif;
#ifdef NMR_ENERGY
                e                           = K2K3.x * dif * dif;
#endif
            }
            else if (rij < R3R4.x)
            {
                df                          = 0.0;
#ifdef NMR_ENERGY
                e                           = 0.0;
#endif
            }
            else if (rij < R3R4.y)
            {
                PMEDouble dif               = rij - R3R4.x;
                df                          = 2.0 * K2K3.y * dif;
#ifdef NMR_ENERGY
                e                           = K2K3.y * dif * dif;
#endif
            }
            else
            {
                PMEDouble dif               = R3R4.y - R3R4.x;
                df                          = 2.0 * K2K3.y * dif;
#ifdef NMR_ENERGY
                e                           = df * (rij - R3R4.y) + K2K3.y * dif * dif;
#endif
            }

            if (cSim.bJar)
            {
                double fold                 = cSim.pNMRJarData[2]; 
                double work                 = cSim.pNMRJarData[3];
                double first                = cSim.pNMRJarData[4];
                double fcurr                = -2.0 * K2K3.x * (rij - R1R2.y);
                if (first == 0.0) {
                    fold                    = -fcurr;
                    cSim.pNMRJarData[4]     = 1.0;
                }
                work                       += 0.5 * (fcurr + fold) * cSim.drjar;
                cSim.pNMRJarData[0]         = R1R2.y;
                cSim.pNMRJarData[1]         = rij;
                cSim.pNMRJarData[2]         = fcurr;
                cSim.pNMRJarData[3]         = work;
            }
            
#ifdef NMR_ENERGY
#ifdef MPI
            if ((atom.x >= cSim.minLocalAtom) && (atom.x < cSim.maxLocalAtom))  
#endif
                sE[threadIdx.x].distance   += lliroundd(ENERGYSCALE * e);
#endif

            df                             *= 1.0 / rij;
            PMEDouble fx                    = df * xij;
            PMEDouble fy                    = df * yij;
            PMEDouble fz                    = df * zij;
            
#if defined(NMR_VIRIAL) || defined(use_SPFP)
            PMEAccumulator ifx              = lliroundd(fx * FORCESCALE);
            PMEAccumulator ify              = lliroundd(fy * FORCESCALE);
            PMEAccumulator ifz              = lliroundd(fz * FORCESCALE);
#endif

#ifdef MPI
            if ((atom.y >= cSim.minLocalAtom) && (atom.y < cSim.maxLocalAtom))
#endif
            {
#if defined(NMR_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom.y], llitoulli(ifx)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom.y], llitoulli(ify));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom.y], llitoulli(ifz));               
#elif defined (NMR_NEIGHBORLIST)                
                cSim.pForceXBuffer[atom.w]  = fx;
                cSim.pForceYBuffer[atom.w]  = fy;
                cSim.pForceZBuffer[atom.w]  = fz;
#else
                cSim.pForceXBuffer[atom.w] += fx;
                cSim.pForceYBuffer[atom.w] += fy;
                cSim.pForceZBuffer[atom.w] += fz; 
#endif
            }
#ifdef MPI
            if ((atom.x >= cSim.minLocalAtom) && (atom.x < cSim.maxLocalAtom))
#endif
            {     
#if defined(NMR_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom.x], llitoulli(-ifx)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom.x], llitoulli(-ify));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom.x], llitoulli(-ifz));   
#elif defined(NMR_NEIGHBORLIST)                   
                cSim.pForceXBuffer[atom.z]  = -fx;
                cSim.pForceYBuffer[atom.z]  = -fy;
                cSim.pForceZBuffer[atom.z]  = -fz;
#else
                cSim.pForceXBuffer[atom.z] += -fx;
                cSim.pForceYBuffer[atom.z] += -fy;
                cSim.pForceZBuffer[atom.z] += -fz;
#endif                
            }
        }
    }

    // Calculate angle restraints
    else if (pos < cSim.NMRAngleOffset)
    {
        pos                                -= cSim.NMRDistanceOffset;
        if (pos < cSim.NMRAngles)
        {
#ifdef NMR_NEIGHBORLIST        
            int4 atom1                      = cSim.pImageNMRAngleID1[pos];
#else
            int4 atom1                      = cSim.pNMRAngleID1[pos];
#endif
#ifdef NMR_NEIGHBORLIST
#if !defined(NMR_VIRIAL) && !defined(use_SPFP)
            int2 atom2                      = cSim.pImageNMRAngleID2[pos];               
#endif
#elif !defined(use_SPFP)
            int2 atom2                      = cSim.pNMRAngleID2[pos];       
#endif
            PMEDouble2 R1R2;
            PMEDouble2 R3R4;
            PMEDouble2 K2K3;
#ifdef NMR_TIMED
            int2 Step                       = cSim.pNMRAngleStep[pos];
            int Inc                         = cSim.pNMRAngleInc[pos];
            // Skip restraint if not active
            if ((step < Step.x) || ((step > Step.y) && (Step.y > 0)))
                goto exit;

            // Read timed data
            PMEDouble2 R1R2Slp              = cSim.pNMRAngleR1R2Slp[pos];
            PMEDouble2 R1R2Int              = cSim.pNMRAngleR1R2Int[pos];
            PMEDouble2 R3R4Slp              = cSim.pNMRAngleR3R4Slp[pos];
            PMEDouble2 R3R4Int              = cSim.pNMRAngleR3R4Int[pos];
            PMEDouble2 K2K3Slp              = cSim.pNMRAngleK2K3Slp[pos];
            PMEDouble2 K2K3Int              = cSim.pNMRAngleK2K3Int[pos];
            
            // Calculate increment
            double dstep                    = step - (double)((step - Step.x) % abs(Inc));

            // Calculate restraint values
            R1R2.x                          = R1R2Slp.x * dstep + R1R2Int.x;
            R1R2.y                          = R1R2Slp.y * dstep + R1R2Int.y;
            R3R4.x                          = R3R4Slp.x * dstep + R3R4Int.x;
            R3R4.y                          = R3R4Slp.y * dstep + R3R4Int.y;
            if (Inc > 0)
            {
                K2K3.x                      = K2K3Slp.x * dstep + K2K3Int.x;
                K2K3.y                      = K2K3Slp.y * dstep + K2K3Int.y;
            }
            else
            {
                int nstepu                  = (step - Step.x) / abs(Inc);
                K2K3.x                      = K2K3Int.x * pow(K2K3Slp.x, nstepu);
                K2K3.y                      = K2K3Int.y * pow(K2K3Slp.y, nstepu);
            }
#else
            R1R2                            = cSim.pNMRAngleR1R2[pos];
            R3R4                            = cSim.pNMRAngleR3R4[pos];
            K2K3                            = cSim.pNMRAngleK2K3[pos];
#endif

#if defined(NODPTEXTURE)
#ifdef NMR_NEIGHBORLIST
            PMEDouble atomIX                = cSim.pImageX[atom1.x];
            PMEDouble atomJX                = cSim.pImageX[atom1.y];
            PMEDouble atomKX                = cSim.pImageX[atom1.z];
            PMEDouble atomIY                = cSim.pImageY[atom1.x];
            PMEDouble atomJY                = cSim.pImageY[atom1.y];
            PMEDouble atomKY                = cSim.pImageY[atom1.z];
            PMEDouble atomIZ                = cSim.pImageZ[atom1.x];
            PMEDouble atomJZ                = cSim.pImageZ[atom1.y];
            PMEDouble atomKZ                = cSim.pImageZ[atom1.z];
#else
            PMEDouble atomIX                = cSim.pAtomX[atom1.x];
            PMEDouble atomJX                = cSim.pAtomX[atom1.y];
            PMEDouble atomKX                = cSim.pAtomX[atom1.z];
            PMEDouble atomIY                = cSim.pAtomY[atom1.x];
            PMEDouble atomJY                = cSim.pAtomY[atom1.y];
            PMEDouble atomKY                = cSim.pAtomY[atom1.z];
            PMEDouble atomIZ                = cSim.pAtomZ[atom1.x];
            PMEDouble atomJZ                = cSim.pAtomZ[atom1.y];
            PMEDouble atomKZ                = cSim.pAtomZ[atom1.z];
#endif
            PMEDouble xij                   = atomIX - atomJX;
            PMEDouble xkj                   = atomKX - atomJX;
            PMEDouble yij                   = atomIY - atomJY;
            PMEDouble ykj                   = atomKY - atomJY;
#else                         
            int2 iatomIX                    = tex1Dfetch(texref, atom1.x);
            int2 iatomJX                    = tex1Dfetch(texref, atom1.y);
            int2 iatomKX                    = tex1Dfetch(texref, atom1.z);
            int2 iatomIY                    = tex1Dfetch(texref, atom1.x + cSim.stride);
            int2 iatomJY                    = tex1Dfetch(texref, atom1.y + cSim.stride);
            int2 iatomKY                    = tex1Dfetch(texref, atom1.z + cSim.stride);
            int2 iatomIZ                    = tex1Dfetch(texref, atom1.x + cSim.stride2);
            int2 iatomJZ                    = tex1Dfetch(texref, atom1.y + cSim.stride2);  
            int2 iatomKZ                    = tex1Dfetch(texref, atom1.z + cSim.stride2);  
            PMEDouble atomIX                = __hiloint2double(iatomIX.y, iatomIX.x);
            PMEDouble atomJX                = __hiloint2double(iatomJX.y, iatomJX.x);
            PMEDouble atomKX                = __hiloint2double(iatomKX.y, iatomKX.x);
            PMEDouble xij                   = atomIX - atomJX; 
            PMEDouble xkj                   = atomKX - atomJX;           
            PMEDouble atomIY                = __hiloint2double(iatomIY.y, iatomIY.x);
            PMEDouble atomJY                = __hiloint2double(iatomJY.y, iatomJY.x);
            PMEDouble atomKY                = __hiloint2double(iatomKY.y, iatomKY.x);
            PMEDouble yij                   = atomIY - atomJY;
            PMEDouble ykj                   = atomKY - atomJY;            
            PMEDouble atomIZ                = __hiloint2double(iatomIZ.y, iatomIZ.x);
            PMEDouble atomJZ                = __hiloint2double(iatomJZ.y, iatomJZ.x);  
            PMEDouble atomKZ                = __hiloint2double(iatomKZ.y, iatomKZ.x);                            
#endif
            PMEDouble zij                   = atomIZ - atomJZ;
            PMEDouble zkj                   = atomKZ - atomJZ;
            PMEDouble rij2                  = xij * xij + yij * yij + zij * zij;
            PMEDouble rkj2                  = xkj * xkj + ykj * ykj + zkj * zkj;
            PMEDouble rij                   = sqrt(rij2);
            PMEDouble rkj                   = sqrt(rkj2);
            PMEDouble rdenom                = rij * rkj;
            PMEDouble cst                   = min(pt999, max(-pt999, (xij * xkj + yij * ykj + zij * zkj) / rdenom));
            PMEDouble theta                 = acos(cst);

            PMEDouble df;
#ifdef NMR_ENERGY
            PMEDouble e;
#endif
            if (theta < R1R2.x)
            {
                PMEDouble dif               = R1R2.x - R1R2.y;
                df                          = 2.0 * K2K3.x * dif;
#ifdef NMR_ENERGY
                e                           = df * (theta - R1R2.x) + K2K3.x * dif * dif;
#endif
            }
            else if (theta < R1R2.y)
            {
                PMEDouble dif               = theta - R1R2.y;
                df                          = 2.0 * K2K3.x * dif;
#ifdef NMR_ENERGY
                e                           = K2K3.x * dif * dif;
#endif
            }
            else if (theta < R3R4.x)
            {
                df                          = 0.0;
#ifdef NMR_ENERGY
                e                           = 0.0;
#endif
            }
            else if (theta < R3R4.y)
            {
                PMEDouble dif               = theta - R3R4.x;
                df                          = 2.0 * K2K3.y * dif;
#ifdef NMR_ENERGY
                e                           = K2K3.y * dif * dif;
#endif
            }
            else
            {
                PMEDouble dif               = R3R4.y - R3R4.x;
                df                          = 2.0 * K2K3.y * dif;
#ifdef NMR_ENERGY
                e                           = df * (theta - R3R4.y) + K2K3.y * dif * dif;
#endif
            }

            if (cSim.bJar)
            {
                double fold                 = cSim.pNMRJarData[2]; 
                double work                 = cSim.pNMRJarData[3];
                double first                = cSim.pNMRJarData[4];
                double fcurr                = -2.0 * K2K3.x * (theta - R1R2.y);
                if (first == 0.0) {
                    fold                    = -fcurr;
                    cSim.pNMRJarData[4]     = 1.0;
                }
                work                       += 0.5 * (fcurr + fold) * cSim.drjar;
                cSim.pNMRJarData[0]         = R1R2.y;
                cSim.pNMRJarData[1]         = theta;
                cSim.pNMRJarData[2]         = fcurr;
                cSim.pNMRJarData[3]         = work;
            }


#ifdef NMR_ENERGY
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom))  
#endif
                sE[threadIdx.x].angle      += lliroundd(ENERGYSCALE * e);
#endif

            PMEDouble snt                   = sin(theta);
            if (abs(snt) < 1.0e-14)
                snt                         = 1.0e-14;
            PMEDouble st                    = -df / snt;
            PMEDouble sth                   = st * cst;
            PMEDouble cik                   = st / rdenom;
            PMEDouble cii                   = sth / rij2;
            PMEDouble ckk                   = sth / rkj2;
            PMEDouble fx1                   = cii * xij - cik * xkj;
            PMEDouble fy1                   = cii * yij - cik * ykj;
            PMEDouble fz1                   = cii * zij - cik * zkj;
            PMEDouble fx2                   = ckk * xkj - cik * xij;
            PMEDouble fy2                   = ckk * ykj - cik * yij;
            PMEDouble fz2                   = ckk * zkj - cik * zij;   

#if defined(NMR_VIRIAL) || defined(use_SPFP)
            PMEAccumulator ifx1             = lliroundd(fx1 * FORCESCALE);
            PMEAccumulator ify1             = lliroundd(fy1 * FORCESCALE);
            PMEAccumulator ifz1             = lliroundd(fz1 * FORCESCALE);
            PMEAccumulator ifx2             = lliroundd(fx2 * FORCESCALE);
            PMEAccumulator ify2             = lliroundd(fy2 * FORCESCALE);
            PMEAccumulator ifz2             = lliroundd(fz2 * FORCESCALE);  
#endif    
                     
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom))
#endif
            {
#if defined(NMR_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.x], llitoulli(ifx1)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.x], llitoulli(ify1));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.x], llitoulli(ifz1));    
#elif defined(NMR_NEIGHBORLIST)            
                cSim.pForceXBuffer[atom1.w] = fx1;
                cSim.pForceYBuffer[atom1.w] = fy1;
                cSim.pForceZBuffer[atom1.w] = fz1;
#else
                cSim.pForceXBuffer[atom1.w]+= fx1;
                cSim.pForceYBuffer[atom1.w]+= fy1;
                cSim.pForceZBuffer[atom1.w]+= fz1;
#endif
            }
#ifdef MPI
            if ((atom1.z >= cSim.minLocalAtom) && (atom1.z < cSim.maxLocalAtom)) 
#endif
            {
#if defined(NMR_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.z], llitoulli(ifx2)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.z], llitoulli(ify2));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.z], llitoulli(ifz2));  
#elif defined (NMR_NEIGHBORLIST)                   
                cSim.pForceXBuffer[atom2.y] = fx2;
                cSim.pForceYBuffer[atom2.y] = fy2;
                cSim.pForceZBuffer[atom2.y] = fz2;
#else
                cSim.pForceXBuffer[atom2.y]+= fx2;
                cSim.pForceYBuffer[atom2.y]+= fy2;
                cSim.pForceZBuffer[atom2.y]+= fz2;
#endif
            }
#ifdef MPI
            if ((atom1.y >= cSim.minLocalAtom) && (atom1.y < cSim.maxLocalAtom)) 
#endif
            {
#if defined(NMR_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.y], llitoulli(-ifx1 - ifx2)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.y], llitoulli(-ify1 - ify2));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.y], llitoulli(-ifz1 - ifz2));                   
#elif defined(NMR_NEIGHBORLIST)              
                cSim.pForceXBuffer[atom2.x] = -fx1 - fx2;
                cSim.pForceYBuffer[atom2.x] = -fy1 - fy2;
                cSim.pForceZBuffer[atom2.x] = -fz1 - fz2;
#else
                cSim.pForceXBuffer[atom2.x]+= -fx1 - fx2;
                cSim.pForceYBuffer[atom2.x]+= -fy1 - fy2;
                cSim.pForceZBuffer[atom2.x]+= -fz1 - fz2;
#endif
            }
        }
    }

    // Calculate torsion restraints
    else if (pos < cSim.NMRTorsionOffset)
    {
        pos                                -= cSim.NMRAngleOffset;
        if (pos < cSim.NMRTorsions)
        {

#ifdef NMR_NEIGHBORLIST
            int4 atom1                      = cSim.pImageNMRTorsionID1[pos];
#if !defined(NMR_VIRIAL) && !defined(use_SPFP)            
            int4 atom2                      = cSim.pImageNMRTorsionID2[pos];
#endif
#else
            int4 atom1                      = cSim.pNMRTorsionID1[pos];
#ifndef use_SPFP
            int4 atom2                      = cSim.pNMRTorsionID2[pos];
#endif
#endif
            PMEDouble2 R1R2;
            PMEDouble2 R3R4;
            PMEDouble2 K2K3;
#ifdef NMR_TIMED
            int2 Step                       = cSim.pNMRTorsionStep[pos];
            int Inc                         = cSim.pNMRTorsionInc[pos];

            // Skip restraint if not active
            if ((step < Step.x) || ((step > Step.y) && (Step.y > 0)))
                goto exit;

            // Read timed data
            PMEDouble2 R1R2Slp              = cSim.pNMRTorsionR1R2Slp[pos];
            PMEDouble2 R1R2Int              = cSim.pNMRTorsionR1R2Int[pos];
            PMEDouble2 R3R4Slp              = cSim.pNMRTorsionR3R4Slp[pos];
            PMEDouble2 R3R4Int              = cSim.pNMRTorsionR3R4Int[pos];
            PMEDouble2 K2K3Slp              = cSim.pNMRTorsionK2K3Slp[pos];
            PMEDouble2 K2K3Int              = cSim.pNMRTorsionK2K3Int[pos];
            
            // Calculate increment
            double dstep                    = step - (double)((step - Step.x) % abs(Inc));

            // Calculate restraint values
            R1R2.x                          = R1R2Slp.x * dstep + R1R2Int.x;
            R1R2.y                          = R1R2Slp.y * dstep + R1R2Int.y;
            R3R4.x                          = R3R4Slp.x * dstep + R3R4Int.x;
            R3R4.y                          = R3R4Slp.y * dstep + R3R4Int.y;
            if (Inc > 0)
            {
                K2K3.x                      = K2K3Slp.x * dstep + K2K3Int.x;
                K2K3.y                      = K2K3Slp.y * dstep + K2K3Int.y;
            }
            else
            {
                int nstepu                  = (step - Step.x) / abs(Inc);
                K2K3.x                      = K2K3Int.x * pow(K2K3Slp.x, nstepu);
                K2K3.y                      = K2K3Int.y * pow(K2K3Slp.y, nstepu);
            }
#else
            R1R2                            = cSim.pNMRTorsionR1R2[pos];
            R3R4                            = cSim.pNMRTorsionR3R4[pos];
            K2K3                            = cSim.pNMRTorsionK2K3[pos];
#endif

#if defined(NODPTEXTURE)
#ifdef NMR_NEIGHBORLIST
            PMEDouble atomIX                = cSim.pImageX[atom1.x];
            PMEDouble atomJX                = cSim.pImageX[atom1.y];    
            PMEDouble atomIY                = cSim.pImageY[atom1.x];
            PMEDouble atomJY                = cSim.pImageY[atom1.y];    
            PMEDouble atomIZ                = cSim.pImageZ[atom1.x];
            PMEDouble atomJZ                = cSim.pImageZ[atom1.y];    
            PMEDouble atomKX                = cSim.pImageX[atom1.z];
            PMEDouble atomKY                = cSim.pImageY[atom1.z];
            PMEDouble atomKZ                = cSim.pImageZ[atom1.z];
            PMEDouble atomLX                = cSim.pImageX[atom1.w];
            PMEDouble atomLY                = cSim.pImageY[atom1.w];
            PMEDouble atomLZ                = cSim.pImageZ[atom1.w];
#else
            PMEDouble atomIX                = cSim.pAtomX[atom1.x];
            PMEDouble atomJX                = cSim.pAtomX[atom1.y];    
            PMEDouble atomIY                = cSim.pAtomY[atom1.x];
            PMEDouble atomJY                = cSim.pAtomY[atom1.y];    
            PMEDouble atomIZ                = cSim.pAtomZ[atom1.x];
            PMEDouble atomJZ                = cSim.pAtomZ[atom1.y];    
            PMEDouble atomKX                = cSim.pAtomX[atom1.z];
            PMEDouble atomKY                = cSim.pAtomY[atom1.z];
            PMEDouble atomKZ                = cSim.pAtomZ[atom1.z];
            PMEDouble atomLX                = cSim.pAtomX[atom1.w];
            PMEDouble atomLY                = cSim.pAtomY[atom1.w];
            PMEDouble atomLZ                = cSim.pAtomZ[atom1.w];
#endif
            PMEDouble xij                   = atomIX - atomJX;
            PMEDouble yij                   = atomIY - atomJY;
            PMEDouble zij                   = atomIZ - atomJZ;
            PMEDouble xkj                   = atomKX - atomJX;
            PMEDouble ykj                   = atomKY - atomJY;
            PMEDouble zkj                   = atomKZ - atomJZ;
            PMEDouble ykl                   = atomKY - atomLY;
#else        
            int2 iatomIX                    = tex1Dfetch(texref, atom1.x);
            int2 iatomJX                    = tex1Dfetch(texref, atom1.y);
            int2 iatomKX                    = tex1Dfetch(texref, atom1.z);
            int2 iatomLX                    = tex1Dfetch(texref, atom1.w);
            int2 iatomIY                    = tex1Dfetch(texref, atom1.x + cSim.stride);
            int2 iatomJY                    = tex1Dfetch(texref, atom1.y + cSim.stride);
            int2 iatomKY                    = tex1Dfetch(texref, atom1.z + cSim.stride);
            int2 iatomLY                    = tex1Dfetch(texref, atom1.w + cSim.stride);
            int2 iatomIZ                    = tex1Dfetch(texref, atom1.x + cSim.stride2);
            int2 iatomJZ                    = tex1Dfetch(texref, atom1.y + cSim.stride2);  
            int2 iatomKZ                    = tex1Dfetch(texref, atom1.z + cSim.stride2);  
            int2 iatomLZ                    = tex1Dfetch(texref, atom1.w + cSim.stride2);
            PMEDouble atomIY                = __hiloint2double(iatomIY.y, iatomIY.x);
            PMEDouble atomJY                = __hiloint2double(iatomJY.y, iatomJY.x);
            PMEDouble yij                   = atomIY - atomJY;  
            PMEDouble atomKZ                = __hiloint2double(iatomKZ.y, iatomKZ.x); 
            PMEDouble atomJZ                = __hiloint2double(iatomJZ.y, iatomJZ.x);               
            PMEDouble zkj                   = atomKZ - atomJZ;  
            PMEDouble atomIZ                = __hiloint2double(iatomIZ.y, iatomIZ.x);
            PMEDouble zij                   = atomIZ - atomJZ;
            PMEDouble atomKX                = __hiloint2double(iatomKX.y, iatomKX.x); 
            PMEDouble atomJX                = __hiloint2double(iatomJX.y, iatomJX.x);   
            PMEDouble xkj                   = atomKX - atomJX;  
            PMEDouble atomIX                = __hiloint2double(iatomIX.y, iatomIX.x);
            PMEDouble xij                   = atomIX - atomJX; 
            PMEDouble atomKY                = __hiloint2double(iatomKY.y, iatomKY.x);       
            PMEDouble ykj                   = atomKY - atomJY; 
            PMEDouble atomLY                = __hiloint2double(iatomLY.y, iatomLY.x);         
            PMEDouble ykl                   = atomKY - atomLY; 
            PMEDouble atomLX                = __hiloint2double(iatomLX.y, iatomLX.x);                     
            PMEDouble atomLZ                = __hiloint2double(iatomLZ.y, iatomLZ.x); 
#endif            
            PMEDouble zkl                   = atomKZ - atomLZ;                                           
            PMEDouble xkl                   = atomKX - atomLX;

            // Calculate ij X jk AND kl X jk:
            PMEDouble dx                    = yij * zkj - zij * ykj;
            PMEDouble dy                    = zij * xkj - xij * zkj;
            PMEDouble dz                    = xij * ykj - yij * xkj;

            PMEDouble gx                    = zkj * ykl - ykj * zkl;
            PMEDouble gy                    = xkj * zkl - zkj * xkl;
            PMEDouble gz                    = ykj * xkl - xkj * ykl;

            // Calculate the magnitudes of above vectors, and their dot product:
            PMEDouble fxi                   = dx * dx + dy * dy + dz * dz + tm24;
            PMEDouble fyi                   = gx * gx + gy * gy + gz * gz + tm24;
            PMEDouble ct                    = dx * gx + dy * gy + dz * gz;

            // Branch if linear dihedral:
            PMEDouble z1                    = rsqrt(fxi);
            PMEDouble z2                    = rsqrt(fyi);
            PMEDouble z11                   = z1 * z1;
            PMEDouble z22                   = z2 * z2;
            PMEDouble z12                   = z1 * z2;
            ct                             *= z1 * z2;
            ct                              = max(-0.9999999999999, min(ct, 0.9999999999999));
            PMEDouble ap                    = acos(ct);

            PMEDouble s                     = xkj * (dz * gy - dy * gz) + ykj * (dx * gz - dz * gx) + zkj * (dy * gx - dx * gy);
            if (s < 0.0) 
                ap = -ap;
            ap                              = PI - ap;
            PMEDouble sphi, cphi;
            faster_sincos(ap, &sphi, &cphi);

            // Translate the value of the torsion (by +- n*360) to bring it as close as
            // possible to one of the two central "cutoff" points (r2,r3). Use this as
            // the value of the torsion in the following comparison.
            PMEDouble apmean = (R1R2.y + R3R4.x) * 0.5;
            if (ap - apmean > PI)
                ap                         -= 2.0 * PI;
            if (apmean - ap > PI)
                ap                         += 2.0 * PI;

            PMEDouble df;
#ifdef NMR_ENERGY
            PMEDouble e;
#endif
            if (ap < R1R2.x)
            {
                PMEDouble dif               = R1R2.x - R1R2.y;
                df                          = 2.0 * K2K3.x * dif;
#ifdef NMR_ENERGY
                e                           = df * (ap - R1R2.x) + K2K3.x * dif * dif;
#endif
            }
            else if (ap < R1R2.y)
            {
                PMEDouble dif               = ap - R1R2.y;
                df                          = 2.0 * K2K3.x * dif;
#ifdef NMR_ENERGY
                e                           = K2K3.x * dif * dif;
#endif
            }
            else if (ap < R3R4.x)
            {
                df                          = 0.0;
#ifdef NMR_ENERGY
                e                           = 0.0;
#endif
            }
            else if (ap < R3R4.y)
            {
                PMEDouble dif               = ap - R3R4.x;
                df                          = 2.0 * K2K3.y * dif;
#ifdef NMR_ENERGY
                e                           = K2K3.y * dif * dif;
#endif
            }
            else
            {
                PMEDouble dif               = R3R4.y - R3R4.x;
                df                          = 2.0 * K2K3.y * dif;
#ifdef NMR_ENERGY
                e                           = df * (ap - R3R4.y) + K2K3.y * dif * dif;
#endif
            }
            df                             *= -1.0 / sphi;

            if (cSim.bJar)
            {
                double fold                 = cSim.pNMRJarData[2]; 
                double work                 = cSim.pNMRJarData[3];
                double first                = cSim.pNMRJarData[4];
                double fcurr                = -2.0 * K2K3.x * (ap - R1R2.y);
                if (first == 0.0) {
                    fold                    = -fcurr;
                    cSim.pNMRJarData[4]     = 1.0;
                }
                work                       += 0.5 * (fcurr + fold) * cSim.drjar;
                cSim.pNMRJarData[0]         = R1R2.y;
                cSim.pNMRJarData[1]         = ap;
                cSim.pNMRJarData[2]         = fcurr;
                cSim.pNMRJarData[3]         = work;
            }


#ifdef NMR_ENERGY
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom))  
#endif
                sE[threadIdx.x].torsion    += lliroundd(ENERGYSCALE * e);
#endif




            // dc = First derivative of cos(phi) w/respect to the cartesian differences t.
            PMEDouble dcdx                  =-gx * z12 - cphi * dx * z11;
            PMEDouble dcdy                  =-gy * z12 - cphi * dy * z11;
            PMEDouble dcdz                  =-gz * z12 - cphi * dz * z11;
            PMEDouble dcgx                  = dx * z12 + cphi * gx * z22;
            PMEDouble dcgy                  = dy * z12 + cphi * gy * z22;
            PMEDouble dcgz                  = dz * z12 + cphi * gz * z22;        
            PMEDouble dr1                   = df * ( dcdz * ykj - dcdy * zkj);
            PMEDouble dr2                   = df * ( dcdx * zkj - dcdz * xkj);
            PMEDouble dr3                   = df * ( dcdy * xkj - dcdx * ykj);
            PMEDouble dr4                   = df * ( dcgz * ykj - dcgy * zkj);
            PMEDouble dr5                   = df * ( dcgx * zkj - dcgz * xkj);
            PMEDouble dr6                   = df * ( dcgy * xkj - dcgx * ykj);
            PMEDouble drx                   = df * (-dcdy * zij + dcdz * yij + dcgy * zkl -  dcgz * ykl);
            PMEDouble dry                   = df * ( dcdx * zij - dcdz * xij - dcgx * zkl +  dcgz * xkl);
            PMEDouble drz                   = df * (-dcdx * yij + dcdy * xij + dcgx * ykl -  dcgy * xkl);  
            
#if defined(NMR_VIRIAL) || defined(use_SPFP)
            PMEAccumulator idr1             = lliroundd(dr1 * FORCESCALE);
            PMEAccumulator idr2             = lliroundd(dr2 * FORCESCALE);
            PMEAccumulator idr3             = lliroundd(dr3 * FORCESCALE);
            PMEAccumulator idr4             = lliroundd(dr4 * FORCESCALE);
            PMEAccumulator idr5             = lliroundd(dr5 * FORCESCALE);
            PMEAccumulator idr6             = lliroundd(dr6 * FORCESCALE);
            PMEAccumulator idrx             = lliroundd(drx * FORCESCALE);
            PMEAccumulator idry             = lliroundd(dry * FORCESCALE);
            PMEAccumulator idrz             = lliroundd(drz * FORCESCALE);
#endif                

#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom))
#endif
            {
#if defined(NMR_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.x], llitoulli(-idr1)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.x], llitoulli(-idr2));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.x], llitoulli(-idr3));  
#elif defined(NMR_NEIGHBORLIST)            
                cSim.pForceXBuffer[atom2.x] = -dr1;
                cSim.pForceYBuffer[atom2.x] = -dr2;
                cSim.pForceZBuffer[atom2.x] = -dr3;
#else
                cSim.pForceXBuffer[atom2.x]+= -dr1;
                cSim.pForceYBuffer[atom2.x]+= -dr2;
                cSim.pForceZBuffer[atom2.x]+= -dr3;
#endif
            }
#ifdef MPI
            if ((atom1.y >= cSim.minLocalAtom) && (atom1.y < cSim.maxLocalAtom)) 
#endif
            {
#if defined(NMR_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.y], llitoulli(-idrx + idr1)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.y], llitoulli(-idry + idr2));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.y], llitoulli(-idrz + idr3));  
#elif defined(NMR_NEIGHBORLIST)   
                cSim.pForceXBuffer[atom2.y] = -drx + dr1;
                cSim.pForceYBuffer[atom2.y] = -dry + dr2;
                cSim.pForceZBuffer[atom2.y] = -drz + dr3;
#else
                cSim.pForceXBuffer[atom2.y]+= -drx + dr1;
                cSim.pForceYBuffer[atom2.y]+= -dry + dr2;
                cSim.pForceZBuffer[atom2.y]+= -drz + dr3;
#endif
            }
#ifdef MPI
            if ((atom1.z >= cSim.minLocalAtom) && (atom1.z < cSim.maxLocalAtom)) 
#endif
            {
#if defined(NMR_VIRIAL) || defined(use_SPFP)          
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.z], llitoulli(idrx + idr4)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.z], llitoulli(idry + idr5));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.z], llitoulli(idrz + idr6));  
#elif defined(NMR_NEIGHBORLIST)  
                cSim.pForceXBuffer[atom2.z] = +drx + dr4;
                cSim.pForceYBuffer[atom2.z] = +dry + dr5;
                cSim.pForceZBuffer[atom2.z] = +drz + dr6;
#else
                cSim.pForceXBuffer[atom2.z]+= +drx + dr4;
                cSim.pForceYBuffer[atom2.z]+= +dry + dr5;
                cSim.pForceZBuffer[atom2.z]+= +drz + dr6;
#endif
            }
#ifdef MPI
            if ((atom1.w >= cSim.minLocalAtom) && (atom1.w < cSim.maxLocalAtom))             
#endif
            {
#if defined(NMR_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.w], llitoulli(-idr4)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.w], llitoulli(-idr5));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.w], llitoulli(-idr6));  
#elif defined(NMR_NEIGHBORLIST)            
                cSim.pForceXBuffer[atom2.w] = -dr4;
                cSim.pForceYBuffer[atom2.w] = -dr5;
                cSim.pForceZBuffer[atom2.w] = -dr6;
#else
                cSim.pForceXBuffer[atom2.w]+= -dr4;
                cSim.pForceYBuffer[atom2.w]+= -dr5;
                cSim.pForceZBuffer[atom2.w]+= -dr6;
#endif
            }           
        }
    }

#ifdef NMR_TIMED
exit:
    {}
#endif

#ifdef NMR_ENERGY
    // Reduce energies if active 
    __syncthreads();
    unsigned int m                          = 1;
    while (m < blockDim.x)
    {
        int p                               = threadIdx.x + m;      
        long long int distance              = ((p < blockDim.x) ? sE[p].distance : 0);
        long long int angle                 = ((p < blockDim.x) ? sE[p].angle : 0);
        long long int torsion               = ((p < blockDim.x) ? sE[p].torsion : 0);
        __syncthreads();
        sE[threadIdx.x].distance           += distance;
        sE[threadIdx.x].angle              += angle;
        sE[threadIdx.x].torsion            += torsion;
        __syncthreads();
        m                                  *= 2;
    }
    if (threadIdx.x == 0)
    {
        if (sE[0].distance != 0)
            atomicAdd(cSim.pENMRDistance, llitoulli(sE[0].distance));
        if (sE[0].angle != 0)
            atomicAdd(cSim.pENMRAngle, llitoulli(sE[0].angle));
        if (sE[0].torsion != 0)
            atomicAdd(cSim.pENMRTorsion, llitoulli(sE[0].torsion));
    }
#endif
}
