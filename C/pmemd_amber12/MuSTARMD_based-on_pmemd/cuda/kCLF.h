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
struct Energy {
    PMEDouble bond;
    PMEDouble angle;
    PMEDouble dihedral;
    PMEDouble el14;
    PMEDouble nb14;
    PMEDouble restraint;
};

struct Virial {
    PMEDouble vir_11;
    PMEDouble vir_22;
    PMEDouble vir_33;
};

#ifdef LOCAL_ENERGY
#if (__CUDA_ARCH__ >= 200)
    __shared__ Energy sE[SM_2X_LOCALFORCES_THREADS_PER_BLOCK];
#else
    __shared__ Energy sE[SM_13_LOCALFORCES_THREADS_PER_BLOCK];
#endif
#endif
#ifdef LOCAL_VIRIAL
#if (__CUDA_ARCH__ >= 200)
    __shared__ Virial sV[SM_2X_LOCALFORCES_THREADS_PER_BLOCK];
#else
    __shared__ Virial sV[SM_13_LOCALFORCES_THREADS_PER_BLOCK];
#endif
#endif
    int pos                                 = blockIdx.x * blockDim.x + threadIdx.x;    
#ifdef LOCAL_ENERGY
    sE[threadIdx.x].bond                    = (PMEDouble)0.0;
    sE[threadIdx.x].angle                   = (PMEDouble)0.0;
    sE[threadIdx.x].dihedral                = (PMEDouble)0.0;
    sE[threadIdx.x].el14                    = (PMEDouble)0.0;
    sE[threadIdx.x].nb14                    = (PMEDouble)0.0;
    sE[threadIdx.x].restraint               = (PMEDouble)0.0;
#endif
#ifdef LOCAL_VIRIAL
    sV[threadIdx.x].vir_11                  = (PMEDouble)0.0;
    sV[threadIdx.x].vir_22                  = (PMEDouble)0.0;
    sV[threadIdx.x].vir_33                  = (PMEDouble)0.0;
#endif    

    // Calculate bond forces
    while (pos < cSim.bondOffset)
    {
        if (pos < cSim.bonds)
        {
#ifdef LOCAL_NEIGHBORLIST           
            int4 atom                       = cSim.pImageBondID[pos];
            PMEDouble2 bond                 = cSim.pBond[pos];       
#else
            int4 atom                       = cSim.pBondID[pos];
            PMEDouble2 bond                 = cSim.pBond[pos];     
#endif  
#ifdef use_SPSP            
            PMEDouble atomIX                = tex1Dfetch(texref, atom.x);
            PMEDouble atomJX                = tex1Dfetch(texref, atom.y);
            PMEDouble atomIY                = tex1Dfetch(texref, atom.x + cSim.stride);
            PMEDouble atomJY                = tex1Dfetch(texref, atom.y + cSim.stride);
            PMEDouble atomIZ                = tex1Dfetch(texref, atom.x + cSim.stride2);
            PMEDouble atomJZ                = tex1Dfetch(texref, atom.y + cSim.stride2);  
            PMEDouble xij                   = atomIX - atomJX;
            PMEDouble yij                   = atomIY - atomJY;            
#elif defined(NODPTEXTURE)
#ifdef LOCAL_NEIGHBORLIST
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

            // Calculate derivative
            PMEDouble rij                   = sqrt(xij * xij + yij * yij + zij * zij);
            PMEDouble da                    = rij - bond.y;
            PMEDouble df                    = bond.x * da;
            PMEDouble dfw                   = (df + df) / rij;
#ifdef LOCAL_ENERGY
#ifdef MPI
            if ((atom.x >= cSim.minLocalAtom) && (atom.x < cSim.maxLocalAtom))  
#endif    
                sE[threadIdx.x].bond       += df * da;
#endif
            
            // Calculate force
            PMEDouble fx                    = dfw * xij;
            PMEDouble fy                    = dfw * yij;
            PMEDouble fz                    = dfw * zij;    

#ifdef LOCAL_NEIGHBORLIST
#ifdef LOCAL_VIRIAL
            unsigned long long int ifx      = fabs(fx) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int ify      = fabs(fy) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int ifz      = fabs(fz) * FORCESCALE + (PMEDouble)0.5;
            if (fx < 0)
                ifx                         = 0ull - ifx;
            if (fy < 0)
                ify                         = 0ull - ify;                    
            if (fz < 0)
                ifz                         = 0ull - ifz;
#endif
#ifdef MPI
            if ((atom.y >= cSim.minLocalAtom) && (atom.y < cSim.maxLocalAtom))
#endif
            {
#ifdef LOCAL_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom.y], ifx); 
                atomicAdd(&cSim.pUllForceY[atom.y], ify);  
                atomicAdd(&cSim.pUllForceZ[atom.y], ifz);               
#else                
                cSim.pForceXBuffer[atom.w]  =  fx;
                cSim.pForceYBuffer[atom.w]  =  fy;
                cSim.pForceZBuffer[atom.w]  =  fz;
#endif
            }
#ifdef MPI
            if ((atom.x >= cSim.minLocalAtom) && (atom.x < cSim.maxLocalAtom))
#endif
            {     
#ifdef LOCAL_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom.x], 0ull - ifx); 
                atomicAdd(&cSim.pUllForceY[atom.x], 0ull - ify);  
                atomicAdd(&cSim.pUllForceZ[atom.x], 0ull - ifz);   
#else                    
                cSim.pForceXBuffer[atom.z]  = -fx;
                cSim.pForceYBuffer[atom.z]  = -fy;
                cSim.pForceZBuffer[atom.z]  = -fz;
#endif                
            }
#else        
#ifdef MPI
            if ((atom.y >= cSim.minLocalAtom) && (atom.y < cSim.maxLocalAtom))
#endif
            {                 
                cSim.pForceXBuffer[atom.w] += fx;
                cSim.pForceYBuffer[atom.w] += fy;
                cSim.pForceZBuffer[atom.w] += fz; 
            }
#ifdef MPI
            if ((atom.x >= cSim.minLocalAtom) && (atom.x < cSim.maxLocalAtom))
#endif            
            {           
                cSim.pForceXBuffer[atom.z] += -fx;
                cSim.pForceYBuffer[atom.z] += -fy;
                cSim.pForceZBuffer[atom.z] += -fz;
            }
#endif
        }
        pos                                += blockDim.x * gridDim.x;
    }
    
    // Calculate bond angle forces
    while (pos < cSim.bondAngleOffset)
    {
        pos                                -= cSim.bondOffset;
        if (pos < cSim.bondAngles)
        {
#ifdef LOCAL_NEIGHBORLIST        
            int4 atom1                      = cSim.pImageBondAngleID1[pos];
#else
            int4 atom1                      = cSim.pBondAngleID1[pos];
#endif
            PMEDouble2 bondAngle            = cSim.pBondAngle[pos];
#ifdef LOCAL_NEIGHBORLIST
#ifndef LOCAL_VIRIAL
            int2 atom2                      = cSim.pImageBondAngleID2[pos];               
#endif
#else
            int2 atom2                      = cSim.pBondAngleID2[pos];       
#endif        
#ifdef use_SPSP
            PMEDouble atomIX                = tex1Dfetch(texref, atom1.x);
            PMEDouble atomJX                = tex1Dfetch(texref, atom1.y);
            PMEDouble atomKX                = tex1Dfetch(texref, atom1.z);
            PMEDouble atomIY                = tex1Dfetch(texref, atom1.x + cSim.stride);
            PMEDouble atomJY                = tex1Dfetch(texref, atom1.y + cSim.stride);
            PMEDouble atomKY                = tex1Dfetch(texref, atom1.z + cSim.stride);
            PMEDouble atomIZ                = tex1Dfetch(texref, atom1.x + cSim.stride2);
            PMEDouble atomJZ                = tex1Dfetch(texref, atom1.y + cSim.stride2);  
            PMEDouble atomKZ                = tex1Dfetch(texref, atom1.z + cSim.stride2);
            PMEDouble xij                   = atomIX - atomJX; 
            PMEDouble xkj                   = atomKX - atomJX;           
            PMEDouble yij                   = atomIY - atomJY;
            PMEDouble ykj                   = atomKY - atomJY;       
#elif defined(NODPTEXTURE)
#ifdef LOCAL_NEIGHBORLIST
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
            PMEDouble rij                   = xij * xij + yij * yij + zij * zij;
            PMEDouble rkj                   = xkj * xkj + ykj * ykj + zkj * zkj;
            PMEDouble rik                   = sqrt(rij * rkj);
            PMEDouble cst                   = min(pt999, max(-pt999, (xij * xkj + yij * ykj + zij * zkj) / rik));
            PMEDouble ant                   = acos(cst);

            // Calculation of the energy and derivative
            PMEDouble da                    = ant - bondAngle.y;
            PMEDouble df                    = bondAngle.x * da;
            PMEDouble dfw                   = -(df + df) / sin(ant);
#ifdef LOCAL_ENERGY    
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom))        
#endif
                sE[threadIdx.x].angle      += df * da;
#endif

            // Calculation of the force
            PMEDouble cik                   = dfw / rik;
            PMEDouble sth                   = dfw * cst;
            PMEDouble cii                   = sth / rij;
            PMEDouble ckk                   = sth / rkj;
            PMEDouble fx1                   = cii * xij - cik * xkj; 
            PMEDouble fy1                   = cii * yij - cik * ykj;
            PMEDouble fz1                   = cii * zij - cik * zkj;
            PMEDouble fx2                   = ckk * xkj - cik * xij;
            PMEDouble fy2                   = ckk * ykj - cik * yij;
            PMEDouble fz2                   = ckk * zkj - cik * zij;            
#ifdef LOCAL_NEIGHBORLIST
#ifdef LOCAL_VIRIAL
             unsigned long long int ifx1    = fabs(fx1) * FORCESCALE + (PMEDouble)0.5;
             unsigned long long int ify1    = fabs(fy1) * FORCESCALE + (PMEDouble)0.5;
             unsigned long long int ifz1    = fabs(fz1) * FORCESCALE + (PMEDouble)0.5;
             if (fx1 < 0.0)
                 ifx1                       = 0ull - ifx1;
             if (fy1 < 0.0)
                 ify1                       = 0ull - ify1;                    
             if (fz1 < 0.0)
                 ifz1                       = 0ull - ifz1;
             unsigned long long int ifx2    = fabs(fx2) * FORCESCALE + (PMEDouble)0.5;
             unsigned long long int ify2    = fabs(fy2) * FORCESCALE + (PMEDouble)0.5;
             unsigned long long int ifz2    = fabs(fz2) * FORCESCALE + (PMEDouble)0.5;
             if (fx2 < 0.0)
                 ifx2                       = 0ull - ifx2;
             if (fy2 < 0.0)
                 ify2                       = 0ull - ify2;                    
             if (fz2 < 0.0)
                 ifz2                       = 0ull - ifz2;     
#endif                             
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom))
#endif
            {
#ifdef LOCAL_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom1.x], ifx1); 
                atomicAdd(&cSim.pUllForceY[atom1.x], ify1);  
                atomicAdd(&cSim.pUllForceZ[atom1.x], ifz1);    
#else            
                cSim.pForceXBuffer[atom1.w] = fx1;
                cSim.pForceYBuffer[atom1.w] = fy1;
                cSim.pForceZBuffer[atom1.w] = fz1;
#endif
            }
#ifdef MPI
            if ((atom1.z >= cSim.minLocalAtom) && (atom1.z < cSim.maxLocalAtom)) 
#endif
            {
#ifdef LOCAL_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom1.z], ifx2); 
                atomicAdd(&cSim.pUllForceY[atom1.z], ify2);  
                atomicAdd(&cSim.pUllForceZ[atom1.z], ifz2);  
#else                   
                cSim.pForceXBuffer[atom2.y] = fx2;
                cSim.pForceYBuffer[atom2.y] = fy2;
                cSim.pForceZBuffer[atom2.y] = fz2;
#endif
            }
#ifdef MPI
            if ((atom1.y >= cSim.minLocalAtom) && (atom1.y < cSim.maxLocalAtom)) 
#endif
            {
#ifdef LOCAL_VIRIAL
                ifx1                        = 0ull - ifx1 - ifx2;
                ify1                        = 0ull - ify1 - ify2;
                ifz1                        = 0ull - ifz1 - ifz2;
                atomicAdd(&cSim.pUllForceX[atom1.y], ifx1); 
                atomicAdd(&cSim.pUllForceY[atom1.y], ify1);  
                atomicAdd(&cSim.pUllForceZ[atom1.y], ifz1);                   
#else            
                fx1                         = -fx1 - fx2;
                fy1                         = -fy1 - fy2;
                fz1                         = -fz1 - fz2;            
                cSim.pForceXBuffer[atom2.x] = fx1;
                cSim.pForceYBuffer[atom2.x] = fy1;
                cSim.pForceZBuffer[atom2.x] = fz1;
#endif
            }                
#else
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom))
#endif
            {
                cSim.pForceXBuffer[atom1.w]+= fx1;
                cSim.pForceYBuffer[atom1.w]+= fy1;
                cSim.pForceZBuffer[atom1.w]+= fz1;
            }
#ifdef MPI
            if ((atom1.z >= cSim.minLocalAtom) && (atom1.z < cSim.maxLocalAtom)) 
#endif            
            {
                cSim.pForceXBuffer[atom2.y]+= fx2;
                cSim.pForceYBuffer[atom2.y]+= fy2;
                cSim.pForceZBuffer[atom2.y]+= fz2;
            }
#ifdef MPI
            if ((atom1.y >= cSim.minLocalAtom) && (atom1.y < cSim.maxLocalAtom)) 
#endif            
            {
                fx1                         = -fx1 - fx2;
                fy1                         = -fy1 - fy2;
                fz1                         = -fz1 - fz2;
                cSim.pForceXBuffer[atom2.x]+= fx1;
                cSim.pForceYBuffer[atom2.x]+= fy1;
                cSim.pForceZBuffer[atom2.x]+= fz1;  
            }    
#endif
        }
        pos                                += cSim.bondOffset + blockDim.x * gridDim.x;
    }
    
    // Calculate dihedral forces
    while (pos < cSim.dihedralOffset)
    {
        pos                                -= cSim.bondAngleOffset;
        if (pos < cSim.dihedrals)
        {
#ifdef LOCAL_NEIGHBORLIST
            int4 atom1                      = cSim.pImageDihedralID1[pos];
#ifndef LOCAL_VIRIAL            
            int4 atom2                      = cSim.pImageDihedralID2[pos];
#endif
#else
            int4 atom1                      = cSim.pDihedralID1[pos];
            int4 atom2                      = cSim.pDihedralID2[pos];
#endif
#ifdef use_SPSP    
            PMEDouble atomIX                = tex1Dfetch(texref, atom1.x);
            PMEDouble atomJX                = tex1Dfetch(texref, atom1.y);
            PMEDouble atomKX                = tex1Dfetch(texref, atom1.z);
            PMEDouble atomLX                = tex1Dfetch(texref, atom1.w);
            PMEDouble atomIY                = tex1Dfetch(texref, atom1.x + cSim.stride);
            PMEDouble atomJY                = tex1Dfetch(texref, atom1.y + cSim.stride);
            PMEDouble atomKY                = tex1Dfetch(texref, atom1.z + cSim.stride);
            PMEDouble atomLY                = tex1Dfetch(texref, atom1.w + cSim.stride);
            PMEDouble atomIZ                = tex1Dfetch(texref, atom1.x + cSim.stride2);
            PMEDouble atomJZ                = tex1Dfetch(texref, atom1.y + cSim.stride2);  
            PMEDouble atomKZ                = tex1Dfetch(texref, atom1.z + cSim.stride2);  
            PMEDouble atomLZ                = tex1Dfetch(texref, atom1.w + cSim.stride2);
            PMEDouble2 dihedral1            = cSim.pDihedral1[pos];
            PMEDouble2 dihedral2            = cSim.pDihedral2[pos];
            PMEDouble dihedral3             = cSim.pDihedral3[pos];
            PMEDouble yij                   = atomIY - atomJY;      
            PMEDouble zkj                   = atomKZ - atomJZ;  
            PMEDouble zij                   = atomIZ - atomJZ;
            PMEDouble xkj                   = atomKX - atomJX;  
            PMEDouble xij                   = atomIX - atomJX; 
            PMEDouble ykj                   = atomKY - atomJY;         
            PMEDouble ykl                   = atomKY - atomLY; 
#elif defined(NODPTEXTURE)
#ifdef LOCAL_NEIGHBORLIST
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
            PMEDouble2 dihedral1            = cSim.pDihedral1[pos];
            PMEDouble2 dihedral2            = cSim.pDihedral2[pos];
            PMEDouble dihedral3             = cSim.pDihedral3[pos];
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
            PMEDouble2 dihedral1            = cSim.pDihedral1[pos];
            PMEDouble2 dihedral2            = cSim.pDihedral2[pos];
            PMEDouble dihedral3             = cSim.pDihedral3[pos];
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

            // Get the normal vector
            PMEDouble dx                    = yij * zkj - zij * ykj;
            PMEDouble dy                    = zij * xkj - xij * zkj;
            PMEDouble dz                    = xij * ykj - yij * xkj;
            PMEDouble gx                    = zkj * ykl - ykj * zkl;
            PMEDouble gy                    = xkj * zkl - zkj * xkl;
            PMEDouble gz                    = ykj * xkl - xkj * ykl;
            PMEDouble fxi                   = sqrt(dx * dx + dy * dy + dz * dz + tm24);
            PMEDouble fyi                   = sqrt(gx * gx + gy * gy + gz * gz + tm24);
            PMEDouble ct                    = dx * gx + dy * gy + dz * gz;

            // Branch if linear dihedral:
            PMEDouble z1                    = (tenm3 <= fxi) ? (one / fxi) : zero;
            PMEDouble z2                    = (tenm3 <= fyi) ? (one / fyi) : zero;
            PMEDouble z12                   = z1 * z2;
            PMEDouble fzi                   = (z12 != zero) ? one : zero;
            PMEDouble s                     = xkj * (dz * gy - dy * gz) + ykj * (dx * gz - dz * gx) + zkj * (dy * gx - dx * gy);
            PMEDouble ap                    = PI - abs(acos(max(-one, min(one, ct * z12)))) * (s >= (PMEDouble)0.0 ? (PMEDouble)1.0 : (PMEDouble)-1.0);
            PMEDouble sphi, cphi;
#ifndef use_SPSP            
            faster_sincos(ap, &sphi, &cphi);
#else
            sincos(ap, &sphi, &cphi);
#endif

            // Calculate the energy and the derivatives with respect to cosphi
            PMEDouble ct0                   = dihedral1.y * ap;
            PMEDouble sinnp, cosnp;
#ifndef use_SPSP
            faster_sincos(ct0, &sinnp, &cosnp);
#else
            sincos(ct0, &sinnp, &cosnp);
#endif
#ifdef LOCAL_ENERGY
            PMEDouble epw                   = (dihedral2.x + cosnp * dihedral2.y + sinnp * dihedral3) * fzi;
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom))        
#endif            
                sE[threadIdx.x].dihedral   += epw;
#endif            
            PMEDouble dums                  = sphi + tm24 * (sphi >= (PMEDouble)0.0 ? (PMEDouble)1.0 : (PMEDouble)-1.0);
            PMEDouble df;
            if (tm06 > abs(dums))
                df                          = fzi * dihedral2.y * (dihedral1.y - dihedral1.x + dihedral1.x * cphi);
            else
                df                          = fzi * dihedral1.y * (dihedral2.y * sinnp - dihedral3 * cosnp) / dums;

	    // AMD If iamd=2,3 a dihedral boost will be added to the dihedral forces, therefore  df=df*fwgtd
	    if((cSim.iamd==2)||(cSim.iamd==3)){
	      PMEDouble fwgtd;
	      fwgtd   = cSim.pAMDfwgtd[0];
	      df *= fwgtd;
	    }
            // Now do torsional first derivatives:

            // Now, set up array dc = 1st der. of cosphi w/respect to cartesian differences:
            PMEDouble z11                   = z1 * z1;
            z12                             = z1 * z2;
            PMEDouble z22                   = z2 * z2;
            PMEDouble dc1                   = -gx * z12 - cphi * dx * z11;
            PMEDouble dc2                   = -gy * z12 - cphi * dy * z11;
            PMEDouble dc3                   = -gz * z12 - cphi * dz * z11;
            PMEDouble dc4                   =  dx * z12 + cphi * gx * z22;
            PMEDouble dc5                   =  dy * z12 + cphi * gy * z22;
            PMEDouble dc6                   =  dz * z12 + cphi * gz * z22;

            // Update the first derivative array:
            PMEDouble dr1                   = df * ( dc3 * ykj - dc2 * zkj);
            PMEDouble dr2                   = df * ( dc1 * zkj - dc3 * xkj);
            PMEDouble dr3                   = df * ( dc2 * xkj - dc1 * ykj);
            PMEDouble dr4                   = df * ( dc6 * ykj - dc5 * zkj);
            PMEDouble dr5                   = df * ( dc4 * zkj - dc6 * xkj);
            PMEDouble dr6                   = df * ( dc5 * xkj - dc4 * ykj);
            PMEDouble drx                   = df * (-dc2 * zij + dc3 * yij + dc5 * zkl - dc6 * ykl);
            PMEDouble dry                   = df * ( dc1 * zij - dc3 * xij - dc4 * zkl + dc6 * xkl);
            PMEDouble drz                   = df * (-dc1 * yij + dc2 * xij + dc4 * ykl - dc5 * xkl);
#ifdef LOCAL_NEIGHBORLIST
#ifdef LOCAL_VIRIAL
            unsigned long long int idr1     = fabs(dr1) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int idr2     = fabs(dr2) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int idr3     = fabs(dr3) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int idr4     = fabs(dr4) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int idr5     = fabs(dr5) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int idr6     = fabs(dr6) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int idrx     = fabs(drx) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int idry     = fabs(dry) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int idrz     = fabs(drz) * FORCESCALE + (PMEDouble)0.5;
            if (dr1 < (PMEDouble)0.0)
                idr1                        = 0ull - idr1;
            if (dr2 < (PMEDouble)0.0)
                idr2                        = 0ull - idr2;                
            if (dr3 < (PMEDouble)0.0)
                idr3                        = 0ull - idr3; 
            if (dr4 < (PMEDouble)0.0)
                idr4                        = 0ull - idr4; 
            if (dr5 < (PMEDouble)0.0)
                idr5                        = 0ull - idr5; 
            if (dr6 < (PMEDouble)0.0)
                idr6                        = 0ull - idr6; 
            if (drx < (PMEDouble)0.0)
                idrx                        = 0ull - idrx;         
            if (dry < (PMEDouble)0.0)
                idry                        = 0ull - idry;         
            if (drz < (PMEDouble)0.0)
                idrz                        = 0ull - idrz;     
#endif                
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom))
#endif
            {
#ifdef LOCAL_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom1.x], 0ull - idr1); 
                atomicAdd(&cSim.pUllForceY[atom1.x], 0ull - idr2);  
                atomicAdd(&cSim.pUllForceZ[atom1.x], 0ull - idr3);  
#else            
                cSim.pForceXBuffer[atom2.x] = -dr1;
                cSim.pForceYBuffer[atom2.x] = -dr2;
                cSim.pForceZBuffer[atom2.x] = -dr3;
#endif
            }
#ifdef MPI
            if ((atom1.y >= cSim.minLocalAtom) && (atom1.y < cSim.maxLocalAtom)) 
#endif
            {
#ifdef LOCAL_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom1.y], 0ull - idrx + idr1); 
                atomicAdd(&cSim.pUllForceY[atom1.y], 0ull - idry + idr2);  
                atomicAdd(&cSim.pUllForceZ[atom1.y], 0ull - idrz + idr3);  
#else   
                cSim.pForceXBuffer[atom2.y] = -drx + dr1;
                cSim.pForceYBuffer[atom2.y] = -dry + dr2;
                cSim.pForceZBuffer[atom2.y] = -drz + dr3;
#endif
            }
#ifdef MPI
            if ((atom1.z >= cSim.minLocalAtom) && (atom1.z < cSim.maxLocalAtom)) 
#endif
            {
#ifdef LOCAL_VIRIAL          
                atomicAdd(&cSim.pUllForceX[atom1.z], idrx + idr4); 
                atomicAdd(&cSim.pUllForceY[atom1.z], idry + idr5);  
                atomicAdd(&cSim.pUllForceZ[atom1.z], idrz + idr6);  
#else  
                cSim.pForceXBuffer[atom2.z] = +drx + dr4;
                cSim.pForceYBuffer[atom2.z] = +dry + dr5;
                cSim.pForceZBuffer[atom2.z] = +drz + dr6;
#endif
            }
#ifdef MPI
            if ((atom1.w >= cSim.minLocalAtom) && (atom1.w < cSim.maxLocalAtom))             
#endif
            {
#ifdef LOCAL_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom1.w], 0ull - idr4); 
                atomicAdd(&cSim.pUllForceY[atom1.w], 0ull - idr5);  
                atomicAdd(&cSim.pUllForceZ[atom1.w], 0ull - idr6);  
#else            
                cSim.pForceXBuffer[atom2.w] = -dr4;
                cSim.pForceYBuffer[atom2.w] = -dr5;
                cSim.pForceZBuffer[atom2.w] = -dr6; 
#endif
            }           
#else
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom))
#endif
            {
                cSim.pForceXBuffer[atom2.x]+= -dr1;
                cSim.pForceYBuffer[atom2.x]+= -dr2;
                cSim.pForceZBuffer[atom2.x]+= -dr3;
            }
#ifdef MPI
            if ((atom1.y >= cSim.minLocalAtom) && (atom1.y < cSim.maxLocalAtom)) 
#endif                
            {
                cSim.pForceXBuffer[atom2.y]+= -drx + dr1;
                cSim.pForceYBuffer[atom2.y]+= -dry + dr2;
                cSim.pForceZBuffer[atom2.y]+= -drz + dr3;
            }
#ifdef MPI
            if ((atom1.z >= cSim.minLocalAtom) && (atom1.z < cSim.maxLocalAtom)) 
#endif             
            {   
                cSim.pForceXBuffer[atom2.z]+= +drx + dr4;
                cSim.pForceYBuffer[atom2.z]+= +dry + dr5;
                cSim.pForceZBuffer[atom2.z]+= +drz + dr6;
            }      
#ifdef MPI
            if ((atom1.w >= cSim.minLocalAtom) && (atom1.w < cSim.maxLocalAtom))             
#endif            
            {   
                cSim.pForceXBuffer[atom2.w]+= -dr4;
                cSim.pForceYBuffer[atom2.w]+= -dr5;
                cSim.pForceZBuffer[atom2.w]+= -dr6;
            }
#endif
        }
        pos                                += cSim.bondAngleOffset + blockDim.x * gridDim.x;
    }
   
   
    // Calculate 1-4 forces
    while (pos < cSim.nb14Offset)
    {
        pos                                -= cSim.dihedralOffset;
        if (pos < cSim.nb14s)
        {
#ifdef LOCAL_NEIGHBORLIST
            int4 atom                       = cSim.pImageNb14ID[pos];
#else
            int4 atom                       = cSim.pNb14ID[pos];
#endif
            PMEDouble2 nb141                = cSim.pNb141[pos];
            PMEDouble2 nb142                = cSim.pNb142[pos];
#ifdef use_SPSP    
            PMEDouble atomIX                = tex1Dfetch(texref, atom.x);
            PMEDouble atomJX                = tex1Dfetch(texref, atom.y);
            PMEDouble atomIY                = tex1Dfetch(texref, atom.x + cSim.stride);
            PMEDouble atomJY                = tex1Dfetch(texref, atom.y + cSim.stride);
            PMEDouble atomIZ                = tex1Dfetch(texref, atom.x + cSim.stride2);
            PMEDouble atomJZ                = tex1Dfetch(texref, atom.y + cSim.stride2);   
            PMEDouble dx                    = atomJX - atomIX;
            PMEDouble dy                    = atomJY - atomIY;            
#elif defined(NODPTEXTURE)
#ifdef LOCAL_NEIGHBORLIST
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
            PMEDouble dx                    = atomJX - atomIX;
            PMEDouble dy                    = atomJY - atomIY;
#else
            int2 iatomIX                    = tex1Dfetch(texref, atom.x);
            int2 iatomJX                    = tex1Dfetch(texref, atom.y);
            int2 iatomIY                    = tex1Dfetch(texref, atom.x + cSim.stride);
            int2 iatomJY                    = tex1Dfetch(texref, atom.y + cSim.stride);
            int2 iatomIZ                    = tex1Dfetch(texref, atom.x + cSim.stride2);
            int2 iatomJZ                    = tex1Dfetch(texref, atom.y + cSim.stride2);   
            PMEDouble atomIX                = __hiloint2double(iatomIX.y, iatomIX.x);
            PMEDouble atomJX                = __hiloint2double(iatomJX.y, iatomJX.x);
            PMEDouble dx                    = atomJX - atomIX;
            PMEDouble atomIY                = __hiloint2double(iatomIY.y, iatomIY.x);
            PMEDouble atomJY                = __hiloint2double(iatomJY.y, iatomJY.x);
            PMEDouble dy                    = atomJY - atomIY;            
            PMEDouble atomIZ                = __hiloint2double(iatomIZ.y, iatomIZ.x);
            PMEDouble atomJZ                = __hiloint2double(iatomJZ.y, iatomJZ.x);
#endif            
            PMEDouble dz                    = atomJZ - atomIZ;           
            PMEDouble r2                    = dx * dx + dy * dy + dz * dz;
            PMEDouble rinv                  = rsqrt(r2);
            PMEDouble r2inv                 = rinv * rinv;
            PMEDouble r6                    = r2inv * r2inv * r2inv;
            PMEDouble g                     = rinv * nb141.x;          
            PMEDouble f6                    = nb142.y * r6;
            PMEDouble f12                   = nb142.x * (r6 * r6);
            PMEDouble df                    = (g + nb141.y * ((PMEDouble)12.0 * f12 - (PMEDouble)6.0 * f6)) * r2inv;
#ifdef LOCAL_ENERGY       
#ifdef MPI
            if ((atom.x >= cSim.minLocalAtom) && (atom.x < cSim.maxLocalAtom))        
#endif     
            {
                sE[threadIdx.x].el14       += g;
                sE[threadIdx.x].nb14       += (f12 - f6) * nb141.y;
            }
#endif
            PMEDouble fx                    = df * dx;
            PMEDouble fy                    = df * dy;
            PMEDouble fz                    = df * dz;
#ifdef LOCAL_VIRIAL
#ifdef MPI
            if ((atom.x >= cSim.minLocalAtom) && (atom.x < cSim.maxLocalAtom))        
#endif     
            {
                sV[threadIdx.x].vir_11     -= fx * dx;
                sV[threadIdx.x].vir_22     -= fy * dy;
                sV[threadIdx.x].vir_33     -= fz * dz;
            }
#endif
#ifdef LOCAL_NEIGHBORLIST           
#ifdef MPI
            if ((atom.y >= cSim.minLocalAtom) && (atom.y < cSim.maxLocalAtom))
#endif
            {
                cSim.pForceXBuffer[atom.w]  =  fx;
                cSim.pForceYBuffer[atom.w]  =  fy;
                cSim.pForceZBuffer[atom.w]  =  fz;
            }
#ifdef MPI
            if ((atom.x >= cSim.minLocalAtom) && (atom.x < cSim.maxLocalAtom))
#endif
            {           
                cSim.pForceXBuffer[atom.z]  = -fx;
                cSim.pForceYBuffer[atom.z]  = -fy;
                cSim.pForceZBuffer[atom.z]  = -fz;
            }
#else
#ifdef MPI
            if ((atom.y >= cSim.minLocalAtom) && (atom.y < cSim.maxLocalAtom))
#endif
            {
                cSim.pForceXBuffer[atom.w] +=  fx;
                cSim.pForceYBuffer[atom.w] +=  fy;
                cSim.pForceZBuffer[atom.w] +=  fz;
            }    
#ifdef MPI
            if ((atom.x >= cSim.minLocalAtom) && (atom.x < cSim.maxLocalAtom))
#endif            
            {
                cSim.pForceXBuffer[atom.z] += -fx;
                cSim.pForceYBuffer[atom.z] += -fy;
                cSim.pForceZBuffer[atom.z] += -fz;
            }
#endif
        }
        pos                                += cSim.dihedralOffset + blockDim.x * gridDim.x;
    }    
    
    // Calculate Constraint forces
    while (pos < cSim.constraintOffset)
    {
        pos                                -= cSim.nb14Offset;
        if (pos < cSim.constraints)
        {
#ifdef LOCAL_NEIGHBORLIST
                int2 atom                   = cSim.pImageConstraintID[pos];           
#else
                int2 atom                   = cSim.pConstraintID[pos];
#endif
#ifdef MPI
            if ((atom.x >= cSim.minLocalAtom) && (atom.x < cSim.maxLocalAtom))
            {
#endif 
                PMEDouble2 constraint1      = cSim.pConstraint1[pos];
                PMEDouble2 constraint2      = cSim.pConstraint2[pos];
#ifdef use_SPSP            
                PMEDouble atomX             = tex1Dfetch(texref, atom.x);          
                PMEDouble atomY             = tex1Dfetch(texref, atom.x + cSim.stride);
                PMEDouble atomZ             = tex1Dfetch(texref, atom.x + cSim.stride2);
                PMEDouble ax                = atomX - constraint1.y;
                PMEDouble ay                = atomY - constraint2.x;            
#elif defined(NODPTEXTURE)
#ifdef LOCAL_NEIGHBORLIST
                PMEDouble atomX             = cSim.pImageX[atom.x];
                PMEDouble atomY             = cSim.pImageY[atom.x];
                PMEDouble atomZ             = cSim.pImageZ[atom.x];  
#else
                PMEDouble atomX             = cSim.pAtomX[atom.x];
                PMEDouble atomY             = cSim.pAtomY[atom.x];
                PMEDouble atomZ             = cSim.pAtomZ[atom.x];  
#endif 
                PMEDouble ax                = atomX - constraint1.y;
                PMEDouble ay                = atomY - constraint2.x;
#else
                int2 iatomX                 = tex1Dfetch(texref, atom.x);          
                int2 iatomY                 = tex1Dfetch(texref, atom.x + cSim.stride);
                int2 iatomZ                 = tex1Dfetch(texref, atom.x + cSim.stride2);
                PMEDouble atomX             = __hiloint2double(iatomX.y, iatomX.x);
                PMEDouble ax                = atomX - constraint1.y;
                PMEDouble atomY             = __hiloint2double(iatomY.y, iatomY.x);
                PMEDouble ay                = atomY - constraint2.x;            
                PMEDouble atomZ             = __hiloint2double(iatomZ.y, iatomZ.x);
#endif            
                PMEDouble az                = atomZ - constraint2.y;
                PMEDouble wx                = constraint1.x * ax;
                PMEDouble wy                = constraint1.x * ay;
                PMEDouble wz                = constraint1.x * az;
#ifdef LOCAL_ENERGY            
                PMEDouble eadd              = wx * ax + wy * ay + wz * az;       
                sE[threadIdx.x].restraint  += eadd;
#endif
#ifdef LOCAL_NEIGHBORLIST
#ifdef LOCAL_VIRIAL
                PMEDouble fx                = (PMEDouble)-2.0 * wx;
                PMEDouble fy                = (PMEDouble)-2.0 * wy;
                PMEDouble fz                = (PMEDouble)-2.0 * wz;
                unsigned long long int ifx  = fabs(fx) * FORCESCALE + (PMEDouble)0.5;
                unsigned long long int ify  = fabs(fy) * FORCESCALE + (PMEDouble)0.5;
                unsigned long long int ifz  = fabs(fz) * FORCESCALE + (PMEDouble)0.5;
                if (fx < 0)
                    ifx                     = 0ull - ifx;
                if (fy < 0)
                    ify                     = 0ull - ify;                    
                if (fz < 0)
                    ifz                     = 0ull - ifz;                
                atomicAdd(&cSim.pUllForceX[atom.x], ifx); 
                atomicAdd(&cSim.pUllForceY[atom.x], ify);  
                atomicAdd(&cSim.pUllForceZ[atom.x], ifz);  
#else
                cSim.pForceXBuffer[atom.y]  = (PMEDouble)-2.0 * wx;
                cSim.pForceYBuffer[atom.y]  = (PMEDouble)-2.0 * wy;
                cSim.pForceZBuffer[atom.y]  = (PMEDouble)-2.0 * wz;
#endif
#else
                cSim.pForceXBuffer[atom.y] += (PMEDouble)-2.0 * wx;
                cSim.pForceYBuffer[atom.y] += (PMEDouble)-2.0 * wy;
                cSim.pForceZBuffer[atom.y] += (PMEDouble)-2.0 * wz;
#endif
#ifdef MPI
            }
#endif
        }
        pos                                += cSim.nb14Offset + blockDim.x * gridDim.x;
    }
    
#if defined(LOCAL_VIRIAL) || defined(LOCAL_ENERGY)    
    __syncthreads();
    unsigned int m                          = 1;
    while (m < blockDim.x)
    {
        int p                               = threadIdx.x + m;
#ifdef LOCAL_ENERGY        
        PMEDouble bond                      = ((p < blockDim.x) ? sE[p].bond : (PMEDouble)0.0);
        PMEDouble angle                     = ((p < blockDim.x) ? sE[p].angle : (PMEDouble)0.0);
        PMEDouble dihedral                  = ((p < blockDim.x) ? sE[p].dihedral : (PMEDouble)0.0);
        PMEDouble el14                      = ((p < blockDim.x) ? sE[p].el14 : (PMEDouble)0.0);
        PMEDouble nb14                      = ((p < blockDim.x) ? sE[p].nb14 : (PMEDouble)0.0);
        PMEDouble restraint                 = ((p < blockDim.x) ? sE[p].restraint : (PMEDouble)0.0); 
#endif
#ifdef LOCAL_VIRIAL
        PMEDouble vir_11                    = ((p < blockDim.x) ? sV[p].vir_11 : (PMEDouble)0.0);
        PMEDouble vir_22                    = ((p < blockDim.x) ? sV[p].vir_22 : (PMEDouble)0.0);
        PMEDouble vir_33                    = ((p < blockDim.x) ? sV[p].vir_33 : (PMEDouble)0.0);
#endif
        __syncthreads();
#ifdef LOCAL_ENERGY
        sE[threadIdx.x].bond               += bond;
        sE[threadIdx.x].angle              += angle;
        sE[threadIdx.x].dihedral           += dihedral;
        sE[threadIdx.x].el14               += el14;
        sE[threadIdx.x].nb14               += nb14;
        sE[threadIdx.x].restraint          += restraint;
#endif
#ifdef LOCAL_VIRIAL
        sV[threadIdx.x].vir_11             += vir_11;
        sV[threadIdx.x].vir_22             += vir_22;
        sV[threadIdx.x].vir_33             += vir_33;
#endif
        __syncthreads();
        m                                  *= 2;
    }
#ifdef LOCAL_ENERGY    
    unsigned long long int val1             = (unsigned long long int)(fabs(sE[threadIdx.x].bond)      * ENERGYSCALE + (PMEDouble)0.5);
    unsigned long long int val2             = (unsigned long long int)(fabs(sE[threadIdx.x].angle)     * ENERGYSCALE + (PMEDouble)0.5);
    unsigned long long int val3             = (unsigned long long int)(fabs(sE[threadIdx.x].dihedral)  * ENERGYSCALE + (PMEDouble)0.5);
    unsigned long long int val4             = (unsigned long long int)(fabs(sE[threadIdx.x].el14)      * ENERGYSCALE + (PMEDouble)0.5);
    unsigned long long int val5             = (unsigned long long int)(fabs(sE[threadIdx.x].nb14)      * ENERGYSCALE + (PMEDouble)0.5);
    unsigned long long int val6             = (unsigned long long int)(fabs(sE[threadIdx.x].restraint) * ENERGYSCALE + (PMEDouble)0.5);
    if (sE[threadIdx.x].bond < (PMEDouble)0.0)
        val1                                = 0ull - val1;
    if (sE[threadIdx.x].angle < (PMEDouble)0.0)
        val2                                = 0ull - val2;
    if (sE[threadIdx.x].dihedral < (PMEDouble)0.0)
        val3                                = 0ull - val3;    
    if (sE[threadIdx.x].el14 < (PMEDouble)0.0)
        val4                                = 0ull - val4;
    if (sE[threadIdx.x].nb14 < (PMEDouble)0.0)
        val5                                = 0ull - val5;
    if (sE[threadIdx.x].restraint < (PMEDouble)0.0)
        val6                                = 0ull - val6;    
#endif
#ifdef LOCAL_VIRIAL
    unsigned long long int val7             = (unsigned long long int)(fabs(sV[threadIdx.x].vir_11)    * ENERGYSCALE + (PMEDouble)0.5);
    unsigned long long int val8             = (unsigned long long int)(fabs(sV[threadIdx.x].vir_22)    * ENERGYSCALE + (PMEDouble)0.5);
    unsigned long long int val9             = (unsigned long long int)(fabs(sV[threadIdx.x].vir_33)    * ENERGYSCALE + (PMEDouble)0.5);
    if (sV[threadIdx.x].vir_11 < (PMEDouble)0.0)
        val7                                = 0ull - val7;
    if (sV[threadIdx.x].vir_22 < (PMEDouble)0.0)
        val8                                = 0ull - val8;
    if (sV[threadIdx.x].vir_33 < (PMEDouble)0.0)
        val9                                = 0ull - val9;        
#endif    
    if (threadIdx.x == 0)
    {
#ifdef LOCAL_ENERGY
        atomicAdd(cSim.pEBond, val1);
        atomicAdd(cSim.pEAngle, val2);
        atomicAdd(cSim.pEDihedral, val3);
        atomicAdd(cSim.pEEL14, val4);
        atomicAdd(cSim.pENB14, val5);
        atomicAdd(cSim.pERestraint, val6); 
#endif
#ifdef LOCAL_VIRIAL
        atomicAdd(cSim.pVirial_11, val7);
        atomicAdd(cSim.pVirial_22, val8);
        atomicAdd(cSim.pVirial_33, val9);
#endif
    }     
#endif     
}

