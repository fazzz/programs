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




#ifdef LOCAL_ENERGY
struct Energy {
    PMEDouble dihedral;
    PMEDouble el14;
    PMEDouble nb14;
#if (__CUDA_ARCH__ < 300) || !defined(LOCAL_VIRIAL)
    PMEDouble bond;
    PMEDouble angle;
    PMEDouble restraint;
#endif
};

#if (__CUDA_ARCH__ < 300) || !defined(LOCAL_VIRIAL)
#define EBOND      sE[threadIdx.x].bond
#define EANGLE     sE[threadIdx.x].angle
#define ERESTRAINT sE[threadIdx.x].restraint
#else
#define EBOND      ebond
#define EANGLE     eangle
#define ERESTRAINT erestraint
#endif

#if (__CUDA_ARCH__ >= 300)
    __shared__ Energy sE[SM_3X_LOCALFORCES_THREADS_PER_BLOCK];
#ifdef LOCAL_VIRIAL
    PMEDouble ebond;
    PMEDouble eangle;
    PMEDouble erestraint;
#endif
#elif (__CUDA_ARCH__ >= 200)
    __shared__ Energy sE[SM_2X_LOCALFORCES_THREADS_PER_BLOCK];
#else
    __shared__ Energy sE[SM_13_LOCALFORCES_THREADS_PER_BLOCK];
#endif
#endif

#ifdef LOCAL_VIRIAL
struct Virial {
    PMEDouble vir_11;
    PMEDouble vir_22;
    PMEDouble vir_33;
};
#if (__CUDA_ARCH__ >= 300)
    __shared__ Virial sV[SM_3X_LOCALFORCES_THREADS_PER_BLOCK];
#elif (__CUDA_ARCH__ >= 200)
    __shared__ Virial sV[SM_2X_LOCALFORCES_THREADS_PER_BLOCK];
#else
    __shared__ Virial sV[SM_13_LOCALFORCES_THREADS_PER_BLOCK];
#endif
#endif

    int pos                                 = blockIdx.x * blockDim.x + threadIdx.x;    
#ifdef LOCAL_ENERGY
    EBOND                                   = (PMEDouble)0.0;
    EANGLE                                  = (PMEDouble)0.0;
    sE[threadIdx.x].dihedral                = (PMEDouble)0.0;
    sE[threadIdx.x].el14                    = (PMEDouble)0.0;
    sE[threadIdx.x].nb14                    = (PMEDouble)0.0;
    ERESTRAINT                              = (PMEDouble)0.0;
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
        
#if defined(NODPTEXTURE)
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
                EBOND                      += df * da;
#endif
            
            // Calculate force
            PMEDouble fx                    = dfw * xij;
            PMEDouble fy                    = dfw * yij;
            PMEDouble fz                    = dfw * zij;    

#if defined(LOCAL_VIRIAL) || defined(use_SPFP)
            PMEAccumulator ifx              = lliroundd(fx * FORCESCALE);
            PMEAccumulator ify              = lliroundd(fy * FORCESCALE);
            PMEAccumulator ifz              = lliroundd(fz * FORCESCALE);
#endif

#ifdef MPI
            if ((atom.y >= cSim.minLocalAtom) && (atom.y < cSim.maxLocalAtom))
#endif
            {
#if defined(LOCAL_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom.y], llitoulli(ifx)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom.y], llitoulli(ify));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom.y], llitoulli(ifz));               
#elif defined (LOCAL_NEIGHBORLIST)                
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
#if defined(LOCAL_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom.x], llitoulli(-ifx)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom.x], llitoulli(-ify));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom.x], llitoulli(-ifz));   
#elif defined(LOCAL_NEIGHBORLIST)                   
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
#if !defined(LOCAL_VIRIAL) && !defined(use_SPFP)
            int2 atom2                      = cSim.pImageBondAngleID2[pos];               
#endif
#elif !defined(use_SPFP)
            int2 atom2                      = cSim.pBondAngleID2[pos];       
#endif        
      
#if defined(NODPTEXTURE)
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
                EANGLE                     += df * da;
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
#if defined(LOCAL_VIRIAL) || defined(use_SPFP)
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
#if defined(LOCAL_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.x], llitoulli(ifx1)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.x], llitoulli(ify1));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.x], llitoulli(ifz1));    
#elif defined(LOCAL_NEIGHBORLIST)            
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
#if defined(LOCAL_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.z], llitoulli(ifx2)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.z], llitoulli(ify2));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.z], llitoulli(ifz2));  
#elif defined (LOCAL_NEIGHBORLIST)                   
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
#if defined(LOCAL_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.y], llitoulli(-ifx1 - ifx2)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.y], llitoulli(-ify1 - ify2));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.y], llitoulli(-ifz1 - ifz2));                   
#elif defined(LOCAL_NEIGHBORLIST)              
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
#if !defined(LOCAL_VIRIAL) && !defined(use_SPFP)            
            int4 atom2                      = cSim.pImageDihedralID2[pos];
#endif
#else
            int4 atom1                      = cSim.pDihedralID1[pos];
#ifndef use_SPFP
            int4 atom2                      = cSim.pDihedralID2[pos];
#endif
#endif

#if defined(NODPTEXTURE)
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
            faster_sincos(ap, &sphi, &cphi);


            // Calculate the energy and the derivatives with respect to cosphi
            PMEDouble ct0                   = dihedral1.y * ap;
            PMEDouble sinnp, cosnp;
            faster_sincos(ct0, &sinnp, &cosnp);
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
	        if ((cSim.iamd == 2) || (cSim.iamd == 3))
            {
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

#if defined(LOCAL_VIRIAL) || defined(use_SPFP)
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
#if defined(LOCAL_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.x], llitoulli(-idr1)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.x], llitoulli(-idr2));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.x], llitoulli(-idr3));  
#elif defined(LOCAL_NEIGHBORLIST)            
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
#if defined(LOCAL_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.y], llitoulli(-idrx + idr1)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.y], llitoulli(-idry + idr2));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.y], llitoulli(-idrz + idr3));  
#elif defined(LOCAL_NEIGHBORLIST)   
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
#if defined(LOCAL_VIRIAL) || defined(use_SPFP)          
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.z], llitoulli(idrx + idr4)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.z], llitoulli(idry + idr5));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.z], llitoulli(idrz + idr6));  
#elif defined(LOCAL_NEIGHBORLIST)  
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
#if defined(LOCAL_VIRIAL) || defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom1.w], llitoulli(-idr4)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom1.w], llitoulli(-idr5));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom1.w], llitoulli(-idr6));  
#elif defined(LOCAL_NEIGHBORLIST)            
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
      
#if defined(NODPTEXTURE)
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

#if defined(use_SPFP)
            PMEAccumulator ifx              = lliroundd(fx * FORCESCALE);
            PMEAccumulator ify              = lliroundd(fy * FORCESCALE);
            PMEAccumulator ifz              = lliroundd(fz * FORCESCALE);
#endif
        
#ifdef MPI
            if ((atom.y >= cSim.minLocalAtom) && (atom.y < cSim.maxLocalAtom))
#endif
            {
#if defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[atom.y], llitoulli(ifx));
                atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[atom.y], llitoulli(ify));
                atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[atom.y], llitoulli(ifz));
#elif defined(LOCAL_NEIGHBORLIST)
                cSim.pForceXBuffer[atom.w]  =  fx;
                cSim.pForceYBuffer[atom.w]  =  fy;
                cSim.pForceZBuffer[atom.w]  =  fz;
#else
                cSim.pForceXBuffer[atom.w] +=  fx;
                cSim.pForceYBuffer[atom.w] +=  fy;
                cSim.pForceZBuffer[atom.w] +=  fz;
#endif
            }
#ifdef MPI
            if ((atom.x >= cSim.minLocalAtom) && (atom.x < cSim.maxLocalAtom))
#endif
            {       
#if defined(use_SPFP)
                atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[atom.x], llitoulli(-ifx));
                atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[atom.x], llitoulli(-ify));
                atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[atom.x], llitoulli(-ifz));
#elif defined(LOCAL_NEIGHBORLIST)    
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
#endif 
            {

                PMEDouble2 constraint1      = cSim.pConstraint1[pos];
                PMEDouble2 constraint2      = cSim.pConstraint2[pos];
           
#if defined(NODPTEXTURE)
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
                ERESTRAINT                 += eadd;
#endif

#if defined(LOCAL_VIRIAL) || defined(use_SPFP)
                PMEAccumulator ifx          = lliroundd(wx * (PMEDouble)-2.0 * FORCESCALE);
                PMEAccumulator ify          = lliroundd(wy * (PMEDouble)-2.0 * FORCESCALE);
                PMEAccumulator ifz          = lliroundd(wz * (PMEDouble)-2.0 * FORCESCALE);
                atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atom.x], llitoulli(ifx)); 
                atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atom.x], llitoulli(ify));  
                atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atom.x], llitoulli(ifz));  
#elif defined(LOCAL_NEIGHBORLIST)
                cSim.pForceXBuffer[atom.y]  = (PMEDouble)-2.0 * wx;
                cSim.pForceYBuffer[atom.y]  = (PMEDouble)-2.0 * wy;
                cSim.pForceZBuffer[atom.y]  = (PMEDouble)-2.0 * wz;
#else
                cSim.pForceXBuffer[atom.y] += (PMEDouble)-2.0 * wx;
                cSim.pForceYBuffer[atom.y] += (PMEDouble)-2.0 * wy;
                cSim.pForceZBuffer[atom.y] += (PMEDouble)-2.0 * wz;
#endif
            }
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
        PMEDouble el14                      = ((p < blockDim.x) ? sE[p].el14 : (PMEDouble)0.0);
        PMEDouble nb14                      = ((p < blockDim.x) ? sE[p].nb14 : (PMEDouble)0.0);
        PMEDouble dihedral                  = ((p < blockDim.x) ? sE[p].dihedral : (PMEDouble)0.0);
#if (__CUDA_ARCH__ < 300) || !defined(LOCAL_VIRIAL)
        PMEDouble bond                      = ((p < blockDim.x) ? sE[p].bond : (PMEDouble)0.0);
        PMEDouble angle                     = ((p < blockDim.x) ? sE[p].angle : (PMEDouble)0.0);
        PMEDouble restraint                 = ((p < blockDim.x) ? sE[p].restraint : (PMEDouble)0.0); 
#endif
#endif
#ifdef LOCAL_VIRIAL
        PMEDouble vir_11                    = ((p < blockDim.x) ? sV[p].vir_11 : (PMEDouble)0.0);
        PMEDouble vir_22                    = ((p < blockDim.x) ? sV[p].vir_22 : (PMEDouble)0.0);
        PMEDouble vir_33                    = ((p < blockDim.x) ? sV[p].vir_33 : (PMEDouble)0.0);
#endif
        __syncthreads();
#ifdef LOCAL_ENERGY
        sE[threadIdx.x].dihedral           += dihedral;
        sE[threadIdx.x].el14               += el14;
        sE[threadIdx.x].nb14               += nb14;
#if (__CUDA_ARCH__ < 300) || !defined(LOCAL_VIRIAL)
        sE[threadIdx.x].bond               += bond;
        sE[threadIdx.x].angle              += angle;
        sE[threadIdx.x].restraint          += restraint;
#endif
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
    unsigned long long int val3             = llitoulli(lliroundd(sE[threadIdx.x].dihedral  * ENERGYSCALE));
    unsigned long long int val4             = llitoulli(lliroundd(sE[threadIdx.x].el14      * ENERGYSCALE));
    unsigned long long int val5             = llitoulli(lliroundd(sE[threadIdx.x].nb14      * ENERGYSCALE));
#if (__CUDA_ARCH__ < 300) || !defined(LOCAL_VIRIAL)
    unsigned long long int val1             = llitoulli(lliroundd(sE[threadIdx.x].bond      * ENERGYSCALE));
    unsigned long long int val2             = llitoulli(lliroundd(sE[threadIdx.x].angle     * ENERGYSCALE));
    unsigned long long int val6             = llitoulli(lliroundd(sE[threadIdx.x].restraint * ENERGYSCALE)); 
#else
    m                                       = 1;
    sE[threadIdx.x].dihedral                = EBOND;
    sE[threadIdx.x].el14                    = EANGLE;
    sE[threadIdx.x].nb14                    = ERESTRAINT;
    __syncthreads();
    while (m < blockDim.x)
    {
        int p                               = threadIdx.x + m;
        PMEDouble bond                      = ((p < blockDim.x) ? sE[p].dihedral : (PMEDouble)0.0);
        PMEDouble angle                     = ((p < blockDim.x) ? sE[p].el14     : (PMEDouble)0.0);
        PMEDouble restraint                 = ((p < blockDim.x) ? sE[p].nb14     : (PMEDouble)0.0);
        __syncthreads();
        sE[threadIdx.x].dihedral           += bond;
        sE[threadIdx.x].el14               += angle;
        sE[threadIdx.x].nb14               += restraint;
        __syncthreads();
        m                                  *= 2;
    }
    unsigned long long int val1             = llitoulli(lliroundd(sE[threadIdx.x].dihedral  * ENERGYSCALE));
    unsigned long long int val2             = llitoulli(lliroundd(sE[threadIdx.x].el14      * ENERGYSCALE));
    unsigned long long int val6             = llitoulli(lliroundd(sE[threadIdx.x].nb14      * ENERGYSCALE)); 
#endif
#endif
#ifdef LOCAL_VIRIAL
    unsigned long long int val7             = llitoulli(lliroundd(sV[threadIdx.x].vir_11    * FORCESCALE));
    unsigned long long int val8             = llitoulli(lliroundd(sV[threadIdx.x].vir_22    * FORCESCALE));
    unsigned long long int val9             = llitoulli(lliroundd(sV[threadIdx.x].vir_33    * FORCESCALE));
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
#undef EBOND
#undef EANGLE
#undef ERESTRAINT
#endif
#ifdef LOCAL_VIRIAL
        atomicAdd(cSim.pVirial_11, val7);
        atomicAdd(cSim.pVirial_22, val8);
        atomicAdd(cSim.pVirial_33, val9);
#endif
    }     
#endif     
}

