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
struct CHARMMEnergy {
    PMEDouble angle_ub;
    PMEDouble imp;
    PMEDouble cmap;
};
#ifdef CHARMM_ENERGY
#if (__CUDA_ARCH__ >= 200)
    __shared__ CHARMMEnergy sE[SM_2X_CHARMMFORCES_THREADS_PER_BLOCK];
#else
    __shared__ CHARMMEnergy sE[SM_13_CHARMMFORCES_THREADS_PER_BLOCK];
#endif
    sE[threadIdx.x].angle_ub                = 0.0;
    sE[threadIdx.x].imp                     = 0.0;
    sE[threadIdx.x].cmap                    = 0.0;
#endif

    unsigned int pos                        = blockIdx.x * blockDim.x + threadIdx.x;    
    // Calculate Urey Bradly angle forces
    while (pos < cSim.UBAngleOffset)
    {
        if (pos < cSim.UBAngles)
        {
#ifdef CHARMM_NEIGHBORLIST           
            int4 atom                       = cSim.pImageUBAngleID[pos];
#else
            int4 atom                       = cSim.pUBAngleID[pos];
#endif
            PMEDouble2 UBAngle              = cSim.pUBAngle[pos];
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
#ifdef CHARMM_NEIGHBORLIST
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
            double da                       = rij - UBAngle.y;
            double df                       = UBAngle.x * da;
#ifdef CHARMM_ENERGY            
#ifdef MPI
            if ((atom.x >= cSim.minLocalAtom) && (atom.x < cSim.maxLocalAtom))  
#endif
                sE[threadIdx.x].angle_ub    += df * da;
#endif
            double dfw                      = (df + df) / rij;
            
            // Calculate force
            PMEDouble fx                    = dfw * xij;
            PMEDouble fy                    = dfw * yij;
            PMEDouble fz                    = dfw * zij;    
#ifdef CHARMM_NEIGHBORLIST
#ifdef CHARMM_VIRIAL
            unsigned long long int ifx      = fabs(fx) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int ify      = fabs(fy) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int ifz      = fabs(fz) * FORCESCALE + (PMEDouble)0.5;
            if (fx < (PMEDouble)0.0)
                ifx                         = 0ull - ifx;
            if (fy < (PMEDouble)0.0)
                ify                         = 0ull - ify;                    
            if (fz < (PMEDouble)0.0)
                ifz                         = 0ull - ifz;
#endif

#ifdef MPI
            if ((atom.y >= cSim.minLocalAtom) && (atom.y < cSim.maxLocalAtom))  
#endif
            {
#ifdef CHARMM_VIRIAL
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
#ifdef CHARMM_VIRIAL
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

    // Calculate improper dihedral forces
    while (pos < cSim.impDihedralOffset)
    {
        pos                                -= cSim.UBAngleOffset;
        if (pos < cSim.impDihedrals)
        {
#ifdef CHARMM_NEIGHBORLIST
            int4 atom1                      = cSim.pImageImpDihedralID1[pos];
#ifndef CHARMM_VIRIAL
            int4 atom2                      = cSim.pImageImpDihedralID2[pos];
#endif
#else
            int4 atom1                      = cSim.pImpDihedralID1[pos];
            int4 atom2                      = cSim.pImpDihedralID2[pos];
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
            PMEDouble2 impDihedral          = cSim.pImpDihedral[pos];
            PMEDouble yij                   = atomIY - atomJY;      
            PMEDouble zkj                   = atomKZ - atomJZ;  
            PMEDouble zij                   = atomIZ - atomJZ;
            PMEDouble xkj                   = atomKX - atomJX;  
            PMEDouble xij                   = atomIX - atomJX; 
            PMEDouble ykj                   = atomKY - atomJY;         
            PMEDouble ylk                   = atomLY - atomKY; 
#elif defined(NODPTEXTURE)
#ifdef CHARMM_NEIGHBORLIST
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
            PMEDouble2 impDihedral          = cSim.pImpDihedral[pos];
            PMEDouble xij                   = atomIX - atomJX;
            PMEDouble yij                   = atomIY - atomJY;
            PMEDouble zij                   = atomIZ - atomJZ;
            PMEDouble xkj                   = atomKX - atomJX;
            PMEDouble ykj                   = atomKY - atomJY;
            PMEDouble zkj                   = atomKZ - atomJZ;
            PMEDouble ylk                   = atomLY - atomKY;
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
            PMEDouble2 impDihedral          = cSim.pImpDihedral[pos];
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
            PMEDouble ylk                   = atomLY - atomKY; 
            PMEDouble atomLX                = __hiloint2double(iatomLX.y, iatomLX.x);                     
            PMEDouble atomLZ                = __hiloint2double(iatomLZ.y, iatomLZ.x); 
#endif            
            PMEDouble zlk                   = atomLZ - atomKZ;                                           
            PMEDouble xlk                   = atomLX - atomKX;
            
            // Calculate phi and quantities required for gradient
            PMEDouble oneOverRKJ            = rsqrt(xkj * xkj + ykj * ykj + zkj * zkj);             
            PMEDouble uxkj                  = xkj * oneOverRKJ;
            PMEDouble uykj                  = ykj * oneOverRKJ;
            PMEDouble uzkj                  = zkj * oneOverRKJ;
            PMEDouble dotIJKJ               = xij * uxkj + yij * uykj + zij * uzkj;  
            PMEDouble upxij                 = xij - dotIJKJ * uxkj;
            PMEDouble upyij                 = yij - dotIJKJ * uykj;
            PMEDouble upzij                 = zij - dotIJKJ * uzkj;
            dotIJKJ                        *= oneOverRKJ;
            PMEDouble oneOverRUIJ           = rsqrt(upxij * upxij + upyij * upyij + upzij * upzij);
            upxij                          *= oneOverRUIJ;
            upyij                          *= oneOverRUIJ;
            upzij                          *= oneOverRUIJ;
            PMEDouble dotLKKJ               = xlk * uxkj + ylk * uykj + zlk * uzkj;
            PMEDouble upxlk                 = xlk - dotLKKJ * uxkj;
            PMEDouble upylk                 = ylk - dotLKKJ * uykj;
            PMEDouble upzlk                 = zlk - dotLKKJ * uzkj;
            dotLKKJ                        *= oneOverRKJ;
            PMEDouble oneOverRULK           = rsqrt(upxlk * upxlk + upylk * upylk + upzlk * upzlk);
            upxlk                          *= oneOverRULK;    
            upylk                          *= oneOverRULK; 
            upzlk                          *= oneOverRULK;    
            PMEDouble dot                   = upxij * upxlk + upyij * upylk + upzij * upzlk;
            PMEDouble cosphi                = min(max(dot, (PMEDouble)-1.0), (PMEDouble)1.0);
            PMEDouble cx                    = upyij * upzlk - upzij * upylk;
            PMEDouble cy                    = upzij * upxlk - upxij * upzlk;
            PMEDouble cz                    = upxij * upylk - upyij * upxlk;
            dot                             = cx * uxkj + cy * uykj + cz * uzkj;
            PMEDouble sinphi                = min(max(dot, (PMEDouble)-1.0), (PMEDouble)1.0);        
            PMEDouble phi                   = acos(cosphi) * (sinphi >= (PMEDouble)0.0 ? (PMEDouble)1.0 : (PMEDouble)-1.0);  
            PMEDouble df                    = (PMEDouble)-2.0 * impDihedral.x * (phi - impDihedral.y);   

#ifdef CHARMM_ENERGY
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom))  
#endif
                sE[threadIdx.x].imp        += impDihedral.x * (phi - impDihedral.y) * (phi - impDihedral.y);
#endif                   
        
            // Calculate gradient
            PMEDouble upxijk                = df * (upyij * uzkj  - upzij * uykj)  * oneOverRUIJ;
            PMEDouble upyijk                = df * (upzij * uxkj  - upxij * uzkj)  * oneOverRUIJ;
            PMEDouble upzijk                = df * (upxij * uykj  - upyij * uxkj)  * oneOverRUIJ;
            PMEDouble upxjkl                = df * (uykj  * upzlk - uzkj  * upylk) * oneOverRULK;
            PMEDouble upyjkl                = df * (uzkj  * upxlk - uxkj  * upzlk) * oneOverRULK;
            PMEDouble upzjkl                = df * (uxkj  * upylk - uykj  * upxlk) * oneOverRULK;
            PMEDouble vx                    = dotIJKJ * upxijk + dotLKKJ * upxjkl;
            PMEDouble vy                    = dotIJKJ * upyijk + dotLKKJ * upyjkl;
            PMEDouble vz                    = dotIJKJ * upzijk + dotLKKJ * upzjkl;     
#ifdef CHARMM_NEIGHBORLIST
#ifdef CHARMM_VIRIAL
            unsigned long long int iupxijk  = fabs(upxijk) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iupyijk  = fabs(upyijk) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iupzijk  = fabs(upzijk) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iupxjkl  = fabs(upxjkl) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iupyjkl  = fabs(upyjkl) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iupzjkl  = fabs(upzjkl) * FORCESCALE + (PMEDouble)0.5;   
            unsigned long long int ivx      = fabs(vx) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int ivy      = fabs(vy) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int ivz      = fabs(vz) * FORCESCALE + (PMEDouble)0.5;  
            if (upxijk < (PMEDouble)0.0)
                iupxijk                     = 0ull - iupxijk;
            if (upyijk < (PMEDouble)0.0)
                iupyijk                     = 0ull - iupyijk;                    
            if (upzijk < (PMEDouble)0.0)
                iupzijk                     = 0ull - iupzijk;   
            if (upxjkl < (PMEDouble)0.0)
                iupxjkl                     = 0ull - iupxjkl;
            if (upyjkl < (PMEDouble)0.0)
                iupyjkl                     = 0ull - iupyjkl;                    
            if (upzjkl < (PMEDouble)0.0)
                iupzjkl                     = 0ull - iupzjkl;   
            if (vx < (PMEDouble)0.0)
                ivx                         = 0ull - ivx;
            if (vy < (PMEDouble)0.0)
                ivy                         = 0ull - ivy;                    
            if (vz < (PMEDouble)0.0)
                ivz                         = 0ull - ivz;                                                                   
#endif
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom))  
#endif
            {
#ifdef CHARMM_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom1.x], iupxijk); 
                atomicAdd(&cSim.pUllForceY[atom1.x], iupyijk);  
                atomicAdd(&cSim.pUllForceZ[atom1.x], iupzijk);       
#else            
                cSim.pForceXBuffer[atom2.x] =  upxijk;
                cSim.pForceYBuffer[atom2.x] =  upyijk;
                cSim.pForceZBuffer[atom2.x] =  upzijk;
#endif
            }
#ifdef MPI
            if ((atom1.y >= cSim.minLocalAtom) && (atom1.y < cSim.maxLocalAtom))  
#endif
            {
#ifdef CHARMM_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom1.y], ivx - iupxijk);
                atomicAdd(&cSim.pUllForceY[atom1.y], ivy - iupyijk);  
                atomicAdd(&cSim.pUllForceZ[atom1.y], ivz - iupzijk);       
#else            
                cSim.pForceXBuffer[atom2.y] =  vx - upxijk;
                cSim.pForceYBuffer[atom2.y] =  vy - upyijk;
                cSim.pForceZBuffer[atom2.y] =  vz - upzijk;
#endif                
            }
#ifdef MPI
            if ((atom1.z >= cSim.minLocalAtom) && (atom1.z < cSim.maxLocalAtom))  
#endif
            {
#ifdef CHARMM_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom1.z], 0ull - ivx - iupxjkl); 
                atomicAdd(&cSim.pUllForceY[atom1.z], 0ull - ivy - iupyjkl);  
                atomicAdd(&cSim.pUllForceZ[atom1.z], 0ull - ivz - iupzjkl);       
#else
                cSim.pForceXBuffer[atom2.z] = -vx - upxjkl;
                cSim.pForceYBuffer[atom2.z] = -vy - upyjkl;
                cSim.pForceZBuffer[atom2.z] = -vz - upzjkl;
#endif
            }
#ifdef MPI
            if ((atom1.w >= cSim.minLocalAtom) && (atom1.w < cSim.maxLocalAtom))  
#endif
            {
#ifdef CHARMM_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom1.w], iupxjkl); 
                atomicAdd(&cSim.pUllForceY[atom1.w], iupyjkl);  
                atomicAdd(&cSim.pUllForceZ[atom1.w], iupzjkl);       
#else
                cSim.pForceXBuffer[atom2.w] =  upxjkl;
                cSim.pForceYBuffer[atom2.w] =  upyjkl;
                cSim.pForceZBuffer[atom2.w] =  upzjkl;
#endif
            }
#else
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom))  
#endif
            {
                cSim.pForceXBuffer[atom2.x]    +=  upxijk;
                cSim.pForceYBuffer[atom2.x]    +=  upyijk;
                cSim.pForceZBuffer[atom2.x]    +=  upzijk;
            }
#ifdef MPI
            if ((atom1.y >= cSim.minLocalAtom) && (atom1.y < cSim.maxLocalAtom))  
#endif            
            {
                cSim.pForceXBuffer[atom2.y]    +=  vx - upxijk;
                cSim.pForceYBuffer[atom2.y]    +=  vy - upyijk;
                cSim.pForceZBuffer[atom2.y]    +=  vz - upzijk;
            }
#ifdef MPI
            if ((atom1.z >= cSim.minLocalAtom) && (atom1.z < cSim.maxLocalAtom))  
#endif            
            {
                cSim.pForceXBuffer[atom2.z]    += -vx - upxjkl;
                cSim.pForceYBuffer[atom2.z]    += -vy - upyjkl;
                cSim.pForceZBuffer[atom2.z]    += -vz - upzjkl;
            }
#ifdef MPI
            if ((atom1.w >= cSim.minLocalAtom) && (atom1.w < cSim.maxLocalAtom))  
#endif            
            {
                cSim.pForceXBuffer[atom2.w]    +=  upxjkl;
                cSim.pForceYBuffer[atom2.w]    +=  upyjkl;
                cSim.pForceZBuffer[atom2.w]    +=  upzjkl;
            }
#endif                      
        }
        pos                                += cSim.UBAngleOffset + blockDim.x * gridDim.x;
    }

    // Calculate cmap forces
    while (pos < cSim.cmapOffset)
    {
        pos                                -= cSim.impDihedralOffset;
        if (pos < cSim.cmaps)
        {
#ifdef CHARMM_NEIGHBORLIST
            int4 atom1                      = cSim.pImageCmapID1[pos];
#else
            int4 atom1                      = cSim.pCmapID1[pos];
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
            PMEDouble yij                   = atomIY - atomJY;      
            PMEDouble zkj                   = atomKZ - atomJZ;  
            PMEDouble zij                   = atomIZ - atomJZ;
            PMEDouble xkj                   = atomKX - atomJX;  
            PMEDouble xij                   = atomIX - atomJX; 
            PMEDouble ykj                   = atomKY - atomJY;         
            PMEDouble xlk                   = atomLX - atomKX;
            PMEDouble ylk                   = atomLY - atomKY;
            PMEDouble zlk                   = atomLZ - atomKZ;
#elif defined(NODPTEXTURE)
#ifdef CHARMM_NEIGHBORLIST
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
            PMEDouble xlk                   = atomLX - atomKX;
            PMEDouble ylk                   = atomLY - atomKY;
            PMEDouble zlk                   = atomLZ - atomKZ;
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
            PMEDouble ylk                   = atomLY - atomKY; 
            PMEDouble atomLX                = __hiloint2double(iatomLX.y, iatomLX.x);                     
            PMEDouble atomLZ                = __hiloint2double(iatomLZ.y, iatomLZ.x); 
            PMEDouble xlk                   = atomLX - atomKX; 
            PMEDouble zlk                   = atomLZ - atomKZ; 
#endif    
            // Calculate phi and quantities required for gradient
            PMEDouble oneOverRKJ            = rsqrt(xkj * xkj + ykj * ykj + zkj * zkj);             
            PMEDouble uxkj                  = xkj * oneOverRKJ;
            PMEDouble uykj                  = ykj * oneOverRKJ;
            PMEDouble uzkj                  = zkj * oneOverRKJ;
            PMEDouble dotIJKJ               = xij * uxkj + yij * uykj + zij * uzkj;  
            PMEDouble upxij                 = xij - dotIJKJ * uxkj;
            PMEDouble upyij                 = yij - dotIJKJ * uykj;
            PMEDouble upzij                 = zij - dotIJKJ * uzkj;
            dotIJKJ                        *= oneOverRKJ;
            PMEDouble oneOverRUIJ           = rsqrt(upxij * upxij + upyij * upyij + upzij * upzij);
            upxij                          *= oneOverRUIJ;
            upyij                          *= oneOverRUIJ;
            upzij                          *= oneOverRUIJ;
            PMEDouble dotLKKJ               = xlk * uxkj + ylk * uykj + zlk * uzkj;
            PMEDouble upxlk                 = xlk - dotLKKJ * uxkj;
            PMEDouble upylk                 = ylk - dotLKKJ * uykj;
            PMEDouble upzlk                 = zlk - dotLKKJ * uzkj;
            dotLKKJ                        *= oneOverRKJ;
            PMEDouble oneOverRULK           = rsqrt(upxlk * upxlk + upylk * upylk + upzlk * upzlk);
            upxlk                          *= oneOverRULK;    
            upylk                          *= oneOverRULK; 
            upzlk                          *= oneOverRULK;    
            PMEDouble dot                   = upxij * upxlk + upyij * upylk + upzij * upzlk;
            PMEDouble cosphi                = min(max(dot, (PMEDouble)-1.0), (PMEDouble)1.0);
            PMEDouble cx                    = upyij * upzlk - upzij * upylk;
            PMEDouble cy                    = upzij * upxlk - upxij * upzlk;
            PMEDouble cz                    = upxij * upylk - upyij * upxlk;
            dot                             = cx * uxkj + cy * uykj + cz * uzkj;
            PMEDouble sinphi                = min(max(dot, (PMEDouble)-1.0), (PMEDouble)1.0);        
            PMEDouble phi                   = acos(cosphi) * (sinphi >= 0.0 ? 1.0 : -1.0);  
            PMEDouble upxijk                = (upyij * uzkj  - upzij * uykj)  * oneOverRUIJ;
            PMEDouble upyijk                = (upzij * uxkj  - upxij * uzkj)  * oneOverRUIJ;
            PMEDouble upzijk                = (upxij * uykj  - upyij * uxkj)  * oneOverRUIJ;
            PMEDouble upxjkl                = (uykj  * upzlk - uzkj  * upylk) * oneOverRULK;
            PMEDouble upyjkl                = (uzkj  * upxlk - uxkj  * upzlk) * oneOverRULK;
            PMEDouble upzjkl                = (uxkj  * upylk - uykj  * upxlk) * oneOverRULK;
        
            // Calculate Psi       
#ifdef CHARMM_NEIGHBORLIST
            int4 atom2                      = cSim.pImageCmapID2[pos];
#else
            int4 atom2                      = cSim.pCmapID2[pos];
#endif 
        
#ifdef use_SPSP    
            PMEDouble atomMX                = tex1Dfetch(texref, atom2.x);
            PMEDouble atomMY                = tex1Dfetch(texref, atom2.x  + cSim.stride);
            PMEDouble atomMZ                = tex1Dfetch(texref, atom2.x  + cSim.stride2);
#elif defined(NODPTEXTURE)
#ifdef CHARMM_NEIGHBORLIST
            PMEDouble atomMX                = cSim.pImageX[atom2.x];
            PMEDouble atomMY                = cSim.pImageY[atom2.x];    
            PMEDouble atomMZ                = cSim.pImageZ[atom2.x];
#else
            PMEDouble atomMX                = cSim.pAtomX[atom2.x];
            PMEDouble atomMY                = cSim.pAtomY[atom2.x];    
            PMEDouble atomMZ                = cSim.pAtomZ[atom2.x];
#endif
#else        
            int2 iatomMX                    = tex1Dfetch(texref, atom2.x);
            int2 iatomMY                    = tex1Dfetch(texref, atom2.x + cSim.stride);
            int2 iatomMZ                    = tex1Dfetch(texref, atom2.x + cSim.stride2);
            PMEDouble atomMX                = __hiloint2double(iatomMX.y, iatomMX.x);
            PMEDouble atomMY                = __hiloint2double(iatomMY.y, iatomMY.x);
            PMEDouble atomMZ                = __hiloint2double(iatomMZ.y, iatomMZ.x); 
#endif        
            PMEDouble oneOverRLK            = rsqrt(xlk * xlk + ylk * ylk + zlk * zlk);  
            PMEDouble uxlk                  = xlk * oneOverRLK;
            PMEDouble uylk                  = ylk * oneOverRLK;
            PMEDouble uzlk                  = zlk * oneOverRLK;
            PMEDouble dotJKLK               = xkj * uxlk + ykj * uylk + zkj * uzlk;  
            PMEDouble upxjk                 = -xkj + dotJKLK * uxlk;
            PMEDouble upyjk                 = -ykj + dotJKLK * uylk;
            PMEDouble upzjk                 = -zkj + dotJKLK * uzlk;
            dotJKLK                        *= oneOverRLK;
            PMEDouble oneOverRUJK           = rsqrt(upxjk * upxjk + upyjk * upyjk + upzjk * upzjk);
            upxjk                          *= oneOverRUJK;
            upyjk                          *= oneOverRUJK;
            upzjk                          *= oneOverRUJK;
            PMEDouble yml                   = atomMY - atomLY;
            PMEDouble zml                   = atomMZ - atomLZ;                                           
            PMEDouble xml                   = atomMX - atomLX;
            PMEDouble dotMLLK               = xml * uxlk + yml * uylk + zml * uzlk;
            PMEDouble upxml                 = xml - dotMLLK * uxlk;
            PMEDouble upyml                 = yml - dotMLLK * uylk;
            PMEDouble upzml                 = zml - dotMLLK * uzlk;
            dotMLLK                        *= oneOverRLK;
            PMEDouble oneOverRUML           = rsqrt(upxml * upxml + upyml * upyml + upzml * upzml);
            upxml                          *= oneOverRUML;    
            upyml                          *= oneOverRUML; 
            upzml                          *= oneOverRUML;    
            dot                             = upxjk * upxml + upyjk * upyml + upzjk * upzml;
            PMEDouble cospsi                = min(max(dot, (PMEDouble)-1.0), (PMEDouble)1.0);
            cx                              = upyjk * upzml - upzjk * upyml;
            cy                              = upzjk * upxml - upxjk * upzml;
            cz                              = upxjk * upyml - upyjk * upxml;
            dot                             = cx * uxlk + cy * uylk + cz * uzlk;
            PMEDouble sinpsi                = min(max(dot, (PMEDouble)-1.0), (PMEDouble)1.0);     
            int cmapType                    = cSim.pCmapType[pos];   
            PMEDouble psi                   = acos(cospsi) * (sinpsi >= 0.0 ? 1.0 : -1.0);  
            PMEDouble upxjkl1               = (upyjk * uzlk  - upzjk * uylk)  * oneOverRUJK;
            PMEDouble upyjkl1               = (upzjk * uxlk  - upxjk * uzlk)  * oneOverRUJK;
            PMEDouble upzjkl1               = (upxjk * uylk  - upyjk * uxlk)  * oneOverRUJK;
            PMEDouble upxklm                = (uylk  * upzml - uzlk  * upyml) * oneOverRUML;
            PMEDouble upyklm                = (uzlk  * upxml - uxlk  * upzml) * oneOverRUML;
            PMEDouble upzklm                = (uxlk  * upyml - uylk  * upxml) * oneOverRUML;
            phi                            += PI;
            psi                            += PI;
            int x                           = phi * (CMAPRESOLUTION / (2.0 * PI));
            int y                           = psi * (CMAPRESOLUTION / (2.0 * PI));
            PMEFloat phifrac               = (phi - x * (2.0 * PI / CMAPRESOLUTION)) * CMAPRESOLUTION / (2.0 * PI);
            PMEFloat psifrac               = (psi - y * (2.0 * PI / CMAPRESOLUTION)) * CMAPRESOLUTION / (2.0 * PI);
            PMEFloat4* pSrc                 = &(cSim.pCmapEnergy[cmapType + y * cSim.cmapRowStride + x * cSim.cmapTermStride]);
            PMEFloat4 E00                   = *pSrc;
            PMEFloat4 E01                   = pSrc[cSim.cmapRowStride];
            PMEFloat4 E10                   = pSrc[cSim.cmapTermStride];
            PMEFloat4 E11                   = pSrc[cSim.cmapRowStride + cSim.cmapTermStride];
#ifdef CHARMM_ENERGY            
            PMEFloat a00                    =                   E00.x;
#endif            
            PMEFloat a10                    =                   E00.y;
            PMEFloat a20                    = -(PMEFloat)3.0 *  E00.x + (PMEFloat)3.0 * E10.x - (PMEFloat)2.0 * E00.y -                 E10.y;
            PMEFloat a30                    =  (PMEFloat)2.0 *  E00.x - (PMEFloat)2.0 * E10.x +                 E00.y +                 E10.y;
            PMEFloat a01                    =                   E00.z;
            PMEFloat a11                    =                   E00.w;
            PMEFloat a21                    = -(PMEFloat)3.0 *  E00.z + (PMEFloat)3.0 * E10.z - (PMEFloat)2.0 * E00.w -                 E10.w;
            PMEFloat a31                    =  (PMEFloat)2.0 *  E00.z - (PMEFloat)2.0 * E10.z +                 E00.w +                 E10.w;
            PMEFloat a02                    = -(PMEFloat)3.0 *  E00.x + (PMEFloat)3.0 * E01.x - (PMEFloat)2.0 * E00.z -                 E01.z;
            PMEFloat a12                    = -(PMEFloat)3.0 *  E00.y + (PMEFloat)3.0 * E01.y - (PMEFloat)2.0 * E00.w -                 E01.w;
            PMEFloat a22                    =  (PMEFloat)9.0 * (E00.x -                 E10.x -                 E01.x +                 E11.x) +
                                               (PMEFloat)6.0 *  E00.y + (PMEFloat)3.0 * E10.y - (PMEFloat)6.0 * E01.y - (PMEFloat)3.0 * E11.y  +
                                               (PMEFloat)6.0 *  E00.z - (PMEFloat)6.0 * E10.z + (PMEFloat)3.0 * E01.z - (PMEFloat)3.0 * E11.z  +
                                               (PMEFloat)4.0 *  E00.w + (PMEFloat)2.0 * E10.w + (PMEFloat)2.0 * E01.w +                 E11.w;            
            PMEFloat a32                    = -(PMEFloat)6.0 * (E00.x -                 E10.x -                 E01.x +                 E11.x) +
                                              -(PMEFloat)3.0 * (E00.y +                 E10.y -                 E01.y -                 E11.y) +                                           
                                              -(PMEFloat)4.0 *  E00.z + (PMEFloat)4.0 * E10.z -                 E01.w -                 E11.w  +
                                              -(PMEFloat)2.0 * (E00.w +                 E10.w +                 E01.z -                 E11.z);                                                      
            PMEFloat a03                    =  (PMEFloat)2.0 *  E00.x - (PMEFloat)2.0 * E01.x +                 E00.z +                 E01.z;
            PMEFloat a13                    =  (PMEFloat)2.0 *  E00.y - (PMEFloat)2.0 * E01.y +                 E00.w +                 E01.w;
            PMEFloat a23                    = -(PMEFloat)6.0 * (E00.x -                 E10.x -                 E01.x +                 E11.x) +
                                              -(PMEFloat)2.0 * (E00.w +                 E10.y +                 E01.w -                 E11.y) +
                                              -(PMEFloat)3.0 * (E00.z -                 E10.z +                 E01.z -                 E11.z) +
                                              -(PMEFloat)4.0 *  E00.y -                 E10.w + (PMEFloat)4.0 * E01.y -                 E11.w; 
            PMEFloat a33                    =  (PMEFloat)4.0 * (E00.x -                 E10.x -                 E01.x +                 E11.x) +
                                               (PMEFloat)2.0 * (E00.y +                 E10.y -                 E01.y -                 E11.y) +
                                               (PMEFloat)2.0 * (E00.z -                 E10.z +                 E01.z -                 E11.z) +
                                                                E00.w +                 E10.w +                 E01.w +                 E11.w;  

#ifdef CHARMM_NEIGHBORLIST
#ifndef CHARMM_VIRIAL
            int2 atom3                      = cSim.pImageCmapID3[pos];
#endif
#else
            int2 atom3                      = cSim.pCmapID3[pos];
#endif
#ifdef CHARMM_ENERGY
            PMEFloat E                      = ((a33 * psifrac + a32) * psifrac + a31) * psifrac + a30;
            E                               = ((a23 * psifrac + a22) * psifrac + a21) * psifrac + a20 + phifrac * E;
            E                               = ((a13 * psifrac + a12) * psifrac + a11) * psifrac + a10 + phifrac * E;
            E                               = ((a03 * psifrac + a02) * psifrac + a01) * psifrac + a00 + phifrac * E;
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom))  
#endif
                sE[threadIdx.x].cmap       += E;
#endif
            PMEFloat dPhi                   = (3.0 * a33 * phifrac + 2.0 * a23) * phifrac + a13;
            dPhi                            = (3.0 * a32 * phifrac + 2.0 * a22) * phifrac + a12 + dPhi * psifrac;
            dPhi                            = (3.0 * a31 * phifrac + 2.0 * a21) * phifrac + a11 + dPhi * psifrac;
            dPhi                            = (3.0 * a30 * phifrac + 2.0 * a20) * phifrac + a10 + dPhi * psifrac;
            PMEFloat dPsi                   = (3.0 * a33 * psifrac + 2.0 * a32) * psifrac + a31;
            dPsi                            = (3.0 * a23 * psifrac + 2.0 * a22) * psifrac + a21 + dPsi * phifrac;
            dPsi                            = (3.0 * a13 * psifrac + 2.0 * a12) * psifrac + a11 + dPsi * phifrac;
            dPsi                            = (3.0 * a03 * psifrac + 2.0 * a02) * psifrac + a01 + dPsi * phifrac;
            dPhi                           *= rad_to_deg_coeff;
            dPsi                           *= rad_to_deg_coeff;
            upxijk                         *= dPhi;
            upyijk                         *= dPhi;
            upzijk                         *= dPhi;
            upxjkl                         *= dPhi;
            upyjkl                         *= dPhi;
            upzjkl                         *= dPhi;
            upxjkl1                        *= dPsi;
            upyjkl1                        *= dPsi;
            upzjkl1                        *= dPsi;
            upxklm                         *= dPsi;
            upyklm                         *= dPsi;
            upzklm                         *= dPsi;
        
      
            // Calculate gradients            
            PMEDouble vx                    =  dotIJKJ * upxijk  + dotLKKJ * upxjkl;
            PMEDouble vy                    =  dotIJKJ * upyijk  + dotLKKJ * upyjkl;
            PMEDouble vz                    =  dotIJKJ * upzijk  + dotLKKJ * upzjkl;
            PMEDouble wx                    = -dotJKLK * upxjkl1 + dotMLLK * upxklm;
            PMEDouble wy                    = -dotJKLK * upyjkl1 + dotMLLK * upyklm;
            PMEDouble wz                    = -dotJKLK * upzjkl1 + dotMLLK * upzklm;
            
#ifdef CHARMM_NEIGHBORLIST
#ifdef CHARMM_VIRIAL
            unsigned long long int iupxijk  = fabs(upxijk) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iupyijk  = fabs(upyijk) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iupzijk  = fabs(upzijk) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iupxjkl  = fabs(upxjkl) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iupyjkl  = fabs(upyjkl) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iupzjkl  = fabs(upzjkl) * FORCESCALE + (PMEDouble)0.5;   
            unsigned long long int iupxjkl1 = fabs(upxjkl1) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iupyjkl1 = fabs(upyjkl1) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iupzjkl1 = fabs(upzjkl1) * FORCESCALE + (PMEDouble)0.5;               
            unsigned long long int iupxklm  = fabs(upxklm) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iupyklm  = fabs(upyklm) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iupzklm  = fabs(upzklm) * FORCESCALE + (PMEDouble)0.5;               
            unsigned long long int ivx      = fabs(vx) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int ivy      = fabs(vy) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int ivz      = fabs(vz) * FORCESCALE + (PMEDouble)0.5;  
            unsigned long long int iwx      = fabs(wx) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iwy      = fabs(wy) * FORCESCALE + (PMEDouble)0.5;
            unsigned long long int iwz      = fabs(wz) * FORCESCALE + (PMEDouble)0.5;                         
            if (upxijk < (PMEDouble)0.0)
                iupxijk                     = 0ull - iupxijk;
            if (upyijk < (PMEDouble)0.0)
                iupyijk                     = 0ull - iupyijk;                    
            if (upzijk < (PMEDouble)0.0)
                iupzijk                     = 0ull - iupzijk;   
            if (upxjkl < (PMEDouble)0.0)
                iupxjkl                     = 0ull - iupxjkl;
            if (upyjkl < (PMEDouble)0.0)
                iupyjkl                     = 0ull - iupyjkl;                    
            if (upzjkl < (PMEDouble)0.0)
                iupzjkl                     = 0ull - iupzjkl;   
            if (upxjkl1 < (PMEDouble)0.0)
                iupxjkl1                    = 0ull - iupxjkl1;
            if (upyjkl1 < (PMEDouble)0.0)
                iupyjkl1                    = 0ull - iupyjkl1;                    
            if (upzjkl1 < (PMEDouble)0.0)
                iupzjkl1                    = 0ull - iupzjkl1;   
            if (upxklm < (PMEDouble)0.0)
                iupxklm                     = 0ull - iupxklm;
            if (upyklm < (PMEDouble)0.0)
                iupyklm                     = 0ull - iupyklm;                    
            if (upzklm < (PMEDouble)0.0)
                iupzklm                     = 0ull - iupzklm;                   
            if (vx < (PMEDouble)0.0)
                ivx                         = 0ull - ivx;
            if (vy < (PMEDouble)0.0)
                ivy                         = 0ull - ivy;                    
            if (vz < (PMEDouble)0.0)
                ivz                         = 0ull - ivz;                                                                   
            if (wx < (PMEDouble)0.0)
                iwx                         = 0ull - iwx;
            if (wy < (PMEDouble)0.0)
                iwy                         = 0ull - iwy;                    
            if (wz < (PMEDouble)0.0)
                iwz                         = 0ull - iwz;     
#endif
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom)) 
#endif
            {
#ifdef CHARMM_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom1.x], 0ull - iupxijk); 
                atomicAdd(&cSim.pUllForceY[atom1.x], 0ull - iupyijk);  
                atomicAdd(&cSim.pUllForceZ[atom1.x], 0ull - iupzijk);       
#else                   
                cSim.pForceXBuffer[atom2.y] =      - upxijk;
                cSim.pForceYBuffer[atom2.y] =      - upyijk;
                cSim.pForceZBuffer[atom2.y] =      - upzijk;
#endif
            }
#ifdef MPI
            if ((atom1.y >= cSim.minLocalAtom) && (atom1.y < cSim.maxLocalAtom)) 
#endif
            {
#ifdef CHARMM_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom1.y], iupxijk - ivx - iupxjkl1);
                atomicAdd(&cSim.pUllForceY[atom1.y], iupyijk - ivy - iupyjkl1);  
                atomicAdd(&cSim.pUllForceZ[atom1.y], iupzijk - ivz - iupzjkl1);       
#else                    
                cSim.pForceXBuffer[atom2.z] =  -vx + upxijk - upxjkl1;
                cSim.pForceYBuffer[atom2.z] =  -vy + upyijk - upyjkl1;
                cSim.pForceZBuffer[atom2.z] =  -vz + upzijk - upzjkl1;
#endif
            }
#ifdef MPI
            if ((atom1.z >= cSim.minLocalAtom) && (atom1.z < cSim.maxLocalAtom)) 
#endif
            {
#ifdef CHARMM_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom1.z], ivx + iupxjkl - iwx + iupxjkl1);
                atomicAdd(&cSim.pUllForceY[atom1.z], ivy + iupyjkl - iwy + iupyjkl1);  
                atomicAdd(&cSim.pUllForceZ[atom1.z], ivz + iupzjkl - iwz + iupzjkl1);       
#else
                cSim.pForceXBuffer[atom2.w] =   vx + upxjkl - wx + upxjkl1;
                cSim.pForceYBuffer[atom2.w] =   vy + upyjkl - wy + upyjkl1;
                cSim.pForceZBuffer[atom2.w] =   vz + upzjkl - wz + upzjkl1;
#endif
            }
#ifdef MPI
            if ((atom1.w >= cSim.minLocalAtom) && (atom1.w < cSim.maxLocalAtom))  
#endif
            {
#ifdef CHARMM_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom1.w], iupxklm + iwx - iupxjkl);
                atomicAdd(&cSim.pUllForceY[atom1.w], iupyklm + iwy - iupyjkl);  
                atomicAdd(&cSim.pUllForceZ[atom1.w], iupzklm + iwz - iupzjkl);       
#else   
                cSim.pForceXBuffer[atom3.x] =       -upxjkl + wx + upxklm;
                cSim.pForceYBuffer[atom3.x] =       -upyjkl + wy + upyklm;
                cSim.pForceZBuffer[atom3.x] =       -upzjkl + wz + upzklm;
#endif
            }
#ifdef MPI
            if ((atom2.x >= cSim.minLocalAtom) && (atom2.x < cSim.maxLocalAtom))  
#endif
            {
#ifdef CHARMM_VIRIAL
                atomicAdd(&cSim.pUllForceX[atom2.x], 0ull - iupxklm); 
                atomicAdd(&cSim.pUllForceY[atom2.x], 0ull - iupyklm);  
                atomicAdd(&cSim.pUllForceZ[atom2.x], 0ull - iupzklm);       
#else            
                cSim.pForceXBuffer[atom3.y] =                    - upxklm;
                cSim.pForceYBuffer[atom3.y] =                    - upyklm;
                cSim.pForceZBuffer[atom3.y] =                    - upzklm;
#endif
            }
#else
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom)) 
#endif
            {
                cSim.pForceXBuffer[atom2.y]+=      - upxijk;
                cSim.pForceYBuffer[atom2.y]+=      - upyijk;
                cSim.pForceZBuffer[atom2.y]+=      - upzijk;
            }
#ifdef MPI
            if ((atom1.y >= cSim.minLocalAtom) && (atom1.y < cSim.maxLocalAtom)) 
#endif            
            {
                cSim.pForceXBuffer[atom2.z]+=  -vx + upxijk - upxjkl1;
                cSim.pForceYBuffer[atom2.z]+=  -vy + upyijk - upyjkl1;
                cSim.pForceZBuffer[atom2.z]+=  -vz + upzijk - upzjkl1;
            }
#ifdef MPI
            if ((atom1.z >= cSim.minLocalAtom) && (atom1.z < cSim.maxLocalAtom)) 
#endif           
            {
                cSim.pForceXBuffer[atom2.w]+=   vx + upxjkl - wx + upxjkl1;
                cSim.pForceYBuffer[atom2.w]+=   vy + upyjkl - wy + upyjkl1;
                cSim.pForceZBuffer[atom2.w]+=   vz + upzjkl - wz + upzjkl1;
            }
#ifdef MPI
            if ((atom1.w >= cSim.minLocalAtom) && (atom1.w < cSim.maxLocalAtom))  
#endif           
            {
                cSim.pForceXBuffer[atom3.x]+=       -upxjkl + wx + upxklm;
                cSim.pForceYBuffer[atom3.x]+=       -upyjkl + wy + upyklm;
                cSim.pForceZBuffer[atom3.x]+=       -upzjkl + wz + upzklm;
            }
#ifdef MPI
            if ((atom2.x >= cSim.minLocalAtom) && (atom2.x < cSim.maxLocalAtom))  
#endif            
            {
                cSim.pForceXBuffer[atom3.y]+=                    - upxklm;
                cSim.pForceYBuffer[atom3.y]+=                    - upyklm;
                cSim.pForceZBuffer[atom3.y]+=                    - upzklm;
            }
#endif                 
        }
        pos                                += cSim.impDihedralOffset + blockDim.x * gridDim.x;
    }


#ifdef CHARMM_ENERGY        
    // Reduce energies
    __syncthreads();
    unsigned int m                          = 1;
    while (m < blockDim.x)
    {
        int p                               = threadIdx.x + m;
        PMEDouble angle_ub                  = ((p < blockDim.x) ? sE[p].angle_ub : (PMEDouble)0.0);
        PMEDouble imp                       = ((p < blockDim.x) ? sE[p].imp      : (PMEDouble)0.0);
        PMEDouble cmap                      = ((p < blockDim.x) ? sE[p].cmap     : (PMEDouble)0.0);
        __syncthreads();
        sE[threadIdx.x].angle_ub           += angle_ub;
        sE[threadIdx.x].imp                += imp;
        sE[threadIdx.x].cmap               += cmap;
        __syncthreads();
        m                                  *= 2;
    }
    unsigned long long int val1             = (unsigned long long int)(fabs(sE[threadIdx.x].angle_ub)  * ENERGYSCALE + (PMEDouble)0.5);
    unsigned long long int val2             = (unsigned long long int)(fabs(sE[threadIdx.x].imp)       * ENERGYSCALE + (PMEDouble)0.5);
    unsigned long long int val3             = (unsigned long long int)(fabs(sE[threadIdx.x].cmap)      * ENERGYSCALE + (PMEDouble)0.5);
    if (sE[threadIdx.x].angle_ub < (PMEDouble)0.0)
        val1                                = 0ull - val1;
    if (sE[threadIdx.x].imp < (PMEDouble)0.0)
        val2                                = 0ull - val2;
    if (sE[threadIdx.x].cmap < (PMEDouble)0.0)
        val3                                = 0ull - val3;    

    if (threadIdx.x == 0)
    {
        atomicAdd(cSim.pEAngle_UB, val1);
        atomicAdd(cSim.pEImp, val2);
        atomicAdd(cSim.pECmap, val3);
    }     
#endif    
}

