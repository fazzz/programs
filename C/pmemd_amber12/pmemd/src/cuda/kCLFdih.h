/***************************************************/
/*                                                 */
/*      AMBER NVIDIA CUDA CPU IMPLEMENTATION       */
/*                 PMEMD VERSION                   */
/*              AMD Dihedral boost                 */
/*                     2011                        */
/*                      by                         */
/*             Romelia Salomon (SDSC)              */
/*                                                 */
/***************************************************/
{
struct Energy {
  PMEDouble dihedral;
};

#if (__CUDA_ARCH__ >= 300)
    __shared__ Energy sE[SM_3X_LOCALFORCES_THREADS_PER_BLOCK];
#elif (__CUDA_ARCH__ >= 200)
    __shared__ Energy sE[SM_2X_LOCALFORCES_THREADS_PER_BLOCK];
#else
    __shared__ Energy sE[SM_13_LOCALFORCES_THREADS_PER_BLOCK];
#endif

    int pos                                 = blockIdx.x * blockDim.x + threadIdx.x;    

    sE[threadIdx.x].dihedral                = (PMEDouble)0.0;

    pos                                    += cSim.bondAngleOffset;
    // Calculate dihedral forces
    while (pos < cSim.dihedralOffset)
    {
        pos                                -= cSim.bondAngleOffset;
        if (pos < cSim.dihedrals)
        {
#ifdef LOCAL_NEIGHBORLIST
            int4 atom1                      = cSim.pImageDihedralID1[pos];
#else
            int4 atom1                      = cSim.pDihedralID1[pos];
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
            faster_sincos2(ap, &sphi, &cphi);

            // Calculate the energy and the derivatives with respect to cosphi
            PMEDouble ct0                   = dihedral1.y * ap;
            PMEDouble sinnp, cosnp;
            faster_sincos2(ct0, &sinnp, &cosnp);
            PMEDouble epw                   = (dihedral2.x + cosnp * dihedral2.y + sinnp * dihedral3) * fzi;
#ifdef MPI
            if ((atom1.x >= cSim.minLocalAtom) && (atom1.x < cSim.maxLocalAtom))        
#endif            
                sE[threadIdx.x].dihedral   += epw;

        }
        pos                                += cSim.bondAngleOffset + blockDim.x * gridDim.x;
    }
   
   
    // Calculate 1-4 forces
    while (pos < cSim.nb14Offset)
    {
        pos                                -= cSim.dihedralOffset;
        if (pos < cSim.nb14s)
        {

        }
        pos                                += cSim.dihedralOffset + blockDim.x * gridDim.x;
    }    
    
    // Calculate Constraint forces
    while (pos < cSim.constraintOffset)
    {
        pos                                -= cSim.nb14Offset;
        if (pos < cSim.constraints)
        {
        }
        pos                                += cSim.nb14Offset + blockDim.x * gridDim.x;
    }
    
    __syncthreads();
    unsigned int m                          = 1;
    while (m < blockDim.x)
    {
        int p                               = threadIdx.x + m;

        PMEDouble dihedral                  = ((p < blockDim.x) ? sE[p].dihedral : (PMEDouble)0.0);
        __syncthreads();

        sE[threadIdx.x].dihedral           += dihedral;

        __syncthreads();
        m                                  *= 2;
    }
    unsigned long long int val3             = llitoulli(lliroundd(sE[threadIdx.x].dihedral  * ENERGYSCALE));
    if (threadIdx.x == 0)
    {
      atomicAdd(cSim.pAMDEDihedral, val3);
    }     
}

