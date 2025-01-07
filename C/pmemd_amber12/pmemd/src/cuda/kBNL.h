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
// #defines: PME_VIRIAL, PME_IS_ORTHOGONAL, PME_ATOMS_PER_WARP 
#define uint unsigned int

struct BNLAtom
{
    float x;
    float y;
    float z;
};

struct BNLWarp
{
#if (__CUDA_ARCH__ < 300)
    uint2 homeCell;
    uint2 cell[NEIGHBORCELLS];
    NLRecord nlRecord;
    uint nlAtoms;
    uint nlpos;
    uint pos;
    uint cellID;
#endif
    uint atomList[GRID];
    uint exclusionMask[GRID];
    NLEntry nlEntry;    
    uint minExclusion;
    uint maxExclusion;
    uint offset;
#if (PME_ATOMS_PER_WARP == 32) && (__CUDA_ARCH__ < 300)
    float x;
    float y;
    float z;
#endif
};

#if (__CUDA_ARCH__ >= 300)
#define POS __shfl(shNlRecord, NEIGHBORCELLS + 2)
#define CELLID shCellID
#define NLPOS __shfl(shNlRecord, NEIGHBORCELLS + 3)
#define NLATOMS shNlatoms
#define MINEXCLUSION psWarp->minExclusion
#define MAXEXCLUSION psWarp->maxExclusion
#define HOMECELLX __shfl(shNlRecord, NEIGHBORCELLS + 4)
#define HOMECELLY __shfl(shNlRecord, NEIGHBORCELLS + 5)
#define CELLX(i) __shfl(shCell.x, i)
#define CELLY(i) __shfl(shCell.y, i)
#define NLRECORDNEIGHBORCELL(i) __shfl(shNlRecord, i)
#define NLRECORDHOMECELL __shfl(shNlRecord, NEIGHBORCELLS)
#define NLRECORDNEIGHBORCELLS __shfl(shNlRecord, NEIGHBORCELLS + 1)
#else
#define POS psWarp->pos
#define CELLID psWarp->cellID
#define NLPOS psWarp->nlpos
#define NLATOMS psWarp->nlAtoms
#define MINEXCLUSION psWarp->minExclusion
#define MAXEXCLUSION psWarp->maxExclusion
#define HOMECELLX psWarp->homeCell.x
#define HOMECELLY psWarp->homeCell.y
#define CELLX(i) psWarp->cell[i].x
#define CELLY(i) psWarp->cell[i].y
#define NLRECORDNEIGHBORCELL(i) psWarp->nlRecord.NL.neighborCell[i]
#define NLRECORDHOMECELL psWarp->nlRecord.NL.homeCell
#define NLRECORDNEIGHBORCELLS psWarp->nlRecord.NL.neighborCells
#endif


#if (__CUDA_ARCH__ >= 300)
uint2 shCell;
uint shCellID;
uint shNlatoms;
volatile uint shNlRecord;
#if (PME_ATOMS_PER_WARP == 32)
const int THREADS_PER_BLOCK = SM_3X_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK;
#else
const int THREADS_PER_BLOCK = SM_3X_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK;
#endif
#elif (__CUDA_ARCH__ >= 200)
#if (PME_ATOMS_PER_WARP == 32)
const int THREADS_PER_BLOCK = SM_2X_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK;
#else
const int THREADS_PER_BLOCK = SM_2X_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK;
#endif
#else
#if (PME_ATOMS_PER_WARP == 32)
const int THREADS_PER_BLOCK = SM_13_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK;
#else
const int THREADS_PER_BLOCK = SM_13_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK;
#endif
#endif
__shared__ volatile BNLWarp sNLWarp[THREADS_PER_BLOCK / GRID];
#if (PME_ATOMS_PER_WARP == 16) && (__CUDA_ARCH__ < 300)
__shared__ volatile BNLAtom sAtom[THREADS_PER_BLOCK];
#endif
#ifdef PME_VIRIAL
__shared__ float sCutPlusSkin2;
__shared__ float sUcellf[9];  
#endif
#if (__CUDA_ARCH__ < 300)   
extern __shared__ uint sExclusion[];
#endif
#ifdef PME_VIRIAL 
    if (threadIdx.x < 9)
        sUcellf[threadIdx.x]                                            = cSim.pNTPData->ucellf[threadIdx.x];
    if (threadIdx.x == 32)
        sCutPlusSkin2                                                   = cSim.pNTPData->cutPlusSkin2;                     
    __syncthreads();   
    float cutPlusSkin2                                                  = sCutPlusSkin2;
#else
    float cutPlusSkin2                                                  = cSim.cutPlusSkin2;
#endif
    unsigned int warp                                                   = threadIdx.x >> GRIDBITS;
    volatile BNLWarp* psWarp                                            = &sNLWarp[warp];
#if (__CUDA_ARCH__ >= 300)
    unsigned int globalWarp                                             = warp + ((blockIdx.x * blockDim.x) >> GRIDBITS);        
    uint* psExclusion                                                   = &cSim.pBNLExclusionBuffer[globalWarp * 768];
#else
    volatile uint* psExclusion                                          = &sExclusion[warp * cSim.NLMaxExclusionsPerWarp];
#endif
#if (__CUDA_ARCH__ >= 300)
    if ((threadIdx.x & GRIDBITSMASK) == (NEIGHBORCELLS + 2))
        shNlRecord                                                      = (blockIdx.x * blockDim.x + threadIdx.x) >> GRIDBITS;
#else
    POS                                                                 = (blockIdx.x * blockDim.x + threadIdx.x) >> GRIDBITS; 
#endif
	unsigned int exclusionMask                                          = cSim.NLAtomsPerWarpMask >> (threadIdx.x & cSim.NLAtomsPerWarpBitsMask);
#if (PME_ATOMS_PER_WARP == 16)
    exclusionMask                                                       = exclusionMask |
                                                                         (exclusionMask << 16);
#endif

#if 0
    if ((blockIdx.x == 0) && (threadIdx.x == 0))
    {
        printf("0x%08x 0x%08x 0x%08x 0x%08x\n", &(psWarp->nlRecord), &(psWarp->nlRecord.NL.homeCell), &(psWarp->nlRecord.NL.neighborCells), &(psWarp->nlRecord.NL.neighborCell));

    }
#endif
  
    while (POS < cSim.NLSize)
    {
        unsigned int tgx                                                = threadIdx.x & GRIDBITSMASK;       
        
        // Read NLRecord information
#if (__CUDA_ARCH__ >= 300)
        unsigned int pos1                                               = POS;
#endif
        if (tgx < 16)
#if (__CUDA_ARCH__ >= 300)
            shNlRecord                                                  = cSim.pNLRecord[pos1].array[tgx];
#else
            psWarp->nlRecord.array[tgx]                                 = cSim.pNLRecord[POS].array[tgx];
#endif
        
        // Calculate Exclusion/neighbor list space required
        int atomOffset                                                  = NLRECORDNEIGHBORCELLS & NLATOMOFFSETMASK;

#if (__CUDA_ARCH__ >= 300)
        uint2 homeCell                                                  = cSim.pNLNonbondCellStartEnd[NLRECORDHOMECELL]; 
        if (tgx == (NEIGHBORCELLS + 4))
            shNlRecord                                                  = homeCell.x;
        if (tgx == (NEIGHBORCELLS + 5))
            shNlRecord                                                  = homeCell.y;
#else
        if (tgx == 0)
        {
            uint2 homeCell                                              = cSim.pNLNonbondCellStartEnd[NLRECORDHOMECELL];   
            psWarp->homeCell.x                                          = homeCell.x;
            psWarp->homeCell.y                                          = homeCell.y;
        }
#endif
        psWarp->nlEntry.NL.xyBufferOffset                               = (((NLRECORDNEIGHBORCELLS >> NLXCELLOFFSETSHIFT) & NLXCELLOFFSETMASK) << NLENTRYXBUFFEROFFSETSHIFT) | atomOffset;
        if ((NLRECORDNEIGHBORCELL(0) & NLCELLTYPEMASK) == 0)
            psWarp->nlEntry.NL.xyBufferOffset                          |= NLENTRYHOMECELLMASK;
        
        int ysize                                                       = max(0, (int)(HOMECELLY - HOMECELLX - atomOffset * cSim.NLAtomsPerWarp));
        if (ysize > 0)
            ysize                                                       = 1 + max(0, ysize - 1) / (cSim.NLAtomsPerWarp * cSim.NLYDivisor);

            
        // Calculate maximum required space
#if (PME_ATOMS_PER_WARP == 16) && (__CUDA_ARCH__ < 300)
        volatile BNLAtom* psAtom                                        = &sAtom[threadIdx.x & GRIDPADDINGMASK];
#endif
        unsigned int cells                                              = NLRECORDNEIGHBORCELLS >> NLCELLCOUNTSHIFT;
        if (tgx < cells)
        {
            uint2 cell                                                  = cSim.pNLNonbondCellStartEnd[NLRECORDNEIGHBORCELL(tgx) >> NLCELLSHIFT];
#if (__CUDA_ARCH__ >= 300)
            shCell.x                                                    = cell.x;
            shCell.y                                                    = cell.y;
#else
            psWarp->cell[tgx].x                                         = cell.x;
            psWarp->cell[tgx].y                                         = cell.y;      
#endif
            psWarp->atomList[tgx]                                       = cell.y - cell.x;
        }
        else
            psWarp->atomList[tgx]                                       = 0;
        
        // Reduce xsize down to thread 0 of each warp
#if (NEIGHBORCELLS > 16)
        if (tgx < 16)
            psWarp->atomList[tgx]                                      += psWarp->atomList[tgx + 16];
#endif        
        if (tgx < 8)
            psWarp->atomList[tgx]                                      += psWarp->atomList[tgx + 8];
        if (tgx < 4)
            psWarp->atomList[tgx]                                      += psWarp->atomList[tgx + 4];         
        if (tgx < 2)
            psWarp->atomList[tgx]                                      += psWarp->atomList[tgx + 2];         
        if (tgx < 1)
            psWarp->atomList[tgx]                                      += psWarp->atomList[tgx + 1];                 
        
#if (__CUDA_ARCH__ >= 300)
        if (tgx == (NEIGHBORCELLS + 2))
#else
        if (tgx == 0)
#endif
        {
            uint totalXSize                                             = ((psWarp->atomList[0] + GRID - 1) >> GRIDBITS);
            uint offset                                                 = atomicAdd(cSim.pNLTotalOffset, totalXSize * ysize * cSim.NLOffsetPerWarp + cSim.NLAtomsPerWarp);
            cSim.pNLOffset[POS]                                         = offset;  
            psWarp->offset                                              = offset;       
        }       
             
        // Generate actual neighbor list entry
        uint ypos                                                       = HOMECELLX + (NLRECORDNEIGHBORCELLS & NLATOMOFFSETMASK) * cSim.NLAtomsPerWarp;
        psWarp->nlEntry.NL.yAtom                                        = ypos;        
        psWarp->nlEntry.NL.yEnd                                         = (HOMECELLY << NLCELLSHIFT) | (NLRECORDHOMECELL & NLCELLTYPEMASK);
#if (__CUDA_ARCH__ >= 300)
        if (tgx == (NEIGHBORCELLS + 3))
            shNlRecord                                                  = 0;
#else
        NLPOS                                                           = 0;
#endif  

        while (ypos < HOMECELLY)
        {
        
            // Calculate y bounds
            uint ymax                                                   = min(ypos + cSim.NLAtomsPerWarp, HOMECELLY);
        
            // Read y atoms
            float xi;
            float yi;
            float zi;
            unsigned int index                                          = ypos + (tgx & (cSim.NLAtomsPerWarpBitsMask));            
            if (index < ymax)
            {         
                PMEFloat2 xy                                            = cSim.pAtomXYSP[index];                     
                zi                                                      = cSim.pAtomZSP[index];                       
                xi                                                      = xy.x;
                yi                                                      = xy.y;
            }
            else
            {
                xi                                                      = (float)10000.0 * index;
                yi                                                      = (float)10000.0 * index;
                zi                                                      = (float)10000.0 * index;
            }       
        
#ifndef PME_IS_ORTHOGONAL   
            // Transform into cartesian space     
#ifdef PME_VIRIAL
            xi                                                          = sUcellf[0] * xi + sUcellf[1] * yi + sUcellf[2] * zi;
            yi                                                          =                   sUcellf[4] * yi + sUcellf[5] * zi;
            zi                                                          =                                     sUcellf[8] * zi;              
#else                
            xi                                                          = cSim.ucell[0][0] * xi + cSim.ucell[0][1] * yi + cSim.ucell[0][2] * zi;
            yi                                                          =                         cSim.ucell[1][1] * yi + cSim.ucell[1][2] * zi;
            zi                                                          =                                                 cSim.ucell[2][2] * zi;
#endif
#endif  

            // Test to see if exclusions fit in shared memory
            if (tgx < ymax - ypos)
            {
                uint atom                                               = cSim.pImageAtom[index];
                uint2 exclusionStartCount                               = cSim.pNLExclusionStartCount[atom];
                psWarp->atomList[tgx]                                   = exclusionStartCount.x;
                psWarp->exclusionMask[tgx]                              = exclusionStartCount.y;
            }
            else
            {
                psWarp->atomList[tgx]                                   = 0;
                psWarp->exclusionMask[tgx]                              = 0;
            }
            uint tempExclusions                                         = psWarp->exclusionMask[tgx];    
    
            
            // Reduce exclusion count
#if (__CUDA_ARCH__ >= 300)
            uint totalExclusions                                        = tempExclusions;
            totalExclusions                                            += __shfl(totalExclusions, tgx + 1);
            totalExclusions                                            += __shfl(totalExclusions, tgx + 2);
            totalExclusions                                            += __shfl(totalExclusions, tgx + 4);
            totalExclusions                                            += __shfl(totalExclusions, tgx + 8);
            totalExclusions                                            += __shfl(totalExclusions, tgx + 16);

            // Exclusions are stored in L1 on SM 3.0 and up so no need to worry about max Exclusions
            {
                totalExclusions                                        = 0; 
#else

            if (tgx < 16)
                psWarp->exclusionMask[tgx]                             += psWarp->exclusionMask[tgx + 16];
            if (tgx < 8)
                psWarp->exclusionMask[tgx]                             += psWarp->exclusionMask[tgx + 8];   
            if (tgx < 4)
                psWarp->exclusionMask[tgx]                             += psWarp->exclusionMask[tgx + 4];
            if (tgx < 2)
                psWarp->exclusionMask[tgx]                             += psWarp->exclusionMask[tgx + 2];
            if (tgx < 1)
                psWarp->exclusionMask[tgx]                             += psWarp->exclusionMask[tgx + 1]; 
            uint totalExclusions                                        = 0; 
            // Cache exclusions if there's room and offset them relative to ypos so that the atoms of adjacent
            // cells increase monotonically (massive reduction in exclusion tests results)]
            if (psWarp->exclusionMask[0] <= cSim.NLMaxExclusionsPerWarp)
            {
#endif
                psWarp->exclusionMask[tgx]                              = tempExclusions;
                unsigned int minExclusion                               = cSim.atoms;
                unsigned int maxExclusion                               = 0;
                uint limit                                              = ymax - ypos;
                for (int i = 0; i < limit; i++)
                {
                    uint start                                          = psWarp->atomList[i];
                    uint count                                          = psWarp->exclusionMask[i];                
                    for (int j = tgx; j < count; j += GRID)
                    {
                        int atom                                        = cSim.pNLExclusionList[start + j];
                        int imageAtom                                   = cSim.pImageAtomLookup[atom];
                        minExclusion                                    = min(minExclusion, imageAtom);
                        maxExclusion                                    = max(maxExclusion, imageAtom);
                        psExclusion[totalExclusions + j]                = (imageAtom << NLEXCLUSIONSHIFT) | i;                        
                    }
                    totalExclusions                                    += count;
                }
#if (__CUDA_ARCH__ >= 500) // Why is this so slow?
                MINEXCLUSION                                            = minExclusion;
                MAXEXCLUSION                                            = maxExclusion;
                MINEXCLUSION                                            = min(__shfl(MINEXCLUSION, tgx + 1), MINEXCLUSION);
                MINEXCLUSION                                            = min(__shfl(MINEXCLUSION, tgx + 2), MINEXCLUSION);
                MINEXCLUSION                                            = min(__shfl(MINEXCLUSION, tgx + 4), MINEXCLUSION);
                MINEXCLUSION                                            = min(__shfl(MINEXCLUSION, tgx + 8), MINEXCLUSION);
                MINEXCLUSION                                            = min(__shfl(MINEXCLUSION, tgx + 16), MINEXCLUSION);
                MAXEXCLUSION                                            = max(__shfl(MAXEXCLUSION, tgx + 1), MAXEXCLUSION);
                MAXEXCLUSION                                            = max(__shfl(MAXEXCLUSION, tgx + 2), MAXEXCLUSION);
                MAXEXCLUSION                                            = max(__shfl(MAXEXCLUSION, tgx + 4), MAXEXCLUSION);
                MAXEXCLUSION                                            = max(__shfl(MAXEXCLUSION, tgx + 8), MAXEXCLUSION);
                MAXEXCLUSION                                            = max(__shfl(MAXEXCLUSION, tgx + 16), MAXEXCLUSION);
#else                
                // Determine min and max exclusions
                psWarp->exclusionMask[tgx]                              = minExclusion;
                psWarp->atomList[tgx]                                   = maxExclusion;
             
                if (tgx < 16)
                {
                    psWarp->exclusionMask[tgx]                          = min(psWarp->exclusionMask[tgx], psWarp->exclusionMask[tgx + 16]);
                    psWarp->atomList[tgx]                               = max(psWarp->atomList[tgx], psWarp->atomList[tgx + 16]);
                }
                if (tgx < 8)
                {
                    psWarp->exclusionMask[tgx]                          = min(psWarp->exclusionMask[tgx], psWarp->exclusionMask[tgx + 8]);
                    psWarp->atomList[tgx]                               = max(psWarp->atomList[tgx], psWarp->atomList[tgx + 8]);
                }
                if (tgx < 4)
                {
                    psWarp->exclusionMask[tgx]                          = min(psWarp->exclusionMask[tgx], psWarp->exclusionMask[tgx + 4]);
                    psWarp->atomList[tgx]                               = max(psWarp->atomList[tgx], psWarp->atomList[tgx + 4]);
                }
                if (tgx < 2)
                {
                    psWarp->exclusionMask[tgx]                          = min(psWarp->exclusionMask[tgx], psWarp->exclusionMask[tgx + 2]);
                    psWarp->atomList[tgx]                               = max(psWarp->atomList[tgx], psWarp->atomList[tgx + 2]);
                }                
                MINEXCLUSION                                            = min(psWarp->exclusionMask[0], psWarp->exclusionMask[1]);
                MAXEXCLUSION                                            = max(psWarp->atomList[0], psWarp->atomList[1]);
#endif
            }
            //else
            //    printf("blown %d\n", psWarp->exclusionMask[0]);

#if 0 
            if (tgx == 1)
                printf("%06d %06d %06d %06d %06d %06d\n", POS, psWarp->ypos, NLPOS, blockIdx.x * blockDim.x + threadIdx.x, MINEXCLUSION, MAXECLUSION);
#endif
            
#if 0
            if (tgx == 0)
                printf("A %06d %06d %06d %06d %06d\n", POS, psWarp->ypos, psWarp->ymax, HOMECELLY, totalExclusions);
#endif  
         
            // Initialize Neighbor List variables for current line of entry
            psWarp->exclusionMask[tgx]                                  = 0;       
            unsigned int cpos                                           = 0;
            unsigned int atoms                                          = 0;
            NLATOMS                                                     = 0;
            uint minAtom                                                = cSim.atoms;
            uint maxAtom                                                = 0;
            unsigned int threadmask                                     = 1 << tgx; 
                                                       
            while (cpos < cells)
            {
                // Check for home cell
                CELLID                                                  = NLRECORDNEIGHBORCELL(cpos) & NLCELLTYPEMASK;
                uint xpos;
                
                // Cell 0 always starts along force matrix diagonal
                if ((cpos == 0) && (psWarp->nlEntry.NL.xyBufferOffset & NLENTRYHOMECELLMASK))
                {                    
                    // Calculate exclusions assuming all atoms are in range of each other
                    for (int i = tgx; i < totalExclusions; i += GRID)
                    {
                        uint atom                                       = psExclusion[i] >> NLEXCLUSIONSHIFT;
                        if ((atom >= ypos) && (atom < ymax))
                        {      
                            unsigned int pos                            = atom - ypos;
                            do
                            {
                                psWarp->atomList[pos]                   = tgx;
                                if (psWarp->atomList[pos] == tgx)
                                    psWarp->exclusionMask[pos]         |= (unsigned int)(1 << (psExclusion[i] & NLEXCLUSIONATOMMASK));                  
                            }
                            while (psWarp->atomList[pos] != tgx);
#if 0                            
                            unsigned int a1             = cSim.pImageAtom[atom];
                            unsigned int a2             = cSim.pImageAtom[psWarp->ypos + (psExclusion[i] & NLEXCLUSIONATOMMASK)];
                            for (int ii = 0; ii < GRID; ii++)
                            {
                                if (ii == tgx)
                                    printf("Ex %06d %06d\n", min(a1, a2), max(a1, a2));
                            }   
#endif                          
                        }
                    }                        

                    // Output exclusion masks
                    if (tgx < cSim.NLAtomsPerWarp) 
                    {
                        unsigned int mask                               = psWarp->exclusionMask[tgx];
                        mask                                            = ((mask >> (1 + tgx)) | (mask << (cSim.NLAtomsPerWarp - tgx - 1))) & cSim.NLAtomsPerWarpMask;
#if 0                        
                        for (int i = 0; i < cSim.NLAtomsPerWarp; i++)
                        if (tgx == i)
                            printf("Ex %06d %06d %06d %06d %1d%1d%1d%1d%1d%1d%1d%1d\n", psWarp->ypos, i, 1 + tgx, cSim.NLAtomsPerWarp - tgx - 1,            
                            (mask&0x1) != 0, (mask&0x2) != 0, (mask&0x4) != 0, (mask &0x8) != 0,
                            (mask&0x10) != 0, (mask&0x20) != 0, (mask&0x40) != 0, (mask &0x80) != 0);
#endif                            
                        cSim.pNLAtomList[psWarp->offset + tgx]          = mask;    
                        psWarp->exclusionMask[tgx]                      = 0;                   

                    }
                    psWarp->offset                                     += cSim.NLAtomsPerWarp;
                    xpos                                                = ypos + cSim.NLAtomsPerWarp;
                   
#if 0       
                    if (tgx == 0)
                    {
                        for (int k = 0; k < cSim.NLAtomsPerWarp; k++)
                        {             
                            printf("%06d %06d 0x%08x\n", psWarp->ypos, k, psWarp->exclusionMask[k]);
                        }
                    }
#endif                       
                }                    
                else
                    xpos                                                = CELLX(cpos);                    
 
                // Read x atoms
                while (xpos < CELLY(cpos))
                {
                    // Calculate number of atoms in this iteration
                    uint xmax                                           = min(xpos + GRID, CELLY(cpos)) - xpos; 

#if 0                      
                    if (tgx == 0)
                        printf("%6d %6d %6d %6d %6d\n", POS, psWarp->ypos, psWarp->ymax, xpos, xmax);
#endif                                    

#if (PME_ATOMS_PER_WARP == 32) || (__CUDA_ARCH__ >= 300)
#define SATOMX sAtomx
#define SATOMY sAtomy
#define SATOMZ sAtomz
                    float sAtomx;
                    float sAtomy;
                    float sAtomz;
#else
#define SATOMX psAtom[tgx].x
#define SATOMY psAtom[tgx].y
#define SATOMZ psAtom[tgx].z
#endif

                    // Read up to GRID atoms
                    if (tgx < xmax)
                    {
                        PMEFloat2 xy                                    = cSim.pAtomXYSP[xpos + tgx];  
                        SATOMZ                                          = cSim.pAtomZSP[xpos + tgx];                      
                        SATOMX                                          = xy.x;  
                        SATOMY                                          = xy.y;                             
                    }
                    else
                    {
                        SATOMX                                          = (float)-10000.0 * tgx;
                        SATOMY                                          = (float)-10000.0 * tgx;
                        SATOMZ                                          = (float)-10000.0 * tgx;                        
                    }
                        
                    // Translate all atoms into a local coordinate system within one unit
                    // cell of the first atom read to avoid PBC handling within inner loops                
                    unsigned int cellID                                 = CELLID;
#if defined(PME_VIRIAL) && defined(PME_IS_ORTHOGONAL)
                    SATOMX                                             += sUcellf[0] * cSim.cellOffset[cellID][0];
                    SATOMY                                             += sUcellf[4] * cSim.cellOffset[cellID][1];
                    SATOMZ                                             += sUcellf[8] * cSim.cellOffset[cellID][2];
#else
                    SATOMX                                             += cSim.cellOffset[cellID][0];
                    SATOMY                                             += cSim.cellOffset[cellID][1];
                    SATOMZ                                             += cSim.cellOffset[cellID][2];
#endif

#ifndef PME_IS_ORTHOGONAL
#ifdef PME_VIRIAL
                    SATOMX                                              = sUcellf[0] * SATOMX + sUcellf[1] * SATOMY + sUcellf[2] * SATOMZ;
                    SATOMY                                              =                       sUcellf[4] * SATOMY + sUcellf[5] * SATOMZ;
                    SATOMZ                                              =                                             sUcellf[8] * SATOMZ;
#else
                    SATOMX                                              = cSim.ucell[0][0] * SATOMX + cSim.ucell[0][1] * SATOMY + cSim.ucell[0][2] * SATOMZ;
                    SATOMY                                              =                             cSim.ucell[1][1] * SATOMY + cSim.ucell[1][2] * SATOMZ;
                    SATOMZ                                              =                                                         cSim.ucell[2][2] * SATOMZ;
#endif
#endif

#undef SATOMX
#undef SATOMY
#undef SATOMZ                     
                   
                    // Check each atom to see if it's in range 
#if (PME_ATOMS_PER_WARP == 32)
                    for (int i = 0; i < xmax; i++)
#else
                    unsigned int offset                                 = tgx >> cSim.NLAtomsPerWarpBits;                   
                    for (int i = 0; i < xmax; i += 2)
#endif
                    {
                        // Calculate distance
#if (PME_ATOMS_PER_WARP == 32)
                        unsigned int pos                                = i;
#else
                        unsigned int pos                                = i + offset;
#endif

#if (__CUDA_ARCH__ >= 300)
                        float ax                                        = __shfl(sAtomx, pos);
                        float ay                                        = __shfl(sAtomy, pos);
                        float az                                        = __shfl(sAtomz, pos);
#endif

                        int pred = 0;
                        if (pos < xmax)
                        {
#if __CUDA_ARCH__ >= 300
                            float dx                                    = xi - ax;
                            float dy                                    = yi - ay;
                            float dz                                    = zi - az;
#elif (PME_ATOMS_PER_WARP == 32)
                            if (pos == tgx)
                            {
                                psWarp->x                               = sAtomx;
                                psWarp->y                               = sAtomy;
                                psWarp->z                               = sAtomz;
                            }
                            float dx                                    = xi - psWarp->x;
                            float dy                                    = yi - psWarp->y;
                            float dz                                    = zi - psWarp->z;
#else
                            float dx                                    = xi - psAtom[pos].x;
                            float dy                                    = yi - psAtom[pos].y;
                            float dz                                    = zi - psAtom[pos].z;
#endif
                            float r2                                    = dx * dx + dy * dy + dz * dz;
                            pred                                        = (r2 < cutPlusSkin2);
                        }
                        
                        // Test if atom is within range of any other atom 
                        pos                                             = i;
                        unsigned int mask                               = cSim.NLAtomsPerWarpMask;
#if (__CUDA_ARCH__ >= 200)
                        unsigned int vote                               = __ballot(pred);
                        while (vote)
                        {          
                            if (vote & mask)
                            {
                                vote                                   &= ~mask;
#else              
                        while (__any(pred))
                        {                   

                            if (__any(pred && (threadmask & mask)))
                            {
                                // Clear predicate bits
                                if (threadmask & mask)
                                    pred                                = 0;     
#endif                 
                                // Add atom to growing neighbor list
                                unsigned int atom                       = xpos + pos;
                                minAtom                                 = min(minAtom, atom);
                                maxAtom                                 = max(maxAtom, atom);
                                psWarp->atomList[atoms]                 = (atom << NLCELLSHIFT) | CELLID;
                                atoms++; 
                                

                                // Output GRID atoms if ready
                                if (atoms == GRID)
                                {
                                    cSim.pNLAtomList[psWarp->offset + tgx]  
                                                                        = psWarp->atomList[tgx];
                                    NLATOMS                            += GRID;
                                    psWarp->offset                     += GRID; 
                                                                    
                                    // Search for y atom exclusions matching any x atom all at once (this should
                                    // reduce exclusion tests by a factor of approximately 100 overall 
                                    // But first rule out skipping exclusion test
                                  
                                    if ((minAtom <= MAXEXCLUSION) && (maxAtom >= MINEXCLUSION))    
                                    { 
                                        uint atom                       = (psWarp->atomList[tgx] >> NLCELLSHIFT);                          
#if (__CUDA_ARCH__ >= 200)
                                        for (int j = 0; j < totalExclusions; j += GRID)
                                        {
                                            int offset                  = j + tgx;
                                            unsigned int eatom          = 0xffffffff;
                                            if (offset < totalExclusions)
                                                eatom                   = psExclusion[offset] >> NLEXCLUSIONSHIFT;
                                            unsigned int vote           = __ballot((eatom >= minAtom) && (eatom <= maxAtom));
                                            while (vote)
                                            {
                                                unsigned int k          = __ffs(vote) - 1;
                                                offset                  = k + j;
                                                eatom                   = psExclusion[offset] >> NLEXCLUSIONSHIFT;
                                                if (atom == eatom)
                                                {
                                                    psWarp->exclusionMask[psExclusion[offset] & NLEXCLUSIONATOMMASK]
                                                                       |= threadmask;
                                                }
                                                vote                   ^= 1 << k;
                                            }
                                        }
#else                                        
                                        for (int j = 0; j < totalExclusions; j++)
                                        {       
                                            if ((psExclusion[j] >> NLEXCLUSIONSHIFT) == atom)
                                            {                                           
                                                atomicOr((unsigned int*)(&psWarp->exclusionMask[psExclusion[j] & NLEXCLUSIONATOMMASK]), threadmask);
#if 0                      
                                                unsigned int a1 = cSim.pImageAtom[psWarp->ypos +  (psExclusion[j] & NLEXCLUSIONATOMMASK)];
                                                unsigned int a2 = cSim.pImageAtom[atom];
                                                for (int i = 0; i < GRID; i++)
                                                    if (tgx == i)
                                                        printf("Ex %06d %06d\n", min(a1, a2), max(a1, a2));
#endif                                                                       
                                            }
                                        } 
#endif
                                    }             
                         
                                    // Output exclusion masks
                                    if (tgx < cSim.NLAtomsPerWarp)
                                    { 
                                        unsigned int emask              = psWarp->exclusionMask[tgx];
                                        cSim.pNLAtomList[psWarp->offset + tgx]
                                                                        = ((emask >> tgx) & exclusionMask) |
                                                                          ((emask << (cSim.NLAtomsPerWarp - tgx)) & ~exclusionMask);                                                                                                                                                                               
                                    }
                                    psWarp->offset                     += cSim.NLAtomsPerWarp;
                                    atoms                               = 0;
                                    psWarp->exclusionMask[tgx]          = 0;
                                    psWarp->atomList[tgx]               = 0;  
                                    minAtom                             = cSim.atoms;
                                    maxAtom                             = 0;                                             
                                }
                            }
                            mask                                      <<= cSim.NLAtomsPerWarp;
                            pos++;
                        }
                    }                    
                    xpos                                               += GRID;                    
                }               

                // Move to next cell           
                cpos++;
            }

            // Output last batch of atoms for this swath
            if (atoms > 0)
            {
                cSim.pNLAtomList[psWarp->offset + tgx]                  = psWarp->atomList[tgx];
                psWarp->offset                                         += GRID; 
                               
                // Search for y atom exclusions matching any x atom all at once (this should
                // reduce exclusion tests by a factor of approximately 100 overall              
                if ((tgx < atoms) && ((minAtom <= MAXEXCLUSION) && (maxAtom >= MINEXCLUSION)))    
                {                                           
                    unsigned int atom                                   = (psWarp->atomList[tgx] >> NLCELLSHIFT);
                    for (int j = 0; j < totalExclusions; j++)
                    {       
                        if ((psExclusion[j] >> NLEXCLUSIONSHIFT) == atom)
                        {                    
                            atomicOr((unsigned int*)(&psWarp->exclusionMask[psExclusion[j] & NLEXCLUSIONATOMMASK]), threadmask);
#if 0                      
                        unsigned int a1 = cSim.pImageAtom[psWarp->ypos +  (psExclusion[j] & NLEXCLUSIONATOMMASK)];
                        unsigned int a2 = cSim.pImageAtom[atom];
                        for (int i = 0; i < GRID; i++)
                        if (tgx == i)
                            printf("Ex %06d %06d\n", min(a1, a2), max(a1, a2));
#endif                       
                        }
                    }
                }                                       

                // Output exclusion masks
                if (tgx < cSim.NLAtomsPerWarp)
                { 
                    unsigned int emask                                  = psWarp->exclusionMask[tgx];
                    cSim.pNLAtomList[psWarp->offset + tgx]
                                                                        = ((emask >> tgx) & exclusionMask) |
                                                                          ((emask << (cSim.NLAtomsPerWarp - tgx)) & ~exclusionMask);                                                                                                                                                                               
                }
                NLATOMS                                                += atoms;
                psWarp->offset                                         += cSim.NLAtomsPerWarp;
                atoms                                                   = 0;
                psWarp->exclusionMask[tgx]                              = 0;
                psWarp->atomList[tgx]                                   = 0;           
            }
            psWarp->nlEntry.NL.xAtoms[NLPOS]                            = NLATOMS;
            
#if 0                      
            unsigned int total = 0;
            for (int i = 0; i < NLRECORDNEIGHBORCELLS >> NLCELLCOUNTSHIFT; i++)
                total += ((CELLY(tgx) - CELLX(tgx) + GRID - 1) >> GRIDBITS) << GRIDBITS;
            if (tgx == 0)
                printf("%06d,%06d,%06d\n", psWarp->ypos, ((NLATOMS + GRID - 1) >> GRIDBITS) << GRIDBITS, total);
#endif             

            // Advance to next swath of atoms
            ypos                                                       += cSim.NLYDivisor * cSim.NLAtomsPerWarp;
#if (__CUDA_ARCH__ >= 300)
            if (tgx == (NEIGHBORCELLS + 3))
                shNlRecord++;
#else
            NLPOS++;
#endif
        }

        // Output neighbor list entry
        cSim.pNLEntry[POS].array[tgx]                                   = psWarp->nlEntry.array[tgx];
        
        //if (tgx == 0)
        //    printf("%06d E %08d\n", POS, psWarp->offset - cSim.pNLOffset[POS]);

        // Advance to next NLRecord entry 
#if (__CUDA_ARCH__ >= 300)
        if (tgx == (NEIGHBORCELLS + 2))
            shNlRecord                                                  = atomicAdd(cSim.pNLPosition, 1);
#else   
        if (tgx == 0)
            POS                                                         = atomicAdd(cSim.pNLPosition, 1);
#endif
    }
}
      
