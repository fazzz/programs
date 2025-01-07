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
    uint2 homeCell;
    uint2 cell[NEIGHBORCELLS];
    uint atomList[GRID];
    uint exclusionMask[GRID];
    NLRecord nlRecord;
    NLEntry nlEntry;    
    uint pos;
    uint cellID;
    uint nlAtoms;
    uint nlpos;
    uint offset;
    int minExclusion;
    int maxExclusion;
#if (PME_ATOMS_PER_WARP == 32)
    float x;
    float y;
    float z;
#endif
};

#if (__CUDA_ARCH__ >= 200)
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
#if (PME_ATOMS_PER_WARP == 16)
__shared__ volatile BNLAtom sAtom[THREADS_PER_BLOCK];
#endif
#ifdef PME_VIRIAL
__shared__ float sCutPlusSkin2;
__shared__ float sUcellf[9];  
#endif    
extern __shared__ uint sExclusion[];
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
    volatile uint* psExclusion                                          = &sExclusion[warp * cSim.NLMaxExclusionsPerWarp];
    psWarp->pos                                                         = (blockIdx.x * blockDim.x + threadIdx.x) >> GRIDBITS; 
	unsigned int exclusionMask                                          = cSim.NLAtomsPerWarpMask >> (threadIdx.x & cSim.NLAtomsPerWarpBitsMask);
#if (PME_ATOMS_PER_WARP == 16)
    exclusionMask                                                       = exclusionMask |
                                                                         (exclusionMask << 16);
#endif
#if 0
    if ((blockIdx.x == 0) && (threadIdx.x < GRID))
        printf("%06d 0x%08x\n", threadIdx.x, exclusionMask);
    return;                                                                      
#endif    
    while (psWarp->pos < cSim.NLSize)
    {
        unsigned int tgx                                                = threadIdx.x & GRIDBITSMASK;       
        
        // Read NLRecord information
        if (tgx < 16)
            psWarp->nlRecord.array[tgx]                                 = cSim.pNLRecord[psWarp->pos].array[tgx];

        
        // Calculate Exclusion/neighbor list space required
        int atomOffset                                                  = psWarp->nlRecord.NL.neighborCells & NLATOMOFFSETMASK;
        if (tgx == 0)
        {
            uint2 homeCell                                              = cSim.pNLNonbondCellStartEnd[psWarp->nlRecord.NL.homeCell];   
            psWarp->homeCell.x                                          = homeCell.x;
            psWarp->homeCell.y                                          = homeCell.y;
            psWarp->nlEntry.NL.xyBufferOffset                           = (((psWarp->nlRecord.NL.neighborCells >> NLXCELLOFFSETSHIFT) & NLXCELLOFFSETMASK) << NLENTRYXBUFFEROFFSETSHIFT) | atomOffset;
            if ((psWarp->nlRecord.NL.neighborCell[0] & NLCELLTYPEMASK) == 0)
                psWarp->nlEntry.NL.xyBufferOffset                      |= NLENTRYHOMECELLMASK;
        }
        
        int ysize                                                       = max(0, (int)(psWarp->homeCell.y - psWarp->homeCell.x - atomOffset * cSim.NLAtomsPerWarp));
        if (ysize > 0)
            ysize                                                       = 1 + max(0, ysize - 1) / (cSim.NLAtomsPerWarp * cSim.NLYDivisor);

            
        // Calculate maximum required space
#if (PME_ATOMS_PER_WARP == 16)
        volatile BNLAtom* psAtom                                        = &sAtom[threadIdx.x & GRIDPADDINGMASK];
#endif
        unsigned int cells                                              = psWarp->nlRecord.NL.neighborCells >> NLCELLCOUNTSHIFT;
        if (tgx < cells)
        {
            uint2 cell                                                  = cSim.pNLNonbondCellStartEnd[psWarp->nlRecord.NL.neighborCell[tgx] >> NLCELLSHIFT]; 
            psWarp->cell[tgx].x                                         = cell.x;
            psWarp->cell[tgx].y                                         = cell.y;      
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
        
        if (tgx == 0)
        {
            uint totalXSize                                             = ((psWarp->atomList[0] + GRID - 1) >> GRIDBITS);
            uint offset                                                 = atomicAdd(cSim.pNLTotalOffset, totalXSize * ysize * cSim.NLOffsetPerWarp + cSim.NLAtomsPerWarp);
            cSim.pNLOffset[psWarp->pos]                                 = offset;  
            psWarp->offset                                              = offset;       
        }       
             
        // Generate actual neighbor list entry
        uint ypos                                                       = psWarp->homeCell.x + (psWarp->nlRecord.NL.neighborCells & NLATOMOFFSETMASK) * cSim.NLAtomsPerWarp;
        psWarp->nlEntry.NL.yAtom                                        = ypos;        
        psWarp->nlEntry.NL.yEnd                                         = (psWarp->homeCell.y << NLCELLSHIFT) | (psWarp->nlRecord.NL.homeCell & NLCELLTYPEMASK);
        psWarp->nlpos                                                   = 0;    

        while (ypos < psWarp->homeCell.y)
        {
        
            // Calculate y bounds
            uint ymax                                                   = min(ypos + cSim.NLAtomsPerWarp, psWarp->homeCell.y);
        
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
             
            // Cache exclusions if there's room and offset them relative to ypos so that the atoms of adjacent
            // cells increase monotonically (massive reduction in exclusion tests results)]
            uint totalExclusions                                        = 0;
            if (psWarp->exclusionMask[0] <= cSim.NLMaxExclusionsPerWarp)
            {
                psWarp->exclusionMask[tgx]                              = tempExclusions;
                int minExclusion                                        = cSim.atoms;
                int maxExclusion                                        = 0;
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
                if (tgx < 1)
                {
                    psWarp->minExclusion                                = min(psWarp->exclusionMask[tgx], psWarp->exclusionMask[tgx + 1]);
                    psWarp->maxExclusion                                = max(psWarp->atomList[tgx], psWarp->atomList[tgx + 1]);
                }
            }

#if 0 
            if (tgx == 1)
                printf("%06d %06d %06d %06d %06d %06d\n", psWarp->pos, psWarp->ypos, psWarp->nlpos, blockIdx.x * blockDim.x + threadIdx.x, psWarp->minExclusion, psWarp->maxExclusion);
#endif
            
#if 0
            if (tgx == 0)
                printf("A %06d %06d %06d %06d %06d\n", psWarp->pos, psWarp->ypos, psWarp->ymax, psWarp->homeCell.y, totalExclusions);
#endif  
         
            // Initialize Neighbor List variables for current line of entry
            psWarp->exclusionMask[tgx]                                  = 0;       
            unsigned int cpos                                           = 0;
            unsigned int atoms                                          = 0;
            psWarp->nlAtoms                                             = 0;
            uint minAtom                                                = cSim.atoms;
            uint maxAtom                                                = 0;
            unsigned int threadmask                                     = 1 << tgx; 
                                                       
            while (cpos < cells)
            {
                // Check for home cell
                psWarp->cellID                                          = psWarp->nlRecord.NL.neighborCell[cpos] & NLCELLTYPEMASK;
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
                    xpos                                                = psWarp->cell[cpos].x;                    
 
                // Read x atoms
                while (xpos < psWarp->cell[cpos].y)
                {
                    // Calculate number of atoms in this iteration
                    uint xmax                                           = min(xpos + GRID, psWarp->cell[cpos].y) - xpos; 

#if 0                      
                    if (tgx == 0)
                        printf("%6d %6d %6d %6d %6d\n", psWarp->pos, psWarp->ypos, psWarp->ymax, xpos, xmax);
#endif                                    

#if (PME_ATOMS_PER_WARP == 32)
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
                    unsigned int cellID                                 = psWarp->cellID;
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
                        int pred = 0;
                        if (pos < xmax)
                        {
#if (PME_ATOMS_PER_WARP == 32)
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
                                psWarp->atomList[atoms]                 = (atom << NLCELLSHIFT) | psWarp->cellID;
                                atoms++; 
                                

                                // Output GRID atoms if ready
                                if (atoms == GRID)
                                {
                                    cSim.pNLAtomList[psWarp->offset + tgx]  
                                                                        = psWarp->atomList[tgx];
                                    psWarp->nlAtoms                    += GRID;
                                    psWarp->offset                     += GRID; 
                                                                    
                                    // Search for y atom exclusions matching any x atom all at once (this should
                                    // reduce exclusion tests by a factor of approximately 100 overall 
                                    // But first rule out skipping exclusion test
                                  
                                    if ((minAtom <= psWarp->maxExclusion) && (maxAtom >= psWarp->minExclusion))    
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
                if ((tgx < atoms) && ((minAtom <= psWarp->maxExclusion) && (maxAtom >= psWarp->minExclusion)))    
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
                psWarp->nlAtoms                                        += atoms;
                psWarp->offset                                         += cSim.NLAtomsPerWarp;
                atoms                                                   = 0;
                psWarp->exclusionMask[tgx]                              = 0;
                psWarp->atomList[tgx]                                   = 0;           
            }
            psWarp->nlEntry.NL.xAtoms[psWarp->nlpos]                    = psWarp->nlAtoms;
            
#if 0                      
            unsigned int total = 0;
            for (int i = 0; i < psWarp->nlRecord.NL.neighborCells >> NLCELLCOUNTSHIFT; i++)
                total += ((psWarp->cell[tgx].y - psWarp->cell[tgx].x + GRID - 1) >> GRIDBITS) << GRIDBITS;
            if (tgx == 0)
                printf("%06d,%06d,%06d\n", psWarp->ypos, ((psWarp->nlAtoms + GRID - 1) >> GRIDBITS) << GRIDBITS, total);
#endif             

            // Advance to next swath of atoms
            ypos                                                       += cSim.NLYDivisor * cSim.NLAtomsPerWarp;
            psWarp->nlpos++;
        }

        // Output neighbor list entry
        cSim.pNLEntry[psWarp->pos].array[tgx]                           = psWarp->nlEntry.array[tgx];
        
        //if (tgx == 0)
        //    printf("%06d E %08d\n", psWarp->pos, psWarp->offset - cSim.pNLOffset[psWarp->pos]);

        // Advance to next NLRecord entry    
        if (tgx == 0)
            psWarp->pos                                                 = atomicAdd(cSim.pNLPosition, 1);
    }
}
      
