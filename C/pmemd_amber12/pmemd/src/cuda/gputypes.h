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

#include "cudaPROC.h"
#ifndef __GPUTYPES_H__
#define __GPUTYPES_H__
#include <stdio.h>
//#include <unistd.h>
#include <sys/types.h>
//#include <sys/ipc.h>
//#include <sys/shm.h>
#include <cuda.h>
#include <cufft.h>
#include <curand.h>
#include <vector_functions.h>
#include <cuda_runtime_api.h>
#include <builtin_types.h>
#include <cstring>
#ifdef MPI
#undef MPI
#include <mpi.h>
#define MPI
#endif

//#define CUDA_P2P

/* Enforce use of CUDA 4.2 or later due to rampant compiler bugs in CUDA 3.0, absence of curand library until 3.2,
   interprocess P2P support starting with CUDA 4.1, and Kepler support starting with CUDA 4.2 */
#if defined(CUDA_VERSION) && (CUDA_VERSION < 4020)
#error "CUDA support requires the use of a 4.2 or later CUDA toolkit. Aborting compilation."
#endif

/* Control of single and double precision use. If neither use_SPFP
   or use_DPDP are defined then the code runs in the default SPDP
   mode which makes mixed use of single and double precision as 
   needed. Defining use_DPDP makes it use double precision throughout
   while defining use_SPFP makes it use single precision for nonbonded
   force components, double-precision for bonded force components and SHAKE, 
   and 24.40 bit fixed-point for force accumulation.

   Note using Double Precision throughout will cause, in some cases large,
   performance degradation.

   This option is now set in the configure script and the original defines
   are left here just for reference.
*/
//#define use_DPDP
//#define use_SPDP
//#define use_SPFP

// Enforce definition of one and only one precision mode
#if !(defined(use_DPDP) && !defined(use_SPDP) && !defined(use_SPFP)) && \
    !(defined(use_SPDP) && !defined(use_DPDP) && !defined(use_SPFP)) && \
    !(defined(use_SPFP) && !defined(use_DPDP) && !defined(use_SPDP))
#error "You must define one and only one precision mode (use_SPFP, use_SPDP, or use_DPDP) to build pmemd.cuda. Aborting compilation."
#endif

#if defined(_MSC_VER)
#define __align(_boundary_size) __declspec(align(_boundary_size))
#else
#define __align(_boundary_size) __attribute__((aligned(_boundary_size)))
#endif

typedef double                 __align(8) aligned_double;
typedef unsigned long int      __align(8) aligned_uli;
typedef long long int          __align(8) aligned_lli;
typedef unsigned long long int __align(8) PMEUllInt;

#ifdef WINDOWS
typedef long long int PMEAccumulator;
typedef double PMEDouble;
#else
typedef long long int  __align(8)  PMEAccumulator;
typedef double  __align(8)  PMEDouble;
#endif
#ifdef MPI
static const MPI_Datatype MPI_PMEDOUBLE = MPI_DOUBLE_PRECISION;
#endif

#if defined(use_DPDP)
#ifdef WINDOWS
typedef double PMEForce;
#else
typedef double             __align(8)  PMEForce;
#endif
typedef double2            __align(16) PMEDouble2;
typedef double4            __align(32) PMEDouble4;
typedef double             __align(8)  PMEFloat;
typedef double2            __align(16) PMEFloat2;
typedef double4            __align(16) PMEFloat4;
typedef cufftDoubleComplex __align(16) PMEComplex;
#ifdef MPI
static const MPI_Datatype MPI_PMEFLOAT = MPI_DOUBLE_PRECISION;
static const MPI_Datatype MPI_PMEACCUMULATOR = MPI_LONG_LONG_INT;
#endif
#elif defined(use_SPFP)
#ifdef WINDOWS
typedef long long int PMEForce;
#else
typedef long long int      __align(8)  PMEForce;
#endif
typedef double2 __align(16) PMEDouble2;
typedef double4 __align(16) PMEDouble4;
typedef float PMEFloat;
typedef float2 __align(8)  PMEFloat2;
typedef float4 __align(16) PMEFloat4;
typedef cufftComplex PMEComplex;
#ifdef MPI
static const MPI_Datatype MPI_PMEFLOAT = MPI_FLOAT;
static const MPI_Datatype MPI_PMEACCUMULATOR = MPI_LONG_LONG_INT;
#endif
#else // use_SPDP
#ifdef WINDOWS
typedef double PMEForce;
#else
typedef double             __align(8)  PMEForce;
#endif
typedef double2 __align(16) PMEDouble2;
typedef double4 __align(16) PMEDouble4;
typedef float PMEFloat;
typedef float2 __align(8)  PMEFloat2;
typedef float4 __align(16) PMEFloat4;
typedef cufftComplex PMEComplex;
#ifdef MPI
static const MPI_Datatype MPI_PMEFLOAT = MPI_FLOAT;
static const MPI_Datatype MPI_PMEACCUMULATOR = MPI_DOUBLE_PRECISION;
#endif
#endif


enum SM_VERSION
{
    SM_10,
    SM_11,
    SM_12,
    SM_13,
    SM_2X,
    SM_3X,
};

// Miscellaneous constants
enum {
    GRID                        = 32,
    GRIDBITS                    = 5,
    GRIDBITSMASK                = (GRID - 1),
    GRIDPADDINGMASK             = 0xffffffff - (GRID - 1),
    CMAPRESOLUTION              = 24,
    CMAPDIMENSION               = CMAPRESOLUTION + 4,
    CMAPSTEPSIZE                = 15,
    NEIGHBORCELLS               = 14,
    CELLIDXBITS                 = 10,
    CELLIDYBITS                 = 10,
    CELLIDZBITS                 = 10,
    CELLIDYSHIFT                = CELLIDXBITS,
    CELLIDZSHIFT                = CELLIDXBITS + CELLIDYBITS,
    CELLIDMASK                  = 0x0000003f,
    CELLHASHXBITS               = 2,
    CELLHASHYBITS               = 2,
    CELLHASHZBITS               = 2,
    CELLHASHBITS                = (CELLHASHXBITS + CELLHASHYBITS + CELLHASHZBITS),
    CELLHASHX                   = (1 << CELLHASHXBITS),
    CELLHASHY                   = (1 << CELLHASHYBITS),
    CELLHASHZ                   = (1 << CELLHASHZBITS),
    CELLHASHXY                  = CELLHASHX * CELLHASHY,
    CELLHASHCELLS               = (1 << CELLHASHBITS),
    NLXENTRIES                  = 29,
    NLCELLSHIFT                 = 8,
    NLCELLTYPEMASK              = ((1 << NLCELLSHIFT) - 1),
    NLCELLCOUNTSHIFT            = 16,
    NLXCELLOFFSETSHIFT          = 8,
    NLXCELLOFFSETMASK           = ((1 << (NLCELLCOUNTSHIFT - NLXCELLOFFSETSHIFT)) - 1),   
    NLATOMOFFSETMASK            = ((1 << NLXCELLOFFSETSHIFT) - 1), 
    
    NLENTRYXBUFFEROFFSETSHIFT   = 8,
    NLENTRYYBUFFEROFFSETMASK    = ((1 << NLENTRYXBUFFEROFFSETSHIFT) - 1), 
    NLENTRYHOMECELLSHIFT        = 16,
    NLENTRYHOMECELLMASK         = (1 << NLENTRYHOMECELLSHIFT),
    NLENTRYXBUFFEROFFSETMASK    = ((1 << (NLENTRYHOMECELLSHIFT - NLENTRYXBUFFEROFFSETSHIFT)) - 1),
    NLEXCLUSIONSHIFT            = 8,
    NLEXCLUSIONATOMMASK         = ((1 << NLEXCLUSIONSHIFT) - 1),
    VIRIALOFFSET                = 17,
    AMDEDIHEDRALOFFSET          = 23,
    ENERGYTERMS                 = 24,
    PADDING                     = 16,
    PADDINGMASK                 = 0xfffffff0
};

struct NLRecordStruct
{
    unsigned int neighborCell[NEIGHBORCELLS];       // Bits 0:7 neighbor cell buffer ID, bits 8:31 neighbor cell numerical ID
    unsigned int homeCell;                          // Homecell numerical ID
    unsigned int neighborCells;                     // Bits 0:7 y atom offset, bits 8:15 x cell offset, bits 16:31 Neighbor cell count
};

union NLRecord
{
    NLRecordStruct NL;
    unsigned int array[NEIGHBORCELLS + 2];
};

struct NLEntryStruct
{
    unsigned int yAtom;                             // First y atom
    unsigned int yEnd;                              // Last Y atom
    unsigned int xyBufferOffset;                    // Bits 0:7 y buffer offset, bits 8:15 x buffer offset, bit 16 cell 0 flag
    unsigned int xAtoms[NLXENTRIES];                // Number of x atoms per row of y atoms
};

union NLEntry
{
    NLEntryStruct NL;
    unsigned int array[NLXENTRIES + 3];
};


// Kernel dimensions - we only support SM 1.3 or better due to double-precision requirements
static const int SM_13_THREADS_PER_BLOCK                        = 256;
static const int SM_2X_THREADS_PER_BLOCK                        = 512;
static const int SM_3X_THREADS_PER_BLOCK                        = 1024;
static const int SM_13_NLCALCULATE_OFFSETS_THREADS_PER_BLOCK    = 240;
static const int SM_2X_NLCALCULATE_OFFSETS_THREADS_PER_BLOCK    = 720;
static const int SM_3X_NLCALCULATE_OFFSETS_THREADS_PER_BLOCK    = 720;
static const int SM_13_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK = 192;
static const int SM_2X_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK = 672;
static const int SM_3X_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK = 640;
static const int SM_13_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK = 160;
static const int SM_2X_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK = 544;
static const int SM_3X_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK = 640;
static const int SM_3X_NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER   = 2;
static const int SM_2X_NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER   = 1;
static const int SM_13_NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER   = 1;
static const int SM_13_LOCALFORCES_THREADS_PER_BLOCK            = 192;
static const int SM_2X_LOCALFORCES_THREADS_PER_BLOCK            = 512;
static const int SM_3X_LOCALFORCES_THREADS_PER_BLOCK            = 1024;
static const int SM_13_CHARMMFORCES_THREADS_PER_BLOCK           = 256;
static const int SM_2X_CHARMMFORCES_THREADS_PER_BLOCK           = 512;
static const int SM_3X_CHARMMFORCES_THREADS_PER_BLOCK           = 1024;
static const int SM_13_NMRFORCES_THREADS_PER_BLOCK              = 256;
static const int SM_2X_NMRFORCES_THREADS_PER_BLOCK              = 512;
static const int SM_3X_NMRFORCES_THREADS_PER_BLOCK              = 512;
static const int SM_13_CLEARFORCES_THREADS_PER_BLOCK            = 512;
static const int SM_2X_CLEARFORCES_THREADS_PER_BLOCK            = 768;
static const int SM_3X_CLEARFORCES_THREADS_PER_BLOCK            = 1024;
static const int SM_13_NLCLEARFORCES_THREADS_PER_BLOCK          = 512;
static const int SM_2X_NLCLEARFORCES_THREADS_PER_BLOCK          = 768;
static const int SM_3X_NLCLEARFORCES_THREADS_PER_BLOCK          = 1024;
static const int SM_13_NLREDUCEFORCES_THREADS_PER_BLOCK         = 512;
static const int SM_2X_NLREDUCEFORCES_THREADS_PER_BLOCK         = 768;
static const int SM_3X_NLREDUCEFORCES_THREADS_PER_BLOCK         = 1024;
static const int SM_13_REDUCEFORCES_THREADS_PER_BLOCK           = 512;
static const int SM_2X_REDUCEFORCES_THREADS_PER_BLOCK           = 768;
static const int SM_3X_REDUCEFORCES_THREADS_PER_BLOCK           = 1024;
static const int SM_13_REDUCEBUFFER_THREADS_PER_BLOCK           = 384;
static const int SM_2X_REDUCEBUFFER_THREADS_PER_BLOCK           = 768;
static const int SM_3X_REDUCEBUFFER_THREADS_PER_BLOCK           = 1024;
static const int SM_13_SHAKE_THREADS_PER_BLOCK                  = 64;
static const int SM_2X_SHAKE_THREADS_PER_BLOCK                  = 64;
static const int SM_3X_SHAKE_THREADS_PER_BLOCK                  = 64;
static const int SM_13_SHAKE_BLOCKS                             = 3;
static const int SM_2X_SHAKE_BLOCKS                             = 8;
static const int SM_3X_SHAKE_BLOCKS                             = 16;
static const int SM_13_UPDATE_THREADS_PER_BLOCK                 = 256;
static const int SM_2X_UPDATE_THREADS_PER_BLOCK                 = 512;
static const int SM_3X_UPDATE_THREADS_PER_BLOCK                 = 1024;
#ifdef use_DPDP
static const int SM_13_GBBORNRADII_THREADS_PER_BLOCK            = 192;
static const int SM_13_GBBORNRADIIIGB78_THREADS_PER_BLOCK       = 128;
static const int SM_2X_GBBORNRADII_THREADS_PER_BLOCK            = 512;
static const int SM_3X_GBBORNRADII_THREADS_PER_BLOCK            = 1024;
static const int SM_3X_GBBORNRADII_BLOCKS_MULTIPLIER            = 1;
static const int SM_2X_GBBORNRADII_BLOCKS_MULTIPLIER            = 1;
static const int SM_13_GBBORNRADII_BLOCKS_MULTIPLIER            = 1;
static const int SM_13_GBNONBONDENERGY1_THREADS_PER_BLOCK       = 192;
static const int SM_2X_GBNONBONDENERGY1_THREADS_PER_BLOCK       = 512;
static const int SM_3X_GBNONBONDENERGY1_THREADS_PER_BLOCK       = 1024;
static const int SM_3X_GBNONBONDENERGY1_BLOCKS_MULTIPLIER       = 1;
static const int SM_2X_GBNONBONDENERGY1_BLOCKS_MULTIPLIER       = 1;
static const int SM_13_GBNONBONDENERGY1_BLOCKS_MULTIPLIER       = 1;
static const int SM_13_GBNONBONDENERGY2_THREADS_PER_BLOCK       = 192;
static const int SM_13_GBNONBONDENERGY2IGB78_THREADS_PER_BLOCK  = 128;
static const int SM_2X_GBNONBONDENERGY2_THREADS_PER_BLOCK       = 512;
static const int SM_3X_GBNONBONDENERGY2_THREADS_PER_BLOCK       = 512;
static const int SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK       = 160;
static const int SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK       = 512;
static const int SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK       = 1024;
static const int SM_3X_GBNONBONDENERGY2_BLOCKS_MULTIPLIER       = 2;
static const int SM_2X_GBNONBONDENERGY2_BLOCKS_MULTIPLIER       = 1;
static const int SM_13_GBNONBONDENERGY2_BLOCKS_MULTIPLIER       = 1;
static const int SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER       = 1;
static const int SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER       = 1;
static const int SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER       = 1;
static const int SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK       = 160;
static const int SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK       = 512;
static const int SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK       = 1024;
static const int SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER       = 1;
static const int SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER       = 1;
static const int SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER       = 1;
#else
static const int SM_13_GBBORNRADII_THREADS_PER_BLOCK            = 384;
static const int SM_13_GBBORNRADIIIGB78_THREADS_PER_BLOCK       = 320;
static const int SM_2X_GBBORNRADII_THREADS_PER_BLOCK            = 768;
static const int SM_3X_GBBORNRADII_THREADS_PER_BLOCK            = 768;
static const int SM_3X_GBBORNRADII_BLOCKS_MULTIPLIER            = 2;
static const int SM_2X_GBBORNRADII_BLOCKS_MULTIPLIER            = 1;
static const int SM_13_GBBORNRADII_BLOCKS_MULTIPLIER            = 1;
static const int SM_13_GBNONBONDENERGY1_THREADS_PER_BLOCK       = 320;
static const int SM_2X_GBNONBONDENERGY1_THREADS_PER_BLOCK       = 640;
static const int SM_3X_GBNONBONDENERGY1_THREADS_PER_BLOCK       = 576;
static const int SM_3X_GBNONBONDENERGY1_BLOCKS_MULTIPLIER       = 2;
static const int SM_2X_GBNONBONDENERGY1_BLOCKS_MULTIPLIER       = 1;
static const int SM_13_GBNONBONDENERGY1_BLOCKS_MULTIPLIER       = 1;
static const int SM_13_GBNONBONDENERGY2_THREADS_PER_BLOCK       = 384;
static const int SM_13_GBNONBONDENERGY2IGB78_THREADS_PER_BLOCK  = 320;
static const int SM_2X_GBNONBONDENERGY2_THREADS_PER_BLOCK       = 768;
static const int SM_3X_GBNONBONDENERGY2_THREADS_PER_BLOCK       = 640;
static const int SM_3X_GBNONBONDENERGY2_BLOCKS_MULTIPLIER       = 2;
static const int SM_2X_GBNONBONDENERGY2_BLOCKS_MULTIPLIER       = 1;
static const int SM_13_GBNONBONDENERGY2_BLOCKS_MULTIPLIER       = 1;
static const int SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK       = 256;
static const int SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK       = 768;
static const int SM_3X_PMENONBONDENERGY_THREADS_PER_BLOCK       = 640;
static const int SM_3X_PMENONBONDENERGY_BLOCKS_MULTIPLIER       = 2;
static const int SM_2X_PMENONBONDENERGY_BLOCKS_MULTIPLIER       = 1;
static const int SM_13_PMENONBONDENERGY_BLOCKS_MULTIPLIER       = 1;
static const int SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK       = 256;
static const int SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK       = 768;
static const int SM_3X_IPSNONBONDENERGY_THREADS_PER_BLOCK       = 640;
static const int SM_3X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER       = 2;
static const int SM_2X_IPSNONBONDENERGY_BLOCKS_MULTIPLIER       = 1;
static const int SM_13_IPSNONBONDENERGY_BLOCKS_MULTIPLIER       = 1;
#endif

static const int SM_13_MAXMOLECULES                             = 420;
static const int SM_2X_MAXMOLECULES                             = 1535;
static const int SM_3X_MAXMOLECULES                             = 1024;
static const int SM_13_MAXPSMOLECULES                           = 660;
static const int SM_2X_MAXPSMOLECULES                           = 2040;
static const int SM_3X_MAXPSMOLECULES                           = 2040;
static const int SM_13_READ_SIZE                                = 64;
static const int SM_2X_READ_SIZE                                = 128;
static const int SM_3X_READ_SIZE                                = 128;
static const int GB_TEXTURE_WIDTH                               = 512;
#ifdef DPDP
static const int MAX_RANDOM_STEPS                               = 64;
#else
static const int MAX_RANDOM_STEPS                               = 128;
#endif

// PME data
#define PI_VAL 3.1415926535897932384626433832795
static __constant__ const PMEDouble PI                                       = (PMEDouble) PI_VAL;
static __constant__ const PMEFloat PI_F                                      = (PMEFloat)  PI_VAL;

// PME data
#define KB_VAL (1.380658 * 6.0221367) / (4.184 * 1000.0)     //Boltzmann's constant in internal units 
static __constant__ const PMEDouble KB                                       = (PMEDouble) KB_VAL;
static __constant__ const PMEFloat KB_F                                      = (PMEFloat)  KB_VAL;

// Compilation flags
static const bool bShadowedOutputBuffers                        = false;    // Turns off sysmem shadowing of really large buffers

struct listdata_rec
{
    int         offset;
    int         cnt;
};

struct bond_rec
{
    int         atm_i;
    int         atm_j;
    int         parm_idx;
};

struct angle_rec
{
    int         atm_i;
    int         atm_j;
    int         atm_k;
    int         parm_idx;
};

struct dihed_rec
{
    int         atm_i;
    int         atm_j;
    int         atm_k;
    int         atm_l;
    int         parm_idx;
};

struct angle_ub_rec
{
    int         atm_i;
    int         atm_j;
    int         parm_idx;
};

struct dihed_imp_rec
{
    int         atm_i;
    int         atm_j;
    int         atm_k;
    int         atm_l;
    int         parm_idx;
};
struct cmap_rec
{
    int         atm_i;
    int         atm_j;
    int         atm_k;
    int         atm_l;
    int         atm_m;
    int         parm_idx;
};

struct shake_bond_rec
{
    int         atm_i;
    int         atm_j;
    double      parm;
};

struct gb_pot_ene_rec
{
    double total;
    double vdw_tot;
    double elec_tot;
    double gb;
    double surf;
    double bond;
    double angle;
    double dihedral;
    double vdw_14;
    double elec_14;
    double restraint;
    double angle_ub;
    double imp;
    double cmap;
};

struct pme_pot_ene_rec
{
    double total;
    double vdw_tot;      // total of dir, recip
    double vdw_dir;
    double vdw_recip;
    double elec_tot;     // total of dir, recip, nb_adjust, self
    double elec_dir;
    double elec_recip;
    double elec_nb_adjust;
    double elec_self;
    double hbond;
    double bond;
    double angle;
    double dihedral;
    double vdw_14;
    double elec_14;
    double restraint;
    double angle_ub;
    double imp;
    double cmap;
};

struct NTPData 
{
    PMEDouble last_recip[9];
    PMEDouble recip[9];
    PMEDouble ucell[9];
    PMEFloat recipf[9];
    PMEFloat ucellf[9];
    PMEFloat one_half_nonbond_skin_squared;
    PMEFloat cutPlusSkin2;
};

struct ep_frame_rec
{
    int extra_pnt[2];
    int ep_cnt;
    int type;
    int parent_atm;
    int frame_atm1;
    int frame_atm2;
    int frame_atm3;
};

#define ESCALE (1ll << 30)
#define FSCALE (1ll << 40)
static __constant__ const PMEDouble ENERGYSCALE                  = (PMEDouble)ESCALE;
static __constant__ const PMEFloat  ENERGYSCALEF                 = (PMEFloat)ESCALE;
static __constant__ const PMEDouble ONEOVERENERGYSCALE           = (PMEDouble)1.0 / (PMEDouble)ESCALE;
static __constant__ const PMEDouble ONEOVERENERGYSCALESQUARED    = (PMEDouble)1.0 / ((PMEDouble)ESCALE * (PMEDouble)ESCALE);
static __constant__ const PMEDouble FORCESCALE                   = (PMEDouble)FSCALE;
static __constant__ const PMEFloat  FORCESCALEF                  = (PMEFloat)FSCALE;
static __constant__ const PMEDouble ONEOVERFORCESCALE            = (PMEDouble)1.0 / (PMEDouble)FSCALE;
static __constant__ const PMEFloat ONEOVERFORCESCALEF            = (PMEDouble)1.0 / (PMEDouble)FSCALE;
static __constant__ const PMEDouble ONEOVERFORCESCALESQUARED     = (PMEDouble)1.0 / ((PMEDouble)FSCALE * (PMEDouble)FSCALE);

struct KineticEnergyRecord
{
    PMEFloat EKE;
    PMEFloat EKPH;
    PMEFloat EKPBS;
};

union KineticEnergy
{
    struct KineticEnergyRecord KE;
    PMEFloat array[3];
};


// Define this to move test issues with GPU RNG.  If whatever issue is being observed
// goes away, something's wrong there as this RNG is identical to the CPU RNG used
// by AMBER.
//#define CPU_RANDOMS
//#define GVERBOSE
//#define MEMTRACKING
//#define SYNCHRONOUS
#ifdef GVERBOSE
#ifndef MEMTRACKING
#define MEMTRACKING
#endif
#ifdef MPI
#define PRINTMETHOD(name) \
{ \
    printf("Method: %s on node %d\n", name, gpu->gpuID); \
    fflush(stdout); \
} 

#ifdef SYNCHRONOUS
#define LAUNCHERROR(s) \
    { \
        printf("Launched %s on node %d\n", s, gpu->gpuID); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
        cudaThreadSynchronize(); \
    }
#else
#define LAUNCHERROR(s) \
    { \
        printf("Launched %s on node %d\n", s, gpu->gpuID); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
    }
#endif
#define LAUNCHERROR_BLOCKING(s) \
    { \
        printf("Launched %s on node %d\n", s, gpu->gpuID); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
        cudaThreadSynchronize(); \
    }
#define LAUNCHERROR_NONBLOCKING(s) \
    { \
        printf("Launched %s on node %d\n", s, gpu->gpuID); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
    }
#else
#define PRINTMETHOD(name) \
{ \
    printf("Method: %s\n", name); \
    fflush(stdout); \
} 
#ifdef SYNCHRONOUS
#define LAUNCHERROR(s) \
    { \
        printf("Launched %s\n", s); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
        cudaThreadSynchronize(); \
    }
#else
#define LAUNCHERROR(s) \
    { \
        printf("Launched %s\n", s); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
    }
#endif
#define LAUNCHERROR_BLOCKING(s) \
    { \
        printf("Launched %s\n", s); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
        cudaThreadSynchronize(); \
    }
#define LAUNCHERROR_NONBLOCKING(s) \
    { \
        printf("Launched %s\n", s); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
    }
#endif
#else
#define PRINTMETHOD(name)
#ifdef SYNCHRONOUS
#define LAUNCHERROR(s) \
{ \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
        cudaThreadSynchronize(); \
    }
#else
#define LAUNCHERROR(s) \
{ \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
    }
#endif
#define LAUNCHERROR_BLOCKING(s) \
{ \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
        cudaThreadSynchronize(); \
    }
#define LAUNCHERROR_NONBLOCKING(s) \
{ \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
    }
#endif

#define RTERROR(status, s) \
    if (status != cudaSuccess) { \
        printf("%s %s\n", s, cudaGetErrorString(status)); \
        cudaThreadExit(); \
        exit(-1); \
    }

struct cudaSimulation {
    int                         grid;                               // Grid size
    int                         gridBits;                           // Grid bit count
    int                         atoms;                              // Total number of atoms
    int                         paddedNumberOfAtoms;                // Atom count padded out to fit grid size
    
    // AMBER parameters
    int                         ntp;                                // AMBER constant pressure setting
    int                         alpb;                               // Analytical Linearized Poisson Boltzmann setting
    int                         igb;                                // Generalized Born overall setting
    PMEDouble                   scnb;                               // 1-4 nonbond scale factor
    PMEDouble                   scee;                               // 1-4 electrostatic scale factor
    PMEFloat                    cut;                                // Nonbond interaction cutoff
    PMEFloat                    cut2;                               // Nonbond interaction cutoff squared
    PMEFloat                    cutPlusSkin;                        // Nonbond interaction cutooff plus skin
    PMEFloat                    cutPlusSkin2;                       // Nonbond interaction cutooff plus skin squared
    PMEDouble                   dielc;                              // Dielectric constant
    PMEDouble                   gamma_ln;                           // Langevin integration parameter
    PMEDouble                   c_ave;                              // Langevin integration parameter
    PMEDouble                   c_implic;                           // Langevin integration parameter
    PMEDouble                   c_explic;                           // Langevin integration parameter
    bool                        bUseVlimit;                         // Use vlimit flag
    PMEDouble                   vlimit;                             // vlimit 
    PMEDouble                   tol;                                // SHAKE tolerance
    PMEDouble                   massH;                              // Hydrogen mass
    aligned_double              invMassH;                           // Inverse hydrogen mass
    PMEFloat                    gb_alpha;                           // Generalized Born parameter
    PMEFloat                    gb_beta;                            // Generalized Born parameter
    PMEFloat                    gb_gamma;                           // Generalized Born parameter
    PMEFloat                    gb_kappa;                           // Generalized Born parameter
    PMEFloat                    gb_kappa_inv;                       // Generalized Born derived parameter    
    PMEFloat                    gb_cutoff;                          // Generalized Born cutoff
    PMEFloat                    gb_fs_max;                          // Generalized Born something or other
    PMEFloat                    gb_neckscale;                       // Generalized Born neck scaling factor
    PMEFloat                    gb_neckcut;                         // Generalized Born neck cutoff
    PMEFloat                    gb_neckoffset;                      // Generalized Born neck offset for LUT index
    PMEFloat                    rgbmax;                             // Generalized Born Born radius cutoff
    PMEFloat                    rgbmax1i;                           // Inverse Generalized Born Radius cutoff
    PMEFloat                    rgbmax2i;                           // Inverse Generalized Born Radius cutoff    
    PMEFloat                    rgbmaxpsmax2;                       // Generalized Born derived quantity
    PMEFloat                    intdiel;                            // Generalized Born interior dielectric
    PMEFloat                    extdiel;                            // Generalized Born exterior dielectric
    PMEFloat                    alpb_alpha;                         // Generalized Born parameter
    PMEFloat                    alpb_beta;                          // Generalized Born derived parameter
    PMEFloat                    extdiel_inv;                        // Generalized Born derived parameter    
    PMEFloat                    intdiel_inv;                        // Generalized Born derived parameter
    PMEFloat                    saltcon;                            // Generalized Born salt conductivity
    PMEFloat                    surften;                            // Generalized Born surface tension
    PMEFloat                    offset;                             // Generalized Born Born radius offset
    PMEFloat                    arad;                               // Generalized Bown atomic radius
    PMEFloat                    one_arad_beta;                      // Generalized Born derived parameters

    // PME parameters
    PMEDouble                   a;                                  // PBC x cell length
    PMEDouble                   b;                                  // PBC y cell length
    PMEDouble                   c;                                  // PBC z cell length
    PMEFloat                    af;                                 // PBC single-precision x cell length
    PMEFloat                    bf;                                 // PBC single-precision y cell length
    PMEFloat                    cf;                                 // PBC single-precision z cell length
    PMEFloat                    alpha;                              // PBC Alpha
    PMEFloat                    beta;                               // PBC Beta
    PMEFloat                    gamma;                              // PBC Gamma 
    PMEFloat                    pbc_box[3];                         // PBC Box
    PMEFloat                    reclng[3];                          // PBC reciprocal cell lengths
    PMEFloat                    cut_factor[3];                      // PBC cut factors
    PMEFloat                    ucell[3][3];                        // PBC cell coordinate system
    aligned_double              recip[3][3];                        // PBC reciprocal cell coordinate system
    PMEFloat                    recipf[3][3];                       // Single precision PBC reciprocal cell coordinate system
    PMEFloat                    uc_volume;                          // PBC cell volume
    PMEFloat                    uc_sphere;                          // PBC bounding sphere
    PMEFloat                    pi_vol_inv;                         // PBC/PME constant
    PMEFloat                    fac;                                // PME Ewald factor
    PMEFloat                    fac2;                               // PME Ewald factor x 2
    bool                        is_orthog;                          // PBC cell orthogonality flag
    int                         nfft1;                              // x fft grid size
    int                         nfft2;                              // y fft grid size
    int                         nfft3;                              // z fft grid size
    int                         nfft1xnfft2;                        // Product of nfft1 and nfft2
    int                         fft_x_dim;                          // x fft dimension
    int                         fft_y_dim;                          // y fft dimension
    int                         fft_z_dim;                          // z fft dimension
    int                         nf1;                                // x scalar sum coefficient
    int                         nf2;                                // y scalar sum coefficient
    int                         nf3;                                // z scalar sum coefficient
    int                         n2Offset;                           // Offset to y prefac data
    int                         n3Offset;                           // Offset to z prefac data
    int                         nSum;                               // Total prefac data
    int                         order;                              // Interpolation order
    int                         XYZStride;                          // Stride of PME buffers
    PMEFloat                    ew_coeff;                           // Ewald coefficient
    PMEFloat                    ew_coeff2;                          // Ewald coefficient squared
    PMEFloat                    negTwoEw_coeffRsqrtPI;              // Switching function constant
    
    // IPS paramers
    bool                        bIPSActive;                         // Flag to indicate IPS is active 
    PMEFloat                    rips;                               // Radius of IPS local region
    PMEFloat                    rips2;                              // rips^2
    PMEFloat                    ripsr;                              // 1/rips
    PMEFloat                    rips2r;                             // 1/rips^2
    PMEFloat                    rips6r;                             // 1/rips^6
    PMEFloat                    rips12r;                            // 1/rips^12
    PMEDouble                   eipssnb;                            // IPS self VDW energy
    PMEDouble                   eipssel;                            // IPS self-electrostatic energy
    PMEDouble                   virips;                             // IPS self-virial energy
    PMEDouble                   EIPSEL;                             // IPS exclusion self electrostatic energy
    PMEDouble                   EIPSNB;                             // IPS exclusion self Nonbond energy

    //  Electrostatic IPS parameters:
    PMEFloat                    aipse0;
    PMEFloat                    aipse1; 
    PMEFloat                    aipse2;
    PMEFloat                    aipse3;
    PMEFloat                    pipsec;
    PMEFloat                    pipse0;
    PMEFloat                    bipse1;
    PMEFloat                    bipse2;
    PMEFloat                    bipse3;
        
    //  Dispersion IPS parameters:
    PMEFloat                    aipsvc0;
    PMEFloat                    aipsvc1;
    PMEFloat                    aipsvc2;
    PMEFloat                    aipsvc3;
    PMEFloat                    pipsvcc;
    PMEFloat                    pipsvc0;
    PMEFloat                    bipsvc1;
    PMEFloat                    bipsvc2;  
    PMEFloat                    bipsvc3; 

    //  Repulsion IPS parameters:
    PMEFloat                    aipsva0;
    PMEFloat                    aipsva1;
    PMEFloat                    aipsva2;
    PMEFloat                    aipsva3;
    PMEFloat                    pipsvac;
    PMEFloat                    pipsva0;
    PMEFloat                    bipsva1;
    PMEFloat                    bipsva2;
    PMEFloat                    bipsva3;             
    
    // AMD parameters:
    int                         iamd;
    int                         iamdlag;
    int                         amd_print_interval;
    int                         AMDNumRecs;
    int                         AMDNumLag;
    PMEDouble                   AMDtboost;
    PMEDouble                   AMDfwgt;
    PMEUllInt*                  pAMDEDihedral;
    PMEDouble*                  pAMDfwgtd;
    PMEDouble                   amd_EthreshP;             
    PMEDouble                   amd_alphaP;             
    PMEDouble                   amd_EthreshD;             
    PMEDouble                   amd_alphaD;             
    PMEDouble                   amd_temp0;             

    // GBSA parameters:
    int                         igbsa;

    // NTP stuff
    NTPData*                    pNTPData;                           // PME NTP mutable values

#ifdef MPI 
    // MPI stuff
    int                         minLocalAtom;                       // First local atom
    int                         maxLocalAtom;                       // Last local atom
    int                         localAtoms;                         // Number of local atoms
    int                         localAtoms3;                        // Number of local atoms x 3 (used for clearing and reducing forces)
    int                         minProcessedAtom;                   // First atom touched locally 
    int                         maxProcessedAtom;                   // Last atom touched locally
    int                         processedAtoms;                     // Number of atoms touched locally
    int                         processedAtoms3;                    // Number of atoms touched locally x 3
    int                         nonLocalAtoms;                      // Number of remote atoms touched locally
    int                         nonLocalAtoms3;                     // Number of remote atoms touched locally x 3
    int                         minReducedAtom;                     // First atom force reduced locally 
    int                         maxReducedAtom;                     // Last atom force reduced locally
    int                         reducedAtoms;                       // Number of atoms force reduced locally
    int                         reducedAtoms3;                      // Number of atoms force reduced locally x 3
#endif
     
        
    // Atom stuff
    PMEDouble*                  pAtomX;                             // Atom X coordinates
    PMEDouble*                  pAtomY;                             // Atom Y coordinates
    PMEDouble*                  pAtomZ;                             // Atom Z coordinates
    PMEFloat2*                  pAtomXYSP;                          // Single Precision Atom X and Y coordinates
    PMEFloat*                   pAtomZSP;                           // Single Precision Atom Z coordinates
    int2*                       pAtomXYFP;                          // Fixed point Atom X and Y coordinates
    int*                        pAtomZFP;                           // Fixed point Atom Z coordinates
    PMEFloat2*                  pAtomSigEps;                        // Atom nonbond parameters
    PMEFloat*                   pAtomS;                             // Atom S Parameter
    PMEFloat*                   pAtomRBorn;                         // Atom Born Radius
    PMEDouble*                  pAtomCharge;                        // Atom charges
    PMEFloat*                   pAtomChargeSP;                      // Single precision atom charges
    PMEDouble*                  pAtomMass;                          // Atom masses
    PMEDouble*                  pAtomInvMass;                       // Atom inverse masses
    PMEDouble*                  pReff;                              // Effective Born Radius
    PMEFloat*                   pPsi;                               // GB intermediate psi value
    PMEFloat*                   pReffSP;                            // Single Precision Effective Born Radius
    PMEFloat*                   pTemp7;                             // Single Precision Born Force
    PMEDouble*                  pForce;                             // Atom forces
    PMEDouble*                  pForceX;                            // Atom X forces
    PMEDouble*                  pForceY;                            // Atom Y forces
    PMEDouble*                  pForceZ;                            // Atom Z forces
    PMEDouble*                  pNBForce;                           // Atom nonbond forces    
    PMEDouble*                  pNBForceX;                          // Atom X nonbond forces
    PMEDouble*                  pNBForceY;                          // Atom Y nonbond forces
    PMEDouble*                  pNBForceZ;                          // Atom Z nonbond forces       
#ifdef MPI
    PMEDouble*                  pTemp7a;                            // Single Precision Born Force
    PMEDouble*                  pReffa;                             // Effective Born Radius
    PMEFloat*                   pPMEForce;                          // PME Force
    PMEDouble*                  pOutForce;                          // Outgoing MPI Forces
    PMEDouble*                  pInForce;                           // Incoming MPI Forces
#endif    
    PMEDouble*                  pVelX;                              // Atom X velocities
    PMEDouble*                  pVelY;                              // Atom Y velocities
    PMEDouble*                  pVelZ;                              // Atom Z velocities
    PMEDouble*                  pLVelX;                             // Atom X last velocities
    PMEDouble*                  pLVelY;                             // Atom Y last velocities
    PMEDouble*                  pLVelZ;                             // Atom Z last velocities
    unsigned int*               pOutputBufferCounter;               // Output buffer counter for bonded interactions
    int                         maxNonbonds;                        // Maximum nonbond interaction buffers
    PMEFloat*                   pXMax;                              // Recentering maximum x
    PMEFloat*                   pYMax;                              // Recentering maximum y
    PMEFloat*                   pZMax;                              // Recentering maximum z
    PMEFloat*                   pXMin;                              // Recentering minimum x
    PMEFloat*                   pYMin;                              // Recentering minimum y
    PMEFloat*                   pZMin;                              // Recentering minimum z
    
    // Neighbor List stuff
    PMEFloat                    skinnb;                             // Input Nonbond skin
    unsigned int*               pImageIndex;                        // Image # for spatial sort
    unsigned int*               pImageIndex2;                       // Image # for spatial sort
    unsigned int*               pImageHash;                         // Image hash for spatial sort
    unsigned int*               pImageHash2;                        // Image hash for spatial sort
    unsigned int*               pImageAtom;                         // Image atom #
    unsigned int*               pImageAtom2;                        // Image atom #
    unsigned int*               pImageAtomLookup;                   // Original atom lookup table
    PMEFloat2*                  pAtomXYSaveSP;                      // Saved atom coordinates from neighbor list generation
    PMEFloat*                   pAtomZSaveSP;                       // Saved atom coordinates from neighbor list generation
    PMEDouble*                  pImageX;                            // Image x coordinates
    PMEDouble*                  pImageY;                            // Image y coordinates
    PMEDouble*                  pImageZ;                            // Image z coordinates
    PMEDouble*                  pImageVelX;                         // Image x velocities
    PMEDouble*                  pImageVelY;                         // Image y velocities
    PMEDouble*                  pImageVelZ;                         // Image z velocities
    PMEDouble*                  pImageLVelX;                        // Image last x velocities
    PMEDouble*                  pImageLVelY;                        // Image last y velocities
    PMEDouble*                  pImageLVelZ;                        // Image last z velocities
    PMEDouble*                  pImageMass;                         // Image masses
    PMEDouble*                  pImageInvMass;                      // Image inverse masses
    PMEDouble*                  pImageCharge;                       // Image charges
    PMEFloat2*                  pImageSigEps;                       // Image sigma/epsilon data for nonbond interactions
    unsigned int*               pImageOutputBuffers;                // Image per atom output buffer count
    unsigned int*               pImageCellID;                       // Image cell ID for calculating local coordinate system
    PMEDouble*                  pImageX2;                           // Image x coordinates
    PMEDouble*                  pImageY2;                           // Image y coordinates
    PMEDouble*                  pImageZ2;                           // Image z coordinates
    PMEDouble*                  pImageVelX2;                        // Image x velocities
    PMEDouble*                  pImageVelY2;                        // Image y velocities
    PMEDouble*                  pImageVelZ2;                        // Image z velocities
    PMEDouble*                  pImageLVelX2;                       // Image last x velocities
    PMEDouble*                  pImageLVelY2;                       // Image last y velocities
    PMEDouble*                  pImageLVelZ2;                       // Image last z velocities
    PMEDouble*                  pImageMass2;                        // Image masses
    PMEDouble*                  pImageInvMass2;                     // Image inverse masses
    PMEDouble*                  pImageCharge2;                      // Image charges
    PMEFloat2*                  pImageSigEps2;                      // Image sigma/epsilon data for nonbond interactions
    unsigned int*               pImageOutputBuffers2;               // Image per atom output buffer count
    unsigned int*               pImageCellID2;                      // Image cell ID for calculating local coordinate system
    int4*                       pImageBondID;                       // Remapped Bond i, j, ibuff, jbuff
    int4*                       pImageBondAngleID1;                 // Remapped Bond Angle i, j, k, ibuff
    int2*                       pImageBondAngleID2;                 // Remapped Bond Angle jbuff, kbuff;    
    int4*                       pImageDihedralID1;                  // Remapped Dihedral i, j, k, l
    int4*                       pImageDihedralID2;                  // Remapped Dihedral ibuff, jbuff, kbuff, lbuff
    int4*                       pImageNb14ID;                       // Remapped 1-4 nonbond i, j, ibuff, jbuff                
    int2*                       pImageConstraintID;                 // Remapped Atom constraint ID 
    int4*                       pImageUBAngleID;                    // Remapped Urey Bradley Angle i, j, ibuff, jbuff
    int4*                       pImageImpDihedralID1;               // Remapped Improper Dihedral i, j, k, l 
    int4*                       pImageImpDihedralID2;               // Remapped Improper Dihedral ibuff, jbuff, kbuff, lbuff
    int4*                       pImageCmapID1;                      // Remapped Cmap i, j, k, l
    int4*                       pImageCmapID2;                      // Remapped Cmap m, ibuff, jbuff, kbuff
    int2*                       pImageCmapID3;                      // Remapped Cmap lbuff, mbuff
    int4*                       pImageNMRDistanceID;                // Remapped NMR Distance i, j, ibuff, kbuff
    int4*                       pImageNMRAngleID1;                  // Remapped NMR Angle i, j, k, ibuff
    int2*                       pImageNMRAngleID2;                  // Remapped NMR Angle jbuff, kbuff
    int4*                       pImageNMRTorsionID1;                // Remapped NMR Torsion i, j, k, l
    int4*                       pImageNMRTorsionID2;                // Remapped NMR Torsion ibuff, jbuff, kbuff, lbuff
    int4*                       pImageShakeID;                      // Remapped Shake Atom ID    
    int4*                       pImageFastShakeID;                  // Remapped Fast Shake Atom ID     
    int*                        pImageSlowShakeID1;                 // Remapped Slow Shake central Atom ID     
    int4*                       pImageSlowShakeID2;                 // Remapped Slow Shake hydrogen Atom IDs
    int4*                       pImageSolventAtomID;                // Remapped solvent molecules/ions
    int*                        pImageSoluteAtomID;                 // Remapped solute atoms
    uint2*                      pNLNonbondCellStartEnd;             // Nonbond cell boundaries pointer
#ifndef use_SPFP
    int*                        pNLChargeGridBufferOffset;          // Precalculated per cell charge grid buffer
    int*                        pNLOddBufferOverlapFlag;            // Precalculated per axis per cell odd buffer overlap indicator
    int                         extraChargeGridBuffers[4];          // extra charge grid buffers per odd axis overlap count
#endif
    int                         maxChargeGridBuffers;               // Total possible charge grid buffers
    int                         cells;                              // Total number of nonbond cells
    int                         xcells;                             // Number of x cells
    int                         ycells;                             // Number of y cells
    int                         zcells;                             // Number of z cells
    int                         xcellsminusone;                     // Number of x cells minus one
    int                         ycellsminusone;                     // Number of y cells minus one
    int                         zcellsminusone;                     // Number of z cells minus one
    int                         xycells;                            // x cells * y cells
    PMEDouble                   xcell;                              // x cell dimension
    PMEDouble                   ycell;                              // y cell dimension
    PMEDouble                   zcell;                              // z cell dimension
    PMEDouble                   minCellX;                           // Minimum x cell coordinate
    PMEDouble                   minCellY;                           // Minimum y cell coordinate
    PMEDouble                   minCellZ;                           // Minimum z cell coordinate
    PMEDouble                   maxCellX;                           // Maximum x cell coordinate
    PMEDouble                   maxCellY;                           // Maximum y cell coordinate
    PMEDouble                   maxCellZ;                           // Maximum z cell coordinate
    bool                        bOddXCells;                         // Odd x cell count
    bool                        bOddYCells;                         // Odd y cell count
    bool                        bOddZCells;                         // Odd z cell count        
    aligned_double              oneOverXcells;                      // Fractional x cell dimension
    aligned_double              oneOverYcells;                      // Fractional y cell dimension
    aligned_double              oneOverZcells;                      // Fractional z cell dimension
    PMEFloat                    oneOverXcellsf;                     // Single precision fractional x cell dimension
    PMEFloat                    oneOverYcellsf;                     // Single precision fractional y cell dimension
    PMEFloat                    oneOverZcellsf;                     // Single precision fractional z cell dimension 
    PMEFloat                    cell;                               // minimum cell dimension
    PMEFloat                    nonbond_skin;                       // Effective nonbond skin
    PMEFloat                    one_half_nonbond_skin_squared;      // Skin test atom movement threshold
    bool*                       pNLbSkinTestFail;                   // Skin test result buffer
    unsigned int*               pNLCellHash;                        // Spatial ordering hash for within cells
    int                         NLWorkUnits;                        // Number of cell to cell interactions
    unsigned int*               pNLTotalOffset;                     // Pointer to total offset    
    int                         NLMaxTotalOffset;                   // Maximum available exclusion masks
    unsigned int*               pNLAtomList;                        // Pointer to neighbor list atoms
    unsigned int*               pNLExclusionList;                   // Pointer to list of nonbond exclusions
    uint2*                      pNLExclusionStartCount;             // Pointer to per-atom exclusions
    unsigned int                NLExclusions;                       // Total number of exclusions   
    unsigned int                NLAtoms;                            // Total number of neighbor list atoms   
    NLRecord*                   pNLRecord;                          // Pointer to neighbor list records
    NLEntry*                    pNLEntry;                           // Active neighbor list
    unsigned int*               pNLOffset;                          // Pointer to atom/exclusion offsets
    unsigned int*               pNLPosition;                        // Position in building neighbor list
    unsigned int                NLBuildWarps;                       // Number of warps in neighbor list build
    unsigned int                NLNonbondWarps;                     // Number of warps in nonbond energy/force calculation
    unsigned int                NLCellBuffers;                      // Number of nonbond cell buffers    
    unsigned int                NLSize;                             // Number of neighbor list records/entries
    unsigned int                NLXDivisor;                         // Neighbor List horizontal atom divisions
    unsigned int                NLYDivisor;                         // Neighbor List vertical atom divisions
    unsigned int                NLYStride;                          // Neighbor List vertical atom stride
    unsigned int                NLEntryTypes;                       // Neighbor List entry types (used to determine output buffers offset)
    unsigned int                NLHomeCellBuffer;                   // Register atom buffer
    unsigned int                NLAtomsPerWarp;                     // Number of atoms to process in each warp's registers
    unsigned int                NLAtomsPerWarpBits;                 // Number of bits in atoms to toprocess in each warp's registers
    unsigned int                NLAtomsPerWarpBitsMask;             // NLAtomsPerWarp - 1
    unsigned int                NLAtomsPerWarpMask;                 // First AtomsPerWarp worth of bits set to 1
    unsigned int                NLOffsetPerWarp;                    // Number of entries needed to process 1 warp iteration's worth of nonbond forces
    unsigned int                NLMaxExclusionsPerWarp;             // Number of atoms to process in each warp's registers
    unsigned int                maxNonbondBuffers;                  // maximum nonbond buffers for NTP simulation
    unsigned int*               pBNLExclusionBuffer;                // Per-warp GMEM exclusion buffer
    
    // GB Data
    PMEFloat2*                  pNeckMaxValPos;                     // GB Born Radii and energy correction data
    PMEFloat*                   pgb_alpha;                          // Pointer to per-atom GB alpha
    PMEFloat*                   pgb_beta;                           // Pointer to per-atom GB beta 
    PMEFloat*                   pgb_gamma;                          // Pointer to per-atom GB gamma
    
    
    // PME data
#ifdef use_SPFP   
    long long int*              plliXYZ_q;                          // Input PME charge grid
#endif
    PMEFloat*                   pXYZ_q;                             // PME charge grid/buffer
    PMEComplex*                 pXYZ_qt;                            // FFTed PME charge grid
    PMEFloat*                   pPrefac1;                           // PME nfft1 pre-factors
    PMEFloat*                   pPrefac2;                           // PME nfft2 pre-factors
    PMEFloat*                   pPrefac3;                           // PME nfft3 pre-factors
    int*                        pIFractX;                           // PME unit cell x coordinates
    int*                        pIFractY;                           // PME unit cell y coordinates
    int*                        pIFractZ;                           // PME unit cell z coordinates   
    PMEFloat4*                  pThetaX;                            // PME x theta spline weights
    PMEFloat4*                  pThetaY;                            // PME y theta spline weights
    PMEFloat4*                  pThetaZ;                            // PME z theta spline weights
    PMEFloat4*                  pDThetaX;                           // PME x theta spline weights
    PMEFloat4*                  pDThetaY;                           // PME y theta spline weights
    PMEFloat4*                  pDThetaZ;                           // PME z theta spline weights
    
    // NTP molecule data
    int                         soluteMolecules;                    // Total solute molecules
    int                         soluteMoleculeStride;               // Total solute molecule stride
    int                         soluteAtoms;                        // Total solute atoms
    int                         soluteAtomsOffset;                  // Used for remapping neighbor list data
    int                         solventMolecules;                   // Total solvent molecules
    int                         solventMoleculeStride;              // Total solvent molecules padded to warp width
    int*                        pSoluteAtomMoleculeID;              // List of solute molecule IDs
    int*                        pSoluteAtomID;                      // List of solute atom IDs
    PMEDouble*                  pSoluteAtomMass;                    // Solute atom masses
    PMEDouble*                  pSoluteCOMX;                        // X Last centers of mass for each solute molecule
    PMEDouble*                  pSoluteCOMY;                        // Y Last centers of mass for each solute molecule
    PMEDouble*                  pSoluteCOMZ;                        // Z Last centers of mass for each solute molecule
    PMEDouble*                  pSoluteDeltaCOMX;                   // X change in center of mass for each solute molecule
    PMEDouble*                  pSoluteDeltaCOMY;                   // Y change in center of mass for each solute molecule
    PMEDouble*                  pSoluteDeltaCOMZ;                   // Z change in center of mass for each solute molecule    
    PMEUllInt*                  pSoluteUllCOMX;                     // X Current center of mass for each solute molecule
    PMEUllInt*                  pSoluteUllCOMY;                     // Y Current center of mass for each solute molecule
    PMEUllInt*                  pSoluteUllCOMZ;                     // Z Current center of mass for each solute molecule  
    PMEUllInt*                  pSoluteUllEKCOMX;                   // Pointer to x component of COM Kinetic energy buffer
    PMEUllInt*                  pSoluteUllEKCOMY;                   // Pointer to x component of COM Kinetic energy buffer
    PMEUllInt*                  pSoluteUllEKCOMZ;                   // Pointer to x component of COM Kinetic energy buffer       
    PMEDouble*                  pSoluteInvMass;                     // Total Inverse mass for each solute molecule
    int4*                       pSolventAtomID;                     // List of solvent molecules/ions of 4 or fewer atoms
    PMEDouble*                  pSolventAtomMass1;                  // First solvent atom mass
    PMEDouble*                  pSolventAtomMass2;                  // Second solvent atom mass
    PMEDouble*                  pSolventAtomMass3;                  // Third solvent atom mass
    PMEDouble*                  pSolventAtomMass4;                  // Fourth solvent atom mass
    PMEDouble*                  pSolventCOMX;                       // X Last centers of mass for each solvent molecule
    PMEDouble*                  pSolventCOMY;                       // Y Last centers of mass for each solvent molecule    
    PMEDouble*                  pSolventCOMZ;                       // Z Last centers of mass for each solvent molecule
    PMEDouble*                  pSolventInvMass;                    // Total inverse mass for eache solvent molecule
    
    // NTP molecule constraint data
    PMEDouble*                  pConstraintAtomX;                   // Pre-centered constraint atom x
    PMEDouble*                  pConstraintAtomY;                   // Pre-centered constraint atom y
    PMEDouble*                  pConstraintAtomZ;                   // Pre-centered constraint atom z
    PMEDouble*                  pConstraintCOMX;                    // Original center of mass X for constraint atom in fractional coordinates
    PMEDouble*                  pConstraintCOMY;                    // Original center of mass Y for constraint atom in fractional coordinates
    PMEDouble*                  pConstraintCOMZ;                    // Original center of mass Z for constraint atom in fractional coordinates 




    // Energy and Virial Buffers
    unsigned long long int*     pEnergyBuffer;                      // Generic energy buffer pointer
    unsigned long long int*     pEELT;                              // Pointer to electrostatic energy
    unsigned long long int*     pEVDW;                              // Pointer to vdw energy
    unsigned long long int*     pEGB;                               // Pointer to Generalized Born energy
    unsigned long long int*     pEBond;                             // Pointer to bond energy
    unsigned long long int*     pEAngle;                            // Pointer to bond angle energy 
    unsigned long long int*     pEDihedral;                         // Pointer to dihedral energy 
    unsigned long long int*     pEEL14;                             // Pointer to 1-4 electrostatic energy 
    unsigned long long int*     pENB14;                             // Pointer to 1-4 vdw energy
    unsigned long long int*     pERestraint;                        // Pointer to restraint energy 
    unsigned long long int*     pEER;                               // Pointer to PME reciprocal electrostatic energy
    unsigned long long int*     pEED;                               // Pointer to PME direct electrostatic energy
    unsigned long long int*     pEAngle_UB;                         // Pointer to CHARMM Urey Bradley energy
    unsigned long long int*     pEImp;                              // Pointer to CHARMM improper dihedral energy
    unsigned long long int*     pECmap;                             // Pointer to CHARMM cmap energy
    unsigned long long int*     pENMRDistance;                      // Pointer to NMR distance energy
    unsigned long long int*     pENMRAngle;                         // Pointer to NMR angle energy
    unsigned long long int*     pENMRTorsion;                       // Pointer to NMR torsion energy
    unsigned long long int*     pVirial;                            // Pointer to PME virial
    unsigned long long int*     pVirial_11;                         // Pointer to PME virial component
    unsigned long long int*     pVirial_22;                         // Pointer to PME virial component
    unsigned long long int*     pVirial_33;                         // Pointer to PME virial component
    unsigned long long int*     pEKCOMX;                            // Pointer to x component of PME center of mass kinetic energy
    unsigned long long int*     pEKCOMY;                            // Pointer to y component of PME center of mass kinetic energy
    unsigned long long int*     pEKCOMZ;                            // Pointer to z component of PME center of mass kinetic energy
    
    // Kinetic Energy Buffers
    KineticEnergy*              pKineticEnergy;                     // Pointer to per-block Kinetic Energy entries
    
    // Random Number stuff
    unsigned int                randomSteps;                        // Number of steps between RNG calls
    unsigned int                randomNumbers;                      // Number of randoms to generate per atom per RNG call
    unsigned int*               pRandomPos;                         // RNG position per block
    PMEDouble*                  pRandom;                            // Pointer to overall RNG buffer
    PMEDouble*                  pRandomX;                           // Pointer to x random numbers
    PMEDouble*                  pRandomY;                           // Pointer to y random numbers
    PMEDouble*                  pRandomZ;                           // Pointer to z random numbers
    
    // Extra points stuff
    int                         EPs;
    int                         EP11s;
    int                         EP12s;
    int                         EP21s;
    int                         EP22s;
    int                         EP11Offset;
    int                         EP12Offset;
    int                         EP21Offset;
    int                         EP22Offset;
    int4*                       pExtraPoint11Frame;
    int*                        pExtraPoint11Index;
    PMEDouble*                  pExtraPoint11X;
    PMEDouble*                  pExtraPoint11Y;
    PMEDouble*                  pExtraPoint11Z;
    int4*                       pExtraPoint12Frame;
    int*                        pExtraPoint12Index;
    PMEDouble*                  pExtraPoint12X;
    PMEDouble*                  pExtraPoint12Y;
    PMEDouble*                  pExtraPoint12Z;
    int4*                       pExtraPoint21Frame;
    int2*                       pExtraPoint21Index;
    PMEDouble*                  pExtraPoint21X1;
    PMEDouble*                  pExtraPoint21Y1;
    PMEDouble*                  pExtraPoint21Z1;    
    PMEDouble*                  pExtraPoint21X2;
    PMEDouble*                  pExtraPoint21Y2;
    PMEDouble*                  pExtraPoint21Z2;        
    int4*                       pExtraPoint22Frame;
    int2*                       pExtraPoint22Index;
    PMEDouble*                  pExtraPoint22X1;
    PMEDouble*                  pExtraPoint22Y1;
    PMEDouble*                  pExtraPoint22Z1;    
    PMEDouble*                  pExtraPoint22X2;
    PMEDouble*                  pExtraPoint22Y2;
    PMEDouble*                  pExtraPoint22Z2;  
    int4*                       pImageExtraPoint11Frame;
    int*                        pImageExtraPoint11Index;
    int4*                       pImageExtraPoint12Frame;
    int*                        pImageExtraPoint12Index;
    int4*                       pImageExtraPoint21Frame;
    int2*                       pImageExtraPoint21Index;    
    int4*                       pImageExtraPoint22Frame;
    int2*                       pImageExtraPoint22Index;         
    
    // Shake constraint stuff
    unsigned int                shakeConstraints;                   // Traditional SHAKE constraints
    unsigned int                shakeOffset;                        // Offset to end of traditional SHAKE constraints
    unsigned int                fastShakeConstraints;               // Fast SHAKE constraints (H2O molecules)
    unsigned int                fastShakeOffset;                    // Offset to end of fast SHAKE constraints
    unsigned int                slowShakeConstraints;               // XH4 (really slow) SHAKE constraints
    unsigned int                slowShakeOffset;                    // Offset to end of slow SHAKE constraints
    int4*                       pShakeID;                           // SHAKE central atom plus up to 3 hydrogens
    double2*                    pShakeParm;                         // SHAKE central atom mass and equilibrium bond length
    int4*                       pFastShakeID;                       // H2O oxygen plus two hydrogens atom ID
    int*                        pSlowShakeID1;                      // Central atom of XH4 Shake constraint 
    int4*                       pSlowShakeID2;                      // XH4 SHAKE constraint hydrogens
    double2*                    pSlowShakeParm;                     // XH4 SHAKE central atom mass and equilibrium bond length
    aligned_double              ra;                                 // Fast SHAKE parameter
    aligned_double              ra_inv;                             // Fast SHAKE parameter
    aligned_double              rb;                                 // Fast SHAKE parameter
    aligned_double              rc;                                 // Fast SHAKE parameter 
    aligned_double              rc2;                                // Fast SHAKE parameter 
    aligned_double              hhhh;                               // Fast SHAKE parameter
    aligned_double              wo_div_wohh;                        // Fast SHAKE parameter
    aligned_double              wh_div_wohh;                        // Fast SHAKE parameter    
    
    // Bonded interaction stuff
    int                         bonds;                              // Total number of bonds
    int                         bondOffset;                         // Offset to end of bondss
    int                         bondAngles;                         // Total number of bond angles
    int                         bondAngleOffset;                    // Offset to end of bond angles
    int                         dihedrals;                          // Total number of dihedrals
    int                         dihedralOffset;                     // Offset to end of dihedrals
    int                         nb14s;                              // Total number of 1-4 nonbond interactions
    int                         nb14Offset;                         // Offset to end of 1-4 nonbond interactions
    int                         constraints;                        // Number of positional constraints
    int                         constraintOffset;                   // Offset to end of positional constraints
    int                         UBAngles;                           // Total number of Urey Bradley angles
    int                         UBAngleOffset;                      // Offset to end of Urey Bradley Angles
    int                         impDihedrals;                       // Total number of improper dihedrals
    int                         impDihedralOffset;                  // Offset to end of improper dihedrals
    int                         cmaps;                              // Total number of cmap terms
    int                         cmapOffset;                         // Offset to end of cmap terms
    PMEDouble2*                 pBond;                              // Bond rk, req
    int4*                       pBondID;                            // Bond i, j, ibuff, jbuff
    PMEDouble2*                 pBondAngle;                         // Bond Angle Kt, teq
    int4*                       pBondAngleID1;                      // Bond Angle i, j, k, ibuff
    int2*                       pBondAngleID2;                      // Bond Angle jbuff, kbuff;
    PMEDouble2*                 pDihedral1;                         // Dihedral Ipn, pn
    PMEDouble2*                 pDihedral2;                         // Dihedral pk, gamc
    PMEDouble*                  pDihedral3;                         // Dihedral gams
    int4*                       pDihedralID1;                       // Dihedral i, j, k, l
    int4*                       pDihedralID2;                       // Dihedral ibuff, jbuff, kbuff, lbuff
    PMEDouble2*                 pNb141;                             // 1-4 nonbond scee, scnb0
    PMEDouble2*                 pNb142;                             // 1-4 nonbond cn1, cn2
    int4*                       pNb14ID;                            // 1-4 nonbond i, j, ibuff, jbuff
    PMEDouble2*                 pConstraint1;                       // Constraint weight and xc
    PMEDouble2*                 pConstraint2;                       // Constraint yc and zc                      
    int2*                       pConstraintID;                      // Atom constraint ID
    
    // CHARMM interaction stuff
    PMEDouble2*                 pUBAngle;                           // Urey Bradley Angle rk, r0
    int4*                       pUBAngleID;                         // Urey Bradley Angle i, j, ibuff, jbuff
    PMEDouble2*                 pImpDihedral;                       // Improper Dihedral pk, phase
    int4*                       pImpDihedralID1;                    // Improper Dihedral i, j, k, l 
    int4*                       pImpDihedralID2;                    // Improper Dihedral ibuff, jbuff, kbuff, lbuff
    int4*                       pCmapID1;                           // Cmap i, j, k, l
    int4*                       pCmapID2;                           // Cmap m, ibuff, jbuff, kbuff
    int2*                       pCmapID3;                           // Cmap lbuff, mbuff
    int*                        pCmapType;                          // Cmap type

    // NMR refinement stuff
    int                         NMRDistances;                       // Number of NMR distance constraints
    int                         NMRDistanceOffset;                  // Offset to end of distance constraints
    int                         NMRAngles;                          // Number of NMR angle constraints
    int                         NMRAngleOffset;                     // Offset to end of angle constraints
    int                         NMRTorsions;                        // Number of NMR torsion constraints
    int                         NMRTorsionOffset;                   // Offset to end of torsion constraints
    bool                        bJar;                               // Determines whether Jarzynski MD is active
    double                      drjar;                              // Jar increment
    double*                     pNMRJarData;                        // Jarzynski accumulated work data
    int4*                       pNMRDistanceID;                     // NMR distance i, j, ibuff, jbuff
    PMEDouble2*                 pNMRDistanceR1R2;                   // NMR distance computed r1, r2
    PMEDouble2*                 pNMRDistanceR3R4;                   // NMR distance computed r3, r4
    PMEDouble2*                 pNMRDistanceK2K3;                   // NMR distance computed k2, k3
    PMEDouble*                  pNMRDistanceAve;                    // NMR distance restraint linear and exponential averaged value
    PMEDouble2*                 pNMRDistanceTgtVal;                 // NMR distance target and value for current step
    int2*                       pNMRDistanceStep;                   // NMR distance first and last step for application of restraint
    int*                        pNMRDistanceInc;                    // NMR distance increment for step weighting
    PMEDouble2*                 pNMRDistanceR1R2Slp;                // NMR distance r1, r2 slope
    PMEDouble2*                 pNMRDistanceR3R4Slp;                // NMR distance r3, r4 slope
    PMEDouble2*                 pNMRDistanceK2K3Slp;                // NMR distance k2, k3 slope
    PMEDouble2*                 pNMRDistanceR1R2Int;                // NMR distance r1, r2 intercept
    PMEDouble2*                 pNMRDistanceR3R4Int;                // NMR distance r3, r4 intercept
    PMEDouble2*                 pNMRDistanceK2K3Int;                // NMR distance k2, k3 intercept
    int4*                       pNMRAngleID1;                       // NMR angle i, j, k, ibuff
    int2*                       pNMRAngleID2;                       // NMR angle jbuff, kbuff
    PMEDouble2*                 pNMRAngleR1R2;                      // NMR angle computed r1, r2
    PMEDouble2*                 pNMRAngleR3R4;                      // NMR angle computed r3, r4
    PMEDouble2*                 pNMRAngleK2K3;                      // NMR angle computed k2, k3
    PMEDouble*                  pNMRAngleAve;                       // NMR angle restraint linear and exponential averaged value
    PMEDouble2*                 pNMRAngleTgtVal;                    // NMR angle target and value for current step
    int2*                       pNMRAngleStep;                      // NMR angle first and last step for application of restraint
    int*                        pNMRAngleInc;                       // NMR angle increment for step weighting
    PMEDouble2*                 pNMRAngleR1R2Slp;                   // NMR angle r1, r2 slope
    PMEDouble2*                 pNMRAngleR3R4Slp;                   // NMR angle r3, r4 slope
    PMEDouble2*                 pNMRAngleK2K3Slp;                   // NMR angle k2, k3 slope
    PMEDouble2*                 pNMRAngleR1R2Int;                   // NMR angle r1, r2 intercept
    PMEDouble2*                 pNMRAngleR3R4Int;                   // NMR angle r3, r4 intercept
    PMEDouble2*                 pNMRAngleK2K3Int;                   // NMR angle k2, k3 intercept
    int4*                       pNMRTorsionID1;                     // NMR torsion i, j, k, l
    int4*                       pNMRTorsionID2;                     // NMR torsion ibuff, jbuff, kbuff, lbuff
    PMEDouble2*                 pNMRTorsionR1R2;                    // NMR torsion computed r1, r2
    PMEDouble2*                 pNMRTorsionR3R4;                    // NMR torsion computed r3, r4
    PMEDouble2*                 pNMRTorsionK2K3;                    // NMR torsion computed k2, k3
    PMEDouble*                  pNMRTorsionAve1;                    // NMR torsion restraint linear and exponential averaged value
    PMEDouble*                  pNMRTorsionAve2;                    // NMR torsion restraint linear and exponential averaged value
    PMEDouble2*                 pNMRTorsionTgtVal;                  // NMR torsion target and value for current step    
    int2*                       pNMRTorsionStep;                    // NMR torsion first and last step for application of restraint
    int*                        pNMRTorsionInc;                     // NMR torsion increment for step weighting
    PMEDouble2*                 pNMRTorsionR1R2Slp;                 // NMR torsion r1, r2 slope
    PMEDouble2*                 pNMRTorsionR3R4Slp;                 // NMR torsion r3, r4 slope
    PMEDouble2*                 pNMRTorsionK2K3Slp;                 // NMR torsion k2, k3 slope
    PMEDouble*                  pNMRTorsionK4Slp;                   // NMR torsion k4 slope
    PMEDouble2*                 pNMRTorsionR1R2Int;                 // NMR torsion r1, r2 intercept
    PMEDouble2*                 pNMRTorsionR3R4Int;                 // NMR torsion r3, r4 intercept
    PMEDouble2*                 pNMRTorsionK2K3Int;                 // NMR torsion k2, k3 intercept    


    // Cmap lookup table
    int                         cmapTermStride;                     // Per term stride 
    int                         cmapRowStride;                      // Per row of terms stride
    PMEFloat4*                  pCmapEnergy;                        // Pointer to Cmap LUT data (E, dPhi, dPsi, dPhi_dPsi)
    
    // Accumulation buffers
#ifdef use_SPFP
    PMEAccumulator*             pForceAccumulator;                  // Bare fixed point force accumulator   
    PMEAccumulator*             pForceXAccumulator;                 // Bare fixed point x force accumulator   
    PMEAccumulator*             pForceYAccumulator;                 // Bare fixed point y force accumulator   
    PMEAccumulator*             pForceZAccumulator;                 // Bare fixed point z force accumulator   
    PMEAccumulator*             pNBForceAccumulator;                // Bare fixed point nonbond force accumulator
    PMEAccumulator*             pNBForceXAccumulator;               // Fixed point nonbond x force accumulator   
    PMEAccumulator*             pNBForceYAccumulator;               // Fixed point nonbond y force accumulator   
    PMEAccumulator*             pNBForceZAccumulator;               // Fixed point nonbond z force accumulator   
    PMEAccumulator*             pReffAccumulator;                   // Effective Born Radius buffer
    PMEAccumulator*             pSumdeijdaAccumulator;              // Atom Sumdeijda buffer
#else    
    PMEDouble*                  pForceBuffer;                       // Bare force buffer pointer (for bonded components)
    PMEDouble*                  pForceXBuffer;                      // Atom X force buffer
    PMEDouble*                  pForceYBuffer;                      // Atom Y force buffer
    PMEDouble*                  pForceZBuffer;                      // Atom Z force buffer    
    PMEDouble*                  pReffBuffer;                        // Effective Born Radius buffer
    PMEDouble*                  pSumdeijdaBuffer;                   // Atom Sumdeijda buffer
#endif
    PMEAccumulator*             pBondedForceAccumulator;            // Bare fixed point bonded force accumulator   
    PMEAccumulator*             pBondedForceXAccumulator;           // Bare fixed point bonded x force accumulator   
    PMEAccumulator*             pBondedForceYAccumulator;           // Bare fixed point bonded y force accumulator   
    PMEAccumulator*             pBondedForceZAccumulator;           // Bare fixed point bonded z force accumulator
    int                         stride;                             // Atom quantity stride
    int                         stride2;                            // Atom quantity 2x stride
    int                         stride3;                            // Atom quantity 3x stride
    int                         stride4;                            // Atom quantity 4x stride
    int                         imageStride;                        // Neighor list index stride (1K-aligned) for Duane Merrill's radix sort
    
    
    // Kernel call stuff
    unsigned int                workUnits;                          // Total work units
    unsigned int                excludedWorkUnits;                  // Total work units with exclusions
    unsigned int*               pExclusion;                         // Exclusion masks                
    unsigned int*               pWorkUnit;                          // Work unit list
    unsigned int*               pGBBRPosition;                      // Generalized Born Born Radius workunit position
    unsigned int*               pGBNB1Position;                     // Generalized Born Nonbond 1 workunit position
    unsigned int*               pGBNB2Position;                     // Generalized Born Nonbond 2 workunit position
    unsigned int                GBTotalWarps[3];                    // Total warps in use for nonbond kernels
    unsigned int                maxForceBuffers;                    // Total Force buffer count
    unsigned int                nonbondForceBuffers;                // Nonbond Force buffer count
    PMEFloat                    cellOffset[NEIGHBORCELLS][3];       // Local coordinate system offsets

};


template <typename T> struct GpuBuffer;


struct _gpuContext {
    SM_VERSION                  sm_version;    
    
    // Memory parameters
    bool                        bECCSupport;                // Flag for ECC support to detect Tesla versus consumer GPU
    aligned_lli                 totalMemory;                // Total memory on GPU
    aligned_lli                 totalCPUMemory;             // Approximate total allocated CPU memory
    aligned_lli                 totalGPUMemory;             // Approximate total allocated CPU memory
        
    // AMBER parameters
    int                         ntt;                        // AMBER Thermostat setting
    int                         ntb;                        // AMBER PBC setting
    int                         ips;                        // AMBER IPS setting
    int                         ntc;                        // AMBER SHAKE setting
    int                         ntf;                        // AMBER force field setting
    int                         ntpr;                       // AMBER status output interval
    int                         ntwe;                       // AMBER energy output interval
    int                         ntr;                        // AMBER Restraint flags

    // Atom stuff
    GpuBuffer<PMEDouble>*       pbAtom;                     // Atom coordinates
    GpuBuffer<PMEFloat2>*       pbAtomXYSP;                 // Single Precision Atom X and Y coordinates    
    GpuBuffer<PMEFloat>*        pbAtomZSP;                  // Single Precision Atom Z coordinate
    GpuBuffer<PMEFloat2>*       pbAtomSigEps;               // Atom nonbond parameters
    GpuBuffer<PMEFloat>*        pbAtomS;                    // Atom scaled Born Radius
    GpuBuffer<PMEFloat>*        pbAtomRBorn;                // Atom Born Radius
    GpuBuffer<PMEDouble>*       pbAtomCharge;               // Atom charges
    GpuBuffer<PMEFloat>*        pbAtomChargeSP;             // Single precision atom charges
    GpuBuffer<PMEDouble>*       pbAtomMass;                 // Atom masses
    GpuBuffer<PMEDouble>*       pbReff;                     // Effective Born Radius
    GpuBuffer<PMEFloat>*        pbReffSP;                   // Single Precision Effective Born Radius
    GpuBuffer<PMEFloat>*        pbTemp7;                    // Single Precision Born Force
    GpuBuffer<PMEFloat>*        pbPsi;                      // GB intermediate psi value
    GpuBuffer<PMEDouble>*       pbForce;                    // Atom forces 
    GpuBuffer<PMEAccumulator>*  pbBondedForce;              // Fixed point bonded forces for NTP runs
#ifdef MPI
    GpuBuffer<PMEDouble>*       pbReffa;                    // Pinned memory unprocessed Born Radius
    GpuBuffer<PMEDouble>*       pbTemp7a;                   // Pinned memory unprocessed Born Force
    GpuBuffer<PMEDouble>*       pbOutForce;                 // Interleaved forces for transmission in MPI mode
    GpuBuffer<PMEDouble>*       pbInForce;                  // Zero copy forces for updating in MPI mode
    GpuBuffer<PMEFloat>*        pbPMEForce;                 // PME forces on node 0 in MPI mode
#endif    
    GpuBuffer<PMEDouble>*       pbVel;                      // Atom velocities
    GpuBuffer<PMEDouble>*       pbLVel;                     // Atom last velocities
    GpuBuffer<unsigned int>*    pbOutputBufferCounter;      // Output buffer counter for bonded interactions
    GpuBuffer<PMEFloat>*        pbCenter;                   // Atom recentering extrema   
    
    // GB Data
    GpuBuffer<PMEFloat2>*       pbNeckMaxValPos;            // GB Born Radii and Energy correction data
    GpuBuffer<PMEFloat>*        pbGBAlphaBetaGamma;         // GB per atom Alpha, Beta, and Gamma
    
    // PME data
    cufftHandle                 forwardPlan;                // PME FFT plan for cuFFT
    cufftHandle                 backwardPlan;               // PME IFFT plan for cuFFT
#ifdef use_SPFP
    GpuBuffer<long long int>*   pblliXYZ_q;                 // SPFP long long int PME charge grid
#endif
    GpuBuffer<PMEFloat>*        pbXYZ_q;                    // PME charge grid/buffer
    GpuBuffer<PMEComplex>*      pbXYZ_qt;                   // FFTed PME charge grid
    GpuBuffer<PMEFloat>*        pbPrefac;                   // PME FFT pre-factors
    GpuBuffer<PMEFloat>*        pbFraction;                 // PME fractional coordinates
    GpuBuffer<int>*             pbIFract;                   // PME Unit cell coordinates
    GpuBuffer<PMEFloat4>*       pbTheta;                    // PME spline coefficients
    GpuBuffer<int2>*            pbTileBoundary;             // PME charge interpolation box boundaries (min/max)

    // Self energy values, calculated on host
    PMEDouble                   ee_plasma;                  // Component of PME self energy
    PMEDouble                   self_energy;                // Total self energy                    
    PMEDouble                   vdw_recip;                  // Reciprocal vdw correction energy       
 
    // Neighbor list stuff
    bool                        bNeighborList;              // Neighbor list activation flag
    bool                        bNewNeighborList;           // Newly generated neighbor list activation flag
    bool                        bNeedNewNeighborList;       // Need newly generated neighbor list activation flag
    bool                        bOddNLCells;                // Flag to determine whether charge grid needs 8 or 16 output buffers
    bool                        bSmallBox;                  // Flag to determine whether any nonbond cell count is 1 or 2 on any axis
    unsigned int                neighborListBits;           // Number of bits in cell list
    GpuBuffer<unsigned int>*    pbImageIndex;               // Image indexing data
    GpuBuffer<PMEFloat2>*       pbAtomXYSaveSP;             // Saved neighbor list coordinates
    GpuBuffer<PMEFloat>*        pbAtomZSaveSP;              // Saved neighbor list coordinates
    GpuBuffer<PMEDouble>*       pbImage;                    // Image coordinates
    GpuBuffer<PMEDouble>*       pbImageVel;                 // Image velocities
    GpuBuffer<PMEDouble>*       pbImageLVel;                // Image last velocities
    GpuBuffer<PMEDouble>*       pbImageMass;                // Image masses
    GpuBuffer<PMEDouble>*       pbImageCharge;              // Image charges
    GpuBuffer<PMEFloat2>*       pbImageSigEps;              // Image sigma/epsilon data for nonbond interactions
    GpuBuffer<unsigned int>*    pbImageOutputBuffers;       // Image output buffer count for all interactions
    GpuBuffer<unsigned int>*    pbImageCellID;              // Image per atom nonbond cell ID
    GpuBuffer<unsigned int>*    pbNLExclusionList;          // Raw exclusion list
    GpuBuffer<uint2>*           pbNLExclusionStartCount;    // Per atom exclusion bounds
    GpuBuffer<uint2>*           pbNLNonbondCellStartEnd;    // Nonbond cells start and end
#ifndef use_SPFP
    GpuBuffer<int>*             pbNLChargeGridBufferOffset; // Precalculated per cell charge grid buffer
    GpuBuffer<int>*             pbNLOddBufferOverlapFlag;   // Per axis grid cell odd buffer overlap flag
#endif
    GpuBuffer<bool>*            pbNLbSkinTestFail;          // Skin test result buffer
    GpuBuffer<unsigned int>*    pbNLCellHash;               // Spatial ordering hash for within cells
    GpuBuffer<NLRecord>*        pbNLRecord;                 // Pointer to neighbor list records
    GpuBuffer<NLEntry>*         pbNLEntry;                  // Active neighbor list
    GpuBuffer<unsigned int>*    pbNLAtomList;               // Pointer to atom list and exclusion data
    GpuBuffer<unsigned int>*    pbNLOffset;                 // Offsets to exclusions/atom list 
    GpuBuffer<unsigned int>*    pbNLTotalOffset;            // Current atom list offset
    GpuBuffer<unsigned int>*    pbNLPosition;               // Position in building neighbor list  
    GpuBuffer<unsigned int>*    pbBNLExclusionBuffer;       // Per-warp GMEM exclusion buffer  
    PMEFloat                    nonbond_skin;               // Nonbond skin thickness

    // Remapped bonded interactions and Shake constraint atom IDs
    GpuBuffer<int4>*            pbImageBondID;              // Remapped Bond i, j, ibuff, jbuff
    GpuBuffer<int4>*            pbImageBondAngleID1;        // Remapped Bond Angle i, j, k, ibuff
    GpuBuffer<int2>*            pbImageBondAngleID2;        // Remapped Bond Angle jbuff, kbuff;
    GpuBuffer<int4>*            pbImageDihedralID1;         // Remapped Dihedral i, j, k, l
    GpuBuffer<int4>*            pbImageDihedralID2;         // Remapped Dihedral ibuff, jbuff, kbuff, lbuff
    GpuBuffer<int4>*            pbImageNb14ID;              // Remapped 1-4 nonbond i, j, ibuff, jbuff                
    GpuBuffer<int2>*            pbImageConstraintID;        // Remapped Atom constraint ID
    GpuBuffer<int4>*            pbImageUBAngleID;           // Remapped Urey Bradley Angle i, j, ibuff, jbuff
    GpuBuffer<int4>*            pbImageImpDihedralID1;      // Remapped Improper Dihedral i, j, k, l 
    GpuBuffer<int4>*            pbImageImpDihedralID2;      // Remapped Improper Dihedral ibuff, jbuff, kbuff, lbuff
    GpuBuffer<int4>*            pbImageCmapID1;             // Remapped Cmap i, j, k, l
    GpuBuffer<int4>*            pbImageCmapID2;             // Remapped Cmap m, ibuff, jbuff, kbuff
    GpuBuffer<int2>*            pbImageCmapID3;             // Remapped Cmap lbuff, mbuff
    GpuBuffer<int4>*            pbImageNMRDistanceID;       // Remapped NMR Distance i, j, ibuff, jbuff
    GpuBuffer<int4>*            pbImageNMRAngleID1;         // Remapped NMR Angle i, j, k, ibuff
    GpuBuffer<int2>*            pbImageNMRAngleID2;         // Remapped NMR Angle jbuff, kbuff;
    GpuBuffer<int4>*            pbImageNMRTorsionID1;       // Remapped NMR Torsion i, j, k, l
    GpuBuffer<int4>*            pbImageNMRTorsionID2;       // Remapped NMR Torsion ibuff, jbuff, kbuff, lbuff
    GpuBuffer<int4>*            pbImageShakeID;             // Remapped traditional SHAKE IDs  
    GpuBuffer<int4>*            pbImageFastShakeID;         // Remapped H2O (Fast) SHAKE IDs  
    GpuBuffer<int>*             pbImageSlowShakeID1;        // Remapped Central atom of XH4 (Slow) Shake constraint 
    GpuBuffer<int4>*            pbImageSlowShakeID2;        // Remapped XH4 (Slow) SHAKE constraint hydrogens
    GpuBuffer<int4>*            pbImageSolventAtomID;       // Remapped solvent molecules/ions
    GpuBuffer<int>*             pbImageSoluteAtomID;        // Remapped solute atoms
    
    // NTP molecule data
    int                         maxSoluteMolecules;         // Maximum solute molecules for most NTP kernels
    int                         maxPSSoluteMolecules;       // Maximum solute molecules for NTP pressure scaling kernels
    GpuBuffer<int>*             pbSoluteAtomID;             // Solute per atom ID 
    GpuBuffer<PMEDouble>*       pbSoluteAtomMass;           // Solute atom masses
    GpuBuffer<PMEDouble>*       pbSolute;                   // Current and last centers of mass, kinetic energy and inverse mass for each solute molecule
    GpuBuffer<PMEUllInt>*       pbUllSolute;                // Current COM and EKCOM in integer form
    GpuBuffer<int4>*            pbSolventAtomID;            // List of solvent molecules/ions of 4 or fewer atoms
    GpuBuffer<PMEDouble>*       pbSolvent;                  // Last centers of mass, atom masses and inverse mass for each solvent molecule
    GpuBuffer<NTPData>*         pbNTPData;                  // NTP mutable values
    
    // NTP constraint molecule data
    GpuBuffer<PMEDouble>*       pbConstraintAtomX;          // Original atom x for constraint
    GpuBuffer<PMEDouble>*       pbConstraintAtomY;          // Original atom y for constraint
    GpuBuffer<PMEDouble>*       pbConstraintAtomZ;          // Original atom z for constraint
    GpuBuffer<PMEDouble>*       pbConstraintCOMX;           // Original x COM for constraint in fractional coordinates
    GpuBuffer<PMEDouble>*       pbConstraintCOMY;           // Original y COM for constraint in fractional coordinates
    GpuBuffer<PMEDouble>*       pbConstraintCOMZ;           // Original z COM for constraint in fractional coordinates
    
    // Radix sort data
    GpuBuffer<unsigned int>*    pbRadixBlockSum;            // Radix sums per thread group
    GpuBuffer<unsigned int>*    pbRadixCounter;             // Radix counters per thread group for scatter
    GpuBuffer<unsigned int>*    pbRadixSum;                 // Pointer to indivudal thread block radix sums

    // Random number stuff
    unsigned int                randomCounter;              // Counter for triggering RNG
    GpuBuffer<unsigned int>*    pbRandomPos;                // RNG position per block
    GpuBuffer<PMEDouble>*       pbRandom;                   // Pointer to overall RNG buffer
    curandGenerator_t           RNG;                        // CURAND RNG 
 
    // Bonded interactions
    bool                        bLocalInteractions;         // Flag indicating presence or absence of local interactions
    bool                        bCharmmInteractions;        // Flag indicating presence or absence of CHARMM interactions
    GpuBuffer<PMEDouble2>*      pbBond;                     // Bond Kr, Req
    GpuBuffer<int4>*            pbBondID;                   // Bond i, j, ibuff, jbuff
    GpuBuffer<PMEDouble2>*      pbBondAngle;                // Bond Angle Ka, Aeq
    GpuBuffer<int4>*            pbBondAngleID1;             // Bond Angle i, j, k, ibuff
    GpuBuffer<int2>*            pbBondAngleID2;             // Bond Angle jbuff, kbuff;
    GpuBuffer<PMEDouble2>*      pbDihedral1;                // Dihedral gmul, pn
    GpuBuffer<PMEDouble2>*      pbDihedral2;                // Dihedral pk, gamc
    GpuBuffer<PMEDouble>*       pbDihedral3;                // Dihedral gams
    GpuBuffer<int4>*            pbDihedralID1;              // Dihedral i, j, k, l
    GpuBuffer<int4>*            pbDihedralID2;              // Dihedral ibuff, jbuff, kbuff, lbuff
    GpuBuffer<PMEDouble2>*      pbNb141;                    // 1-4 nonbond scee, scnb0
    GpuBuffer<PMEDouble2>*      pbNb142;                    // 1-4 nonbond cn1, cn2
    GpuBuffer<int4>*            pbNb14ID;                   // 1-4 nonbond i, j, ibuff, jbuff
    GpuBuffer<PMEDouble2>*      pbConstraint1;              // Constraint weight and xc
    GpuBuffer<PMEDouble2>*      pbConstraint2;              // Constraint yc and zc                      
    GpuBuffer<int2>*            pbConstraintID;             // Atom constraint ID
    GpuBuffer<PMEDouble2>*      pbUBAngle;                  // Urey Bradley Angle rk, r0
    GpuBuffer<int4>*            pbUBAngleID;                // Urey Bradley Angle i, j, ibuff, jbuff
    GpuBuffer<PMEDouble2>*      pbImpDihedral;              // Improper Dihedral pk, phase
    GpuBuffer<int4>*            pbImpDihedralID1;           // Improper Dihedral i, j, k, l 
    GpuBuffer<int4>*            pbImpDihedralID2;           // Improper Dihedral ibuff, jbuff, kbuff, lbuff
    GpuBuffer<int4>*            pbCmapID1;                  // Cmap i, j, k, l
    GpuBuffer<int4>*            pbCmapID2;                  // Cmap m, ibuff, jbuff, kbuff
    GpuBuffer<int2>*            pbCmapID3;                  // Cmap lbuff, mbuff
    GpuBuffer<int>*             pbCmapType;                 // Cmap type
    GpuBuffer<PMEFloat4>*       pbCmapEnergy;               // Pointer to Cmap LUT data (E, dPhi, dPsi, dPhi_dPsi)
    
    // NMR stuff
    bool                        bNMRInteractions;           // Flag indicating presence or absence of NMR interactions
    int                         NMRnstep;                   // Imported NMR variable for time-dependent restraints
    GpuBuffer<double>*          pbNMRJarData;               // Jarzynski accumulated work data
    GpuBuffer<int4>*            pbNMRDistanceID;            // NMR distance i, j, ibuff, jbuff
    GpuBuffer<PMEDouble2>*      pbNMRDistanceR1R2;          // NMR distance computed r1, r2
    GpuBuffer<PMEDouble2>*      pbNMRDistanceR3R4;          // NMR distance computed r3, r4
    GpuBuffer<PMEDouble2>*      pbNMRDistanceK2K3;          // NMR distance computed k2, k3
    GpuBuffer<PMEDouble>*       pbNMRDistanceK4;            // NMR distance computed k4
    GpuBuffer<PMEDouble>*       pbNMRDistanceAve;           // NMR distance restraint linear and exponential averaged value
    GpuBuffer<PMEDouble2>*      pbNMRDistanceTgtVal;        // NMR distance target and actual value for current step
    GpuBuffer<int2>*            pbNMRDistanceStep;          // NMR distance first and last step for application of restraint
    GpuBuffer<int>*             pbNMRDistanceInc;           // NMR distance increment for step weighting
    GpuBuffer<PMEDouble2>*      pbNMRDistanceR1R2Slp;       // NMR distance r1, r2 slope
    GpuBuffer<PMEDouble2>*      pbNMRDistanceR3R4Slp;       // NMR distance r3, r4 slope
    GpuBuffer<PMEDouble2>*      pbNMRDistanceK2K3Slp;       // NMR distance k2, k3 slope
    GpuBuffer<PMEDouble>*       pbNMRDistanceK4Slp;         // NMR distance k4 slope
    GpuBuffer<PMEDouble2>*      pbNMRDistanceR1R2Int;       // NMR distance r1, r2 intercept
    GpuBuffer<PMEDouble2>*      pbNMRDistanceR3R4Int;       // NMR distance r3, r4 intercept
    GpuBuffer<PMEDouble2>*      pbNMRDistanceK2K3Int;       // NMR distance k2, k3 intercept
    GpuBuffer<PMEDouble>*       pbNMRDistanceK4Int;         // NMR distance k4 intercept   
    GpuBuffer<int4>*            pbNMRAngleID1;              // NMR angle i, j, k, ibuff
    GpuBuffer<int2>*            pbNMRAngleID2;              // NMR angle jbuff, kbuff
    GpuBuffer<PMEDouble2>*      pbNMRAngleR1R2;             // NMR angle computed r1, r2
    GpuBuffer<PMEDouble2>*      pbNMRAngleR3R4;             // NMR angle computed r3, r4
    GpuBuffer<PMEDouble2>*      pbNMRAngleK2K3;             // NMR angle computed k2, k3
    GpuBuffer<PMEDouble>*       pbNMRAngleK4;               // NMR angle computed k4
    GpuBuffer<PMEDouble>*       pbNMRAngleAve;              // NMR angle restraint linear and exponential averaged value
    GpuBuffer<PMEDouble2>*      pbNMRAngleTgtVal;           // NMR angle target and actual value for current step
    GpuBuffer<int2>*            pbNMRAngleStep;             // NMR angle first and last step for application of restraint
    GpuBuffer<int>*             pbNMRAngleInc;              // NMR angle increment for step weighting
    GpuBuffer<PMEDouble2>*      pbNMRAngleR1R2Slp;          // NMR angle r1, r2 slope
    GpuBuffer<PMEDouble2>*      pbNMRAngleR3R4Slp;          // NMR angle r3, r4 slope
    GpuBuffer<PMEDouble2>*      pbNMRAngleK2K3Slp;          // NMR angle k2, k3 slope
    GpuBuffer<PMEDouble>*       pbNMRAngleK4Slp;            // NMR angle k4 slope
    GpuBuffer<PMEDouble2>*      pbNMRAngleR1R2Int;          // NMR angle r1, r2 intercept
    GpuBuffer<PMEDouble2>*      pbNMRAngleR3R4Int;          // NMR angle r3, r4 intercept
    GpuBuffer<PMEDouble2>*      pbNMRAngleK2K3Int;          // NMR angle k2, k3 intercept
    GpuBuffer<PMEDouble>*       pbNMRAngleK4Int;            // NMR angle k4 intercept   
    GpuBuffer<int4>*            pbNMRTorsionID1;            // NMR torsion i, j, k, l
    GpuBuffer<int4>*            pbNMRTorsionID2;            // NMR torsion ibuff, jbuff, kbuff, lbuff
    GpuBuffer<PMEDouble2>*      pbNMRTorsionR1R2;           // NMR torsion computed r1, r2
    GpuBuffer<PMEDouble2>*      pbNMRTorsionR3R4;           // NMR torsion computed r3, r4
    GpuBuffer<PMEDouble2>*      pbNMRTorsionK2K3;           // NMR torsion computed k2, k3
    GpuBuffer<PMEDouble>*       pbNMRTorsionK4;             // NMR torsion computed k4
    GpuBuffer<PMEDouble>*       pbNMRTorsionAve1;           // NMR torsion restraint linear and exponential averaged value
    GpuBuffer<PMEDouble>*       pbNMRTorsionAve2;           // NMR torsion restraint linear and exponential averaged value
    GpuBuffer<PMEDouble2>*      pbNMRTorsionTgtVal;         // NMR torsion target and actual value for current step
    GpuBuffer<int2>*            pbNMRTorsionStep;           // NMR torsion first and last step for application of restraint
    GpuBuffer<int>*             pbNMRTorsionInc;            // NMR torsion increment for step weighting
    GpuBuffer<PMEDouble2>*      pbNMRTorsionR1R2Slp;        // NMR torsion r1, r2 slope
    GpuBuffer<PMEDouble2>*      pbNMRTorsionR3R4Slp;        // NMR torsion r3, r4 slope
    GpuBuffer<PMEDouble2>*      pbNMRTorsionK2K3Slp;        // NMR torsion k2, k3 slope
    GpuBuffer<PMEDouble>*       pbNMRTorsionK4Slp;          // NMR torsion k4 slope
    GpuBuffer<PMEDouble2>*      pbNMRTorsionR1R2Int;        // NMR torsion r1, r2 intercept
    GpuBuffer<PMEDouble2>*      pbNMRTorsionR3R4Int;        // NMR torsion r3, r4 intercept
    GpuBuffer<PMEDouble2>*      pbNMRTorsionK2K3Int;        // NMR torsion k2, k3 intercept
    GpuBuffer<PMEDouble>*       pbNMRTorsionK4Int;          // NMR torsion k4 intercept           

#ifdef use_SPFP
    // Force accumulators for atomic ops
    GpuBuffer<PMEAccumulator>*  pbReffAccumulator;          // Effective Born Radius accumulator
    GpuBuffer<PMEAccumulator>*  pbSumdeijdaAccumulator;     // Sumdeijda accumulator
    GpuBuffer<PMEAccumulator>*  pbForceAccumulator;         // Force accumulators
#else    
    // Force output buffers
    GpuBuffer<PMEDouble>*       pbReffBuffer;               // Effective Born Radius output buffer
    GpuBuffer<PMEDouble>*       pbSumdeijdaBuffer;          // Sumdeijda output buffer
    GpuBuffer<PMEDouble>*       pbForceBuffer;              // Force output buffer
#endif

    GpuBuffer<unsigned long long int>*   pbEnergyBuffer;    // Energy accumulation buffer
    GpuBuffer<KineticEnergy>*   pbKineticEnergyBuffer;      // Kinetic energy accumulation buffer
    
    // AMD buffers
    PMEDouble*                  pAmdWeightsAndEnergy;       // AMD
    GpuBuffer<PMEDouble>*       pbAMDfwgtd;

    // Extra points data
    GpuBuffer<int4>*            pbExtraPoint11Frame;        // Type 1 single extra point parent_atm, atm1, atm2, atm3
    GpuBuffer<int>*             pbExtraPoint11Index;        // Type 1 single extra point EP
    GpuBuffer<PMEDouble>*       pbExtraPoint11;             // Type 1 single extra point coordinates
    GpuBuffer<int4>*            pbExtraPoint12Frame;        // Type 2 single extra point parent_atm, atm1, atm2, atm3
    GpuBuffer<int>*             pbExtraPoint12Index;        // Type 2 single extra point EP
    GpuBuffer<PMEDouble>*       pbExtraPoint12;             // Type 2 single extra point coordinates
    GpuBuffer<int4>*            pbExtraPoint21Frame;        // Type 1 dual extra point parent_atm, atm1, atm2, atm3
    GpuBuffer<int2>*            pbExtraPoint21Index;        // type 1 dual extra point EP1, EP2
    GpuBuffer<PMEDouble>*       pbExtraPoint21;             // Type 1 dual extra point coordinates
    GpuBuffer<int4>*            pbExtraPoint22Frame;        // Type 2 dual extra point parent_atm, atm1, atm2, atm3
    GpuBuffer<int2>*            pbExtraPoint22Index;        // Type 2 dual extra point EP1. EP2
    GpuBuffer<PMEDouble>*       pbExtraPoint22;             // Type 2 dual extra point coordinates
    GpuBuffer<int4>*            pbImageExtraPoint11Frame;   // Remapped Type 1 single extra point parent_atm, atm1, atm2, atm3
    GpuBuffer<int>*             pbImageExtraPoint11Index;   // Remapped Type 1 single extra point EP
    GpuBuffer<int4>*            pbImageExtraPoint12Frame;   // Remapped Type 2 single extra point parent_atm, atm1, atm2, atm3
    GpuBuffer<int>*             pbImageExtraPoint12Index;   // Remapped Type 2 single extra point EP
    GpuBuffer<int4>*            pbImageExtraPoint21Frame;   // Remapped Type 1 dual extra point parent_atm, atm1, atm2, atm3
    GpuBuffer<int2>*            pbImageExtraPoint21Index;   // Remapped type 1 dual extra point EP1, EP2
    GpuBuffer<int4>*            pbImageExtraPoint22Frame;   // Remapped Type 2 dual extra point parent_atm, atm1, atm2, atm3
    GpuBuffer<int2>*            pbImageExtraPoint22Index;   // Remapped Type 2 dual extra point EP1. EP2
    
    // SHAKE constraints
    GpuBuffer<int4>*            pbShakeID;                  // SHAKE central Atom plus up to 3 hydrogens
    GpuBuffer<double2>*         pbShakeParm;                // SHAKE central atom mass and equilibrium bond length
    GpuBuffer<int4>*            pbFastShakeID;              // H2O oxygen plus two hydrogens atom ID
    GpuBuffer<int>*             pbSlowShakeID1;             // Central atom of XH4 Shake constraint 
    GpuBuffer<int4>*            pbSlowShakeID2;             // XH4 SHAKE constraint hydrogens
    GpuBuffer<double2>*         pbSlowShakeParm;            // XH4 SHAKE central atom mass and equilibrium bond length

    
    // Launch parameters
    int                         gpu_device_id;              // allow -1 for default = choose gpu based on memory.
    unsigned int                blocks;
    unsigned int                threadsPerBlock;
    unsigned int                clearForcesThreadsPerBlock; 
    unsigned int                NLClearForcesThreadsPerBlock; 
    unsigned int                reduceForcesThreadsPerBlock; 
    unsigned int                NLReduceForcesThreadsPerBlock; 
    unsigned int                reduceBufferThreadsPerBlock;
    unsigned int                localForcesBlocks;          
    unsigned int                localForcesThreadsPerBlock;
    unsigned int                CHARMMForcesThreadsPerBlock;
    unsigned int                NMRForcesThreadsPerBlock;
    unsigned int                GBBornRadiiThreadsPerBlock;
    unsigned int                GBBornRadiiIGB78ThreadsPerBlock;
    unsigned int                GBNonbondEnergy1ThreadsPerBlock;
    unsigned int                GBNonbondEnergy2ThreadsPerBlock;
    unsigned int                GBNonbondEnergy2IGB78ThreadsPerBlock;
    unsigned int                PMENonbondEnergyThreadsPerBlock; 
    unsigned int                IPSNonbondEnergyThreadsPerBlock;
    unsigned int                updateThreadsPerBlock;
    unsigned int                shakeThreadsPerBlock;
    unsigned int                GBBornRadiiBlocks;
    unsigned int                GBNonbondEnergy1Blocks;
    unsigned int                GBNonbondEnergy2Blocks;
    unsigned int                PMENonbondBlocks;
    unsigned int                IPSNonbondBlocks;
    unsigned int                BNLBlocks;
    unsigned int                NLCalculateOffsetsThreadsPerBlock;
    unsigned int                NLBuildNeighborList32ThreadsPerBlock;
    unsigned int                NLBuildNeighborList16ThreadsPerBlock;
    unsigned int                readSize;
    
    // Nonbond Kernel stuff
    GpuBuffer<unsigned int>*    pbExclusion;                // Exclusion masks                
    GpuBuffer<unsigned int>*    pbWorkUnit;                 // Work unit list
    GpuBuffer<unsigned int>*    pbGBPosition;               // Nonbond kernel worunit positions
    
    cudaSimulation              sim;                        // All simulation data destined for GPU constant RAM


#ifdef MPI
    // Single-node multi-gpu parameters
    PMEDouble*                  pSharedForceBuffer0;        // Shared memory Force Buffer 0    
    PMEDouble*                  pSharedForceBuffer1;        // Shared memory Force Buffer 1    
    PMEDouble*                  pSharedForceBuffer2;        // Shared memory Force Buffer 2    
    PMEFloat*                   pSharedPMEForceBuffer0;     // Shared memory PME Force Buffer 0
    PMEFloat*                   pSharedPMEForceBuffer1;     // Shared memory PME Force Buffer 1
    PMEFloat*                   pSharedPMEForceBuffer2;     // Shared memory PME Force Buffer 2
    bool                        bSingleNode;                // Flag to indicate MPI run is all on one node
#ifdef CUDA_P2P
    bool                        bP2P;                       // Flag to trigger P2P force transmission
    PMEDouble**                 pP2PReffa;                  // Pointers to per-process Born Radius Buffers
    PMEDouble**                 pP2PTemp7a;                 // Pointers to per-process Temp7a Buffers 
    PMEDouble**                 pP2PInForce;                // Pointers to per-process Force Buffers;
    cudaIpcMemHandle_t*         pP2PReffaHandle;            // Pointers to per-process Born Radius Buffers
    cudaIpcMemHandle_t*         pP2PTemp7aHandle;           // Pointers to per-process Temp7a Buffers 
    cudaIpcMemHandle_t*         pP2PInForceHandle;          // Pointers to per-process Force Buffers;
#endif

    // Multi-gpu parameters
    int                         nGpus;                      // Number of GPUs involved in calculation
    int                         gpuID;                      // Local GPU ID
    MPI_Comm                    comm;                       // MPI Communicator for all collective MPI calls
    int                         minLocalCell;               // First nonbond cell owned by local GPU
    int                         maxLocalCell;               // End of nonbond cells owned by local GPU
    int                         minProcessedCell;           // First cell requiring force reduction
    int                         maxProcessedCell;           // Last cell requiring force reduction
    int*                        pMinLocalCell;              // First cell owned by each node
    int*                        pMaxLocalCell;              // Last cell owned by each node
    int*                        pMinLocalAtom;              // First atom owned by each node
    int*                        pMaxLocalAtom;              // Last atom owned by each node    
    int*                        pMinProcessedCell;          // First cell in each node requiring force reduction
    int*                        pMaxProcessedCell;          // Last cell in each node requiring force reduction
    int*                        pMinProcessedAtom;          // First atom in each node requiring force reduction
    int*                        pMaxProcessedAtom;          // Last atom in each node requiring force reduction
    int*                        pAllGathervRecvCountSoA;    // Length of first allGatherv per node
    int*                        pAllGathervRecvDisplSoA;    // Displacement of first allGatherv per node
    int*                        pAllGathervRecvCountAoS;    // Length of first allGatherv per node
    int*                        pAllGathervRecvDisplAoS;    // Displacement of first allGatherv per node    
    PMEFloat*                   pPMEForceData;              // Buffer for incoming reciprocal force data
    int*                        pPMEStart;                  // Starts for PME data on each node
    int*                        pPMELength;                 // Length of PME data on each node
    int                         forceSendNodes;             // Number of nodes to which force data is written
    int                         forceSendOffset;            // Offset for double-buffered force messages
    int*                        pForceSendNode;             // List of nodes to send force data to
    int*                        pForceSendMinCell;          // Minimum atom cell containing force data to send
    int*                        pForceSendMaxCell;          // Maximum atom cell containing force data to send   
    int*                        pOutForceSendStart;         // Beginning of force data to send to each node
    int*                        pForceSendStart;            // Minimum atom containing force data to send
    int*                        pForceSendLength;           // Maximum atom containing force data to send        
    int                         forceReceiveFirstCell;      // First cell that receives force data
    int                         forceReceiveLastCell;       // Last cell that receives force data          
    int                         forceReceiveFirstAtom;      // First atom that receives force data
    int                         forceReceiveLastAtom;       // Last atom that receives force data              
    void*                       pForceData;                 // Incoming force data buffer;
    PMEDouble*                  pForceData0;                // Incoming force data buffer 0
    PMEDouble*                  pForceData1;                // Incoming force data buffer 1 
    MPI_Win                     MPIPMEForceWindow;          // MPI window for distributing reciprocal sum forces
    MPI_Group                   MPIPMEForceSendGroup;       // List of nodes that will receive force data from this node
    MPI_Group                   MPIPMEForceReceiveGroup;    // List of nodes that send force data to this node        
#endif

    // Methods
    _gpuContext();
    ~_gpuContext();
};

typedef struct _gpuContext *gpuContext;

template <typename T>
struct GpuBuffer
{
    unsigned int    _length;
    bool            _bSysMem;
    bool            _bPinned;
    T*              _pSysData;
    T*              _pDevData;
    GpuBuffer(int length, bool bSysMem = true, bool bPinned = false);
    GpuBuffer(unsigned int length, bool bSysMem = true, bool bPinned = false);
    virtual ~GpuBuffer();
    void Allocate();
    void Deallocate();
    void Upload(T* pBuff = NULL);
    void Download(T * pBuff = NULL);
};

#ifdef GPU_CPP
static gpuContext gpu = NULL;

template <typename T>
GpuBuffer<T>::GpuBuffer(unsigned int length, bool bSysMem, bool bPinned) : _length(length), _bSysMem(bSysMem), _bPinned(bPinned), _pSysData(NULL), _pDevData(NULL)
{
    Allocate();   
}

template <typename T>
GpuBuffer<T>::GpuBuffer(int length, bool bSysMem, bool bPinned) : _length(length), _bSysMem(bSysMem), _bPinned(bPinned), _pSysData(NULL), _pDevData(NULL)
{
    Allocate();   
}

template <typename T>
GpuBuffer<T>::~GpuBuffer()
{
    Deallocate();
}

template <typename T>
void GpuBuffer<T>::Allocate()
{
#ifdef MEMTRACKING
    printf("Allocating %d bytes of GPU memory", _length * sizeof(T));
    if (!_bSysMem && !_bPinned)
        printf(", unshadowed");
    else if (_bPinned)
        printf(", pinned");
    printf("\n");   
#endif
    cudaError_t status;
    if (_bPinned)
    {
        status = cudaHostAlloc((void **)&_pSysData, _length * sizeof(T), cudaHostAllocMapped);
        RTERROR(status, "cudaHostalloc GpuBuffer::Allocate failed");
        gpu->totalCPUMemory                +=  _length * sizeof(T);
        gpu->totalGPUMemory                +=  _length * sizeof(T);
        status = cudaHostGetDevicePointer((void **)&_pDevData, (void *)_pSysData, 0);
        RTERROR(status, "cudaGetDevicePointer GpuBuffer::failed to get device pointer");
        memset(_pSysData, 0, _length * sizeof(T));
    }
    else 
    {
        if (_bSysMem)
        {
            _pSysData =     new T[_length];
            gpu->totalCPUMemory            +=  _length * sizeof(T);
            memset(_pSysData, 0, _length * sizeof(T));
        }

        status = cudaMalloc((void **) &_pDevData, _length * sizeof(T));
        gpu->totalGPUMemory                +=  _length * sizeof(T);
        RTERROR(status, "cudaMalloc GpuBuffer::Allocate failed");
        status = cudaMemset((void *) _pDevData, 0, _length * sizeof(T));
        RTERROR(status, "cudaMemset GpuBuffer::Allocate failed");
    }
#ifdef MEMTRACKING
    printf("Mem++: %lld %lld\n", gpu->totalGPUMemory, gpu->totalCPUMemory);     
#endif
}

template <typename T>
void GpuBuffer<T>::Deallocate()
{
    cudaError_t status;
    if (_bPinned)
    {
        status = cudaFreeHost(_pSysData);
        gpu->totalCPUMemory                -=  _length * sizeof(T);
        gpu->totalGPUMemory                -=  _length * sizeof(T);        
    }
    else
    {
        if (_bSysMem)
        {
            delete[] _pSysData;
            gpu->totalCPUMemory            -=  _length * sizeof(T);   
        }         
        status = cudaFree(_pDevData);
        gpu->totalGPUMemory                -=  _length * sizeof(T);
    }
    RTERROR(status, "cudaFree GpuBuffer::Deallocate failed");   
    _pSysData = NULL;
    _pDevData = NULL;
#ifdef MEMTRACKING    
    printf("Mem--: %lld %lld\n", gpu->totalGPUMemory, gpu->totalCPUMemory);     
#endif
}
#endif

template <typename T>
void GpuBuffer<T>::Upload(T* pBuff)
{
    if (pBuff)
    {
        cudaError_t status;
        status = cudaMemcpy(_pDevData, pBuff, _length * sizeof(T), cudaMemcpyHostToDevice);
        RTERROR(status, "cudaMemcpy GpuBuffer::Upload failed");
    }
    else if (_bSysMem)
    {
        cudaError_t status;
        status = cudaMemcpy(_pDevData, _pSysData, _length * sizeof(T), cudaMemcpyHostToDevice);
        RTERROR(status, "cudaMemcpy GpuBuffer::Upload failed");
    }
}

template <typename T>
void GpuBuffer<T>::Download(T* pBuff)
{
    if (pBuff)
    {
        cudaError_t status;
        status = cudaMemcpy(pBuff, _pDevData, _length * sizeof(T), cudaMemcpyDeviceToHost);
        RTERROR(status, "cudaMemcpy GpuBuffer::Download failed");
    }
    else if (_bSysMem)
    {
        cudaError_t status;
        status = cudaMemcpy(_pSysData, _pDevData, _length * sizeof(T), cudaMemcpyDeviceToHost);
        RTERROR(status, "cudaMemcpy GpuBuffer::Download failed");
    }
}

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
#define printf(f,...)
#endif
#endif


