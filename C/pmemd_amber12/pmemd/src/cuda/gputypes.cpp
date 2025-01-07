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

#include "gputypes.h"
#include <math.h>

void clearCudaSimulation(cudaSimulation& sim)
{
    sim.atoms                   = 0;
    sim.grid                    = GRID;
    sim.gridBits                = GRIDBITS;
    sim.ntp                     = 0;
    sim.alpb                    = 0;
    sim.igb                     = 0;
    sim.scnb                    = (1.0 / 2.0);
    sim.scee                    = (1.0 / 1.2);
    sim.cut                     = 8.0;
    sim.cut2                    = sim.cut * sim.cut;
    sim.skinnb                  = 2.0f;
    sim.dielc                   = 1.0;
    sim.tol                     = 0.0001;
    sim.bUseVlimit              = false;
    sim.vlimit                  = 0.0;
    sim.massH                   = 1.008;
    sim.gb_alpha                = 0.0;
    sim.gb_beta                 = 0.0;
    sim.gb_gamma                = 0.0; 
    sim.gb_fs_max               = 2.0;
    sim.rgbmax                  = 25.0;
    sim.intdiel                 = 1.0;
    sim.extdiel                 = 78.5;
    sim.alpb_alpha              = 0.571412;
    sim.saltcon                 = 0.0;
    sim.surften                 = 0.005;
    sim.offset                  = 0.09;
    sim.gb_neckscale            = 0.361825;
    sim.gb_neckcut              = 2.8;
    sim.gb_neckoffset           = 1.0 - 0.09;
    sim.arad                    = 15.0;
    sim.a                       = 100.0;
    sim.b                       = 100.0;
    sim.c                       = 100.0;
    sim.alpha                   = 90.0;
    sim.beta                    = 90.0;
    sim.gamma                   = 90.0;
    sim.pi_vol_inv              = 0.0;
    sim.fac                     = 0.0;
    sim.fac2                    = 0.0;
    sim.is_orthog               = true;
    sim.paddedNumberOfAtoms     = 0;
    sim.pAtomX                  = NULL;
    sim.pAtomY                  = NULL;
    sim.pAtomXYSP               = NULL;
    sim.pAtomZ                  = NULL;
    sim.pAtomZSP                = NULL;
    sim.pAtomSigEps             = NULL;
    sim.pAtomRBorn              = NULL;
    sim.pAtomS                  = NULL;
    sim.pgb_alpha               = NULL;
    sim.pgb_beta                = NULL;
    sim.pgb_gamma               = NULL;
    sim.pAtomInvMass            = NULL;
    sim.pAtomMass               = NULL;
    sim.pAtomCharge             = NULL;
    sim.pAtomChargeSP           = NULL;
    sim.pReff                   = NULL;
    sim.pReffSP                 = NULL;
    sim.pPsi                    = NULL;
    sim.pTemp7                  = NULL;
    sim.pForce                  = NULL;
    sim.pForceX                 = NULL;
    sim.pForceY                 = NULL;
    sim.pForceZ                 = NULL;
    sim.pBondedForceAccumulator = NULL;
    sim.pBondedForceXAccumulator= NULL;
    sim.pBondedForceYAccumulator= NULL;
    sim.pBondedForceZAccumulator= NULL;
    sim.pNBForce                = NULL;
    sim.pNBForceX               = NULL;
    sim.pNBForceY               = NULL;
    sim.pNBForceZ               = NULL;   
#ifdef MPI    
    sim.pInForce                = NULL;
    sim.pOutForce               = NULL;
    sim.pPMEForce               = NULL;  
    sim.pTemp7a                 = NULL;
    sim.pReffa                  = NULL;   
#endif
    sim.pVelX                   = NULL;
    sim.pVelY                   = NULL;
    sim.pVelZ                   = NULL;
    sim.pLVelX                  = NULL;
    sim.pLVelY                  = NULL;
    sim.pLVelZ                  = NULL;
    sim.randomSteps             = 0;
    sim.randomNumbers           = 0;
    sim.pRandomPos              = NULL;
    sim.pRandom                 = NULL;
    sim.pRandomX                = NULL;
    sim.pRandomY                = NULL;
    sim.pRandomZ                = NULL;
    sim.shakeConstraints        = 0;
    sim.shakeOffset             = 0;
    sim.fastShakeConstraints    = 0;
    sim.fastShakeOffset         = 0;
    sim.slowShakeConstraints    = 0;
    sim.slowShakeOffset         = 0;
    sim.pShakeID                = NULL;
    sim.pShakeParm              = NULL;
    sim.pFastShakeID            = NULL;
    sim.pSlowShakeID1           = NULL;
    sim.pSlowShakeID2           = NULL;
    sim.pSlowShakeParm          = NULL;
    sim.pImageShakeID           = NULL;
    sim.pImageFastShakeID       = NULL;
    sim.pImageSlowShakeID1      = NULL;
    sim.pImageSlowShakeID2      = NULL;
    sim.bonds                   = 0;
    sim.bondAngles              = 0;
    sim.bondAngleOffset         = 0;
    sim.dihedrals               = 0;
    sim.dihedralOffset          = 0;   
    sim.nb14s                   = 0;
    sim.nb14Offset              = 0;
    sim.constraints             = 0;   
    sim.constraintOffset        = 0;
    sim.UBAngles                = 0;
    sim.UBAngleOffset           = 0;
    sim.impDihedrals            = 0;
    sim.impDihedralOffset       = 0;
    sim.cmaps                   = 0;
    sim.cmapOffset              = 0;
    sim.pBond                   = NULL;
    sim.pBondID                 = NULL;
    sim.pBondAngle              = NULL;
    sim.pBondAngleID1           = NULL;
    sim.pBondAngleID2           = NULL;
    sim.pDihedral1              = NULL;
    sim.pDihedral2              = NULL;
    sim.pDihedral3              = NULL;
    sim.pDihedralID1            = NULL;
    sim.pDihedralID2            = NULL;
    sim.pNb141                  = NULL;
    sim.pNb142                  = NULL;
    sim.pNb14ID                 = NULL;
    sim.pConstraint1            = NULL;
    sim.pConstraint2            = NULL;
    sim.pConstraintID           = NULL;
    sim.pUBAngle                = NULL;
    sim.pUBAngleID              = NULL;
    sim.pImpDihedral            = NULL;    
    sim.pImpDihedralID1         = NULL; 
    sim.pImpDihedralID2         = NULL;
    sim.pCmapID1                = NULL;
    sim.pCmapID2                = NULL;
    sim.pCmapID3                = NULL;
    sim.pCmapType               = NULL;
    sim.pCmapEnergy             = NULL;


    sim.NMRDistances            = 0;
    sim.NMRDistanceOffset       = 0;
    sim.NMRAngles               = 0;
    sim.NMRAngleOffset          = 0;
    sim.NMRTorsions             = 0;
    sim.NMRTorsionOffset        = 0;
    sim.bJar                    = false;
    sim.drjar                   = 0.0; 
    sim.pNMRJarData             = NULL;    
    sim.pNMRDistanceID          = NULL;
    sim.pNMRDistanceR1R2        = NULL;
    sim.pNMRDistanceR3R4        = NULL;
    sim.pNMRDistanceK2K3        = NULL;
    sim.pNMRDistanceAve         = NULL;
    sim.pNMRDistanceTgtVal      = NULL;
    sim.pNMRDistanceStep        = NULL; 
    sim.pNMRDistanceInc         = NULL; 
    sim.pNMRDistanceR1R2Slp     = NULL;
    sim.pNMRDistanceR3R4Slp     = NULL;
    sim.pNMRDistanceK2K3Slp     = NULL;
    sim.pNMRDistanceR1R2Int     = NULL;
    sim.pNMRDistanceR3R4Int     = NULL;
    sim.pNMRDistanceK2K3Int     = NULL;
    sim.pNMRAngleID1            = NULL;
    sim.pNMRAngleID2            = NULL;
    sim.pNMRAngleR1R2           = NULL;
    sim.pNMRAngleR3R4           = NULL;
    sim.pNMRAngleK2K3           = NULL;
    sim.pNMRAngleAve            = NULL;
    sim.pNMRAngleTgtVal         = NULL;
    sim.pNMRAngleStep           = NULL;
    sim.pNMRAngleInc            = NULL;
    sim.pNMRAngleR1R2Slp        = NULL;
    sim.pNMRAngleR3R4Slp        = NULL;
    sim.pNMRAngleK2K3Slp        = NULL;
    sim.pNMRAngleR1R2Int        = NULL; 
    sim.pNMRAngleR3R4Int        = NULL; 
    sim.pNMRAngleK2K3Int        = NULL;  
    sim.pNMRTorsionID1          = NULL;
    sim.pNMRTorsionID2          = NULL;
    sim.pNMRTorsionR1R2         = NULL;
    sim.pNMRTorsionR3R4         = NULL;
    sim.pNMRTorsionK2K3         = NULL;
    sim.pNMRTorsionAve1         = NULL;
    sim.pNMRTorsionAve2         = NULL;
    sim.pNMRAngleTgtVal         = NULL;
    sim.pNMRTorsionStep         = NULL;
    sim.pNMRTorsionInc          = NULL;
    sim.pNMRTorsionR1R2Slp      = NULL;
    sim.pNMRTorsionR3R4Slp      = NULL;
    sim.pNMRTorsionK2K3Slp      = NULL;
    sim.pNMRTorsionR1R2Int      = NULL;
    sim.pNMRTorsionR3R4Int      = NULL;
    sim.pNMRTorsionK2K3Int      = NULL; 
    sim.EPs                     = 0;
    sim.EP11s                   = 0;
    sim.EP12s                   = 0;
    sim.EP21s                   = 0;
    sim.EP22s                   = 0;
    sim.EP11Offset              = 0;
    sim.EP12Offset              = 0;
    sim.EP21Offset              = 0;
    sim.EP22Offset              = 0;
    sim.pExtraPoint11Frame      = NULL;
    sim.pExtraPoint11Index      = NULL;
    sim.pExtraPoint11X          = NULL;
    sim.pExtraPoint11Y          = NULL;
    sim.pExtraPoint11Z          = NULL;
    sim.pExtraPoint12Frame      = NULL;
    sim.pExtraPoint12Index      = NULL;
    sim.pExtraPoint12X          = NULL;
    sim.pExtraPoint12Y          = NULL;
    sim.pExtraPoint12Z          = NULL;
    sim.pExtraPoint21Frame      = NULL;
    sim.pExtraPoint21Index      = NULL;
    sim.pExtraPoint21X1         = NULL;
    sim.pExtraPoint21Y1         = NULL;
    sim.pExtraPoint21Z1         = NULL;    
    sim.pExtraPoint21X2         = NULL;
    sim.pExtraPoint21Y2         = NULL;
    sim.pExtraPoint21Z2         = NULL;        
    sim.pExtraPoint22Frame      = NULL;
    sim.pExtraPoint22Index      = NULL;
    sim.pExtraPoint22X1         = NULL;
    sim.pExtraPoint22Y1         = NULL;
    sim.pExtraPoint22Z1         = NULL;    
    sim.pExtraPoint22X2         = NULL;
    sim.pExtraPoint22Y2         = NULL;
    sim.pExtraPoint22Z2         = NULL;     
    sim.pImageExtraPoint11Frame = NULL;
    sim.pImageExtraPoint11Index = NULL;
    sim.pImageExtraPoint12Frame = NULL;
    sim.pImageExtraPoint12Index = NULL;
    sim.pImageExtraPoint21Frame = NULL;
    sim.pImageExtraPoint21Index = NULL;    
    sim.pImageExtraPoint22Frame = NULL;
    sim.pImageExtraPoint22Index = NULL;
#ifdef use_SPFP
    sim.pForceAccumulator       = NULL;
    sim.pForceXAccumulator      = NULL;
    sim.pForceYAccumulator      = NULL;
    sim.pForceZAccumulator      = NULL;
    sim.pNBForceAccumulator     = NULL;
    sim.pNBForceXAccumulator    = NULL;
    sim.pNBForceYAccumulator    = NULL;
    sim.pNBForceZAccumulator    = NULL;          
    sim.pReffAccumulator        = NULL;
    sim.pSumdeijdaAccumulator   = NULL;  
#else         
    sim.pForceBuffer            = NULL;
    sim.pForceXBuffer           = NULL;
    sim.pForceYBuffer           = NULL;
    sim.pForceZBuffer           = NULL;          
    sim.pReffBuffer             = NULL;
    sim.pSumdeijdaBuffer        = NULL;
#endif
    sim.pEnergyBuffer           = NULL; 
    sim.pNeckMaxValPos          = NULL;
    sim.stride                  = 0;
    sim.stride2                 = 0;
    sim.stride3                 = 0;
    sim.stride4                 = 0;
    sim.imageStride             = 0;     
    sim.ew_coeff                = (PMEFloat)1.0;
    sim.ew_coeff2               = sim.ew_coeff * sim.ew_coeff;
    sim.negTwoEw_coeffRsqrtPI   = (PMEFloat)(-2.0 * sim.ew_coeff / sqrt(PI));
    sim.nfft1                   = 0;
    sim.nfft2                   = 0;
    sim.nfft3                   = 0;
    sim.fft_x_dim               = 0;
    sim.fft_y_dim               = 0;
    sim.fft_z_dim               = 0;
    sim.order                   = 4;
    sim.pPrefac1                = NULL;
    sim.pPrefac2                = NULL;
    sim.pPrefac3                = NULL;
    sim.pIFractX                = NULL;
    sim.pIFractY                = NULL;
    sim.pIFractZ                = NULL;
    sim.pThetaX                 = NULL;
    sim.pThetaY                 = NULL;
    sim.pThetaZ                 = NULL;
    sim.pDThetaX                = NULL;
    sim.pDThetaY                = NULL;
    sim.pDThetaZ                = NULL;
#ifdef use_SPFP
    sim.plliXYZ_q               = NULL;
#endif
    sim.pXYZ_q                  = NULL;
    sim.pXYZ_qt                 = NULL;
    
    sim.cells                   = 0;
    sim.cell                    = (PMEFloat)1.0;
    sim.xcell                   = (PMEDouble)1.0;
    sim.ycell                   = (PMEDouble)1.0;
    sim.zcell                   = (PMEDouble)1.0;
    sim.minCellX                = (PMEDouble)-2.0;
    sim.minCellY                = (PMEDouble)-2.0;
    sim.minCellZ                = (PMEDouble)-2.0;
    sim.maxCellX                = (PMEDouble)2.0;
    sim.maxCellY                = (PMEDouble)2.0;
    sim.maxCellZ                = (PMEDouble)2.0;    
    sim.oneOverXcells           = (PMEFloat)1.0;
    sim.oneOverYcells           = (PMEFloat)1.0;
    sim.oneOverZcells           = (PMEFloat)1.0;    
    sim.xcells                  = 0;
    sim.ycells                  = 0;
    sim.zcells                  = 0;
    sim.pNLbSkinTestFail        = NULL;
    sim.pNLCellHash             = NULL;
    sim.pImageIndex             = NULL;
    sim.pImageIndex2            = NULL;
    sim.pImageAtom              = NULL;
    sim.pImageAtom2             = NULL;
    sim.pImageAtomLookup        = NULL;
    sim.pImageHash              = NULL;
    sim.pAtomXYSaveSP           = NULL;
    sim.pAtomZSaveSP            = NULL;
    sim.pImageX                 = NULL;
    sim.pImageY                 = NULL;
    sim.pImageZ                 = NULL; 
    sim.pImageVelX              = NULL;
    sim.pImageVelY              = NULL;
    sim.pImageVelZ              = NULL;
    sim.pImageLVelX             = NULL;
    sim.pImageLVelY             = NULL;
    sim.pImageLVelZ             = NULL;
    sim.pImageMass              = NULL;
    sim.pImageCharge            = NULL;
    sim.pImageSigEps            = NULL; 
    sim.pImageOutputBuffers     = NULL;
    sim.pImageCellID            = NULL;    
    sim.pImageX2                = NULL;
    sim.pImageY2                = NULL;
    sim.pImageZ2                = NULL; 
    sim.pImageVelX2             = NULL;
    sim.pImageVelY2             = NULL;
    sim.pImageVelZ2             = NULL;
    sim.pImageLVelX2            = NULL;
    sim.pImageLVelY2            = NULL;
    sim.pImageLVelZ2            = NULL;
    sim.pImageMass2             = NULL;
    sim.pImageCharge2           = NULL;
    sim.pImageSigEps2           = NULL; 
    sim.pImageOutputBuffers2    = NULL;
    sim.pImageCellID2           = NULL;
    sim.pImageBondID            = NULL;
    sim.pImageBondAngleID1      = NULL;
    sim.pImageBondAngleID2      = NULL;
    sim.pImageDihedralID1       = NULL;
    sim.pImageDihedralID2       = NULL;
    sim.pImageNb14ID            = NULL;
    sim.pImageConstraintID      = NULL;
    sim.pImageUBAngleID         = NULL;
    sim.pImageImpDihedralID1    = NULL;
    sim.pImageImpDihedralID2    = NULL;
    sim.pImageCmapID1           = NULL;
    sim.pImageCmapID2           = NULL;
    sim.pImageCmapID3           = NULL;
    sim.pImageNMRDistanceID     = NULL;
    sim.pImageNMRAngleID1       = NULL;
    sim.pImageNMRAngleID2       = NULL;
    sim.pImageNMRTorsionID1     = NULL;
    sim.pImageNMRTorsionID2     = NULL; 
    sim.pImageSoluteAtomID      = NULL;
    sim.pImageSolventAtomID     = NULL;
    sim.pBNLExclusionBuffer     = NULL;
    sim.pNLNonbondCellStartEnd  = NULL;
    sim.pNLExclusionStartCount  = NULL;
    sim.pNLExclusionList        = NULL; 
    sim.pNLAtomList             = NULL;  
    sim.pNLTotalOffset          = NULL; 
    sim.maxNonbonds             = 0;
    sim.maxNonbondBuffers       = 0;
    sim.soluteMolecules         = 0;
    sim.soluteAtoms             = 0;
    sim.solventMolecules        = 0;
    sim.solventMoleculeStride   = 0;
    sim.pSoluteAtomMoleculeID   = NULL;
    sim.pSoluteAtomID           = NULL;
    sim.pSoluteAtomMass         = NULL;
    sim.pSoluteCOMX             = NULL;
    sim.pSoluteCOMY             = NULL;
    sim.pSoluteCOMZ             = NULL;
    sim.pSoluteDeltaCOMX        = NULL;
    sim.pSoluteDeltaCOMY        = NULL;
    sim.pSoluteDeltaCOMZ        = NULL;  
    sim.pSoluteUllCOMX          = NULL;
    sim.pSoluteUllCOMY          = NULL;
    sim.pSoluteUllCOMZ          = NULL;      
    sim.pSoluteUllEKCOMX        = NULL;
    sim.pSoluteUllEKCOMY        = NULL;
    sim.pSoluteUllEKCOMZ        = NULL;          
    sim.pSoluteInvMass          = NULL;
    sim.pSolventAtomID          = NULL;
    sim.pSolventAtomMass1       = NULL;
    sim.pSolventAtomMass2       = NULL;
    sim.pSolventAtomMass3       = NULL;
    sim.pSolventAtomMass4       = NULL;
    sim.pSolventCOMX            = NULL;
    sim.pSolventCOMY            = NULL;
    sim.pSolventCOMZ            = NULL;
    sim.pSolventInvMass         = NULL;
    sim.pKineticEnergy          = NULL; 
    sim.pConstraintAtomX        = NULL;
    sim.pConstraintAtomY        = NULL;
    sim.pConstraintAtomZ        = NULL;
    sim.pConstraintCOMX         = NULL;
    sim.pConstraintCOMY         = NULL;
    sim.pConstraintCOMZ         = NULL;

    // General IPS parameters
    sim.bIPSActive              = false;
    sim.rips                    = (PMEFloat)8.0;
    sim.eipssnb                 = 0.0;
    sim.eipssel                 = 0.0;
    sim.virips                  = 0.0;
    
    //  Electrostatic IPS parameters:
    sim.aipse0                  = -35.0 / 16.0;
    sim.aipse1                  =  35.0 / 16.0;
    sim.aipse2                  = -21.0 / 16.0;
    sim.aipse3                  =   5.0 / 16.0;
    sim.pipsec                  =   1.0 + sim.aipse0 + sim.aipse1 + sim.aipse2 + sim.aipse3;
    sim.pipse0                  =   sim.aipse0 - sim.pipsec; 
    sim.bipse1                  =   2.0 * sim.aipse1;
    sim.bipse2                  =   4.0 * sim.aipse2;
    sim.bipse3                  =   6.0 * sim.aipse3;
        
    //  Dispersion IPS parameters:
    sim.aipsvc0                 =   7.0 / 16.0;
    sim.aipsvc1                 =   9.0 / 14.0;
    sim.aipsvc2                 =  -3.0 / 28.0;
    sim.aipsvc3                 =   6.0 / 7.0;
    sim.pipsvcc                 =   1.0 + sim.aipsvc0 + sim.aipsvc1 + sim.aipsvc2 + sim.aipsvc3;
    sim.pipsvc0                 =   sim.aipsvc0 - sim.pipsvcc;
    sim.bipsvc1                 =   2.0 * sim.aipsvc1;
    sim.bipsvc2                 =   4.0 * sim.aipsvc2;  
    sim.bipsvc3                 =   6.0 * sim.aipsvc3; 

    //  Repulsion IPS parameters:
    sim.aipsva0                 =  5.0 / 787.0;
    sim.aipsva1                 =  9.0 /  26.0;
    sim.aipsva2                 = -3.0 /  13.0;
    sim.aipsva3                 = 27.0 /  26.0;
    sim.pipsvac                 =  1.0 + sim.aipsva0 + sim.aipsva1 + sim.aipsva2 + sim.aipsva3;
    sim.pipsva0                 =  sim.aipsva0 - sim.pipsvac;
    sim.bipsva1                 =  4.0 * sim.aipsva1;
    sim.bipsva2                 =  8.0 * sim.aipsva2;
    sim.bipsva3                 = 12.0 * sim.aipsva3;

    // Derived parameters
    sim.rips2                   = sim.rips * sim.rips;
    sim.ripsr                   = (PMEFloat)1.0 / (sim.rips);
    sim.rips2r                  = (PMEFloat)1.0 / sim.rips2;
    sim.rips6r                  = sim.rips2r * sim.rips2r * sim.rips2r;
    sim.rips12r                 = sim.rips6r * sim.rips6r;              
    sim.rgbmax1i                = 1.0 / sim.rgbmax;
    sim.rgbmax2i                = sim.rgbmax1i * sim.rgbmax1i;
    sim.rgbmaxpsmax2            = (sim.rgbmax + sim.gb_fs_max) * (sim.rgbmax + sim.gb_fs_max);
    sim.alpb_beta               = sim.alpb_alpha * (sim.intdiel / sim.extdiel);
    sim.one_arad_beta           = sim.alpb_beta / sim.arad;
    if (sim.saltcon >= 0) 
    {
        sim.gb_kappa            = 0.73 * sqrt(0.10806 * sim.saltcon);
        if (sim.gb_kappa > 0.0)
            sim.gb_kappa_inv    = 1.0 / sim.gb_kappa;
    }
    else
    {
        sim.gb_kappa            = 0.0;
        sim.gb_kappa_inv        = 0.0;
    }
    if (sim.alpb == 0)
    {
        sim.intdiel_inv         = 1.0 / sim.intdiel;
        sim.extdiel_inv         = 1.0 / sim.extdiel;
    }
    else
    {
        sim.extdiel_inv         = 1.0 / (sim.extdiel * (1.0 + sim.alpb_beta));
        sim.intdiel_inv         = 1.0 / (sim.intdiel * (1.0 + sim.alpb_beta));
    }
    sim.invMassH                = 1.0 / sim.massH;
    
    sim.ucell[0][0]             = sim.a;
    sim.ucell[0][1]             = 0.0;
    sim.ucell[0][2]             = 0.0;
    sim.ucell[1][0]             = 0.0;
    sim.ucell[1][1]             = sim.b;
    sim.ucell[1][2]             = 0.0;
    sim.ucell[2][0]             = 0.0;
    sim.ucell[2][1]             = 0.0;
    sim.ucell[2][2]             = sim.c;
    sim.recip[0][0]             = 1.0 / sim.a;
    sim.recip[0][1]             = 0.0;
    sim.recip[0][2]             = 0.0;
    sim.recip[1][0]             = 0.0;
    sim.recip[1][1]             = 1.0 / sim.b;
    sim.recip[1][2]             = 0.0;
    sim.recip[2][0]             = 0.0;
    sim.recip[2][1]             = 0.0;
    sim.recip[2][2]             = 1.0 / sim.c;
    sim.pbc_box[0]              = sim.a;
    sim.pbc_box[1]              = sim.b;
    sim.pbc_box[2]              = sim.c;
    sim.reclng[0]               = 1.0 / sim.a;
    sim.reclng[1]               = 1.0 / sim.b;
    sim.reclng[2]               = 1.0 / sim.c;
    sim.uc_volume               = sim.a * sim.b * sim.c;
    sim.uc_sphere               = sqrt(sim.a * sim.a + sim.b * sim.b + sim.c * sim.c);
    sim.NLMaxTotalOffset        = 0;   
    sim.NLWorkUnits             = 0;   
     
    // Set up AMD parameters
    sim.iamd                    = 0;
    sim.iamdlag                 = 1;
    sim.amd_print_interval      = 0;
    sim.amd_EthreshP            = 0.0;
    sim.amd_alphaP              = 0.0;
    sim.amd_EthreshD            = 0.0;
    sim.amd_alphaD              = 0.0;
    sim.amd_temp0               = 0.0;
    sim.AMDNumLag               = 0;
    sim.AMDNumRecs              = -1;
    sim.AMDtboost               = 0.0;
    sim.pAMDfwgtd               = NULL;
    sim.AMDfwgt                 = 1.0;
    sim.pAMDEDihedral           = NULL;

    // Set up GBSA parameters
    sim.igbsa                    = 0;
}

_gpuContext::_gpuContext() :
bECCSupport(false),
totalCPUMemory(0),
totalGPUMemory(0),
#ifdef MPI
nGpus(1),
gpuID(0),
#endif
sm_version(SM_10),
ntb(1),
ips(0),
ntc(1),
ntf(1),
ntr(0),
forwardPlan(0),
backwardPlan(0),
pbAtom(NULL),
pbAtomXYSP(NULL),
pbAtomZSP(NULL),
pbAtomSigEps(NULL),
pbAtomCharge(NULL),
pbAtomChargeSP(NULL),
pbAtomRBorn(NULL),
pbAtomS(NULL),
pbAtomMass(NULL),
pbCenter(NULL),
pbForce(NULL),
#ifndef use_SPFP
pbBondedForce(NULL),
#endif
#ifdef MPI
pbReffa(NULL),
pbTemp7a(NULL),
pbInForce(NULL),
pbOutForce(NULL),
pbPMEForce(NULL),
pPMEStart(NULL),
pPMELength(NULL),
pMinLocalCell(NULL),
pMaxLocalCell(NULL),
pMinLocalAtom(NULL),
pMaxLocalAtom(NULL),
pMinProcessedCell(NULL),
pMaxProcessedCell(NULL),
pMinProcessedAtom(NULL),
pMaxProcessedAtom(NULL),
pAllGathervRecvCountAoS(NULL),
pAllGathervRecvDisplAoS(NULL),
pAllGathervRecvCountSoA(NULL),
pAllGathervRecvDisplSoA(NULL),
pForceSendNode(NULL),
pForceSendMinCell(NULL),
pForceSendMaxCell(NULL),
pOutForceSendStart(NULL),
pForceSendStart(NULL),
pForceSendLength(NULL),     
pPMEForceData(NULL),
pForceData(NULL),
bSingleNode(false),
#ifdef CUDA_P2P
bP2P(false),
#endif
pSharedForceBuffer0(NULL),
pSharedForceBuffer1(NULL),
pSharedForceBuffer2(NULL),
pSharedPMEForceBuffer0(NULL),
pSharedPMEForceBuffer1(NULL),
pSharedPMEForceBuffer2(NULL),    
#ifdef CUDA_P2P
pP2PReffa(NULL),
pP2PTemp7a(NULL), 
pP2PInForce(NULL), 
pP2PReffaHandle(NULL),
pP2PTemp7aHandle(NULL), 
pP2PInForceHandle(NULL), 
#endif
#endif
pbVel(NULL),
pbLVel(NULL),
bLocalInteractions(true),
bCharmmInteractions(false),
bNMRInteractions(false),
NMRnstep(0),
pbBond(NULL),
pbBondID(NULL),
pbBondAngle(NULL),
pbBondAngleID1(NULL),
pbBondAngleID2(NULL),
pbDihedral1(NULL),
pbDihedral2(NULL),
pbDihedral3(NULL),
pbDihedralID1(NULL),
pbDihedralID2(NULL),
pbNb141(NULL),
pbNb142(NULL),
pbNb14ID(NULL),
pbConstraint1(NULL),
pbConstraint2(NULL),
pbConstraintID(NULL),
pbUBAngle(NULL),
pbUBAngleID(NULL),
pbImpDihedral(NULL),
pbImpDihedralID1(NULL), 
pbImpDihedralID2(NULL),
pbCmapID1(NULL),
pbCmapID2(NULL),
pbCmapID3(NULL),
pbCmapType(NULL),
pbCmapEnergy(NULL),
pbNMRJarData(NULL),
pbNMRDistanceID(NULL),
pbNMRDistanceR1R2(NULL),
pbNMRDistanceR3R4(NULL),
pbNMRDistanceK2K3(NULL),
pbNMRDistanceK4(NULL),
pbNMRDistanceAve(NULL),
pbNMRDistanceTgtVal(NULL),
pbNMRDistanceStep(NULL), 
pbNMRDistanceInc(NULL),
pbNMRDistanceR1R2Slp(NULL),
pbNMRDistanceR3R4Slp(NULL),
pbNMRDistanceK2K3Slp(NULL),
pbNMRDistanceK4Slp(NULL),
pbNMRDistanceR1R2Int(NULL),
pbNMRDistanceR3R4Int(NULL),
pbNMRDistanceK2K3Int(NULL),
pbNMRDistanceK4Int(NULL),
pbNMRAngleID1(NULL),
pbNMRAngleID2(NULL),
pbNMRAngleR1R2(NULL),
pbNMRAngleR3R4(NULL),
pbNMRAngleK2K3(NULL),
pbNMRAngleK4(NULL),
pbNMRAngleAve(NULL),
pbNMRAngleTgtVal(NULL),
pbNMRAngleStep(NULL),
pbNMRAngleInc(NULL),
pbNMRAngleR1R2Slp(NULL),
pbNMRAngleR3R4Slp(NULL),
pbNMRAngleK2K3Slp(NULL),
pbNMRAngleK4Slp(NULL),
pbNMRAngleR1R2Int(NULL), 
pbNMRAngleR3R4Int(NULL), 
pbNMRAngleK2K3Int(NULL), 
pbNMRAngleK4Int(NULL), 
pbNMRTorsionID1(NULL),
pbNMRTorsionID2(NULL),
pbNMRTorsionR1R2(NULL),
pbNMRTorsionR3R4(NULL),
pbNMRTorsionK2K3(NULL),
pbNMRTorsionK4(NULL),
pbNMRTorsionAve1(NULL),
pbNMRTorsionAve2(NULL),
pbNMRTorsionTgtVal(NULL),
pbNMRTorsionStep(NULL),
pbNMRTorsionInc(NULL),
pbNMRTorsionR1R2Slp(NULL),
pbNMRTorsionR3R4Slp(NULL),
pbNMRTorsionK2K3Slp(NULL),
pbNMRTorsionK4Slp(NULL),
pbNMRTorsionR1R2Int(NULL),
pbNMRTorsionR3R4Int(NULL),
pbNMRTorsionK2K3Int(NULL),     
pbNMRTorsionK4Int(NULL),     
pbOutputBufferCounter(NULL),
pbReff(NULL),
pbReffSP(NULL),
pbPsi(NULL),
pbTemp7(NULL),
#ifdef use_SPFP
pbReffAccumulator(NULL),
pbSumdeijdaAccumulator(NULL),
pbForceAccumulator(NULL),
#else
pbReffBuffer(NULL),
pbSumdeijdaBuffer(NULL),
pbForceBuffer(NULL),
#endif
pbEnergyBuffer(NULL),
pbKineticEnergyBuffer(NULL),
pbWorkUnit(NULL),
pbExclusion(NULL),
pbGBPosition(NULL),
pbNeckMaxValPos(NULL),
pbGBAlphaBetaGamma(NULL),
randomCounter(0),
pbRandomPos(NULL),
pbRandom(NULL),
pbShakeID(NULL),
pbShakeParm(NULL),
pbFastShakeID(NULL),
pbSlowShakeID1(NULL),
pbSlowShakeID2(NULL),
pbSlowShakeParm(NULL),
pbExtraPoint11Frame(NULL),
pbExtraPoint11Index(NULL),
pbExtraPoint11(NULL),
pbExtraPoint12Frame(NULL),
pbExtraPoint12Index(NULL),
pbExtraPoint12(NULL),
pbExtraPoint21Frame(NULL),
pbExtraPoint21Index(NULL),
pbExtraPoint21(NULL),
pbExtraPoint22Frame(NULL),
pbExtraPoint22Index(NULL),
pbExtraPoint22(NULL),
pbImageExtraPoint11Frame(NULL),
pbImageExtraPoint11Index(NULL),
pbImageExtraPoint12Frame(NULL),
pbImageExtraPoint12Index(NULL),
pbImageExtraPoint21Frame(NULL),
pbImageExtraPoint21Index(NULL),
pbImageExtraPoint22Frame(NULL),
pbImageExtraPoint22Index(NULL),
pbPrefac(NULL),
pbIFract(NULL),
pbTheta(NULL),
#ifdef use_SPFP
pblliXYZ_q(NULL),
#endif
pbXYZ_q(NULL),
pbXYZ_qt(NULL),
bNeighborList(false),
bNeedNewNeighborList(true),
bNewNeighborList(false),
bSmallBox(false),
bOddNLCells(false),
neighborListBits(32),
pbAtomXYSaveSP(NULL),
pbAtomZSaveSP(NULL),
pbImage(NULL),
pbImageIndex(NULL),
pbImageVel(NULL),
pbImageLVel(NULL),
pbImageMass(NULL),
pbImageCharge(NULL),
pbImageSigEps(NULL),
pbImageOutputBuffers(NULL),
pbImageCellID(NULL),
pbImageBondID(NULL),
pbImageBondAngleID1(NULL),      
pbImageBondAngleID2(NULL),      
pbImageDihedralID1(NULL),
pbImageDihedralID2(NULL),
pbImageNb14ID(NULL),  
pbImageConstraintID(NULL),
pbImageUBAngleID(NULL),
pbImageImpDihedralID1(NULL), 
pbImageImpDihedralID2(NULL),
pbImageCmapID1(NULL),
pbImageCmapID2(NULL),
pbImageCmapID3(NULL),
pbImageNMRDistanceID(NULL),
pbImageNMRAngleID1(NULL),      
pbImageNMRAngleID2(NULL),      
pbImageNMRTorsionID1(NULL),
pbImageNMRTorsionID2(NULL),
pbImageShakeID(NULL),
pbImageFastShakeID(NULL),
pbImageSlowShakeID1(NULL),
pbImageSlowShakeID2(NULL),
pbImageSolventAtomID(NULL),
pbImageSoluteAtomID(NULL),
pbBNLExclusionBuffer(NULL),
pbNLExclusionList(NULL),
pbNLExclusionStartCount(NULL),
pbNLAtomList(NULL),
pbNLOffset(NULL),
pbNLTotalOffset(NULL),
pbNLPosition(NULL),
pbNLNonbondCellStartEnd(NULL),
#ifndef use_SPFP
pbNLChargeGridBufferOffset(NULL),
pbNLOddBufferOverlapFlag(NULL),
#endif
pbNLbSkinTestFail(NULL),
pbNLCellHash(NULL),
maxSoluteMolecules(0),
maxPSSoluteMolecules(0),
pbSoluteAtomID(NULL),
pbSoluteAtomMass(NULL),
pbSolute(NULL),
pbUllSolute(NULL),
pbSolventAtomID(NULL),
pbSolvent(NULL),
pbNTPData(NULL),
pbConstraintAtomX(NULL),
pbConstraintAtomY(NULL),
pbConstraintAtomZ(NULL),
pbConstraintCOMX(NULL),
pbConstraintCOMY(NULL),
pbConstraintCOMZ(NULL),
ee_plasma(0.0),
self_energy(0.0),
vdw_recip(0.0),
// AMD buffers
pbAMDfwgtd(NULL),
pAmdWeightsAndEnergy(NULL)
{
   clearCudaSimulation(sim);
}

_gpuContext::~_gpuContext()
{
    // Delete Atom data
    delete pbAtom;
    delete pbAtomXYSP;
    delete pbAtomZSP;    
    delete pbAtomSigEps;
    delete pbAtomRBorn;
    delete pbAtomS;
    delete pbAtomCharge;
    delete pbAtomChargeSP;
    delete pbAtomMass;
    delete pbReff;
    delete pbReffSP;
    delete pbPsi;
    delete pbTemp7;
    delete pbForce;
#ifndef use_SPFP
    delete pbBondedForce;
#endif
#ifdef MPI
    delete[] pMinLocalCell;
    delete[] pMaxLocalCell;
    delete[] pMinProcessedCell;
    delete[] pMaxProcessedCell;
    delete[] pMinProcessedAtom;
    delete[] pMaxProcessedAtom;
    delete[] pAllGathervRecvCountAoS;
    delete[] pAllGathervRecvDisplAoS;
    delete[] pAllGathervRecvCountSoA;
    delete[] pAllGathervRecvDisplSoA;
    delete pbInForce;
    delete pbOutForce;
    delete pbPMEForce;
    delete pbTemp7a;
    delete pbReffa;
#ifdef CUDA_P2P
    delete[] pP2PReffa;
    delete[] pP2PTemp7a; 
    delete[] pP2PInForce; 
    delete[] pP2PReffaHandle;
    delete[] pP2PTemp7aHandle; 
    delete[] pP2PInForceHandle;
#endif
    delete[] pPMEStart;
    delete[] pPMELength;
    delete[] pForceSendNode;
    delete[] pForceSendMinCell;
    delete[] pForceSendMaxCell;
    delete[] pOutForceSendStart;  
    delete[] pForceSendStart;  
    delete[] pForceSendLength;      
#endif
    delete pbVel;
    delete pbLVel;
    delete pbCenter;
    delete pbOutputBufferCounter;
    
    // Delete PME stuff
    delete pbPrefac;
    delete pbIFract;
    delete pbTheta;
#ifdef use_SPFP
    delete pblliXYZ_q;
#endif
    delete pbXYZ_q;
    delete pbXYZ_qt;   
    cufftDestroy(forwardPlan);
    cufftDestroy(backwardPlan);
    
    // Delete neighbor list stuff
    delete pbAtomXYSaveSP;
    delete pbAtomZSaveSP;
    delete pbImage;
    delete pbImageIndex;
    delete pbImageVel;
    delete pbImageLVel;
    delete pbImageMass;
    delete pbImageCharge;
    delete pbImageSigEps;
    delete pbImageOutputBuffers;
    delete pbImageCellID;
    delete pbImageBondID;
    delete pbImageBondAngleID1;
    delete pbImageBondAngleID2;
    delete pbImageDihedralID1;
    delete pbImageDihedralID2;
    delete pbImageNb14ID;                
    delete pbImageConstraintID;
    delete pbImageUBAngleID;
    delete pbImageImpDihedralID1; 
    delete pbImageImpDihedralID2;
    delete pbImageCmapID1;
    delete pbImageCmapID2;
    delete pbImageCmapID3;
    delete pbImageShakeID;
    delete pbImageFastShakeID;
    delete pbImageSlowShakeID1;
    delete pbImageSlowShakeID2;
    delete pbImageSolventAtomID;
    delete pbImageSoluteAtomID;
    delete pbNLNonbondCellStartEnd;
#ifndef use_SPFP
    delete pbNLChargeGridBufferOffset;
    delete pbNLOddBufferOverlapFlag;
#endif
    delete pbBNLExclusionBuffer;
    delete pbNLExclusionList;
    delete pbNLExclusionStartCount;
    delete pbNLAtomList;
    delete pbNLOffset;
    delete pbNLTotalOffset;
    delete pbNLPosition;
    delete pbNLbSkinTestFail;
    delete pbNLCellHash;

    // Delete NTP stuff
    delete pbSoluteAtomID;
    delete pbSoluteAtomMass;
    delete pbSolute;
    delete pbUllSolute;
    delete pbSolventAtomID;
    delete pbSolvent;
    delete pbNTPData;
    delete pbConstraintAtomX;
    delete pbConstraintAtomY;
    delete pbConstraintAtomZ;
    delete pbConstraintCOMX;
    delete pbConstraintCOMY;
    delete pbConstraintCOMZ;
    
    // Delete GB stuff
    delete pbWorkUnit;
    delete pbExclusion;
    delete pbGBPosition;
    delete pbNeckMaxValPos;
    delete pbGBAlphaBetaGamma;
    
    // Delete random number stuff
    delete pbRandomPos;
    delete pbRandom;

    // Delete bonded parameter data
    delete pbBond;
    delete pbBondID;
    delete pbBondAngle;
    delete pbBondAngleID1;
    delete pbBondAngleID2;
    delete pbDihedral1;
    delete pbDihedral2;
    delete pbDihedral3;
    delete pbDihedralID1;
    delete pbDihedralID2;
    delete pbNb141;
    delete pbNb142;
    delete pbNb14ID;
    delete pbConstraint1;
    delete pbConstraint2;
    delete pbConstraintID;
    delete pbUBAngle;
    delete pbUBAngleID;
    delete pbImpDihedral;
    delete pbImpDihedralID1; 
    delete pbImpDihedralID2;
    delete pbCmapID1;
    delete pbCmapID2;
    delete pbCmapID3;
    delete pbCmapType;
    delete pbCmapEnergy;
    
    // Delete NMR stuff
    delete pbNMRJarData;
    delete pbNMRDistanceID;
    delete pbNMRDistanceR1R2;
    delete pbNMRDistanceR3R4;
    delete pbNMRDistanceK2K3;
    delete pbNMRDistanceAve;
    delete pbNMRDistanceTgtVal;
    delete pbNMRDistanceStep; 
    delete pbNMRDistanceInc; 
    delete pbNMRDistanceR1R2Slp;
    delete pbNMRDistanceR3R4Slp;
    delete pbNMRDistanceK2K3Slp;
    delete pbNMRDistanceR1R2Int;
    delete pbNMRDistanceR3R4Int;
    delete pbNMRDistanceK2K3Int;
    delete pbNMRAngleID1;
    delete pbNMRAngleID2;
    delete pbNMRAngleR1R2;
    delete pbNMRAngleR3R4;
    delete pbNMRAngleK2K3;
    delete pbNMRAngleAve;
    delete pbNMRAngleTgtVal;
    delete pbNMRAngleStep;
    delete pbNMRAngleInc;
    delete pbNMRAngleR1R2Slp;
    delete pbNMRAngleR3R4Slp;
    delete pbNMRAngleK2K3Slp;
    delete pbNMRAngleR1R2Int; 
    delete pbNMRAngleR3R4Int; 
    delete pbNMRAngleK2K3Int;   
    delete pbNMRTorsionID1;
    delete pbNMRTorsionID2;
    delete pbNMRTorsionR1R2;
    delete pbNMRTorsionR3R4;
    delete pbNMRTorsionK2K3;
    delete pbNMRTorsionAve1;
    delete pbNMRTorsionAve2;
    delete pbNMRTorsionTgtVal;
    delete pbNMRTorsionStep;
    delete pbNMRTorsionInc;
    delete pbNMRTorsionR1R2Slp;
    delete pbNMRTorsionR3R4Slp;
    delete pbNMRTorsionK2K3Slp;
    delete pbNMRTorsionR1R2Int;
    delete pbNMRTorsionR3R4Int;
    delete pbNMRTorsionK2K3Int;    
    delete pbImageNMRDistanceID;
    delete pbImageNMRAngleID1;
    delete pbImageNMRAngleID2;
    delete pbImageNMRTorsionID1;
    delete pbImageNMRTorsionID2;
    
    // Delete Shake constraint data
    delete pbShakeID;
    delete pbShakeParm;
    delete pbFastShakeID;
    delete pbSlowShakeID1;
    delete pbSlowShakeID2;
    delete pbSlowShakeParm;

    // Delete extra points data
    delete pbExtraPoint11Frame;
    delete pbExtraPoint11Index;
    delete pbExtraPoint11;
    delete pbExtraPoint12Frame;
    delete pbExtraPoint12Index;
    delete pbExtraPoint12;
    delete pbExtraPoint21Frame;
    delete pbExtraPoint21Index;
    delete pbExtraPoint21;  
    delete pbExtraPoint22Frame;
    delete pbExtraPoint22Index;
    delete pbExtraPoint22;   
    delete pbImageExtraPoint11Frame;
    delete pbImageExtraPoint11Index;
    delete pbImageExtraPoint12Frame;
    delete pbImageExtraPoint12Index;
    delete pbImageExtraPoint21Frame;
    delete pbImageExtraPoint21Index;
    delete pbImageExtraPoint22Frame;
    delete pbImageExtraPoint22Index;
    
    // Delete output and/or accumulator buffers
#ifdef use_SPFP
    delete pbReffAccumulator;
    delete pbSumdeijdaAccumulator;
    delete pbForceAccumulator; 
#else
    delete pbReffBuffer;
    delete pbSumdeijdaBuffer;
    delete pbForceBuffer; 
#endif
    delete pbEnergyBuffer;
    delete pbKineticEnergyBuffer;

    // Delete AMD buffers
    delete  pAmdWeightsAndEnergy;       // AMD
    delete  pbAMDfwgtd;
}
