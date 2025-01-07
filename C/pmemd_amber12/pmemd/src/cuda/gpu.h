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

#ifndef __GPU_H__
#define __GPU_H__
#include "gputypes.h"
// F95 interface
extern "C" void gpu_init_(void);
extern "C" void gpu_shutdown_(void);
extern "C" void gpu_setup_system_(int* atoms, double* tol, int* ntf, int* ntb, int* ips, int* ntp, int* ntt, int* vrand);
extern "C" void gpu_upload_crd_(double atm_crd[][3]);
extern "C" void gpu_download_crd_(double atm_crd[][3]);
extern "C" void gpu_upload_charges_(double charge[]);
extern "C" void gpu_download_charges_(double charge[]);
extern "C" void gpu_upload_fs_(double fs[]);
extern "C" void gpu_download_fs_(double fs[]);
extern "C" void gpu_upload_rborn_(double rborn[]);
extern "C" void gpu_download_rborn_(double rborn[]);
extern "C" void gpu_upload_reff_(double reff[]);
extern "C" void gpu_download_reff_(double reff[]);
extern "C" void gpu_upload_frc_(double atm_frc[][3]);
extern "C" void gpu_upload_frc_add_(double atm_frc[][3]);
extern "C" void gpu_download_frc_(double atm_frc[][3]);
extern "C" void gpu_upload_vel_(double atm_vel[][3]);
extern "C" void gpu_download_vel_(double atm_vel[][3]);
extern "C" void gpu_upload_last_vel_(double atm_last_vel[][3]);
extern "C" void gpu_download_last_vel_(double atm_last_vel[][3]);
extern "C" void gpu_bonds_setup_(int* cit_nbona, bond_rec cit_a_bond[], int* cit_nbonh, bond_rec cit_h_bond[], double gbl_req[], double gbl_rk[]);
extern "C" void gpu_angles_setup_(int* angle_cnt, int* ntheth, angle_rec cit_angle[], double gbl_teq[], double gbl_tk[]);
extern "C" void gpu_dihedrals_setup_(int* dihed_cnt, int* nhpih, dihed_rec cit_dihed[], int gbl_ipn[], double gbl_pn[], double gbl_pk[], double gbl_gamc[], double gbl_gams[]);
extern "C" void gpu_nb14_setup_(int* cit_nb14_cnt, int cit_nb14[][3], double gbl_one_scee[], double gbl_one_scnb[], int* ntypes, int iac[], int ico[], double cn1[], double cn2[], double cn114[], double cn214[]);
extern "C" void gpu_molecule_list_setup(int* molecules, listdata_rec listdata[]);
extern "C" void gpu_create_outputbuffers_();
extern "C" void gpu_amd_setup_(int* iamd, int* iamdlag, int* ntwx, double* EthreshP, double* alphaP, double* EthreshD, double* alphaD, double* temp0);
extern "C" void gpu_calculate_and_apply_amd_weights_(double* pot_ene_tot, double* dih_ene_tot, double* num_amd_lag);
extern "C" void gpu_download_amd_weights_(double amd_weights_and_energy[]);
extern "C" void gpu_gbsa_setup_();
extern "C" void gpu_gbsa_frc_add_(double atm_frc[][3]);
#ifdef MPI
extern "C" void gpu_calculate_amd_dihedral_weight_(double* totdih);
extern "C" void gpu_calculate_amd_dihedral_energy_(double* totdih);
extern "C" void gpu_calculate_gb_amd_dihedral_energy_(double* totdih);
#else
extern "C" void gpu_calculate_amd_dihedral_energy_weight_();
extern "C" void gpu_calculate_gb_amd_dihedral_energy_weight_();
#endif

#ifdef CPU_RANDOMS
void cpu_amrset(int seed);
void cpu_kRandom(gpuContext gpu);
#endif

// Local interface
extern "C" void gpuCopyConstants();


// Kernel interfaces
extern "C" void kCalculateGBBornRadiiInitKernels(gpuContext gpu);
extern "C" void kCalculateGBNonbondEnergy1InitKernels(gpuContext gpu);
extern "C" void kCalculateGBNonbondEnergy2InitKernels(gpuContext gpu);
extern "C" void kNeighborListInitKernels(gpuContext gpu);
extern "C" void kCalculatePMENonbondEnergyInitKernels(gpuContext gpu);
extern "C" void kPMEInterpolationInitKernels(gpuContext gpu);  
extern "C" void kCalculateLocalForcesInitKernels(gpuContext gpu);
extern "C" void kShakeInitKernels(gpuContext gpu);
extern "C" void SetkForcesUpdateSim(gpuContext gpu);
extern "C" void GetkForcesUpdateSim(gpuContext gpu);
extern "C" void SetkCalculateLocalForcesSim(gpuContext gpu);
extern "C" void GetkCalculateLocalForcesSim(gpuContext gpu);
extern "C" void SetkCalculateGBBornRadiiSim(gpuContext gpu);
extern "C" void GetkCalculateGBBornRadiiSim(gpuContext gpu);
extern "C" void SetkCalculateGBNonbondEnergy1Sim(gpuContext gpu);
extern "C" void GetkCalculateGBNonbondEnergy1Sim(gpuContext gpu);
extern "C" void SetkCalculateGBNonbondEnergy2Sim(gpuContext gpu);
extern "C" void GetkCalculateGBNonbondEnergy2Sim(gpuContext gpu);
extern "C" void SetkShakeSim(gpuContext gpu);
extern "C" void GetkShakeSim(gpuContext gpu);
extern "C" void SetkNeighborListSim(gpuContext gpu);
extern "C" void GetkNeighborListSim(gpuContext gpu);
extern "C" void SetkPMEInterpolationSim(gpuContext gpu);
extern "C" void GetkPMEInterpolationSim(gpuContext gpu);
extern "C" void SetkCalculatePMENonbondEnergySim(gpuContext gpu);
extern "C" void GetkCalculatePMENonbondEnergySim(gpuContext gpu);
extern "C" void SetkCalculateAMDWeightsSim(gpuContext gpu);
extern "C" void GetkCalculateAMDWeightsSim(gpuContext gpu);

// Kernels
extern "C" void SetNLClearForcesKernel(gpuContext gpu);
extern "C" void kClearForces(gpuContext gpu);
extern "C" void kClearNBForces(gpuContext gpu);
extern "C" void SetNLReduceForcesKernel(gpuContext gpu);
extern "C" void kReduceForces(gpuContext gpu);
extern "C" void kReduceNBForces(gpuContext gpu);
extern "C" void kOrientForces(gpuContext gpu);
extern "C" void kLocalToGlobal(gpuContext gpu);
extern "C" void kCalculateLocalForces(gpuContext gpu);
extern "C" void kCalculateLocalEnergy(gpuContext gpu);
extern "C" void kCalculateCHARMMForces(gpuContext gpu);
extern "C" void kCalculateCHARMMEnergy(gpuContext gpu);
extern "C" void kCalculateNMRForces(gpuContext gpu);
extern "C" void kCalculateNMREnergy(gpuContext gpu);
extern "C" void kCalculateGBBornRadii(gpuContext gpu);
extern "C" void kReduceGBBornRadii(gpuContext gpu);
extern "C" void kProcessGBBornRadii(gpuContext gpu);
extern "C" void kClearGBBuffers(gpuContext gpu);
extern "C" void kCalculateGBNonbondEnergy1(gpuContext gpu);
extern "C" void kCalculateGBNonbondForces1(gpuContext gpu);
extern "C" void kReduceGBTemp7(gpuContext gpu);
extern "C" void kProcessGBTemp7(gpuContext gpu);
extern "C" void kReduceGBTemp7Energy(gpuContext gpu);
extern "C" void kProcessGBTemp7Energy(gpuContext gpu);
extern "C" void kCalculateGBNonbondEnergy2(gpuContext gpu);
extern "C" void kUpdate(gpuContext gpu, PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln);
extern "C" void kShake(gpuContext gpu);
extern "C" void kFastShake(gpuContext gpu);
extern "C" void kCalculateKineticEnergy(gpuContext gpu, PMEFloat c_ave);
extern "C" void kCalculateCOM(gpuContext gpu);
extern "C" void kCalculateSoluteCOM(gpuContext gpu);
extern "C" void kReduceSoluteCOM(gpuContext gpu);
extern "C" void kCalculateCOMKineticEnergy(gpuContext gpu);
extern "C" void kReduceCOMKineticEnergy(gpuContext gpu);
extern "C" void kCalculateMolecularVirial(gpuContext gpu);
extern "C" void kCalculateSoluteCOM(gpuContext gpu);
extern "C" void kPressureScaleCoordinates(gpuContext gpu);
extern "C" void kPressurScaleConstraints(gpuContext gpu);
extern "C" void kPressureScaleConstraintCoordinates(gpuContext gpu);
extern "C" void kResetVelocities(gpuContext gpu, double temp, double half_dtx);
extern "C" void kRecalculateVelocities(gpuContext gpu, PMEDouble dtx_inv);
extern "C" void kScaleVelocities(gpuContext gpu, PMEDouble scale);
extern "C" void kRecenter_Molecule(gpuContext gpu);
extern "C" void kRandom(gpuContext gpu);
extern "C" void kPMEGetGridWeights(gpuContext gpu);
extern "C" void kPMEClearChargeGrid(gpuContext gpu);
extern "C" void kPMEFillChargeGrid(gpuContext gpu);
extern "C" void kPMEConvertChargeGrid(gpuContext gpu);
extern "C" void kPMEClearChargeGridBuffer(gpuContext gpu);
extern "C" void kPMEFillChargeGridBuffer(gpuContext gpu);
extern "C" void kPMEReduceChargeGridBuffer(gpuContext gpu);
extern "C" void kPMEScalarSumRC(gpuContext gpu, PMEDouble ewaldcof, PMEDouble vol);
extern "C" void kPMEScalarSumRCEnergy(gpuContext gpu, PMEDouble ewaldcof, PMEDouble vol);
extern "C" void kPMEC2CScalarSumRC(gpuContext gpu, PMEDouble ewaldcof, PMEDouble vol);
extern "C" void kPMEGradSum(gpuContext gpu);
extern "C" void kNLGenerateSpatialHash(gpuContext gpu);
extern "C" void kNLInitRadixSort(gpuContext gpu);
extern "C" void kNLDeleteRadixSort(gpuContext gpu);
extern "C" void kNLRadixSort(gpuContext gpu);
extern "C" void kNLRemapLocalInteractions(gpuContext gpu);
extern "C" void kNLRemapImage(gpuContext gpu);
extern "C" void kNLCalculateOffsets(gpuContext gpu);
extern "C" void kNLCalculateCellCoordinates(gpuContext gpu);
extern "C" void kNLBuildNeighborList(gpuContext gpu);
extern "C" void kNLClearCellBoundaries(gpuContext gpu);
extern "C" void kNLCalculateCellBoundaries(gpuContext gpu);
extern "C" void kCalculatePMENonbondEnergy(gpuContext gpu);
extern "C" void kCalculatePMENonbondForces(gpuContext gpu);
extern "C" void kCalculatePMELocalForces(gpuContext gpu);
extern "C" void kCalculatePMELocalEnergy(gpuContext gpu);
extern "C" void kCalculateIPSNonbondEnergy(gpuContext gpu);
extern "C" void kCalculateIPSNonbondForces(gpuContext gpu);
extern "C" void kNLSkinTest(gpuContext gpu);
extern "C" void kCalculateAMDWeights(gpuContext gpu);
extern "C" void kAMDCalcWeightAndScaleForces(gpuContext gpu, PMEDouble pot_ene_tot, PMEDouble dih_ene_tot, PMEDouble fwgt);
extern "C" void kCalculateAmdDihedralEnergy(gpuContext gpu);
#ifdef MPI
extern "C" void kTransposeForces(gpuContext gpu);
#endif
#endif


