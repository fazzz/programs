# Platform-specific info should be found in config.h

include ../../../config.h

CU_OBJS = cuda_info.o gpu.o gputypes.o \
	kForcesUpdate.o kCalculateLocalForces.o kCalculateGBBornRadii.o \
	kCalculatePMENonbondEnergy.o kCalculateGBNonbondEnergy1.o kNLRadixSort.o \
        kCalculateGBNonbondEnergy2.o kShake.o kNeighborList.o kPMEInterpolation.o \
        kCalculateAMDWeights.o

.SUFFIXES: .F90 .cpp .cu .o

# CUDA compilations flags
kForcesUpdate.cu.CU_FLAGS                   =
kCalculateLocalForces.cu.CU_FLAGS           =
kCalculateGBBornRadii.cu.CU_FLAGS           =
kCalculateGBNonbondEnergy1.cu.CU_FLAGS      =
kCalculateGBNonbondEnergy2.cu.CU_FLAGS      =
kRandom.cu.CU_FLAGS                         =
kShake.cu.CU_FLAGS                          =
kPMEInterpolation.cu.CU_FLAGS               =
kNeighborList.cu.CU_FLAGS                   =
kCalculatePMENonbondEnergy.cu.CU_FLAGS      =
kCalculateAMDWeights.cu.CU_FLAGS            =

cuda.a: $(CU_OBJS) 
	ar rvs cuda.a $(CU_OBJS)

.F90.o:
	$(PMEMD_F90) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(PMEMD_CU_INCLUDES) -c $*.F90

.cpp.o:
	$(PMEMD_CC) $(PMEMD_COPTFLAGS) $(PMEMD_CU_DEFINES) $(PMEMD_CU_INCLUDES) -c $*.cpp

.cu.o:  gputypes.h
	$(NVCC) $(PMEMD_CU_DEFINES) $(PMEMD_CU_INCLUDES) $($?.CU_FLAGS) -c $?

clean:
	rm -f *.o *.linkinfo cuda.a *.mod

