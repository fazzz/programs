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
#include <stdio.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h>

#ifdef WINDOWS
#include "winfunctions.h"
#endif

using namespace std;

#define GPU_CPP
#include "gputypes.h"
#include "gpu.h"
#undef GPU_CPP

#ifdef MPI
extern "C" void gpu_startup_(int* gpuID, int* nGpus, int* comm_number) 
#else
extern "C" void gpu_startup_(void)
#endif
{
    gpu = new _gpuContext;
PRINTMETHOD("gpu_startup"); 
#ifdef MPI
    gpu->gpuID                          = *gpuID;
    gpu->nGpus                          = *nGpus;
    gpu->comm                           = MPI_COMM_WORLD;
#ifdef GVERBOSE    
    printf("Node %d of %d nodes\n", gpu->gpuID, gpu->nGpus);
#endif

    // Create local communicator from comm number
    MPI_Comm_split(MPI_COMM_WORLD, *comm_number, *gpuID, &gpu->comm);
#endif

}

extern "C" void gpu_set_device_(int* device)
{
PRINTMETHOD("gpu_set_device"); 
    gpu->gpu_device_id                              = *device;
}

extern "C" void gpu_init_(void)
{

PRINTMETHOD("gpu_init");
    int LRFSize = 0;
    int SMCount = 0;
    int SMMajor = 0;
    int SMMinor = 0;
    
#ifdef MPI
    if(getenv("CUDA_PROFILE") != 0) 
    {
      char profile_log[80];
      if(getenv("CUDA_PROFILE_LOG")) {
        sprintf(profile_log, "../%s%d", getenv("CUDA_PROFILE_LOG"), gpu->gpuID);
      } else {
        sprintf(profile_log, "../cu%d.log", gpu->gpuID);
      }
#ifndef WINDOWS
	  setenv("CUDA_PROFILE_LOG", profile_log, 1);
#else
	  char env[100];
	  strcpy(env, "CUDA_PROFILE_LOG");
	  strcat(env, profile_log);
	  _putenv(env);
#endif
    }
#endif

    int device	                                    = -1;
    int gpuCount                                    = 0;
    cudaError_t status;
    cudaDeviceProp deviceProp;
    status = cudaGetDeviceCount(&gpuCount);
    RTERROR(status, "cudaGetDeviceCount failed");
    if (gpuCount == 0)
    {
        printf("No CUDA-capable devices found, exiting.\n");
        cudaThreadExit();
        exit(-1);
    }

#ifdef MPI
    // Grab node names from all other processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int length;
    char myName[MPI_MAX_PROCESSOR_NAME + 1];
    char* pName                                     = new char[world_size * (MPI_MAX_PROCESSOR_NAME + 1)];
    int* pNameCount                                 = new int[world_size];
    int* pNameDisp                                  = new int[world_size];
    MPI_Get_processor_name(myName, &length);
    strcpy(&pName[world_rank * (MPI_MAX_PROCESSOR_NAME + 1)], myName); 
    for (int i = 0; i < world_size; i++)
    {
        pNameCount[i]                               = MPI_MAX_PROCESSOR_NAME + 1;
        pNameDisp[i]                                = i * (MPI_MAX_PROCESSOR_NAME + 1);
    }
    MPI_Allgatherv(myName, MPI_MAX_PROCESSOR_NAME + 1, MPI_CHAR, pName, pNameCount, pNameDisp, 
            MPI_CHAR, MPI_COMM_WORLD);

    // Test for single node run
    bool bSingleNode = true;
    for (int i = 0; i < gpu->nGpus; i++)
    {
        if (!strcmp(&pName[i * (MPI_MAX_PROCESSOR_NAME + 1)], myName))
            bSingleNode                             = false;
    }
#endif

    // If the device id is -1 this means it was left at the default so we 
    // or the user specified -1 on the command line meaning autoselect GPU.
    // choose CUDA device with the most memory (flesh out later for multi-gpu)
    // otherwise we use the device id specified on the command line.
    if ( gpu->gpu_device_id == -1 ) 
    {
#ifdef MPI       
        // Check for duplicate processes on current node
        int localCount                              = 0;
        int offset                                  = 1;

        for (int i = 0; i < world_size; i++)
        {
            if (!strcmp(&pName[i * (MPI_MAX_PROCESSOR_NAME + 1)], myName))
            {
                localCount++;
                if (i < world_rank)
                    offset++;
            }
        }
        
        if (localCount > 1)
        {
            // Choose nth gpu that can run AMBER
            int pos = 0;
            while (offset > 0)
            {
                cudaGetDeviceProperties(&deviceProp, pos);
                if (deviceProp.canMapHostMemory && 
                    ((deviceProp.major >= 2) || ((deviceProp.major == 1) && (deviceProp.minor == 3))))          
                {
                    device                          = pos;
                    offset--;
                }   
                pos++;
                if (pos == gpuCount)
                    pos                             = 0;         
            }       
#ifdef GVERBOSE            
            printf("Node %d running on device %d\n", gpu->gpuID, device);
#endif
        }  
        else
#endif
        {
            // Choose GPU with the most memory
            size_t maxMem                           = 0;            
            for (int i = 0; i < gpuCount; i++)
            {
                cudaGetDeviceProperties(&deviceProp, i);

                if (
#ifdef MPI
                    deviceProp.canMapHostMemory && 
#endif
                    ((deviceProp.major >= 2) || ((deviceProp.major == 1) && (deviceProp.minor == 3))) &&
                     (deviceProp.totalGlobalMem >= maxMem))          
                {
                    maxMem                          = deviceProp.totalGlobalMem;
                    device                          = i;
                }
            }
        
        }           
    }
    else
    {
        cudaGetDeviceProperties(&deviceProp, gpu->gpu_device_id);
#ifdef MPI
        if (deviceProp.canMapHostMemory && (deviceProp.major >= 2) || ((deviceProp.major == 1) && (deviceProp.minor == 3)))
#else        
        if ((deviceProp.major >= 2) || ((deviceProp.major == 1) && (deviceProp.minor == 3)))
#endif
            device = gpu->gpu_device_id;
        else
        {
#ifdef MPI
            printf("Selected GPU does not support both zero-copy and double-precision, exiting.\n");
#else        
            printf("Selected GPU lacks double-precision support, exiting.\n");
#endif            
            cudaThreadExit();
            exit(-1);
        }
    }

#ifdef MPI
    // Release list of gpu names
    delete[] pName;
    delete[] pNameCount;
    delete[] pNameDisp;
#endif

    if (device == -1)
    {
#ifdef MPI
        printf("No zero-copy and double-precision capable gpu located, exiting.\n");
#else    
        printf("No double-precision capable gpu located, exiting.\n");
#endif        
        cudaThreadExit();
        exit(-1);
    }

#if defined(MPI) && defined(CUDA_P2P)
    // Test for universal P2P access
    if (bSingleNode)
    {
        int* pDevice                                = new int[gpu->nGpus];
        pDevice[gpu->gpuID]                         = device;
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, pDevice, sizeof(int), MPI_BYTE, gpu->comm);
        for (int i = 0; i < gpu->nGpus; i++)
        {     
            if (i != gpu->gpuID)
            {
                int canAccessPeer;
                cudaError_t status                  = cudaDeviceCanAccessPeer(&canAccessPeer, device, pDevice[i]);
                RTERROR(status, "cudaDeviceCanAccessPeer");
                if (canAccessPeer == 0)
                    bSingleNode                     = false;
            }
        }
        delete[] pDevice;
    }
    gpu->bSingleNode                                = bSingleNode;
    gpu->bP2P                                       = bSingleNode;
#else
    // Activate zero-copy
    cudaSetDeviceFlags(cudaDeviceMapHost);
#endif
    
    status = cudaSetDevice(device); 
    RTERROR(status, "Error setting CUDA device");  
    cudaThreadSynchronize();

    // Determine kernel call configuration and grab desired additional GPU properties
    cudaGetDeviceProperties(&deviceProp, device);
    gpu->bECCSupport                                = deviceProp.ECCEnabled;
    gpu->totalMemory                                = deviceProp.totalGlobalMem;

#ifdef GVERBOSE
    double memsize = (double)deviceProp.totalGlobalMem / (1024.0 * 1024.0);
    printf("Using GPU %d, %s, SM %d.%d, %.1f MBytes of memory\n", device, deviceProp.name, deviceProp.major, deviceProp.minor, memsize);
#endif

    // Store GPU Device ID for later use
    gpu->gpu_device_id                              = device;
    
    gpu->blocks                                     = deviceProp.multiProcessorCount;
    gpu->GBBornRadiiBlocks                          = gpu->blocks;
    gpu->GBNonbondEnergy1Blocks                     = gpu->blocks;
    gpu->GBNonbondEnergy2Blocks                     = gpu->blocks; 

    // Determine SM version
    unsigned int blocksPerSM;
    if (deviceProp.major == 1)
    {
        switch (deviceProp.minor)
        {
        case 0:
        case 1:
        case 2:
        case 5:
            printf("GPU SM revision must be 1.3 or better.\n");
            gpu_shutdown_();
            exit(-1);
            break;
        }

        gpu->sm_version                             = SM_13;
        gpu->threadsPerBlock                        = SM_13_THREADS_PER_BLOCK;
        gpu->NLCalculateOffsetsThreadsPerBlock      = SM_13_NLCALCULATE_OFFSETS_THREADS_PER_BLOCK;
        gpu->NLBuildNeighborList32ThreadsPerBlock   = SM_13_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK;
        gpu->NLBuildNeighborList16ThreadsPerBlock   = SM_13_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK;
        gpu->localForcesThreadsPerBlock             = SM_13_LOCALFORCES_THREADS_PER_BLOCK;
        gpu->CHARMMForcesThreadsPerBlock            = SM_13_CHARMMFORCES_THREADS_PER_BLOCK;
        gpu->clearForcesThreadsPerBlock             = SM_13_CLEARFORCES_THREADS_PER_BLOCK;
        gpu->NLClearForcesThreadsPerBlock           = SM_13_NLCLEARFORCES_THREADS_PER_BLOCK;
        gpu->reduceForcesThreadsPerBlock            = SM_13_REDUCEFORCES_THREADS_PER_BLOCK;
        gpu->NLReduceForcesThreadsPerBlock          = SM_13_NLREDUCEFORCES_THREADS_PER_BLOCK;
        gpu->reduceBufferThreadsPerBlock            = SM_13_REDUCEBUFFER_THREADS_PER_BLOCK;
        gpu->GBBornRadiiThreadsPerBlock             = SM_13_GBBORNRADII_THREADS_PER_BLOCK;
        gpu->GBBornRadiiIGB78ThreadsPerBlock        = SM_13_GBBORNRADIIIGB78_THREADS_PER_BLOCK;
        gpu->GBNonbondEnergy1ThreadsPerBlock        = SM_13_GBNONBONDENERGY1_THREADS_PER_BLOCK;
        gpu->GBNonbondEnergy2ThreadsPerBlock        = SM_13_GBNONBONDENERGY2_THREADS_PER_BLOCK;
        gpu->GBNonbondEnergy2IGB78ThreadsPerBlock   = SM_13_GBNONBONDENERGY2IGB78_THREADS_PER_BLOCK;        
        gpu->PMENonbondEnergyThreadsPerBlock        = SM_13_PMENONBONDENERGY_THREADS_PER_BLOCK;
        gpu->IPSNonbondEnergyThreadsPerBlock        = SM_13_IPSNONBONDENERGY_THREADS_PER_BLOCK;
        gpu->updateThreadsPerBlock                  = SM_13_UPDATE_THREADS_PER_BLOCK;
        gpu->shakeThreadsPerBlock                   = SM_13_SHAKE_THREADS_PER_BLOCK;
        gpu->readSize                               = SM_13_READ_SIZE;
        gpu->maxSoluteMolecules                     = SM_13_MAXMOLECULES;
        gpu->maxPSSoluteMolecules                   = SM_13_MAXPSMOLECULES;
    }
    else if (deviceProp.major == 2)
    {
        gpu->sm_version                             = SM_2X;
        gpu->threadsPerBlock                        = SM_2X_THREADS_PER_BLOCK;
        gpu->NLCalculateOffsetsThreadsPerBlock      = SM_2X_NLCALCULATE_OFFSETS_THREADS_PER_BLOCK;
        gpu->NLBuildNeighborList32ThreadsPerBlock   = SM_2X_NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK;
        gpu->NLBuildNeighborList16ThreadsPerBlock   = SM_2X_NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK;
        gpu->localForcesThreadsPerBlock             = SM_2X_LOCALFORCES_THREADS_PER_BLOCK;
        gpu->CHARMMForcesThreadsPerBlock            = SM_2X_CHARMMFORCES_THREADS_PER_BLOCK;
        gpu->clearForcesThreadsPerBlock             = SM_2X_CLEARFORCES_THREADS_PER_BLOCK;
        gpu->NLClearForcesThreadsPerBlock           = SM_2X_NLCLEARFORCES_THREADS_PER_BLOCK;
        gpu->reduceForcesThreadsPerBlock            = SM_2X_REDUCEFORCES_THREADS_PER_BLOCK;
        gpu->NLReduceForcesThreadsPerBlock          = SM_2X_NLREDUCEFORCES_THREADS_PER_BLOCK;
        gpu->reduceBufferThreadsPerBlock            = SM_2X_REDUCEBUFFER_THREADS_PER_BLOCK;
        gpu->GBBornRadiiThreadsPerBlock             = SM_2X_GBBORNRADII_THREADS_PER_BLOCK;
        gpu->GBBornRadiiIGB78ThreadsPerBlock        = SM_2X_GBBORNRADII_THREADS_PER_BLOCK;        
        gpu->GBNonbondEnergy1ThreadsPerBlock        = SM_2X_GBNONBONDENERGY1_THREADS_PER_BLOCK;
        gpu->GBNonbondEnergy2ThreadsPerBlock        = SM_2X_GBNONBONDENERGY2_THREADS_PER_BLOCK;
        gpu->GBNonbondEnergy2IGB78ThreadsPerBlock   = SM_2X_GBNONBONDENERGY2_THREADS_PER_BLOCK;        
        gpu->PMENonbondEnergyThreadsPerBlock        = SM_2X_PMENONBONDENERGY_THREADS_PER_BLOCK;
        gpu->IPSNonbondEnergyThreadsPerBlock        = SM_2X_IPSNONBONDENERGY_THREADS_PER_BLOCK;
        gpu->updateThreadsPerBlock                  = SM_2X_UPDATE_THREADS_PER_BLOCK;
        gpu->shakeThreadsPerBlock                   = SM_2X_SHAKE_THREADS_PER_BLOCK;
        gpu->readSize                               = SM_2X_READ_SIZE;
        gpu->maxSoluteMolecules                     = SM_2X_MAXMOLECULES;
        gpu->maxPSSoluteMolecules                   = SM_2X_MAXPSMOLECULES;
    }
    else
    {
        printf("Currently only GPU SM Revisions 1.3 and 2.X are supported.\n");
        gpu_shutdown_();
        exit(-1);
    }
    gpu->bNeighborList                              = false;
#ifdef GVERBOSE
    printf("LRF size is %d\n", deviceProp.regsPerBlock);
#endif  
    return;
}

extern "C" void gpu_shutdown_(void)
{
PRINTMETHOD("gpu_shutdown");
#ifdef MPI
    if (gpu->bNeighborList)
    {
        MPI_Win_free(&gpu->MPIPMEForceWindow); 
        MPI_Free_mem(gpu->pForceData);   
    }
#ifdef CUDA_P2P
    if (gpu->bP2P)
    {
        for (int i = 0; i < gpu->nGpus; i++)
        {
            if (i != gpu->gpuID)
            {
                cudaError_t status;
                status                              = cudaIpcCloseMemHandle((void*)gpu->pP2PReffa[i]);
                RTERROR(status, "cudaIpcCloseMemHandle failed on gpu->pbReffa");
                status                              = cudaIpcCloseMemHandle((void*)gpu->pP2PTemp7a[i]);
                RTERROR(status, "cudaIpcCloseMemHandle failed on gpu->pbTemp7a");
                status                              = cudaIpcCloseMemHandle((void*)gpu->pP2PInForce[i]);
                RTERROR(status, "cudaIpcCloseMemHandle failed on gpu->pbInForce");
            }
        }

        // Make sure all processes reach here before deleting GPU memory
        MPI_Barrier(gpu->comm);
    }
#endif
#endif
    delete gpu;
    cudaThreadExit();
    return;
}

extern "C" void gpu_get_device_info_(int* gpu_dev_count, int* gpu_dev_id, int* gpu_dev_mem, int* gpu_num_proc, double* gpu_core_freq, char* gpu_name, int* name_len)
{
PRINTMETHOD("gpu_get_device_info");
    cudaError_t status;
    cudaDeviceProp deviceProp;
    size_t device_mem;

    status = cudaGetDeviceCount(gpu_dev_count);
    RTERROR(status, "cudaGetDeviceCount failed");

    *gpu_dev_id = gpu->gpu_device_id;

    cudaGetDeviceProperties(&deviceProp, *gpu_dev_id);
    device_mem = deviceProp.totalGlobalMem/(1024*1024);
    *gpu_dev_mem = (int )device_mem;

    *gpu_num_proc = (int )deviceProp.multiProcessorCount;
    *gpu_core_freq = (double )(deviceProp.clockRate * 1e-6f);

    strcpy(gpu_name,deviceProp.name);
    *name_len=strlen(deviceProp.name);
}

#ifdef MPI
extern "C" void gpu_get_slave_device_info_(int* taskID, int* gpu_dev_count, int* gpu_dev_id, int* gpu_dev_mem, int* gpu_num_proc, double* gpu_core_freq, char* gpu_name, int* name_len)
{
PRINTMETHOD("gpu_get_slave_device_info");
    MPI_Status status;
    MPI_Recv(gpu_dev_count, 1, MPI_INT, *taskID, 0, gpu->comm, &status);
    MPI_Recv(gpu_dev_id, 1, MPI_INT, *taskID, 0, gpu->comm, &status);
    MPI_Recv(gpu_dev_mem, 1, MPI_INT, *taskID, 0, gpu->comm, &status);
    MPI_Recv(gpu_num_proc, 1, MPI_INT, *taskID, 0, gpu->comm, &status);
    MPI_Recv(gpu_core_freq, 1, MPI_DOUBLE, *taskID, 0, gpu->comm, &status);
    MPI_Recv(gpu_name, 80, MPI_CHAR, *taskID, 0, gpu->comm, &status);
    MPI_Recv(name_len, 1, MPI_INT, *taskID, 0, gpu->comm, &status);
}

extern "C" void gpu_send_slave_device_info_()
{
PRINTMETHOD("gpu_send_slave_device_info");
    int gpu_dev_count;
    int gpu_dev_id;
    int gpu_dev_mem;
    int gpu_num_proc;
    double gpu_core_freq;
    char gpu_name[81];
    int name_len;
    gpu_get_device_info_(&gpu_dev_count, &gpu_dev_id, &gpu_dev_mem, &gpu_num_proc, &gpu_core_freq, gpu_name, &name_len);
    MPI_Send(&gpu_dev_count, 1, MPI_INT, 0, 0, gpu->comm);
    MPI_Send(&gpu_dev_id, 1, MPI_INT, 0, 0, gpu->comm);
    MPI_Send(&gpu_dev_mem, 1, MPI_INT, 0, 0, gpu->comm);
    MPI_Send(&gpu_num_proc, 1, MPI_INT, 0, 0, gpu->comm);
    MPI_Send(&gpu_core_freq, 1, MPI_DOUBLE, 0, 0, gpu->comm);
    MPI_Send(gpu_name, 80, MPI_CHAR, 0, 0, gpu->comm);
    MPI_Send(&name_len, 1, MPI_INT, 0, 0, gpu->comm);    
}
#endif

// Returns KB of memory in use on CPU and GPU
extern "C" void gpu_get_memory_info_(int* gpumemory, int* cpumemory)
{
    *gpumemory                          = (int)(gpu->totalGPUMemory / 1024ll);
    *cpumemory                          = (int)(gpu->totalCPUMemory / 1024ll);
    return;
}

extern "C" void gpu_setup_system_(int* atoms, double* tol, int* ntf, int* ntb, int* ntp, int* ntt, int* vrand)
{
PRINTMETHOD("gpu_setup_system");

    // Grab simulation parameters;
    gpu->ntf                            = *ntf;
    gpu->ntb                            = *ntb;
    gpu->ntt                            = *ntt;
    gpu->sim.ntp                        = *ntp;

    // Allocate system based on atom count
    gpu->sim.atoms                      = *atoms;
    gpu->sim.paddedNumberOfAtoms        = ((*atoms + gpu->sim.grid - 1) >> gpu->sim.gridBits) << gpu->sim.gridBits;

    // Calculate stride to insure buffers begin on safe texture boundaries
    if (gpu->sm_version < SM_2X)
    {
#ifndef use_SPSP
        gpu->sim.stride                 = ((*atoms + 31) >> 5) << 5;
#else
        gpu->sim.stride                 = ((*atoms + 63) >> 6) << 6;
#endif
    }
    else
    {
#ifndef use_SPSP
        gpu->sim.stride                 = ((*atoms + 63) >> 6) << 6;
#else
        gpu->sim.stride                 = ((*atoms + 127) >> 7) << 7;
#endif    
    }
    gpu->sim.stride2                    = 2 * gpu->sim.stride;
    gpu->sim.stride3                    = 3 * gpu->sim.stride;
    gpu->sim.stride4                    = 4 * gpu->sim.stride;
    gpu->sim.tol                        = *tol;
#ifdef GVERBOSE
    printf("I see %d atoms, padded to %d\n", *atoms, gpu->sim.paddedNumberOfAtoms);
#endif

    // Determine number of randoms to generate
    if ((gpu->ntt == 3) || (gpu->ntt == 2))
    {
        gpu->sim.randomSteps            = MAX_RANDOM_STEPS;                                
        if (gpu->totalMemory < aligned_uli(2u * 1024u * 1024u * 1024u )) // - 1 here to avoid overflowing a signed int. = 2GB.
            gpu->sim.randomSteps       /= 2;
        if (*atoms >= 131072)
            gpu->sim.randomSteps       /= 2;
        if (*atoms >= 262144)
            gpu->sim.randomSteps       /= 2;
        if (*atoms >= 524288)
            gpu->sim.randomSteps       /= 2;
        gpu->sim.randomNumbers          = gpu->sim.paddedNumberOfAtoms * 3 * gpu->sim.randomSteps;
    }
    else
    {
        gpu->sim.randomNumbers          = 1;
    }

    // Clear any previous stuff
    delete gpu->pbAtom;
    delete gpu->pbAtomXYSP;
    delete gpu->pbAtomZSP;
    delete gpu->pbAtomSigEps;
    delete gpu->pbAtomRBorn;
    delete gpu->pbAtomS;
    delete gpu->pbAtomCharge;
    delete gpu->pbAtomChargeSP;
    delete gpu->pbAtomMass;
    delete gpu->pbReff;
    gpu->pbReff                         = NULL;
    delete gpu->pbReffSP;
    delete gpu->pbPsi;
    delete gpu->pbTemp7;
    delete gpu->pbForce;
    delete gpu->pbUllForce;
#ifdef MPI
    delete gpu->pbInForce;
    delete gpu->pbOutForce;
    delete gpu->pbPMEForce;
    delete gpu->pbTemp7a;
    gpu->pbTemp7a                       = NULL;
    delete gpu->pbReffa;
    gpu->pbReffa                        = NULL;
#ifdef CUDA_P2P
    delete[] gpu->pP2PReffa;
    gpu->pP2PReffa                      = NULL;
    delete[] gpu->pP2PTemp7a;
    gpu->pP2PTemp7a                     = NULL;
    delete[] gpu->pP2PInForce;
#endif
#endif    
    delete gpu->pbVel;
    delete gpu->pbLVel;
    delete gpu->pbCenter;
    delete gpu->pbOutputBufferCounter;
    delete gpu->pbRandomPos;
    delete gpu->pbRandom;
    gpu->pbAtom                         = new GpuBuffer<PMEDouble>(gpu->sim.stride3);
    gpu->pbAtomXYSP                     = new GpuBuffer<PMEFloat2>(gpu->sim.stride);
    gpu->pbAtomZSP                      = new GpuBuffer<PMEFloat>(gpu->sim.stride);
    gpu->pbAtomCharge                   = new GpuBuffer<PMEDouble>(gpu->sim.stride);
    gpu->pbAtomChargeSP                 = new GpuBuffer<PMEFloat>(gpu->sim.stride);
    gpu->pbAtomSigEps                   = new GpuBuffer<PMEFloat2>(gpu->sim.stride);
    gpu->pbAtomRBorn                    = new GpuBuffer<PMEFloat>(gpu->sim.stride);
    gpu->pbAtomS                        = new GpuBuffer<PMEFloat>(gpu->sim.stride);
    gpu->pbAtomMass                     = new GpuBuffer<PMEDouble>(gpu->sim.stride2);
    gpu->pbReff                         = new GpuBuffer<PMEDouble>(gpu->sim.stride);
    if (gpu->ntb == 0)
    {
        gpu->pbReffSP                   = new GpuBuffer<PMEFloat>(gpu->sim.stride);
        gpu->pbPsi                      = new GpuBuffer<PMEFloat>(gpu->sim.stride);
        gpu->pbTemp7                    = new GpuBuffer<PMEFloat>(gpu->sim.stride);
#ifdef MPI
#ifdef CUDA_P2P
        if (gpu->bP2P)
        {
            gpu->pbReffa                = new GpuBuffer<PMEDouble>(gpu->sim.stride);
            gpu->pbTemp7a               = new GpuBuffer<PMEDouble>(gpu->sim.stride);
            gpu->pbInForce              = new GpuBuffer<PMEDouble>(gpu->sim.stride3);
            if (gpu->bSingleNode)
            {
                gpu->pP2PReffaHandle    = new cudaIpcMemHandle_t[gpu->nGpus];
                gpu->pP2PTemp7aHandle   = new cudaIpcMemHandle_t[gpu->nGpus];
                gpu->pP2PInForceHandle  = new cudaIpcMemHandle_t[gpu->nGpus];
                gpu->pP2PReffa          = new PMEDouble*[gpu->nGpus];
                gpu->pP2PTemp7a         = new PMEDouble*[gpu->nGpus];
                gpu->pP2PInForce        = new PMEDouble*[gpu->nGpus];
                cudaError_t status;
                status                  = cudaIpcGetMemHandle(&(gpu->pP2PReffaHandle[gpu->gpuID]), gpu->pbReffa->_pDevData);
                RTERROR(status, "cudaIpcGetMemHandle failed on gpu->pbReffa");
                status                  = cudaIpcGetMemHandle(&(gpu->pP2PTemp7aHandle[gpu->gpuID]), gpu->pbTemp7a->_pDevData);
                RTERROR(status, "cudaIpcGetMemHandle failed on gpu->pbTemp7a");
                status                  = cudaIpcGetMemHandle(&(gpu->pP2PInForceHandle[gpu->gpuID]), gpu->pbInForce->_pDevData);
                RTERROR(status, "cudaIpcGetMemHandle failed on gpu->pbInForce");
                MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, gpu->pP2PReffaHandle, sizeof(cudaIpcMemHandle_t), MPI_BYTE, gpu->comm);
                MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, gpu->pP2PTemp7aHandle, sizeof(cudaIpcMemHandle_t), MPI_BYTE, gpu->comm);
                MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, gpu->pP2PInForceHandle, sizeof(cudaIpcMemHandle_t), MPI_BYTE, gpu->comm);
                for (int i = 0; i < gpu->nGpus; i++)
                {
                    if (i != gpu->gpuID)
                    {
                        status          = cudaIpcOpenMemHandle((void**)&(gpu->pP2PReffa[i]), gpu->pP2PReffaHandle[i], cudaIpcMemLazyEnablePeerAccess);
                        RTERROR(status, "cudaIpcOpenMemHandle failed on gpu->pbReffaHandle");
                        status          = cudaIpcOpenMemHandle((void**)&(gpu->pP2PTemp7a[i]), gpu->pP2PTemp7aHandle[i], cudaIpcMemLazyEnablePeerAccess);
                        RTERROR(status, "cudaIpcOpenMemHandle failed on gpu->pbTemp7aHandle");
                        status          = cudaIpcOpenMemHandle((void**)&(gpu->pP2PInForce[i]), gpu->pP2PInForceHandle[i], cudaIpcMemLazyEnablePeerAccess);
                        RTERROR(status, "cudaIpcOpenMemHandle failed on gpu->pbInForceHandle");
                    }
                }
            }
        }
        else
#endif
        {
            gpu->pbReffa                = new GpuBuffer<PMEDouble>(gpu->sim.stride, false, true);
            gpu->pbTemp7a               = new GpuBuffer<PMEDouble>(gpu->sim.stride, false, true);
            gpu->pbInForce              = new GpuBuffer<PMEDouble>(gpu->sim.stride3, false, true);
        }
    }
    else 
    {
        gpu->pbPMEForce                 = new GpuBuffer<PMEFloat>(gpu->sim.stride3, false, true);
        gpu->pbInForce                  = new GpuBuffer<PMEDouble>(gpu->sim.stride3, false, true);
    } 
    gpu->pbOutForce                     = new GpuBuffer<PMEDouble>(gpu->sim.stride3);
#else
    }
#endif
    if (gpu->sim.ntp > 0)
    {
        gpu->pbForce                    = new GpuBuffer<PMEDouble>(gpu->sim.stride * 6);
        gpu->pbUllForce                 = new GpuBuffer<unsigned long long>(gpu->sim.stride3);   
    }
    else
        gpu->pbForce                    = new GpuBuffer<PMEDouble>(gpu->sim.stride3);

    gpu->pbVel                          = new GpuBuffer<PMEDouble>(gpu->sim.stride3);
    gpu->pbLVel                         = new GpuBuffer<PMEDouble>(gpu->sim.stride3);
    gpu->pbCenter                       = new GpuBuffer<PMEFloat>(gpu->blocks * 6);
    gpu->pbOutputBufferCounter          = new GpuBuffer<unsigned int>(gpu->sim.stride);
    gpu->pbRandomPos                    = new GpuBuffer<unsigned int>(gpu->blocks);
    gpu->pbRandom                       = new GpuBuffer<PMEDouble>(gpu->sim.randomNumbers, bShadowedOutputBuffers);
    gpu->sim.pAtomX                     = gpu->pbAtom->_pDevData;
    gpu->sim.pAtomY                     = gpu->pbAtom->_pDevData + gpu->sim.stride;
    gpu->sim.pAtomXYSP                  = gpu->pbAtomXYSP->_pDevData;
    gpu->sim.pAtomZ                     = gpu->pbAtom->_pDevData + gpu->sim.stride2;
    gpu->sim.pAtomZSP                   = gpu->pbAtomZSP->_pDevData;
    gpu->sim.pAtomSigEps                = gpu->pbAtomSigEps->_pDevData;
    gpu->sim.pAtomRBorn                 = gpu->pbAtomRBorn->_pDevData;
    gpu->sim.pAtomS                     = gpu->pbAtomS->_pDevData;
    gpu->sim.pAtomCharge                = gpu->pbAtomCharge->_pDevData;
    gpu->sim.pAtomChargeSP              = gpu->pbAtomChargeSP->_pDevData;
    gpu->sim.pAtomMass                  = gpu->pbAtomMass->_pDevData;
    gpu->sim.pAtomInvMass               = gpu->pbAtomMass->_pDevData + gpu->sim.stride;
    gpu->sim.pForce                     = gpu->pbForce->_pDevData;
    gpu->sim.pForceX                    = gpu->pbForce->_pDevData;
    gpu->sim.pForceY                    = gpu->pbForce->_pDevData + gpu->sim.stride;
    gpu->sim.pForceZ                    = gpu->pbForce->_pDevData + gpu->sim.stride2;
    if (gpu->ntb == 0)
    {
        gpu->sim.pReff                  = gpu->pbReff->_pDevData;
        gpu->sim.pReffSP                = gpu->pbReffSP->_pDevData;
        gpu->sim.pPsi                   = gpu->pbPsi->_pDevData;
        gpu->sim.pTemp7                 = gpu->pbTemp7->_pDevData;
#ifdef MPI
        gpu->sim.pReffa                 = gpu->pbReffa->_pDevData;
        gpu->sim.pTemp7a                = gpu->pbTemp7a->_pDevData;
    }
    else
    {
        gpu->sim.pPMEForce              = gpu->pbPMEForce->_pDevData;
    }
#else
    }   
#endif    
#ifdef MPI
    gpu->sim.pInForce               = gpu->pbInForce->_pDevData;
    gpu->sim.pOutForce              = gpu->pbOutForce->_pDevData;
#endif

    if (gpu->sim.ntp > 0)
    {
        gpu->sim.pNBForce               = gpu->pbForce->_pDevData + gpu->sim.stride3;
        gpu->sim.pNBForceX              = gpu->pbForce->_pDevData + gpu->sim.stride3;
        gpu->sim.pNBForceY              = gpu->pbForce->_pDevData + gpu->sim.stride * 4;
        gpu->sim.pNBForceZ              = gpu->pbForce->_pDevData + gpu->sim.stride * 5;
        gpu->sim.pUllForce              = (unsigned long long int*)gpu->pbUllForce->_pDevData;
        gpu->sim.pUllForceX             = (unsigned long long int*)gpu->pbUllForce->_pDevData;
        gpu->sim.pUllForceY             = (unsigned long long int*)gpu->pbUllForce->_pDevData + gpu->sim.stride;
        gpu->sim.pUllForceZ             = (unsigned long long int*)gpu->pbUllForce->_pDevData + gpu->sim.stride2;    
    }
    gpu->sim.pVelX                      = gpu->pbVel->_pDevData;
    gpu->sim.pVelY                      = gpu->pbVel->_pDevData + gpu->sim.stride;
    gpu->sim.pVelZ                      = gpu->pbVel->_pDevData + gpu->sim.stride2;
    gpu->sim.pLVelX                     = gpu->pbLVel->_pDevData;
    gpu->sim.pLVelY                     = gpu->pbLVel->_pDevData + gpu->sim.stride;
    gpu->sim.pLVelZ                     = gpu->pbLVel->_pDevData + gpu->sim.stride2;
    gpu->sim.pXMin                      = gpu->pbCenter->_pDevData + 0 * gpu->blocks;
    gpu->sim.pYMin                      = gpu->pbCenter->_pDevData + 1 * gpu->blocks;
    gpu->sim.pZMin                      = gpu->pbCenter->_pDevData + 2 * gpu->blocks;
    gpu->sim.pXMax                      = gpu->pbCenter->_pDevData + 3 * gpu->blocks;
    gpu->sim.pYMax                      = gpu->pbCenter->_pDevData + 4 * gpu->blocks;
    gpu->sim.pZMax                      = gpu->pbCenter->_pDevData + 5 * gpu->blocks;
    gpu->sim.pOutputBufferCounter       = gpu->pbOutputBufferCounter->_pDevData;
    gpu->sim.pRandomPos                 = gpu->pbRandomPos->_pDevData;
    gpu->sim.pRandom                    = gpu->pbRandom->_pDevData;
    gpu->sim.pRandomX                   = gpu->sim.pRandom;
    gpu->sim.pRandomY                   = gpu->sim.pRandomX + gpu->sim.randomSteps * gpu->sim.paddedNumberOfAtoms;
    gpu->sim.pRandomZ                   = gpu->sim.pRandomY + gpu->sim.randomSteps * gpu->sim.paddedNumberOfAtoms;
    
    // Add additional force output for PME reciprocal forces
    if (gpu->bNeighborList)
    {
        int buffers = 1;
        for (int i = 0; i < gpu->sim.atoms; i++)
            gpu->pbOutputBufferCounter->_pSysData[i] = buffers;
    }
    gpuCopyConstants();
}

extern "C" void gpu_upload_crd_(double atm_crd[][3])
{
PRINTMETHOD("gpu_upload_crd");
    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL))
    {
        if (gpu->bNewNeighborList)
        {
            gpu->pbImageIndex->Download();
            gpu->bNewNeighborList                           = false;
        }
        gpu->pbImage->Download();
        unsigned int* pImageAtomLookup                      = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        PMEDouble *pCrd                                     = gpu->pbImage->_pSysData;
        if (gpu->sim.pImageX != gpu->pbImage->_pDevData)
            pCrd                                            = gpu->pbImage->_pSysData + gpu->sim.stride3;
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            int i1                                          = pImageAtomLookup[i];
            pCrd[i1]                                        = atm_crd[i][0];
            pCrd[i1 + gpu->sim.stride]                      = atm_crd[i][1];
            pCrd[i1 + gpu->sim.stride2]                     = atm_crd[i][2];
        }
        gpu->pbImage->Upload();
    }
    else 
    {
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            gpu->pbAtom->_pSysData[i]                       = atm_crd[i][0];
            gpu->pbAtom->_pSysData[i + gpu->sim.stride]     = atm_crd[i][1];
            gpu->pbAtom->_pSysData[i + gpu->sim.stride2]    = atm_crd[i][2];
            gpu->pbAtomXYSP->_pSysData[i].x                 = atm_crd[i][0];
            gpu->pbAtomXYSP->_pSysData[i].y                 = atm_crd[i][1];
            gpu->pbAtomZSP->_pSysData[i]                    = atm_crd[i][2];
        }
        for (int i = gpu->sim.atoms; i < gpu->sim.stride; i++)
        {
            gpu->pbAtom->_pSysData[i]                       = 9999990000.0 + i * 2000.0;
            gpu->pbAtom->_pSysData[i + gpu->sim.stride]     = 9999990000.0 + i * 2000.0;
            gpu->pbAtom->_pSysData[i + gpu->sim.stride2]    = 9999990000.0 + i * 2000.0;
            gpu->pbAtomXYSP->_pSysData[i].x                 = 9999990000.0f + i * 2000.0;
            gpu->pbAtomXYSP->_pSysData[i].y                 = 9999990000.0f + i * 2000.0;
            gpu->pbAtomZSP->_pSysData[i]                    = 9999990000.0f + i * 2000.0;
        }
        gpu->pbAtom->Upload();
        gpu->pbAtomXYSP->Upload();
        gpu->pbAtomZSP->Upload();   
    }
    

}

extern "C" void gpu_download_crd_(double atm_crd[][3])
{
PRINTMETHOD("gpu_download_crd");
    gpu->pbAtom->Download();
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        atm_crd[i][0]                   = gpu->pbAtom->_pSysData[i];
        atm_crd[i][1]                   = gpu->pbAtom->_pSysData[i + gpu->sim.stride];
        atm_crd[i][2]                   = gpu->pbAtom->_pSysData[i + gpu->sim.stride2];
    } 
    
    if (gpu->bNeighborList)
    {
        if (gpu->bNewNeighborList)
        {
            gpu->pbImageIndex->Download();
            gpu->bNewNeighborList                       = false;
        }
        gpu->pbImage->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        PMEDouble *pCrd                                 = gpu->pbImage->_pSysData;
        if (gpu->sim.pImageX != gpu->pbImage->_pDevData)
            pCrd                                        = gpu->pbImage->_pSysData + gpu->sim.stride3;
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            int i1                                      = pImageAtomLookup[i];
            atm_crd[i][0]                               = pCrd[i1];
            atm_crd[i][1]                               = pCrd[i1 + gpu->sim.stride];
            atm_crd[i][2]                               = pCrd[i1 + gpu->sim.stride2];
        }
    }
    else
    {   
        gpu->pbVel->Download();
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            atm_crd[i][0]                               = gpu->pbAtom->_pSysData[i];
            atm_crd[i][1]                               = gpu->pbAtom->_pSysData[i + gpu->sim.stride];
            atm_crd[i][2]                               = gpu->pbAtom->_pSysData[i + gpu->sim.stride2];
        }
    }     
}

extern "C" void gpu_upload_charges_(double charge[])
{
PRINTMETHOD("gpu_upload_charges");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        gpu->pbAtomCharge->_pSysData[i]     = charge[i];
        gpu->pbAtomChargeSP->_pSysData[i]   = charge[i];
    }
    for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++)
    {
        gpu->pbAtomCharge->_pSysData[i]     = (PMEDouble)0.0;
        gpu->pbAtomChargeSP->_pSysData[i]   = (PMEFloat)0.0;
    }
    gpu->pbAtomCharge->Upload();   
    gpu->pbAtomChargeSP->Upload();       
}

extern "C" void gpu_download_charges_(double charge[])
{
PRINTMETHOD("gpu_download_charges");
    gpu->pbAtomCharge->Download();  
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        charge[i]                           = gpu->pbAtomCharge->_pSysData[i];
    } 
}

extern "C" void gpu_upload_sigeps_(double sig[], double eps[])
{
PRINTMETHOD("gpu_upload_sigeps");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        gpu->pbAtomSigEps->_pSysData[i].x   = sig[i];
        gpu->pbAtomSigEps->_pSysData[i].y   = eps[i];
    }
    for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++)
    {
        gpu->pbAtomSigEps->_pSysData[i].x   = (PMEFloat)0.0;
        gpu->pbAtomSigEps->_pSysData[i].y   = (PMEFloat)0.0;
    }
    gpu->pbAtomSigEps->Upload();       
}

extern "C" void gpu_download_sigeps_(double sig[], double eps[])
{
PRINTMETHOD("gpu_download_sigeps");
    gpu->pbAtomSigEps->Download();  
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        sig[i]                              = gpu->pbAtomSigEps->_pSysData[i].x;
        eps[i]                              = gpu->pbAtomSigEps->_pSysData[i].y;
    } 
}

extern "C" void gpu_upload_fs_(double fs[])
{
PRINTMETHOD("gpu_upload_fs");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        gpu->pbAtomS->_pSysData[i]      = fs[i];
    }
    for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++)
    {
        gpu->pbAtomS->_pSysData[i]      = (PMEFloat)0.0;
    }
    gpu->pbAtomS->Upload();       
}

extern "C" void gpu_download_fs_(double fs[])
{
PRINTMETHOD("gpu_download_fs");
    gpu->pbAtomS->Download();  
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        fs[i]                           = gpu->pbAtomS->_pSysData[i];
    } 
}

extern "C" void gpu_upload_rborn_(double rborn[])
{
PRINTMETHOD("gpu_upload_rborn");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        gpu->pbAtomRBorn->_pSysData[i]  = rborn[i];
    }
    for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++)
    {
        gpu->pbAtomRBorn->_pSysData[i]  = (PMEFloat)0.0001;
        gpu->pbReff->_pSysData[i]       = (PMEFloat)0.0001;
        gpu->pbReffSP->_pSysData[i]     = (PMEFloat)0.0001;
        gpu->pbTemp7->_pSysData[i]      = (PMEFloat)0.0;
    }
    gpu->pbAtomRBorn->Upload();  
    gpu->pbReff->Upload();
    gpu->pbReffSP->Upload();     
    gpu->pbTemp7->Upload();
}

extern "C" void gpu_download_rborn_(double rborn[])
{
PRINTMETHOD("gpu_download_rborn");
    gpu->pbAtomRBorn->Download();  
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        rborn[i]                        = gpu->pbAtomRBorn->_pSysData[i];
    } 
}

extern "C" void gpu_upload_masses_(double mass[])
{
PRINTMETHOD("gpu_upload_masses");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        gpu->pbAtomMass->_pSysData[i]                       = mass[i];
        if (mass[i] > 0.0)
            gpu->pbAtomMass->_pSysData[i + gpu->sim.stride] = (PMEDouble)(1.0 / mass[i]);
        else
            gpu->pbAtomMass->_pSysData[i + gpu->sim.stride] = (PMEDouble)0.0;
    }
    for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++)
    {
        gpu->pbAtomMass->_pSysData[i]                       = (PMEDouble)0.0;
        gpu->pbAtomMass->_pSysData[i + gpu->sim.stride]     = (PMEDouble)10000000.0;
    }
    gpu->pbAtomMass->Upload();       
}

extern "C" void gpu_download_masses_(double mass[])
{
PRINTMETHOD("gpu_download_masses");
    gpu->pbAtomMass->Download();  
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        mass[i]                         = gpu->pbAtomMass->_pSysData[i];
    } 
}


extern "C" void gpu_upload_reff_(double reff[])
{
PRINTMETHOD("gpu_upload_reff");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        gpu->pbReff->_pSysData[i]      = reff[i];
    }
    for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++)
    {
        gpu->pbReff->_pSysData[i]      = (PMEFloat)0.0001;
    }
    gpu->pbReff->Upload();       
}

extern "C" void gpu_download_reff_(double reff[])
{
PRINTMETHOD("gpu_download_reff");
    gpu->pbReff->Download();  
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        reff[i]                        = gpu->pbReff->_pSysData[i];
    } 
}

extern "C" void gpu_upload_frc_(double atm_frc[][3])
{
PRINTMETHOD("gpu_upload_frc");
#ifdef MPI
    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL))
    {
        if (gpu->bNewNeighborList)
        {
            gpu->pbImageIndex->Download();
            gpu->bNewNeighborList                       = false;
        }
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        PMEDouble *pForce                               = gpu->pbInForce->_pSysData;
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            int i1                                      = pImageAtomLookup[i];
            pForce[i1 * 3]                              = atm_frc[i][0];
            pForce[i1 * 3 + 1]                          = atm_frc[i][1];
            pForce[i1 * 3 + 2]                          = atm_frc[i][2];
        }
        gpu->pbInForce->Upload();
    }
    else
    {
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            gpu->pbInForce->_pSysData[i * 3]                = atm_frc[i][0];
            gpu->pbInForce->_pSysData[i * 3 + 1]            = atm_frc[i][1];
            gpu->pbInForce->_pSysData[i * 3 + 2]            = atm_frc[i][2];
        }
    }
#else
    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL))
    {
        if (gpu->bNewNeighborList)
        {
            gpu->pbImageIndex->Download();
            gpu->bNewNeighborList                       = false;
        }
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        PMEDouble *pForce                               = gpu->pbForce->_pSysData;
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            int i1                                      = pImageAtomLookup[i];
            pForce[i1]                                  = atm_frc[i][0];
            pForce[i1 + gpu->sim.stride]                = atm_frc[i][1];
            pForce[i1 + gpu->sim.stride2]               = atm_frc[i][2];
        }
        gpu->pbForce->Upload();
    }
    else
    {
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            gpu->pbForce->_pSysData[i]                      = atm_frc[i][0];
            gpu->pbForce->_pSysData[i + gpu->sim.stride]    = atm_frc[i][1];
            gpu->pbForce->_pSysData[i + gpu->sim.stride2]   = atm_frc[i][2];
        }
        gpu->pbForce->Upload();
    }
#endif
}

extern "C" void gpu_download_frc_(double atm_frc[][3])
{
PRINTMETHOD("gpu_download_frc");
#ifdef MPI
    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL))
    {
        if (gpu->bNewNeighborList)
        {
            gpu->pbImageIndex->Download();
            gpu->bNewNeighborList                       = false;
        }
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        PMEDouble *pForce                               = gpu->pbInForce->_pSysData;
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            int i1                                      = pImageAtomLookup[i];
            atm_frc[i][0]                               = pForce[i1 * 3];
            atm_frc[i][1]                               = pForce[i1 * 3 + 1];
            atm_frc[i][2]                               = pForce[i1 * 3 + 2];
        }
    }
    else
    {
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            atm_frc[i][0]                               = gpu->pbInForce->_pSysData[i * 3];
            atm_frc[i][1]                               = gpu->pbInForce->_pSysData[i * 3 + 1];
            atm_frc[i][2]                               = gpu->pbInForce->_pSysData[i * 3 + 2];
        } 
    }
#else
    gpu->pbForce->Download();
    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL))
    {
        if (gpu->bNewNeighborList)
        {
            gpu->pbImageIndex->Download();
            gpu->bNewNeighborList                       = false;
        }
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        PMEDouble *pForce                               = gpu->pbForce->_pSysData;
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            int i1                                      = pImageAtomLookup[i];
            atm_frc[i][0]                               = pForce[i1];
            atm_frc[i][1]                               = pForce[i1 + gpu->sim.stride];
            atm_frc[i][2]                               = pForce[i1 + gpu->sim.stride2];
        }
    }
    else
    {
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            atm_frc[i][0]                               = gpu->pbForce->_pSysData[i];
            atm_frc[i][1]                               = gpu->pbForce->_pSysData[i + gpu->sim.stride];
            atm_frc[i][2]                               = gpu->pbForce->_pSysData[i + gpu->sim.stride2];
        } 
    }
#endif
}

extern "C" void gpu_upload_vel_(double atm_vel[][3])
{
PRINTMETHOD("gpu_upload_vel");
    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL))
    {
        if (gpu->bNewNeighborList)
        {
            gpu->pbImageIndex->Download();
            gpu->bNewNeighborList                       = false;
        }
        gpu->pbImageVel->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        PMEDouble *pVel                                 = gpu->pbImageVel->_pSysData;
        if (gpu->sim.pImageVelX != gpu->pbImageVel->_pDevData)
            pVel                                        = gpu->pbImageVel->_pSysData + gpu->sim.stride3;
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            int i1                                      = pImageAtomLookup[i];
            pVel[i1]                                    = atm_vel[i][0];
            pVel[i1 + gpu->sim.stride]                  = atm_vel[i][1];
            pVel[i1 + gpu->sim.stride2]                 = atm_vel[i][2];
        }
        gpu->pbImageVel->Upload();
    }
    else
    {
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            gpu->pbVel->_pSysData[i]                    = atm_vel[i][0];
            gpu->pbVel->_pSysData[i + gpu->sim.stride]  = atm_vel[i][1];
            gpu->pbVel->_pSysData[i + gpu->sim.stride2] = atm_vel[i][2];
        }
        gpu->pbVel->Upload();
    } 
}

extern "C" void gpu_download_vel_(double atm_vel[][3])
{
PRINTMETHOD("gpu_download_vel");
    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL))
    {
        if (gpu->bNewNeighborList)
        {
            gpu->pbImageIndex->Download();
            gpu->bNewNeighborList                       = false;
        }
        gpu->pbImageVel->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        PMEDouble *pVel                                 = gpu->pbImageVel->_pSysData;
        if (gpu->sim.pImageVelX != gpu->pbImageVel->_pDevData)
            pVel                                        = gpu->pbImageVel->_pSysData + gpu->sim.stride3;
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            int i1                                      = pImageAtomLookup[i];
            atm_vel[i][0]                               = pVel[i1];
            atm_vel[i][1]                               = pVel[i1 + gpu->sim.stride];
            atm_vel[i][2]                               = pVel[i1 + gpu->sim.stride2];
        }
    }
    else
    {   
        gpu->pbVel->Download();
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            atm_vel[i][0]                               = gpu->pbVel->_pSysData[i];
            atm_vel[i][1]                               = gpu->pbVel->_pSysData[i + gpu->sim.stride];
            atm_vel[i][2]                               = gpu->pbVel->_pSysData[i + gpu->sim.stride2];
        }
    } 
}

extern "C" void gpu_upload_last_vel_(double atm_last_vel[][3])
{
PRINTMETHOD("gpu_upload_last_vel");
    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL))
    {
        if (gpu->bNewNeighborList)
        {
            gpu->pbImageIndex->Download();
            gpu->bNewNeighborList                       = false;
        }
        gpu->pbImageLVel->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        PMEDouble *pLVel                                = gpu->pbImageLVel->_pSysData;
        if (gpu->sim.pImageLVelX != gpu->pbImageLVel->_pDevData)
            pLVel                                       = gpu->pbImageLVel->_pSysData + gpu->sim.stride3;
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            int i1                                      = pImageAtomLookup[i];
            pLVel[i1]                                   = atm_last_vel[i][0];
            pLVel[i1 + gpu->sim.stride]                 = atm_last_vel[i][1];
            pLVel[i1 + gpu->sim.stride2]                = atm_last_vel[i][2];
        }
        gpu->pbImageLVel->Upload();
    }
    else
    {
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            gpu->pbLVel->_pSysData[i]                   = atm_last_vel[i][0];
            gpu->pbLVel->_pSysData[i + gpu->sim.stride] = atm_last_vel[i][1];
            gpu->pbLVel->_pSysData[i + gpu->sim.stride2]= atm_last_vel[i][2];
        }  
        gpu->pbLVel->Upload();
    } 
}

extern "C" void gpu_download_last_vel_(double atm_last_vel[][3])
{
PRINTMETHOD("gpu_download_last_vel");
    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL))
    {
        if (gpu->bNewNeighborList)
        {
            gpu->pbImageIndex->Download();
            gpu->bNewNeighborList                       = false;
        }
        gpu->pbImageLVel->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        PMEDouble *pLVel                                = gpu->pbImageLVel->_pSysData;
        if (gpu->sim.pImageLVelX != gpu->pbImageLVel->_pDevData)
            pLVel                                       = gpu->pbImageLVel->_pSysData + gpu->sim.stride3;
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            int i1                                      = pImageAtomLookup[i];
            atm_last_vel[i][0]                          = pLVel[i1];
            atm_last_vel[i][1]                          = pLVel[i1 + gpu->sim.stride];
            atm_last_vel[i][2]                          = pLVel[i1 + gpu->sim.stride2];
        }
    }
    else
    {
        gpu->pbLVel->Download();
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            atm_last_vel[i][0]                          = gpu->pbLVel->_pSysData[i];
            atm_last_vel[i][1]                          = gpu->pbLVel->_pSysData[i + gpu->sim.stride];
            atm_last_vel[i][2]                          = gpu->pbLVel->_pSysData[i + gpu->sim.stride2];
        }
    } 
}

extern "C" void gpu_bonds_setup_(int* cit_nbona, bond_rec cit_a_bond[], int* cit_nbonh, bond_rec cit_h_bond[], double gbl_req[], double gbl_rk[])
{
PRINTMETHOD("gpu_bonds_setup");    
    // Count non-zero bonds
    int bonds                                               = 0;
    if (gpu->ntf < 3)
    {
        for (int i = 0; i < *cit_nbona; i++)
            if (gbl_rk[cit_a_bond[i].parm_idx - 1])
                bonds++;
    }
    if (gpu->ntf < 2)
    {
        for (int i = 0; i < *cit_nbonh; i++)
            if (gbl_rk[cit_h_bond[i].parm_idx - 1])
                bonds++;
    }
    
    // Allocate/reallocate GPU bond buffers
    delete gpu->pbBond;
    delete gpu->pbBondID;                

#ifdef GVERBOSE
    printf("%d bonds, %d active\n", *cit_nbona + *cit_nbonh, bonds); 
#endif    
    gpu->pbBond                                             = new GpuBuffer<PMEDouble2>(bonds);
    gpu->pbBondID                                           = new GpuBuffer<int4>(bonds);
    
    bonds                                                   = 0;
    if (gpu->ntf < 3)
    {
        for (int i = 0; i < *cit_nbona; i++)
        {
            if (gbl_rk[cit_a_bond[i].parm_idx - 1])
            {          
                gpu->pbBond->_pSysData[bonds].x             = gbl_rk[cit_a_bond[i].parm_idx - 1];
                gpu->pbBond->_pSysData[bonds].y             = gbl_req[cit_a_bond[i].parm_idx - 1];
                gpu->pbBondID->_pSysData[bonds].x           = abs(cit_a_bond[i].atm_i) - 1;
                gpu->pbBondID->_pSysData[bonds].y           = abs(cit_a_bond[i].atm_j) - 1;
                gpu->pbBondID->_pSysData[bonds].z           = abs(cit_a_bond[i].atm_i) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cit_a_bond[i].atm_i) - 1]++;
                gpu->pbBondID->_pSysData[bonds].w           = abs(cit_a_bond[i].atm_j) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cit_a_bond[i].atm_j) - 1]++;
                bonds++;                
            }
        }
    }
    if (gpu->ntf < 2)
    {
        for (int i = 0; i < *cit_nbonh; i++)
        {
            if (gbl_rk[cit_h_bond[i].parm_idx - 1])
            {        
                gpu->pbBond->_pSysData[bonds].x             = gbl_rk[cit_h_bond[i].parm_idx - 1];
                gpu->pbBond->_pSysData[bonds].y             = gbl_req[cit_h_bond[i].parm_idx - 1];
                gpu->pbBondID->_pSysData[bonds].x           = abs(cit_h_bond[i].atm_i) - 1;
                gpu->pbBondID->_pSysData[bonds].y           = abs(cit_h_bond[i].atm_j) - 1;
                gpu->pbBondID->_pSysData[bonds].z           = abs(cit_h_bond[i].atm_i) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cit_h_bond[i].atm_i) - 1]++;
                gpu->pbBondID->_pSysData[bonds].w           = abs(cit_h_bond[i].atm_j) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cit_h_bond[i].atm_j) - 1]++;
                bonds++;             
            }
        }
    }
    gpu->pbBond->Upload();
    gpu->pbBondID->Upload();

    // Set constants
    gpu->sim.bonds                                          = bonds;
    gpu->sim.pBond                                          = gpu->pbBond->_pDevData;
    gpu->sim.pBondID                                        = gpu->pbBondID->_pDevData;
    gpuCopyConstants();
}

extern "C" void gpu_angles_setup_(int* angle_cnt, int* ntheth, angle_rec cit_angle[], double gbl_teq[], double gbl_tk[])
{
PRINTMETHOD("gpu_angles_setup");

    // Allocate/reallocate GPU bond angle buffers
    delete gpu->pbBondAngle;
    delete gpu->pbBondAngleID1;
    delete gpu->pbBondAngleID2;      
    int angles                                              = 0;  
    if (gpu->ntf < 4)
        angles                                              = *angle_cnt;
    else if (gpu->ntf < 5)
        angles                                              = *angle_cnt - *ntheth;
     
#ifdef GVERBOSE
    printf("%d bond angles, %d active\n", *angle_cnt, angles);
#endif        
    gpu->pbBondAngle                                        = new GpuBuffer<PMEDouble2>(angles);
    gpu->pbBondAngleID1                                     = new GpuBuffer<int4>(angles);
    gpu->pbBondAngleID2                                     = new GpuBuffer<int2>(angles);

    // Copy bond angles
    angles = 0;
    if (gpu->ntf < 5)
    {
        for (int i = 0; i < *angle_cnt; i++)
        {
            if ((gpu->ntf < 4) || (i >= *ntheth))
            { 
                gpu->pbBondAngle->_pSysData[angles].x       = gbl_tk[cit_angle[i].parm_idx -1];
                gpu->pbBondAngle->_pSysData[angles].y       = gbl_teq[cit_angle[i].parm_idx - 1];
                gpu->pbBondAngleID1->_pSysData[angles].x    = abs(cit_angle[i].atm_i) - 1;
                gpu->pbBondAngleID1->_pSysData[angles].y    = abs(cit_angle[i].atm_j) - 1;
                gpu->pbBondAngleID1->_pSysData[angles].z    = abs(cit_angle[i].atm_k) - 1;
                gpu->pbBondAngleID1->_pSysData[angles].w    = abs(cit_angle[i].atm_i) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cit_angle[i].atm_i) - 1]++;
                gpu->pbBondAngleID2->_pSysData[angles].x    = abs(cit_angle[i].atm_j) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cit_angle[i].atm_j) - 1]++;
                gpu->pbBondAngleID2->_pSysData[angles].y    = abs(cit_angle[i].atm_k) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cit_angle[i].atm_k) - 1]++;
                angles++;
            }
        }
    }
    gpu->pbBondAngle->Upload();
    gpu->pbBondAngleID1->Upload();
    gpu->pbBondAngleID2->Upload();

    // Set constants
    gpu->sim.bondAngles                                     = angles;
    gpu->sim.pBondAngle                                     = gpu->pbBondAngle->_pDevData;
    gpu->sim.pBondAngleID1                                  = gpu->pbBondAngleID1->_pDevData;
    gpu->sim.pBondAngleID2                                  = gpu->pbBondAngleID2->_pDevData;
    gpuCopyConstants();
}

static PMEDouble gmul[] = {0.0, 2.0, 0.0, 4.0, 0.0, 6.0, 0.0, 8.0, 0.0, 10.0};

extern "C" void gpu_dihedrals_setup_(int* dihed_cnt, int* nphih, dihed_rec cit_dihed[], int gbl_ipn[], double gbl_pn[], double gbl_pk[], double gbl_gamc[], double gbl_gams[])
{
PRINTMETHOD("gpu_dihedrals_setup");

    // Allocate/reallocate GPU dihedral buffers
    delete gpu->pbDihedral1;
    delete gpu->pbDihedral2;
    delete gpu->pbDihedral3;
    delete gpu->pbDihedralID1;
    delete gpu->pbDihedralID2;   
    int dihedrals                                           = 0;  
    if (gpu->ntf < 6)
        dihedrals                                           = *dihed_cnt;
    else if (gpu->ntf < 7)
        dihedrals                                           = *dihed_cnt - *nphih;
                       
#ifdef GVERBOSE
    printf("%d dihedrals, %d active\n", *dihed_cnt, dihedrals);
#endif                
    gpu->pbDihedral1                                        = new GpuBuffer<PMEDouble2>(dihedrals);
    gpu->pbDihedral2                                        = new GpuBuffer<PMEDouble2>(dihedrals);
    gpu->pbDihedral3                                        = new GpuBuffer<PMEDouble>(dihedrals);
    gpu->pbDihedralID1                                      = new GpuBuffer<int4>(dihedrals);
    gpu->pbDihedralID2                                      = new GpuBuffer<int4>(dihedrals);

    // Copy dihedrals
    dihedrals                                               = 0;
    if (gpu->ntf < 7)
    {
        for (int i = 0; i < *dihed_cnt; i++)
        {
            if ((gpu->ntf < 6) || (i >= *nphih))
            {
                gpu->pbDihedral1->_pSysData[dihedrals].x    = gmul[gbl_ipn[cit_dihed[i].parm_idx -1]];
                gpu->pbDihedral1->_pSysData[dihedrals].y    = gbl_pn[cit_dihed[i].parm_idx - 1];
                gpu->pbDihedral2->_pSysData[dihedrals].x    = gbl_pk[cit_dihed[i].parm_idx - 1];
                gpu->pbDihedral2->_pSysData[dihedrals].y    = gbl_gamc[cit_dihed[i].parm_idx - 1];
                gpu->pbDihedral3->_pSysData[dihedrals]      = gbl_gams[cit_dihed[i].parm_idx - 1];
                gpu->pbDihedralID1->_pSysData[dihedrals].x  = abs(cit_dihed[i].atm_i) - 1;
                gpu->pbDihedralID1->_pSysData[dihedrals].y  = abs(cit_dihed[i].atm_j) - 1;
                gpu->pbDihedralID1->_pSysData[dihedrals].z  = abs(cit_dihed[i].atm_k) - 1;
                gpu->pbDihedralID1->_pSysData[dihedrals].w  = abs(cit_dihed[i].atm_l) - 1;
                gpu->pbDihedralID2->_pSysData[dihedrals].x  = abs(cit_dihed[i].atm_i) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cit_dihed[i].atm_i) - 1]++;
                gpu->pbDihedralID2->_pSysData[dihedrals].y  = abs(cit_dihed[i].atm_j) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cit_dihed[i].atm_j) - 1]++;
                gpu->pbDihedralID2->_pSysData[dihedrals].z  = abs(cit_dihed[i].atm_k) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cit_dihed[i].atm_k) - 1]++;
                gpu->pbDihedralID2->_pSysData[dihedrals].w  = abs(cit_dihed[i].atm_l) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cit_dihed[i].atm_l) - 1]++;
                dihedrals++;
            }
        }
    }

    gpu->pbDihedral1->Upload();
    gpu->pbDihedral2->Upload();
    gpu->pbDihedral3->Upload();   
    gpu->pbDihedralID1->Upload();   
    gpu->pbDihedralID2->Upload();
 
    // Set constants
    gpu->sim.dihedrals                                      = dihedrals;
    gpu->sim.pDihedral1                                     = gpu->pbDihedral1->_pDevData;
    gpu->sim.pDihedral2                                     = gpu->pbDihedral2->_pDevData;
    gpu->sim.pDihedral3                                     = gpu->pbDihedral3->_pDevData;
    gpu->sim.pDihedralID1                                   = gpu->pbDihedralID1->_pDevData;
    gpu->sim.pDihedralID2                                   = gpu->pbDihedralID2->_pDevData;

    gpuCopyConstants();
}

extern "C" void gpu_nb14_setup_(int* nb14_cnt, int cit_nb14[][3], double gbl_one_scee[], double gbl_one_scnb[], int* ntypes, int iac[], int ico[], double cn1[], double cn2[], double cn114[], double cn214[])
{
PRINTMETHOD("gpu_nb14_setup");

    // Allocate/reallocate GPU nb14 buffers
    int nb14s                                       = 0;
    if (gpu->ntf < 8)
        nb14s                                       = *nb14_cnt;
    else
        nb14s                                       = 0;
    delete gpu->pbNb141;
    delete gpu->pbNb142;
    delete gpu->pbNb14ID;        

#ifdef GVERBOSE
    printf("%d 1-4 nonbonds, %d active\n", *nb14_cnt, nb14s);
#endif
    gpu->pbNb141                                    = new GpuBuffer<PMEDouble2>(nb14s);
    gpu->pbNb142                                    = new GpuBuffer<PMEDouble2>(nb14s);
    gpu->pbNb14ID                                   = new GpuBuffer<int4>(nb14s);
    
#ifdef GVERBOSE
    printf("%d atom types\n", *ntypes);
#endif    
    
    // Re-order output buffers to put 1-4s at beginning if necessary to calculate molecular virial
    if (gpu->sim.ntp > 0)
    {
        int* pNB14s                                 = new int[gpu->sim.atoms];
        for (int i = 0; i < gpu->sim.atoms; i++)        
            pNB14s[i]                               = 0;

        int maxNonbonds                             = 0;

        for (int i = 0; i < nb14s; i++)
        {
            pNB14s[abs(cit_nb14[i][0]) - 1]++;
            pNB14s[abs(cit_nb14[i][1]) - 1]++;
        }
        for (int i = 0; i < gpu->sim.atoms; i++)
            if (maxNonbonds < pNB14s[i])
                maxNonbonds                          = pNB14s[i];
        delete[] pNB14s;
        maxNonbonds++;

#ifdef GVERBOSE
        printf("Max nb14s = %d\n", maxNonbonds);
#endif
        
        gpu->sim.maxNonbonds                        = maxNonbonds;
    }    
    
    // Generate sigma and epsilon from atom types
    double* sigma                                   = new double[*ntypes];
    double* epsilon                                 = new double[*ntypes];
    for (int i = 0; i < *ntypes; i++)
    {
        int nbtype                                  = ico[*ntypes * i + i] - 1;
        if (nbtype >= 0)
        {
            double c1                               = cn1[nbtype];
            double c2                               = cn2[nbtype];
            double sig                              = pow(c1 / c2, 1.0 / 6.0);
            double eps                              = (c2 * c2) / (4.0 * c1);
            sigma[i]                                = sig / 2.0;
            epsilon[i]                              = sqrt(eps);
        }
        else
        {
            sigma[i]                                = 0.5;
            epsilon[i]                              = 0.0;
        }
    }
#if 0
    printf("%d nonbond types\n", *ntypes);
    exit(-1);
#endif
    int* pOutputBufferCounter = NULL;
    if (gpu->sim.ntp > 0)
    {
        pOutputBufferCounter                        = new int[gpu->sim.atoms];
        for (int i = 0; i < gpu->sim.atoms; i++)
            pOutputBufferCounter[i]                 = 1; 
    }
    
    // Copy atomic nonbond parameters
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        gpu->pbAtomSigEps->_pSysData[i].x           = sigma[iac[i] - 1];
        gpu->pbAtomSigEps->_pSysData[i].y           = (PMEFloat)(2.0 * epsilon[iac[i] - 1]);    
    }
    for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++)
    {
        gpu->pbAtomSigEps->_pSysData[i].x           = (PMEFloat)0.0;
        gpu->pbAtomSigEps->_pSysData[i].y           = (PMEFloat)0.0;
    }
    gpu->pbAtomSigEps->Upload();
    
    // Copy 1-4 interactions
    if (gpu->ntf < 8)
    {
        for (int i = 0; i < nb14s; i++)
        {
            int parm_idx                            = cit_nb14[i][2] - 1;
         //   printf("%5d: %5d - %5d  %5d %9.4f %9.4f\n", i, cit_nb14[i][0], cit_nb14[i][1], gbl_one_scee[parm_idx], gbl_one_scnb[parm_idx]);
            gpu->pbNb141->_pSysData[i].x            = gbl_one_scee[parm_idx] * gpu->pbAtomCharge->_pSysData[abs(cit_nb14[i][0]) - 1] * gpu->pbAtomCharge->_pSysData[abs(cit_nb14[i][1]) - 1];
            gpu->pbNb141->_pSysData[i].y            = gbl_one_scnb[parm_idx];
            int nbtype                              = ico[*ntypes * (iac[abs(cit_nb14[i][0]) - 1] - 1) + (iac[abs(cit_nb14[i][1]) - 1] - 1)] - 1;
            if (nbtype >= 0)
            {
                gpu->pbNb142->_pSysData[i].x        = cn114[nbtype];
                gpu->pbNb142->_pSysData[i].y        = cn214[nbtype];   
            }
            else
            {
                gpu->pbNb142->_pSysData[i].x        = (PMEDouble)0.0;
                gpu->pbNb142->_pSysData[i].y        = (PMEDouble)0.0;   
            }  
            gpu->pbNb14ID->_pSysData[i].x           = abs(cit_nb14[i][0]) - 1;
            gpu->pbNb14ID->_pSysData[i].y           = abs(cit_nb14[i][1]) - 1;
            
            // Special case constant pressure simulation
            if (gpu->sim.ntp > 0)
            {
                gpu->pbNb14ID->_pSysData[i].z       = abs(cit_nb14[i][0]) - 1 + gpu->sim.stride3 * pOutputBufferCounter[abs(cit_nb14[i][0]) - 1]++;
                gpu->pbNb14ID->_pSysData[i].w       = abs(cit_nb14[i][1]) - 1 + gpu->sim.stride3 * pOutputBufferCounter[abs(cit_nb14[i][1]) - 1]++;
            }
            else
            {
                gpu->pbNb14ID->_pSysData[i].z       = abs(cit_nb14[i][0]) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cit_nb14[i][0]) - 1]++;
                gpu->pbNb14ID->_pSysData[i].w       = abs(cit_nb14[i][1]) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cit_nb14[i][1]) - 1]++;
            }
        }
    }
    delete[] pOutputBufferCounter;
    gpu->pbNb141->Upload();
    gpu->pbNb142->Upload();
    gpu->pbNb14ID->Upload();
 
    
    // Set constants
    gpu->sim.nb14s                                  = nb14s;
    gpu->sim.pNb141                                 = gpu->pbNb141->_pDevData;
    gpu->sim.pNb142                                 = gpu->pbNb142->_pDevData;
    gpu->sim.pNb14ID                                = gpu->pbNb14ID->_pDevData;
    gpuCopyConstants();

    
    // Delete temporary arrays
    delete[] sigma;
    delete[] epsilon;
}

extern "C" void gpu_init_extra_pnts_nb14_(int* frame_cnt, ep_frame_rec* ep_frames, double ep_lcl_crd[][2][3])
{
PRINTMETHOD("gpu_init_extra_pnts_nb14");
    int EP11s           = 0;
    int EP12s           = 0;
    int EP21s           = 0;
    int EP22s           = 0;    

    // Count each type of extra point frame
    for (int i = 0; i < *frame_cnt; i++)
    {
        if (ep_frames[i].ep_cnt == 1)
        {
            if (ep_frames[i].type == 1)
                EP11s++;
            else if (ep_frames[i].type == 2)
                EP12s++;
        }
        else if (ep_frames[i].ep_cnt == 2)
        {
            if (ep_frames[i].type == 1)
                EP21s++;
            else if (ep_frames[i].type == 2)
                EP22s++;        
        }
    }
    int EP11Stride                           = ((EP11s + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
    int EP12Stride                           = ((EP12s + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
    int EP21Stride                           = ((EP21s + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
    int EP22Stride                           = ((EP22s + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;

    // Delete existing extra point frames
    delete gpu->pbExtraPoint11Frame;
    delete gpu->pbExtraPoint11Index;
    delete gpu->pbExtraPoint11;
    delete gpu->pbExtraPoint12Frame;
    delete gpu->pbExtraPoint12Index;
    delete gpu->pbExtraPoint12;
    delete gpu->pbExtraPoint21Frame;
    delete gpu->pbExtraPoint21Index;
    delete gpu->pbExtraPoint21;  
    delete gpu->pbExtraPoint22Frame;
    delete gpu->pbExtraPoint22Index;
    delete gpu->pbExtraPoint22;   
    
    // Allocate extra point frames
    gpu->pbExtraPoint11Frame                = new GpuBuffer<int4>(EP11s);
    gpu->pbExtraPoint11Index                = new GpuBuffer<int>(EP11s);
    gpu->pbExtraPoint11                     = new GpuBuffer<PMEDouble>(3 * EP11Stride);
    gpu->pbExtraPoint12Frame                = new GpuBuffer<int4>(EP12s);
    gpu->pbExtraPoint12Index                = new GpuBuffer<int>(EP12s);
    gpu->pbExtraPoint12                     = new GpuBuffer<PMEDouble>(3 * EP12Stride);
    gpu->pbExtraPoint21Frame                = new GpuBuffer<int4>(EP21s);
    gpu->pbExtraPoint21Index                = new GpuBuffer<int2>(EP21s);
    gpu->pbExtraPoint21                     = new GpuBuffer<PMEDouble>(6 * EP21Stride);   
    gpu->pbExtraPoint22Frame                = new GpuBuffer<int4>(EP22s);
    gpu->pbExtraPoint22Index                = new GpuBuffer<int2>(EP22s);
    gpu->pbExtraPoint22                     = new GpuBuffer<PMEDouble>(6 * EP22Stride);    
    
    // Copy Extra point data
    EP11s                                                               = 0;
    EP12s                                                               = 0;
    EP21s                                                               = 0;
    EP22s                                                               = 0;
    for (int i =0; i < *frame_cnt; i++)
    {
        if (ep_frames[i].ep_cnt == 1)
        {
            if (ep_frames[i].type == 1)
            {
                gpu->pbExtraPoint11Index->_pSysData[EP11s]              = ep_frames[i].extra_pnt[0] - 1;
                gpu->pbExtraPoint11Frame->_pSysData[EP11s].x            = ep_frames[i].parent_atm - 1;
                gpu->pbExtraPoint11Frame->_pSysData[EP11s].y            = ep_frames[i].frame_atm1 - 1;
                gpu->pbExtraPoint11Frame->_pSysData[EP11s].z            = ep_frames[i].frame_atm2 - 1;
                gpu->pbExtraPoint11Frame->_pSysData[EP11s].w            = ep_frames[i].frame_atm3 - 1;
                gpu->pbExtraPoint11->_pSysData[EP11s]                   = ep_lcl_crd[i][0][0];
                gpu->pbExtraPoint11->_pSysData[EP11s + EP11Stride]      = ep_lcl_crd[i][0][1];
                gpu->pbExtraPoint11->_pSysData[EP11s + EP11Stride * 2]  = ep_lcl_crd[i][0][2];
                EP11s++;
            }
            else if (ep_frames[i].type == 2)
            {
                gpu->pbExtraPoint12Index->_pSysData[EP12s]              = ep_frames[i].extra_pnt[0] - 1;
                gpu->pbExtraPoint12Frame->_pSysData[EP12s].x            = ep_frames[i].parent_atm - 1;
                gpu->pbExtraPoint12Frame->_pSysData[EP12s].y            = ep_frames[i].frame_atm1 - 1;
                gpu->pbExtraPoint12Frame->_pSysData[EP12s].z            = ep_frames[i].frame_atm2 - 1;
                gpu->pbExtraPoint12Frame->_pSysData[EP12s].w            = ep_frames[i].frame_atm3 - 1;
                gpu->pbExtraPoint12->_pSysData[EP12s]                   = ep_lcl_crd[i][0][0];
                gpu->pbExtraPoint12->_pSysData[EP12s + EP12Stride]      = ep_lcl_crd[i][0][1];
                gpu->pbExtraPoint12->_pSysData[EP12s + EP12Stride * 2]  = ep_lcl_crd[i][0][2];                 
                EP12s++;               
            }
        }
        else if (ep_frames[i].ep_cnt == 2)
        {
            if (ep_frames[i].type == 1)
            {
                gpu->pbExtraPoint21Index->_pSysData[EP21s].x            = ep_frames[i].extra_pnt[0] - 1;
                gpu->pbExtraPoint21Index->_pSysData[EP21s].y            = ep_frames[i].extra_pnt[1] - 1;
                gpu->pbExtraPoint21Frame->_pSysData[EP21s].x            = ep_frames[i].parent_atm - 1;
                gpu->pbExtraPoint21Frame->_pSysData[EP21s].y            = ep_frames[i].frame_atm1 - 1;
                gpu->pbExtraPoint21Frame->_pSysData[EP21s].z            = ep_frames[i].frame_atm2 - 1;
                gpu->pbExtraPoint21Frame->_pSysData[EP21s].w            = ep_frames[i].frame_atm3 - 1;
                gpu->pbExtraPoint21->_pSysData[EP21s]                   = ep_lcl_crd[i][0][0];
                gpu->pbExtraPoint21->_pSysData[EP21s + EP21Stride]      = ep_lcl_crd[i][0][1];
                gpu->pbExtraPoint21->_pSysData[EP21s + EP21Stride * 2]  = ep_lcl_crd[i][0][2];    
                gpu->pbExtraPoint21->_pSysData[EP21s + EP21Stride * 3]  = ep_lcl_crd[i][1][0];
                gpu->pbExtraPoint21->_pSysData[EP21s + EP21Stride * 4]  = ep_lcl_crd[i][1][1];
                gpu->pbExtraPoint21->_pSysData[EP21s + EP21Stride * 5]  = ep_lcl_crd[i][1][2];            
                EP21s++;
            }
            else if (ep_frames[i].type == 2)
            {
                gpu->pbExtraPoint22Index->_pSysData[EP22s].x            = ep_frames[i].extra_pnt[0] - 1;
                gpu->pbExtraPoint22Index->_pSysData[EP22s].y            = ep_frames[i].extra_pnt[1] - 1;
                gpu->pbExtraPoint22Frame->_pSysData[EP22s].x            = ep_frames[i].parent_atm - 1;
                gpu->pbExtraPoint22Frame->_pSysData[EP22s].y            = ep_frames[i].frame_atm1 - 1;
                gpu->pbExtraPoint22Frame->_pSysData[EP22s].z            = ep_frames[i].frame_atm2 - 1;
                gpu->pbExtraPoint22Frame->_pSysData[EP22s].w            = ep_frames[i].frame_atm3 - 1;
                gpu->pbExtraPoint22->_pSysData[EP22s]                   = ep_lcl_crd[i][0][0];
                gpu->pbExtraPoint22->_pSysData[EP22s + EP22Stride]      = ep_lcl_crd[i][0][1];
                gpu->pbExtraPoint22->_pSysData[EP22s + EP22Stride * 2]  = ep_lcl_crd[i][0][2];    
                gpu->pbExtraPoint22->_pSysData[EP22s + EP22Stride * 3]  = ep_lcl_crd[i][1][0];
                gpu->pbExtraPoint22->_pSysData[EP22s + EP22Stride * 4]  = ep_lcl_crd[i][1][1];
                gpu->pbExtraPoint22->_pSysData[EP22s + EP22Stride * 5]  = ep_lcl_crd[i][1][2];                        
            
                EP22s++;        
            }
        }    
    }
        
    // Upload constants and data
    gpu->pbExtraPoint11Frame->Upload();
    gpu->pbExtraPoint11Index->Upload();
    gpu->pbExtraPoint11->Upload();
    gpu->pbExtraPoint12Frame->Upload();
    gpu->pbExtraPoint12Index->Upload();
    gpu->pbExtraPoint12->Upload();
    gpu->pbExtraPoint21Frame->Upload();
    gpu->pbExtraPoint21Index->Upload();
    gpu->pbExtraPoint21->Upload();  
    gpu->pbExtraPoint22Frame->Upload();
    gpu->pbExtraPoint22Index->Upload();
    gpu->pbExtraPoint22->Upload();   
    gpu->sim.EPs                            = *frame_cnt;
    gpu->sim.EP11s                          = EP11s;
    gpu->sim.EP12s                          = EP12s;
    gpu->sim.EP21s                          = EP21s;
    gpu->sim.EP22s                          = EP22s;
    gpu->sim.EP11Offset                     =                       EP11Stride;
    gpu->sim.EP12Offset                     = gpu->sim.EP11Offset + EP12Stride;
    gpu->sim.EP21Offset                     = gpu->sim.EP12Offset + EP21Stride;
    gpu->sim.EP22Offset                     = gpu->sim.EP21Offset + EP22Stride;
    gpu->sim.pExtraPoint11Frame             = gpu->pbExtraPoint11Frame->_pDevData;
    gpu->sim.pExtraPoint11Index             = gpu->pbExtraPoint11Index->_pDevData;
    gpu->sim.pExtraPoint11X                 = gpu->pbExtraPoint11->_pDevData;
    gpu->sim.pExtraPoint11Y                 = gpu->pbExtraPoint11->_pDevData + EP11Stride;
    gpu->sim.pExtraPoint11Z                 = gpu->pbExtraPoint11->_pDevData + EP11Stride * 2;
    gpu->sim.pExtraPoint12Frame             = gpu->pbExtraPoint12Frame->_pDevData;
    gpu->sim.pExtraPoint12Index             = gpu->pbExtraPoint12Index->_pDevData;
    gpu->sim.pExtraPoint12X                 = gpu->pbExtraPoint12->_pDevData;
    gpu->sim.pExtraPoint12Y                 = gpu->pbExtraPoint12->_pDevData + EP12Stride;
    gpu->sim.pExtraPoint12Z                 = gpu->pbExtraPoint12->_pDevData + EP12Stride * 2;
    gpu->sim.pExtraPoint21Frame             = gpu->pbExtraPoint21Frame->_pDevData;
    gpu->sim.pExtraPoint21Index             = gpu->pbExtraPoint21Index->_pDevData;
    gpu->sim.pExtraPoint21X1                = gpu->pbExtraPoint21->_pDevData;
    gpu->sim.pExtraPoint21Y1                = gpu->pbExtraPoint21->_pDevData + EP21Stride;
    gpu->sim.pExtraPoint21Z1                = gpu->pbExtraPoint21->_pDevData + EP21Stride * 2;    
    gpu->sim.pExtraPoint21X2                = gpu->pbExtraPoint21->_pDevData + EP21Stride * 3;
    gpu->sim.pExtraPoint21Y2                = gpu->pbExtraPoint21->_pDevData + EP21Stride * 4;
    gpu->sim.pExtraPoint21Z2                = gpu->pbExtraPoint21->_pDevData + EP21Stride * 5;        
    gpu->sim.pExtraPoint22Frame             = gpu->pbExtraPoint22Frame->_pDevData;
    gpu->sim.pExtraPoint22Index             = gpu->pbExtraPoint22Index->_pDevData;
    gpu->sim.pExtraPoint22X1                = gpu->pbExtraPoint22->_pDevData;
    gpu->sim.pExtraPoint22Y1                = gpu->pbExtraPoint22->_pDevData + EP22Stride;
    gpu->sim.pExtraPoint22Z1                = gpu->pbExtraPoint22->_pDevData + EP22Stride * 2;    
    gpu->sim.pExtraPoint22X2                = gpu->pbExtraPoint22->_pDevData + EP22Stride * 3;
    gpu->sim.pExtraPoint22Y2                = gpu->pbExtraPoint22->_pDevData + EP22Stride * 4;
    gpu->sim.pExtraPoint22Z2                = gpu->pbExtraPoint22->_pDevData + EP22Stride * 5;           
    gpuCopyConstants();
#if 0    
    printf("%d Extra Point Frames, %d EP11, %d EP12, %d EP21 %d EP22\n", *frame_cnt, EP11s, EP12s, EP21s, EP22s);
    printf("Extra Point Offsets: %d EP11, %d EP12, %d EP21 %d EP22\n", gpu->sim.EP11Offset, gpu->sim.EP12Offset, gpu->sim.EP21Offset, gpu->sim.EP22Offset);
    exit(-1);
#endif
}

extern "C" void gpu_constraints_setup_(int* natc, int* atm_jrc, double* atm_weight, double* atm_xc)
{
PRINTMETHOD("gpu_constraints_setup");
    
    // Delete old constraints
    delete gpu->pbConstraint1;
    delete gpu->pbConstraint2;
    delete gpu->pbConstraintID;

    // Punt if any constraint pointer is NULL
    if ((atm_jrc == NULL) || (atm_weight == NULL) || (atm_xc == NULL))
        return;

    int constraints = *natc;
#ifdef GVERBOSE
    printf("%d constraints\n", constraints);
#endif   

    // Allocate/reallocate constraint buffers 
    gpu->pbConstraint1                              = new GpuBuffer<PMEDouble2>(constraints);
    gpu->pbConstraint2                              = new GpuBuffer<PMEDouble2>(constraints);
    gpu->pbConstraintID                             = new GpuBuffer<int2>(constraints);
    
    // Copy constraints
    for (int i = 0; i < constraints; i++)
    {
        int j = atm_jrc[i] - 1;
        gpu->pbConstraint1->_pSysData[i].x          = atm_weight[i];
        gpu->pbConstraint1->_pSysData[i].y          = atm_xc[3 * j];
        gpu->pbConstraint2->_pSysData[i].x          = atm_xc[3 * j + 1];
        gpu->pbConstraint2->_pSysData[i].y          = atm_xc[3 * j + 2];
        gpu->pbConstraintID->_pSysData[i].x         = j;
        gpu->pbConstraintID->_pSysData[i].y         = j + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[j]++;
    }
    gpu->pbConstraint1->Upload();
    gpu->pbConstraint2->Upload();
    gpu->pbConstraintID->Upload();

    // Set constants
    gpu->sim.constraints                            = constraints;
    gpu->sim.pConstraint1                           = gpu->pbConstraint1->_pDevData;
    gpu->sim.pConstraint2                           = gpu->pbConstraint2->_pDevData;
    gpu->sim.pConstraintID                          = gpu->pbConstraintID->_pDevData;
    gpuCopyConstants(); 
}

extern "C" void gpu_angles_ub_setup_(int* angle_ub_cnt, angle_ub_rec angle_ub[], double ub_r0[], double ub_rk[])
{
PRINTMETHOD("gpu_angles_ub_setup");

    delete gpu->pbUBAngle;
    delete gpu->pbUBAngleID;


#ifdef GVERBOSE
    printf("%d Urey Bradley angles\n", *angle_ub_cnt);
#endif
    int UBAngles                                    = *angle_ub_cnt;

    // Allocate/reallocate Urey Bradley angles
    gpu->pbUBAngle                                  = new GpuBuffer<PMEDouble2>(UBAngles);
    gpu->pbUBAngleID                                = new GpuBuffer<int4>(UBAngles);

    // Add Urey Bradley angles
    UBAngles                                        = 0;
    for (int i = 0; i < *angle_ub_cnt; i++)
    {
        gpu->pbUBAngle->_pSysData[UBAngles].x       = ub_rk[angle_ub[i].parm_idx - 1];
        gpu->pbUBAngle->_pSysData[UBAngles].y       = ub_r0[angle_ub[i].parm_idx - 1];
        gpu->pbUBAngleID->_pSysData[UBAngles].x     = abs(angle_ub[i].atm_i) - 1;
        gpu->pbUBAngleID->_pSysData[UBAngles].y     = abs(angle_ub[i].atm_j) - 1;
        gpu->pbUBAngleID->_pSysData[UBAngles].z     = abs(angle_ub[i].atm_i) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(angle_ub[i].atm_i) - 1]++;
        gpu->pbUBAngleID->_pSysData[UBAngles].w     = abs(angle_ub[i].atm_j) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(angle_ub[i].atm_j) - 1]++;
        UBAngles++;
    }
    gpu->pbUBAngle->Upload();
    gpu->pbUBAngleID->Upload();
    
    // Set constants
    gpu->sim.UBAngles                               = UBAngles;
    gpu->sim.pUBAngle                               = gpu->pbUBAngle->_pDevData;
    gpu->sim.pUBAngleID                             = gpu->pbUBAngleID->_pDevData;
    gpuCopyConstants();            
}
extern "C" void  gpu_dihedrals_imp_setup_(int* dihed_imp_cnt, dihed_imp_rec dihed_imp[], double imp_pk[], double imp_phase[])
{
PRINTMETHOD("gpu_dihedrals_imp_setup");
    // Delete old improper dihedrals
    delete gpu->pbImpDihedral;
    delete gpu->pbImpDihedralID1;
    delete gpu->pbImpDihedralID2;
    
#ifdef GVERBOSE
    printf("%d improper dihedrals\n", *dihed_imp_cnt);
#endif
    int impDihedrals                                = *dihed_imp_cnt;

    // Allocate/reallocate improper dihedrals
    gpu->pbImpDihedral                              = new GpuBuffer<PMEDouble2>(impDihedrals);
    gpu->pbImpDihedralID1                           = new GpuBuffer<int4>(impDihedrals);
    gpu->pbImpDihedralID2                           = new GpuBuffer<int4>(impDihedrals);
    
    // Add improper dihedrals
    impDihedrals                                    = 0;  
    for (int i = 0; i < *dihed_imp_cnt; i++)
    {
        gpu->pbImpDihedral->_pSysData[impDihedrals].x       = imp_pk[dihed_imp[i].parm_idx - 1];
        gpu->pbImpDihedral->_pSysData[impDihedrals].y       = imp_phase[dihed_imp[i].parm_idx - 1];
        gpu->pbImpDihedralID1->_pSysData[impDihedrals].x    = abs(dihed_imp[i].atm_i) - 1;
        gpu->pbImpDihedralID1->_pSysData[impDihedrals].y    = abs(dihed_imp[i].atm_j) - 1;
        gpu->pbImpDihedralID1->_pSysData[impDihedrals].z    = abs(dihed_imp[i].atm_k) - 1;
        gpu->pbImpDihedralID1->_pSysData[impDihedrals].w    = abs(dihed_imp[i].atm_l) - 1;
        gpu->pbImpDihedralID2->_pSysData[impDihedrals].x    = abs(dihed_imp[i].atm_i) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(dihed_imp[i].atm_i) - 1]++;
        gpu->pbImpDihedralID2->_pSysData[impDihedrals].y    = abs(dihed_imp[i].atm_j) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(dihed_imp[i].atm_j) - 1]++;
        gpu->pbImpDihedralID2->_pSysData[impDihedrals].z    = abs(dihed_imp[i].atm_k) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(dihed_imp[i].atm_k) - 1]++;
        gpu->pbImpDihedralID2->_pSysData[impDihedrals].w    = abs(dihed_imp[i].atm_l) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(dihed_imp[i].atm_l) - 1]++;
        impDihedrals++;
    }
    gpu->pbImpDihedral->Upload();
    gpu->pbImpDihedralID1->Upload();
    gpu->pbImpDihedralID2->Upload();
    
    // Set constants
    gpu->sim.impDihedrals                           = impDihedrals;
    gpu->sim.pImpDihedral                           = gpu->pbImpDihedral->_pDevData;
    gpu->sim.pImpDihedralID1                        = gpu->pbImpDihedralID1->_pDevData;
    gpu->sim.pImpDihedralID2                        = gpu->pbImpDihedralID2->_pDevData;
    gpuCopyConstants();            
}

extern "C" void gpu_cmap_setup_(int* cmap_cnt, cmap_rec cmap[], int* cmap_type_count, double cmap_grid[], double cmap_dPhi[], double cmap_dPsi[], double cmap_dPhi_dPsi[]) 
{
PRINTMETHOD("gpu_cmap_setup");

    // Delete old cmap data
    delete gpu->pbCmapID1;
    delete gpu->pbCmapID2;
    delete gpu->pbCmapID3;
    delete gpu->pbCmapType;
    delete gpu->pbCmapEnergy;

#ifdef GVERBOSE
    printf("%d cmaps\n", *cmap_cnt);
#endif
    int cmaps                                       = *cmap_cnt;
    
    // Allocate/reallocate cmaps and cmap data
    gpu->pbCmapID1                                  = new GpuBuffer<int4>(cmaps);
    gpu->pbCmapID2                                  = new GpuBuffer<int4>(cmaps);
    gpu->pbCmapID3                                  = new GpuBuffer<int2>(cmaps);
    gpu->pbCmapType                                 = new GpuBuffer<int>(cmaps);
    
    // Add cmaps
    cmaps                                           = 0;
    for (int i = 0; i < *cmap_cnt; i++)
    {
        gpu->pbCmapID1->_pSysData[cmaps].x          = abs(cmap[i].atm_i) - 1;
        gpu->pbCmapID1->_pSysData[cmaps].y          = abs(cmap[i].atm_j) - 1;
        gpu->pbCmapID1->_pSysData[cmaps].z          = abs(cmap[i].atm_k) - 1;
        gpu->pbCmapID1->_pSysData[cmaps].w          = abs(cmap[i].atm_l) - 1;
        gpu->pbCmapID2->_pSysData[cmaps].x          = abs(cmap[i].atm_m) - 1;
        gpu->pbCmapID2->_pSysData[cmaps].y          = abs(cmap[i].atm_i) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cmap[i].atm_i) - 1]++;
        gpu->pbCmapID2->_pSysData[cmaps].z          = abs(cmap[i].atm_j) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cmap[i].atm_j) - 1]++;
        gpu->pbCmapID2->_pSysData[cmaps].w          = abs(cmap[i].atm_k) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cmap[i].atm_k) - 1]++;
        gpu->pbCmapID3->_pSysData[cmaps].x          = abs(cmap[i].atm_l) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cmap[i].atm_l) - 1]++;
        gpu->pbCmapID3->_pSysData[cmaps].y          = abs(cmap[i].atm_m) - 1 + gpu->sim.stride3 * gpu->pbOutputBufferCounter->_pSysData[abs(cmap[i].atm_m) - 1]++;        
        gpu->pbCmapType->_pSysData[cmaps]           = cmap[i].parm_idx - 1;
        cmaps++;
    }
    gpu->pbCmapID1->Upload();
    gpu->pbCmapID2->Upload();
    gpu->pbCmapID3->Upload();
    gpu->pbCmapType->Upload();
    
    
    // Allocate energy map
    int termReads                                   = (*cmap_type_count * sizeof(PMEDouble) + gpu->readSize - 1) / gpu->readSize;
    int cmapTermStride                              = (termReads * gpu->readSize) / sizeof(PMEDouble);   
    int cmapRowStride                               = cmapTermStride * CMAPDIMENSION;
    gpu->pbCmapEnergy                               = new GpuBuffer<PMEFloat4>(cmapTermStride * CMAPDIMENSION * CMAPDIMENSION); 
    
    // Copy cmap terms
    for (int i = 0; i < *cmap_type_count; i++)
    {
        PMEFloat4* pCmap                            = gpu->pbCmapEnergy->_pSysData + i + cmapRowStride + cmapTermStride;
        for (int x = -1; x < (int)CMAPRESOLUTION + 3; x++)
        {
            int ix = x;
            if (ix < 0)
                ix                                 += CMAPRESOLUTION;
            else if (ix >= (int)CMAPRESOLUTION)
                ix                                 -= CMAPRESOLUTION;
            for (int y = -1; y < (int)CMAPRESOLUTION + 3; y++)
            {
                int iy = y;
                if (iy < 0)
                    iy                             += CMAPRESOLUTION;
                else if (iy >= (int)CMAPRESOLUTION)
                    iy                             -= CMAPRESOLUTION;
                pCmap[y * cmapRowStride + 
                      x * cmapTermStride].x         = cmap_grid[i + iy * *cmap_type_count + ix * CMAPRESOLUTION * *cmap_type_count]; 
                pCmap[y * cmapRowStride + 
                      x * cmapTermStride].y         = cmap_dPhi[i + iy * *cmap_type_count + ix * CMAPRESOLUTION * *cmap_type_count] * CMAPSTEPSIZE; 
                pCmap[y * cmapRowStride + 
                      x * cmapTermStride].z         = cmap_dPsi[i + iy * *cmap_type_count + ix * CMAPRESOLUTION * *cmap_type_count] * CMAPSTEPSIZE; 
                pCmap[y * cmapRowStride + 
                      x * cmapTermStride].w         = cmap_dPhi_dPsi[i + iy * *cmap_type_count + ix * CMAPRESOLUTION * *cmap_type_count] * CMAPSTEPSIZE * CMAPSTEPSIZE;       
            }
        }
    }
    gpu->pbCmapEnergy->Upload();
    
    // Set constants
    gpu->sim.cmaps                                  = cmaps;
    gpu->sim.pCmapID1                               = gpu->pbCmapID1->_pDevData;
    gpu->sim.pCmapID2                               = gpu->pbCmapID2->_pDevData;
    gpu->sim.pCmapID3                               = gpu->pbCmapID3->_pDevData;
    gpu->sim.pCmapType                              = gpu->pbCmapType->_pDevData;
    gpu->sim.cmapTermStride                         = cmapTermStride;
    gpu->sim.cmapRowStride                          = cmapRowStride;
    gpu->sim.pCmapEnergy                            = gpu->pbCmapEnergy->_pDevData + cmapTermStride + cmapRowStride;              
    gpuCopyConstants();     
}

extern "C" void gpu_shake_setup_(double atm_mass[], int* my_nonfastwat_bond_cnt, shake_bond_rec my_nonfastwat_bond_dat[], int* my_fastwat_res_cnt, int* iorwat, int my_fastwat_res_lst[])
{
PRINTMETHOD("gpu_shake_setup");
    // Delete any existing Shake constraints
    delete gpu->pbShakeID;
    delete gpu->pbShakeParm;
    delete gpu->pbFastShakeID;
    delete gpu->pbSlowShakeID1;
    delete gpu->pbSlowShakeID2;
    delete gpu->pbSlowShakeParm;

#ifdef GVERBOSE    
    printf("%d nonfast SHAKE constraints, %d fast SHAKE constraints\n", *my_nonfastwat_bond_cnt, *my_fastwat_res_cnt);
#endif
    
    // Determine number of H-bond networks in traditional Shake constraint list
    int shakeConstraints                                = 0;
    int slowShakeConstraints                            = 0;
    bool* shakeMap                                      = new bool[*my_nonfastwat_bond_cnt];
    int4* pShakeID                                      = new int4[*my_nonfastwat_bond_cnt];
    int* pSlowShakeID1                                  = new int[*my_nonfastwat_bond_cnt];
    int4* pSlowShakeID2                                 = new int4[*my_nonfastwat_bond_cnt];
    double2* pShakeParm                                 = new double2[*my_nonfastwat_bond_cnt];
    double2* pSlowShakeParm                             = new double2[*my_nonfastwat_bond_cnt];

    for (int i = 0; i < *my_nonfastwat_bond_cnt; i++)
    {
        shakeMap[i]                                     = false;
    }
    
    bool fail = false;
    for (int i = 0; i < *my_nonfastwat_bond_cnt; i++)
    {
        if (!shakeMap[i])
        {
            shakeMap[i]                                 = true;
            unsigned int atm_i;
            unsigned int atm_j;
            double mass_i;
            double mass_j;
            double r                                    = my_nonfastwat_bond_dat[i].parm;
        
            // Determine central atom in network
            if (atm_mass[my_nonfastwat_bond_dat[i].atm_i - 1] > atm_mass[my_nonfastwat_bond_dat[i].atm_j - 1])
            {
                atm_i                                   = my_nonfastwat_bond_dat[i].atm_i;
                atm_j                                   = my_nonfastwat_bond_dat[i].atm_j;
                mass_i                                  = atm_mass[my_nonfastwat_bond_dat[i].atm_i - 1];
                mass_j                                  = atm_mass[my_nonfastwat_bond_dat[i].atm_j - 1];
            }
            else           
            {
                atm_i                                   = my_nonfastwat_bond_dat[i].atm_j;
                atm_j                                   = my_nonfastwat_bond_dat[i].atm_i;
                mass_i                                  = atm_mass[my_nonfastwat_bond_dat[i].atm_j - 1];
                mass_j                                  = atm_mass[my_nonfastwat_bond_dat[i].atm_i - 1];
            }
            unsigned int atm_i_cnt                      = 0;
            int atom_jj[4];
            double mass_jj[4];
            double r_jj[4];
            
            for (int j = 0; j < 3; j++)
            {
                atom_jj[j]                              = 0;
                mass_jj[j]                              = 0.0;
                r_jj[j]                                 = -10.0;
            }
            
            for (int j = i + 1; j <  *my_nonfastwat_bond_cnt; j++)
            {
                if (!shakeMap[j])
                {
                    // Check for additional members of network
                    if (my_nonfastwat_bond_dat[j].atm_i == atm_i)
                    {
                        atom_jj[atm_i_cnt]              = my_nonfastwat_bond_dat[j].atm_j;
                        mass_jj[atm_i_cnt]              = atm_mass[my_nonfastwat_bond_dat[j].atm_j - 1];
                        r_jj[atm_i_cnt]                 = my_nonfastwat_bond_dat[j].parm;
                        shakeMap[j]                     = true;
                    }
                    else if (my_nonfastwat_bond_dat[j].atm_j == atm_i)
                    {
                        atom_jj[atm_i_cnt]              = my_nonfastwat_bond_dat[j].atm_i;
                        mass_jj[atm_i_cnt]              = atm_mass[my_nonfastwat_bond_dat[j].atm_i - 1];
                        r_jj[atm_i_cnt]                 = my_nonfastwat_bond_dat[j].parm;
                        shakeMap[j]                     = true;
                    }
                    
                    // Check for inconsistent network
                    if ((my_nonfastwat_bond_dat[j].atm_i == atm_j) || (my_nonfastwat_bond_dat[j].atm_j == atm_j))
                    {
                        printf("Hydrogen atom %d appears to have multiple bonds to atoms %d and %d which is illegal for SHAKEH.\n", atm_j, atm_i, my_nonfastwat_bond_dat[j].atm_j);
                        fail                            = true;
                    }
                    
                    if (shakeMap[j])
                    {
                        atm_i_cnt++;
                        if (atm_i_cnt == 4)
                        {
                            printf("Too many hydrogens for a hydrogen network, exiting.\n");
                            gpu_shutdown_();
                            exit(-1);
                        }
                    }
                }
            }

#if 0       
            printf("%5d: %16.7f %16.7f ", shakeConstraints, mass_i, mass_j);
            for (int ik = 0; ik < atm_i_cnt; ik++)
                printf("%16.7f ", mass_jj[ik]);
            printf("\n");
#endif       

#if 0       
            printf("%5d: %16.7f ", shakeConstraints, r);
            for (int ik = 0; ik < atm_i_cnt; ik++)
                printf("%16.7f ", r_jj[ik]);
            printf("\n");
#endif    
       
            // Construct constraint
            if (atm_i_cnt < 3)
            {
                pShakeID[shakeConstraints].x            = atm_i - 1;
                pShakeID[shakeConstraints].y            = atm_j - 1;
                pShakeID[shakeConstraints].z            = atom_jj[0] - 1;
                pShakeID[shakeConstraints].w            = atom_jj[1] - 1;
                pShakeParm[shakeConstraints].x          = 1.0 / mass_i;
                pShakeParm[shakeConstraints].y          = r;
            
                // Check for consistency
                //printf("%5d: %9d %9d %9d %9d\n", shakeConstraints, atm_i - 1, atm_j - 1, atom_jj[0] - 1, atom_jj[1] - 1);
                //printf("%5d: %9.4f %9.4f %9.4f %9.4f\n", shakeConstraints, mass_i, mass_j, mass_jj[0], mass_jj[1]); 
                //printf("%5d: %9.4f %9.4f %9.4f\n", shakeConstraints, r, r_jj[0], r_jj[1]); 
                shakeConstraints++;
            }
            else
            {
                pSlowShakeID1[slowShakeConstraints]     = atm_i - 1;
                pSlowShakeID2[slowShakeConstraints].x   = atm_j - 1;
                pSlowShakeID2[slowShakeConstraints].y   = atom_jj[0] - 1;
                pSlowShakeID2[slowShakeConstraints].z   = atom_jj[1] - 1;
                pSlowShakeID2[slowShakeConstraints].w   = atom_jj[2] - 1;
                pSlowShakeParm[slowShakeConstraints].x  = 1.0 / mass_i;
                pSlowShakeParm[slowShakeConstraints].y  = r;
                slowShakeConstraints++;
            }
        }
        
        // Exit if one or more illegal SHAKEH networks detected
        if (fail)
        {
            printf("Exiting due to the presence of inconsistent SHAKEH hydrogen clusters.\n");
            gpu_shutdown_();
            exit(-1);
        }
    }

    // Allocate and copy standard SHAKE constraints to gpu RAM
    gpu->sim.shakeConstraints                       = shakeConstraints;
    gpu->pbShakeID                                  = new GpuBuffer<int4>(shakeConstraints);
    gpu->pbShakeParm                                = new GpuBuffer<double2>(shakeConstraints);
    gpu->sim.pShakeID                               = gpu->pbShakeID->_pDevData;
    gpu->sim.pShakeParm                             = gpu->pbShakeParm->_pDevData;
    for (int i = 0; i < shakeConstraints; i++)
    {
        gpu->pbShakeID->_pSysData[i]                = pShakeID[i];
        gpu->pbShakeParm->_pSysData[i]              = pShakeParm[i];
    }
    gpu->pbShakeID->Upload();
    gpu->pbShakeParm->Upload();

    // Allocate and copy rarely occurring sp3-hybridized SHAKE constraints to GPU RAM
    gpu->sim.slowShakeConstraints                   = slowShakeConstraints;
    gpu->pbSlowShakeID1                             = new GpuBuffer<int>(slowShakeConstraints);
    gpu->pbSlowShakeID2                             = new GpuBuffer<int4>(slowShakeConstraints);
    gpu->pbSlowShakeParm                            = new GpuBuffer<double2>(slowShakeConstraints);
    gpu->sim.pSlowShakeID1                          = gpu->pbSlowShakeID1->_pDevData;
    gpu->sim.pSlowShakeID2                          = gpu->pbSlowShakeID2->_pDevData;
    gpu->sim.pSlowShakeParm                         = gpu->pbSlowShakeParm->_pDevData;
    for (int i = 0; i < slowShakeConstraints; i++)
    {
        gpu->pbSlowShakeID1->_pSysData[i]           = pSlowShakeID1[i];
        gpu->pbSlowShakeID2->_pSysData[i]           = pSlowShakeID2[i];
        gpu->pbSlowShakeParm->_pSysData[i]          = pSlowShakeParm[i];
    }
    gpu->pbSlowShakeID1->Upload();
    gpu->pbSlowShakeID2->Upload();
    gpu->pbSlowShakeParm->Upload();
    
    int ind1, ind2, ind3;
    switch (*iorwat)
    {
        case 1:
            ind1                                    = 0;
            ind2                                    = 1;
            ind3                                    = 2;
            break;
            
        case 2:
            ind1                                    = 1;
            ind2                                    = 2;
            ind3                                    = 0;
            break;
            
        default:
            ind1                                    = 2;
            ind2                                    = 0;
            ind3                                    = 1;
            break;
    }
    
    // Set up fast SHAKE constraints
    gpu->pbFastShakeID                              = new GpuBuffer<int4>(*my_fastwat_res_cnt);
    for (int i = 0; i < *my_fastwat_res_cnt; i++)
    {
        int atm_1                                   = my_fastwat_res_lst[i] + ind1 - 1;
        int atm_2                                   = my_fastwat_res_lst[i] + ind2 - 1;
        int atm_3                                   = my_fastwat_res_lst[i] + ind3 - 1;
        gpu->pbFastShakeID->_pSysData[i].x          = atm_1;
        gpu->pbFastShakeID->_pSysData[i].y          = atm_2;
        gpu->pbFastShakeID->_pSysData[i].z          = atm_3;
        gpu->pbFastShakeID->_pSysData[i].w          = -1;
    }
    gpu->sim.fastShakeConstraints                   = *my_fastwat_res_cnt;
    gpu->pbFastShakeID->Upload();
    gpu->sim.pFastShakeID                           = gpu->pbFastShakeID->_pDevData;

    // Calculate SHAKE offsets
    gpu->sim.shakeOffset                            = (((gpu->sim.shakeConstraints          + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
    gpu->sim.fastShakeOffset                        = (((gpu->sim.fastShakeConstraints      + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
    gpu->sim.slowShakeOffset                        = (((gpu->sim.slowShakeConstraints      + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
    gpuCopyConstants();
    
    // Release temporary space
    delete[] shakeMap;
    delete[] pShakeID;
    delete[] pShakeParm;
    delete[] pSlowShakeID1;
    delete[] pSlowShakeID2;
    delete[] pSlowShakeParm;
}

// Set up fast SHAKE parameters
extern "C" void gpu_get_water_distances_(double* rbtarg)
{
    gpu->sim.ra                                     = rbtarg[0];
    gpu->sim.ra_inv                                 = 1.0 / gpu->sim.ra;
    gpu->sim.rb                                     = rbtarg[1];
    gpu->sim.rc                                     = rbtarg[2];
    gpu->sim.rc2                                    = rbtarg[3];
    gpu->sim.hhhh                                   = rbtarg[4];
    gpu->sim.wo_div_wohh                            = rbtarg[5] / rbtarg[7];
    gpu->sim.wh_div_wohh                            = rbtarg[6] / rbtarg[7];
    gpuCopyConstants();
}

extern "C" void gpu_create_outputbuffers_()
{
PRINTMETHOD("gpu_create_outputbuffers");

    // Determine size needed based on presence or absence of neighbor list
    int nonbondForceBuffers;
    int maxForceBuffers;
    if (gpu->bNeighborList) 
    {
        nonbondForceBuffers                         = max((int)gpu->sim.NLCellBuffers, gpu->sim.maxNonbonds);
        maxForceBuffers                             = nonbondForceBuffers;
    }
    else
    {
        nonbondForceBuffers                         = gpu->sim.paddedNumberOfAtoms / GRID;
        maxForceBuffers                             = 2 * nonbondForceBuffers;
    }
    
    if (gpu->sim.ntp == 0)
    {
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            if (gpu->pbOutputBufferCounter->_pSysData[i] > maxForceBuffers)
            {
                maxForceBuffers                         = gpu->pbOutputBufferCounter->_pSysData[i];
            }   
        }
    }
    gpu->sim.maxForceBuffers                        = maxForceBuffers;
    gpu->sim.nonbondForceBuffers                    = nonbondForceBuffers;
    gpu->sim.maxNonbondBuffers                      = nonbondForceBuffers * gpu->sim.stride3;
 
    // delete any existing output buffers
    delete gpu->pbForceBuffer;
    delete gpu->pbReffBuffer;
    delete gpu->pbSumdeijdaBuffer;
    gpu->pbReffBuffer                               = NULL;
    gpu->pbSumdeijdaBuffer                          = NULL;
    delete gpu->pbXYZ_q;
    delete gpu->pbXYZ_qt;
    gpu->pbXYZ_q                                    = NULL;
    gpu->pbXYZ_qt                                   = NULL;
    delete gpu->pbEnergyBuffer;
    delete gpu->pbKineticEnergyBuffer;
    
    // Only allocate GB buffers if they're needed
    if (gpu->ntb == 0)
    {
        // Allocate Born Radius buffer
        gpu->pbReffBuffer                           = new GpuBuffer<PMEDouble>(nonbondForceBuffers * gpu->sim.stride, bShadowedOutputBuffers);
        gpu->sim.pReffBuffer                        = gpu->pbReffBuffer->_pDevData;       
    
        // Allocate Born Force buffer
        gpu->pbSumdeijdaBuffer                      = new GpuBuffer<PMEDouble>(nonbondForceBuffers * gpu->sim.stride, bShadowedOutputBuffers);
        gpu->sim.pSumdeijdaBuffer                   = gpu->pbSumdeijdaBuffer->_pDevData;    
    }
    else
    {
        // Allocate PME buffers if needed     
        gpu->pbXYZ_q                                = new GpuBuffer<PMEFloat>(gpu->sim.maxChargeGridBuffers * gpu->sim.XYZStride, bShadowedOutputBuffers);
        gpu->pbXYZ_qt                               = new GpuBuffer<PMEComplex>(gpu->sim.fft_x_dim * gpu->sim.fft_y_dim * gpu->sim.fft_z_dim);
        gpu->sim.pXYZ_q                             = gpu->pbXYZ_q->_pDevData;
        gpu->sim.plXYZ_q                            = (int*)(gpu->pbXYZ_q->_pDevData);
        gpu->sim.pXYZ_qt                            = gpu->pbXYZ_qt->_pDevData;
    }
    
    
    
    // Allocate new force buffers
    gpu->pbForceBuffer                              = new GpuBuffer<PMEDouble>(gpu->sim.stride3 * maxForceBuffers, bShadowedOutputBuffers);
    gpu->sim.pForceBuffer                           = gpu->pbForceBuffer->_pDevData;
    gpu->sim.pForceXBuffer                          = gpu->sim.pForceBuffer;
    gpu->sim.pForceYBuffer                          = gpu->sim.pForceBuffer + gpu->sim.stride;
    gpu->sim.pForceZBuffer                          = gpu->sim.pForceBuffer + gpu->sim.stride2;
    

    
    // Allocate energy buffer
    gpu->pbEnergyBuffer                             = new GpuBuffer<unsigned long long int>(ENERGYTERMS);
    gpu->sim.pEnergyBuffer                          = gpu->pbEnergyBuffer->_pDevData;
    gpu->sim.pEELT                                  = gpu->sim.pEnergyBuffer + 0;
    gpu->sim.pEVDW                                  = gpu->sim.pEnergyBuffer + 1;
    gpu->sim.pEGB                                   = gpu->sim.pEnergyBuffer + 2;
    gpu->sim.pEBond                                 = gpu->sim.pEnergyBuffer + 3;
    gpu->sim.pEAngle                                = gpu->sim.pEnergyBuffer + 4;
    gpu->sim.pEDihedral                             = gpu->sim.pEnergyBuffer + 5;
    gpu->sim.pEEL14                                 = gpu->sim.pEnergyBuffer + 6;
    gpu->sim.pENB14                                 = gpu->sim.pEnergyBuffer + 7;
    gpu->sim.pERestraint                            = gpu->sim.pEnergyBuffer + 8;
    gpu->sim.pEER                                   = gpu->sim.pEnergyBuffer + 9;
    gpu->sim.pEED                                   = gpu->sim.pEnergyBuffer + 10;
    gpu->sim.pEAngle_UB                             = gpu->sim.pEnergyBuffer + 11;
    gpu->sim.pEImp                                  = gpu->sim.pEnergyBuffer + 12;
    gpu->sim.pECmap                                 = gpu->sim.pEnergyBuffer + 13;
    gpu->sim.pVirial                                = gpu->sim.pEnergyBuffer + VIRIALOFFSET;
    gpu->sim.pVirial_11                             = gpu->sim.pEnergyBuffer + VIRIALOFFSET;
    gpu->sim.pVirial_22                             = gpu->sim.pEnergyBuffer + VIRIALOFFSET + 1;
    gpu->sim.pVirial_33                             = gpu->sim.pEnergyBuffer + VIRIALOFFSET + 2;
    gpu->sim.pEKCOMX                                = gpu->sim.pEnergyBuffer + VIRIALOFFSET + 3;
    gpu->sim.pEKCOMY                                = gpu->sim.pEnergyBuffer + VIRIALOFFSET + 4;
    gpu->sim.pEKCOMZ                                = gpu->sim.pEnergyBuffer + VIRIALOFFSET + 5;
    
    // Allocate kinetic energy buffer
    gpu->pbKineticEnergyBuffer                      = new GpuBuffer<KineticEnergy>(gpu->blocks, false, true);
    gpu->sim.pKineticEnergy                         = gpu->pbKineticEnergyBuffer->_pDevData;
    
    // Calculate offsets
    gpu->sim.bondOffset                             =                                 (((gpu->sim.bonds          + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
    gpu->sim.bondAngleOffset                        = gpu->sim.bondOffset           + (((gpu->sim.bondAngles     + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
    gpu->sim.dihedralOffset                         = gpu->sim.bondAngleOffset      + (((gpu->sim.dihedrals      + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
    gpu->sim.nb14Offset                             = gpu->sim.dihedralOffset       + (((gpu->sim.nb14s          + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
    gpu->sim.constraintOffset                       = gpu->sim.nb14Offset           + (((gpu->sim.constraints    + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
    gpuCopyConstants();
    
    if (!gpu->bNeighborList)
        kClearGBBuffers(gpu);
    
    // Check for zero local interations and skip local force kernel in that case.
    if (gpu->sim.bonds + gpu->sim.bondAngles + gpu->sim.dihedrals + gpu->sim.nb14s + gpu->sim.constraints > 0)
        gpu->bLocalInteractions                     = true;
    else
        gpu->bLocalInteractions                     = false;       
        
        
    // Calculate charmm offsets
    gpu->sim.UBAngleOffset                          =                                 (((gpu->sim.UBAngles       + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
    gpu->sim.impDihedralOffset                      = gpu->sim.UBAngleOffset        + (((gpu->sim.impDihedrals   + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
    gpu->sim.cmapOffset                             = gpu->sim.impDihedralOffset    + (((gpu->sim.cmaps          + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);

    if (gpu->sim.UBAngles + gpu->sim.impDihedrals + gpu->sim.cmaps > 0)
        gpu->bCharmmInteractions                    = true;
    else
        gpu->bCharmmInteractions                    = false;   

#if 0   
    printf("%d %d\n", gpu->sim.bonds, gpu->sim.bondOffset);
    printf("%d %d\n", gpu->sim.bondAngles, gpu->sim.bondAngleOffset);
    printf("%d %d\n", gpu->sim.dihedrals, gpu->sim.dihedralOffset);
    printf("%d %d\n", gpu->sim.nb14s, gpu->sim.nb14Offset);
    printf("%d %d\n", gpu->sim.constraints, gpu->sim.constraintOffset);
    printf("%d %d\n", gpu->sim.UBAngles, gpu->sim.UBAngleOffset);
    printf("%d %d\n", gpu->sim.impDihedrals, gpu->sim.impDihedralOffset);
    printf("%d %d\n", gpu->sim.cmaps, gpu->sim.cmapOffset);
    exit(-1);
#endif
}


extern "C" void gpuCopyConstants()
{
PRINTMETHOD("gpuCopyConstants");
    SetkForcesUpdateSim(gpu);
    SetkShakeSim(gpu);
    SetkCalculateLocalForcesSim(gpu);
    SetkCalculateAMDWeightsSim(gpu);
    if (gpu->ntb == 0)
    {
        SetkCalculateGBBornRadiiSim(gpu);
        SetkCalculateGBNonbondEnergy1Sim(gpu);
        SetkCalculateGBNonbondEnergy2Sim(gpu);
    }
    else
    {
        SetkPMEInterpolationSim(gpu);
        SetkNeighborListSim(gpu);
        SetkCalculatePMENonbondEnergySim(gpu);
    }
}

extern "C" void gpu_init_pbc_(double* a, double* b, double* c, double* alpha, double* beta, double* gamma, double* uc_volume, double* uc_sphere, double* max_cutoff, double pbc_box[], double reclng[], double cut_factor[], double ucell[][3], double recip[][3])
{
PRINTMETHOD("gpu_init_pbc");

    gpu->bNeighborList                      = true;
    gpu->bNeedNewNeighborList               = true;
    gpu->bNewNeighborList                   = false;

    // Check for orthogonal system
    if ((*alpha == 90.0) && (*beta  == 90.0) && (*gamma == 90.0))    
        gpu->sim.is_orthog                  = true;
    else
        gpu->sim.is_orthog                  = false;

    gpu->sim.af                             = (PMEFloat)*a;
    gpu->sim.bf                             = (PMEFloat)*b;
    gpu->sim.cf                             = (PMEFloat)*c;
    gpu->sim.a                              = (PMEDouble)gpu->sim.af;
    gpu->sim.b                              = (PMEDouble)gpu->sim.bf;
    gpu->sim.c                              = (PMEDouble)gpu->sim.cf;
    gpu->sim.alpha                          = *alpha * PI / 180.0;
    gpu->sim.beta                           = *beta * PI / 180.0;
    gpu->sim.gamma                          = *gamma * PI / 180.0;

	for (int i = 0; i < 3; i++)
    {
        gpu->sim.pbc_box[i]                 = pbc_box[i];
        gpu->sim.reclng[i]                  = reclng[i];
        gpu->sim.cut_factor[i]              = cut_factor[i]; 
        for (int j = 0; j < 3; j++)
        {
            gpu->sim.ucell[i][j]            = ucell[j][i];
            gpu->sim.recipf[i][j]           = recip[j][i];
			gpu->sim.recip[i][j]            = gpu->sim.recipf[i][j];
        }
    }
                


#if 0
    printf("%14.6f %14.6f %14.6f\n", gpu->sim.ucell[0][0], gpu->sim.ucell[0][1], gpu->sim.ucell[0][2]);
    printf("%14.6f %14.6f %14.6f\n", gpu->sim.ucell[1][0], gpu->sim.ucell[1][1], gpu->sim.ucell[1][2]);
    printf("%14.6f %14.6f %14.6f\n\n", gpu->sim.ucell[2][0], gpu->sim.ucell[2][1], gpu->sim.ucell[2][2]);
    
    printf("%14.6f %14.6f %14.6f\n", gpu->sim.recip[0][0], gpu->sim.recip[0][1], gpu->sim.recip[0][2]);
    printf("%14.6f %14.6f %14.6f\n", gpu->sim.recip[1][0], gpu->sim.recip[1][1], gpu->sim.recip[1][2]);
    printf("%14.6f %14.6f %14.6f\n\n", gpu->sim.recip[2][0], gpu->sim.recip[2][1], gpu->sim.recip[2][2]);    
    
    printf("%14.6f %14.6f %14.6f\n", *alpha, *beta, *gamma);
#endif
    
    gpu->sim.uc_volume                      = *uc_volume;
    gpu->sim.uc_sphere                      = *uc_sphere;
    gpu->sim.pi_vol_inv                     = 1.0 / (PI * gpu->sim.uc_volume);

    gpuCopyConstants();
    return;
}

static const int cellOffset[][3] =    {  
                                        { 0,  0,  0},
                                        { 1,  0,  0},
                                        {-1,  1,  0},
                                        { 0,  1,  0},
                                        { 1,  1,  0},
                                        {-1, -1,  1},
                                        { 0, -1,  1},
                                        { 1, -1,  1},
                                        {-1,  0,  1},
                                        { 0,  0,  1},
                                        { 1,  0,  1},
                                        {-1,  1,  1},
                                        { 0,  1,  1},
                                        { 1,  1,  1}
};
                                 
                                 
static const unsigned int cellHash[CELLHASHCELLS] = {
                        34, 35, 60, 61,
                        33, 32, 63, 62,
                        30, 31,  0,  1,
                        29, 28,  3,  2,
                        
                        37, 36, 59, 58,
                        38, 39, 56, 57,
                        25, 24,  7,  6,
                        26, 27,  4,  5,
                        
                        42, 43, 54, 53,
                        41, 40, 55, 52,
                        20, 23,  8,  9,
                        21, 22, 11, 10,
                        
                        45, 44, 49, 50,
                        46, 47, 48, 51,
                        19, 16, 15, 14,
                        18, 17, 12, 13,
};                         
                                 
extern "C" void gpu_neighbor_list_setup_(int numex[], int natex[], double* vdw_cutoff, double* skinnb)
{
PRINTMETHOD("gpu_neighbor_list_setup");

    // Delete any existing neighbor list data
    delete gpu->pbImageIndex;
    delete gpu->pbAtomXYSaveSP;
    delete gpu->pbAtomZSaveSP;
    delete gpu->pbImage;
    delete gpu->pbImageVel;
    delete gpu->pbImageLVel;
    delete gpu->pbImageMass;
    delete gpu->pbImageCharge;
    delete gpu->pbImageSigEps;
    delete gpu->pbImageOutputBuffers;
    delete gpu->pbImageCellID;
    delete gpu->pbNLCellHash;
    delete gpu->pbNLNonbondCellStartEnd;
    delete gpu->pbNLChargeGridBufferOffset;
    delete gpu->pbNLOddBufferOverlapFlag;
    delete gpu->pbNLExclusionList;
    delete gpu->pbNLExclusionStartCount;
    delete gpu->pbNLAtomList;
    delete gpu->pbNLOffset;
    delete gpu->pbNLTotalOffset;
    delete gpu->pbNLPosition;
    delete gpu->pbNLbSkinTestFail;
    delete gpu->pbImageBondID;                                      
    delete gpu->pbImageBondAngleID1;
    delete gpu->pbImageBondAngleID2;
    delete gpu->pbImageDihedralID1;
    delete gpu->pbImageDihedralID2;
    delete gpu->pbImageNb14ID;         
    delete gpu->pbImageUBAngleID;
    delete gpu->pbImageImpDihedralID1;
    delete gpu->pbImageImpDihedralID2;
    delete gpu->pbImageCmapID1;
    delete gpu->pbImageCmapID2;
    delete gpu->pbImageCmapID3;
    delete gpu->pbImageConstraintID;
    delete gpu->pbImageShakeID;
    delete gpu->pbImageFastShakeID;
    delete gpu->pbImageSlowShakeID1;
    delete gpu->pbImageSlowShakeID2;
    delete gpu->pbImageSoluteAtomID;
    delete gpu->pbImageSolventAtomID;
    delete gpu->pbImageExtraPoint11Frame;
    delete gpu->pbImageExtraPoint11Index;
    delete gpu->pbImageExtraPoint12Frame;
    delete gpu->pbImageExtraPoint12Index;    
    delete gpu->pbImageExtraPoint21Frame;
    delete gpu->pbImageExtraPoint21Index;
    delete gpu->pbImageExtraPoint22Frame;
    delete gpu->pbImageExtraPoint22Index;    
#ifdef MPI
    delete[] gpu->pMinLocalCell;
    delete[] gpu->pMaxLocalCell;
    delete[] gpu->pMinLocalAtom;
    delete[] gpu->pMaxLocalAtom;
    delete[] gpu->pMinProcessedCell;
    delete[] gpu->pMaxProcessedCell;
    delete[] gpu->pMinProcessedAtom;
    delete[] gpu->pMaxProcessedAtom;
    delete[] gpu->pPMEStart;
    delete[] gpu->pPMELength;
    MPI_Free_mem(gpu->pForceData);
    delete[] gpu->pAllGathervRecvCountSoA;
    delete[] gpu->pAllGathervRecvDisplSoA;
    delete[] gpu->pForceSendNode;
    delete[] gpu->pForceSendMinCell;
    delete[] gpu->pForceSendMaxCell;
    delete[] gpu->pOutForceSendStart;
    delete[] gpu->pForceSendStart;
    delete[] gpu->pForceSendLength;     
#endif    
    
    // Copy constants
    gpu->sim.cut                                            = *vdw_cutoff;
    gpu->sim.cut2                                           = (*vdw_cutoff) * (*vdw_cutoff);
    gpu->sim.skinnb                                         = *skinnb;
    gpu->sim.cutPlusSkin                                    = *vdw_cutoff + *skinnb;
    gpu->sim.cutPlusSkin2                                   = gpu->sim.cutPlusSkin * gpu->sim.cutPlusSkin;    

    // Allocate skin test buffer
    gpu->pbNLbSkinTestFail                                  = new GpuBuffer<bool>(1, false, true);
    
    // Allocate cell hash
    gpu->pbNLCellHash                                       = new GpuBuffer<unsigned int>(CELLHASHCELLS);
    
    // Fill cell hash array
    for (int i = 0; i < CELLHASHCELLS; i++)
        gpu->pbNLCellHash->_pSysData[i]                     = cellHash[i];
    gpu->pbNLCellHash->Upload();

    // Allocate new exclusion lists
    gpu->pbNLExclusionStartCount                            = new GpuBuffer<uint2>(gpu->sim.atoms);
  
    // Double count exclusions
    unsigned int* pExclusionCount                           = new unsigned int[gpu->sim.atoms];   
    memset(pExclusionCount, 0, gpu->sim.atoms * sizeof(unsigned int));    
 
    int offset                                              = 0;
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        for (int j = 0; j < numex[i]; j++)
        {
            if (natex[offset + j] > 0)
            {
                pExclusionCount[i]++;
                pExclusionCount[natex[offset + j] - 1]++;
            }
        }                
        offset                                             += numex[i];
    }

    // Count total exclusions
    int exclusions = 0;
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        gpu->pbNLExclusionStartCount->_pSysData[i].x        = exclusions;
        gpu->pbNLExclusionStartCount->_pSysData[i].y        = 0;
        exclusions                                         += ((pExclusionCount[i] + (GRID - 1)) >> GRIDBITS) << GRIDBITS;
    }
    gpu->pbNLExclusionList                                  = new GpuBuffer<unsigned int>(exclusions);
    memset(gpu->pbNLExclusionList->_pSysData, 0xff, exclusions * sizeof(unsigned int));
    offset                                                  = 0;
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        int pos_i                                           = gpu->pbNLExclusionStartCount->_pSysData[i].x + 
                                                              gpu->pbNLExclusionStartCount->_pSysData[i].y;
        for (int j = 0; j < numex[i]; j++)
        {
            if (natex[offset + j] > 0)
            {
                gpu->pbNLExclusionList->_pSysData[pos_i++]  = natex[offset + j] - 1;
                gpu->pbNLExclusionStartCount->_pSysData[i].y++;
                int pos_j                                   = gpu->pbNLExclusionStartCount->_pSysData[natex[offset + j] - 1].x + 
                                                              gpu->pbNLExclusionStartCount->_pSysData[natex[offset + j] - 1].y++;
                gpu->pbNLExclusionList->_pSysData[pos_j]    = i;
            }
        }
        offset                                             += numex[i];
    }
    delete[] pExclusionCount;
#if 0   
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        printf("%6d %6d %6d ", i, gpu->pbNLExclusionStartCount->_pSysData[i].x, gpu->pbNLExclusionStartCount->_pSysData[i].y);
        for (int j = 0; j < gpu->pbNLExclusionStartCount->_pSysData[i].y; j++)
            printf("%6d ", gpu->pbNLExclusionList->_pSysData[gpu->pbNLExclusionStartCount->_pSysData[i].x + j]);
        printf("\n");
    }
    exit(-1);
#endif 
    gpu->pbNLExclusionStartCount->Upload();    
    gpu->pbNLExclusionList->Upload();
    gpu->sim.pNLExclusionList                               = gpu->pbNLExclusionList->_pDevData;
    gpu->sim.pNLExclusionStartCount                         = gpu->pbNLExclusionStartCount->_pDevData;

    // Calculate nonbond cell size
    PMEFloat cell                                           = gpu->sim.cut + gpu->sim.skinnb;
    gpu->sim.xcells                                         = int((gpu->sim.a) / (cell * gpu->sim.cut_factor[0]));
    if (gpu->sim.xcells < 1)
        gpu->sim.xcells                                     = 1;
    if ((gpu->sim.xcells > 20) && (gpu->sim.xcells & 1))
        gpu->sim.xcells--;        
    gpu->sim.ycells                                         = int((gpu->sim.b) / (cell * gpu->sim.cut_factor[1]));
    if (gpu->sim.ycells < 1)
        gpu->sim.ycells                                     = 1;
    if ((gpu->sim.ycells > 20) && (gpu->sim.ycells & 1))
        gpu->sim.ycells--;        
    gpu->sim.zcells                                         = int((gpu->sim.c) / (cell * gpu->sim.cut_factor[2]));
    if (gpu->sim.zcells < 1)
        gpu->sim.zcells                                     = 1; 
    if ((gpu->sim.zcells > 20) && (gpu->sim.zcells & 1))
        gpu->sim.zcells--;        
    gpu->sim.xcellsminusone                                 = gpu->sim.xcells - 1;
    gpu->sim.ycellsminusone                                 = gpu->sim.ycells - 1;
    gpu->sim.zcellsminusone                                 = gpu->sim.zcells - 1;
    if (gpu->sim.xcells & 0x1)
        gpu->sim.bOddXCells                                 = true;
    else
        gpu->sim.bOddXCells                                 = false;
     if (gpu->sim.ycells & 0x1)
        gpu->sim.bOddYCells                                 = true;
    else
        gpu->sim.bOddYCells                                 = false;
    if (gpu->sim.zcells & 0x1)
        gpu->sim.bOddZCells                                 = true;
    else
        gpu->sim.bOddZCells                                 = false;     
    if (gpu->sim.bOddXCells || gpu->sim.bOddYCells || gpu->sim.bOddZCells)
        gpu->bOddNLCells                                    = true;
    else
        gpu->bOddNLCells                                    = false;  
    gpu->sim.xycells                                        = gpu->sim.xcells * gpu->sim.ycells;
    gpu->sim.xcell                                          = (PMEFloat)(gpu->sim.a / (PMEDouble)gpu->sim.xcells);
    gpu->sim.ycell                                          = (PMEFloat)(gpu->sim.b / (PMEDouble)gpu->sim.ycells);
    gpu->sim.zcell                                          = (PMEFloat)(gpu->sim.c / (PMEDouble)gpu->sim.zcells);   
    gpu->sim.minCellX                                       = -0.5 / gpu->sim.xcells;
    gpu->sim.minCellY                                       = -0.5 / gpu->sim.ycells;
    gpu->sim.minCellZ                                       = -0.5 / gpu->sim.zcells;  
    gpu->sim.maxCellX                                       = 1.5 / gpu->sim.xcells;
    gpu->sim.maxCellY                                       = 1.5 / gpu->sim.ycells;
    gpu->sim.maxCellZ                                       = 1.5 / gpu->sim.zcells;      
    gpu->sim.oneOverXcellsf                                 = (PMEFloat)1.0 / (PMEFloat)(gpu->sim.xcells);
    gpu->sim.oneOverYcellsf                                 = (PMEFloat)1.0 / (PMEFloat)(gpu->sim.ycells);
    gpu->sim.oneOverZcellsf                                 = (PMEFloat)1.0 / (PMEFloat)(gpu->sim.zcells);      
    gpu->sim.oneOverXcells                                  = gpu->sim.oneOverXcellsf;
    gpu->sim.oneOverYcells                                  = gpu->sim.oneOverYcellsf;
    gpu->sim.oneOverZcells                                  = gpu->sim.oneOverZcellsf;  
    
    // Test for small nonbond cell count on any dimension
    if ((gpu->sim.xcells <= 2) || (gpu->sim.ycells <= 2) || (gpu->sim.zcells <= 2))
        gpu->bSmallBox                                      = true;
    else
        gpu->bSmallBox                                      = false;

    gpu->sim.cell                                           = gpu->sim.xcell;
    if (gpu->sim.ycell > gpu->sim.cell)
        gpu->sim.cell                                       = gpu->sim.ycell;
    if (gpu->sim.zcell > gpu->sim.cell)
        gpu->sim.cell                                       = gpu->sim.zcell;
    PMEFloat xskin                                          = (gpu->sim.xcell - gpu->sim.cut) / gpu->sim.cut_factor[0];
    PMEFloat yskin                                          = (gpu->sim.ycell - gpu->sim.cut) / gpu->sim.cut_factor[1];
    PMEFloat zskin                                          = (gpu->sim.zcell - gpu->sim.cut) / gpu->sim.cut_factor[2];
    gpu->sim.nonbond_skin                                   = xskin;
    if (yskin < gpu->sim.nonbond_skin)
        gpu->sim.nonbond_skin                               = yskin;
    if (zskin < gpu->sim.nonbond_skin)
        gpu->sim.nonbond_skin                               = zskin;
    gpu->sim.one_half_nonbond_skin_squared                  = (gpu->sim.nonbond_skin * gpu->sim.nonbond_skin) * 0.25;
    gpu->sim.cutPlusSkin                                    = gpu->sim.cut + gpu->sim.nonbond_skin;
    gpu->sim.cutPlusSkin2                                   = gpu->sim.cutPlusSkin * gpu->sim.cutPlusSkin;
    gpu->sim.cells                                          = gpu->sim.xcells * gpu->sim.ycells * gpu->sim.zcells;
    gpu->pbNLNonbondCellStartEnd                            = new GpuBuffer<uint2>(gpu->sim.cells);
    gpu->pbNLChargeGridBufferOffset                         = new GpuBuffer<int>(gpu->sim.cells);
    gpu->neighborListBits                                   = int(log((double)gpu->sim.cells) / log((double)2) + 1 + CELLHASHBITS);
        
    // Calculate local cell offsets
    for (int i = 0; i < NEIGHBORCELLS; i++)
    {
        if ((gpu->sim.is_orthog) && (gpu->sim.ntp == 0))
        {
            gpu->sim.cellOffset[i][0]                       = (PMEDouble)cellOffset[i][0] * gpu->sim.xcell;
            gpu->sim.cellOffset[i][1]                       = (PMEDouble)cellOffset[i][1] * gpu->sim.ycell;
            gpu->sim.cellOffset[i][2]                       = (PMEDouble)cellOffset[i][2] * gpu->sim.zcell;         
        }
        else
        {
            gpu->sim.cellOffset[i][0]                       = (PMEDouble)cellOffset[i][0] / (PMEDouble)gpu->sim.xcells;
            gpu->sim.cellOffset[i][1]                       = (PMEDouble)cellOffset[i][1] / (PMEDouble)gpu->sim.ycells;
            gpu->sim.cellOffset[i][2]                       = (PMEDouble)cellOffset[i][2] / (PMEDouble)gpu->sim.zcells;       
        }
       // printf("%3d %16.9f %16.9f %16.9f\n", i, gpu->sim.cellOffset[i][0], gpu->sim.cellOffset[i][1], gpu->sim.cellOffset[i][2]);
    }
    //exit(-1);


#if 0
    printf("cut %16.7f, skin %16.7f\n", gpu->sim.cut, gpu->sim.skinnb);
    printf("%16.7f %16.7f %16.7f %16.7f\n", gpu->sim.a, gpu->sim.b, gpu->sim.c, gpu->sim.nonbond_skin);
    printf("%16.7f %16.7f %16.7f %16.7f\n", gpu->sim.xcell, gpu->sim.ycell, gpu->sim.zcell, gpu->sim.cell);
    printf("%16d %16d %16d %16d\n", gpu->sim.xcells, gpu->sim.ycells, gpu->sim.zcells, gpu->sim.cells);
    printf("%16.7f %16.7f %16.7f\n", gpu->sim.oneOverXcells, gpu->sim.oneOverYcells, gpu->sim.oneOverZcells);
    printf("%16d %16d %16d\n", gpu->sim.nfft1, gpu->sim.nfft2, gpu->sim.nfft3);
    exit(-1);
#endif  

    // Calculate per axis odd axis counter for charge grid buffer clearing and reduction
    gpu->pbNLOddBufferOverlapFlag                           = new GpuBuffer<int>(gpu->sim.nfft1 + gpu->sim.nfft2 + gpu->sim.nfft3);
    int* pOddXBufferOverlapFlag                             = gpu->pbNLOddBufferOverlapFlag->_pSysData;
    int* pOddYBufferOverlapFlag                             = gpu->pbNLOddBufferOverlapFlag->_pSysData + gpu->sim.nfft1;
    int* pOddZBufferOverlapFlag                             = gpu->pbNLOddBufferOverlapFlag->_pSysData + gpu->sim.nfft1 + gpu->sim.nfft2;
    for (int i = 0; i < gpu->sim.nfft1 + gpu->sim.nfft2 + gpu->sim.nfft3; i++)
        gpu->pbNLOddBufferOverlapFlag->_pSysData[i]         = 0;
 
    if (gpu->sim.bOddXCells)
    {
        int ix                                              = (((gpu->sim.a * gpu->sim.xcellsminusone) / gpu->sim.xcells - gpu->sim.nonbond_skin / 2) / gpu->sim.a) * gpu->sim.nfft1 - 3 - 1;
        int ixend                                           = (gpu->sim.nonbond_skin / 2) / gpu->sim.a * gpu->sim.nfft1 + 1;
        //printf("x: %d %d\n", ix, ixend);
        while (ix != ixend)
        {
            pOddXBufferOverlapFlag[ix]                      = 1;
            ix++;
            if (ix >= gpu->sim.nfft1)
                ix                                         -= gpu->sim.nfft1;
        }
    }
    
    if (gpu->sim.bOddYCells)
    {
        int iy                                              = (((gpu->sim.b * gpu->sim.ycellsminusone) / gpu->sim.ycells - gpu->sim.nonbond_skin / 2) / gpu->sim.b) * gpu->sim.nfft2 - 3 - 1;
        int iyend                                           = (gpu->sim.nonbond_skin / 2) / gpu->sim.b * gpu->sim.nfft2 + 1;
        //printf("y: %d %d\n", iy, iyend);
        while (iy != iyend)
        {
            pOddYBufferOverlapFlag[iy]                      = 1;
            iy++;
            if (iy >= gpu->sim.nfft2)
                iy                                         -= gpu->sim.nfft2;
        }
    }
    
    if (gpu->sim.bOddZCells)
    {
        int iz                                              = (((gpu->sim.c * gpu->sim.zcellsminusone) / gpu->sim.zcells - gpu->sim.nonbond_skin / 2) / gpu->sim.c) * gpu->sim.nfft3 - 3 - 1;
        int izend                                           = (gpu->sim.nonbond_skin / 2) / gpu->sim.c * gpu->sim.nfft3 + 1;
        //printf("z: %d %d\n", iz, izend);
        while (iz != izend)
        {
            pOddZBufferOverlapFlag[iz]                      = 1;
            iz++;
            if (iz >= gpu->sim.nfft3)
                iz                                         -= gpu->sim.nfft3;
        }
    }
    gpu->pbNLOddBufferOverlapFlag->Upload();

    // Calculate total charge grid buffers
    int totalOddAxes = 0;
    if (gpu->sim.bOddXCells)
        totalOddAxes++;
    if (gpu->sim.bOddYCells)
        totalOddAxes++;
    if (gpu->sim.bOddZCells)
        totalOddAxes++;
    int singleAxisBuffers                                   = totalOddAxes * 4;
    int doubleAxisBuffers                                   = 0;
    int tripleAxisBuffers                                   = 0;
    switch (totalOddAxes)
    {
        case 2:
            doubleAxisBuffers                               = 2;
            break;
            
        case 3:
            doubleAxisBuffers                               = 6;
            tripleAxisBuffers                               = 1;
    }
    int singleAxisBufferOffset                              = 8;
    int doubleAxisBufferOffset                              = singleAxisBufferOffset + singleAxisBuffers;
    int tripleAxisBufferOffset                              = doubleAxisBufferOffset + doubleAxisBuffers;
    int totalChargeGridBuffers                              = 8 + singleAxisBuffers + doubleAxisBuffers + tripleAxisBuffers;
    int singleBuffer                                        = singleAxisBufferOffset;
    int doubleBuffer                                        = doubleAxisBufferOffset;
    
    // Calculate extra charge grid buffers when clearing or reducing a charge grid cell with odd nonbond cell overlap
    gpu->sim.extraChargeGridBuffers[0]                      = 0;
    gpu->sim.extraChargeGridBuffers[1]                      = singleAxisBuffers;
    gpu->sim.extraChargeGridBuffers[2]                      = singleAxisBuffers + doubleAxisBuffers;
    gpu->sim.extraChargeGridBuffers[3]                      = singleAxisBuffers + doubleAxisBuffers + tripleAxisBuffers;
    gpu->sim.maxChargeGridBuffers                           = 8 + singleAxisBuffers + doubleAxisBuffers + tripleAxisBuffers;

#if 0    
    printf("%d %d %d %d\n", gpu->sim.extraChargeGridBuffers[0], gpu->sim.extraChargeGridBuffers[1], gpu->sim.extraChargeGridBuffers[2], gpu->sim.extraChargeGridBuffers[3]);
    exit(-1);
#endif
    
    int bufferOffset[8][8];
    // Normal buffer offsets at zero
    bufferOffset[0][0]                                      = 0;
    bufferOffset[0][1]                                      = 1;
    bufferOffset[0][2]                                      = 2;
    bufferOffset[0][3]                                      = 3;
    bufferOffset[0][4]                                      = 4;
    bufferOffset[0][5]                                      = 5;
    bufferOffset[0][6]                                      = 6;
    bufferOffset[0][7]                                      = 7;
    
    // Odd x axis
    bufferOffset[1][0]                                      = singleBuffer;
    bufferOffset[1][1]                                      = singleBuffer;
    bufferOffset[1][2]                                      = singleBuffer + 1;
    bufferOffset[1][3]                                      = singleBuffer + 1;
    bufferOffset[1][4]                                      = singleBuffer + 2;
    bufferOffset[1][5]                                      = singleBuffer + 2;
    bufferOffset[1][6]                                      = singleBuffer + 3;
    bufferOffset[1][7]                                      = singleBuffer + 3;
    if (gpu->sim.bOddXCells)
        singleBuffer                                       += 4;
    
    // Odd y axis
    bufferOffset[2][0]                                      = singleBuffer;
    bufferOffset[2][1]                                      = singleBuffer + 1;
    bufferOffset[2][2]                                      = singleBuffer;
    bufferOffset[2][3]                                      = singleBuffer + 1;
    bufferOffset[2][4]                                      = singleBuffer + 2;
    bufferOffset[2][5]                                      = singleBuffer + 3;
    bufferOffset[2][6]                                      = singleBuffer + 2;
    bufferOffset[2][7]                                      = singleBuffer + 3;
    if (gpu->sim.bOddYCells)
        singleBuffer                                       += 4;
    
    // Odd x and y axis
    bufferOffset[3][0]                                      = doubleBuffer;
    bufferOffset[3][1]                                      = doubleBuffer;
    bufferOffset[3][2]                                      = doubleBuffer;
    bufferOffset[3][3]                                      = doubleBuffer;
    bufferOffset[3][4]                                      = doubleBuffer + 1;
    bufferOffset[3][5]                                      = doubleBuffer + 1;
    bufferOffset[3][6]                                      = doubleBuffer + 1;
    bufferOffset[3][7]                                      = doubleBuffer + 1;
    if (gpu->sim.bOddXCells && gpu->sim.bOddYCells)
        doubleBuffer                                       += 2;
    
    // Odd z axis
    bufferOffset[4][0]                                      = singleBuffer;
    bufferOffset[4][1]                                      = singleBuffer + 1;
    bufferOffset[4][2]                                      = singleBuffer + 2;
    bufferOffset[4][3]                                      = singleBuffer + 3;
    bufferOffset[4][4]                                      = singleBuffer;
    bufferOffset[4][5]                                      = singleBuffer + 1;
    bufferOffset[4][6]                                      = singleBuffer + 2;
    bufferOffset[4][7]                                      = singleBuffer + 3;
    if (gpu->sim.bOddZCells)
        singleBuffer                                       += 4;
        
    // Odd x and z axis
    bufferOffset[5][0]                                      = doubleBuffer;
    bufferOffset[5][1]                                      = doubleBuffer;
    bufferOffset[5][2]                                      = doubleBuffer + 1;
    bufferOffset[5][3]                                      = doubleBuffer + 1;
    bufferOffset[5][4]                                      = doubleBuffer;
    bufferOffset[5][5]                                      = doubleBuffer;
    bufferOffset[5][6]                                      = doubleBuffer + 1;
    bufferOffset[5][7]                                      = doubleBuffer + 1;
    if (gpu->sim.bOddXCells && gpu->sim.bOddZCells)
        doubleBuffer                                       += 2;
        
    // Odd y and z axis
    bufferOffset[6][0]                                      = doubleBuffer;
    bufferOffset[6][1]                                      = doubleBuffer + 1;
    bufferOffset[6][2]                                      = doubleBuffer;
    bufferOffset[6][3]                                      = doubleBuffer + 1;
    bufferOffset[6][4]                                      = doubleBuffer;
    bufferOffset[6][5]                                      = doubleBuffer + 1;
    bufferOffset[6][6]                                      = doubleBuffer;
    bufferOffset[6][7]                                      = doubleBuffer + 1;
    if (gpu->sim.bOddYCells && gpu->sim.bOddZCells)
        doubleBuffer                                       += 2;
    
    // Only one triple axis buffer
    bufferOffset[7][0]                                      = tripleAxisBufferOffset;
    bufferOffset[7][1]                                      = tripleAxisBufferOffset;
    bufferOffset[7][2]                                      = tripleAxisBufferOffset;
    bufferOffset[7][3]                                      = tripleAxisBufferOffset;
    bufferOffset[7][4]                                      = tripleAxisBufferOffset;
    bufferOffset[7][5]                                      = tripleAxisBufferOffset;
    bufferOffset[7][6]                                      = tripleAxisBufferOffset;
    bufferOffset[7][7]                                      = tripleAxisBufferOffset;


#ifdef MPI
    bool* bCellReceive                                      = new bool[gpu->sim.cells];
    bool* bCellSend                                         = new bool[gpu->sim.cells];
    bool* bSendGroup                                        = new bool[gpu->nGpus];
    bool* bReceiveGroup                                     = new bool[gpu->nGpus];
    
    int* cellOwner                                          = new int[gpu->sim.cells];
    
    // Really braindead but surprisingly accurate load-balancing
    unsigned int cellweight                                 = 10000 / gpu->nGpus;
    unsigned int pmeweight                                  = 2540;
    if (gpu->nGpus >= 4)
        pmeweight                                           = 2540;
    else if (gpu->sim.atoms > 60000)
        pmeweight                                           = 2390;
    unsigned int totalweight;
    if (pmeweight < cellweight)
    {
        totalweight                                         = gpu->nGpus * cellweight - pmeweight;
    }
    else
    {
        totalweight                                         = (gpu->nGpus - 1) * cellweight;
    }

    // Determine global cell allocation
    gpu->pMinLocalCell                                      = new int[gpu->nGpus];
    gpu->pMaxLocalCell                                      = new int[gpu->nGpus];
    gpu->pMinLocalAtom                                      = new int[gpu->nGpus];
    gpu->pMaxLocalAtom                                      = new int[gpu->nGpus];    
    gpu->pMinProcessedCell                                  = new int[gpu->nGpus];
    gpu->pMaxProcessedCell                                  = new int[gpu->nGpus];
    gpu->pMinProcessedAtom                                  = new int[gpu->nGpus];
    gpu->pMaxProcessedAtom                                  = new int[gpu->nGpus];
    gpu->pAllGathervRecvCountAoS                            = new int[gpu->nGpus];
    gpu->pAllGathervRecvDisplAoS                            = new int[gpu->nGpus];
    memset(bCellReceive, 0, gpu->sim.cells * sizeof(bool)); 
    memset(bCellSend, 0, gpu->sim.cells * sizeof(bool));
    
    for (int i = 0; i < gpu->nGpus; i++)
    {
        int minCell                                         = 0;
        int maxCell                                         = 0;
        // Hardcode writes from gpu 0 on gpu 0 and reads from it on all other nodes for PME, skip for IPS
        if (gpu->sim.bIPSActive)
        {
            minCell                                         = (i * gpu->sim.cells) / gpu->nGpus;
            maxCell                                         = ((i + 1) * gpu->sim.cells) / gpu->nGpus;
        }
        else
        {
            if (gpu->gpuID != 0)
            {
                bSendGroup[i]                               = false;
                bReceiveGroup[i]                            = (i == 0);
            }
            else 
            {
                bSendGroup[i]                               = true;
                bReceiveGroup[i]                            = false;            
            }


            if (i == 0)
            {
                if (pmeweight < cellweight)
                    maxCell                                     = (cellweight - pmeweight) * gpu->sim.cells / totalweight;
            }
            else if (pmeweight < cellweight)
            {
                minCell                                         = (i * cellweight - pmeweight) * gpu->sim.cells / totalweight;
                maxCell                                         = ((i + 1) * cellweight - pmeweight) * gpu->sim.cells / totalweight;
            }
            else
            {
                minCell                                         = ((i - 1) * cellweight) * gpu->sim.cells / totalweight;
                maxCell                                         = (i * cellweight) * gpu->sim.cells / totalweight;    
            }
        }
        for (int j = minCell; j < maxCell; j++)
        {
            cellOwner[j]                                    = i; 
        }
        gpu->pMinLocalCell[i]                               = minCell;
        gpu->pMaxLocalCell[i]                               = maxCell;
        gpu->pMinLocalAtom[i]                               = 0;
        gpu->pMaxLocalAtom[i]                               = 0;
        gpu->pMinProcessedCell[i]                           = minCell;
        gpu->pMaxProcessedCell[i]                           = maxCell;
        gpu->pMinProcessedAtom[i]                           = 0;
        gpu->pMaxProcessedAtom[i]                           = 0;
        gpu->pAllGathervRecvCountAoS[i]                     = 0;
        gpu->pAllGathervRecvDisplAoS[i]                     = 0;
    }
    
    // Determine local cell allocation
    gpu->minLocalCell                                       = gpu->pMinProcessedCell[gpu->gpuID];
    gpu->maxLocalCell                                       = max(0, gpu->pMaxProcessedCell[gpu->gpuID]);
    gpu->minProcessedCell                                   = gpu->minLocalCell;
    gpu->maxProcessedCell                                   = gpu->maxLocalCell;
#ifdef GVERBOSE    
    printf("Node %d %d %d\n", gpu->gpuID, gpu->minLocalCell, gpu->maxLocalCell);
#endif
#endif

    // Build work unit list
    int workUnits                                           = 0;
    NLRecord* pNLCellRecord                                 = new NLRecord[gpu->sim.cells];
    unsigned int homeCells                                  = 0;
    
    int offUnits                                            = 0;
    for (int i = 0; i < gpu->sim.cells; i++)
    {
        int x                                               = i % gpu->sim.xcells;
        int y                                               = ((i - x) / gpu->sim.xcells) % gpu->sim.ycells;
        int z                                               = (i - x - y * gpu->sim.xcells) / (gpu->sim.xcells * gpu->sim.ycells);
              
        // Clear neighbor list data
        pNLCellRecord[i].NL.homeCell                        = i;
        pNLCellRecord[i].NL.neighborCells                   = 0;
              
        // Calculate buffer offset;
        int buffer                                          = 0;
        if (gpu->sim.bOddXCells && (x == gpu->sim.xcells - 1))
            buffer                                         |= 1;
        if (gpu->sim.bOddYCells && (y == gpu->sim.ycells - 1))
            buffer                                         |= 2;
        if (gpu->sim.bOddZCells && (z == gpu->sim.zcells - 1))
            buffer                                         |= 4;
        int offset                                          = (4 * (z & 0x1) + 2 * (y & 0x1) + (x & 0x1));
        gpu->pbNLChargeGridBufferOffset->_pSysData[i]       = bufferOffset[buffer][offset] * gpu->sim.XYZStride;
        
        for (int j = 0; j < NEIGHBORCELLS; j++)
        {
            int x1                                          = x + cellOffset[j][0];
            int y1                                          = y + cellOffset[j][1];
            int z1                                          = z + cellOffset[j][2];
            if (x1 < 0)
                x1                                         += gpu->sim.xcells;
            if (y1 < 0)
                y1                                         += gpu->sim.ycells;
            if (z1 < 0)
                z1                                         += gpu->sim.zcells;      
            if (x1 >= gpu->sim.xcells)
                x1                                         -= gpu->sim.xcells;
            if (y1 >= gpu->sim.ycells)
                y1                                         -= gpu->sim.ycells;
            if (z1 >= gpu->sim.zcells)
                z1                                         -= gpu->sim.zcells;  
            int cell1                                       = (z1 * gpu->sim.ycells + y1) * gpu->sim.xcells + x1;  

#ifdef MPI
            int cell2                                       = cell1;
            if ((z1 == 0) && (z != 0))           
                cell2                                      += gpu->sim.cells;
            if ((i <= cell1) || (z1 != z))
            {
                if (cell2 < gpu->pMinProcessedCell[cellOwner[i]])
                    gpu->pMinProcessedCell[cellOwner[i]]    = cell2;
                if (cell2 >= gpu->pMaxProcessedCell[cellOwner[i]])
                    gpu->pMaxProcessedCell[cellOwner[i]]    = cell2 + 1; 
            }
            else
            {
                if (i < gpu->pMinProcessedCell[cellOwner[cell1]])
                    gpu->pMinProcessedCell[cellOwner[cell1]]= cell2;
                if (i >= gpu->pMaxProcessedCell[cellOwner[cell1]])
                    gpu->pMaxProcessedCell[cellOwner[cell1]]= cell2 + 1; 
            }
            
            
			if (((cellOwner[i] == gpu->gpuID) && ((i <= cell1) || (z1 != z))) ||
			    ((cellOwner[cell1] == gpu->gpuID) && (i > cell1) && (z1 == z)))
			
			{
#endif
                workUnits++;
                
                // Add to neighbor list records
                if (pNLCellRecord[i].NL.neighborCells == 0)
                    homeCells++;
                pNLCellRecord[i].NL.neighborCell[pNLCellRecord[i].NL.neighborCells++]    = (cell1 << NLCELLSHIFT) | j;
#ifdef MPI
                // Add to cell write flags if updating non-local cell               
                if (cellOwner[i] != gpu->gpuID)
                {
                    bCellSend[i]                            = true;
                    bSendGroup[cellOwner[i]]                = true;
                    offUnits++;
                }
                if (cellOwner[cell1] != gpu->gpuID)
                {
                    bCellSend[cell1]                        = true;   
                    bSendGroup[cellOwner[cell1]]            = true;
                    offUnits++;
                }                 
			}
			else 
			{
			    // Update cell read flags
			    if (cellOwner[cell1] == gpu->gpuID)
			    {
    			    bCellReceive[cell1]                     = true;
    			    bReceiveGroup[cellOwner[i]]             = true;
    			}
    			if (cellOwner[i] == gpu->gpuID)
    			{
    			    bCellReceive[i]                         = true;
    			    bReceiveGroup[cellOwner[cell1]]         = true;
    			}
			}
#endif 
        }
    }

    //printf("%d home cells\n", homeCells);
    // Determine cell divisors for neighbor list generation
    unsigned int totalWarps                                 = (gpu->PMENonbondEnergyThreadsPerBlock * gpu->blocks) / GRID;
    unsigned int xDivisor                                   = 1;
    unsigned int yDivisor                                   = 1;
    unsigned int atomsPerCell                               = (gpu->sim.atoms + GRID) / gpu->sim.cells;
    unsigned int atomsPerWarp                               = 32;
    unsigned int maxExclusionsPerWarp                       = 570; // 295 (768 threads), 322 (704 threads), 338 (672 threads) 
                                                                   // 510 (606 threads), 538 (574 threads), 570 (544 threads), 606 (512 threads), 646 (480 threads)

    // Make sure GPU is busy by adjusting Y divisor
    if (homeCells < 8)
        yDivisor                                            = 8;
    else if (homeCells < totalWarps / 8)
        yDivisor                                            = 4;
    else if (homeCells < totalWarps * 3)
        yDivisor                                            = 3;        
    else if (homeCells < totalWarps * 6)
        yDivisor                                            = 2; 
    else
        yDivisor                                            = 1;
    
   
    // Adjust X divisor further for really really small cell counts and high multi-GPU node counts
    if (homeCells < 3 * totalWarps / 5)
    {
        if (homeCells < 9)
            xDivisor                                        = 14;
        else if (homeCells < 27)
            xDivisor                                        = 7;
        else if (homeCells < totalWarps / 2)
            xDivisor                                        = 5;
        else
            xDivisor                                        = 2;
    }

    // Adjust warp atoms based on cell volume
    if (atomsPerCell < 256)
        atomsPerWarp                                        = 16;

    
    // Determine exclusions per warp
    if (atomsPerWarp == 32)
        maxExclusionsPerWarp                                = 570;
    else
        maxExclusionsPerWarp                                = 338;
       
        
    // Allocate and fill buffer cell lookup table
    gpu->sim.NLCellBuffers                                  = yDivisor * NEIGHBORCELLS + xDivisor;
    gpu->sim.NLHomeCellBuffer                               = yDivisor * NEIGHBORCELLS;
 
    // Count neighbor list records
    unsigned int NLRecords                                  = 0;

    for (int i = 0; i < gpu->sim.cells; i++)
    {
        if (pNLCellRecord[i].NL.neighborCells != 0)
            NLRecords                                          += yDivisor * min(pNLCellRecord[i].NL.neighborCells, xDivisor);
    }
    //printf("%d neighbor list records %d %d %d %d %d\n", NLRecords, gpu->sim.xcells, gpu->sim.ycells, gpu->sim.zcells, yDivisor, xDivisor);
    if (atomsPerWarp == 32)
        gpu->sim.NLBuildWarps                               = (gpu->blocks * gpu->NLBuildNeighborList32ThreadsPerBlock) / GRID;
    else
        gpu->sim.NLBuildWarps                               = (gpu->blocks * gpu->NLBuildNeighborList16ThreadsPerBlock) / GRID;
    gpu->sim.NLNonbondWarps                                 = totalWarps;
   
    gpu->pbNLRecord                                         = new GpuBuffer<NLRecord>(NLRecords);
    gpu->pbNLEntry                                          = new GpuBuffer<NLEntry>(NLRecords);
    gpu->pbNLOffset                                         = new GpuBuffer<unsigned int>(NLRecords);
    gpu->sim.NLSize                                         = NLRecords;
    gpu->sim.NLXDivisor                                     = xDivisor;    
    gpu->sim.NLYDivisor                                     = yDivisor;
    SetNLClearForcesKernel(gpu);    
    SetNLReduceForcesKernel(gpu);
    gpu->sim.NLEntryTypes                                   = yDivisor; // Placeholder, not redundant upon the return of xDivisor
    gpu->sim.NLYStride                                      = yDivisor * atomsPerWarp;
    gpu->sim.NLAtomsPerWarp                                 = atomsPerWarp;
    gpu->sim.NLAtomsPerWarpBits                             = ffs(atomsPerWarp) - 1;
    gpu->sim.NLAtomsPerWarpBitsMask                         = (1 << gpu->sim.NLAtomsPerWarpBits) - 1;
    gpu->sim.NLAtomsPerWarpMask                             = (unsigned long)((1ull << (1ull << (gpu->sim.NLAtomsPerWarpBits))) - 1ull);
    gpu->sim.NLMaxExclusionsPerWarp                         = maxExclusionsPerWarp;
    NLRecord* pNLRecord                                     = gpu->pbNLRecord->_pSysData;
    gpu->sim.pNLRecord                                      = gpu->pbNLRecord->_pDevData;
    gpu->sim.pNLEntry                                       = gpu->pbNLEntry->_pDevData;
    gpu->sim.pNLOffset                                      = gpu->pbNLOffset->_pDevData;
    NLRecords                                               = 0;

    // Build Neighbor list records from raw neighbor cell data
    for (int i = 0; i < gpu->sim.cells; i++)
    {
        if (pNLCellRecord[i].NL.neighborCells != 0)
        {
            int divisor                                    = min(xDivisor, pNLCellRecord[i].NL.neighborCells);
            for (int y = 0; y < yDivisor; y++)
            {

                for (int x = 0; x < divisor; x++)
                {
                    int start                               = (x       * pNLCellRecord[i].NL.neighborCells) / divisor;
                    int end                                 = ((x + 1) * pNLCellRecord[i].NL.neighborCells) / divisor;
                    if (end > start)
                    {
                        pNLRecord[NLRecords].NL.homeCell    = i;                    
                        pNLRecord[NLRecords].NL.neighborCells  
                                                            = ((end - start) << NLCELLCOUNTSHIFT) | (x << NLXCELLOFFSETSHIFT) | y;  
                      //  printf("A %06d %06d %06d\n", NLRecords, y, x);              
                        for (int k = start; k < end; k++)
                            pNLRecord[NLRecords].NL.neighborCell[k - start] 
                                                            = pNLCellRecord[i].NL.neighborCell[k];
                       // printf("%06d %06d %06d %06d %06d\n", NLRecords, y, x, start, end);
#if 0
                printf("%06d %3d %3d: ", pNLRecord[NLRecords].NL.homeCell,
                                         pNLRecord[NLRecords].NL.neighborCells >> NLCELLCOUNTSHIFT,
                                        (pNLRecord[NLRecords].NL.neighborCells & NLATOMOFFSETMASK) * gpu->sim.NLAtomsPerWarp
                                              );
                for (int l = 0; l < (pNLRecord[NLRecords].NL.neighborCells >> NLCELLCOUNTSHIFT); l++)
                    printf("%4d %3d ", pNLRecord[NLRecords].NL.neighborCell[l] >> NLCELLSHIFT,
                                       pNLRecord[NLRecords].NL.neighborCell[l] & NLBUFFERMASK);
                printf("\n");                           
#endif
                        NLRecords++;  
                    }
                }                  
            }
        }
    }
    gpu->pbNLRecord->Upload();
//    printf("%06d %06d %06d %06d\n", gpu->sim.NLAtomsPerWarp, gpu->sim.NLAtomsPerWarpBits, yDivisor, xDivisor);
//    exit(-1); 
    

    
   // printf("WU %d %d %d\n", gpu->gpuID, workUnits, offUnits);
   // printf("WU %d %d\n", workUnits, offUnits); //3024, 5636 : 378, 1033 : 300, 770 - 100, 589 | 3024 820, 4738
    
    // Map excessive exclusion mask space
    unsigned int warpsPerCell                                       = (atomsPerCell + GRID - 1) / GRID;
    unsigned int multiplier                                         = GRID / gpu->sim.NLAtomsPerWarp;
    unsigned int offsetPerWarp                                      = GRID + gpu->sim.NLAtomsPerWarp;
    unsigned int maxTotalOffset                                     = 2 * workUnits * multiplier * warpsPerCell * warpsPerCell * offsetPerWarp;
    gpu->pbNLAtomList                                               = new GpuBuffer<unsigned int>(maxTotalOffset, bShadowedOutputBuffers); 
    gpu->sim.NLMaxTotalOffset                                       = maxTotalOffset;
    gpu->sim.NLOffsetPerWarp                                        = offsetPerWarp;
    gpu->pbNLTotalOffset                                            = new GpuBuffer<unsigned int>(1); 
    gpu->sim.pNLTotalOffset                                         = gpu->pbNLTotalOffset->_pDevData;
    gpu->pbNLPosition                                               = new GpuBuffer<unsigned int>(1); 
    gpu->sim.pNLPosition                                            = gpu->pbNLPosition->_pDevData;    
    gpu->sim.pNLAtomList                                            = gpu->pbNLAtomList->_pDevData;
    
#ifdef MPI
#ifdef GVERBOSE
    printf("Node %d, %d work units\n", gpu->gpuID, workUnits);
#endif
    for (int i = 0; i < gpu->nGpus; i++)    
    {
#ifdef GVERBOSE    
        printf("Node %d Before %d, min/max %d %d\n", gpu->gpuID, i, gpu->pMinProcessedCell[i], gpu->pMaxProcessedCell[i]);
#endif
        if (gpu->pMaxProcessedCell[i] - gpu->pMinProcessedCell[i] >= gpu->sim.cells)
        {
            gpu->pMinProcessedCell[i]                               = 0;
            gpu->pMaxProcessedCell[i]                               = gpu->sim.cells;
        }    
    }
    gpu->minProcessedCell                                           = gpu->pMinProcessedCell[gpu->gpuID];
    gpu->maxProcessedCell                                           = gpu->pMaxProcessedCell[gpu->gpuID];
#ifdef GVERBOSE
    printf("Node %d, After min/max %d %d\n", gpu->gpuID, gpu->minProcessedCell, gpu->maxProcessedCell);
#endif
    MPI_Alloc_mem(gpu->sim.stride3 * (2 * sizeof(PMEDouble) + sizeof(PMEFloat)), MPI_INFO_NULL, &gpu->pForceData);  
    MPI_Win_create(gpu->pForceData, gpu->sim.stride3 * (2 * sizeof(PMEDouble) + sizeof(PMEFloat)), 1, MPI_INFO_NULL, gpu->comm, &gpu->MPIPMEForceWindow);
    gpu->pForceData0                                                = (PMEDouble*)gpu->pForceData;
    gpu->pForceData1                                                = (PMEDouble*)gpu->pForceData + gpu->sim.stride3;
    gpu->forceSendOffset                                            = 0;
    gpu->pPMEForceData                                              = (PMEFloat*)(gpu->pForceData1 + gpu->sim.stride3);
    gpu->pPMEStart                                                  = new int[gpu->nGpus];
    gpu->pPMELength                                                 = new int[gpu->nGpus];
    memset(gpu->pForceData, 0, gpu->sim.stride3 * (2 * sizeof(PMEDouble) + sizeof(PMEFloat))); 
    MPI_Win_fence(0, gpu->MPIPMEForceWindow);    
    
    // Create read/write groups
    int* receiveGroup                                               = new int[gpu->nGpus];
    int* sendGroup                                                  = new int[gpu->nGpus];
  
    int receiveGroupSize                                            = 0;
    int sendGroupSize                                               = 0;
    for (int i = 0; i < gpu->nGpus; i++)
    {
        if (i != gpu->gpuID)
        {
            if (bReceiveGroup[i])
                receiveGroup[receiveGroupSize++]                   = i;
            if (bSendGroup[i])
                sendGroup[sendGroupSize++]                         = i;
        }
    }

    // Create force groups and open up initial window     
    MPI_Group worldGroup;
    MPI_Comm_group(gpu->comm, &worldGroup); 
    MPI_Group_incl(worldGroup, receiveGroupSize, receiveGroup, &gpu->MPIPMEForceReceiveGroup);
    MPI_Group_incl(worldGroup, sendGroupSize, sendGroup, &gpu->MPIPMEForceSendGroup);
    MPI_Group_free(&worldGroup);   
    
    // Allocate force send data
    gpu->pForceSendNode                                             = new int[sendGroupSize];
    gpu->pForceSendMinCell                                          = new int[sendGroupSize];
    gpu->pForceSendMaxCell                                          = new int[sendGroupSize];
    gpu->pOutForceSendStart                                         = new int[sendGroupSize];
    gpu->pForceSendStart                                            = new int[sendGroupSize];
    gpu->pForceSendLength                                           = new int[sendGroupSize];
    
   
    // Determine cell receive bounds
    gpu->forceReceiveFirstCell                                      = gpu->sim.cells;
    gpu->forceReceiveLastCell                                       = -1;
    gpu->forceReceiveFirstAtom                                      = 0;
    gpu->forceReceiveLastAtom                                       = 0;
#ifdef GVERBOSE
    printf("M %d %d %d\n", gpu->gpuID, gpu->minLocalCell, gpu->maxLocalCell);
#endif
    for (int i = gpu->minLocalCell; i < gpu->maxLocalCell; i++)
    {
        if (bCellReceive[i])
        {
            if (i < gpu->forceReceiveFirstCell)
                gpu->forceReceiveFirstCell                          = i;
            if (i > gpu->forceReceiveLastCell)
                gpu->forceReceiveLastCell                           = i;                
        }
    } 
   
    // Create list of writes
    int forceSendNodes                                              = 0;
    if (gpu->minLocalCell != gpu->maxLocalCell)
    {

        for (int i = 0; i < sendGroupSize; i++)
        {
            int minCell                                             = gpu->sim.cells;
            int maxCell                                             = -1;
            int node                                                = sendGroup[i];
            
            // Look for max and min cells to write
            for (int j = gpu->pMinLocalCell[node]; j < gpu->pMaxLocalCell[node]; j++)
            {
                if (bCellSend[j])
                {
                    if (j < minCell)
                        minCell                                     = j;
                    if (j > maxCell)
                        maxCell                                     = j;
                }
            }
            if (maxCell >= minCell)
            {
                gpu->pForceSendNode[forceSendNodes]                 = sendGroup[i];
                gpu->pForceSendMinCell[forceSendNodes]              = minCell;
                gpu->pForceSendMaxCell[forceSendNodes]              = maxCell;
                gpu->pOutForceSendStart[forceSendNodes]             = 0;
                gpu->pForceSendStart[forceSendNodes]                = 0;
                gpu->pForceSendLength[forceSendNodes]               = 0;
                forceSendNodes++;
#ifdef GVERBOSE
                printf("S %d %d %d %d\n", gpu->gpuID, sendGroup[i], minCell, maxCell);
#endif
            }
        }
    }
    gpu->forceSendNodes                                             = forceSendNodes;
    
    delete[] receiveGroup;
    delete[] sendGroup;
    delete[] bCellReceive;
    delete[] bCellSend;          
    delete[] bReceiveGroup;
    delete[] bSendGroup;
    delete[] cellOwner;        
#endif    

    // Determine block scheduling    
    gpu->sim.NLWorkUnits                                            = workUnits;
    gpu->pbNLChargeGridBufferOffset->Upload();


    // Allocate remapped local interactions
    gpu->pbImageBondID                                      = new GpuBuffer<int4>(gpu->sim.bonds);             
    gpu->pbImageBondAngleID1                                = new GpuBuffer<int4>(gpu->sim.bondAngles);
    gpu->pbImageBondAngleID2                                = new GpuBuffer<int2>(gpu->sim.bondAngles);
    gpu->pbImageDihedralID1                                 = new GpuBuffer<int4>(gpu->sim.dihedrals);
    gpu->pbImageDihedralID2                                 = new GpuBuffer<int4>(gpu->sim.dihedrals);
    gpu->pbImageNb14ID                                      = new GpuBuffer<int4>(gpu->sim.nb14s);
    gpu->pbImageConstraintID                                = new GpuBuffer<int2>(gpu->sim.constraints);
    gpu->pbImageUBAngleID                                   = new GpuBuffer<int4>(gpu->sim.UBAngles);   
    gpu->pbImageImpDihedralID1                              = new GpuBuffer<int4>(gpu->sim.impDihedrals);   
    gpu->pbImageImpDihedralID2                              = new GpuBuffer<int4>(gpu->sim.impDihedrals);   
    gpu->pbImageCmapID1                                     = new GpuBuffer<int4>(gpu->sim.cmaps);   
    gpu->pbImageCmapID2                                     = new GpuBuffer<int4>(gpu->sim.cmaps);   
    gpu->pbImageCmapID3                                     = new GpuBuffer<int2>(gpu->sim.cmaps);            
    gpu->pbImageShakeID                                     = new GpuBuffer<int4>(gpu->sim.shakeConstraints);
    gpu->pbImageFastShakeID                                 = new GpuBuffer<int4>(gpu->sim.fastShakeConstraints);
    gpu->pbImageSlowShakeID1                                = new GpuBuffer<int>(gpu->sim.slowShakeConstraints);
    gpu->pbImageSlowShakeID2                                = new GpuBuffer<int4>(gpu->sim.slowShakeConstraints);
    if ((gpu->sim.solventMolecules > 0) || (gpu->sim.soluteAtoms > 0))
    {
        gpu->pbImageSolventAtomID                           = new GpuBuffer<int4>(gpu->sim.solventMolecules);
        gpu->pbImageSoluteAtomID                            = new GpuBuffer<int>(gpu->sim.soluteAtoms);
    }

    // Allocate new lists
    gpu->pbImageIndex                                       = new GpuBuffer<unsigned int>(7 * gpu->sim.stride);
    gpu->pbAtomXYSaveSP                                     = new GpuBuffer<PMEFloat2>(gpu->sim.stride);
    gpu->pbAtomZSaveSP                                      = new GpuBuffer<PMEFloat>(gpu->sim.stride);
    gpu->pbImage                                            = new GpuBuffer<PMEDouble>(6 * gpu->sim.stride);
    gpu->pbImageVel                                         = new GpuBuffer<PMEDouble>(6 * gpu->sim.stride);
    gpu->pbImageLVel                                        = new GpuBuffer<PMEDouble>(6 * gpu->sim.stride);
    gpu->pbImageMass                                        = new GpuBuffer<PMEDouble>(4 * gpu->sim.stride);
    gpu->pbImageCharge                                      = new GpuBuffer<PMEDouble>(2 * gpu->sim.stride);  
    gpu->pbImageSigEps                                      = new GpuBuffer<PMEFloat2>(2 * gpu->sim.stride);
    gpu->pbImageOutputBuffers                               = new GpuBuffer<unsigned int>(2 * gpu->sim.stride);
    gpu->pbImageCellID                                      = new GpuBuffer<unsigned int>(2 * gpu->sim.stride);
    
    // Allocate new extra point arrays
    gpu->pbImageExtraPoint11Frame                           = new GpuBuffer<int4>(gpu->sim.EP11s);
    gpu->pbImageExtraPoint11Index                           = new GpuBuffer<int>(gpu->sim.EP11s);
    gpu->pbImageExtraPoint12Frame                           = new GpuBuffer<int4>(gpu->sim.EP12s);
    gpu->pbImageExtraPoint12Index                           = new GpuBuffer<int>(gpu->sim.EP12s);    
    gpu->pbImageExtraPoint21Frame                           = new GpuBuffer<int4>(gpu->sim.EP21s);
    gpu->pbImageExtraPoint21Index                           = new GpuBuffer<int2>(gpu->sim.EP21s);
    gpu->pbImageExtraPoint22Frame                           = new GpuBuffer<int4>(gpu->sim.EP22s);
    gpu->pbImageExtraPoint22Index                           = new GpuBuffer<int2>(gpu->sim.EP22s);    
  
    // Copy data
    for (int i = 0; i < gpu->sim.stride; i++)
    {
        gpu->pbImageIndex->_pSysData[i]                             = i;
        gpu->pbImageIndex->_pSysData[i + gpu->sim.stride]           = i;
        gpu->pbImageIndex->_pSysData[i + gpu->sim.stride2]          = i;
        gpu->pbImageIndex->_pSysData[i + gpu->sim.stride3]          = 0;
        gpu->pbImageIndex->_pSysData[i + gpu->sim.stride4]          = i;
        gpu->pbImageIndex->_pSysData[i + gpu->sim.stride * 5]       = i;
        gpu->pbImageIndex->_pSysData[i + gpu->sim.stride * 6]       = i;
        gpu->pbImage->_pSysData[i]                                  = gpu->pbAtom->_pSysData[i];
        gpu->pbImage->_pSysData[i + gpu->sim.stride]                = gpu->pbAtom->_pSysData[i + gpu->sim.stride];
        gpu->pbImage->_pSysData[i + gpu->sim.stride2]               = gpu->pbAtom->_pSysData[i + gpu->sim.stride2];
        gpu->pbAtomXYSaveSP->_pSysData[i].x                         = gpu->pbAtom->_pSysData[i];
        gpu->pbAtomXYSaveSP->_pSysData[i].y                         = gpu->pbAtom->_pSysData[i + gpu->sim.stride];
        gpu->pbAtomZSaveSP->_pSysData[i]                            = gpu->pbAtom->_pSysData[i + gpu->sim.stride2];
        gpu->pbImageVel->_pSysData[i]                               = gpu->pbVel->_pSysData[i];
        gpu->pbImageVel->_pSysData[i + gpu->sim.stride]             = gpu->pbVel->_pSysData[i + gpu->sim.stride];
        gpu->pbImageVel->_pSysData[i + gpu->sim.stride2]            = gpu->pbVel->_pSysData[i + gpu->sim.stride2];
        gpu->pbImageLVel->_pSysData[i]                              = gpu->pbLVel->_pSysData[i];
        gpu->pbImageLVel->_pSysData[i + gpu->sim.stride]            = gpu->pbLVel->_pSysData[i + gpu->sim.stride];
        gpu->pbImageLVel->_pSysData[i + gpu->sim.stride2]           = gpu->pbLVel->_pSysData[i + gpu->sim.stride2];
        gpu->pbImageMass->_pSysData[i]                              = gpu->pbAtomMass->_pSysData[i];
        gpu->pbImageMass->_pSysData[i + gpu->sim.stride]            = gpu->pbAtomMass->_pSysData[i + gpu->sim.stride];
        gpu->pbImageCharge->_pSysData[i]                            = gpu->pbAtomCharge->_pSysData[i];
        gpu->pbImageSigEps->_pSysData[i]                            = gpu->pbAtomSigEps->_pSysData[i];
        gpu->pbImageOutputBuffers->_pSysData[i]                     = max(gpu->sim.NLCellBuffers, gpu->pbOutputBufferCounter->_pSysData[i]) * gpu->sim.stride3;
        gpu->pbImageOutputBuffers->_pSysData[i + gpu->sim.stride]   = max(gpu->sim.NLCellBuffers, gpu->pbOutputBufferCounter->_pSysData[i]) * gpu->sim.stride3;
        gpu->pbImageCellID->_pSysData[i]                            = 0;
    }
    gpu->pbImageIndex->Upload();
    gpu->pbImage->Upload();
    gpu->pbImageVel->Upload();
    gpu->pbImageLVel->Upload();
    gpu->pbImageMass->Upload();
    gpu->pbImageCharge->Upload();
    gpu->pbImageSigEps->Upload();
    gpu->pbImageOutputBuffers->Upload();
    gpu->pbImageCellID->Upload();

    // Copy bonded interactions
    for (int i = 0; i < gpu->sim.bonds; i++)
    {
        gpu->pbImageBondID->_pSysData[i]                    = gpu->pbBondID->_pSysData[i];
        gpu->pbBondID->_pSysData[i].z                      -= gpu->pbBondID->_pSysData[i].x;
        gpu->pbBondID->_pSysData[i].w                      -= gpu->pbBondID->_pSysData[i].y;
    }
    for (int i = 0; i < gpu->sim.bondAngles; i++)
    {
        gpu->pbImageBondAngleID1->_pSysData[i]              = gpu->pbBondAngleID1->_pSysData[i];
        gpu->pbImageBondAngleID2->_pSysData[i]              = gpu->pbBondAngleID2->_pSysData[i];
        gpu->pbBondAngleID1->_pSysData[i].w                -= gpu->pbBondAngleID1->_pSysData[i].x;
        gpu->pbBondAngleID2->_pSysData[i].x                -= gpu->pbBondAngleID1->_pSysData[i].y;
        gpu->pbBondAngleID2->_pSysData[i].y                -= gpu->pbBondAngleID1->_pSysData[i].z;        
    }
    for (int i = 0; i < gpu->sim.dihedrals; i++)
    {
        gpu->pbImageDihedralID1->_pSysData[i]               = gpu->pbDihedralID1->_pSysData[i];
        gpu->pbImageDihedralID2->_pSysData[i]               = gpu->pbDihedralID2->_pSysData[i];
        gpu->pbDihedralID2->_pSysData[i].x                 -= gpu->pbDihedralID1->_pSysData[i].x;
        gpu->pbDihedralID2->_pSysData[i].y                 -= gpu->pbDihedralID1->_pSysData[i].y;
        gpu->pbDihedralID2->_pSysData[i].z                 -= gpu->pbDihedralID1->_pSysData[i].z;
        gpu->pbDihedralID2->_pSysData[i].w                 -= gpu->pbDihedralID1->_pSysData[i].w;
    }
    for (int i = 0; i < gpu->sim.nb14s; i++)
    {
        gpu->pbImageNb14ID->_pSysData[i]                    = gpu->pbNb14ID->_pSysData[i];
        gpu->pbNb14ID->_pSysData[i].z                      -= gpu->pbNb14ID->_pSysData[i].x;
        gpu->pbNb14ID->_pSysData[i].w                      -= gpu->pbNb14ID->_pSysData[i].y;
    }
    for (int i = 0; i < gpu->sim.constraints; i++)
    {
        gpu->pbImageConstraintID->_pSysData[i]              = gpu->pbConstraintID->_pSysData[i];
        gpu->pbConstraintID->_pSysData[i].y                -= gpu->pbConstraintID->_pSysData[i].x;
    }
    for (int i = 0; i < gpu->sim.UBAngles; i++)
    {
        gpu->pbImageUBAngleID->_pSysData[i]                 = gpu->pbUBAngleID->_pSysData[i];
        gpu->pbUBAngleID->_pSysData[i].z                   -= gpu->pbUBAngleID->_pSysData[i].x;
        gpu->pbUBAngleID->_pSysData[i].w                   -= gpu->pbUBAngleID->_pSysData[i].y;
    }
    for (int i = 0; i < gpu->sim.impDihedrals; i++)
    {
        gpu->pbImageImpDihedralID1->_pSysData[i]            = gpu->pbImpDihedralID1->_pSysData[i];
        gpu->pbImageImpDihedralID2->_pSysData[i]            = gpu->pbImpDihedralID2->_pSysData[i];
        gpu->pbImpDihedralID2->_pSysData[i].x              -= gpu->pbImpDihedralID1->_pSysData[i].x;
        gpu->pbImpDihedralID2->_pSysData[i].y              -= gpu->pbImpDihedralID1->_pSysData[i].y;
        gpu->pbImpDihedralID2->_pSysData[i].z              -= gpu->pbImpDihedralID1->_pSysData[i].z;
        gpu->pbImpDihedralID2->_pSysData[i].w              -= gpu->pbImpDihedralID1->_pSysData[i].w;
    }
    for (int i = 0; i < gpu->sim.cmaps; i++)
    {
        gpu->pbImageCmapID1->_pSysData[i]                   = gpu->pbCmapID1->_pSysData[i];
        gpu->pbImageCmapID2->_pSysData[i]                   = gpu->pbCmapID2->_pSysData[i];
        gpu->pbImageCmapID3->_pSysData[i]                   = gpu->pbCmapID3->_pSysData[i];
        gpu->pbCmapID2->_pSysData[i].y                     -= gpu->pbCmapID1->_pSysData[i].x;
        gpu->pbCmapID2->_pSysData[i].z                     -= gpu->pbCmapID1->_pSysData[i].y;
        gpu->pbCmapID2->_pSysData[i].w                     -= gpu->pbCmapID1->_pSysData[i].z;
        gpu->pbCmapID3->_pSysData[i].x                     -= gpu->pbCmapID1->_pSysData[i].w;
        gpu->pbCmapID3->_pSysData[i].y                     -= gpu->pbCmapID2->_pSysData[i].x;
    }
    for (int i = 0; i < gpu->sim.shakeConstraints; i++)
    {
        gpu->pbImageShakeID->_pSysData[i]                   = gpu->pbShakeID->_pSysData[i];
    }
    for (int i = 0; i < gpu->sim.fastShakeConstraints; i++)
    {
        gpu->pbImageFastShakeID->_pSysData[i]               = gpu->pbFastShakeID->_pSysData[i];
    }
    for (int i = 0; i < gpu->sim.slowShakeConstraints; i++)
    {
        gpu->pbImageSlowShakeID1->_pSysData[i]              = gpu->pbSlowShakeID1->_pSysData[i];
        gpu->pbImageSlowShakeID2->_pSysData[i]              = gpu->pbSlowShakeID2->_pSysData[i];
    }
    for (int i = 0; i < gpu->sim.solventMolecules; i++)
    {
        gpu->pbImageSolventAtomID->_pSysData[i]             = gpu->pbSolventAtomID->_pSysData[i];
    }
    for (int i = 0; i < gpu->sim.soluteAtoms; i++)
    {
        gpu->pbImageSoluteAtomID->_pSysData[i]              = gpu->pbSoluteAtomID->_pSysData[i];
    }
    for (int i = 0; i < gpu->sim.EP11s; i++)
    {
        gpu->pbImageExtraPoint11Frame->_pSysData[i]         = gpu->pbExtraPoint11Frame->_pSysData[i];
        gpu->pbImageExtraPoint11Index->_pSysData[i]         = gpu->pbExtraPoint11Index->_pSysData[i];
    }        
    for (int i = 0; i < gpu->sim.EP12s; i++)
    {
        gpu->pbImageExtraPoint12Frame->_pSysData[i]         = gpu->pbExtraPoint12Frame->_pSysData[i];
        gpu->pbImageExtraPoint12Index->_pSysData[i]         = gpu->pbExtraPoint12Index->_pSysData[i];
    }     
    for (int i = 0; i < gpu->sim.EP21s; i++)
    {
        gpu->pbImageExtraPoint21Frame->_pSysData[i]         = gpu->pbExtraPoint21Frame->_pSysData[i];
        gpu->pbImageExtraPoint21Index->_pSysData[i]         = gpu->pbExtraPoint21Index->_pSysData[i];
    }             
    for (int i = 0; i < gpu->sim.EP22s; i++)
    {
        gpu->pbImageExtraPoint22Frame->_pSysData[i]         = gpu->pbExtraPoint22Frame->_pSysData[i];
        gpu->pbImageExtraPoint22Index->_pSysData[i]         = gpu->pbExtraPoint22Index->_pSysData[i];
    }     
        
    if (gpu->pbBondID)
    {
        gpu->pbBondID->Upload();
        gpu->pbImageBondID->Upload(); 
    }
    if (gpu->pbBondAngleID1)
    {
        gpu->pbBondAngleID1->Upload();                                    
        gpu->pbImageBondAngleID1->Upload();
        gpu->pbBondAngleID2->Upload();
        gpu->pbImageBondAngleID2->Upload();
    }
    if (gpu->pbDihedralID1)
    {
        gpu->pbImageDihedralID1->Upload();
        gpu->pbDihedralID2->Upload();
        gpu->pbImageDihedralID2->Upload();
    }
    if (gpu->pbNb14ID)
    {
        gpu->pbNb14ID->Upload();    
        gpu->pbImageNb14ID->Upload(); 
    }
    if (gpu->pbConstraintID)
    {
        gpu->pbConstraintID->Upload();       
        gpu->pbImageConstraintID->Upload();
    }
    if (gpu->pbUBAngle)
    {
        gpu->pbUBAngleID->Upload();
        gpu->pbImageUBAngleID->Upload();
    }
    if (gpu->pbImpDihedralID1)
    {
        gpu->pbImageImpDihedralID1->Upload();
        gpu->pbImpDihedralID2->Upload();
        gpu->pbImageImpDihedralID2->Upload();
    }
    if (gpu->pbCmapID1)
    {
        gpu->pbImageCmapID1->Upload();
        gpu->pbCmapID2->Upload();
        gpu->pbImageCmapID2->Upload();
        gpu->pbCmapID3->Upload();
        gpu->pbImageCmapID3->Upload();
    }
    gpu->pbImageShakeID->Upload();
    gpu->pbImageFastShakeID->Upload();
    gpu->pbImageSlowShakeID1->Upload();
    gpu->pbImageSlowShakeID2->Upload();
    if ((gpu->sim.solventMolecules > 0) || (gpu->sim.soluteAtoms > 0))
    { 
        gpu->pbImageSolventAtomID->Upload();
        gpu->pbImageSoluteAtomID->Upload();
    }
    if (gpu->sim.EP11s > 0)
    {
        gpu->pbImageExtraPoint11Frame->Upload();
        gpu->pbImageExtraPoint11Index->Upload();
    }
    if (gpu->sim.EP12s > 0)
    {
        gpu->pbImageExtraPoint12Frame->Upload();
        gpu->pbImageExtraPoint12Index->Upload();
    }
    if (gpu->sim.EP21s > 0)
    {
        gpu->pbImageExtraPoint21Frame->Upload();
        gpu->pbImageExtraPoint21Index->Upload();
    }
    if (gpu->sim.EP22s > 0)
    {
        gpu->pbImageExtraPoint22Frame->Upload();
        gpu->pbImageExtraPoint22Index->Upload();
    }
    
#if 0    
    printf("%d %d %d %d %d %d %d\n", gpu->sim.bonds, gpu->sim.bondAngles, gpu->sim.dihedrals, gpu->sim.nb14s,
        gpu->sim.constraints, gpu->sim.shakeConstraints, gpu->sim.fastShakeConstraints);
    
    exit(-1);
#endif

     
    // Set up pointers
    gpu->sim.pImageIndex                                    = gpu->pbImageIndex->_pDevData;
    gpu->sim.pImageAtom                                     = gpu->pbImageIndex->_pDevData + gpu->sim.stride;
    gpu->sim.pImageAtomLookup                               = gpu->pbImageIndex->_pDevData + gpu->sim.stride2;
    gpu->sim.pImageHash                                     = gpu->pbImageIndex->_pDevData + gpu->sim.stride3;
    gpu->sim.pImageIndex2                                   = gpu->pbImageIndex->_pDevData + gpu->sim.stride4;
    gpu->sim.pImageAtom2                                    = gpu->pbImageIndex->_pDevData + gpu->sim.stride * 5;
    gpu->sim.pImageHash2                                    = gpu->pbImageIndex->_pDevData + gpu->sim.stride * 6;
    gpu->sim.pAtomXYSaveSP                                  = gpu->pbAtomXYSaveSP->_pDevData;
    gpu->sim.pAtomZSaveSP                                   = gpu->pbAtomZSaveSP->_pDevData;
    gpu->sim.pImageX                                        = gpu->pbImage->_pDevData;
    gpu->sim.pImageY                                        = gpu->pbImage->_pDevData + gpu->sim.stride;
    gpu->sim.pImageZ                                        = gpu->pbImage->_pDevData + gpu->sim.stride2;
    gpu->sim.pImageX2                                       = gpu->pbImage->_pDevData + gpu->sim.stride3;
    gpu->sim.pImageY2                                       = gpu->pbImage->_pDevData + gpu->sim.stride4;
    gpu->sim.pImageZ2                                       = gpu->pbImage->_pDevData + gpu->sim.stride * 5;
    gpu->sim.pImageVelX                                     = gpu->pbImageVel->_pDevData;
    gpu->sim.pImageVelY                                     = gpu->pbImageVel->_pDevData + gpu->sim.stride;
    gpu->sim.pImageVelZ                                     = gpu->pbImageVel->_pDevData + gpu->sim.stride2;
    gpu->sim.pImageVelX2                                    = gpu->pbImageVel->_pDevData + gpu->sim.stride3;
    gpu->sim.pImageVelY2                                    = gpu->pbImageVel->_pDevData + gpu->sim.stride4;
    gpu->sim.pImageVelZ2                                    = gpu->pbImageVel->_pDevData + gpu->sim.stride * 5;
    gpu->sim.pImageLVelX                                    = gpu->pbImageLVel->_pDevData;
    gpu->sim.pImageLVelY                                    = gpu->pbImageLVel->_pDevData + gpu->sim.stride;
    gpu->sim.pImageLVelZ                                    = gpu->pbImageLVel->_pDevData + gpu->sim.stride2;
    gpu->sim.pImageLVelX2                                   = gpu->pbImageLVel->_pDevData + gpu->sim.stride3;
    gpu->sim.pImageLVelY2                                   = gpu->pbImageLVel->_pDevData + gpu->sim.stride4;
    gpu->sim.pImageLVelZ2                                   = gpu->pbImageLVel->_pDevData + gpu->sim.stride * 5;
    gpu->sim.pImageMass                                     = gpu->pbImageMass->_pDevData;
    gpu->sim.pImageInvMass                                  = gpu->pbImageMass->_pDevData + gpu->sim.stride;
    gpu->sim.pImageMass2                                    = gpu->pbImageMass->_pDevData + gpu->sim.stride2;
    gpu->sim.pImageInvMass2                                 = gpu->pbImageMass->_pDevData + gpu->sim.stride3;
    gpu->sim.pImageCharge                                   = gpu->pbImageCharge->_pDevData;
    gpu->sim.pImageCharge2                                  = gpu->pbImageCharge->_pDevData + gpu->sim.stride;
    gpu->sim.pImageSigEps                                   = gpu->pbImageSigEps->_pDevData;
    gpu->sim.pImageSigEps2                                  = gpu->pbImageSigEps->_pDevData + gpu->sim.stride;
    gpu->sim.pImageOutputBuffers                            = gpu->pbImageOutputBuffers->_pDevData;
    gpu->sim.pImageOutputBuffers2                           = gpu->pbImageOutputBuffers->_pDevData + gpu->sim.stride;
    gpu->sim.pImageCellID                                   = gpu->pbImageCellID->_pDevData;
    gpu->sim.pImageCellID2                                  = gpu->pbImageCellID->_pDevData + gpu->sim.stride;

    if ((gpu->sim.solventMolecules > 0) || (gpu->sim.soluteAtoms > 0))
    {
        gpu->sim.pImageSolventAtomID                        = gpu->pbImageSolventAtomID->_pDevData;
        gpu->sim.pImageSoluteAtomID                         = gpu->pbImageSoluteAtomID->_pDevData;
    }

    gpu->sim.pImageBondID                                   = gpu->pbImageBondID->_pDevData;
    gpu->sim.pImageBondAngleID1                             = gpu->pbImageBondAngleID1->_pDevData;
    gpu->sim.pImageBondAngleID2                             = gpu->pbImageBondAngleID2->_pDevData;
    gpu->sim.pImageDihedralID1                              = gpu->pbImageDihedralID1->_pDevData;
    gpu->sim.pImageDihedralID2                              = gpu->pbImageDihedralID2->_pDevData;
    gpu->sim.pImageUBAngleID                                = gpu->pbImageUBAngleID->_pDevData;
    gpu->sim.pImageImpDihedralID1                           = gpu->pbImageImpDihedralID1->_pDevData;
    gpu->sim.pImageImpDihedralID2                           = gpu->pbImageImpDihedralID2->_pDevData;
    gpu->sim.pImageCmapID1                                  = gpu->pbImageCmapID1->_pDevData;
    gpu->sim.pImageCmapID2                                  = gpu->pbImageCmapID2->_pDevData;
    gpu->sim.pImageCmapID3                                  = gpu->pbImageCmapID3->_pDevData;
    gpu->sim.pImageNb14ID                                   = gpu->pbImageNb14ID->_pDevData;
    gpu->sim.pImageConstraintID                             = gpu->pbImageConstraintID->_pDevData;
    gpu->sim.pImageShakeID                                  = gpu->pbImageShakeID->_pDevData;
    gpu->sim.pImageFastShakeID                              = gpu->pbImageFastShakeID->_pDevData;
    gpu->sim.pImageSlowShakeID1                             = gpu->pbImageSlowShakeID1->_pDevData;
    gpu->sim.pImageSlowShakeID2                             = gpu->pbImageSlowShakeID2->_pDevData;
    gpu->sim.pImageExtraPoint11Frame                        = gpu->pbImageExtraPoint11Frame->_pDevData;
    gpu->sim.pImageExtraPoint11Index                        = gpu->pbImageExtraPoint11Index->_pDevData;
    gpu->sim.pImageExtraPoint12Frame                        = gpu->pbImageExtraPoint12Frame->_pDevData;
    gpu->sim.pImageExtraPoint12Index                        = gpu->pbImageExtraPoint12Index->_pDevData;
    gpu->sim.pImageExtraPoint21Frame                        = gpu->pbImageExtraPoint21Frame->_pDevData;
    gpu->sim.pImageExtraPoint21Index                        = gpu->pbImageExtraPoint21Index->_pDevData;
    gpu->sim.pImageExtraPoint22Frame                        = gpu->pbImageExtraPoint22Frame->_pDevData;
    gpu->sim.pImageExtraPoint22Index                        = gpu->pbImageExtraPoint22Index->_pDevData;
    gpu->sim.pNLNonbondCellStartEnd                         = gpu->pbNLNonbondCellStartEnd->_pDevData;
    gpu->sim.pNLChargeGridBufferOffset                      = gpu->pbNLChargeGridBufferOffset->_pDevData;
    gpu->sim.pNLOddBufferOverlapFlag                        = gpu->pbNLOddBufferOverlapFlag->_pDevData;
    gpu->sim.pNLAtomList                                    = gpu->pbNLAtomList->_pDevData;
    gpu->sim.pNLbSkinTestFail                               = gpu->pbNLbSkinTestFail->_pDevData;
    gpu->sim.pNLCellHash                                    = gpu->pbNLCellHash->_pDevData;
    
    // Set up radix sort
    kNLInitRadixSort(gpu);
                     
    gpuCopyConstants();    
    
    gpu_create_outputbuffers_();   
}

extern "C" void gpu_skin_test_()
{
PRINTMETHOD("gpu_skin_test");
     kNLSkinTest(gpu);
}

extern "C" void gpu_build_neighbor_list_()
{
PRINTMETHOD("gpu_build_neighbor_list");
    if (gpu->bNeedNewNeighborList)
    {  
        gpu->bNeedNewNeighborList = false;
        gpu->bNewNeighborList = true;
        kNLGenerateSpatialHash(gpu);
#if 0
        cudaThreadSynchronize();
        gpu->pbImageIndex->Download();
        gpu->pbImage->Download();
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            unsigned int index = gpu->pbImageIndex->_pSysData[i];
            printf("%5d %5d 0x%08x %16.7f %16.7f %16.7f\n", i, index,
                   gpu->pbImageIndex->_pSysData[i + gpu->sim.stride3],
                   gpu->pbImage->_pSysData[i],
                   gpu->pbImage->_pSysData[i + gpu->sim.stride],
                   gpu->pbImage->_pSysData[i + gpu->sim.stride2]
                   );
        }
        
        exit(-1);
#endif      

#if 0 
        // Test hash
        gpu->pbImage->Download();
        gpu->pbImageCellID->Download();
        gpu->pbImageCellID->Download();
        gpu->pbImageIndex->Download();
        unsigned int* pImageIndex = gpu->sim.pImageIndex - gpu->pbImageIndex->_pDevData + gpu->pbImageIndex->_pSysData;
        PMEDouble* pImageX = gpu->pbImage->_pSysData;
        unsigned int* pImageHash = gpu->sim.pImageHash - gpu->pbImageIndex->_pDevData + gpu->pbImageIndex->_pSysData;
        unsigned int* pImageCellID = gpu->pbImageCellID->_pSysData;
        if (gpu->sim.pImageX != gpu->pbImage->_pDevData)
        {
            pImageX = gpu->pbImage->_pSysData + gpu->sim.stride3;
            pImageCellID = gpu->pbImageCellID->_pSysData + gpu->sim.stride;
        }
        PMEDouble* pImageY = pImageX + gpu->sim.stride;
        PMEDouble* pImageZ = pImageX + gpu->sim.stride2;
        for (int i = 0; i < gpu->sim.stride; i++)
        {
            printf("B: %7d %7d 0x%08x %16.7f %16.7f %16.7f 0x%08x\n", i,
                   pImageIndex[i], pImageHash[i], pImageX[i], pImageY[i], pImageZ[i], pImageCellID[i]);
        }    
        printf("\n\n\n");       
#endif
      
        kNLRadixSort(gpu);

#if 0
        cudaThreadSynchronize();
        gpu->pbImage->Download();
        gpu->pbImageCellID->Download();
        gpu->pbImageCellID->Download();
        gpu->pbImageIndex->Download();
        unsigned int* pImageIndex = gpu->sim.pImageIndex - gpu->pbImageIndex->_pDevData + gpu->pbImageIndex->_pSysData;
        PMEDouble* pImageX = gpu->pbImage->_pSysData;
        unsigned int* pImageHash = gpu->sim.pImageHash - gpu->pbImageIndex->_pDevData + gpu->pbImageIndex->_pSysData;
        unsigned int* pImageCellID = gpu->pbImageCellID->_pSysData;
        unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);

        if (gpu->sim.pImageX != gpu->pbImage->_pDevData)
        {
            pImageX = gpu->pbImage->_pSysData + gpu->sim.stride3;
            pImageCellID = gpu->pbImageCellID->_pSysData + gpu->sim.stride;
        }
        PMEDouble* pImageY = pImageX + gpu->sim.stride;
        PMEDouble* pImageZ = pImageX + gpu->sim.stride2;
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            printf("A: %7d %7d %7d 0x%08x %16.7f %16.7f %16.7f 0x%08x\n", i,
                   pImageIndex[i], pImageAtomLookup[i], pImageHash[i], pImageX[i], pImageY[i], pImageZ[i], pImageCellID[i]);
        }
        exit(-1);    
#endif


        kNLRemapImage(gpu);
#if 0        
        gpu->pbImageCellID->Download();
        for (int i = 0; i < gpu->sim.stride; i++)
        {   
            printf("%6d 0x%08x 0x%08x\n", i, gpu->pbImageCellID->_pSysData[i + gpu->sim.stride], gpu->pbImageCellID->_pSysData[i]);
        }  
        exit(-1);
#endif                
        kNLClearCellBoundaries(gpu);
        kNLCalculateCellBoundaries(gpu);
#ifdef MPI        
        // Regenerate force send/receives
        gpu->pbNLNonbondCellStartEnd->Download();  

        
        // Calculate overall limits and PME scatter info
        for (int i = 0; i < gpu->nGpus; i++)
        {
            if (gpu->pMinLocalCell[i] != gpu->pMaxLocalCell[i])  
            {    
                // We have to detect crazy users who have large areas of vacuum in their simulation.
                // Technically speaking, I suspect we ought to crash the simulation with a nasty warning message
                // but in the interest of world peace, here's a complicated loop to handle this (1 or more empty nonbond cells).
                int j                                   = gpu->pMinLocalCell[i];
                while ((j < gpu->pMaxLocalCell[i]) && (gpu->pbNLNonbondCellStartEnd->_pSysData[j].x == gpu->pbNLNonbondCellStartEnd->_pSysData[j].y))
                    j++;
                if (j < gpu->pMaxLocalCell[i])  
                {
                    gpu->pMinLocalAtom[i]                   = gpu->pbNLNonbondCellStartEnd->_pSysData[j].x;
                    int k                                   = gpu->pMaxLocalCell[i] - 1;

                    // We'll never get here unless there's at least one nonempty cell so no need to check cell bounds
                    while (gpu->pbNLNonbondCellStartEnd->_pSysData[k].x == gpu->pbNLNonbondCellStartEnd->_pSysData[k].y)
                        k--;
                    gpu->pMaxLocalAtom[i]                   = gpu->pbNLNonbondCellStartEnd->_pSysData[k].y;
                }
                else
                {
                    gpu->pMinLocalAtom[i]                   = 0;
                    gpu->pMaxLocalAtom[i]                   = 0;
                }
            }
            else
            {
                gpu->pMinLocalAtom[i]                       = 0;
                gpu->pMaxLocalAtom[i]                       = 0;
            }     
            gpu->pPMEStart[i]                               = gpu->pMinLocalAtom[i] * 3;
            gpu->pPMELength[i]                              = (gpu->pMaxLocalAtom[i] - gpu->pMinLocalAtom[i]) * 3;   
            gpu->pAllGathervRecvCountAoS[i]                 = (gpu->pMaxLocalAtom[i] - gpu->pMinLocalAtom[i]) * 3;
            gpu->pAllGathervRecvDisplAoS[i]                 = gpu->pMinLocalAtom[i] * 3; 
        }
  
        // Calculate local limits
        if (gpu->minLocalCell != gpu->maxLocalCell)  
        {    
            gpu->sim.minLocalAtom                           = gpu->pMinLocalAtom[gpu->gpuID];
            gpu->sim.maxLocalAtom                           = gpu->pMaxLocalAtom[gpu->gpuID];
            int j                                           = gpu->minProcessedCell;
            int jlimit                                      = gpu->maxProcessedCell % gpu->sim.cells;
            while ((j != jlimit) && (gpu->pbNLNonbondCellStartEnd->_pSysData[j].x == gpu->pbNLNonbondCellStartEnd->_pSysData[j].y))
            {
                j++;
                if (j == gpu->sim.cells)
                    j                                       = 0;
            }
            if (j != jlimit)  
            {
                gpu->sim.minProcessedAtom                   = gpu->pbNLNonbondCellStartEnd->_pSysData[j].x;
                int k                                       = (gpu->maxProcessedCell - 1) % gpu->sim.cells;

                // Once again we'll never get here unless there's at least one nonempty cell so no need to check cell bounds
                while (gpu->pbNLNonbondCellStartEnd->_pSysData[k].x == gpu->pbNLNonbondCellStartEnd->_pSysData[k].y)
                {
                    k--;
                    if (k == -1)
                        k                                   = gpu->sim.cells - 1;
                }       
                gpu->sim.maxProcessedAtom                   = gpu->pbNLNonbondCellStartEnd->_pSysData[k].y;
                gpu->sim.processedAtoms                     = gpu->sim.maxProcessedAtom - gpu->sim.minProcessedAtom;
                if (gpu->sim.processedAtoms <= 0)
                {
                    gpu->sim.processedAtoms                += gpu->sim.atoms;
                }
            }
            else
            {          
                gpu->sim.minProcessedAtom                   = 0;
                gpu->sim.maxProcessedAtom                   = 0;
                gpu->sim.processedAtoms                     = 0;
            }
            gpu->sim.minReducedAtom                         = gpu->sim.minProcessedAtom & (0xffffffff ^ GRIDBITSMASK);
            gpu->sim.maxReducedAtom                         = (gpu->sim.maxProcessedAtom + (GRID - 1)) & (0xffffffff ^ GRIDBITSMASK);   
            gpu->sim.reducedAtoms                           = gpu->sim.maxReducedAtom - gpu->sim.minReducedAtom;
            if (gpu->sim.reducedAtoms <= 0)
            {
                gpu->sim.reducedAtoms                      += gpu->sim.paddedNumberOfAtoms;
            }                     
        }
        else
        {
            gpu->sim.minLocalAtom                           = 0;
            gpu->sim.maxLocalAtom                           = 0;
            gpu->sim.minProcessedAtom                       = 0;
            gpu->sim.maxProcessedAtom                       = 0;  
            gpu->sim.processedAtoms                         = 0;   
            gpu->sim.minReducedAtom                         = 0;
            gpu->sim.maxReducedAtom                         = 0;  
            gpu->sim.reducedAtoms                           = 0;   
        }
        gpu->sim.localAtoms                                 = gpu->sim.maxLocalAtom - gpu->sim.minLocalAtom;
        gpu->sim.localAtoms3                                = gpu->sim.localAtoms * 3;
        gpu->sim.processedAtoms3                            = gpu->sim.processedAtoms * 3;
        gpu->sim.reducedAtoms3                              = gpu->sim.reducedAtoms * 3;
        // printf("N %d | %7d %7d | %7d %7d | %7d %7d\n", gpu->gpuID, gpu->sim.minLocalAtom, gpu->sim.maxLocalAtom, gpu->minProcessedCell, gpu->maxProcessedCell, gpu->sim.minReducedAtom, gpu->sim.maxReducedAtom);

        // Calculate force send limits
        for (int i = 0; i < gpu->forceSendNodes; i++)
        {                                 
            int node                                        = gpu->pForceSendNode[i];   
            int j                                           = gpu->pForceSendMinCell[i];
            while ((j <= gpu->pForceSendMaxCell[i]) && (gpu->pbNLNonbondCellStartEnd->_pSysData[j].x == gpu->pbNLNonbondCellStartEnd->_pSysData[j].y))
                j++;
            if (j <= gpu->pForceSendMaxCell[i])
            {
                int offset                                  = gpu->pbNLNonbondCellStartEnd->_pSysData[j].x - gpu->sim.minReducedAtom;
                if (offset < 0)
                    offset                                 += gpu->sim.paddedNumberOfAtoms;
                gpu->pOutForceSendStart[i]                  = offset * 3;
                gpu->pForceSendStart[i]                     = gpu->pbNLNonbondCellStartEnd->_pSysData[j].x * 3;
                int k                                       = gpu->pForceSendMaxCell[i];
                while (gpu->pbNLNonbondCellStartEnd->_pSysData[k].x == gpu->pbNLNonbondCellStartEnd->_pSysData[k].y)
                    k--;
                gpu->pForceSendLength[i]                    = gpu->pbNLNonbondCellStartEnd->_pSysData[k].y * 3 - gpu->pForceSendStart[i];            
                // printf("FS %d %d %d %d %d %d %d\n", gpu->gpuID, node, gpu->pOutForceSendStart[i], gpu->pForceSendStart[i], gpu->pForceSendLength[i], gpu->sim.minReducedAtom, gpu->sim.maxReducedAtom);
            }
            else
            {
                gpu->pOutForceSendStart[i]                  = 0;
                gpu->pForceSendStart[i]                     = 0;
                gpu->pForceSendLength[i]                    = 0;
            }
        }     
        
        // Calculate force receive limits
        if (gpu->sim.minLocalAtom != gpu->sim.maxLocalAtom)
        {
            int j                                           = gpu->forceReceiveFirstCell;
            while ((j <= gpu->forceReceiveLastCell) && (gpu->pbNLNonbondCellStartEnd->_pSysData[j].x == gpu->pbNLNonbondCellStartEnd->_pSysData[j].y))
                j++;
            if (j <= gpu->forceReceiveLastCell)
            {
                gpu->forceReceiveFirstAtom                  = gpu->pbNLNonbondCellStartEnd->_pSysData[j].x;
                int k                                       = gpu->forceReceiveLastCell;
                while (gpu->pbNLNonbondCellStartEnd->_pSysData[k].x == gpu->pbNLNonbondCellStartEnd->_pSysData[k].y)
                    k--;
                gpu->forceReceiveLastAtom                   = gpu->pbNLNonbondCellStartEnd->_pSysData[k].y;
            }
            else 
            {
                gpu->forceReceiveFirstAtom                  = 0;
                gpu->forceReceiveLastAtom                   = 0;  
            }
        }
        else
        {
            gpu->forceReceiveFirstAtom                      = 0;
            gpu->forceReceiveLastAtom                       = 0;        
        }
        // printf("FR %d %d %d\n", gpu->gpuID, gpu->forceReceiveFirstAtom, gpu->forceReceiveLastAtom);
        
#if 0        
        printf("Node %d, %d to %d local, %d to %d processed, %d total processed, %d %d\n", gpu->gpuID,
        gpu->sim.minLocalAtom,
        gpu->sim.maxLocalAtom,
        gpu->sim.minProcessedAtom,
        gpu->sim.maxProcessedAtom,
        gpu->sim.processedAtoms,
        gpu->minLocalCell,
        gpu->maxLocalCell);
#endif                    
#endif
        gpuCopyConstants();
        kNLCalculateCellCoordinates(gpu);  
        kNLBuildNeighborList(gpu);  
        kNLRemapLocalInteractions(gpu);

#if 0
        {
        gpu->pbImageIndex->Download();
        unsigned int* pImageIndex = gpu->sim.pImageIndex - gpu->pbImageIndex->_pDevData + gpu->pbImageIndex->_pSysData;
        unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);

        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            printf("%6d: %8d %8d\n", i, pImageIndex[i], pImageAtomLookup[i]);
        }
        }
#endif
#if 0
        gpu->pbNLNonbondCellStartEnd->Download();
        for (int i = 0; i < gpu->sim.cells; i++)
        {
            printf("%6d: %8d %8d %8d\n", i, gpu->pbNLNonbondCellStartEnd->_pSysData[i].x, gpu->pbNLNonbondCellStartEnd->_pSysData[i].y,
                gpu->pbNLNonbondCellStartEnd->_pSysData[i].y - gpu->pbNLNonbondCellStartEnd->_pSysData[i].x);
        }
        exit(-1);
#endif    

#if 0
        gpu->pbImageIndex->Download();
        gpu->pbImage->Download();
        for (int i = 0; i < gpu->sim.stride; i++)
        {
            unsigned int index = gpu->pbImageIndex->_pSysData[i + 5 * gpu->sim.stride];
            printf("%5d %5d %16.7f %16.7f %16.7f\n", i, index,
                   gpu->pbImage->_pSysData[i + gpu->sim.stride3],
                   gpu->pbImage->_pSysData[i + gpu->sim.stride * 4],
                   gpu->pbImage->_pSysData[i + gpu->sim.stride * 5]
                   );
        }
        exit(-1);
#endif

#ifdef MPI
        if (gpu->sim.ntp > 0)
            kClearNBForces(gpu);
            
        // Clear force receive buffer again since dimensions have potentially changed
        if (gpu->sim.minLocalAtom != gpu->sim.maxLocalAtom)   
            memset(gpu->pForceData0 + gpu->forceReceiveFirstAtom * 3, 0, (gpu->forceReceiveLastAtom - gpu->forceReceiveFirstAtom) * 3 * sizeof(PMEDouble));     
        MPI_Win_fence(0, gpu->MPIPMEForceWindow);             
#endif
    }  
}

extern "C" void gpu_molecule_list_setup_(int* molecules, listdata_rec listdata[])
{
PRINTMETHOD("gpu_molecule_list_setup");
   
    // Delete previous molecule list
    delete gpu->pbSoluteAtomID;
    delete gpu->pbSoluteAtomMass;
    delete gpu->pbSolute;
    delete gpu->pbUllSolute;
    delete gpu->pbSolventAtomID;
    delete gpu->pbSolvent;
    
    // Count molecule parameters
    int soluteAtoms                             = 0;
    int soluteMolecules                         = 0;
    int solventMolecules                        = 0;
    for (int i = 0; i < *molecules; i++)
    {
        // Distinguish between solute and solvet
        if (listdata[i].cnt < 5)
        {
            solventMolecules++;
        }
        else
        {
            int offset                          = ((listdata[i].cnt + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
            soluteAtoms                        += offset;
            soluteMolecules++;
        }
    }
    int soluteMoleculeStride                    = ((soluteMolecules + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
    int solventMoleculeStride                   = ((solventMolecules + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
    
    // Allocate solvent/solute data
    gpu->pbSoluteAtomID                                     = new GpuBuffer<int>(2 * soluteAtoms);
    gpu->pbSoluteAtomMass                                   = new GpuBuffer<PMEDouble>(soluteAtoms);
    if (soluteMolecules <= gpu->maxPSSoluteMolecules)
    {
        gpu->pbSolute                                       = new GpuBuffer<PMEDouble>(4 * soluteMoleculeStride);
    }
    else
    {
        gpu->pbSolute                                       = new GpuBuffer<PMEDouble>(7 * soluteMoleculeStride);    
    }
    gpu->pbUllSolute                                        = new GpuBuffer<PMEUllInt>(6 * soluteMoleculeStride);
    gpu->pbSolventAtomID                                    = new GpuBuffer<int4>(solventMoleculeStride);
    gpu->pbSolvent                                          = new GpuBuffer<PMEDouble>(8 * solventMoleculeStride);
    
    // Fill in data
    int* pSoluteAtomID                                      = gpu->pbSoluteAtomID->_pSysData;
    int* pSoluteAtomMoleculeID                              = gpu->pbSoluteAtomID->_pSysData + soluteAtoms;
    PMEDouble* pSoluteAtomMass                              = gpu->pbSoluteAtomMass->_pSysData;
    PMEDouble* pSoluteInvMass                               = gpu->pbSolute->_pSysData + soluteMoleculeStride * 3;
    int4* pSolventAtomID                                    = gpu->pbSolventAtomID->_pSysData;
    PMEDouble* pSolventMass1                                = gpu->pbSolvent->_pSysData;
    PMEDouble* pSolventMass2                                = gpu->pbSolvent->_pSysData + solventMoleculeStride;
    PMEDouble* pSolventMass3                                = gpu->pbSolvent->_pSysData + solventMoleculeStride * 2;
    PMEDouble* pSolventMass4                                = gpu->pbSolvent->_pSysData + solventMoleculeStride * 3;
    PMEDouble* pSolventInvMass                              = gpu->pbSolvent->_pSysData + solventMoleculeStride * 7;
    soluteAtoms                                             = 0;
    soluteMolecules                                         = 0;
    solventMolecules                                        = 0;    
    for (int i = 0; i < *molecules; i++)
    {
        // Distinguish between solute and solvent
        double totalMass                                    = 0.0;
        if (listdata[i].cnt < 5)
        {
            pSolventAtomID[solventMolecules].y              = -1;
            pSolventAtomID[solventMolecules].z              = -1;
            pSolventAtomID[solventMolecules].w              = -1;
            pSolventMass2[solventMolecules]                 = 0.0;
            pSolventMass3[solventMolecules]                 = 0.0;
            pSolventMass4[solventMolecules]                 = 0.0;
            double totalMass                                = 0.0;
            switch (listdata[i].cnt)
            {
                case 4:
                    pSolventAtomID[solventMolecules].w      = listdata[i].offset + 3;
                    pSolventMass4[solventMolecules]         = gpu->pbAtomMass->_pSysData[listdata[i].offset + 3];   
                    totalMass                              += gpu->pbAtomMass->_pSysData[listdata[i].offset + 3]; 
                
                case 3:
                    pSolventAtomID[solventMolecules].z      = listdata[i].offset + 2;
                    pSolventMass3[solventMolecules]         = gpu->pbAtomMass->_pSysData[listdata[i].offset + 2];   
                    totalMass                              += gpu->pbAtomMass->_pSysData[listdata[i].offset + 2]; 
                
                case 2:
                    pSolventAtomID[solventMolecules].y      = listdata[i].offset + 1;
                    pSolventMass2[solventMolecules]         = gpu->pbAtomMass->_pSysData[listdata[i].offset + 1];   
                    totalMass                              += gpu->pbAtomMass->_pSysData[listdata[i].offset + 1]; 
                                
                case 1:
                    pSolventAtomID[solventMolecules].x      = listdata[i].offset;
                    pSolventMass1[solventMolecules]         = gpu->pbAtomMass->_pSysData[listdata[i].offset];   
                    totalMass                              += gpu->pbAtomMass->_pSysData[listdata[i].offset];          
            }
            pSolventInvMass[solventMolecules]               = 1.0 / totalMass;
            solventMolecules++;
        }
        else
        {
            double totalMass                                = 0.0;
            for (int j = 0; j < listdata[i].cnt; j++)
            {
                pSoluteAtomID[soluteAtoms + j]              = listdata[i].offset + j;
                pSoluteAtomMoleculeID[soluteAtoms + j]      = soluteMolecules;
                pSoluteAtomMass[soluteAtoms + j]            = gpu->pbAtomMass->_pSysData[listdata[i].offset + j];
                totalMass                                  += gpu->pbAtomMass->_pSysData[listdata[i].offset + j];
            }        
            int offset                                      = ((listdata[i].cnt + (GRID - 1)) >> GRIDBITS) << GRIDBITS;
            for (int j = listdata[i].cnt; j < offset; j++)
            {
                pSoluteAtomID[soluteAtoms + j]              = -1;
                pSoluteAtomMoleculeID[soluteAtoms + j]      = soluteMolecules;
                pSoluteAtomMass[soluteAtoms + j]            = 0.0;
            }
            soluteAtoms                                    += offset;
            pSoluteInvMass[soluteMolecules]                 = 1.0 / totalMass;    
            soluteMolecules++;
        }
    }

    // Upload data
    gpu->pbSoluteAtomID->Upload();
    gpu->pbSoluteAtomMass->Upload();
    gpu->pbSolute->Upload();
    gpu->pbUllSolute->Upload();
    gpu->pbSolventAtomID->Upload();
    gpu->pbSolvent->Upload();
    
    // Set up constant pointers
    gpu->sim.soluteMolecules                    = soluteMolecules;
    gpu->sim.soluteMoleculeStride               = soluteMoleculeStride;
    gpu->sim.soluteAtoms                        = soluteAtoms;
    gpu->sim.solventMolecules                   = solventMolecules;
    gpu->sim.solventMoleculeStride              = solventMoleculeStride;

    gpu->sim.pSoluteCOMX                        = gpu->pbSolute->_pDevData;
    gpu->sim.pSoluteCOMY                        = gpu->pbSolute->_pDevData + soluteMoleculeStride;
    gpu->sim.pSoluteCOMZ                        = gpu->pbSolute->_pDevData + soluteMoleculeStride * 2;
    gpu->sim.pSoluteInvMass                     = gpu->pbSolute->_pDevData + soluteMoleculeStride * 3;
    gpu->sim.pSoluteDeltaCOMX                   = gpu->pbSolute->_pDevData + soluteMoleculeStride * 4;
    gpu->sim.pSoluteDeltaCOMY                   = gpu->pbSolute->_pDevData + soluteMoleculeStride * 5;
    gpu->sim.pSoluteDeltaCOMZ                   = gpu->pbSolute->_pDevData + soluteMoleculeStride * 6;
    gpu->sim.pSoluteUllCOMX                     = gpu->pbUllSolute->_pDevData;
    gpu->sim.pSoluteUllCOMY                     = gpu->pbUllSolute->_pDevData + soluteMoleculeStride;
    gpu->sim.pSoluteUllCOMZ                     = gpu->pbUllSolute->_pDevData + soluteMoleculeStride * 2;
    gpu->sim.pSoluteUllEKCOMX                   = gpu->pbUllSolute->_pDevData + soluteMoleculeStride * 3;
    gpu->sim.pSoluteUllEKCOMY                   = gpu->pbUllSolute->_pDevData + soluteMoleculeStride * 4;
    gpu->sim.pSoluteUllEKCOMZ                   = gpu->pbUllSolute->_pDevData + soluteMoleculeStride * 5;    
    gpu->sim.pSoluteAtomID                      = gpu->pbSoluteAtomID->_pDevData;
    gpu->sim.pSoluteAtomMoleculeID              = gpu->pbSoluteAtomID->_pDevData + soluteAtoms;
    gpu->sim.pSoluteAtomMass                    = gpu->pbSoluteAtomMass->_pDevData;
    
    gpu->sim.pSolventAtomMass1                  = gpu->pbSolvent->_pDevData;
    gpu->sim.pSolventAtomMass2                  = gpu->pbSolvent->_pDevData + solventMoleculeStride;
    gpu->sim.pSolventAtomMass3                  = gpu->pbSolvent->_pDevData + solventMoleculeStride * 2;
    gpu->sim.pSolventAtomMass4                  = gpu->pbSolvent->_pDevData + solventMoleculeStride * 3;
    gpu->sim.pSolventCOMX                       = gpu->pbSolvent->_pDevData + solventMoleculeStride * 4;
    gpu->sim.pSolventCOMY                       = gpu->pbSolvent->_pDevData + solventMoleculeStride * 5;
    gpu->sim.pSolventCOMZ                       = gpu->pbSolvent->_pDevData + solventMoleculeStride * 6;
    gpu->sim.pSolventInvMass                    = gpu->pbSolvent->_pDevData + solventMoleculeStride * 7;
    gpu->sim.pSolventAtomID                     = gpu->pbSolventAtomID->_pDevData;
    
    gpuCopyConstants();
    
#ifdef GVERBOSE
    printf("Found %d molecules, %d solvent, %d solute, and %d effective solute atoms\n", *molecules, solventMolecules, soluteMolecules, soluteAtoms);
#endif
}

extern "C" void gpu_constraint_molecule_list_setup_(int* molecules, listdata_rec listdata[], double atm_xc[])
{
    PRINTMETHOD("gpu_constraint_molecule_list_setup");
    // Delete previous constained molecule list
    delete gpu->pbConstraintSoluteID;
    delete gpu->pbConstraintSoluteAtom;
    delete gpu->pbConstraintSolute;
    delete gpu->pbConstraintUllSolute;
    delete gpu->pbConstraintSolventAtoms;
    delete gpu->pbConstraintSolventAtom;
    delete gpu->pbConstraintSolventConstraint;

    // Allocate and clear molecule IDs for NTP constaints (assumes regular molecule list is already generated!)
    int* pSoluteMoleculeAtoms                               = new int[gpu->sim.soluteMolecules];   
    int* pSoluteMoleculeID                                  = new int[gpu->sim.atoms];
    int* pSolventMoleculeID                                 = new int[gpu->sim.atoms];
    int* pMoleculeID                                        = new int[gpu->sim.atoms];
    int* pAtomConstraintID                                  = new int[gpu->sim.atoms];
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        pMoleculeID[i]                                      = -1;
        pSoluteMoleculeID[i]                                = -1;
        pSolventMoleculeID[i]                               = -1;
        pAtomConstraintID[i]                                = -1;
    }

    // Count molecule parameters and assign molecule IDs to all atoms
    int soluteAtoms                                         = 0;
    int soluteMolecules                                     = 0;
    int solventMolecules                                    = 0;
    for (int i = 0; i < *molecules; i++)
    {
        // Distinguish between solute and solvet
        if (listdata[i].cnt < 5)
        {
            for (int j = listdata[i].offset; j < listdata[i].offset + listdata[i].cnt; j++)
            {
                pSolventMoleculeID[j]                       = solventMolecules;
                pMoleculeID[j]                              = i;
            }            
            solventMolecules++;
        }
        else
        {
            int offset                                      = ((listdata[i].cnt + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
            soluteAtoms                                    += offset;
            for (int j = listdata[i].offset; j < listdata[i].offset + listdata[i].cnt; j++)
            {
                pSoluteMoleculeID[j]                        = soluteMolecules;    
                pMoleculeID[j]                              = i;
            }
            pSoluteMoleculeAtoms[soluteMolecules]           = listdata[i].cnt;                             
            soluteMolecules++;
        }
    }
    
   
    // Now locate and count molecules with constraints
    bool* pSoluteMoleculeFlag                               = new bool[*molecules];
    bool* pSolventMoleculeFlag                              = new bool[*molecules];
    for (int i = 0; i < *molecules; i++)
    {
        pSoluteMoleculeFlag[i]                              = false;
        pSolventMoleculeFlag[i]                             = false;
    }            
    int constraintSoluteAtoms                               = 0;
    int constraintSoluteMolecules                           = 0;
    int constraintSolventMolecules                          = 0;
    for (int i = 0; i < gpu->sim.constraints; i++)
    {
        int2 atom                                           = gpu->pbConstraintID->_pSysData[i];
        int soluteMoleculeID                                = pSoluteMoleculeID[atom.x];
        int solventMoleculeID                               = pSolventMoleculeID[atom.x];
        int moleculeID                                      = pMoleculeID[atom.x];
        pAtomConstraintID[atom.x]                           = i;
        if ((soluteMoleculeID != -1) && !pSoluteMoleculeFlag[moleculeID])
        {
            pSoluteMoleculeFlag[moleculeID]                 = true;
            int offset                                      = ((pSoluteMoleculeAtoms[soluteMoleculeID] + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
            constraintSoluteAtoms                          += offset;
            constraintSoluteMolecules++;                
        }
        else if ((solventMoleculeID != -1) && !pSolventMoleculeFlag[moleculeID])
        {
            pSolventMoleculeFlag[moleculeID]                = true;              
            constraintSolventMolecules++;  
        }
    }
    int constraintSoluteMoleculeStride                      = ((constraintSoluteMolecules + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;  
    int constraintSolventMoleculeStride                     = ((constraintSolventMolecules + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;         

    // Generate constrained solute molecules
    if (constraintSoluteMolecules > 0)
    {
        gpu->pbConstraintSoluteID                           = new GpuBuffer<int>(2 * constraintSoluteAtoms);
        gpu->pbConstraintSoluteAtom                         = new GpuBuffer<PMEDouble>(4 * constraintSoluteAtoms);
        if (constraintSoluteMolecules <= gpu->maxPSSoluteMolecules) 
        {
            gpu->pbConstraintSolute                             = new GpuBuffer<PMEDouble>(4 * constraintSoluteMoleculeStride);
        }
        else
        {
            gpu->pbConstraintSolute                             = new GpuBuffer<PMEDouble>(7 * constraintSoluteMoleculeStride);
        }
        gpu->pbConstraintUllSolute                          = new GpuBuffer<PMEUllInt>(3 * constraintSoluteMoleculeStride);
        int* pConstraintSoluteMoleculeID                    = gpu->pbConstraintSoluteID->_pSysData;
        int* pConstraintSoluteConstraint                    = gpu->pbConstraintSoluteID->_pSysData + constraintSoluteAtoms;
        PMEDouble* pConstraintSoluteAtomX                   = gpu->pbConstraintSoluteAtom->_pSysData;
        PMEDouble* pConstraintSoluteAtomY                   = gpu->pbConstraintSoluteAtom->_pSysData + constraintSoluteAtoms;
        PMEDouble* pConstraintSoluteAtomZ                   = gpu->pbConstraintSoluteAtom->_pSysData + constraintSoluteAtoms * 2;
        PMEDouble* pConstraintSoluteAtomMass                = gpu->pbConstraintSoluteAtom->_pSysData + constraintSoluteAtoms * 3;
        PMEDouble* pConstraintSoluteInvMass                 = gpu->pbConstraintSolute->_pSysData + constraintSoluteMoleculeStride * 3;      
        constraintSoluteMolecules                           = 0;
        constraintSoluteAtoms                               = 0;
        for (int i = 0; i < *molecules; i++)
        {
            if (pSoluteMoleculeFlag[i])
            {
                double totalMass                                            = 0.0;
                int offset                                                  = ((listdata[i].cnt + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
                for (int j = 0; j < listdata[i].cnt; j++)
                {
                    pConstraintSoluteMoleculeID[constraintSoluteAtoms + j]  = constraintSoluteMolecules;
                    pConstraintSoluteConstraint[constraintSoluteAtoms + j]  = pAtomConstraintID[listdata[i].offset + j];
                    pConstraintSoluteAtomMass[constraintSoluteAtoms + j]    = gpu->pbAtomMass->_pSysData[listdata[i].offset + j];
                    pConstraintSoluteAtomX[constraintSoluteAtoms + j]       = atm_xc[(listdata[i].offset + j) * 3];
                    pConstraintSoluteAtomY[constraintSoluteAtoms + j]       = atm_xc[(listdata[i].offset + j) * 3 + 1];
                    pConstraintSoluteAtomZ[constraintSoluteAtoms + j]       = atm_xc[(listdata[i].offset + j) * 3 + 2];
                    totalMass                                              += gpu->pbAtomMass->_pSysData[listdata[i].offset + j];
                }
                for (int j = listdata[i].cnt; j < offset; j++)
                {
                    pConstraintSoluteMoleculeID[constraintSoluteAtoms + j]  = constraintSoluteMolecules;                
                    pConstraintSoluteConstraint[constraintSoluteAtoms + j]  = -1;
                    pConstraintSoluteAtomX[constraintSoluteAtoms + j]       = (PMEDouble)0.0;
                    pConstraintSoluteAtomY[constraintSoluteAtoms + j]       = (PMEDouble)0.0;
                    pConstraintSoluteAtomZ[constraintSoluteAtoms + j]       = (PMEDouble)0.0;
                    pConstraintSoluteMoleculeID[constraintSoluteAtoms + j]  = constraintSoluteMolecules;
                    pConstraintSoluteAtomMass[constraintSoluteAtoms + j]    = (PMEDouble)0.0;
                }                
                constraintSoluteAtoms                                      += offset;
                pConstraintSoluteInvMass[constraintSoluteMolecules]         = 1.0 / totalMass;
                constraintSoluteMolecules++;
            }   
        }   
        
        // Upload data
        gpu->pbConstraintSoluteID->Upload();
        gpu->pbConstraintSoluteAtom->Upload();
        gpu->pbConstraintSolute->Upload();
        gpu->pbConstraintUllSolute->Upload(); 
        
        // Set up pointers
        gpu->sim.pConstraintSoluteAtomMoleculeID                            = gpu->pbConstraintSoluteID->_pDevData;
        gpu->sim.pConstraintSoluteAtomConstraintID                          = gpu->pbConstraintSoluteID->_pDevData + constraintSoluteAtoms;
        gpu->sim.pConstraintSoluteAtomX                                     = gpu->pbConstraintSoluteAtom->_pDevData;
        gpu->sim.pConstraintSoluteAtomY                                     = gpu->pbConstraintSoluteAtom->_pDevData + constraintSoluteAtoms;
        gpu->sim.pConstraintSoluteAtomZ                                     = gpu->pbConstraintSoluteAtom->_pDevData + constraintSoluteAtoms * 2;
        gpu->sim.pConstraintSoluteAtomMass                                  = gpu->pbConstraintSoluteAtom->_pDevData + constraintSoluteAtoms * 3;    
        gpu->sim.pConstraintSoluteCOMX                                      = gpu->pbConstraintSolute->_pDevData;
        gpu->sim.pConstraintSoluteCOMY                                      = gpu->pbConstraintSolute->_pDevData + constraintSoluteMoleculeStride;
        gpu->sim.pConstraintSoluteCOMZ                                      = gpu->pbConstraintSolute->_pDevData + constraintSoluteMoleculeStride * 2;
        gpu->sim.pConstraintSoluteInvMass                                   = gpu->pbConstraintSolute->_pDevData + constraintSoluteMoleculeStride * 3;
        gpu->sim.pConstraintSoluteDeltaCOMX                                 = gpu->pbConstraintSolute->_pDevData + constraintSoluteMoleculeStride * 4;
        gpu->sim.pConstraintSoluteDeltaCOMY                                 = gpu->pbConstraintSolute->_pDevData + constraintSoluteMoleculeStride * 5;
        gpu->sim.pConstraintSoluteDeltaCOMZ                                 = gpu->pbConstraintSolute->_pDevData + constraintSoluteMoleculeStride * 6;
        gpu->sim.pConstraintSoluteUllCOMX                                   = gpu->pbConstraintUllSolute->_pDevData;
        gpu->sim.pConstraintSoluteUllCOMY                                   = gpu->pbConstraintUllSolute->_pDevData + constraintSoluteMoleculeStride;
        gpu->sim.pConstraintSoluteUllCOMZ                                   = gpu->pbConstraintUllSolute->_pDevData + constraintSoluteMoleculeStride * 2;            
    }

    // Generate constrained solvent molecules
    if (constraintSolventMolecules > 0)
    {   
        gpu->pbConstraintSolventAtoms                       = new GpuBuffer<int>(constraintSolventMoleculeStride);
        gpu->pbConstraintSolventAtom                        = new GpuBuffer<PMEDouble>(20 * constraintSolventMoleculeStride);
        gpu->pbConstraintSolventConstraint                  = new GpuBuffer<int4>(constraintSolventMoleculeStride);
        int* pConstraintSolventAtoms                        = gpu->pbConstraintSolventAtoms->_pSysData;
        PMEDouble* pConstraintSolventAtomMass1              = gpu->pbConstraintSolventAtom->_pSysData;
        PMEDouble* pConstraintSolventAtomX1                 = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride;
        PMEDouble* pConstraintSolventAtomY1                 = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride * 2;
        PMEDouble* pConstraintSolventAtomZ1                 = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride * 3;
        PMEDouble* pConstraintSolventAtomMass2              = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride * 4;
        PMEDouble* pConstraintSolventAtomX2                 = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride * 5;
        PMEDouble* pConstraintSolventAtomY2                 = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride * 6;
        PMEDouble* pConstraintSolventAtomZ2                 = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride * 7;        
        PMEDouble* pConstraintSolventAtomMass3              = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride * 8;
        PMEDouble* pConstraintSolventAtomX3                 = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride * 9;
        PMEDouble* pConstraintSolventAtomY3                 = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride * 10;
        PMEDouble* pConstraintSolventAtomZ3                 = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride * 11;
        PMEDouble* pConstraintSolventAtomMass4              = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride * 12;
        PMEDouble* pConstraintSolventAtomX4                 = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride * 13;
        PMEDouble* pConstraintSolventAtomY4                 = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride * 14;
        PMEDouble* pConstraintSolventAtomZ4                 = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride * 15;
        PMEDouble* pConstraintSolventInvMass                = gpu->pbConstraintSolventAtom->_pSysData + constraintSolventMoleculeStride * 16;
        int4* pConstraintSolventConstraint                  = gpu->pbConstraintSolventConstraint->_pSysData;
        
        constraintSolventMolecules                          = 0;
        for (int i = 0; i < *molecules; i++)
        {
            if (pSolventMoleculeFlag[i])
            {
                double totalMass                                                    = 0.0;
                pConstraintSolventAtoms[constraintSolventMolecules]                 = listdata[i].cnt;
                pConstraintSolventAtomMass2[constraintSolventMolecules]             = 0.0;
                pConstraintSolventAtomX2[constraintSolventMolecules]                = 0.0;
                pConstraintSolventAtomY2[constraintSolventMolecules]                = 0.0;
                pConstraintSolventAtomZ2[constraintSolventMolecules]                = 0.0;
                pConstraintSolventAtomMass3[constraintSolventMolecules]             = 0.0;
                pConstraintSolventAtomX3[constraintSolventMolecules]                = 0.0;
                pConstraintSolventAtomY3[constraintSolventMolecules]                = 0.0;
                pConstraintSolventAtomZ3[constraintSolventMolecules]                = 0.0;
                pConstraintSolventAtomMass4[constraintSolventMolecules]             = 0.0;
                pConstraintSolventAtomX4[constraintSolventMolecules]                = 0.0;
                pConstraintSolventAtomY4[constraintSolventMolecules]                = 0.0;
                pConstraintSolventAtomZ4[constraintSolventMolecules]                = 0.0;
                pConstraintSolventConstraint[constraintSolventMolecules].y          = -1;
                pConstraintSolventConstraint[constraintSolventMolecules].z          = -1;
                pConstraintSolventConstraint[constraintSolventMolecules].w          = -1;
                
                switch (listdata[i].cnt)
                {
                    case 4:
                        pConstraintSolventAtomMass4[constraintSolventMolecules]     = gpu->pbAtomMass->_pSysData[listdata[i].offset + 3]; 
                        pConstraintSolventAtomX4[constraintSolventMolecules]        = atm_xc[(listdata[i].offset + 3) * 3];
                        pConstraintSolventAtomY4[constraintSolventMolecules]        = atm_xc[(listdata[i].offset + 3) * 3 + 1];
                        pConstraintSolventAtomZ4[constraintSolventMolecules]        = atm_xc[(listdata[i].offset + 3) * 3 + 2];
                        pConstraintSolventConstraint[constraintSolventMolecules].w  = pAtomConstraintID[listdata[i].offset + 3];
                        totalMass                                                  += gpu->pbAtomMass->_pSysData[listdata[i].offset + 3]; 
                
                    case 3:
                        pConstraintSolventAtomMass3[constraintSolventMolecules]     = gpu->pbAtomMass->_pSysData[listdata[i].offset + 2]; 
                        pConstraintSolventAtomX3[constraintSolventMolecules]        = atm_xc[(listdata[i].offset + 2) * 3];
                        pConstraintSolventAtomY3[constraintSolventMolecules]        = atm_xc[(listdata[i].offset + 2) * 3 + 1];
                        pConstraintSolventAtomZ3[constraintSolventMolecules]        = atm_xc[(listdata[i].offset + 2) * 3 + 2];
                        pConstraintSolventConstraint[constraintSolventMolecules].z  = pAtomConstraintID[listdata[i].offset + 2];
                        totalMass                                                  += gpu->pbAtomMass->_pSysData[listdata[i].offset + 2];  
                
                    case 2:
                        pConstraintSolventAtomMass2[constraintSolventMolecules]     = gpu->pbAtomMass->_pSysData[listdata[i].offset + 1]; 
                        pConstraintSolventAtomX2[constraintSolventMolecules]        = atm_xc[(listdata[i].offset + 1) * 3];
                        pConstraintSolventAtomY2[constraintSolventMolecules]        = atm_xc[(listdata[i].offset + 1) * 3 + 1];
                        pConstraintSolventAtomZ2[constraintSolventMolecules]        = atm_xc[(listdata[i].offset + 1) * 3 + 2];
                        pConstraintSolventConstraint[constraintSolventMolecules].y  = pAtomConstraintID[listdata[i].offset + 1];
                        totalMass                                                  += gpu->pbAtomMass->_pSysData[listdata[i].offset + 1];  
                                
                    case 1:
                        pConstraintSolventAtomMass1[constraintSolventMolecules]     = gpu->pbAtomMass->_pSysData[listdata[i].offset]; 
                        pConstraintSolventAtomX1[constraintSolventMolecules]        = atm_xc[listdata[i].offset * 3];
                        pConstraintSolventAtomY1[constraintSolventMolecules]        = atm_xc[listdata[i].offset * 3 + 1];
                        pConstraintSolventAtomZ1[constraintSolventMolecules]        = atm_xc[listdata[i].offset * 3 + 2];
                        pConstraintSolventConstraint[constraintSolventMolecules].x  = pAtomConstraintID[listdata[i].offset];
                        totalMass                                                  += gpu->pbAtomMass->_pSysData[listdata[i].offset];      
                } 
                pConstraintSolventInvMass[constraintSolventMolecules]               = 1.0 / totalMass;
                constraintSolventMolecules++;
            }
        }
        
        // Upload Data
        gpu->pbConstraintSolventAtoms->Upload();
        gpu->pbConstraintSolventAtom->Upload();
        gpu->pbConstraintSolventConstraint->Upload();
        
        // Set up pointers
        gpu->sim.pConstraintSolventAtoms                                            = gpu->pbConstraintSolventAtoms->_pDevData;
        gpu->sim.pConstraintSolventConstraint                                       = gpu->pbConstraintSolventConstraint->_pDevData;
        gpu->sim.pConstraintSolventAtomMass1                                        = gpu->pbConstraintSolventAtom->_pDevData;
        gpu->sim.pConstraintSolventAtomX1                                           = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride;
        gpu->sim.pConstraintSolventAtomY1                                           = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 2;
        gpu->sim.pConstraintSolventAtomZ1                                           = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 3;
        gpu->sim.pConstraintSolventAtomMass2                                        = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 4;
        gpu->sim.pConstraintSolventAtomX2                                           = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 5;
        gpu->sim.pConstraintSolventAtomY2                                           = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 6;
        gpu->sim.pConstraintSolventAtomZ2                                           = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 7;
        gpu->sim.pConstraintSolventAtomMass3                                        = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 8;
        gpu->sim.pConstraintSolventAtomX3                                           = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 9;
        gpu->sim.pConstraintSolventAtomY3                                           = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 10;
        gpu->sim.pConstraintSolventAtomZ3                                           = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 11;
        gpu->sim.pConstraintSolventAtomMass4                                        = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 12;
        gpu->sim.pConstraintSolventAtomX4                                           = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 13;
        gpu->sim.pConstraintSolventAtomY4                                           = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 14;
        gpu->sim.pConstraintSolventAtomZ4                                           = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 15;
        gpu->sim.pConstraintSolventInvMass                                          = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 16;
        gpu->sim.pConstraintSolventCOMX                                             = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 17;
        gpu->sim.pConstraintSolventCOMY                                             = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 18;    
        gpu->sim.pConstraintSolventCOMZ                                             = gpu->pbConstraintSolventAtom->_pDevData + constraintSolventMoleculeStride * 19;
    }
    
    // Set up constraint counters
    gpu->sim.constraintSoluteMolecules                                  = constraintSoluteMolecules;
    gpu->sim.constraintSoluteMoleculeStride                             = constraintSoluteMoleculeStride;
    gpu->sim.constraintSoluteAtoms                                      = constraintSoluteAtoms;
    gpu->sim.constraintSolventMolecules                                 = constraintSolventMolecules;

    gpuCopyConstants();

    // Free temporary arrays    
    delete[] pSoluteMoleculeFlag;
    delete[] pSolventMoleculeFlag;
    delete[] pMoleculeID;
    delete[] pSoluteMoleculeID;
    delete[] pSolventMoleculeID;
    delete[] pAtomConstraintID;
    delete[] pSoluteMoleculeAtoms;
    
#ifdef GVERBOSE
    printf("Found %d constrained solvent molecules, %d solute molecules, and %d effective solute atoms\n", constraintSolventMolecules, constraintSoluteMolecules, constraintSoluteAtoms);
#endif    
}

extern "C" void gpu_ntp_setup_()
{
PRINTMETHOD("gpu_pme_ntp_setup");
    gpu->pbNTPData                              = new GpuBuffer<NTPData>(1);
    NTPData* pNTPData                           = gpu->pbNTPData->_pSysData;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
               pNTPData->ucell[i * 3 + j]       = gpu->sim.ucell[i][j];
               pNTPData->ucellf[i * 3 + j]      = gpu->sim.ucell[i][j];
               pNTPData->recip[i * 3 + j]       = gpu->sim.recip[i][j];
               pNTPData->recipf[i * 3 + j]      = gpu->sim.recipf[i][j];
               pNTPData->last_recip[i * 3 + j]  = gpu->sim.recip[i][j];
        }
    }
    pNTPData->one_half_nonbond_skin_squared     = gpu->sim.one_half_nonbond_skin_squared;
    pNTPData->cutPlusSkin2                      = (gpu->sim.cut + gpu->sim.nonbond_skin) * (gpu->sim.cut + gpu->sim.nonbond_skin);
    gpu->sim.pNTPData                           = gpu->pbNTPData->_pDevData;
    gpu->pbNTPData->Upload();
    
    gpuCopyConstants();
    
   
    // Calculate initial COM
    kCalculateCOM(gpu);
    kReduceSoluteCOM(gpu);
    if (gpu->sim.constraints > 0)
    {
        kCalculateConstraintsCOM(gpu);
        kReduceSoluteConstraintsCOM(gpu);  
    }
#if 0
    gpu->pbSolute->Download();
    gpu->pbSolvent->Download();
    for (int i = 0; i < gpu->sim.soluteMolecules; i++)
    {
        printf("%6d %20.10f %20.10f %20.10f\n", i, gpu->pbSolute->_pSysData[i], gpu->pbSolute->_pSysData[i + gpu->sim.soluteMoleculeStride], gpu->pbSolute->_pSysData[i + 2 * gpu->sim.soluteMoleculeStride]);
    }
    int solventMoleculeStride                   = ((gpu->sim.solventMolecules + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
    printf("%d\n", solventMoleculeStride);
    for (int i = 0; i < gpu->sim.solventMolecules; i++)
    {
        printf("%6d %20.10f %20.10f %20.10f\n", i, gpu->pbSolvent->_pSysData[i + 4 * solventMoleculeStride], gpu->pbSolvent->_pSysData[i + 5 * solventMoleculeStride], gpu->pbSolvent->_pSysData[i + 6 * solventMoleculeStride]);
    }
    exit(-1);
#endif
}

extern "C" void gpu_pme_alltasks_setup_(int* nfft1, int* nfft2, int* nfft3, double* prefac1, double* prefac2, double* prefac3, double* ew_coeff, int* ips)
{
PRINTMETHOD("gpu_pme_alltasks_setup");

    // Delete any existing PME state
    delete gpu->pbPrefac;
    delete gpu->pbIFract;
    delete gpu->pbTheta;
    
    // Determine whether running IPS or PME
    if (*ips != 0)
        gpu->sim.bIPSActive                     = true;
    else
        gpu->sim.bIPSActive                     = false;
    
    // Allocate GPU data
    int n1                                      = ((*nfft1 + 1) + PADDING) & PADDINGMASK;
    int n2                                      = ((*nfft2 + 1) + PADDING) & PADDINGMASK;
    int n3                                      = ((*nfft3 + 1) + PADDING) & PADDINGMASK;
    gpu->sim.n2Offset                           = n1;
    gpu->sim.n3Offset                           = n1 + n2;
    gpu->sim.nSum                               = n1 + n2 + n3;
    gpu->sim.nfft1                              = *nfft1;
    gpu->sim.nfft2                              = *nfft2;
    gpu->sim.nfft3                              = *nfft3;
    gpu->sim.nfft1xnfft2                        = gpu->sim.nfft1 * gpu->sim.nfft2;
    gpu->sim.fft_x_dim                          = *nfft1 / 2 + 1;
    gpu->sim.fft_y_dim                          = *nfft2;
    gpu->sim.fft_z_dim                          = *nfft3;
    gpu->sim.nf1                                = *nfft1 / 2;
    if (2 * gpu->sim.nf1 < *nfft1) 
        gpu->sim.nf1++;
    gpu->sim.nf2                                = *nfft2 / 2;
    if (2 * gpu->sim.nf2 < *nfft2) 
        gpu->sim.nf2++;
    gpu->sim.nf3                                = *nfft3 / 2;
    if (2 * gpu->sim.nf3 < *nfft3) 
        gpu->sim.nf3++;    
        
    // Set up PME charge buffer
    gpu->sim.ew_coeff                           = *ew_coeff;
    gpu->sim.ew_coeff2                          = (*ew_coeff) * (*ew_coeff);
    gpu->sim.negTwoEw_coeffRsqrtPI              = -2.0 * gpu->sim.ew_coeff / sqrt(PI);
    gpu->sim.fac                                = (PI * PI) / ((*ew_coeff) * (*ew_coeff));
    gpu->sim.fac2                               = 2.0 * gpu->sim.fac;
    gpu->sim.XYZStride                          = ((2 * gpu->sim.fft_x_dim * gpu->sim.fft_y_dim * gpu->sim.fft_z_dim) + 31) & 0xffffffe0;
    gpu->pbPrefac                               = new GpuBuffer<PMEFloat>(n1 + n2 + n3);
    gpu->pbIFract                               = new GpuBuffer<int>(3 * gpu->sim.stride);
    gpu->pbTheta                                = new GpuBuffer<PMEFloat4>(6 * gpu->sim.stride);
      
    // Copy Prefac data
    for (int i = 0; i < *nfft1; i++)
        gpu->pbPrefac->_pSysData[i + 1]             = prefac1[i];
    for (int i = 0; i < *nfft2; i++)
        gpu->pbPrefac->_pSysData[i + 1 + n1]        = prefac2[i];        
    for (int i = 0; i < *nfft3; i++)
        gpu->pbPrefac->_pSysData[i + 1 + n1 + n2]   = prefac3[i];   
    gpu->pbPrefac->Upload();    
    
    // Set up pointers
    gpu->sim.pPrefac1                           = gpu->pbPrefac->_pDevData;
    gpu->sim.pPrefac2                           = gpu->pbPrefac->_pDevData + n1;
    gpu->sim.pPrefac3                           = gpu->pbPrefac->_pDevData + n1 + n2;
    gpu->sim.pIFractX                           = gpu->pbIFract->_pDevData;
    gpu->sim.pIFractY                           = gpu->pbIFract->_pDevData + gpu->sim.stride;
    gpu->sim.pIFractZ                           = gpu->pbIFract->_pDevData + gpu->sim.stride2;
    gpu->sim.pThetaX                            = gpu->pbTheta->_pDevData;
    gpu->sim.pThetaY                            = gpu->pbTheta->_pDevData + gpu->sim.stride;
    gpu->sim.pThetaZ                            = gpu->pbTheta->_pDevData + gpu->sim.stride2;
    gpu->sim.pDThetaX                           = gpu->pbTheta->_pDevData + gpu->sim.stride3;
    gpu->sim.pDThetaY                           = gpu->pbTheta->_pDevData + gpu->sim.stride4;
    gpu->sim.pDThetaZ                           = gpu->pbTheta->_pDevData + gpu->sim.stride * 5;

    // Set up PME kernel cache preferences
    PMEInitKernels(gpu);    

    // Set up FFT plans
#ifdef use_DPDP
    cufftPlan3d(&(gpu->forwardPlan), gpu->sim.nfft3, gpu->sim.nfft2, gpu->sim.nfft1, CUFFT_D2Z);
    cufftPlan3d(&(gpu->backwardPlan),  gpu->sim.nfft3, gpu->sim.nfft2, gpu->sim.nfft1, CUFFT_Z2D);
#else
    cufftPlan3d(&(gpu->forwardPlan), gpu->sim.nfft3, gpu->sim.nfft2, gpu->sim.nfft1, CUFFT_R2C);
    cufftPlan3d(&(gpu->backwardPlan),  gpu->sim.nfft3, gpu->sim.nfft2, gpu->sim.nfft1, CUFFT_C2R);
#endif
    gpuCopyConstants();
    
    return;
}

extern "C" void gpu_get_self(double* ee_plasma, double* ene)
{
    if ((gpu->ee_plasma != *ee_plasma) || (*ene != gpu->self_energy))
    {
        gpu->ee_plasma                          = *ee_plasma;
        gpu->self_energy                        = *ene;
        gpuCopyConstants();
    }
}

//  Electrostatic IPS parameters:
static const double aipse0                      = -35.0 / 16.0;
static const double aipse1                      =  35.0 / 16.0;
static const double aipse2                      = -21.0 / 16.0;
static const double aipse3                      =   5.0 / 16.0;
static const double pipsec                      =   1.0 + aipse0 + aipse1 + aipse2 + aipse3;
static const double pipse0                      =   aipse0 - pipsec; 
static const double bipse1                      =   2.0 * aipse1;
static const double bipse2                      =   4.0 * aipse2;
static const double bipse3                      =   6.0 * aipse3;
    
//  Dispersion IPS parameters:
static const double aipsvc0                     =   7.0 / 16.0;
static const double aipsvc1                     =   9.0 / 14.0;
static const double aipsvc2                     =  -3.0 / 28.0;
static const double aipsvc3                     =   6.0 / 7.0;
static const double pipsvcc                     =   1.0 + aipsvc0 + aipsvc1 + aipsvc2 + aipsvc3;
static const double pipsvc0                     =   aipsvc0 - pipsvcc;
static const double bipsvc1                     =   2.0 * aipsvc1;
static const double bipsvc2                     =   4.0 * aipsvc2;  
static const double bipsvc3                     =   6.0 * aipsvc3; 

//  Repulsion IPS parameters:
static const double aipsva0                     =  5.0 / 787.0;
static const double aipsva1                     =  9.0 /  26.0;
static const double aipsva2                     = -3.0 /  13.0;
static const double aipsva3                     = 27.0 /  26.0;
static const double pipsvac                     =  1.0 + aipsva0 + aipsva1 + aipsva2 + aipsva3;
static const double pipsva0                     =  aipsva0 - pipsvac;
static const double bipsva1                     =  4.0 * aipsva1;
static const double bipsva2                     =  8.0 * aipsva2;
static const double bipsva3                     = 12.0 * aipsva3;

extern "C" void gpu_ips_setup_(double* rips, double atm_qterm[], int* ntypes, int iac[], int ico[], double cn1[], double cn2[], double* eipssnb, double* eipssel, double* virips)
{
PRINTMETHOD("gpu_ips_setup");
      
    gpu->sim.rips                               = *rips;
    
    // Copied data
    gpu->sim.eipssnb                            = *eipssnb;
    gpu->sim.eipssel                            = *eipssel;
    gpu->sim.virips                             = 0.5 * *virips;
    
    // Derived coefficients
    gpu->sim.rips2                              = *rips * *rips;
    gpu->sim.ripsr                              = 1.0 / (*rips);
    gpu->sim.rips2r                             = 1.0 / (*rips * *rips);
    gpu->sim.rips6r                             = pow(1.0 / *rips, 6);
    gpu->sim.rips12r                            = pow(1.0 / *rips, 12); 
    
    // Premultiplied coefficients
    gpu->sim.aipse0                             = aipse0 / *rips;
    gpu->sim.aipse1                             = aipse1 / pow(*rips, 3);
    gpu->sim.aipse2                             = aipse2 / pow(*rips, 5);
    gpu->sim.aipse3                             = aipse3 / pow(*rips, 7);
    gpu->sim.bipse1                             = bipse1 / pow(*rips, 3);
    gpu->sim.bipse2                             = bipse2 / pow(*rips, 5);
    gpu->sim.bipse3                             = bipse3 / pow(*rips, 7);
    gpu->sim.pipsec                             = pipsec / *rips;
    
    gpu->sim.aipsvc0                            = aipsvc0 / pow(*rips, 6);
    gpu->sim.aipsvc1                            = aipsvc1 / pow(*rips, 8);
    gpu->sim.aipsvc2                            = aipsvc2 / pow(*rips, 10);
    gpu->sim.aipsvc3                            = aipsvc3 / pow(*rips, 12);
    gpu->sim.bipsvc1                            = bipsvc1 / pow(*rips, 8);
    gpu->sim.bipsvc2                            = bipsvc2 / pow(*rips, 10);
    gpu->sim.bipsvc3                            = bipsvc3 / pow(*rips, 12);
    gpu->sim.pipsvcc                            = pipsvcc / pow(*rips, 6);
    
    gpu->sim.aipsva0                            = aipsva0 / pow(*rips, 12);
    gpu->sim.aipsva1                            = aipsva1 / pow(*rips, 16);
    gpu->sim.aipsva2                            = aipsva2 / pow(*rips, 20);
    gpu->sim.aipsva3                            = aipsva3 / pow(*rips, 24);
    gpu->sim.bipsva1                            = bipsva1 / pow(*rips, 16);
    gpu->sim.bipsva2                            = bipsva2 / pow(*rips, 20);
    gpu->sim.bipsva3                            = bipsva3 / pow(*rips, 24);
    gpu->sim.pipsvac                            = pipsvac / pow(*rips, 12);
    
    
    // Calculate system energies
    double EEL                                  = 0.0;
    double ENB                                  = 0.0;
    double ripsr                                = 1.0 / *rips;
    double rips6r                               = pow(ripsr, 6);
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        double qi                               = atm_qterm[i];
        double qiqj                             = qi * qi; 
        EEL                                    += 0.5 * qiqj * pipse0 * ripsr;        
        int j                                   = iac[i] - 1;
        int nbtype                              = ico[*ntypes * j + j] - 1;
        double A                                = cn1[nbtype];
        double B                                = cn2[nbtype];          
        ENB                                    += 0.5 * (A * pipsva0 * rips6r - B * pipsvc0) * rips6r;
    }
    gpu->sim.EIPSNB                             = ENB;
    gpu->sim.EIPSEL                             = EEL;
    gpuCopyConstants();
}

extern "C" void gpu_ips_update_(double* eipssnb, double* eipssel, double* virips) 
{
PRINTMETHOD("gpu_ips_update");
    gpu->sim.eipssnb                            = *eipssnb;
    gpu->sim.eipssel                            = *eipssel;
    gpu->sim.virips                             = 0.5 * *virips;
}


extern "C" void gpu_get_grid_weights_()
{
PRINTMETHOD("gpu_get_grid_weights");

    kPMEGetGridWeights(gpu);

#if 0
    cudaThreadSynchronize();
    gpu->pbIFract->Download();  
    gpu->pbImageIndex->Download();
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        printf("%5d: %9d %9d %9d\n", i, gpu->pbIFract->_pSysData[i], gpu->pbIFract->_pSysData[i + gpu->sim.stride], gpu->pbIFract->_pSysData[i + gpu->sim.stride2]);
    }
    exit(-1);
#endif         
    
   
#if 0
    gpu->pbTheta->Download();
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        printf("%5d: %13.7f %13.7f %13.7f %13.7f\n", i,
            gpu->pbTheta->_pSysData[i].x,
            gpu->pbTheta->_pSysData[i].y,
            gpu->pbTheta->_pSysData[i].z,
            gpu->pbTheta->_pSysData[i].w
              );
        printf("%5d: %13.7f %13.7f %13.7f %13.7f\n", i,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride].x,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride].y,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride].z,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride].w
              );
        printf("%5d: %13.7f %13.7f %13.7f %13.7f\n", i,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 2].x,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 2].y,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 2].z,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 2].w
              );
            
    }
    exit(-1);
#endif
#if 0
    gpu->pbTheta->Download();
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        printf("%5d: %13.7f %13.7f %13.7f %13.7f\n", i,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 3].x,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 3].y,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 3].z,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 3].w
              );
        printf("%5d: %13.7f %13.7f %13.7f %13.7f\n", i,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 4].x,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 4].y,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 4].z,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 4].w
              );
        printf("%5d: %13.7f %13.7f %13.7f %13.7f\n", i,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 5].x,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 5].y,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 5].z,
            gpu->pbTheta->_pSysData[i + gpu->sim.stride * 5].w
              );
            
    }
    exit(-1);
#endif

#if 0   
    gpu->pbIFract->Download();
    for (int i = 0; i < gpu->sim.atoms; i++)
        printf("%5d: %9d %9d %9d\n", i, gpu->pbIFract->_pSysData[i], gpu->pbIFract->_pSysData[i + gpu->sim.stride],
        gpu->pbIFract->_pSysData[i + gpu->sim.stride2]);
    exit(-1);
#endif
    
}

extern "C" void gpu_fill_charge_grid_()
{
PRINTMETHOD("gpu_fill_charge_grid");
    gpu_build_neighbor_list_();
    kPMEGetGridWeights(gpu);
    
#if 0
    cudaThreadSynchronize();
    gpu->pbIFract->Download();  
    gpu->pbImageIndex->Download();
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        for (int j = 0; j < gpu->sim.atoms; j++)
            if (gpu->pbImageIndex->_pSysData[j + gpu->sim.stride2] == i)
                printf("%5d: %9d %9d %9d\n", i, gpu->pbIFract->_pSysData[j], gpu->pbIFract->_pSysData[j + gpu->sim.stride], gpu->pbIFract->_pSysData[j + gpu->sim.stride2]);
    }
    exit(-1);
#endif        
    kPMEClearChargeGridBuffer(gpu);  
    kPMEFillChargeGridBuffer(gpu);
    kPMEReduceChargeGridBuffer(gpu);  
    
#if 0
    gpu->pbXYZ_q->Download();
    for (int i = 0; i < gpu->sim.nfft1; i++)
    {
        for (int j = 0; j < gpu->sim.nfft2; j++)
        {
            for (int k = 0; k < gpu->sim.nfft3; k++)
            {
                printf("%3d %3d %3d %32.15f\n", i, j, k, gpu->pbXYZ_q->_pSysData[(k * gpu->sim.nfft2 + j) * gpu->sim.nfft1 + i]); 
            }
        }
    }
    exit(-1);
#endif
}

extern "C" void gpu_fft3drc_forward_()
{

    PRINTMETHOD("gpu_fft3drc_forward");
#ifdef use_DPDP
    cufftExecD2Z(gpu->forwardPlan, gpu->sim.pXYZ_q, gpu->sim.pXYZ_qt);
#else
    cufftExecR2C(gpu->forwardPlan, gpu->sim.pXYZ_q, gpu->sim.pXYZ_qt);
#endif
#if 0 
    gpu->pbXYZ_qt->Download();
    PMEFloat2* pFloat2 = gpu->pbXYZ_qt->_pSysData;       
    for (int i = 0; i < gpu->sim.fft_x_dim; i++)
    {
        for (int j = 0; j < gpu->sim.fft_y_dim; j++)
        {
            for (int k = 0; k < gpu->sim.fft_z_dim; k++)
            {
                printf("%3d %3d %3d: %20.13f %20.13f\n", i, j, k, 
                    pFloat2[(k * gpu->sim.fft_y_dim + j) * gpu->sim.fft_x_dim + i].x,
                    pFloat2[(k * gpu->sim.fft_y_dim + j) * gpu->sim.fft_x_dim + i].y);
            }
        }
    }
    exit(-1);
#endif
}

extern "C" void gpu_scalar_sumrc_(double* ewaldcof, double* vol)
{
PRINTMETHOD("gpu_scalar_sumrc");
    kPMEScalarSumRC(gpu, *ewaldcof, *vol);

#if 0
    gpu->pbXYZ_qt->Download();
    PMEDouble2* pDouble2 = gpu->pbXYZ_qt->_pSysData;   
    for (int i = 0; i < gpu->sim.fft_x_dim; i++)
    {
        for (int j = 0; j < gpu->sim.fft_y_dim; j++)
        {
            for (int k = 0; k < gpu->sim.fft_z_dim; k++)
            {
                printf("%3d %3d %3d %40.33f %40.33f\n", i, j, k, 
                    pDouble2[(k * gpu->sim.fft_y_dim + j) * gpu->sim.fft_x_dim + i].x,
                    pDouble2[(k * gpu->sim.fft_y_dim + j) * gpu->sim.fft_x_dim + i].y);
            }
        }
    }
    exit(-1);
#endif    
}


extern "C" void gpu_fft3drc_back_()
{
    PRINTMETHOD("gpu_fft3drc_back");
#ifdef use_DPDP
    cufftExecZ2D(gpu->backwardPlan, gpu->sim.pXYZ_qt, gpu->sim.pXYZ_q);
#else
    cufftExecC2R(gpu->backwardPlan, gpu->sim.pXYZ_qt, gpu->sim.pXYZ_q);    
#endif
#if 0
    gpu->pbXYZ_q->Download();
    for (int i = 0; i < gpu->sim.nfft1; i++)
    {
        for (int j = 0; j < gpu->sim.nfft2; j++)
        {
            for (int k = 0; k < gpu->sim.nfft3; k++)
            {
                printf("%3d %3d %3d %32.15f\n", i, j, k, gpu->pbXYZ_q->_pSysData[(k * gpu->sim.nfft2 + j) * gpu->sim.nfft1 + i]); 
            }
        }
    }
    exit(-1);   
#endif
    
#if 0    
    // Check virial and energies
    gpu->pbEnergyBuffer->Download();
    for (int i = 0; i < 12; i++)
    {
        unsigned long long int val  = gpu->pbVirialBuffer->_pSysData[i];
        PMEDouble dval;
        if (val >= 0x8000000000000000ull)
        {
            dval                    = -(PMEDouble)(val ^ 0xffffffffffffffffull) / ENERGYSCALE;
        }
        else
        {
            dval                    = (PMEDouble)val / ENERGYSCALE;
        }
        printf("%d %20.13lf\n", i, dval);
    }
    unsigned long long int val  = gpu->pbEnergyBuffer->_pSysData[9];
    PMEDouble dval;
    if (val >= 0x8000000000000000ull)
    {
        dval                    = -(PMEDouble)(val ^ 0xffffffffffffffffull) / ENERGYSCALE;
    }
    else
    {
        dval                    = (PMEDouble)val / ENERGYSCALE;
    }
    printf("Energy %20.13lf\n", dval);
    exit(-1);
#endif
}

extern "C" void gpu_grad_sum_()
{
PRINTMETHOD("gpu_grad_sum");
    kPMEGradSum(gpu);

#if 0
    gpu->pbForceBuffer->Download();
    gpu->pbImageIndex->Download();
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        for (int j = 0; j < gpu->sim.atoms; j++)
        {
            if (gpu->pbImageIndex->_pSysData[j + gpu->sim.stride * 5] == i)
                printf("%6d %20.15f %20.15f %20.15f\n", i, gpu->pbForceBuffer->_pSysData[j], gpu->pbForceBuffer->_pSysData[j + gpu->sim.stride], gpu->pbForceBuffer->_pSysData[j + gpu->sim.stride2]);
        }
    }
    exit(-1);
#endif    
}

extern "C" void gpu_self_(double* ee_plasma, double* ene)
{
    gpu->ee_plasma                              = *ee_plasma;
    gpu->self_energy                            = *ene;
}
extern "C" void gpu_vdw_correction_(double* ene)
{
    gpu->vdw_recip                              = *ene;
}

extern "C" void print_virial(char* message)
{
    static PMEDouble oldvirial[3] = {0.0, 0.0, 0.0};
    gpu->pbEnergyBuffer->Download();
    PMEDouble virial[ENERGYTERMS];
    for (int i = 0; i < 3; i++)
    {
        unsigned long long int val  = gpu->pbEnergyBuffer->_pSysData[i + VIRIALOFFSET];
        if (val >= 0x8000000000000000ull)
        {
            virial[i]               = -(PMEDouble)(val ^ 0xffffffffffffffffull) / ENERGYSCALE;
        }
        else
        {
            virial[i]               = (PMEDouble)val / ENERGYSCALE;
        }
    }
    printf("%20s %20.10f %20.10f %20.10f\n", message, virial[0] - oldvirial[0], virial[1] - oldvirial[1], virial[2] - oldvirial[2]);
    oldvirial[0] = virial[0];
    oldvirial[1] = virial[1];
    oldvirial[2] = virial[2];
    
}

#ifdef MPI
extern "C" void gpu_download_partial_forces()
{
PRINTMETHOD("gpu_download_partial_forces");
    cudaError_t status;
    PMEDouble *pDev                 = gpu->pbOutForce->_pDevData;
    PMEDouble *pSys                 = gpu->pbOutForce->_pSysData;
    //printf("%d hello %d %d 0x%016llx 0x%016llx\n", gpu->gpuID, gpu->sim.minReducedAtom, gpu->sim.reducedAtoms3, pDev, pSys);
    status = cudaMemcpy(pSys, pDev, gpu->sim.reducedAtoms3 * sizeof(PMEDouble), cudaMemcpyDeviceToHost); 
    RTERROR(status, "gpu_download_partial_forces: download failed");
}
#endif
extern "C" void gpu_pme_ene_(double* ewaldcof, double* vol, pme_pot_ene_rec* pEnergy, double virial[3], double ekcmt[3])
{
PRINTMETHOD("gpu_pme_ene");
    // Rebuild neighbor list
    gpu_build_neighbor_list_();

    // Clear forces        
    kClearForces(gpu);     

#ifdef MPI
    if (gpu->sim.minLocalAtom != gpu->sim.maxLocalAtom)   
        memset(gpu->pForceData1 + gpu->forceReceiveFirstAtom * 3, 0, (gpu->forceReceiveLastAtom - gpu->forceReceiveFirstAtom) * 3 * sizeof(PMEDouble));       
#endif
    

#if 0
    {
        gpu->pbImageIndex->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);        
        PMEDouble *pCrd                                 = gpu->pbImage->_pSysData; 
        gpu->pbImage->Download();     
     
        if (gpu->sim.pImageX != gpu->pbImage->_pDevData)
           pCrd                                        = gpu->pbImage->_pSysData + gpu->sim.stride3;

        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            unsigned int j = pImageAtomLookup[i];
            printf("%6d,%32.15f,%32.15f,%32.15f\n", i, pCrd[j], pCrd[j + gpu->sim.stride], pCrd[j + gpu->sim.stride2]);
        }
    }
#endif        


    if (gpu->ntf != 8)
    {

        // Local energy
#ifdef MPI
        if (gpu->sim.minLocalAtom != gpu->sim.maxLocalAtom)
        {
#endif
            kCalculateLocalEnergy(gpu);
            kCalculateCHARMMEnergy(gpu);
#ifdef MPI
        }
#endif
#ifdef MPI
        if (gpu->gpuID == 0)
#endif        
        {

            // Reciprocal energy
            kPMEGetGridWeights(gpu);               
            kPMEClearChargeGridBuffer(gpu);  
            kPMEFillChargeGridBuffer(gpu);
            kPMEReduceChargeGridBuffer(gpu);   
#if 0
    gpu->pbXYZ_q->Download();
    for (int i = 0; i < gpu->sim.nfft1; i++)
    {
        for (int j = 0; j < gpu->sim.nfft2; j++)
        {
            for (int k = 0; k < gpu->sim.nfft3; k++)
            {
                printf("%3d %3d %3d %20.14f\n", i, j, k, gpu->pbXYZ_q->_pSysData[(k * gpu->sim.nfft2 + j) * gpu->sim.nfft1 + i]); 
            }
        }
    }
    exit(-1);   
#endif            
            
#ifdef use_DPDP
            cufftExecD2Z(gpu->forwardPlan, gpu->sim.pXYZ_q, gpu->sim.pXYZ_qt);
#else
            cufftExecR2C(gpu->forwardPlan, gpu->sim.pXYZ_q, gpu->sim.pXYZ_qt);
#endif

#if 0
    gpu->pbXYZ_qt->Download();
    
    for (int i = 0; i < gpu->sim.fft_x_dim; i++)
    {
        for (int j = 0; j < gpu->sim.fft_y_dim; j++)
        {
            for (int k = 0; k < gpu->sim.fft_z_dim; k++)
            {
                printf("%3d %3d %3d %20.14f %20.14f\n", i, j, k, gpu->pbXYZ_qt->_pSysData[(k * gpu->sim.fft_y_dim + j) * gpu->sim.fft_x_dim + i].x, gpu->pbXYZ_qt->_pSysData[(k * gpu->sim.fft_y_dim + j) * gpu->sim.fft_x_dim + i].y); 
            }
        }
    }
    exit(-1);   
#endif            


            kPMEScalarSumRCEnergy(gpu, *ewaldcof, *vol);    
            
#if 0
    gpu->pbXYZ_qt->Download();
    
    for (int i = 0; i < gpu->sim.fft_x_dim; i++)
    {
        for (int j = 0; j < gpu->sim.fft_y_dim; j++)
        {
            for (int k = 0; k < gpu->sim.fft_z_dim; k++)
            {
                printf("%3d %3d %3d %20.14f %20.14f\n", i, j, k, gpu->pbXYZ_qt->_pSysData[(k * gpu->sim.fft_y_dim + j) * gpu->sim.fft_x_dim + i].x, gpu->pbXYZ_qt->_pSysData[(k * gpu->sim.fft_y_dim + j) * gpu->sim.fft_x_dim + i].y); 
            }
        }
    }
    exit(-1);   
#endif                 
            
#ifdef use_DPDP    
            cufftExecZ2D(gpu->backwardPlan, gpu->sim.pXYZ_qt, gpu->sim.pXYZ_q);
#else
            cufftExecC2R(gpu->backwardPlan, gpu->sim.pXYZ_qt, gpu->sim.pXYZ_q);
#endif

#if 0
    gpu->pbXYZ_q->Download();
    for (int i = 0; i < gpu->sim.nfft1; i++)
    {
        for (int j = 0; j < gpu->sim.nfft2; j++)
        {
            for (int k = 0; k < gpu->sim.nfft3; k++)
            {
                printf("%3d %3d %3d %20.14f\n", i, j, k, gpu->pbXYZ_q->_pSysData[(k * gpu->sim.nfft2 + j) * gpu->sim.nfft1 + i]); 
            }
        }
    }
    exit(-1);   
#endif       

            kPMEGradSum(gpu);        
        }
#ifdef MPI        
        else
        {
            kPMEGetGridWeights(gpu); 
        }
#endif

        // Direct energy and force reduction
#ifdef MPI        
        if (gpu->sim.minLocalAtom != gpu->sim.maxLocalAtom)
        {                
            kCalculatePMENonbondEnergy(gpu);    
            kReduceForces(gpu);
            if (gpu->sim.EPs > 0)
                kOrientForces(gpu);            
            kTransposeForces(gpu);      
        }
        
        // Distribute reciprocal forces
        if (gpu->gpuID == 0)
        {         
            for (int i = 1; i < gpu->nGpus; i++)
            {
                MPI_Put(&gpu->pbPMEForce->_pSysData[gpu->pPMEStart[i]], gpu->pPMELength[i], MPI_PMEFLOAT, i, 
                     gpu->pPMEStart[i] * sizeof(PMEFloat) + 2 * gpu->sim.stride3 * sizeof(PMEDouble), 
                     gpu->pPMELength[i], MPI_PMEFLOAT, gpu->MPIPMEForceWindow);    
            }
        }
        if (gpu->sim.minLocalAtom != gpu->sim.maxLocalAtom)
        {
            cudaThreadSynchronize();
            gpu_download_partial_forces();
        }

#else
        kCalculatePMENonbondEnergy(gpu);    
        kReduceForces(gpu); 
        if (gpu->sim.EPs > 0)
            kOrientForces(gpu);        
#if 0      
        gpu->pbForce->Download();
        for (int i = 0; i < gpu->sim.atoms; i++)
            printf("%06d %18.12f %18.12f %18.12f\n", i, gpu->pbForce->_pSysData[i], gpu->pbForce->_pSysData[i + gpu->sim.stride], gpu->pbForce->_pSysData[i + gpu->sim.stride2]);
        exit(-1);        
#endif

#if 0
        gpu->pbForceBuffer->Download();     
        for (int i = 0; i < gpu->sim.atoms; i++)
            for (int j = 0; j < gpu->sim.NLCellBuffers; j++)
                printf("%06d %06d %18.12f %18.12f %18.12f\n", i, j, gpu->pbForceBuffer->_pSysData[i + j * gpu->sim.stride3], gpu->pbForceBuffer->_pSysData[i + gpu->sim.stride + j * gpu->sim.stride3], gpu->pbForceBuffer->_pSysData[i + gpu->sim.stride2 + j * gpu->sim.stride3]);
        exit(-1);        
#endif
#endif   
    }
    
#if 0
        gpu->pbForce->Download();
        gpu->pbImageIndex->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            unsigned int j = pImageAtomLookup[i];
            printf("%6d,%32.15f,%32.15f,%32.15f\n", i, 
                gpu->pbForce->_pSysData[j],
                gpu->pbForce->_pSysData[gpu->sim.stride + j],
                gpu->pbForce->_pSysData[gpu->sim.stride2 + j]
            );
            
        }
#endif    
        
#ifdef MPI 
    // Transmit non-local forces to other nodes   
    for (int i = 0; i < gpu->forceSendNodes; i++)
    {
        //printf("MPA %d %d %d %d %d\n", gpu->gpuID, gpu->pForceSendNode[i], gpu->pOutForceSendStart[i], gpu->pForceSendStart[i], gpu->pForceSendLength[i]);
        MPI_Accumulate(&gpu->pbOutForce->_pSysData[gpu->pOutForceSendStart[i]], gpu->pForceSendLength[i], 
        MPI_PMEDOUBLE, gpu->pForceSendNode[i], 
        (gpu->pForceSendStart[i] + gpu->forceSendOffset) * sizeof(PMEDouble), gpu->pForceSendLength[i], 
        MPI_PMEDOUBLE, MPI_SUM, gpu->MPIPMEForceWindow);
    }
#endif    

#if 0
        gpu->pbForce->Download();
        gpu->pbImageIndex->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            unsigned int j = pImageAtomLookup[i];
            printf("%6d,%32.15f,%32.15f,%32.15f\n", i, 
                gpu->pbForce->_pSysData[j],
                gpu->pbForce->_pSysData[gpu->sim.stride + j],
                gpu->pbForce->_pSysData[gpu->sim.stride2 + j]
            );
            
        }
#endif    


        
    if (gpu->sim.ntp > 0)
    {
#ifdef MPI    
        if (gpu->gpuID == 1)
#endif
        {   
            kCalculateCOMKineticEnergy(gpu);
            kReduceCOMKineticEnergy(gpu);
        }
        kCalculateMolecularVirial(gpu);
   //     print_virial("molecular");
#if 0
        gpu->pbSolute->Download();
        gpu->pbSolvent->Download();
        double* pCOMX = gpu->pbSolute->_pSysData + 3 * 32;
        double* pCOMY = gpu->pbSolute->_pSysData + 4 * 32;
        double* pCOMZ = gpu->pbSolute->_pSysData + 5 * 32;
        for (int i = 0; i < gpu->sim.soluteMolecules; i++)
            printf("%6d %20.10f %20.10f %20.10f\n", i, pCOMX[i], pCOMY[i], pCOMZ[i]);
        pCOMX = gpu->pbSolvent->_pSysData + 4 * 7040;
        pCOMY = gpu->pbSolvent->_pSysData + 5 * 7040;
        pCOMZ = gpu->pbSolvent->_pSysData + 6 * 7040;
        for (int i = 0; i < gpu->sim.solventMolecules; i++)
            printf("%6d %20.10f %20.10f %20.10f\n", i + gpu->sim.soluteMolecules, pCOMX[i], pCOMY[i], pCOMZ[i]);
        exit(-1);
#endif
    }
#if 0
        gpu->pbForce->Download();
        gpu->pbImageIndex->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            unsigned int j = pImageAtomLookup[i];
            printf("%6d,%32.15f,%32.15f,%32.15f\n", i, 
                gpu->pbForce->_pSysData[j],
                gpu->pbForce->_pSysData[gpu->sim.stride + j],
                gpu->pbForce->_pSysData[gpu->sim.stride2 + j]
            );
            
        }
#endif    
    
    gpu->pbEnergyBuffer->Download(); 
    PMEDouble energy[ENERGYTERMS];
    pEnergy->total                  = 0.0;
    for (int i = 0; i < ENERGYTERMS; i++)
    {
        unsigned long long int val  = gpu->pbEnergyBuffer->_pSysData[i];
        if (val >= 0x8000000000000000ull)
        {
            energy[i]               = -(PMEDouble)(val ^ 0xffffffffffffffffull) / ENERGYSCALE;
        }
        else
        {
            energy[i]               = (PMEDouble)val / ENERGYSCALE;
        }
        if (i < VIRIALOFFSET)
            pEnergy->total         += energy[i];
        //printf("%6d %16.7f\n", i, energy[i]);
        //printf("%06d %6d %16.7f\n", gpu->gpuID, i, energy[i]);
    }
#ifdef MPI    
    if (gpu->gpuID == 0)
#endif
    pEnergy->total                 += gpu->vdw_recip + gpu->self_energy;
#ifdef MPI
    if (gpu->gpuID > 0)
        pEnergy->vdw_tot            = energy[1];
    else
#endif
        pEnergy->vdw_tot            = energy[1] + gpu->vdw_recip;
    pEnergy->vdw_dir                = energy[1];
    pEnergy->vdw_recip              = gpu->vdw_recip;
#ifdef MPI
    if (gpu->gpuID > 0)    
        pEnergy->elec_tot           = energy[9] + energy[10];
    else
#endif
        pEnergy->elec_tot           = energy[9] + energy[10] + gpu->self_energy;
    pEnergy->elec_dir               = energy[10];
    pEnergy->elec_recip             = energy[9];
    pEnergy->elec_nb_adjust         = 0.0;
    pEnergy->hbond                  = 0.0;
    pEnergy->bond                   = energy[3];
    pEnergy->angle                  = energy[4];
    pEnergy->dihedral               = energy[5];
    pEnergy->vdw_14                 = energy[7];
    pEnergy->elec_14                = energy[6];
    pEnergy->restraint              = energy[8];
    pEnergy->angle_ub               = energy[11];
    pEnergy->imp                    = energy[12];
    pEnergy->cmap                   = energy[13];
    
    // printf("E: %20.10f %20.10f %20.10f\n", energy[9], energy[10], energy[1]);

    // Grab virial if needed
    if (gpu->sim.ntp > 0)
    {

        virial[0]                   = 0.5 * energy[VIRIALOFFSET + 0];
        virial[1]                   = 0.5 * energy[VIRIALOFFSET + 1];
        virial[2]                   = 0.5 * energy[VIRIALOFFSET + 2];
#ifdef MPI
        if (gpu->gpuID == 0)
#endif 
        {      
            virial[0]              -= 0.5 * (gpu->ee_plasma + 2.0 * gpu->vdw_recip);
            virial[1]              -= 0.5 * (gpu->ee_plasma + 2.0 * gpu->vdw_recip);
            virial[2]              -= 0.5 * (gpu->ee_plasma + 2.0 * gpu->vdw_recip);
        }
        ekcmt[0]                    = energy[VIRIALOFFSET + 3];
        ekcmt[1]                    = energy[VIRIALOFFSET + 4];
        ekcmt[2]                    = energy[VIRIALOFFSET + 5];

       //printf("VE%3d: %20.10f %20.10f %20.10f | %20.10f %20.10f %20.10f\n", gpu->gpuID, energy[VIRIALOFFSET + 0], energy[VIRIALOFFSET + 1], energy[VIRIALOFFSET + 2], ekcmt[0], ekcmt[1], ekcmt[2]);
       //printf("VE %20.10f %20.10f %20.10f | %20.10f %20.10f %20.10f\n", energy[VIRIALOFFSET + 0], energy[VIRIALOFFSET + 1], energy[VIRIALOFFSET + 2], ekcmt[0], ekcmt[1], ekcmt[2]);
       //printf("E %20.10f %20.10f %20.10f\n", pEnergy->elec_dir, pEnergy->elec_recip, pEnergy->vdw_dir);
    }
    

#if 0
    gpu->pbImageIndex->Download();
    unsigned int* pImageAtomLookup                      = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
    for (int i = 0; i < 32; i++)
    {
        printf("%6d %6d\n", i, pImageAtomLookup[i]);
    }
#endif
}

extern "C" void gpu_pme_force_(double* ewaldcof, double* vol, double virial[3], double ekcmt[3])
{
PRINTMETHOD("gpu_pme_force");
    // Rebuild neighbor list
    gpu_build_neighbor_list_(); 

    // Clear forces
#ifdef MPI    
    if (gpu->sim.minLocalAtom != gpu->sim.maxLocalAtom)
#endif
        kClearForces(gpu); 
#ifdef MPI
    else if (gpu->sim.ntp > 0)
       cudaMemset(gpu->sim.pEnergyBuffer, 0, 8 * ENERGYTERMS);
#endif

#ifdef MPI
    if (gpu->sim.minLocalAtom != gpu->sim.maxLocalAtom)   
        memset(gpu->pForceData1 + gpu->forceReceiveFirstAtom * 3, 0, (gpu->forceReceiveLastAtom - gpu->forceReceiveFirstAtom) * 3 * sizeof(PMEDouble));       
#endif

    
    if (gpu->ntf != 8)
    {

        // Local forces
#ifdef MPI
        if (gpu->sim.minLocalAtom != gpu->sim.maxLocalAtom)
        {
#endif
            kCalculateLocalForces(gpu);
            kCalculateCHARMMForces(gpu);
#ifdef MPI
        }
#endif
#ifdef MPI
        if (gpu->gpuID == 0)
#endif        
        {  
            // Reciprocal forces
            kPMEGetGridWeights(gpu);          
            kPMEClearChargeGridBuffer(gpu);  
            kPMEFillChargeGridBuffer(gpu);
            kPMEReduceChargeGridBuffer(gpu);            
#ifdef use_DPDP   
            cufftExecD2Z(gpu->forwardPlan, gpu->sim.pXYZ_q, gpu->sim.pXYZ_qt);
#else
            cufftExecR2C(gpu->forwardPlan, gpu->sim.pXYZ_q, gpu->sim.pXYZ_qt);
#endif 

            kPMEScalarSumRC(gpu, *ewaldcof, *vol);
    
#ifdef use_DPDP
            cufftExecZ2D(gpu->backwardPlan, gpu->sim.pXYZ_qt, gpu->sim.pXYZ_q);
#else
            cufftExecC2R(gpu->backwardPlan, gpu->sim.pXYZ_qt, gpu->sim.pXYZ_q);
#endif
            kPMEGradSum(gpu);         
        }
#ifdef MPI
        else
        {
            kPMEGetGridWeights(gpu);
        }
#endif

        // Direct energy and force reduction
#ifdef MPI    
        if (gpu->sim.minLocalAtom != gpu->sim.maxLocalAtom)
        {
            kCalculatePMENonbondForces(gpu);                  
            kReduceForces(gpu);    
            if (gpu->sim.EPs > 0)
                kOrientForces(gpu);                
            kTransposeForces(gpu);                   
        }

        // Distribute reciprocal forces
        if (gpu->gpuID == 0)
        {
            for (int i = 1; i < gpu->nGpus; i++)
            {
                MPI_Put(&gpu->pbPMEForce->_pSysData[gpu->pPMEStart[i]], gpu->pPMELength[i], MPI_PMEFLOAT, i, 
                     gpu->pPMEStart[i] * sizeof(PMEFloat) + 2 * gpu->sim.stride3 * sizeof(PMEDouble), 
                     gpu->pPMELength[i], MPI_PMEFLOAT, gpu->MPIPMEForceWindow);                                     
            }
        }
        
        if (gpu->sim.minLocalAtom != gpu->sim.maxLocalAtom)
        {
            cudaThreadSynchronize();
            gpu_download_partial_forces();
        }       
#else
        kCalculatePMENonbondForces(gpu);    
        kReduceForces(gpu); 
        if (gpu->sim.EPs > 0)
            kOrientForces(gpu);    
#endif   
    }
    


#ifdef MPI 
    // Transmit non-local forces to other nodes   
    for (int i = 0; i < gpu->forceSendNodes; i++)
    {
       // printf("MPA %d %d %d %d\n", gpu->gpuID, gpu->pForceSendNode[i], gpu->pForceSendStart[i], gpu->pForceSendLength[i]);
        MPI_Accumulate(&gpu->pbOutForce->_pSysData[gpu->pOutForceSendStart[i]], gpu->pForceSendLength[i], 
        MPI_PMEDOUBLE, gpu->pForceSendNode[i], 
        (gpu->pForceSendStart[i] + gpu->forceSendOffset) * sizeof(PMEDouble), gpu->pForceSendLength[i], 
        MPI_PMEDOUBLE, MPI_SUM, gpu->MPIPMEForceWindow);
    }
#endif    
   
    // Grab virial and ekcmt if needed
    if (gpu->sim.ntp > 0)
    {
#ifdef MPI
        if (gpu->gpuID == 1)
        {
#endif
            kCalculateCOMKineticEnergy(gpu);
            kReduceCOMKineticEnergy(gpu);
#ifdef MPI
        }
#endif
        kCalculateMolecularVirial(gpu);    

        gpu->pbEnergyBuffer->Download();
        for(int i = 0; i < 6; i++)
        {
            unsigned long long int val  = gpu->pbEnergyBuffer->_pSysData[VIRIALOFFSET + i];
            PMEDouble dval;
            if (val >= 0x8000000000000000ull)
            {
                dval                    = -(PMEDouble)(val ^ 0xffffffffffffffffull) / ENERGYSCALE;
            }
            else
            {
                dval                    = (PMEDouble)val / ENERGYSCALE;
            }        
        
            if (i < 3)
                virial[i]               = 0.5 * dval;
            else
                ekcmt[i - 3]            = dval;
        }
#ifdef MPI
        if (gpu->gpuID == 0)
#endif 
        {      
            virial[0]                  -= 0.5 * (gpu->ee_plasma + 2.0 * gpu->vdw_recip);
            virial[1]                  -= 0.5 * (gpu->ee_plasma + 2.0 * gpu->vdw_recip);
            virial[2]                  -= 0.5 * (gpu->ee_plasma + 2.0 * gpu->vdw_recip);
        }


        //printf("VF%3d: %20.10f %20.10f %20.10f | %20.10f %20.10f %20.10f\n", gpu->gpuID, virial[0], virial[1], virial[2], ekcmt[0], ekcmt[1], ekcmt[2]);
    }    
}

#ifdef MPI
#ifdef CUDA_MPI
extern "C" void P2P_Broadcast(void* pSrc, void** pDst, int length, int displacement, size_t size)
{
    int offset                          = (displacement * size) / sizeof(char);
    int sendlength                      = length * size;
    char* pSrcChar                      = (char*)pSrc + offset;
    for (int i = 0; i < gpu->nGpus; i++)
    {
        if (i != gpu->gpuID)
        {
            char* pDstChar              = (char*)pDst[i] + offset;
            cudaError_t status          = cudaMemcpy((void*)pSrcChar, (void*)pDstChar, sendlength, cudaMemcpyDeviceToDevice);
            RTERROR(status, "P2P_Broadcast failed.\n");
        }
    }
    cudaDeviceSynchronize();
    MPI_Barrier(gpu->comm);
}
#endif

extern "C" void gpu_gather_pme_forces_()
{
PRINTMETHOD("gpu_gather_pme_forces");
    static int step = 0;
    // Initialize force vectors
    //printf("Node %d, mina/maxa %d/%d\n", gpu->gpuID, gpu->sim.minLocalAtom, gpu->sim.maxLocalAtom);


    // Wait for all local forces to arrive
    MPI_Win_fence(0, gpu->MPIPMEForceWindow);
#if 0
        {
        char buff[512];
        sprintf(buff, "../npforces%d_%d.out", gpu->gpuID, step);
        FILE* fp = fopen(buff, "w");
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbOutForce->_pSysData[i * 3],
                                                            gpu->pbOutForce->_pSysData[i * 3 + 1], 
                                                            gpu->pbOutForce->_pSysData[i * 3 + 2]);
            fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pForceData0[i * 3],
                                                            gpu->pForceData0[i * 3 + 1], 
                                                            gpu->pForceData0[i * 3 + 2]);
            if (gpu->gpuID != 0)
                fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pPMEForceData[i * 3],
                                                            gpu->pPMEForceData[i * 3 + 1], 
                                                            gpu->pPMEForceData[i * 3 + 2]);   
            else
                fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbPMEForce->_pSysData[i * 3],
                                                            gpu->pbPMEForce->_pSysData[i * 3 + 1], 
                                                            gpu->pbPMEForce->_pSysData[i * 3 + 2]);                                                                        
                                                      
        }
        fclose(fp);
        }
#endif                


#if 0
        {
        char buff[512];
        sprintf(buff, "../npforces%d_%d.out", gpu->gpuID, step);
        FILE* fp = fopen(buff, "w");
        gpu->pbForce->Download();
        for (int i = gpu->sim.minReducedAtom; i != gpu->sim.maxReducedAtom; i = ((i + 1) % gpu->sim.paddedNumberOfAtoms))
        {
            fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbForce->_pSysData[i],
                                                            gpu->pbForce->_pSysData[i + gpu->sim.stride], 
                                                            gpu->pbForce->_pSysData[i + gpu->sim.stride2]);
            int j                                       = i - gpu->sim.minReducedAtom;
            if (j < 0)
                j                                      += gpu->sim.paddedNumberOfAtoms;
            fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbOutForce->_pSysData[j * 3],
                                                            gpu->pbOutForce->_pSysData[j * 3 + 1], 
                                                            gpu->pbOutForce->_pSysData[j * 3 + 2]);                                                            
        }
        fclose(fp);
        }
#endif                

        // Copy PME and received force data
        int atoms                                       = gpu->forceReceiveFirstAtom - gpu->sim.minLocalAtom;
        int length                                      = atoms >> 3;       
        PMEFloat* pPMEForce                             = ((gpu->gpuID != 0) ? gpu->pPMEForceData : gpu->pbPMEForce->_pSysData) + gpu->sim.minLocalAtom * 3;                
        PMEDouble* pForce                               = gpu->pbOutForce->_pSysData + (gpu->sim.minLocalAtom - gpu->sim.minReducedAtom) * 3;
        PMEDouble* pInForce                             = gpu->pbInForce->_pSysData + gpu->sim.minLocalAtom * 3;    
        //printf("X %d %d %d %d %d\n", gpu->gpuID, atoms, length, gpu->forceReceiveFirstAtom, gpu->sim.minLocalAtom);
        int i;
        for (i = 0; i < length; i++)
        {
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;             
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;             
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;                                     
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;             
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;             
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;             
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;                                     
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;        
        }
        i                                             <<= 3;
        for (;i < atoms; i++)
        {
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;             
        }                                       
        
        atoms                                           = gpu->forceReceiveLastAtom - gpu->forceReceiveFirstAtom;
        length                                          = atoms >> 3;            
        PMEDouble* pReceivedForce                       = gpu->pForceData0 + gpu->forceReceiveFirstAtom * 3;
        //printf("Y %d %d %d %d %d\n", gpu->gpuID, atoms, length, gpu->forceReceiveLastAtom, gpu->forceReceiveFirstAtom);        
        for (i = 0; i < length; i++)
        {
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;   
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;                          
        }
        i                                             <<= 3;
        for (;i < atoms; i++)
        {
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++ + *pReceivedForce++;
        }                              
        
        atoms                                           = gpu->sim.maxLocalAtom - gpu->forceReceiveLastAtom;
        length                                          = atoms >> 3;
        //printf("Z %d %d %d %d %d\n", gpu->gpuID, atoms, length, gpu->sim.maxLocalAtom, gpu->forceReceiveLastAtom);                
        for (i = 0; i < length; i++)
        {
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++; 
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;                                                
        }
        i                                             <<= 3;
        for (;i < atoms; i++)
        {
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
            *pInForce++                                 = *pForce++ + *pPMEForce++;
        }                              
        
        

        // Update force pointers and offset
        gpu->forceSendOffset                            = gpu->sim.stride3 - gpu->forceSendOffset;
        pForce                                          = gpu->pForceData0;
        gpu->pForceData0                                = gpu->pForceData1;
        gpu->pForceData1                                = pForce;
    
     

#if 0
        {
        char buff[512];
        sprintf(buff, "../nbforces%d_%d.out", gpu->gpuID, step);
        FILE* fp = fopen(buff, "w");
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            if ((i >= gpu->sim.minLocalAtom) && (i < gpu->sim.maxLocalAtom))
                fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbInForce->_pSysData[i],
                                                            gpu->pbInForce->_pSysData[i + gpu->sim.stride], 
                                                            gpu->pbInForce->_pSysData[i + gpu->sim.stride2]);   
            else
                fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, 0.0,
                                                            0.0, 
                                                            0.0);                                                       
        }
        fclose(fp);
        }
#endif        
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_PMEDOUBLE, gpu->pbInForce->_pSysData, gpu->pAllGathervRecvCountAoS, gpu->pAllGathervRecvDisplAoS, MPI_PMEDOUBLE, gpu->comm);
       
#if 0
    {
    char buff[512];
    sprintf(buff, "../naforces%d_%d.out", gpu->gpuID, step);
    FILE* fp = fopen(buff, "w");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbInForce->_pSysData[i * 3],
                                                        gpu->pbInForce->_pSysData[i * 3 + 1], 
                                                        gpu->pbInForce->_pSysData[i * 3 + 2]);             
    }
    fclose(fp);
  //  exit(-1);
    }
#endif    
    step++;
}
#endif






extern "C" void gpu_ips_ene_(double* vol, pme_pot_ene_rec* pEnergy, double virial[3], double ekcmt[3])
{
PRINTMETHOD("gpu_ips_ene");
    // Rebuild neighbor list
    gpu_build_neighbor_list_();

    // Clear forces        
    kClearForces(gpu);     

#ifdef MPI
    memset(gpu->pForceData1 + gpu->forceReceiveFirstAtom * 3, 0, (gpu->forceReceiveLastAtom - gpu->forceReceiveFirstAtom) * 3 * sizeof(PMEDouble));       
#endif
    

#if 0
    {
        gpu->pbImageIndex->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);        
        PMEDouble *pCrd                                 = gpu->pbImage->_pSysData; 
        gpu->pbImage->Download();     
     
        if (gpu->sim.pImageX != gpu->pbImage->_pDevData)
           pCrd                                        = gpu->pbImage->_pSysData + gpu->sim.stride3;

        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            unsigned int j = pImageAtomLookup[i];
            printf("%6d,%32.15f,%32.15f,%32.15f\n", i, pCrd[j], pCrd[j + gpu->sim.stride], pCrd[j + gpu->sim.stride2]);
        }
    }
#endif        


    if (gpu->ntf != 8)
    {

        // Local energy
        kCalculateLocalEnergy(gpu);
        kCalculateCHARMMEnergy(gpu);
        kNLCalculateCellCoordinates(gpu);   
        kCalculateIPSNonbondEnergy(gpu);    
        kReduceForces(gpu);
        if (gpu->sim.EPs > 0)
            kOrientForces(gpu);        
#ifdef MPI             
        kTransposeForces(gpu);         
        cudaThreadSynchronize();
        gpu_download_partial_forces();
#endif
    }
    
#if 0
        gpu->pbForce->Download();
        gpu->pbImageIndex->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            unsigned int j = pImageAtomLookup[i];
            printf("%6d,%32.15f,%32.15f,%32.15f\n", i, 
                gpu->pbForce->_pSysData[j],
                gpu->pbForce->_pSysData[gpu->sim.stride + j],
                gpu->pbForce->_pSysData[gpu->sim.stride2 + j]
            );
            
        }
#endif    
        
        
#if 0
    {
    char buff[512];
    gpu->pbForce->Download();
    sprintf(buff, "../naforces.out");
    FILE* fp = fopen(buff, "w");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbForce->_pSysData[i],
                                                        gpu->pbForce->_pSysData[i + gpu->sim.stride], 
                                                        gpu->pbForce->_pSysData[i + gpu->sim.stride2]);             
    }
    fclose(fp);
    }
#endif    

        
        
#ifdef MPI 
    // Transmit non-local forces to other nodes   
    for (int i = 0; i < gpu->forceSendNodes; i++)
    {
        //printf("MPA %d %d %d %d %d\n", gpu->gpuID, gpu->pForceSendNode[i], gpu->pOutForceSendStart[i], gpu->pForceSendStart[i], gpu->pForceSendLength[i]);
        MPI_Accumulate(&gpu->pbOutForce->_pSysData[gpu->pOutForceSendStart[i]], gpu->pForceSendLength[i], 
        MPI_PMEDOUBLE, gpu->pForceSendNode[i], 
        (gpu->pForceSendStart[i] + gpu->forceSendOffset) * sizeof(PMEDouble), gpu->pForceSendLength[i], 
        MPI_PMEDOUBLE, MPI_SUM, gpu->MPIPMEForceWindow);
    }
#endif    

#if 0
        gpu->pbForce->Download();
        gpu->pbImageIndex->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            unsigned int j = pImageAtomLookup[i];
            printf("%6d,%32.15f,%32.15f,%32.15f\n", i, 
                gpu->pbForce->_pSysData[j],
                gpu->pbForce->_pSysData[gpu->sim.stride + j],
                gpu->pbForce->_pSysData[gpu->sim.stride2 + j]
            );
            
        }
#endif    


        
    if (gpu->sim.ntp > 0)
    {
#ifdef MPI    
        if (gpu->gpuID == 1)
#endif
        {   
            kCalculateCOMKineticEnergy(gpu);
            kReduceCOMKineticEnergy(gpu);
        }
        kCalculateMolecularVirial(gpu);
   //     print_virial("molecular");
#if 0
        gpu->pbSolute->Download();
        gpu->pbSolvent->Download();
        double* pCOMX = gpu->pbSolute->_pSysData + 3 * 32;
        double* pCOMY = gpu->pbSolute->_pSysData + 4 * 32;
        double* pCOMZ = gpu->pbSolute->_pSysData + 5 * 32;
        for (int i = 0; i < gpu->sim.soluteMolecules; i++)
            printf("%6d %20.10f %20.10f %20.10f\n", i, pCOMX[i], pCOMY[i], pCOMZ[i]);
        pCOMX = gpu->pbSolvent->_pSysData + 4 * 7040;
        pCOMY = gpu->pbSolvent->_pSysData + 5 * 7040;
        pCOMZ = gpu->pbSolvent->_pSysData + 6 * 7040;
        for (int i = 0; i < gpu->sim.solventMolecules; i++)
            printf("%6d %20.10f %20.10f %20.10f\n", i + gpu->sim.soluteMolecules, pCOMX[i], pCOMY[i], pCOMZ[i]);
        exit(-1);
#endif
    }
#if 0
        gpu->pbForce->Download();
        gpu->pbImageIndex->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            unsigned int j = pImageAtomLookup[i];
            printf("%6d,%32.15f,%32.15f,%32.15f\n", i, 
                gpu->pbForce->_pSysData[j],
                gpu->pbForce->_pSysData[gpu->sim.stride + j],
                gpu->pbForce->_pSysData[gpu->sim.stride2 + j]
            );
            
        }
#endif    
    
    gpu->pbEnergyBuffer->Download(); 
    PMEDouble energy[ENERGYTERMS];
    pEnergy->total                  = 0.0;
    for (int i = 0; i < ENERGYTERMS; i++)
    {
        unsigned long long int val  = gpu->pbEnergyBuffer->_pSysData[i];
        if (val >= 0x8000000000000000ull)
        {
            energy[i]               = -(PMEDouble)(val ^ 0xffffffffffffffffull) / ENERGYSCALE;
        }
        else
        {
            energy[i]               = (PMEDouble)val / ENERGYSCALE;
        }
        if (i < VIRIALOFFSET)
            pEnergy->total         += energy[i];
        //printf("%6d %16.7f\n", i, energy[i]);
        //printf("%06d %6d %16.7f\n", gpu->gpuID, i, energy[i]);
    }
#ifdef MPI    
    if (gpu->gpuID == 0)
#endif
    pEnergy->total                 += gpu->sim.EIPSEL + gpu->sim.EIPSNB + gpu->sim.eipssnb + gpu->sim.eipssel;
#ifdef MPI
    if (gpu->gpuID > 0)
        pEnergy->vdw_tot            = energy[1];
    else
#endif
    pEnergy->vdw_tot                = energy[1] + gpu->sim.EIPSNB + + gpu->sim.eipssnb;
    pEnergy->vdw_dir                = pEnergy->vdw_tot;
#ifdef MPI
    if (gpu->gpuID > 0)    
        pEnergy->elec_tot           = energy[10];
    else
#endif
    pEnergy->elec_tot               = energy[10] + gpu->sim.EIPSEL + gpu->sim.eipssel;
    pEnergy->elec_dir               = pEnergy->elec_tot;
    pEnergy->elec_recip             = 0.0;
    pEnergy->elec_nb_adjust         = 0.0;
    pEnergy->hbond                  = 0.0;
    pEnergy->bond                   = energy[3];
    pEnergy->angle                  = energy[4];
    pEnergy->dihedral               = energy[5];
    pEnergy->vdw_14                 = energy[7];
    pEnergy->elec_14                = energy[6];
    pEnergy->restraint              = energy[8];
    pEnergy->angle_ub               = energy[11];
    pEnergy->imp                    = energy[12];
    pEnergy->cmap                   = energy[13];
    
    // printf("E: %20.10f %20.10f %20.10f\n", energy[9], energy[10], energy[1]);

    // Grab virial if needed
    if (gpu->sim.ntp > 0)
    {

        virial[0]                   = 0.5 * energy[VIRIALOFFSET + 0];
        virial[1]                   = 0.5 * energy[VIRIALOFFSET + 1];
        virial[2]                   = 0.5 * energy[VIRIALOFFSET + 2];
        ekcmt[0]                    = energy[VIRIALOFFSET + 3];
        ekcmt[1]                    = energy[VIRIALOFFSET + 4];
        ekcmt[2]                    = energy[VIRIALOFFSET + 5];

#ifdef MPI
        if (gpu->gpuID == 0)
#endif 
        {      
            virial[0]              += gpu->sim.virips;
            virial[1]              += gpu->sim.virips;
            virial[2]              += gpu->sim.virips;
        }        

       //printf("VE%3d: %20.10f %20.10f %20.10f | %20.10f %20.10f %20.10f\n", gpu->gpuID, energy[VIRIALOFFSET + 0], energy[VIRIALOFFSET + 1], energy[VIRIALOFFSET + 2], ekcmt[0], ekcmt[1], ekcmt[2]);
       //printf("VE%3d %20.10f %20.10f %20.10f | %20.10f %20.10f %20.10f\n", gpu->gpuID, virial[0], virial[1], virial[2], ekcmt[0], ekcmt[1], ekcmt[2]);
       //printf("E %20.10f %20.10f %20.10f\n", pEnergy->elec_dir, pEnergy->elec_recip, pEnergy->vdw_dir);
    }
    

#if 0
    gpu->pbImageIndex->Download();
    unsigned int* pImageAtomLookup                      = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
    for (int i = 0; i < 32; i++)
    {
        printf("%6d %6d\n", i, pImageAtomLookup[i]);
    }
#endif
}

extern "C" void gpu_ips_force_(double* vol, double virial[3], double ekcmt[3])
{
PRINTMETHOD("gpu_ips_force");
    // Rebuild neighbor list
    gpu_build_neighbor_list_(); 

    // Clear forces
    kClearForces(gpu); 

#ifdef MPI  
    memset(gpu->pForceData1 + gpu->forceReceiveFirstAtom * 3, 0, (gpu->forceReceiveLastAtom - gpu->forceReceiveFirstAtom) * 3 * sizeof(PMEDouble));       
#endif

    
    if (gpu->ntf != 8)
    {

        // Local forces
        kCalculateLocalForces(gpu);
        kCalculateCHARMMForces(gpu);
        kNLCalculateCellCoordinates(gpu);  

        // Direct energy and force reduction
        kCalculateIPSNonbondForces(gpu);                  
        kReduceForces(gpu);    
        if (gpu->sim.EPs > 0)
            kOrientForces(gpu);  
#ifdef MPI    
        kTransposeForces(gpu);                   
        cudaThreadSynchronize();
        gpu_download_partial_forces();
#endif             
    }
    
    
#ifdef MPI 
    // Transmit non-local forces to other nodes   
    for (int i = 0; i < gpu->forceSendNodes; i++)
    {
        //printf("MPA %d %d %d %d %d\n", gpu->gpuID, gpu->pForceSendNode[i], gpu->pOutForceSendStart[i], gpu->pForceSendStart[i], gpu->pForceSendLength[i]);
        MPI_Accumulate(&gpu->pbOutForce->_pSysData[gpu->pOutForceSendStart[i]], gpu->pForceSendLength[i], 
        MPI_PMEDOUBLE, gpu->pForceSendNode[i], 
        (gpu->pForceSendStart[i] + gpu->forceSendOffset) * sizeof(PMEDouble), gpu->pForceSendLength[i], 
        MPI_PMEDOUBLE, MPI_SUM, gpu->MPIPMEForceWindow);
    }
#endif     
   
    // Grab virial and ekcmt if needed
    if (gpu->sim.ntp > 0)
    {
#ifdef MPI
        if (gpu->gpuID == 1)
#endif
        {
            kCalculateCOMKineticEnergy(gpu);
            kReduceCOMKineticEnergy(gpu);
        }
        kCalculateMolecularVirial(gpu);    

        gpu->pbEnergyBuffer->Download();
        for(int i = 0; i < 6; i++)
        {
            unsigned long long int val  = gpu->pbEnergyBuffer->_pSysData[VIRIALOFFSET + i];
            PMEDouble dval;
            if (val >= 0x8000000000000000ull)
            {
                dval                    = -(PMEDouble)(val ^ 0xffffffffffffffffull) / ENERGYSCALE;
            }
            else
            {
                dval                    = (PMEDouble)val / ENERGYSCALE;
            }        
        
            if (i < 3)
                virial[i]               = 0.5 * dval;
            else
                ekcmt[i - 3]            = dval;
        }
        
#ifdef MPI
        if (gpu->gpuID == 0)
#endif 
        {      
            virial[0]                  += gpu->sim.virips;
            virial[1]                  += gpu->sim.virips;
            virial[2]                  += gpu->sim.virips;
        }   

       //printf("V%3d %20.10f %20.10f %20.10f | %20.10f %20.10f %20.10f| %20.10f\n", gpu->gpuID, virial[0], virial[1], virial[2], ekcmt[0], ekcmt[1], ekcmt[2], gpu->sim.virips);
       //printf("V %20.10f %20.10f %20.10f | %20.10f %20.10f %20.10f| %20.10f\n", virial[0], virial[1], virial[2], ekcmt[0], ekcmt[1], ekcmt[2], gpu->sim.virips);
    }    
}











#ifdef MPI
extern "C" void gpu_gather_ips_forces_()
{
PRINTMETHOD("gpu_gather_ips_forces");
    static int step = 0;
    // Initialize force vectors
    //printf("Node %d, mina/maxa %d/%d\n", gpu->gpuID, gpu->sim.minLocalAtom, gpu->sim.maxLocalAtom);


    // Wait for all local forces to arrive
    MPI_Win_fence(0, gpu->MPIPMEForceWindow);            
#if 0
        {
        char buff[512];
        sprintf(buff, "../npforces%d_%d.out", gpu->gpuID, step);
        FILE* fp = fopen(buff, "w");
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbOutForce->_pSysData[i * 3],
                                                            gpu->pbOutForce->_pSysData[i * 3 + 1], 
                                                            gpu->pbOutForce->_pSysData[i * 3 + 2]);
            fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pForceData0[i * 3],
                                                            gpu->pForceData0[i * 3 + 1], 
                                                            gpu->pForceData0[i * 3 + 2]);                                                                     
        }
        fclose(fp);
        }
#endif                


#if 0
        {
        char buff[512];
        sprintf(buff, "../npforces%d_%d.out", gpu->gpuID, step);
        FILE* fp = fopen(buff, "w");
        gpu->pbForce->Download();
        for (int i = gpu->sim.minReducedAtom; i != gpu->sim.maxReducedAtom; i = ((i + 1) % gpu->sim.paddedNumberOfAtoms))
        {
            fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbForce->_pSysData[i],
                                                            gpu->pbForce->_pSysData[i + gpu->sim.stride], 
                                                            gpu->pbForce->_pSysData[i + gpu->sim.stride2]);
            int j                                       = i - gpu->sim.minReducedAtom;
            if (j < 0)
                j                                      += gpu->sim.paddedNumberOfAtoms;
            fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbOutForce->_pSysData[j * 3],
                                                            gpu->pbOutForce->_pSysData[j * 3 + 1], 
                                                            gpu->pbOutForce->_pSysData[j * 3 + 2]);                                                            
        }
        fclose(fp);
        }
#endif                

        // Copy received force data
        int atoms                                       = gpu->forceReceiveFirstAtom - gpu->sim.minLocalAtom;
        int length                                      = atoms >> 3;                    
        PMEDouble* pForce                               = gpu->pbOutForce->_pSysData + (gpu->sim.minLocalAtom - gpu->sim.minReducedAtom) * 3;
        PMEDouble* pInForce                             = gpu->pbInForce->_pSysData + gpu->sim.minLocalAtom * 3;    
        //printf("X %d %d %d %d %d\n", gpu->gpuID, atoms, length, gpu->forceReceiveFirstAtom, gpu->sim.minLocalAtom);
        int i;
        for (i = 0; i < length; i++)
        {
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;             
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;             
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;                                     
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;             
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;             
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;             
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;                                     
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;        
        }
        i                                             <<= 3;
        for (;i < atoms; i++)
        {
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;             
        }                                       
        
        atoms                                           = gpu->forceReceiveLastAtom - gpu->forceReceiveFirstAtom;
        length                                          = atoms >> 3;            
        PMEDouble* pReceivedForce                       = gpu->pForceData0 + gpu->forceReceiveFirstAtom * 3;
        //printf("Y %d %d %d %d %d\n", gpu->gpuID, atoms, length, gpu->forceReceiveLastAtom, gpu->forceReceiveFirstAtom);        
        for (i = 0; i < length; i++)
        {
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;   
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;                          
        }
        i                                             <<= 3;
        for (;i < atoms; i++)
        {
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
            *pInForce++                                 = *pForce++ + *pReceivedForce++;
        }                              
        
        atoms                                           = gpu->sim.maxLocalAtom - gpu->forceReceiveLastAtom;
        length                                          = atoms >> 3;
        //printf("Z %d %d %d %d %d\n", gpu->gpuID, atoms, length, gpu->sim.maxLocalAtom, gpu->forceReceiveLastAtom);                
        for (i = 0; i < length; i++)
        {
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++; 
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;                                                
        }
        i                                             <<= 3;
        for (;i < atoms; i++)
        {
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
            *pInForce++                                 = *pForce++;
        }                                     
       

        // Update force pointers and offset
        gpu->forceSendOffset                            = gpu->sim.stride3 - gpu->forceSendOffset;
        pForce                                          = gpu->pForceData0;
        gpu->pForceData0                                = gpu->pForceData1;
        gpu->pForceData1                                = pForce;
    
     

#if 0
        {
        char buff[512];
        sprintf(buff, "../nbforces%d_%d.out", gpu->gpuID, step);
        FILE* fp = fopen(buff, "w");
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            if ((i >= gpu->sim.minLocalAtom) && (i < gpu->sim.maxLocalAtom))
                fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbInForce->_pSysData[i * 3],
                                                            gpu->pbInForce->_pSysData[i * 3 + 1], 
                                                            gpu->pbInForce->_pSysData[i * 3 + 2]);   
            else
                fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, 0.0,
                                                            0.0, 
                                                            0.0);                                                       
        }
        fclose(fp);
        }
#endif        
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_PMEDOUBLE, gpu->pbInForce->_pSysData, gpu->pAllGathervRecvCountAoS, gpu->pAllGathervRecvDisplAoS, MPI_PMEDOUBLE, gpu->comm);
       
#if 0
    {
    char buff[512];
    sprintf(buff, "../naforces%d_%d.out", gpu->gpuID, step);
    FILE* fp = fopen(buff, "w");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbInForce->_pSysData[i * 3],
                                                        gpu->pbInForce->_pSysData[i * 3 + 1], 
                                                        gpu->pbInForce->_pSysData[i * 3 + 2]);             
    }
    fclose(fp);
    }
#endif    
    step++;
}
#endif

extern "C" void gpu_pressure_scale_(double factor[])
{
    NTPData* pNTPData                                       = gpu->pbNTPData->_pSysData;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            pNTPData->last_recip[i * 3 + j]                 = pNTPData->recip[i * 3 + j];
            pNTPData->recip[i * 3 + j]                     /= factor[i];
            pNTPData->ucell[i * 3 + j]                     *= factor[i];
            pNTPData->recipf[i * 3 + j]                     = pNTPData->recip[i * 3 + j];
            pNTPData->ucellf[i * 3 + j]                     = pNTPData->ucell[i * 3 + j];
        }
    }
   
        
    // Recalculate a, b, c
    double a                                                = pNTPData->ucell[0];
    double b                                                = sqrt(pNTPData->ucell[1] * pNTPData->ucell[1] +
                                                                   pNTPData->ucell[4] * pNTPData->ucell[4]);
    double c                                                = sqrt(pNTPData->ucell[2] * pNTPData->ucell[2] + 
                                                                   pNTPData->ucell[5] * pNTPData->ucell[5] + 
                                                                   pNTPData->ucell[8] * pNTPData->ucell[8]);
                                                                   
  //printf("MOL: %d %d %d\n", gpu->sim.soluteMolecules, gpu->maxSoluteMolecules, gpu->maxPSSoluteMolecules);
  // printf("U1: %20.10f %20.10f %20.10f\n", pNTPData->ucell[0], pNTPData->ucell[1], pNTPData->ucell[2]);
  // printf("U2: %20.10f %20.10f %20.10f\n", pNTPData->ucell[3], pNTPData->ucell[4], pNTPData->ucell[5]);
  // printf("U3: %20.10f %20.10f %20.10f\n", pNTPData->ucell[6], pNTPData->ucell[7], pNTPData->ucell[8]);
  // printf("C1: %20.10f %20.10f %20.10f\n", a, b, c);
  // printf("C2: %20.10f %20.10f %20.10f\n", gpu->sim.a, gpu->sim.b, gpu->sim.c);        
        
    // Recalculate nonbond skin
    PMEFloat skin                                           = (a / gpu->sim.xcells - gpu->sim.cut) / gpu->sim.cut_factor[0];
    //printf("%f\n", skin);
    PMEFloat yskin                                          = (b / gpu->sim.ycells - gpu->sim.cut) / gpu->sim.cut_factor[1];
    if (yskin < skin)
        skin                                                = yskin;
  //  printf("%f\n", skin);    
    PMEFloat zskin                                          = (c / gpu->sim.zcells - gpu->sim.cut) / gpu->sim.cut_factor[2];
    if (zskin < skin)
        skin                                                = zskin;
  //  printf("%f\n", skin);    
    
  //  printf("C4: %20.10f %20.10f\n", gpu->sim.nonbond_skin, skin); 
  //  printf("C5: %20.10f %20.10f %20.10f\n", gpu->sim.cut_factor[0], gpu->sim.cut_factor[1], gpu->sim.cut_factor[2]);
    skin *= 0.99;
    pNTPData->one_half_nonbond_skin_squared                 = 0.25 * skin * skin;
    pNTPData->cutPlusSkin2                                  = (gpu->sim.cut + skin) * (gpu->sim.cut + skin);
 //   printf("Skin %f\n", pNTPData->one_half_nonbond_skin_squared);
    gpu->pbNTPData->Upload();    

    kCalculateSoluteCOM(gpu);
    kPressureScaleCoordinates(gpu);
    if (gpu->sim.constraints > 0)
    {
        kCalculateSoluteConstraintsCOM(gpu);
        kReduceSoluteConstraintsCOM(gpu);
        kPressureScaleConstraintCoordinates(gpu);
    }  
}

extern "C" void gpu_final_gb_setup_(unsigned int* igb, unsigned int* alpb, double* saltcon, double* rgbmax, double *gb_neckcut, double* gb_fs_max, double atm_gb_radii[], double atm_gb_fs[], double neckMaxVal[][21], double neckMaxPos[][21])
{
PRINTMETHOD("gpu_final_gb_setup");
    static const PMEDouble ta                   = 1.0 / 3.0;
    static const PMEDouble tb                   = 2.0 / 5.0;
    static const PMEDouble tc                   = 3.0 / 7.0;
    static const PMEDouble td                   = 4.0 / 9.0;
    static const PMEDouble tdd                  = 5.0 / 11.0;
    static const PMEDouble te                   = 4.0 / 3.0;
    static const PMEDouble tf                   = 12.0 / 5.0;
    static const PMEDouble tg                   = 24.0 / 7.0;
    static const PMEDouble th                   = 40.0 / 9.0;
    static const PMEDouble thh                  = 60.0 / 11.0;
    static const PMEDouble alpb_alpha           = 0.571412;


    // Count atom types
    unsigned int ntypes                         = 1;
    PMEDouble* gb_radii                         = new PMEDouble[gpu->sim.atoms];
    PMEDouble* gb_fs                            = new PMEDouble[gpu->sim.atoms];
    unsigned int* gb_type                       = new unsigned int[gpu->sim.atoms];
    gpu->sim.igb                                = *igb;
    gpu->sim.alpb                               = *alpb;
    gpu->sim.saltcon                            = *saltcon;
    gpu->sim.rgbmax                             = *rgbmax;
    
    // Set up simulation parameters
    gpu->sim.gb_alpha                           = 0.0;
    gpu->sim.gb_beta                            = 0.0;
    gpu->sim.gb_gamma                           = 0.0;
    gpu->sim.gb_fs_max                          = 0.0;
    gpu->sim.gb_kappa                           = 0.0;
    gpu->sim.gb_kappa_inv                       = 0.0;
    gpu->sim.extdiel_inv                        = 0.0;
    gpu->sim.intdiel_inv                        = 0.0;
    gpu->sim.offset                             = 0.09;
    if (gpu->sim.igb == 1)
    {
        gpu->sim.gb_alpha                       = 1.0;
        gpu->sim.gb_beta                        = 0.0;
        gpu->sim.gb_gamma                       = 0.0;
    }
    if (gpu->sim.igb == 2)
    {
        // Use our best guesses for Onufriev/Case GB  (GB^OBC I):
        gpu->sim.gb_alpha                       = 0.8;
        gpu->sim.gb_beta                        = 0.0;
        gpu->sim.gb_gamma                       = 2.909125;
    }
    else if (gpu->sim.igb == 5)
    {
        // Use our second best guesses for Onufriev/Case GB (GB^OBC II):
        gpu->sim.gb_alpha                       = 1.0;
        gpu->sim.gb_beta                        = 0.80;
        gpu->sim.gb_gamma                       = 4.8510;
    }
    else if (gpu->sim.igb == 7)
    {
        gpu->sim.gb_alpha                       = 1.09511284;
        gpu->sim.gb_beta                        = 1.90792938;
        gpu->sim.gb_gamma                       = 2.50798245;
        gpu->sim.gb_neckscale                   = 0.361825;
        gpu->sim.gb_neckoffset                  = 1.0 - gpu->sim.offset;        
    }
    else if (gpu->sim.igb == 8)
    {
        gpu->sim.gb_neckscale                   = 0.826836;
        gpu->sim.offset                         = 0.195141;
        gpu->sim.gb_neckoffset                  = 1.0 - gpu->sim.offset;
    }
    gpu->sim.gb_neckcut                         = *gb_neckcut + 2.0 * gpu->sim.offset;
    
    
    if ((gpu->sim.igb == 7) || (gpu->sim.igb == 8))
    {
        delete gpu->pbNeckMaxValPos;
        gpu->pbNeckMaxValPos                    = new GpuBuffer<PMEFloat2>(441);
        for (int i = 0; i < 21; i++)
        {
            for (int j = 0; j < 21; j++)
            {
                gpu->pbNeckMaxValPos->_pSysData[i * 21 + j].x
                                                = neckMaxVal[j][i];
                gpu->pbNeckMaxValPos->_pSysData[i * 21 + j].y
                                                = neckMaxPos[j][i];
            }
        }
        gpu->pbNeckMaxValPos->Upload();
        gpu->sim.pNeckMaxValPos                 = gpu->pbNeckMaxValPos->_pDevData;
    }
    
    if (gpu->sim.saltcon >= 0.0)
    {
        gpu->sim.gb_kappa = 0.73 * sqrt(0.108060 * gpu->sim.saltcon);
    }
        
    if (gpu->sim.alpb == 0)
    {
        // Standard Still's GB
        gpu->sim.extdiel_inv                    = 1.0 / gpu->sim.extdiel;
        gpu->sim.intdiel_inv                    = 1.0 / gpu->sim.intdiel;
    }
    else
    {
        // Sigalov Onufriev ALPB (epsilon-dependent GB):
        gpu->sim.alpb_beta                      = alpb_alpha * (gpu->sim.intdiel / gpu->sim.extdiel);
        gpu->sim.extdiel_inv                    = 1.0 / (gpu->sim.extdiel * (1.0 + gpu->sim.alpb_beta));
        gpu->sim.intdiel_inv                    = 1.0 / (gpu->sim.intdiel * (1.0 + gpu->sim.alpb_beta));
        gpu->sim.one_arad_beta                  = gpu->sim.alpb_beta / gpu->sim.arad;
        if (gpu->sim.gb_kappa != 0.0) 
            gpu->sim.gb_kappa_inv               = 1.0 / gpu->sim.gb_kappa;
    }        
        
        
    
    
    gb_radii[0]                                 = atm_gb_radii[0];
    gb_fs[0]                                    = atm_gb_fs[0];
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        int j;
        for (j = 0; j < ntypes; j++)
        {
            if (((PMEDouble)atm_gb_radii[i] == gb_radii[j]) && ((PMEDouble)atm_gb_fs[i] == gb_fs[j]))
                break;
        }
        gb_type[i]                              = j;
        
        // Add new type if discovered
        if (j == ntypes)
        {
            gb_radii[j]                         = atm_gb_radii[i];
            gb_fs[j]                            = atm_gb_fs[i];
            if (gpu->sim.gb_fs_max < gb_fs[j])
                gpu->sim.gb_fs_max              = gb_fs[j];
            ntypes++;
        }
    }
#ifdef GVERBOSE
    printf("%d GB atom types encountered.\n", ntypes);
    for (int j = 0; j < ntypes; j++)
        printf("%5d %16.8f %16.8f\n", j, gb_radii[j], gb_fs[j]);
#endif
    gpu->sim.rgbmax1i                           = 1.0 / gpu->sim.rgbmax;
    gpu->sim.rgbmax2i                           = gpu->sim.rgbmax1i * gpu->sim.rgbmax1i;
    gpu->sim.rgbmaxpsmax2                       = (gpu->sim.rgbmax + gpu->sim.gb_fs_max) * (gpu->sim.rgbmax + gpu->sim.gb_fs_max);
    gpuCopyConstants();
        
        
    // Allocate and generate look-up table
#if 0
    GpuBuffer<PMEFloat>* psGBBRTexture          = new GpuBuffer<PMEFloat>(ntypes * ntypes * GB_TEXTURE_WIDTH);
    PMEFloat* pBRData                           = psGBBRTexture->_pSysData;
    GpuBuffer<PMEFloat>* psGBNBTexture          = new GpuBuffer<PMEFloat>(ntypes * ntypes * GB_TEXTURE_WIDTH);
    PMEFloat* pNBData                           = psGBNBTexture->_pSysData;
    for (int i = 0; i < ntypes; i++)
    {
        double ri                               = gb_radii[i] - gpu->sim.offset;
        double ri1i                             = 1.0 / ri;
        double si                               = gb_fs[i];
        double si2                              = si * si;
        int ni                                  = (int)((gb_radii[i] - 1.0) * 20.0 + 0.5);
        for (int j = 0; j < ntypes; j++)
        {
            double rj                           = gb_radii[j] - gpu->sim.offset;
            double rj1i                         = 1.0 / rj;
            double sj                           = gb_fs[j];
            double sj2                          = sj * sj;
            int nj                              = (int)((gb_radii[j] - 1.0) * 20.0 + 0.5);
            for (unsigned int k = 0; k < GB_TEXTURE_WIDTH; k++)
            {
                double dij                      = k * (gpu->sim.rgbmax + gpu->sim.gb_fs_max) / GB_TEXTURE_WIDTH;
                if (dij < 0.00001)      // Avoid divide by zero
                    dij = 0.00001;
                double r2                       = dij * dij;
                double dij1i                    = 1.0 / dij;
                double dij2i                    = dij1i * dij1i;
                double dij3i                    = dij2i * dij1i;
                if (dij > gpu->sim.rgbmax + sj)
                {
                    pBRData[k]                  = 0.0;
                    pNBData[k]                  = 0.0;
                }
                else if (dij > gpu->sim.rgbmax - sj)
                {
                    double uij                  = 1.0 / (dij - sj);
                    pBRData[k]                  = 0.125 * dij1i * (1.0 + 2.0 * dij * uij + 
                                                  gpu->sim.rgbmax2i * (r2 - 4.0 * gpu->sim.rgbmax * dij - sj2) + 
                                                  2.0 * log((dij - sj) * gpu->sim.rgbmax1i));
                    double temp1                = 1.0 / (dij - sj);
                    pNBData[k]                  = 0.125 * dij3i * ((r2 + sj2) * 
                                                  (temp1 * temp1 - gpu->sim.rgbmax2i) - 2.0 * log(gpu->sim.rgbmax * temp1));
                }
                else if (dij > 4.0 * sj)
                {            
                    double dij2i                = dij1i * dij1i;
                    double tmpsd                = sj2 * dij2i;
                    double dumbo                = ta + tmpsd *  (tb + tmpsd * (tc + tmpsd * (td + tmpsd * tdd)));
                    pBRData[k]                  = tmpsd * sj * dij2i * dumbo;
                    dumbo                       = te + tmpsd * (tf + tmpsd * (tg + tmpsd * (th + tmpsd * thh)));
                    pNBData[k]                  = tmpsd * sj * dij2i * dij2i * dumbo;
                }
                else if (dij > ri + sj)
                {
                    double v2                   = 1.0 / (dij + sj);
                    double v4                   = log(v2 * (dij - sj));
                    pBRData[k]                  = 0.5 * (sj / (r2 - sj2) + 0.5 * dij1i * v4);
                    v2                          = 1.0 / (r2 - sj2);
                    v4                          = log((dij - sj) / (dij + sj));
                    pNBData[k]                  = v2 * sj * (-0.5 * dij2i + v2) +
                                                  0.25 * dij3i * v4;
                }
                else if (dij > fabs(ri - sj))
                {
                    double v2                   = 1.0 / (dij + sj);
                    double v4                   = log(v2 * ri);
                    double theta                = 0.5 * ri1i * dij1i * (r2 + ri * ri - sj2);
                    pBRData[k]                  = 0.25 * (ri1i * (2.0 - theta) - 
                                                  v2 + dij1i * v4);
                    v4                          = log(ri / (dij + sj));
                    pNBData[k]                  = -0.25 * (-0.5 * (r2 - ri * ri + sj2) *
                                                  dij3i * ri1i * ri1i + dij1i * v2 *
                                                  (v2 - dij1i) - dij3i * v4);                   
                }        
                else if (ri < sj)
                {
                    double v2                   = 1.0 / (dij + sj);
                    double v4                   = log(v2 * (sj - dij));
                    pBRData[k]                  = 0.5 * (sj / (r2 - sj2) + 2.0 * ri1i + 
                                                  0.5 * dij1i * v4);
                    v2                          = 1.0 / (r2 - sj2);
                    v4                          = log ((sj - dij) / (dij + sj));
                    pNBData[k]                  = -0.5 * (sj * dij2i * v2 -
                                                  2.0 * sj * v2 * v2 - 
                                                  0.5 * dij3i * v4);
                }
                else
                {
                    pBRData[k]                  = 0.0;   
                    pNBData[k]                  = 0.0;        
                }
                // printf("%4d %16.8f %16.8f\n", k, pBRData[k], pNBData[k]); 
            }      
            pBRData                                += GB_TEXTURE_WIDTH;
        }
   
    }
#endif           
    delete[] gb_radii;
    delete[] gb_fs;
    delete[] gb_type;

}



extern "C" void gpu_gb_igb8_setup_(double gb_alpha[], double gb_beta[], double gb_gamma[])
{
PRINTMETHOD("gpu_gb_igb8_setup");

    // Delete existing alpha, beta, and gamma arrays
    delete gpu->pbGBAlphaBetaGamma;
    
    // Allocate new arrays
    gpu->pbGBAlphaBetaGamma                     = new GpuBuffer<PMEFloat>(gpu->sim.stride3);
    gpu->sim.pgb_alpha                          = gpu->pbGBAlphaBetaGamma->_pDevData;
    gpu->sim.pgb_beta                           = gpu->sim.pgb_alpha + gpu->sim.stride;
    gpu->sim.pgb_gamma                          = gpu->sim.pgb_alpha + gpu->sim.stride2;
    
    // Fill arrays
    PMEFloat* pAlpha                            = gpu->pbGBAlphaBetaGamma->_pSysData;
    PMEFloat* pBeta                             = gpu->pbGBAlphaBetaGamma->_pSysData + gpu->sim.stride;
    PMEFloat* pGamma                            = gpu->pbGBAlphaBetaGamma->_pSysData + gpu->sim.stride2;
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        pAlpha[i]                               = gb_alpha[i];
        pBeta[i]                                = gb_beta[i];
        pGamma[i]                               = gb_gamma[i];
    }
    gpu->pbGBAlphaBetaGamma->Upload();
    gpuCopyConstants();
}


extern "C" void gpu_build_threadblock_work_list_(int numex[], int natex[])
{
PRINTMETHOD("gpu_build_threadblock_work_list");
    const unsigned int dim                      = (gpu->sim.paddedNumberOfAtoms + (GRID - 1)) / GRID;
    const unsigned int cells                    = dim * (dim + 1) / 2;
  
    // Delete existing data
    delete gpu->pbWorkUnit;
    delete gpu->pbExclusion;
    delete gpu->pbGBPosition;
#ifdef MPI    
    delete[] gpu->pMinLocalAtom;
    delete[] gpu->pMaxLocalAtom;
    delete[] gpu->pAllGathervRecvCountAoS;
    delete[] gpu->pAllGathervRecvDisplAoS;
    delete[] gpu->pAllGathervRecvCountSoA;
    delete[] gpu->pAllGathervRecvDisplSoA;
#endif
    
    // Determine number of work units with exclusions so we can place them at the head of the list
    
    // Generate symmetric exclusion list based on assymmetric -1 padded list
    // which is capped at 2x the previous list size, but is likely smaller due to
    // superfluous -1 entries
    unsigned int totalExclusions                = 0; 
    unsigned int* pExclCount                    = new unsigned int[gpu->sim.paddedNumberOfAtoms + 1];
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        totalExclusions                        += numex[i];
        pExclCount[i]                           = 0;
    }
    for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++)
        pExclCount[i]                           = 0;
       
        

    unsigned int* pExclList                     = new unsigned int[totalExclusions * 2];
    unsigned int* pExclOffset                   = new unsigned int[gpu->sim.paddedNumberOfAtoms + 1];
    // Count Exclusions and filter out -1 entries

    unsigned int offset                         = 0;
    for (unsigned int x = 0; x < gpu->sim.atoms; x++)
    {
        for (unsigned int j = offset; j < offset + numex[x]; j++)
        {
            if (natex[j] > 0)
            {
                unsigned int y                  = natex[j] - 1;
                pExclCount[x]++;
                pExclCount[y]++;
            }
        }
        offset                                 += numex[x]; 
    } 
    
    // Calculate symmetric exclusion offsets
    unsigned int* pExclCounter                  = new unsigned int[gpu->sim.paddedNumberOfAtoms + 1];
    offset                                      = 0;
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        pExclOffset[i]                          = offset;
        pExclCounter[i]                         = offset;
        offset                                 += pExclCount[i];
    }
    for (int i = gpu->sim.atoms; i <= gpu->sim.paddedNumberOfAtoms; i++)
    {
        pExclOffset[i]                          = offset;
        pExclCounter[i]                         = offset;
    }
    
    // Now regenerate exclusions
    offset                                      = 0;
    for (int x = 0; x < gpu->sim.atoms; x++)
    {
        for (int j = offset; j < offset + numex[x]; j++)
        {
            if (natex[j] > 0)
            {
                unsigned int y                  = natex[j] - 1;
                pExclList[pExclCounter[x]++]    = y;
                pExclList[pExclCounter[y]++]    = x;
            }
        }
        offset                                 += numex[x];
    }
    
#if 0
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        printf("%05d: %05d ", i, pExclCount[i]);
     //   for (int j = pExclOffset[i]; j < pExclOffset[i + 1]; j++)
     //       printf("%5d ", pExclList[j]);
        printf("\n");
    }
    exit(-1);
#endif
    
    // Partition nonbond calculation amongst nodes
#ifdef MPI
    unsigned int cellCounter                    = 0;
    unsigned int minDim                         = (gpu->gpuID * dim) / gpu->nGpus;
    unsigned int maxDim                         = ((gpu->gpuID + 1) * dim) / gpu->nGpus;
    gpu->sim.minLocalAtom                       = minDim * GRID;
    gpu->sim.maxLocalAtom                       = min(maxDim * GRID, (unsigned int)gpu->sim.atoms);

    {
        gpu->pMinLocalAtom                      = new int[gpu->nGpus];
        gpu->pMaxLocalAtom                      = new int[gpu->nGpus];
        gpu->pAllGathervRecvCountAoS            = new int[gpu->nGpus];
        gpu->pAllGathervRecvDisplAoS            = new int[gpu->nGpus];
        gpu->pAllGathervRecvCountSoA            = new int[gpu->nGpus];
        gpu->pAllGathervRecvDisplSoA            = new int[gpu->nGpus];

        
        for (int i = 0; i < gpu->nGpus; i++)
        {
            unsigned int minDim                 = (i * dim) / gpu->nGpus;
            unsigned int maxDim                 = ((i + 1) * dim) / gpu->nGpus;
            gpu->pMinLocalAtom[i]               = minDim * GRID;
            gpu->pMaxLocalAtom[i]               = min(maxDim * GRID, (unsigned int)gpu->sim.atoms);    
            gpu->pAllGathervRecvCountAoS[i]     = (gpu->pMaxLocalAtom[i] - gpu->pMinLocalAtom[i]) * 3;
            gpu->pAllGathervRecvDisplAoS[i]     = gpu->pMinLocalAtom[i] * 3;    
            gpu->pAllGathervRecvCountSoA[i]     = gpu->pMaxLocalAtom[i] - gpu->pMinLocalAtom[i];
            gpu->pAllGathervRecvDisplSoA[i]     = gpu->pMinLocalAtom[i];    
        }
    }
    gpu->sim.minReducedAtom                     = GRID * ((gpu->gpuID * dim) / gpu->nGpus);
    gpu->sim.maxReducedAtom                     = GRID * (((gpu->gpuID + 1) * dim) / gpu->nGpus);
    gpu->sim.reducedAtoms                       = gpu->sim.maxReducedAtom - gpu->sim.minReducedAtom;
    gpu->sim.reducedAtoms3                      = gpu->sim.reducedAtoms * 3;
#endif
 
    unsigned int excludedWorkUnits              = 0;   
    for (unsigned int y = 0; y < dim; y++)
    {
        for (unsigned int x = y; x < dim; x++)
        {
#ifdef MPI
            if (((x >= minDim) && (x < maxDim)) || ((y >= minDim) && (y < maxDim))) 
            {
#endif
                unsigned int xstart             = x * GRID;
                unsigned int ystart             = y * GRID;
                unsigned int xend               = xstart + GRID;   
                unsigned int yend               = ystart + GRID;
                bool excluded                   = false;
                for (int i = pExclOffset[xstart]; i < pExclOffset[xend]; i++)
                {
                    if ((pExclList[i] >= ystart) && (pExclList[i] < yend))
                    {
                        excluded                = true;
                        break;
                    }
                }
                if (excluded)
                    excludedWorkUnits++;
#ifdef MPI
                cellCounter++;
            }
#endif
        }
        
    }
#ifdef MPI
    gpu->sim.workUnits                          = cellCounter;
#else   
    gpu->sim.workUnits                          = cells;
#endif


    // Decrease thread count for extra small molecules to spread computation
    // across entire chip   
    int balancedWorkBlock                       = GRID * (1 + gpu->sim.workUnits / gpu->blocks);
    int activeWorkUnits                         = (gpu->blocks * gpu->GBBornRadiiThreadsPerBlock) / GRID;
    if (activeWorkUnits > gpu->sim.workUnits + gpu->blocks)
    {
        gpu->GBBornRadiiThreadsPerBlock         = balancedWorkBlock;
    }
    activeWorkUnits                             = (gpu->blocks * gpu->GBNonbondEnergy1ThreadsPerBlock) / GRID;
    if (activeWorkUnits > gpu->sim.workUnits + gpu->blocks)
    {
        gpu->GBNonbondEnergy1ThreadsPerBlock    = balancedWorkBlock;
    }
    activeWorkUnits                             = (gpu->blocks * gpu->GBNonbondEnergy2ThreadsPerBlock) / GRID;
    if (activeWorkUnits > gpu->sim.workUnits + gpu->blocks)
    {
        gpu->GBNonbondEnergy2ThreadsPerBlock    = balancedWorkBlock;
    }    


    // Build work unit list
    gpu->sim.excludedWorkUnits                  = excludedWorkUnits;
    unsigned int excludedOffset                 = 0;
    unsigned int unexcludedOffset               = excludedWorkUnits;
    unsigned int exclusions                     = 0;
     
    GpuBuffer<unsigned int>* pbExclusion        = new GpuBuffer<unsigned int>(GRID * excludedWorkUnits);
    unsigned int* pExclusion                    = pbExclusion->_pSysData;
    GpuBuffer<unsigned int>* pbWorkUnit         = new GpuBuffer<unsigned int>(gpu->sim.workUnits);
    unsigned int* pWorkUnit                     = pbWorkUnit->_pSysData; 
    gpu->pbGBPosition                           = new GpuBuffer<unsigned int>(3);   


    for (unsigned int y = 0; y < dim; y++)
    {
        for (unsigned int x = y; x < dim; x++)
        {       
#ifdef MPI
            if (((x >= minDim) && (x < maxDim)) || ((y >= minDim) && (y < maxDim))) 
            {
#endif
                // Check for exclusions
                unsigned int xstart             = x * GRID;
                unsigned int ystart             = y * GRID;
                unsigned int xend               = xstart + GRID;   
                unsigned int yend               = ystart + GRID;
                bool excluded                   = false;
                for (int i = pExclOffset[xstart]; i < pExclOffset[xend]; i++)
                {
                    if ((pExclList[i] >= ystart) && (pExclList[i] < yend))
                    {
                        excluded                = true;
                        break;
                    }
                }

                // Add exclusions if present
                if (excluded)
                {
            
                    // Create exclusion masks
                    unsigned int excl[GRID];
                    for (int i = 0; i < GRID; i++)
                        excl[i]                 = 0xffffffff;
                    
                    for (int i = 0; i < GRID; i++)
                    {
                        unsigned int x          = xstart + i;
                        for (int j = pExclOffset[x]; j < pExclOffset[x + 1]; j++)
                        {
                            unsigned int y = pExclList[j];
                            if ((y >= ystart) && (y < yend))
                            {
                                excl[i] ^= 1 << (y - ystart);
                            }
                        }
                    }
                
                    // Skip padded atoms
                    if (xend > gpu->sim.atoms)
                    {
                        for (int i = gpu->sim.atoms - xstart; i < GRID; i++)
                            excl[i]             = 0;
                    }
                    if (yend > gpu->sim.atoms)
                    {
                        unsigned int offset = yend - gpu->sim.atoms;
                        for (int i = 0; i < GRID; i++)
                            excl[i]             = (excl[i] << offset) >> offset;
                    }
                
                    // Post-process exclusion masks
                    for (int i = 0; i < GRID; i++)
                    {
                
                        unsigned int offset = i;
                        if (xstart == ystart)
                            offset++;
                        excl[i]                 = (excl[i] >> offset) | (excl[i] << (GRID - offset));
                        pExclusion[exclusions++]= excl[i];
                    }
                    pWorkUnit[excludedOffset]   = (x << 17) | (y << 2);
                    pWorkUnit[excludedOffset]  |= 0x1;
                    excludedOffset++;
                }
                else
                {
                    pWorkUnit[unexcludedOffset] = (x << 17) | (y << 2);
                    unexcludedOffset++;
                }
#ifdef MPI
            }
#endif            
        }
    }    
    
    // Delete temporary data
    delete[] pExclOffset;
    delete[] pExclCount;
    delete[] pExclCounter;
    delete[] pExclList;

    // Set up GPU pointers and constants
    gpu->pbWorkUnit                             = pbWorkUnit;
    gpu->sim.pWorkUnit                          = pbWorkUnit->_pDevData;
    gpu->pbExclusion                            = pbExclusion;
    gpu->sim.pExclusion                         = pbExclusion->_pDevData;
    gpu->sim.pGBBRPosition                      = gpu->pbGBPosition->_pDevData;
    gpu->sim.pGBNB1Position                     = gpu->pbGBPosition->_pDevData + 1;
    gpu->sim.pGBNB2Position                     = gpu->pbGBPosition->_pDevData + 2;
    gpu->sim.GBTotalWarps[1]                    = (gpu->GBNonbondEnergy1ThreadsPerBlock * gpu->blocks) / GRID;
    
    // Special case thread counts for C1060/GT2xx
    if ((gpu->sm_version < SM_2X) && ((gpu->sim.igb == 7) || (gpu->sim.igb == 8)))
    {
        gpu->sim.GBTotalWarps[0]                = (gpu->GBBornRadiiIGB78ThreadsPerBlock * gpu->blocks) / GRID;
        gpu->sim.GBTotalWarps[2]                = (gpu->GBNonbondEnergy2IGB78ThreadsPerBlock * gpu->blocks) / GRID;    
    }
    else 
    {
        gpu->sim.GBTotalWarps[0]                = (gpu->GBBornRadiiThreadsPerBlock * gpu->blocks) / GRID;
        gpu->sim.GBTotalWarps[2]                = (gpu->GBNonbondEnergy2ThreadsPerBlock * gpu->blocks) / GRID;
    }
    pbWorkUnit->Upload();
    pbExclusion->Upload();
    gpuCopyConstants();
    return;
}

extern "C" void gpu_final_setup_()
{
PRINTMETHOD("gpu_final_setup");
    
}


extern "C" void gpu_update_(double* dt, double* temp0, double* gamma_ln)
{
#if 0
    {
        FILE* fp = fopen("../forces.out", "w");
        gpu->pbForce->Download();
        PMEDouble *pCrd                                 = gpu->pbForce->_pSysData;
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, pCrd[i], pCrd[i + gpu->sim.stride], pCrd[i + gpu->sim.stride2]);
        }
        fclose(fp);
    }
#endif

PRINTMETHOD("gpu_update");
#if 0
    {
        gpu->pbImageIndex->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);        
        PMEDouble *pCrd                                 = gpu->pbImage->_pSysData; 
        gpu->pbImage->Download();     
     
        if (gpu->sim.pImageX != gpu->pbImage->_pDevData)
           pCrd                                        = gpu->pbImage->_pSysData + gpu->sim.stride3;

        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            unsigned int j = pImageAtomLookup[i];
            printf("%6d %30.20f %30.20f %30.20f\n", i, pCrd[j], pCrd[j + gpu->sim.stride], pCrd[j + gpu->sim.stride2]);
        }
    }
    exit(-1);
#endif  
    kUpdate(gpu, *dt, *temp0, *gamma_ln);
    
#if 0
    {
        char buff[512];
        PMEDouble *pCrd                                 = gpu->pbImage->_pSysData; 
        gpu->pbImage->Download();     
#ifdef MPI  
        sprintf(buff, "../npos%d.out", gpu->gpuID);
#else
        sprintf(buff, "../pos1.out");
#endif
        
        if (gpu->sim.pImageX != gpu->pbImage->_pDevData)
           pCrd                                        = gpu->pbImage->_pSysData + gpu->sim.stride3;

        FILE* fp = fopen(buff, "w");
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            fprintf(fp, "%6d %30.20f %30.20f %30.20f\n", i, pCrd[i], pCrd[i + gpu->sim.stride], pCrd[i + gpu->sim.stride2]);
        }
        fclose(fp);
    }
#endif
#if 0
    {
        char buff[512];
#ifdef MPI  
        sprintf(buff, "../nopos%d.out", gpu->gpuID);
#else
        sprintf(buff, "../opos1.out");
#endif
        gpu->pbForce->Download();  
        PMEDouble *pCrd                                = gpu->pbForce->_pSysData;    
        FILE* fp = fopen(buff, "w");
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            fprintf(fp, "%6d %16.7f %16.7f %16.7f\n", i, pCrd[i], pCrd[i + gpu->sim.stride], pCrd[i + gpu->sim.stride2]);
        }
        fclose(fp);
    }
#endif
#if 0
    {
        gpu->pbImageIndex->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);        
        PMEDouble *pCrd                                 = gpu->pbImage->_pSysData; 
        gpu->pbImage->Download();     
     
        if (gpu->sim.pImageX != gpu->pbImage->_pDevData)
           pCrd                                        = gpu->pbImage->_pSysData + gpu->sim.stride3;

        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            unsigned int j = pImageAtomLookup[i];
            printf("%6d %30.20f %30.20f %30.20f\n", i, pCrd[j], pCrd[j + gpu->sim.stride], pCrd[j + gpu->sim.stride2]);
        }
    }
    exit(-1);
#endif  

}

extern "C" void gpu_shake_()
{
PRINTMETHOD("gpu_shake");
    kShake(gpu);
    
#if 0
    {
        char buff[512];
#ifdef MPI  
        PMEDouble *pCrd                                 = gpu->pbExchangeData->_pSysData;
        sprintf(buff, "../nspos%d.out", gpu->gpuID);
#else
        sprintf(buff, "../spos1.out");
        gpu->pbImage->Download();  
        PMEDouble *pCrd                                 = gpu->pbImage->_pSysData;    
        if (gpu->sim.pImageX != gpu->pbImage->_pDevData)
           pCrd                                        = gpu->pbImage->_pSysData + gpu->sim.stride3;
#endif
        FILE* fp = fopen(buff, "w");
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            fprintf(fp, "%6d %16.7f %16.7f %16.7f\n", i, pCrd[i], pCrd[i + gpu->sim.stride], pCrd[i + gpu->sim.stride2]);
        }
        fclose(fp);
    }
#endif
}

extern "C" void gpu_vrand_reset_velocities_(double* temp, double* half_dtx)
{

PRINTMETHOD("gpu_vrand_reset_velocities");
    kResetVelocities(gpu, *temp, *half_dtx);   
}


extern "C" void gpu_recalculate_velocities_(double* dtx_inv)
{
PRINTMETHOD("gpu_recalculate_velocities");
    kRecalculateVelocities(gpu, *dtx_inv);
    
#if 0
    {
        gpu->pbImageIndex->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);        
        PMEDouble *pVel                                 = gpu->pbImageVel->_pSysData; 
        gpu->pbImageVel->Download();     
     
        if (gpu->sim.pImageVelX != gpu->pbImageVel->_pDevData)
           pVel                                         = gpu->pbImageVel->_pSysData + gpu->sim.stride3;

        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            unsigned int j = pImageAtomLookup[i];
            printf("%6d,%32.15f,%32.15f,%32.15f\n", i, pVel[j], pVel[j + gpu->sim.stride], pVel[j + gpu->sim.stride2]);
        }
    }
//    exit(-1);
#endif        
}

extern "C" void gpu_scale_velocities_(double* scale)
{
PRINTMETHOD("gpu_scale_velocities");
    kScaleVelocities(gpu, *scale);
}

extern "C" void gpu_calculate_kinetic_energy_(double* c_ave, double* eke, double* ekph, double* ekpbs)
{
PRINTMETHOD("gpu_calculate_kinetic_energy");
    kCalculateKineticEnergy(gpu, *c_ave);
    
    PMEFloat E[3];
    E[0]                                        = 0.0f;     // EKE
    E[1]                                        = 0.0f;     // EKPH
    E[2]                                        = 0.0f;     // EKPBS
    for (int i = 0; i < gpu->blocks; i++)
    {
        E[0]                                   += gpu->pbKineticEnergyBuffer->_pSysData[i].KE.EKE;
        E[1]                                   += gpu->pbKineticEnergyBuffer->_pSysData[i].KE.EKPH;
        E[2]                                   += gpu->pbKineticEnergyBuffer->_pSysData[i].KE.EKPBS;
    }

    *eke                                        = E[0];
    *ekph                                       = E[1];
    *ekpbs                                      = E[2];


#if 0
    {
	// Calculate bounding box to test Minimum Image Convention
	gpu->pbImage->Download();
	PMEDouble xmin = 99999999.0;
	PMEDouble ymin = 99999999.0;
    PMEDouble zmin = 99999999.0;
	PMEDouble xmax = -99999999.0;
	PMEDouble ymax = -99999999.0;
    PMEDouble zmax = -99999999.0;
	PMEDouble* pImageX = gpu->pbImage->_pSysData;
	if (gpu->sim.pImageX != gpu->pbImage->_pDevData)
		pImageX += gpu->sim.stride3;
	PMEDouble* pImageY = pImageX + gpu->sim.stride;
	PMEDouble* pImageZ = pImageX + gpu->sim.stride2;

	for (int i = 0; i < gpu->sim.atoms; i++)
	{
		if (pImageX[i] < xmin)
			xmin = pImageX[i];
		if (pImageX[i] > xmax)
			xmax = pImageX[i];
		if (pImageY[i] < ymin)
			ymin = pImageY[i];
		if (pImageY[i] > ymax)
			ymax = pImageY[i];
		if (pImageZ[i] < zmin)
			zmin = pImageZ[i];
		if (pImageZ[i] > zmax)
			zmax = pImageZ[i];
	}
	printf("Dbounds: %16.8f %16.8f %16.8f\n", xmax - xmin, ymax - ymin, zmax - zmin);
	fflush(stdout);
    }
#endif


#if 0
    {
    // Calculate bounding box to test Minimum Image Convention
	gpu->pbAtomXYSP->Download();
    gpu->pbAtomZSP->Download();
	PMEFloat xmin = 99999999.0;
	PMEFloat ymin = 99999999.0;
    PMEFloat zmin = 99999999.0;
	PMEFloat xmax = -99999999.0;
	PMEFloat ymax = -99999999.0;
    PMEFloat zmax = -99999999.0;
	PMEFloat2* pAtomXY = gpu->pbAtomXYSP->_pSysData;
	PMEFloat* pAtomZ = gpu->pbAtomZSP->_pSysData;

	for (int i = 0; i < gpu->sim.atoms; i++)
	{
		if (pAtomXY[i].x < xmin)
			xmin = pAtomXY[i].x;
		if (pAtomXY[i].x > xmax)
			xmax = pAtomXY[i].x;
		if (pAtomXY[i].y < ymin)
			ymin = pAtomXY[i].y;
		if (pAtomXY[i].y > ymax)
			ymax = pAtomXY[i].y;
		if (pAtomZ[i] < zmin)
			zmin = pAtomZ[i];
		if (pAtomZ[i] > zmax)
			zmax = pAtomZ[i];
	}
	printf("Sbounds: %16.8f %16.8f %16.8f\n", xmax - xmin, ymax - ymin, zmax - zmin);
	fflush(stdout);
    }
#endif

}

extern "C" void gpu_recenter_molecule_()
{
PRINTMETHOD("gpu_recenter_molecule");
    kRecenter_Molecule(gpu);
}

extern "C" void gpu_amrset_(int* seed)
{
PRINTMETHOD("gpu_amrset");
#ifdef CPU_RANDOMS
    cpu_amrset(*seed);
    cpu_kRandom(gpu);
#else
    curandCreateGenerator(&gpu->RNG, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gpu->RNG, *seed);
    kRandom(gpu);
#endif
    gpu->randomCounter = 0;
}

extern "C" void gpu_dump_float_(float* pFloat)
{
    FILE *fp = fopen("float.txt", "w");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        fprintf(fp, "%32.15f\n", pFloat[i]);
        printf("%32.15f\n", pFloat[i]);
    }
    fclose(fp);
    exit(-1);
    
}

extern "C" void gpu_dump_grid_weights_(int* atoms, int* map, double* theta)
{
    for (int i = 1; i <= *atoms; i++)
    {
        for (int j = 0; j < *atoms; j++)
        {
            if (map[j] == i)
            {
                printf("%5d: %13.7f %13.7f %13.7f %13.7f\n", i - 1, theta[j * 12],     theta[j * 12 + 1], theta[j * 12 + 2], theta[j * 12 + 3]);
                printf("%5d: %13.7f %13.7f %13.7f %13.7f\n", i - 1, theta[j * 12 + 4], theta[j * 12 + 5], theta[j * 12 + 6], theta[j * 12 + 7]);
                printf("%5d: %13.7f %13.7f %13.7f %13.7f\n", i - 1, theta[j * 12 + 8], theta[j * 12 + 9], theta[j * 12 + 10], theta[j * 12 + 11]);
            }
        }
    }
    exit(-1);
}

extern "C" void gpu_dump_img_double_vector_(int* count, int* img, double* pDouble)
{
    for (int i = 1; i <= *count; i++)
    {
        for (int j = 0; j < *count; j++)
        {
            if (img[j] == i)
            {
                printf("%5d: %13.7f %13.7f %13.7f\n", i - 1, pDouble[j * 3], pDouble[j * 3 + 1], pDouble[j * 3 + 2]);
            }
        }
    }



    exit(-1);
}

extern "C" void gpu_dump_img_int_vector_(int* img, int* pInt)
{
    for (int i = 1; i <= gpu->sim.atoms; i++)
    {
        for (int j = 0; j < gpu->sim.atoms; j++)
        {
            if (img[j] == i)
            {
                printf("%5d: %9d %9d %9d\n", i - 1, pInt[j * 3], pInt[j * 3 + 1], pInt[j * 3 + 2]);
            }
        }
    }



}

extern "C" void gpu_dump_complex_grid_(double zxy_qt[], int *fft_x_dim, int *fft_y_dim, int *fft_z_dim)
{
    for (int i = 0; i < *fft_x_dim; i++)
    {
        for (int j = 0; j < *fft_y_dim; j++)
        {
            for (int k = 0; k < *fft_z_dim; k++)
            {
                printf("%3d %3d %3d: %40.33f %40.33f\n", i, j, k, 
                    zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 *(*fft_z_dim) + 2 * k],
                    zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 *(*fft_z_dim) + 2 * k + 1]); 
            }
        }
    }
    exit(-1);
}

extern "C" void gpu_load_complex_grid_(double zxy_qt[], int *fft_x_dim, int *fft_y_dim, int *fft_z_dim)
{
    gpu->pbXYZ_qt->Download();
    double maxerror = 0.0;
    for (int i = 0; i < *fft_x_dim; i++)
    {
        for (int j = 0; j < *fft_y_dim; j++)
        {
            for (int k = 0; k < *fft_z_dim; k++)
            {
                double error;
                error = fabs((zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 * (*fft_z_dim) + 2 * k] - 
                    gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].x) / 
                    gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].x);
                if (error > maxerror)
                {
                    printf("%3d %3d %3d: %40.30f %40.30f\n", i, j, k, 
                    zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 *(*fft_z_dim) + 2 * k],
                    gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].x);
                    maxerror = error;
                }
                error = fabs((zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 * (*fft_z_dim) + 2 * k + 1] - 
                    gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].y) / 
                    gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].y);
                if (error > maxerror)
                {
                    maxerror = error;  
                    printf("%3d %3d %3d: %40.30f %40.30f\n", i, j, k, 
                    zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 *(*fft_z_dim) + 2 * k + 1],
                    gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].y);
                    maxerror = error;                    
                }  
#if 0
                printf("%3d %3d %3d: %20.13f %20.13f\n", i, j, k, 
                    zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 *(*fft_z_dim) + 2 * k] - 
                    gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].x ,
                    zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 *(*fft_z_dim) + 2 * k + 1] -
                    gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].y ); 
#endif
#if 0                    
                gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].x  = 
                zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 *(*fft_z_dim) + 2 * k];
                gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].y  = 
                zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 *(*fft_z_dim) + 2 * k + 1];
#else 
                zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 *(*fft_z_dim) + 2 * k];
                gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].x; 
                zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 *(*fft_z_dim) + 2 * k + 1] =
                gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].y;
                
#endif                
            }
        }
    }
    printf("%40.20lf\n", maxerror);
  //  gpu->pbXYZ_qt->Upload();
     exit(-1);
}


extern "C" void gpu_dump_grid_(double xyz_q[], int *fft_x_dim, int *fft_y_dim, int *fft_z_dim)
{
    double dmin = 99999.0;
    double dmax = -99999.0;
    int xstride = 2 * (*fft_x_dim / 2 + 1);
    for (int i = 0; i < *fft_x_dim; i++)
    {
        for (int j = 0; j < *fft_y_dim; j++)
        {
            for (int k = 0; k < *fft_z_dim; k++)
            {
                printf("%3d %3d %3d %32.15f\n", i, j, k, xyz_q[(k * *fft_y_dim + j) * xstride + i]); 
                double v = xyz_q[(k * *fft_y_dim + j) * xstride + i];
                if (v > dmax)
                    dmax = v;
                if (v < dmin)
                    dmin = v;
            }
        }
    }
    printf("%f %f\n", dmin, dmax);
    exit(-1);
}

extern "C" void gpu_dump_int_vector_(int* pInt, int* atoms)
{
#if 1
    FILE *fp = fopen("int.txt", "w");
    for (int i = 0; i < *atoms; i++)
    {
        fprintf(fp, "%5d: %9d %9d %9d\n", i, pInt[i * 3], pInt[i * 3 + 1], pInt[i * 3 + 2]);
        printf("%5d: %9d %9d %9d\n", i, pInt[i * 3], pInt[i * 3 + 1], pInt[i * 3 + 2]);
    }
    fclose(fp);
    exit(-1);
#endif
}

extern "C" void gpu_dump_mapped_int_vector_(int* pInt, int* pMap, int* atoms)
{
    FILE *fp = fopen("int.txt", "w");
    for (int i = 1; i <= *atoms; i++)
    {
        for (int j = 0; j < *atoms; j++)
        {
            if (pMap[j] == i)
            {
                fprintf(fp, "%5d: %9d %9d %9d\n", i - 1, pInt[j * 3], pInt[j * 3 + 1], pInt[j * 3 + 2]);
                printf("%5d: %9d %9d %9d\n", i - 1, pInt[j * 3], pInt[j * 3 + 1], pInt[j * 3 + 2]);
            }
        }
    }
    fclose(fp);
    exit(-1);
}

extern "C" void gpu_dump_mapped_double_vector_(double* pDouble, int* pMap, int* atoms)
{
    for (int i = 1; i <= *atoms; i++)
    {
        for (int j = 0; j < *atoms; j++)
        {
            if (pMap[j] == i)
            {
                printf("%5d: %20.15f %20.15f %20.15f\n", i - 1, pDouble[j * 3], pDouble[j * 3 + 1], pDouble[j * 3 + 2]);
            }
        }
    }
    exit(-1);
}

extern "C" void gpu_dump_int_(int* pInt)
{
    FILE *fp = fopen("int.txt", "w");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        fprintf(fp, "%5d: %9d\n", i, pInt[i]);
        printf("%5d: %9d\n", i, pInt[i]);
    }
    fclose(fp);
    exit(-1);
}


extern "C" void gpu_dump_double_(double* pDouble)
{
    FILE *fp = fopen("double.txt", "w");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        fprintf(fp, "%32.15f\n", pDouble[i]);
        printf("%32.15f\n", pDouble[i]);
    }
    fclose(fp);
    exit(-1);
    
}

extern "C" void gpu_dump_double_vector_(int* atoms, double* pDouble)
{
 //   FILE *fp = fopen("double.txt", "w");
    for (int i = 0; i < *atoms; i++)
    {
   //     fprintf(fp, "%5d: %32.15f %32.15f %32.15f\n", i, pDouble[i * 3], pDouble[i * 3 + 1], pDouble[i * 3 + 2]);
        printf("%6d,%32.15f,%32.15f,%32.15f\n", i, pDouble[i * 3], pDouble[i * 3 + 1], pDouble[i * 3 + 2]);
    }
 //   fclose(fp);
 //   exit(-1);
    
}



extern "C" void gpu_dump_dval_(double* pDouble)
{
        printf("C %32.15f\n", *pDouble);    
}

#ifdef MPI
void gpu_gather_gb_reffa()
{
PRINTMETHOD("gpu_gather_gb_reffa"); 
#ifdef CUDA_P2P
    if (gpu->bSingleNode)
    {
        P2P_Broadcast(gpu->pbReffa->_pDevData, (void**)gpu->pP2PReffa, gpu->pAllGathervRecvCountSoA[gpu->gpuID], gpu->pAllGathervRecvDisplSoA[gpu->gpuID], sizeof(PMEDouble));
    }
    else
#endif
    {
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_PMEDOUBLE, gpu->pbReffa->_pSysData, gpu->pAllGathervRecvCountSoA, gpu->pAllGathervRecvDisplSoA, MPI_PMEDOUBLE, gpu->comm);
    }
}

void gpu_gather_gb_temp7a()
{
PRINTMETHOD("gpu_gather_gb_temp7a");
#ifdef CUDA_P2P
    if (gpu->bSingleNode)
    {
        P2P_Broadcast(gpu->pbTemp7a->_pDevData, (void**)gpu->pP2PTemp7a, gpu->pAllGathervRecvCountSoA[gpu->gpuID], gpu->pAllGathervRecvDisplSoA[gpu->gpuID], sizeof(PMEDouble));
    }
    else
#endif
    {
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_PMEDOUBLE, gpu->pbTemp7a->_pSysData, gpu->pAllGathervRecvCountSoA, gpu->pAllGathervRecvDisplSoA, MPI_PMEDOUBLE, gpu->comm);
    }
}

void gpu_gather_gb_forces()
{
PRINTMETHOD("gpu_gather_gb_forces"); 
#if 0
        {
        static int step = 0;
        char buff[512];
        sprintf(buff, "../nbforces%d_%d.out", step, gpu->gpuID);
        FILE* fp = fopen(buff, "w");
        fprintf(fp, "%06d %06d\n", gpu->sim.minReducedAtom, gpu->sim.maxReducedAtom);
        for (int i = 0; i < gpu->sim.atoms; i++)
        {

                fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbOutForce->_pSysData[i * 3],
                                                            gpu->pbOutForce->_pSysData[i * 3 + 1], 
                                                            gpu->pbOutForce->_pSysData[i * 3 + 2]);                                                       
        }
        fclose(fp);
        step++;
        }
#endif
#ifdef CUDA_P2P   
    if (gpu->bSingleNode)
    {
        P2P_Broadcast(gpu->pbInForce->_pDevData, (void**)gpu->pP2PInForce, gpu->pAllGathervRecvCountAoS[gpu->gpuID], gpu->pAllGathervRecvDisplAoS[gpu->gpuID], sizeof(PMEDouble));
    }
    else
#endif
    {      
        MPI_Allgatherv(gpu->pbOutForce->_pSysData, gpu->pAllGathervRecvCountAoS[gpu->gpuID], MPI_PMEDOUBLE, gpu->pbInForce->_pSysData, gpu->pAllGathervRecvCountAoS, gpu->pAllGathervRecvDisplAoS, MPI_PMEDOUBLE, gpu->comm);
    }
#if 0
        {
        static int step = 0;
        char buff[512];
        sprintf(buff, "../naforces%d_%d.out", step, gpu->gpuID);
        FILE* fp = fopen(buff, "w");
        fprintf(fp, "%06d %06d\n", gpu->sim.minReducedAtom, gpu->sim.maxReducedAtom);
        for (int i = 0; i < gpu->sim.atoms; i++)
        {

                fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbInForce->_pSysData[i * 3],
                                                            gpu->pbInForce->_pSysData[i * 3 + 1], 
                                                            gpu->pbInForce->_pSysData[i * 3 + 2]);                                                       
        }
        fclose(fp);
        step++;
        }
    {
    static int step = 0;
    char buff[512];
    sprintf(buff, "../nreff%d_%d.out", step, gpu->gpuID);
    gpu->pbReff->Download();
    FILE* fp = fopen(buff, "w");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        fprintf(fp, "%6d %20.10f\n", i, gpu->pbReff->_pSysData[i]);
    }
    fclose(fp);
    step++;
    }
    {
    static int step = 0;
    char buff[512];
    sprintf(buff, "../ntemp7%d_%d.out", step, gpu->gpuID);
    gpu->pbTemp7->Download();
    FILE* fp = fopen(buff, "w");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        fprintf(fp, "%6d %20.10f\n", i, gpu->pbTemp7->_pSysData[i]);
    }
    fclose(fp);
    step++;
    }    
#endif            
}
#endif


extern "C" void gpu_gb_ene_(gb_pot_ene_rec* pEnergy)
{
PRINTMETHOD("gpu_gb_ene"); 
    kClearForces(gpu);   
    if (gpu->ntf != 8)
    {   
        kCalculateGBBornRadii(gpu);
        kReduceGBBornRadii(gpu);
#ifdef MPI        
        gpu_gather_gb_reffa();
        kProcessGBBornRadii(gpu);
#endif
        kCalculateGBNonbondEnergy1(gpu);     
        kReduceGBTemp7Energy(gpu);
#ifdef MPI        
        gpu_gather_gb_temp7a();
        if (gpu->gpuID == 0)
            kProcessGBTemp7Energy(gpu);
        else
            kProcessGBTemp7(gpu);        
#endif
        kCalculateGBNonbondEnergy2(gpu);
        kCalculateLocalEnergy(gpu);
        kCalculateCHARMMEnergy(gpu);
        kReduceForces(gpu);   
#ifdef MPI
        kTransposeForces(gpu);
#endif        
    }
#ifdef MPI
    cudaThreadSynchronize();
    gpu_download_partial_forces();
    gpu_gather_gb_forces(); 
#endif          


#if 0
#ifndef MPI
    {
    static int step = 0;
    char buff[512];
    sprintf(buff, "../nreff%d.out", step);
    gpu->pbReff->Download();
    FILE* fp = fopen(buff, "w");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        fprintf(fp, "%6d %20.10f\n", i, gpu->pbReff->_pSysData[i]);
    }
    fclose(fp);
    step++;
    }
    {
    static int step = 0;
    char buff[512];
    sprintf(buff, "../ntemp7%d.out", step);
    gpu->pbTemp7->Download();
    FILE* fp = fopen(buff, "w");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        fprintf(fp, "%6d %20.10f\n", i, gpu->pbTemp7->_pSysData[i]);
    }
    fclose(fp);
    step++;
    }    
    {
    static int step = 0;
    char buff[512];
    sprintf(buff, "../nforces%d.out", step);
    gpu->pbForce->Download();
    FILE* fp = fopen(buff, "w");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbForce->_pSysData[i],
                                                        gpu->pbForce->_pSysData[i + gpu->sim.stride], 
                                                        gpu->pbForce->_pSysData[i + gpu->sim.stride2]);
    }
    fclose(fp);
    step++;
    }
#endif
#endif    
   
    gpu->pbEnergyBuffer->Download();
    PMEDouble energy[ENERGYTERMS];
    pEnergy->total                  = 0.0;
    for (int i = 0; i < ENERGYTERMS; i++)
    {
        unsigned long long int val  = gpu->pbEnergyBuffer->_pSysData[i];
        if (val >= 0x8000000000000000ull)
        {
            energy[i]               = -(PMEDouble)(val ^ 0xffffffffffffffffull) / ENERGYSCALE;
        }
        else
        {
            energy[i]               = (PMEDouble)val / ENERGYSCALE;
        }
        pEnergy->total             += energy[i];
    }
    pEnergy->vdw_tot                = energy[1];
    pEnergy->elec_tot               = energy[0];
    pEnergy->gb                     = energy[2];
    pEnergy->bond                   = energy[3];
    pEnergy->angle                  = energy[4];
    pEnergy->dihedral               = energy[5];
    pEnergy->vdw_14                 = energy[7];
    pEnergy->elec_14                = energy[6];
    pEnergy->restraint              = energy[8];
    pEnergy->angle_ub               = energy[11];
    pEnergy->imp                    = energy[12];
    pEnergy->cmap                   = energy[13];
    
#if 0
    gpu->pbForce->Download();
    PMEDouble xsum = 0.0;
    PMEDouble ysum = 0.0;
    PMEDouble zsum = 0.0;
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        printf("%5d: %32.15f %32.15f %32.15f\n", i, gpu->pbForce->_pSysData[i], gpu->pbForce->_pSysData[i + gpu->sim.stride], gpu->pbForce->_pSysData[i + gpu->sim.stride2]);
        xsum += gpu->pbForce->_pSysData[i];
        ysum += gpu->pbForce->_pSysData[i + gpu->sim.stride];
        zsum += gpu->pbForce->_pSysData[i + gpu->sim.stride2];
    } 
    printf("Sum: %32.15f %32.15f %32.15f\n", xsum, ysum, zsum);
    exit(-1);
#endif  
}

extern "C" void gpu_clear_forces_()
{
PRINTMETHOD("gpu_clear_forces"); 
    kClearForces(gpu);
}   

extern "C" void gpu_gb_forces_()
{
PRINTMETHOD("gpu_gb_forces"); 
    kClearForces(gpu);   
    if (gpu->ntf != 8)
    {
        kCalculateGBBornRadii(gpu);
        kReduceGBBornRadii(gpu);
#ifdef MPI        
        gpu_gather_gb_reffa();
        kProcessGBBornRadii(gpu);
#endif        
        kCalculateGBNonbondForces1(gpu);
        kReduceGBTemp7(gpu);
#ifdef MPI        
        gpu_gather_gb_temp7a();
        kProcessGBTemp7(gpu);
#endif        
        kCalculateGBNonbondEnergy2(gpu);   
        kCalculateLocalForces(gpu);
        kCalculateCHARMMForces(gpu);        
        kReduceForces(gpu);
#ifdef MPI
        kTransposeForces(gpu);
#endif 
    }
#ifdef MPI       
    cudaThreadSynchronize();
    gpu_download_partial_forces();
    gpu_gather_gb_forces();
#endif
#if 0
#ifndef MPI
    {
    char buff[512];
    sprintf(buff, "../naforces.out");
    gpu->pbForce->Download();
    FILE* fp = fopen(buff, "w");
    for (int i = 0; i < gpu->sim.atoms; i++)
    {
        fprintf(fp, "%6d %20.10f %20.10f %20.10f\n", i, gpu->pbForce->_pSysData[i],
                                                        gpu->pbForce->_pSysData[i + gpu->sim.stride], 
                                                        gpu->pbForce->_pSysData[i + gpu->sim.stride2]);
    }
    fclose(fp);
    }
#endif
#endif  
}

extern "C" void gpu_get_nb_energy_()
{
PRINTMETHOD("gpu_get_nb_energy"); 
    kCalculatePMENonbondForces(gpu);
}

extern "C" void gpu_local_to_global_()
{
PRINTMETHOD("gpu_local_to_global");
    kLocalToGlobal(gpu);
    
#if 0
    {
        gpu->pbImageIndex->Download();
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);        
        PMEDouble *pCrd                                 = gpu->pbImage->_pSysData; 
        gpu->pbImage->Download();     
     
        if (gpu->sim.pImageX != gpu->pbImage->_pDevData)
           pCrd                                        = gpu->pbImage->_pSysData + gpu->sim.stride3;

        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            unsigned int j = pImageAtomLookup[i];
            printf("%6d,%32.15f,%32.15f,%32.15f\n", i, pCrd[j], pCrd[j + gpu->sim.stride], pCrd[j + gpu->sim.stride2]);
        }
    }
//    exit(-1);
#endif    
}

extern "C" void gpu_amd_setup_(int* iamd, int* iamdlag, int* ntwx,double* EthreshP, double* alphaP, double* EthreshD, double* alphaD, double* temp0)
{
PRINTMETHOD("gpu_amd_setup");
    // Determine what type of AMD is being used
    gpu->sim.iamd                    = *iamd;
    gpu->sim.iamdlag                    = *iamdlag;
    
    // Allocate GPU data
    // Set up AMD parameters
    gpu->sim.amd_print_interval                    = *ntwx;
    gpu->sim.amd_EthreshP                            = *EthreshP;
    gpu->sim.amd_alphaP                           = *alphaP;
    gpu->sim.amd_EthreshD                             = *EthreshD;
    gpu->sim.amd_alphaD                            = *alphaD;
    gpu->sim.amd_temp0                             = *temp0;
    gpu->pbAmdWeightsAndEnergy                     = new GpuBuffer<PMEDouble>(*ntwx * 6, false, true);
    gpu->sim.pAmdWeightsAndEnergy                  = gpu->pbAmdWeightsAndEnergy->_pDevData;

    gpu->sim.pAmdNumLag                           = 0;

    gpu->sim.pAmdNumRecs                           = -1;

    gpu->sim.pAMDtboost                            = 0.0;

    gpu->pbAMDfwgtd                                = new GpuBuffer<PMEDouble>(1);
    gpu->sim.pAMDfwgtd                             = gpu->pbAMDfwgtd->_pDevData;

    gpu->sim.pAMDfwgt                             = 1.0;

    gpu->pbAMDEDihedral                            = new GpuBuffer<PMEUllInt>(1);
    gpu->sim.pAMDEDihedral                         = gpu->pbAMDEDihedral->_pDevData;
    gpuCopyConstants();
    return;
}



extern "C" void gpu_download_amd_weights_(double amd_weights_and_energy[][6])
{
  PRINTMETHOD("gpu_download_amd_weights");

  for (int i = 0; i < gpu->sim.amd_print_interval; i++)
    {
      amd_weights_and_energy[i][0]                   = gpu->pbAmdWeightsAndEnergy->_pSysData[i*6 + 0];
      amd_weights_and_energy[i][1]                   = gpu->pbAmdWeightsAndEnergy->_pSysData[i*6 + 1];
      amd_weights_and_energy[i][2]                   = gpu->pbAmdWeightsAndEnergy->_pSysData[i*6 + 2];
      amd_weights_and_energy[i][3]                   = gpu->pbAmdWeightsAndEnergy->_pSysData[i*6 + 3];
      amd_weights_and_energy[i][4]                   = gpu->pbAmdWeightsAndEnergy->_pSysData[i*6 + 4];
      amd_weights_and_energy[i][5]                   = gpu->pbAmdWeightsAndEnergy->_pSysData[i*6 + 5];
      
    } 
}

#ifdef MPI
extern "C" void gpu_calculate_amd_dihedral_weight_(double* totdih){
  PRINTMETHOD("gpu_calculate_amd_dihedral_weight");

  //calculate AMD weight for dihedral boost (tboost)
    PMEDouble EthreshD   = gpu->sim.amd_EthreshD; 
    PMEDouble alphaD  = gpu->sim.amd_alphaD;
    PMEDouble tboost = 0.0;
    PMEDouble fwgtd   = 1.0;
    PMEDouble temp0 = gpu->sim.amd_temp0;

    if((*totdih) <= EthreshD) {
      if(gpu->sim.pAmdNumLag == 0){
	tboost = ((EthreshD - (*totdih))*(EthreshD - (*totdih)))*1000.0/ 
	  ((alphaD + (EthreshD - (*totdih)))*temp0*1.987);
	fwgtd = (alphaD*alphaD)/((alphaD + EthreshD - (*totdih))*(alphaD + EthreshD - (*totdih)));
      }
    }

    gpu->sim.pAMDtboost = tboost;
    gpu->pbAMDfwgtd->_pSysData[0] = fwgtd;
    gpu->pbAMDfwgtd->Upload();
    //gpuCopyConstants();
}
extern "C" void gpu_calculate_amd_dihedral_energy_(double* totdih){
  PRINTMETHOD("gpu_calculate_amd_dihedral_weight");

  // Rebuild neighbor list 
  gpu_build_neighbor_list_();
  // Local energy
  if (gpu->sim.minLocalAtom != gpu->sim.maxLocalAtom){
    //calculate AMD energy for dihedral boost (tboost)
    kCalculateAmdDihedralEnergy(gpu);
    gpu->pbAMDEDihedral->Download();
    PMEDouble totdih2;
    if (gpu->pbAMDEDihedral->_pSysData[0] >= 0x8000000000000000ull)
      {
	totdih2               = -(PMEDouble)(gpu->pbAMDEDihedral->_pSysData[0] ^ 0xffffffffffffffffull) / ENERGYSCALE;
      }
    else
      {
	totdih2               = (PMEDouble)gpu->pbAMDEDihedral->_pSysData[0] / ENERGYSCALE;
      }
    *totdih = totdih2;
  }
  else{
    *totdih = 0.0;
  }
}
#else
extern "C" void gpu_calculate_amd_dihedral_energy_weight_(){
  PRINTMETHOD("gpu_calculate_amd_dihedral_energy_weight");

  PMEDouble totdih = 0.0;
  if(gpu->sim.pAmdNumLag == 0){
    // Rebuild neighbor list 
    gpu_build_neighbor_list_();
    //calculate AMD calculate energy for dihedral boost (tboost)
    kCalculateAmdDihedralEnergy(gpu);
    
    gpu->pbAMDEDihedral->Download();
    if (gpu->pbAMDEDihedral->_pSysData[0] >= 0x8000000000000000ull)
      {
	totdih               = -(PMEDouble)(gpu->pbAMDEDihedral->_pSysData[0] ^ 0xffffffffffffffffull) / ENERGYSCALE;
      }
    else
      {
	totdih               = (PMEDouble)gpu->pbAMDEDihedral->_pSysData[0] / ENERGYSCALE;
      }
  }
  //calculate AMD weight for dihedral boost (tboost)
  PMEDouble EthreshD   = gpu->sim.amd_EthreshD; 
  PMEDouble alphaD  = gpu->sim.amd_alphaD;
  PMEDouble tboost = 0.0;
  PMEDouble fwgtd   = 1.0;
  PMEDouble temp0 = gpu->sim.amd_temp0;
  
  if((totdih) <= EthreshD) {
    if(gpu->sim.pAmdNumLag == 0){
      tboost = ((EthreshD - (totdih))*(EthreshD - (totdih)))*1000.0/ 
	((alphaD + (EthreshD - (totdih)))*temp0*1.987);
      fwgtd = (alphaD*alphaD)/((alphaD + EthreshD - (totdih))*(alphaD + EthreshD - (totdih)));
    }
  }
  
  gpu->sim.pAMDtboost = tboost;
  gpu->pbAMDfwgtd->_pSysData[0] = fwgtd;
  gpu->pbAMDfwgtd->Upload();
  
}
#endif

extern "C" void gpu_calculate_and_apply_amd_weights_(double* pot_ene_tot, double* dih_ene_tot, double* num_amd_lag)
{
  
  PRINTMETHOD("gpu_calculate_and_apply_amd_weights");
  PMEDouble tboost = 0.0;
  PMEDouble fwgtd = 1.0;
  PMEDouble fwgt   = 1.0;
  PMEDouble tboostall = 0.0;
  
  if(gpu->sim.pAmdNumLag == 0){
    
    //calculate AMD weight, seting dihedral boost (tboost) to zero for now
    PMEDouble EthreshP  =  gpu->sim.amd_EthreshP;     
    PMEDouble alphaP = gpu->sim.amd_alphaP;   
    PMEDouble temp0 = gpu->sim.amd_temp0;   
    
    if ((gpu->sim.iamd == 2)||(gpu->sim.iamd == 3)){
      tboost = gpu->sim.pAMDtboost;
      fwgtd   = gpu->pbAMDfwgtd->_pSysData[0];
    }
    PMEDouble totalenergy = *pot_ene_tot + (tboost * (temp0)*1.987)/1000.0;
    
    if (((gpu->sim.iamd == 1)||(gpu->sim.iamd == 3)) && (totalenergy <= EthreshP)) {
      tboostall = ((EthreshP - totalenergy)*(EthreshP - totalenergy))*1000.0;
      tboostall = tboostall/((alphaP + (EthreshP - totalenergy))*temp0*1.987);
      fwgt = (alphaP*alphaP)/((alphaP + EthreshP - totalenergy)*(alphaP + EthreshP - totalenergy));
    }  
    
    //calculate AMD weight
    
#ifdef MPI
    for (int i = 0; i < gpu->sim.atoms; i++){
      gpu->pbInForce->_pSysData[i * 3]       *= fwgt;
      gpu->pbInForce->_pSysData[i * 3 + 1]   *= fwgt;
      gpu->pbInForce->_pSysData[i * 3 + 2]   *= fwgt;
    }
    
#else
    kAMDCalcWeightAndScaleForces(gpu, *pot_ene_tot, *dih_ene_tot, fwgt);
#endif   
    
  }
  
  
  // Output sum if thread 0
  int numrecs = gpu->sim.pAmdNumRecs * 6;
  
  if(numrecs >= 0){
    gpu->pbAmdWeightsAndEnergy->_pSysData[numrecs + 0] = *pot_ene_tot; 
    gpu->pbAmdWeightsAndEnergy->_pSysData[numrecs + 1] = *dih_ene_tot;
    gpu->pbAmdWeightsAndEnergy->_pSysData[numrecs + 2] = fwgt;
    gpu->pbAmdWeightsAndEnergy->_pSysData[numrecs + 3] = fwgtd;
    gpu->pbAmdWeightsAndEnergy->_pSysData[numrecs + 4] = tboostall;
    gpu->pbAmdWeightsAndEnergy->_pSysData[numrecs + 5] = tboost;
  }

 if(gpu->sim.pAmdNumLag == gpu->sim.iamdlag){
    gpu->sim.pAmdNumLag = 0;
  }
  else{
    gpu->sim.pAmdNumLag = gpu->sim.pAmdNumLag +1;
  } 
  
  
  gpu->sim.pAmdNumRecs = gpu->sim.pAmdNumRecs +1;
  if (gpu->sim.pAmdNumRecs  >= gpu->sim.amd_print_interval){
    //set pAmdNumRecs back to 0 
    gpu->sim.pAmdNumRecs = 0 ; 
  }
  
  *num_amd_lag = gpu->sim.pAmdNumLag;  
}

//AMD functions for gb, no need to do neighbor lists here
#ifdef MPI
extern "C" void gpu_calculate_gb_amd_dihedral_energy_(double* totdih){
  PRINTMETHOD("gpu_calculate_amd_dihedral_weight");

  // Local energy
  if (gpu->sim.minLocalAtom != gpu->sim.maxLocalAtom){
    //calculate AMD energy for dihedral boost (tboost)
    kCalculateAmdDihedralEnergy(gpu);
    gpu->pbAMDEDihedral->Download();
    PMEDouble totdih2;
    if (gpu->pbAMDEDihedral->_pSysData[0] >= 0x8000000000000000ull)
      {
	totdih2               = -(PMEDouble)(gpu->pbAMDEDihedral->_pSysData[0] ^ 0xffffffffffffffffull) / ENERGYSCALE;
      }
    else
      {
	totdih2               = (PMEDouble)gpu->pbAMDEDihedral->_pSysData[0] / ENERGYSCALE;
      }
    *totdih = totdih2;
  }
  else{
    *totdih = 0.0;
  }
}
#else
extern "C" void gpu_calculate_gb_amd_dihedral_energy_weight_(){
  PRINTMETHOD("gpu_calculate_amd_dihedral_energy_weight");

  PMEDouble totdih = 0.0;
  if(gpu->sim.pAmdNumLag == 0){
    //calculate AMD calculate energy for dihedral boost (tboost)
    kCalculateAmdDihedralEnergy(gpu);
    
    gpu->pbAMDEDihedral->Download();
    if (gpu->pbAMDEDihedral->_pSysData[0] >= 0x8000000000000000ull)
      {
	totdih               = -(PMEDouble)(gpu->pbAMDEDihedral->_pSysData[0] ^ 0xffffffffffffffffull) / ENERGYSCALE;
      }
    else
      {
	totdih               = (PMEDouble)gpu->pbAMDEDihedral->_pSysData[0] / ENERGYSCALE;
      }
  }
  //calculate AMD weight for dihedral boost (tboost)
  PMEDouble EthreshD   = gpu->sim.amd_EthreshD; 
  PMEDouble alphaD  = gpu->sim.amd_alphaD;
  PMEDouble tboost = 0.0;
  PMEDouble fwgtd   = 1.0;
  PMEDouble temp0 = gpu->sim.amd_temp0;
  
  if((totdih) <= EthreshD) {
    if(gpu->sim.pAmdNumLag == 0){
      tboost = ((EthreshD - (totdih))*(EthreshD - (totdih)))*1000.0/ 
	((alphaD + (EthreshD - (totdih)))*temp0*1.987);
      fwgtd = (alphaD*alphaD)/((alphaD + EthreshD - (totdih))*(alphaD + EthreshD - (totdih)));
    }
  }
  
  gpu->sim.pAMDtboost = tboost;
  gpu->pbAMDfwgtd->_pSysData[0] = fwgtd;
  gpu->pbAMDfwgtd->Upload();
  
}
#endif
