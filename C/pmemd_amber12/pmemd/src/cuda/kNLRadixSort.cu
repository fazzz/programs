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

#include <cuda.h>
#include "gpu.h"
#include <radix_sort/enactor.cuh>
#include <util/multi_buffer.cuh>

static b40c::radix_sort::Enactor* pSortEnactor = NULL;
static b40c::util::MultiBuffer<2, unsigned int, unsigned int>* pDeviceStorage = NULL;


extern "C" void kNLDeleteRadixSort(gpuContext gpu)
{
    if (pSortEnactor)
        delete pSortEnactor;
    pSortEnactor                                = NULL; 
    if (pDeviceStorage)
        delete pDeviceStorage;
    pDeviceStorage                              = NULL;
}

extern "C" void EnactSort(int bits, int sort_atoms)
{
	switch (bits)
	{
	    case 9:
	    case 10:
	        (*pSortEnactor).Sort<b40c::radix_sort::SMALL_PROBLEM, 10, 0>(*pDeviceStorage, sort_atoms);
	        break;
	    case 11:
	    case 12:	    
	        (*pSortEnactor).Sort<b40c::radix_sort::SMALL_PROBLEM, 12, 0>(*pDeviceStorage, sort_atoms);
	        break;
	    case 13:
	    case 14:
	    case 15:
	    case 16:
	        (*pSortEnactor).Sort<b40c::radix_sort::SMALL_PROBLEM, 16, 0>(*pDeviceStorage, sort_atoms);
	        break;
	    case 17:
	    case 18:
	    case 19:
	    case 20:
	        (*pSortEnactor).Sort<b40c::radix_sort::SMALL_PROBLEM, 20, 0>(*pDeviceStorage, sort_atoms);
	        break;
	    case 21:
	    case 22:
	    case 23:
	    case 24:
	        (*pSortEnactor).Sort<b40c::radix_sort::SMALL_PROBLEM, 24, 0>(*pDeviceStorage, sort_atoms);
	        break;	        	        
	    case 25:
	    case 26:
	    case 27:
	    case 28:
	        (*pSortEnactor).Sort<b40c::radix_sort::SMALL_PROBLEM, 28, 0>(*pDeviceStorage, sort_atoms);
	        break;
	    case 29:
	    case 30:
	        (*pSortEnactor).Sort<b40c::radix_sort::SMALL_PROBLEM, 30, 0>(*pDeviceStorage, sort_atoms);
	        break;
	    case 31:
	    case 32:
	        (*pSortEnactor).Sort<b40c::radix_sort::SMALL_PROBLEM, 32, 0>(*pDeviceStorage, sort_atoms);
	        break;	        
	}
}


extern "C" void kNLInitRadixSort(gpuContext gpu)
{


    // Delete old Radix sort
    cudaDeviceSynchronize();
    kNLDeleteRadixSort(gpu);
    cudaDeviceSynchronize();
   
    // Create new sort
    pDeviceStorage                              = new b40c::util::MultiBuffer<2, unsigned int, unsigned int>();
    pSortEnactor                                = new b40c::radix_sort::Enactor();
    cudaDeviceSynchronize();
    pDeviceStorage->d_keys[pDeviceStorage->selector]                   = gpu->sim.pImageHash;
	pDeviceStorage->d_values[pDeviceStorage->selector]                 = gpu->sim.pImageIndex;
    pDeviceStorage->d_keys[pDeviceStorage->selector ^ 1]               = gpu->sim.pImageHash2;
	pDeviceStorage->d_values[pDeviceStorage->selector ^ 1]             = gpu->sim.pImageIndex2;
    cudaDeviceSynchronize();
    EnactSort(gpu->neighborListBits, gpu->sim.atoms);
    cudaDeviceSynchronize();
}




extern "C" void kNLRadixSort(gpuContext gpu)
{
    pDeviceStorage->d_keys[pDeviceStorage->selector]                   = gpu->sim.pImageHash;
	pDeviceStorage->d_values[pDeviceStorage->selector]                 = gpu->sim.pImageIndex;
    pDeviceStorage->d_keys[pDeviceStorage->selector ^ 1]               = gpu->sim.pImageHash2;
	pDeviceStorage->d_values[pDeviceStorage->selector ^ 1]             = gpu->sim.pImageIndex2;

    EnactSort(gpu->neighborListBits, gpu->sim.atoms);

    gpu->sim.pImageHash             = pDeviceStorage->d_keys[pDeviceStorage->selector];
	gpu->sim.pImageIndex            = pDeviceStorage->d_values[pDeviceStorage->selector];
    gpu->sim.pImageHash2            = pDeviceStorage->d_keys[pDeviceStorage->selector ^ 1];
	gpu->sim.pImageIndex2           = pDeviceStorage->d_values[pDeviceStorage->selector ^ 1];
}

