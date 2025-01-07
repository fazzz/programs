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
#include <radixsort_single_grid.cu>
#include <radixsort_early_exit.cu>		

using namespace b40c;
static SingleGridRadixSortingEnactor<unsigned int, unsigned int>* pSortEnactor = NULL;
static MultiCtaRadixSortStorage<unsigned int, unsigned int>* pDeviceStorage = NULL;
extern "C" void kNLDeleteRadixSort(gpuContext gpu)
{
    if (pSortEnactor)
        delete pSortEnactor;
    pSortEnactor                                = NULL; 
    if (pDeviceStorage)
        delete pDeviceStorage;
    pDeviceStorage                              = NULL;
}

extern "C" void EnactSort(int bits)
{
	switch (bits)
	{
	    case 9:
	    case 10:
	    case 11:
	    case 12:	    
	        (*pSortEnactor).EnactSort<12>(*pDeviceStorage);
	        break;
	        
	    case 13:
	    case 14:
	    case 15:
	    case 16:
	        (*pSortEnactor).EnactSort<16>(*pDeviceStorage);
	        break;

	    case 17:
	    case 18:
	    case 19:
	    case 20:
	        (*pSortEnactor).EnactSort<20>(*pDeviceStorage);
	        break;

	    case 21:
	    case 22:
	    case 23:
	    case 24:
	        (*pSortEnactor).EnactSort<24>(*pDeviceStorage);
	        break;	        	        

	    case 25:
	    case 26:
	    case 27:
	    case 28:
	        (*pSortEnactor).EnactSort<28>(*pDeviceStorage);
	        break;

	    case 29:
	    case 30:
	    case 31:
	    case 32:
	        (*pSortEnactor).EnactSort<32>(*pDeviceStorage);
	        break;	        
	}
}


extern "C" void kNLInitRadixSort(gpuContext gpu)
{
    // Delete old Radix sort
    kNLDeleteRadixSort(gpu);
   
    // Create new sort
    pDeviceStorage                              = new MultiCtaRadixSortStorage<unsigned int, unsigned int>(gpu->sim.atoms);
    pSortEnactor                                = new SingleGridRadixSortingEnactor<unsigned int, unsigned int>;
    pDeviceStorage->d_keys[0]                   = gpu->sim.pImageHash;
	pDeviceStorage->d_values[0]                 = gpu->sim.pImageIndex;
    pDeviceStorage->d_keys[1]                   = gpu->sim.pImageHash2;
	pDeviceStorage->d_values[1]                 = gpu->sim.pImageIndex2;	
    EnactSort(gpu->neighborListBits);
}



extern "C" void kNLRadixSort(gpuContext gpu)
{
    pDeviceStorage->d_keys[0]                   = gpu->sim.pImageHash;
	pDeviceStorage->d_values[0]                 = gpu->sim.pImageIndex;
    pDeviceStorage->d_keys[1]                   = gpu->sim.pImageHash2;
	pDeviceStorage->d_values[1]                 = gpu->sim.pImageIndex2;	
    EnactSort(gpu->neighborListBits);
    gpu->sim.pImageHash                         = pDeviceStorage->d_keys[pDeviceStorage->selector];           
	gpu->sim.pImageIndex                        = pDeviceStorage->d_values[pDeviceStorage->selector];
    gpu->sim.pImageHash2                        = pDeviceStorage->d_keys[1 - pDeviceStorage->selector];           
	gpu->sim.pImageIndex2                       = pDeviceStorage->d_values[1 - pDeviceStorage->selector];
}

