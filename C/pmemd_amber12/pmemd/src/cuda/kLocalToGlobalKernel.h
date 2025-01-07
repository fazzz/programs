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
    unsigned int pos                        = blockIdx.x * blockDim.x + threadIdx.x;
#define EP_ONEPOINT
    while (pos < cSim.EP11Offset)
    {    
#define EP_TYPE1
        if (pos < cSim.EP11s)
        {
#include "kLocalToGlobalLoop.h"
        }
#undef EP_TYPE1
        pos                                += blockDim.x * gridDim.x;
    }
    
    while (pos < cSim.EP12Offset)
    {
        pos                                -= cSim.EP11Offset;
#define EP_TYPE2
        if (pos < cSim.EP12s)
        {    
#include "kLocalToGlobalLoop.h"
        }
#undef EP_TYPE2
        pos                                += cSim.EP11Offset + blockDim.x * gridDim.x;
    }
#undef EP_ONEPOINT

#define EP_TWOPOINTS
    while (pos < cSim.EP21Offset)
    {  
        pos                                -= cSim.EP12Offset;      
#define EP_TYPE1
        if (pos < cSim.EP21s)
        {
#include "kLocalToGlobalLoop.h"
        }
#undef EP_TYPE1
        pos                                += cSim.EP12Offset + blockDim.x * gridDim.x;
    }
    
    while (pos < cSim.EP22Offset)
    {
        pos                                -= cSim.EP21Offset;
#define EP_TYPE2
        if (pos < cSim.EP22s)
        {    
#include "kLocalToGlobalLoop.h"
        }
#undef EP_TYPE2
        pos                                += cSim.EP21Offset + blockDim.x * gridDim.x;
    }
#undef EP_TWOPOINTS   
}
