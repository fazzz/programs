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

static __constant__ cudaSimulation cSim;
struct Atom 
{
    double invMassI;
    double xpl;
    double ypl;
    double zpl;
    double xil;
    double yil;
    double zil;
};

// Texture reference for double-precision coordinates (disguised as int2 to work around HW limitations)
#ifndef use_SPSP
texture<int2, 1, cudaReadModeElementType> texref;
#else
texture<float, 1, cudaReadModeElementType> texref;
#endif

void SetkShakeSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

void GetkShakeSim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_SHAKE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_SHAKE_THREADS_PER_BLOCK, 1)
#endif
kShake_kernel()
{
#if (__CUDA_ARCH__ >= 200)
    __shared__ Atom sA[SM_2X_SHAKE_THREADS_PER_BLOCK];
#else
    __shared__ Atom sA[SM_13_SHAKE_THREADS_PER_BLOCK];
#endif
    Atom* psA                                   = &sA[threadIdx.x];
    unsigned int pos                            = blockIdx.x * blockDim.x + threadIdx.x;   
    while (pos < cSim.shakeConstraints)
    {
    
        // Read SHAKE network data
        int4 shakeID                            = cSim.pShakeID[pos];
        double2 shakeParm                       = cSim.pShakeParm[pos];
        
        // Read SHAKE network components
#ifdef use_SPSP        
        double xi                               = tex1Dfetch(texref, shakeID.x);
        double yi                               = tex1Dfetch(texref, shakeID.x + cSim.stride);
        double zi                               = tex1Dfetch(texref, shakeID.x + cSim.stride2);
        double xij                              = tex1Dfetch(texref, shakeID.y);
        double yij                              = tex1Dfetch(texref, shakeID.y + cSim.stride);
        double zij                              = tex1Dfetch(texref, shakeID.y + cSim.stride2);
#elif defined(NODPTEXTURE) && (__CUDA_ARCH__ < 200)
        double xi                               = cSim.pForceX[shakeID.x];
        double yi                               = cSim.pForceY[shakeID.x];
        double zi                               = cSim.pForceZ[shakeID.x];
        double xij                              = cSim.pForceX[shakeID.y];
        double yij                              = cSim.pForceY[shakeID.y];
        double zij                              = cSim.pForceZ[shakeID.y];
#else        
        int2 ixi                                = tex1Dfetch(texref, shakeID.x);
        int2 iyi                                = tex1Dfetch(texref, shakeID.x + cSim.stride);
        int2 izi                                = tex1Dfetch(texref, shakeID.x + cSim.stride2);
        int2 ixij                               = tex1Dfetch(texref, shakeID.y);
        int2 iyij                               = tex1Dfetch(texref, shakeID.y + cSim.stride);
        int2 izij                               = tex1Dfetch(texref, shakeID.y + cSim.stride2);
#endif                  
        double xpi                              = cSim.pAtomX[shakeID.x];
        double ypi                              = cSim.pAtomY[shakeID.x];
        double zpi                              = cSim.pAtomZ[shakeID.x];
        double xpj                              = cSim.pAtomX[shakeID.y];
        double ypj                              = cSim.pAtomY[shakeID.y];
        double zpj                              = cSim.pAtomZ[shakeID.y];
#if !defined(use_SPSP) && (!defined(NODPTEXTURE) || (__CUDA_ARCH__ >= 200))
        double xi                               = __hiloint2double(ixi.y, ixi.x);
        double yi                               = __hiloint2double(iyi.y, iyi.x);
        double zi                               = __hiloint2double(izi.y, izi.x);
        double xij                              = __hiloint2double(ixij.y, ixij.x);
        double yij                              = __hiloint2double(iyij.y, iyij.x);
        double zij                              = __hiloint2double(izij.y, izij.x); 
#endif               
        psA->invMassI                           = shakeParm.x;
        double toler                            = shakeParm.y;
        
        
        // Optionally read 2nd hydrogen
        double xpk, ypk, zpk, xik, yik, zik;
        if (shakeID.z != -1)
        {
#ifdef use_SPSP
            xik                                 = tex1Dfetch(texref, shakeID.z);
            yik                                 = tex1Dfetch(texref, shakeID.z + cSim.stride);
            zik                                 = tex1Dfetch(texref, shakeID.z + cSim.stride2); 
#elif defined(NODPTEXTURE) && (__CUDA_ARCH__ < 200)
            xik                                 = cSim.pForceX[shakeID.z];
            yik                                 = cSim.pForceY[shakeID.z];
            zik                                 = cSim.pForceZ[shakeID.z];    
#else
            int2 ixik                           = tex1Dfetch(texref, shakeID.z);
            int2 iyik                           = tex1Dfetch(texref, shakeID.z + cSim.stride);
            int2 izik                           = tex1Dfetch(texref, shakeID.z + cSim.stride2);
#endif 
            xpk                                 = cSim.pAtomX[shakeID.z];
            ypk                                 = cSim.pAtomY[shakeID.z];
            zpk                                 = cSim.pAtomZ[shakeID.z];
#if !defined(use_SPSP) && (!defined(NODPTEXTURE) || (__CUDA_ARCH__ >= 200))
            xik                                 = __hiloint2double(ixik.y, ixik.x);
            yik                                 = __hiloint2double(iyik.y, iyik.x);
            zik                                 = __hiloint2double(izik.y, izik.x);  
#endif                  
        }
        
        // Optionally read 3rd hydrogen into shared memory
        if (shakeID.w != -1)
        {
#ifdef use_SPSP
            psA->xil                            = tex1Dfetch(texref, shakeID.w);
            psA->yil                            = tex1Dfetch(texref, shakeID.w + cSim.stride);
            psA->zil                            = tex1Dfetch(texref, shakeID.w + cSim.stride2); 
#elif defined(NODPTEXTURE) && (__CUDA_ARCH__ < 200)     
            psA->xil                            = cSim.pForceX[shakeID.w];
            psA->yil                            = cSim.pForceY[shakeID.w];
            psA->zil                            = cSim.pForceZ[shakeID.w];       
#else            
            int2 ixil                           = tex1Dfetch(texref, shakeID.w);
            int2 iyil                           = tex1Dfetch(texref, shakeID.w + cSim.stride);
            int2 izil                           = tex1Dfetch(texref, shakeID.w + cSim.stride2); 
#endif               
            psA->xpl                            = cSim.pAtomX[shakeID.w];
            psA->ypl                            = cSim.pAtomY[shakeID.w];
            psA->zpl                            = cSim.pAtomZ[shakeID.w]; 
#if !defined(use_SPSP) && (!defined(NODPTEXTURE) || (__CUDA_ARCH__ >= 200))
            psA->xil                            = __hiloint2double(ixil.y, ixil.x);
            psA->yil                            = __hiloint2double(iyil.y, iyil.x);
            psA->zil                            = __hiloint2double(izil.y, izil.x);   
#endif          
        }
        
        // Calculate unchanging quantities
        xij                                     = xi - xij;
        yij                                     = yi - yij;
        zij                                     = zi - zij;
        
        if (shakeID.z != -1)
        {
            xik                                 = xi - xik;
            yik                                 = yi - yik;
            zik                                 = zi - zik;
        }        
         
        if (shakeID.w != -1)
        {
            psA->xil                            = xi - psA->xil;
            psA->yil                            = yi - psA->yil;
            psA->zil                            = zi - psA->zil;
        }      
       
        bool done                               = false;
        for (int i = 0; i < 3000; i++)
        {
            done = true;
            
            // Calculate nominal distance squared
            double xpxx                         = xpi - xpj;
            double ypxx                         = ypi - ypj;
            double zpxx                         = zpi - zpj;
            double rpxx2                        = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
      
            // Apply correction
            double diff                         = toler - rpxx2;
            if (abs(diff) >= toler * cSim.tol)
            {
                done                            = false;
               
                // Shake resetting of coordinate is done here
                double rrpr                     = xij * xpxx + yij * ypxx + zij * zpxx;     
                if (rrpr >= toler * 1.0e-06)
                {
                
                    double acor                 = diff / (rrpr * 2.0 * (psA->invMassI + cSim.invMassH));
                    double h                    = xij * acor;
                    xpi                        += h * psA->invMassI;
                    xpj                        -= h * cSim.invMassH;
                    h                           = yij * acor;
                    ypi                        += h * psA->invMassI;
                    ypj                        -= h * cSim.invMassH;
                    h                           = zij * acor;
                    zpi                        += h * psA->invMassI;
                    zpj                        -= h * cSim.invMassH;             
                }
            }
      
            // Second bond if present
            if (shakeID.z != -1)
            {
                xpxx                            = xpi - xpk;
                ypxx                            = ypi - ypk;
                zpxx                            = zpi - zpk;
                rpxx2                           = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
      
                // Apply correction
                diff                            = toler - rpxx2;
                if (abs(diff) >= toler * cSim.tol)
                {
                    done                        = false;
               
                    // Shake resetting of coordinate is done here
                    double rrpr              = xik * xpxx + yik * ypxx + zik * zpxx;     
                    if (rrpr >= toler * 1.0e-06)
                    {
                
                        double acor          = diff / (rrpr * 2.0 * (psA->invMassI + cSim.invMassH));
                        double h             = xik * acor;
                        xpi                    += h * psA->invMassI;
                        xpk                    -= h * cSim.invMassH;
                        h                       = yik * acor;
                        ypi                    += h * psA->invMassI;
                        ypk                    -= h * cSim.invMassH;
                        h                       = zik * acor;
                        zpi                    += h * psA->invMassI;
                        zpk                    -= h * cSim.invMassH;             
                    }
                }
            }
            
            // Third bond if present
            if (shakeID.w != -1)
            {
                xpxx                            = xpi - psA->xpl;
                ypxx                            = ypi - psA->ypl;
                zpxx                            = zpi - psA->zpl;
                rpxx2                           = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
      
                // Apply correction
                diff                            = toler - rpxx2;
                if (abs(diff) >= toler * cSim.tol)
                {
                    done                        = false;
               
                    // Shake resetting of coordinate is done here
                    double rrpr              = psA->xil * xpxx + psA->yil * ypxx + psA->zil * zpxx;     
                    if (rrpr >= toler * 1.0e-06)
                    {
                
                        double acor             = diff / (rrpr * 2.0 * (psA->invMassI + cSim.invMassH));
                        double h                = psA->xil * acor;
                        xpi                    += h * psA->invMassI;
                        psA->xpl               -= h * cSim.invMassH;
                        h                       = psA->yil * acor;
                        ypi                    += h * psA->invMassI;
                        psA->ypl               -= h * cSim.invMassH;
                        h                       = psA->zil * acor;
                        zpi                    += h * psA->invMassI;
                        psA->zpl               -= h * cSim.invMassH;             
                    }
                }
            }
            
            
            // Check for convergence
            if (done)
                break;
        }
      
        // Write out results if converged, but there's no really good
        // way to indicate failure so we'll let the simulation heading
        // off to Neptune do that for us.  Wish there were a better way,
        // but until the CPU needs something from the GPU, those are the
        // the breaks.  I guess, technically, we could just set a flag to NOP
        // the simulation from here and then carry that result through upon
        // the next ntpr, ntwc, or ntwx update, but I leave that up to you 
        // guys to implement that (or not). 
        if (done)
        {
            cSim.pAtomX[shakeID.x]              = xpi;
            cSim.pAtomY[shakeID.x]              = ypi;
            cSim.pAtomZ[shakeID.x]              = zpi;
            PMEFloat2 xyi                       = {xpi, ypi};
            cSim.pAtomXYSP[shakeID.x]           = xyi;
            cSim.pAtomZSP[shakeID.x]            = zpi;
            
            cSim.pAtomX[shakeID.y]              = xpj;
            cSim.pAtomY[shakeID.y]              = ypj;
            cSim.pAtomZ[shakeID.y]              = zpj;
            PMEFloat2 xyj                       = {xpj, ypj};
            cSim.pAtomXYSP[shakeID.y]           = xyj;
            cSim.pAtomZSP[shakeID.y]            = zpj;

            if (shakeID.z != -1)
            {
                cSim.pAtomX[shakeID.z]          = xpk;
                cSim.pAtomY[shakeID.z]          = ypk;
                cSim.pAtomZ[shakeID.z]          = zpk;
                PMEFloat2 xyk                   = {xpk, ypk};
                cSim.pAtomXYSP[shakeID.z]       = xyk;
                cSim.pAtomZSP[shakeID.z]        = zpk;
            }
        
            if (shakeID.w != -1)
            {
                cSim.pAtomX[shakeID.w]          = psA->xpl;
                cSim.pAtomY[shakeID.w]          = psA->ypl;
                cSim.pAtomZ[shakeID.w]          = psA->zpl;
                PMEFloat2 xyl                   = {psA->xpl, psA->ypl};
                cSim.pAtomXYSP[shakeID.w]       = xyl;
                cSim.pAtomZSP[shakeID.w]        = psA->zpl;
            }
        
        }
        
        pos                                    += gridDim.x * blockDim.x;
    }

    if (cSim.fastShakeConstraints > 0)
    {    
        while (pos < cSim.shakeOffset)
        {
            pos                                += gridDim.x * blockDim.x;
        }
        pos                                    -= cSim.shakeOffset;

        while (pos < cSim.fastShakeConstraints)
        {
            
            // Read atom data
            int4 shakeID                        = cSim.pFastShakeID[pos];
#ifdef use_SPSP
            double x1                           = tex1Dfetch(texref, shakeID.x);
            double y1                           = tex1Dfetch(texref, shakeID.x + cSim.stride);
            double z1                           = tex1Dfetch(texref, shakeID.x + cSim.stride2);
            double x2                           = tex1Dfetch(texref, shakeID.y);
            double y2                           = tex1Dfetch(texref, shakeID.y + cSim.stride);
            double z2                           = tex1Dfetch(texref, shakeID.y + cSim.stride2);
            double x3                           = tex1Dfetch(texref, shakeID.z);
            double y3                           = tex1Dfetch(texref, shakeID.z + cSim.stride);
            double z3                           = tex1Dfetch(texref, shakeID.z + cSim.stride2);  
#elif defined(NODPTEXTURE) && (__CUDA_ARCH__ < 200)
            double x1                           = cSim.pForceX[shakeID.x];
            double y1                           = cSim.pForceY[shakeID.x];
            double z1                           = cSim.pForceZ[shakeID.x]; 
            double x2                           = cSim.pForceX[shakeID.y];
            double y2                           = cSim.pForceY[shakeID.y];
            double z2                           = cSim.pForceZ[shakeID.y]; 
            double x3                           = cSim.pForceX[shakeID.z];
            double y3                           = cSim.pForceY[shakeID.z];
            double z3                           = cSim.pForceZ[shakeID.z]; 
#else           
            int2 ix1                            = tex1Dfetch(texref, shakeID.x);
            int2 iy1                            = tex1Dfetch(texref, shakeID.x + cSim.stride);
            int2 iz1                            = tex1Dfetch(texref, shakeID.x + cSim.stride2);
            int2 ix2                            = tex1Dfetch(texref, shakeID.y);
            int2 iy2                            = tex1Dfetch(texref, shakeID.y + cSim.stride);
            int2 iz2                            = tex1Dfetch(texref, shakeID.y + cSim.stride2);
            int2 ix3                            = tex1Dfetch(texref, shakeID.z);
            int2 iy3                            = tex1Dfetch(texref, shakeID.z + cSim.stride);
            int2 iz3                            = tex1Dfetch(texref, shakeID.z + cSim.stride2);  
#endif
            double xp1                          = cSim.pAtomX[shakeID.x];
            double yp1                          = cSim.pAtomY[shakeID.x];
            double zp1                          = cSim.pAtomZ[shakeID.x];
            double xp2                          = cSim.pAtomX[shakeID.y];
            double yp2                          = cSim.pAtomY[shakeID.y];
            double zp2                          = cSim.pAtomZ[shakeID.y];
            double xp3                          = cSim.pAtomX[shakeID.z];
            double yp3                          = cSim.pAtomY[shakeID.z];
            double zp3                          = cSim.pAtomZ[shakeID.z];
#if !defined(use_SPSP) && (!defined(NODPTEXTURE) || (__CUDA_ARCH__ >= 200))
            double x1                           = __hiloint2double(ix1.y, ix1.x);
            double y1                           = __hiloint2double(iy1.y, iy1.x);
            double z1                           = __hiloint2double(iz1.y, iz1.x);
            double x2                           = __hiloint2double(ix2.y, ix2.x);
            double y2                           = __hiloint2double(iy2.y, iy2.x);
            double z2                           = __hiloint2double(iz2.y, iz2.x);
            double x3                           = __hiloint2double(ix3.y, ix3.x);
            double y3                           = __hiloint2double(iy3.y, iy3.x);
            double z3                           = __hiloint2double(iz3.y, iz3.x);
#endif
    
            // Step1  A1_prime:
            double xb0                          = x2 - x1;
            double yb0                          = y2 - y1;
            double zb0                          = z2 - z1;
            double xc0                          = x3 - x1;
            double yc0                          = y3 - y1;
            double zc0                          = z3 - z1;

            psA->xpl                            = xp1 * cSim.wo_div_wohh + (xp2 + xp3) * cSim.wh_div_wohh;
            psA->ypl                            = yp1 * cSim.wo_div_wohh + (yp2 + yp3) * cSim.wh_div_wohh;
            psA->zpl                            = zp1 * cSim.wo_div_wohh + (zp2 + zp3) * cSim.wh_div_wohh;

            double xa1                          = xp1 - psA->xpl;
            double ya1                          = yp1 - psA->ypl;
            double za1                          = zp1 - psA->zpl;
            double xb1                          = xp2 - psA->xpl;
            double yb1                          = yp2 - psA->ypl;
            double zb1                          = zp2 - psA->zpl;
            double xc1                          = xp3 - psA->xpl;
            double yc1                          = yp3 - psA->ypl;
            double zc1                          = zp3 - psA->zpl;

            double xakszd                       = yb0 * zc0 - zb0 * yc0;
            double yakszd                       = zb0 * xc0 - xb0 * zc0;
            double zakszd                       = xb0 * yc0 - yb0 * xc0;
            double xaksxd                       = ya1 * zakszd - za1 * yakszd;
            double yaksxd                       = za1 * xakszd - xa1 * zakszd;
            double zaksxd                       = xa1 * yakszd - ya1 * xakszd;
            double xaksyd                       = yakszd * zaksxd - zakszd * yaksxd;
            double yaksyd                       = zakszd * xaksxd - xakszd * zaksxd;
            double zaksyd                       = xakszd * yaksxd - yakszd * xaksxd;

            double axlng_inv                    = rsqrt(xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd);
            double aylng_inv                    = rsqrt(xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd);
            double azlng_inv                    = rsqrt(xakszd * xakszd + yakszd * yakszd + zakszd * zakszd);

            double trns11                       = xaksxd * axlng_inv;
            double trns21                       = yaksxd * axlng_inv;
            double trns31                       = zaksxd * axlng_inv;
            double trns12                       = xaksyd * aylng_inv;
            double trns22                       = yaksyd * aylng_inv;
            double trns32                       = zaksyd * aylng_inv;
            double trns13                       = xakszd * azlng_inv;
            double trns23                       = yakszd * azlng_inv;
            double trns33                       = zakszd * azlng_inv;

            double xb0d                         = trns11 * xb0 + trns21 * yb0 + trns31 * zb0;
            double yb0d                         = trns12 * xb0 + trns22 * yb0 + trns32 * zb0;
            double xc0d                         = trns11 * xc0 + trns21 * yc0 + trns31 * zc0;
            double yc0d                         = trns12 * xc0 + trns22 * yc0 + trns32 * zc0;
            double za1d                         = trns13 * xa1 + trns23 * ya1 + trns33 * za1;
            double xb1d                         = trns11 * xb1 + trns21 * yb1 + trns31 * zb1;
            double yb1d                         = trns12 * xb1 + trns22 * yb1 + trns32 * zb1;
            double zb1d                         = trns13 * xb1 + trns23 * yb1 + trns33 * zb1;
            double xc1d                         = trns11 * xc1 + trns21 * yc1 + trns31 * zc1;
            double yc1d                         = trns12 * xc1 + trns22 * yc1 + trns32 * zc1;
            double zc1d                         = trns13 * xc1 + trns23 * yc1 + trns33 * zc1;

            // Step2  A2_prime:
            double sinphi                       = za1d * cSim.ra_inv;
            double cosphi                       = sqrt(1.0 - sinphi * sinphi);
            double sinpsi                       = (zb1d - zc1d) / (cSim.rc2 * cosphi);
            double cospsi                       = sqrt(1.0 - sinpsi * sinpsi);
 
            double ya2d                         =  cSim.ra * cosphi;
            double xb2d                         = -cSim.rc * cospsi;
            double yb2d                         = -cSim.rb * cosphi - cSim.rc * sinpsi * sinphi;
            double yc2d                         = -cSim.rb * cosphi + cSim.rc * sinpsi * sinphi;
            xb2d                                = -0.5 * sqrt(cSim.hhhh - (yb2d - yc2d) * (yb2d - yc2d) - (zb1d - zc1d) * (zb1d - zc1d));

            // Step3  al,be,ga:

            double alpa                         = (xb2d * (xb0d-xc0d) + yb0d * yb2d + yc0d * yc2d);
            double beta                         = (xb2d * (yc0d-yb0d) + xb0d * yb2d + xc0d * yc2d);
            double gama                         = xb0d * yb1d - xb1d * yb0d + xc0d * yc1d - xc1d * yc0d;

            double al2be2                       =  alpa * alpa + beta * beta;
            double sinthe                       = (alpa * gama - beta * sqrt(al2be2 - gama * gama)) / al2be2;

            // Step4  A3_prime:

            double costhe                       =  sqrt(1.0 - sinthe * sinthe);
            double xa3d                         = -ya2d * sinthe;
            double ya3d                         =  ya2d * costhe;
            double za3d                         =  za1d;
            double xb3d                         =  xb2d * costhe - yb2d * sinthe;
            double yb3d                         =  xb2d * sinthe + yb2d * costhe;
            double zb3d                         =  zb1d;
            double xc3d                         = -xb2d * costhe - yc2d * sinthe;
            double yc3d                         = -xb2d * sinthe + yc2d * costhe;
            double zc3d                         =  zc1d;

            // Step5  A3:
            cSim.pAtomX[shakeID.x]              = psA->xpl + trns11 * xa3d + trns12 * ya3d + trns13 * za3d;
            cSim.pAtomY[shakeID.x]              = psA->ypl + trns21 * xa3d + trns22 * ya3d + trns23 * za3d;
            cSim.pAtomZ[shakeID.x]              = psA->zpl + trns31 * xa3d + trns32 * ya3d + trns33 * za3d;
            cSim.pAtomX[shakeID.y]              = psA->xpl + trns11 * xb3d + trns12 * yb3d + trns13 * zb3d;
            cSim.pAtomY[shakeID.y]              = psA->ypl + trns21 * xb3d + trns22 * yb3d + trns23 * zb3d;
            cSim.pAtomZ[shakeID.y]              = psA->zpl + trns31 * xb3d + trns32 * yb3d + trns33 * zb3d;
            cSim.pAtomX[shakeID.z]              = psA->xpl + trns11 * xc3d + trns12 * yc3d + trns13 * zc3d;
            cSim.pAtomY[shakeID.z]              = psA->ypl + trns21 * xc3d + trns22 * yc3d + trns23 * zc3d;
            cSim.pAtomZ[shakeID.z]              = psA->zpl + trns31 * xc3d + trns32 * yc3d + trns33 * zc3d;

            pos                                += gridDim.x * blockDim.x;                                     
        }
    }

    if (cSim.slowShakeConstraints > 0)
    {    
        while (pos < cSim.fastShakeOffset)
        {
            pos                                += gridDim.x * blockDim.x;
        }
        pos                                    -= cSim.fastShakeOffset;

        while (pos < cSim.slowShakeConstraints)
        {
            int shakeID1;
		    int4 shakeID2;
            double toler;

            // Read SHAKE network data
            Atom* psA                           = &sA[threadIdx.x];
            shakeID1                            = cSim.pSlowShakeID1[pos];
            shakeID2                            = cSim.pSlowShakeID2[pos];
            double2 shakeParm                   = cSim.pSlowShakeParm[pos];
        
            // Read SHAKE network components
#ifdef use_SPSP        
            double xi                           = tex1Dfetch(texref, shakeID1);
            double yi                           = tex1Dfetch(texref, shakeID1 + cSim.stride);
            double zi                           = tex1Dfetch(texref, shakeID1 + cSim.stride2);
            double xij                          = tex1Dfetch(texref, shakeID2.x);
            double yij                          = tex1Dfetch(texref, shakeID2.x + cSim.stride);
            double zij                          = tex1Dfetch(texref, shakeID2.x + cSim.stride2);
            double xik                          = tex1Dfetch(texref, shakeID2.y);
            double yik                          = tex1Dfetch(texref, shakeID2.y + cSim.stride);
            double zik                          = tex1Dfetch(texref, shakeID2.y + cSim.stride2); 
            psA->xil                            = tex1Dfetch(texref, shakeID2.z);
            psA->yil                            = tex1Dfetch(texref, shakeID2.z + cSim.stride);
            psA->zil                            = tex1Dfetch(texref, shakeID2.z + cSim.stride2); 
            double xim                          = tex1Dfetch(texref, shakeID2.w);
            double yim                          = tex1Dfetch(texref, shakeID2.w + cSim.stride);
            double zim                          = tex1Dfetch(texref, shakeID2.w + cSim.stride2); 
#elif defined(NODPTEXTURE) && (__CUDA_ARCH__ < 200)
            double xi                           = cSim.pForceX[shakeID1];
            double yi                           = cSim.pForceY[shakeID1];
            double zi                           = cSim.pForceZ[shakeID1];
            double xij                          = cSim.pForceX[shakeID2.x];
            double yij                          = cSim.pForceY[shakeID2.x];
            double zij                          = cSim.pForceZ[shakeID2.x];
            double xik                          = cSim.pForceX[shakeID2.y];
            double yik                          = cSim.pForceY[shakeID2.y];
            double zik                          = cSim.pForceZ[shakeID2.y];  
            psA->xil                            = cSim.pForceX[shakeID2.z];
            psA->yil                            = cSim.pForceY[shakeID2.z];
            psA->zil                            = cSim.pForceZ[shakeID2.z]; 
            double xim                          = cSim.pForceX[shakeID2.w];
            double yim                          = cSim.pForceY[shakeID2.w];
            double zim                          = cSim.pForceZ[shakeID2.w];       
#else        
            int2 ixi                            = tex1Dfetch(texref, shakeID1);
            int2 iyi                            = tex1Dfetch(texref, shakeID1 + cSim.stride);
            int2 izi                            = tex1Dfetch(texref, shakeID1 + cSim.stride2);
            int2 ixij                           = tex1Dfetch(texref, shakeID2.x);
            int2 iyij                           = tex1Dfetch(texref, shakeID2.x + cSim.stride);
            int2 izij                           = tex1Dfetch(texref, shakeID2.x + cSim.stride2);
            int2 ixik                           = tex1Dfetch(texref, shakeID2.y);
            int2 iyik                           = tex1Dfetch(texref, shakeID2.y + cSim.stride);
            int2 izik                           = tex1Dfetch(texref, shakeID2.y + cSim.stride2);
            int2 ixil                           = tex1Dfetch(texref, shakeID2.z);
            int2 iyil                           = tex1Dfetch(texref, shakeID2.z + cSim.stride);
            int2 izil                           = tex1Dfetch(texref, shakeID2.z + cSim.stride2); 
            int2 ixim                           = tex1Dfetch(texref, shakeID2.w);
            int2 iyim                           = tex1Dfetch(texref, shakeID2.w + cSim.stride);
            int2 izim                           = tex1Dfetch(texref, shakeID2.w + cSim.stride2); 
#endif 
            double xpi                          = cSim.pAtomX[shakeID1];
            double ypi                          = cSim.pAtomY[shakeID1];
            double zpi                          = cSim.pAtomZ[shakeID1];
            double xpj                          = cSim.pAtomX[shakeID2.x];
            double ypj                          = cSim.pAtomY[shakeID2.x];
            double zpj                          = cSim.pAtomZ[shakeID2.x];
            double xpk                          = cSim.pAtomX[shakeID2.y];
            double ypk                          = cSim.pAtomY[shakeID2.y];
            double zpk                          = cSim.pAtomZ[shakeID2.y];
            psA->xpl                            = cSim.pAtomX[shakeID2.z];
            psA->ypl                            = cSim.pAtomY[shakeID2.z];
            psA->zpl                            = cSim.pAtomZ[shakeID2.z];
            double xpm                          = cSim.pAtomX[shakeID2.w];
            double ypm                          = cSim.pAtomY[shakeID2.w];
            double zpm                          = cSim.pAtomZ[shakeID2.w];
#if !defined(use_SPSP) && (!defined(NODPTEXTURE) || (__CUDA_ARCH__ >= 200))
            double xi                           = __hiloint2double(ixi.y, ixi.x);
            double yi                           = __hiloint2double(iyi.y, iyi.x);
            double zi                           = __hiloint2double(izi.y, izi.x);
            double xij                          = __hiloint2double(ixij.y, ixij.x);
            double yij                          = __hiloint2double(iyij.y, iyij.x);
            double zij                          = __hiloint2double(izij.y, izij.x);
            double xik                          = __hiloint2double(ixik.y, ixik.x);
            double yik                          = __hiloint2double(iyik.y, iyik.x);
            double zik                          = __hiloint2double(izik.y, izik.x);   
            psA->xil                            = __hiloint2double(ixil.y, ixil.x);
            psA->yil                            = __hiloint2double(iyil.y, iyil.x);
            psA->zil                            = __hiloint2double(izil.y, izil.x);   
            double xim                          = __hiloint2double(ixim.y, ixim.x);
            double yim                          = __hiloint2double(iyim.y, iyim.x);
            double zim                          = __hiloint2double(izim.y, izim.x);   
#endif                       
            psA->invMassI                       = shakeParm.x;
            toler                               = shakeParm.y;
            
            // Calculate unchanging quantities
            xij                                 = xi - xij;
            yij                                 = yi - yij;
            zij                                 = zi - zij;
            xik                                 = xi - xik;
            yik                                 = yi - yik;
            zik                                 = zi - zik; 
            psA->xil                            = xi - psA->xil;
            psA->yil                            = yi - psA->yil;
            psA->zil                            = zi - psA->zil;
            xim                                 = xi - xim;
            yim                                 = yi - yim;
            zim                                 = zi - zim;               
           
            bool done                           = false;
            for (int i = 0; i < 3000; i++)
            {
                done = true;
                
                // Calculate nominal distance squared
                double xpxx                     = xpi - xpj;
                double ypxx                     = ypi - ypj;
                double zpxx                     = zpi - zpj;
                double rpxx2                    = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
          
                // Apply correction to first hydrogen
                double diff                     = toler - rpxx2;
                if (abs(diff) >= toler * cSim.tol)
                {
                    done                        = false;
                   
                    // Shake resetting of coordinate is done here
                    double rrpr                 = xij * xpxx + yij * ypxx + zij * zpxx;     
                    if (rrpr >= toler * 1.0e-06)
                    {
                    
                        double acor             = diff / (rrpr * (double)2.0 * (psA->invMassI + cSim.invMassH));
                        double h                = xij * acor;
                        xpi                    += h * psA->invMassI;
                        xpj                    -= h * cSim.invMassH;
                        h                       = yij * acor;
                        ypi                    += h * psA->invMassI;
                        ypj                    -= h * cSim.invMassH;
                        h                       = zij * acor;
                        zpi                    += h * psA->invMassI;
                        zpj                    -= h * cSim.invMassH;             
                    }
                }
          
     
                xpxx                            = xpi - xpk;
                ypxx                            = ypi - ypk;
                zpxx                            = zpi - zpk;
                rpxx2                           = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
          
                // Apply correction to second hydrogen
                diff                            = toler - rpxx2;
                if (abs(diff) >= toler * cSim.tol)
                {
                    done                        = false;
               
                    // Shake resetting of coordinate is done here
                    double rrpr                 = xik * xpxx + yik * ypxx + zik * zpxx;     
                    if (rrpr >= toler * 1.0e-06)
                    {
                    
                        double acor             = diff / (rrpr * 2.0 * (psA->invMassI + cSim.invMassH));
                        double h                = xik * acor;
                        xpi                    += h * psA->invMassI;
                        xpk                    -= h * cSim.invMassH;
                        h                       = yik * acor;
                        ypi                    += h * psA->invMassI;
                        ypk                    -= h * cSim.invMassH;
                        h                       = zik * acor;
                        zpi                    += h * psA->invMassI;
                        zpk                    -= h * cSim.invMassH;             
                    }
                }
                
 
                xpxx                            = xpi - psA->xpl;
                ypxx                            = ypi - psA->ypl;
                zpxx                            = zpi - psA->zpl;        
                rpxx2                           = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
          
                // Apply correction to third hydrogen
                diff                            = toler - rpxx2;
                if (abs(diff) >= toler * cSim.tol)
                {
                    done                        = false;
                  
                    // Shake resetting of coordinate is done here
                    double rrpr                 = psA->xil * xpxx + psA->yil * ypxx + psA->zil * zpxx;     
                    if (rrpr >= toler * 1.0e-06)
                    {             
                        double acor             = diff / (rrpr * (double)2.0 * (psA->invMassI + cSim.invMassH));
                        double h                = psA->xil * acor;
                        xpi                    += h * psA->invMassI;
                        psA->xpl               -= h * cSim.invMassH;
                        h                       = psA->yil * acor;
                        ypi                    += h * psA->invMassI;
                        psA->ypl               -= h * cSim.invMassH;
                        h                       = psA->zil * acor;
                        zpi                    += h * psA->invMassI;
                        psA->zpl               -= h * cSim.invMassH;             
                    }
                }

                xpxx                            = xpi - xpm;
                ypxx                            = ypi - ypm;
                zpxx                            = zpi - zpm;        
                rpxx2                           = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;

                // Apply correction to third hydrogen
                diff                            = toler - rpxx2;
                if (abs(diff) >= toler * cSim.tol)
                {
                    done                        = false;
                  
                    // Shake resetting of coordinate is done here
                    double rrpr                 = xim * xpxx + yim * ypxx + zim * zpxx;     
                    if (rrpr >= toler * 1.0e-06)
                    {             
                        double acor             = diff / (rrpr * (double)2.0 * (psA->invMassI + cSim.invMassH));
                        double h                = xim * acor;
                        xpi                    += h * psA->invMassI;
                        xpm                    -= h * cSim.invMassH;
                        h                       = yim * acor;
                        ypi                    += h * psA->invMassI;
                        ypm                    -= h * cSim.invMassH;
                        h                       = zim * acor;
                        zpi                    += h * psA->invMassI;
                        zpm                    -= h * cSim.invMassH;             
                    }
                }

                
                
                // Check for convergence
                if (done)
                    break;
            }
          
            // Write out results if converged, but there's no really good
            // way to indicate failure so we'll let the simulation heading
            // off to Neptune do that for us.  Wish there were a better way,
            // but until the CPU needs something from the GPU, those are the
            // the breaks.  I guess, technically, we could just set a flag to NOP
            // the simulation from here and then carry that result through upon
            // the next ntpr, ntwc, or ntwx update, but I leave that up to you 
            // guys to implement that (or not). 
            if (done)
            {
                cSim.pAtomX[shakeID1]           = xpi;
                cSim.pAtomY[shakeID1]           = ypi;
                cSim.pAtomZ[shakeID1]           = zpi;
                PMEFloat2 xyi                   = {xpi, ypi};
                cSim.pAtomXYSP[shakeID1]        = xyi;
                cSim.pAtomZSP[shakeID1]         = zpi;
                cSim.pAtomX[shakeID2.x]         = xpj;
                cSim.pAtomY[shakeID2.x]         = ypj;
                cSim.pAtomZ[shakeID2.x]         = zpj;
                PMEFloat2 xyj                   = {xpj, ypj};
                cSim.pAtomXYSP[shakeID2.x]      = xyj;
                cSim.pAtomZSP[shakeID2.x]       = zpj;
                cSim.pAtomX[shakeID2.y]         = xpk;
                cSim.pAtomY[shakeID2.y]         = ypk;
                cSim.pAtomZ[shakeID2.y]         = zpk;
                PMEFloat2 xyk                   = {xpk, ypk};
                cSim.pAtomXYSP[shakeID2.y]      = xyk;
                cSim.pAtomZSP[shakeID2.y]       = zpk;
                cSim.pAtomX[shakeID2.z]         = psA->xpl;
                cSim.pAtomY[shakeID2.z]         = psA->ypl;
                cSim.pAtomZ[shakeID2.z]         = psA->zpl;
                PMEFloat2 xyl                   = {psA->xpl, psA->ypl};
                cSim.pAtomXYSP[shakeID2.z]      = xyl;
                cSim.pAtomZSP[shakeID2.z]       = psA->zpl;
                cSim.pAtomX[shakeID2.w]         = xpm;
                cSim.pAtomY[shakeID2.w]         = ypm;
                cSim.pAtomZ[shakeID2.w]         = zpm;
                PMEFloat2 xym                   = {xpm, ypm};
                cSim.pAtomXYSP[shakeID2.w]      = xym;
                cSim.pAtomZSP[shakeID2.w]       = zpm;
            }        

    
            pos                                += gridDim.x * blockDim.x;                                 
        }
    }


}

#if (__CUDA_ARCH__ >= 200)
struct PMEAtom 
{
    double invMassI;
    double xpl;
    double ypl;
    double zpl;
    double xil;
    double yil;
    double zil;
    double toler;
    int4 shakeID;
    double dummy1;
    double dummy2;
};


struct PMEFastAtom
{
    double xcom;
    double ycom;
    double zcom;
    double trns11;
    double trns12;
    double trns13;
    double trns21;
    double trns22;
    double trns23;
    double trns31;
    double trns32;
    double trns33;
};
#else
struct PMEAtom 
{
    double invMassI;
    double xpl;
    double ypl;
    double zpl;
    double xil;
    double yil;
    double zil;
    double dummy1;
    double dummy2;
    double dummy3;
};


struct PMEFastAtom
{
    double xcom;
    double ycom;
    double zcom;
    double trns11;
    double trns12;
    double trns13;
    double trns21;
    double trns22;
    double trns23;
    double trns31;
};
#endif

__global__ void 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_SHAKE_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(SM_13_SHAKE_THREADS_PER_BLOCK, 1)
#endif
kPMEShake_kernel()
{
#if (__CUDA_ARCH__ >= 200)
__shared__ PMEAtom sA[SM_2X_SHAKE_THREADS_PER_BLOCK];
#else
__shared__ PMEAtom sA[SM_13_SHAKE_THREADS_PER_BLOCK];
#endif

    unsigned int pos                            = blockIdx.x * blockDim.x + threadIdx.x;   
    while (pos < cSim.shakeConstraints)
    {
        PMEAtom* psA                            = &sA[threadIdx.x];
#if (__CUDA_ARCH__ >= 200)
#define TOLER psA->toler
#define SHAKEID psA->shakeID
#define SHAKEIDX psA->shakeID.x
#define SHAKEIDY psA->shakeID.y
#define SHAKEIDZ psA->shakeID.z
#define SHAKEIDW psA->shakeID.w
#else
		int4 shakeID;
        double toler;
#define TOLER toler
#define SHAKEID shakeID
#define SHAKEIDX shakeID.x
#define SHAKEIDY shakeID.y
#define SHAKEIDZ shakeID.z
#define SHAKEIDW shakeID.w
#endif
        // Read SHAKE network data
        SHAKEID                            		= cSim.pImageShakeID[pos];
        double2 shakeParm                       = cSim.pShakeParm[pos];
        
        // Read SHAKE network components
#ifdef use_SPSP        
        double xi                               = tex1Dfetch(texref, SHAKEIDX);
        double yi                               = tex1Dfetch(texref, SHAKEIDX + cSim.stride);
        double zi                               = tex1Dfetch(texref, SHAKEIDX + cSim.stride2);
        double xij                              = tex1Dfetch(texref, SHAKEIDY);
        double yij                              = tex1Dfetch(texref, SHAKEIDY + cSim.stride);
        double zij                              = tex1Dfetch(texref, SHAKEIDY + cSim.stride2);
#elif defined(NODPTEXTURE) && (__CUDA_ARCH__ < 200)
        double xi                               = cSim.pForceX[SHAKEIDX];
        double yi                               = cSim.pForceY[SHAKEIDX];
        double zi                               = cSim.pForceZ[SHAKEIDX];
        double xij                              = cSim.pForceX[SHAKEIDY];
        double yij                              = cSim.pForceY[SHAKEIDY];
        double zij                              = cSim.pForceZ[SHAKEIDY];
#else        
        int2 ixi                                = tex1Dfetch(texref, SHAKEIDX);
        int2 iyi                                = tex1Dfetch(texref, SHAKEIDX + cSim.stride);
        int2 izi                                = tex1Dfetch(texref, SHAKEIDX + cSim.stride2);
        int2 ixij                               = tex1Dfetch(texref, SHAKEIDY);
        int2 iyij                               = tex1Dfetch(texref, SHAKEIDY + cSim.stride);
        int2 izij                               = tex1Dfetch(texref, SHAKEIDY + cSim.stride2);
#endif 
        double xpi                              = cSim.pImageX[SHAKEIDX];
        double ypi                              = cSim.pImageY[SHAKEIDX];
        double zpi                              = cSim.pImageZ[SHAKEIDX];
        double xpj                              = cSim.pImageX[SHAKEIDY];
        double ypj                              = cSim.pImageY[SHAKEIDY];
        double zpj                              = cSim.pImageZ[SHAKEIDY];
#if !defined(use_SPSP) && (!defined(NODPTEXTURE) || (__CUDA_ARCH__ >= 200))
        double xi                               = __hiloint2double(ixi.y, ixi.x);
        double yi                               = __hiloint2double(iyi.y, iyi.x);
        double zi                               = __hiloint2double(izi.y, izi.x);
        double xij                              = __hiloint2double(ixij.y, ixij.x);
        double yij                              = __hiloint2double(iyij.y, iyij.x);
        double zij                              = __hiloint2double(izij.y, izij.x); 
#endif                       
        psA->invMassI                           = shakeParm.x;
        TOLER                                   = shakeParm.y;
        
        
        // Optionally read 2nd hydrogen
        double xpk, ypk, zpk, xik, yik, zik;
        if (SHAKEIDZ != -1)
        {
#ifdef use_SPSP
            xik                                 = tex1Dfetch(texref, SHAKEIDZ);
            yik                                 = tex1Dfetch(texref, SHAKEIDZ + cSim.stride);
            zik                                 = tex1Dfetch(texref, SHAKEIDZ + cSim.stride2); 
#elif defined(NODPTEXTURE) && (__CUDA_ARCH__ < 200)
            xik                                 = cSim.pForceX[SHAKEIDZ];
            yik                                 = cSim.pForceY[SHAKEIDZ];
            zik                                 = cSim.pForceZ[SHAKEIDZ];    
#else
            int2 ixik                           = tex1Dfetch(texref, SHAKEIDZ);
            int2 iyik                           = tex1Dfetch(texref, SHAKEIDZ + cSim.stride);
            int2 izik                           = tex1Dfetch(texref, SHAKEIDZ + cSim.stride2);
#endif 
            xpk                                 = cSim.pImageX[SHAKEIDZ];
            ypk                                 = cSim.pImageY[SHAKEIDZ];
            zpk                                 = cSim.pImageZ[SHAKEIDZ];
#if !defined(use_SPSP) && (!defined(NODPTEXTURE) || (__CUDA_ARCH__ >= 200))
            xik                                 = __hiloint2double(ixik.y, ixik.x);
            yik                                 = __hiloint2double(iyik.y, iyik.x);
            zik                                 = __hiloint2double(izik.y, izik.x);  
#endif                  
        }
        
        // Optionally read 3rd hydrogen into shared memory
        if (SHAKEIDW != -1)
        {
#ifdef use_SPSP
            psA->xil                            = tex1Dfetch(texref, SHAKEIDW);
            psA->yil                            = tex1Dfetch(texref, SHAKEIDW + cSim.stride);
            psA->zil                            = tex1Dfetch(texref, SHAKEIDW + cSim.stride2); 
#elif defined(NODPTEXTURE) && (__CUDA_ARCH__ < 200)
            psA->xil                            = cSim.pForceX[SHAKEIDW];
            psA->yil                            = cSim.pForceY[SHAKEIDW];
            psA->zil                            = cSim.pForceZ[SHAKEIDW];      
#else            
            int2 ixil                           = tex1Dfetch(texref, SHAKEIDW);
            int2 iyil                           = tex1Dfetch(texref, SHAKEIDW + cSim.stride);
            int2 izil                           = tex1Dfetch(texref, SHAKEIDW + cSim.stride2); 
#endif          
            psA->xpl                            = cSim.pImageX[SHAKEIDW];
            psA->ypl                            = cSim.pImageY[SHAKEIDW];
            psA->zpl                            = cSim.pImageZ[SHAKEIDW];
#if !defined(use_SPSP) && (!defined(NODPTEXTURE) || (__CUDA_ARCH__ >= 200))
            psA->xil                            = __hiloint2double(ixil.y, ixil.x);
            psA->yil                            = __hiloint2double(iyil.y, iyil.x);
            psA->zil                            = __hiloint2double(izil.y, izil.x);   
#endif          
        }
        
        // Calculate unchanging quantities
        xij                                     = xi - xij;
        yij                                     = yi - yij;
        zij                                     = zi - zij;
        
        if (SHAKEIDZ != -1)
        {
            xik                                 = xi - xik;
            yik                                 = yi - yik;
            zik                                 = zi - zik; 
        }        
         
        if (SHAKEIDW != -1)
        {
            psA->xil                            = xi - psA->xil;
            psA->yil                            = yi - psA->yil;
            psA->zil                            = zi - psA->zil;
        }      
       
        bool done                               = false;
        for (int i = 0; i < 3000; i++)
        {
            done = true;
            
            // Calculate nominal distance squared
            double xpxx                         = xpi - xpj;
            double ypxx                         = ypi - ypj;
            double zpxx                         = zpi - zpj;
            double rpxx2                        = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
      
            // Apply correction
            double diff                         = TOLER - rpxx2;
            if (abs(diff) >= TOLER * cSim.tol)
            {
                done                            = false;
               
                // Shake resetting of coordinate is done here
                double rrpr                     = xij * xpxx + yij * ypxx + zij * zpxx;     
                if (rrpr >= TOLER * 1.0e-06)
                {
                
                    double acor                 = diff / (rrpr * (double)2.0 * (psA->invMassI + cSim.invMassH));
                    double h                    = xij * acor;
                    xpi                        += h * psA->invMassI;
                    xpj                        -= h * cSim.invMassH;
                    h                           = yij * acor;
                    ypi                        += h * psA->invMassI;
                    ypj                        -= h * cSim.invMassH;
                    h                           = zij * acor;
                    zpi                        += h * psA->invMassI;
                    zpj                        -= h * cSim.invMassH;             
                }
            }
      
            // Second bond if present
            if (SHAKEIDZ != -1)
            {
                xpxx                            = xpi - xpk;
                ypxx                            = ypi - ypk;
                zpxx                            = zpi - zpk;
                rpxx2                           = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
      
                // Apply correction
                diff                            = TOLER - rpxx2;
                if (abs(diff) >= TOLER * cSim.tol)
                {
                    done                        = false;
               
                    // Shake resetting of coordinate is done here
                    double rrpr                 = xik * xpxx + yik * ypxx + zik * zpxx;     
                    if (rrpr >= TOLER * 1.0e-06)
                    {
                
                        double acor             = diff / (rrpr * 2.0 * (psA->invMassI + cSim.invMassH));
                        double h                = xik * acor;
                        xpi                    += h * psA->invMassI;
                        xpk                    -= h * cSim.invMassH;
                        h                       = yik * acor;
                        ypi                    += h * psA->invMassI;
                        ypk                    -= h * cSim.invMassH;
                        h                       = zik * acor;
                        zpi                    += h * psA->invMassI;
                        zpk                    -= h * cSim.invMassH;             
                    }
                }
            }
            
            // Third bond if present
            if (SHAKEIDW != -1)
            {
                xpxx                            = xpi - psA->xpl;
                ypxx                            = ypi - psA->ypl;
                zpxx                            = zpi - psA->zpl;        
                rpxx2                           = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
      
                // Apply correction
                diff                            = TOLER - rpxx2;
                if (abs(diff) >= TOLER * cSim.tol)
                {
                    done                        = false;
               
                    // Shake resetting of coordinate is done here
                    double rrpr                 = psA->xil * xpxx + psA->yil * ypxx + psA->zil * zpxx;     
                    if (rrpr >= TOLER * 1.0e-06)
                    {
                
                        double acor             = diff / (rrpr * (double)2.0 * (psA->invMassI + cSim.invMassH));
                        double h                = psA->xil * acor;
                        xpi                    += h * psA->invMassI;
                        psA->xpl               -= h * cSim.invMassH;
                        h                       = psA->yil * acor;
                        ypi                    += h * psA->invMassI;
                        psA->ypl               -= h * cSim.invMassH;
                        h                       = psA->zil * acor;
                        zpi                    += h * psA->invMassI;
                        psA->zpl               -= h * cSim.invMassH;             
                    }
                }
            }
            
            
            // Check for convergence
            if (done)
                break;
        }
      
        // Write out results if converged, but there's no really good
        // way to indicate failure so we'll let the simulation heading
        // off to Neptune do that for us.  Wish there were a better way,
        // but until the CPU needs something from the GPU, those are the
        // the breaks.  I guess, technically, we could just set a flag to NOP
        // the simulation from here and then carry that result through upon
        // the next ntpr, ntwc, or ntwx update, but I leave that up to you 
        // guys to implement that (or not). 
        if (done)
        {
            cSim.pImageX[SHAKEIDX]              = xpi;
            cSim.pImageY[SHAKEIDX]              = ypi;
            cSim.pImageZ[SHAKEIDX]              = zpi;
            cSim.pImageX[SHAKEIDY]              = xpj;
            cSim.pImageY[SHAKEIDY]              = ypj;
            cSim.pImageZ[SHAKEIDY]              = zpj;

            if (SHAKEIDZ != -1)
            {
                cSim.pImageX[SHAKEIDZ]          = xpk;
                cSim.pImageY[SHAKEIDZ]          = ypk;
                cSim.pImageZ[SHAKEIDZ]          = zpk;
            }

            if (SHAKEIDW != -1)
            {
                cSim.pImageX[SHAKEIDW]          = psA->xpl;
                cSim.pImageY[SHAKEIDW]          = psA->ypl;
                cSim.pImageZ[SHAKEIDW]          = psA->zpl;
            }
        }        
        pos                                    += gridDim.x * blockDim.x;
    }

#undef TOLER
#undef SHAKEID
#undef SHAKEIDX
#undef SHAKEIDY
#undef SHAKEIDZ
#undef SHAKEIDW
#if (__CUDA_ARCH__ >= 200)
#define TRNS32 psA->trns32
#define TRNS33 psA->trns33
#else
#define TRNS32 trns32
#define TRNS33 trns33
#endif
    if (cSim.fastShakeConstraints > 0)
    {    
        while (pos < cSim.shakeOffset)
        {
            pos                                += gridDim.x * blockDim.x;
        }
        pos                                    -= cSim.shakeOffset;

        while (pos < cSim.fastShakeConstraints)
        {
            PMEFastAtom* psA                    = (PMEFastAtom*)&sA[threadIdx.x];
            
            // Read atom data
            int4 shakeID                        = cSim.pImageFastShakeID[pos];
#ifdef use_SPSP
            double x1                           = tex1Dfetch(texref, shakeID.x);
            double y1                           = tex1Dfetch(texref, shakeID.x + cSim.stride);
            double z1                           = tex1Dfetch(texref, shakeID.x + cSim.stride2);
            double x2                           = tex1Dfetch(texref, shakeID.y);
            double y2                           = tex1Dfetch(texref, shakeID.y + cSim.stride);
            double z2                           = tex1Dfetch(texref, shakeID.y + cSim.stride2);
            double x3                           = tex1Dfetch(texref, shakeID.z);
            double y3                           = tex1Dfetch(texref, shakeID.z + cSim.stride);
            double z3                           = tex1Dfetch(texref, shakeID.z + cSim.stride2);  
#elif defined(NODPTEXTURE) && (__CUDA_ARCH__ < 200)
            double x1                           = cSim.pForceX[shakeID.x];
            double y1                           = cSim.pForceY[shakeID.x];
            double z1                           = cSim.pForceZ[shakeID.x]; 
            double x2                           = cSim.pForceX[shakeID.y];
            double y2                           = cSim.pForceY[shakeID.y];
            double z2                           = cSim.pForceZ[shakeID.y]; 
            double x3                           = cSim.pForceX[shakeID.z];
            double y3                           = cSim.pForceY[shakeID.z];
            double z3                           = cSim.pForceZ[shakeID.z]; 
#else           
            int2 ix1                            = tex1Dfetch(texref, shakeID.x);
            int2 iy1                            = tex1Dfetch(texref, shakeID.x + cSim.stride);
            int2 iz1                            = tex1Dfetch(texref, shakeID.x + cSim.stride2);
            int2 ix2                            = tex1Dfetch(texref, shakeID.y);
            int2 iy2                            = tex1Dfetch(texref, shakeID.y + cSim.stride);
            int2 iz2                            = tex1Dfetch(texref, shakeID.y + cSim.stride2);
            int2 ix3                            = tex1Dfetch(texref, shakeID.z);
            int2 iy3                            = tex1Dfetch(texref, shakeID.z + cSim.stride);
            int2 iz3                            = tex1Dfetch(texref, shakeID.z + cSim.stride2);  
#endif
            double xp1                          = cSim.pImageX[shakeID.x];
            double yp1                          = cSim.pImageY[shakeID.x];
            double zp1                          = cSim.pImageZ[shakeID.x];
            double xp2                          = cSim.pImageX[shakeID.y];
            double yp2                          = cSim.pImageY[shakeID.y];
            double zp2                          = cSim.pImageZ[shakeID.y];
            double xp3                          = cSim.pImageX[shakeID.z];
            double yp3                          = cSim.pImageY[shakeID.z];
            double zp3                          = cSim.pImageZ[shakeID.z];
#if !defined(use_SPSP) && (!defined(NODPTEXTURE) || (__CUDA_ARCH__ >= 200))
            double x1                           = __hiloint2double(ix1.y, ix1.x);
            double y1                           = __hiloint2double(iy1.y, iy1.x);
            double z1                           = __hiloint2double(iz1.y, iz1.x);
            double x2                           = __hiloint2double(ix2.y, ix2.x);
            double y2                           = __hiloint2double(iy2.y, iy2.x);
            double z2                           = __hiloint2double(iz2.y, iz2.x);
            double x3                           = __hiloint2double(ix3.y, ix3.x);
            double y3                           = __hiloint2double(iy3.y, iy3.x);
            double z3                           = __hiloint2double(iz3.y, iz3.x);
#endif

            // Step1  A1_prime:
            double xb0                          = x2 - x1;
            double yb0                          = y2 - y1;
            double zb0                          = z2 - z1;
            double xc0                          = x3 - x1;
            double yc0                          = y3 - y1;
            double zc0                          = z3 - z1;
            psA->xcom                           = xp1 * cSim.wo_div_wohh + (xp2 + xp3) * cSim.wh_div_wohh;
            psA->ycom                           = yp1 * cSim.wo_div_wohh + (yp2 + yp3) * cSim.wh_div_wohh;
            psA->zcom                           = zp1 * cSim.wo_div_wohh + (zp2 + zp3) * cSim.wh_div_wohh;

            double xa1                          = xp1 - psA->xcom;
            double ya1                          = yp1 - psA->ycom;
            double za1                          = zp1 - psA->zcom;
            double xb1                          = xp2 - psA->xcom;
            double yb1                          = yp2 - psA->ycom;
            double zb1                          = zp2 - psA->zcom;
            double xc1                          = xp3 - psA->xcom;
            double yc1                          = yp3 - psA->ycom;
            double zc1                          = zp3 - psA->zcom;
            double xakszd                       = yb0 * zc0 - zb0 * yc0;
            double yakszd                       = zb0 * xc0 - xb0 * zc0;
            double zakszd                       = xb0 * yc0 - yb0 * xc0;
            double xaksxd                       = ya1 * zakszd - za1 * yakszd;
            double yaksxd                       = za1 * xakszd - xa1 * zakszd;
            double zaksxd                       = xa1 * yakszd - ya1 * xakszd;
            double xaksyd                       = yakszd * zaksxd - zakszd * yaksxd;
            double yaksyd                       = zakszd * xaksxd - xakszd * zaksxd;
            double zaksyd                       = xakszd * yaksxd - yakszd * xaksxd;

            double axlng_inv                    = rsqrt(xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd);
            double aylng_inv                    = rsqrt(xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd);
            double azlng_inv                    = rsqrt(xakszd * xakszd + yakszd * yakszd + zakszd * zakszd);

            psA->trns11                         = xaksxd * axlng_inv;
            psA->trns21                         = yaksxd * axlng_inv;
            psA->trns31                         = zaksxd * axlng_inv;
            psA->trns12                         = xaksyd * aylng_inv;
            psA->trns22                         = yaksyd * aylng_inv;
#if (__CUDA_ARCH__ < 200)
			double trns32;
#endif
            TRNS32                              = zaksyd * aylng_inv;
            psA->trns13                         = xakszd * azlng_inv;
            psA->trns23                         = yakszd * azlng_inv;
#if (__CUDA_ARCH__ < 200)
			double trns33;
#endif
            TRNS33                              = zakszd * azlng_inv;

            double xb0d                         = psA->trns11 * xb0 + psA->trns21 * yb0 + psA->trns31 * zb0;
            double yb0d                         = psA->trns12 * xb0 + psA->trns22 * yb0 + TRNS32      * zb0;
            double xc0d                         = psA->trns11 * xc0 + psA->trns21 * yc0 + psA->trns31 * zc0;
            double yc0d                         = psA->trns12 * xc0 + psA->trns22 * yc0 + TRNS32      * zc0;
            double za1d                         = psA->trns13 * xa1 + psA->trns23 * ya1 + TRNS33      * za1;
            double xb1d                         = psA->trns11 * xb1 + psA->trns21 * yb1 + psA->trns31 * zb1;
            double yb1d                         = psA->trns12 * xb1 + psA->trns22 * yb1 + TRNS32      * zb1;
            double zb1d                         = psA->trns13 * xb1 + psA->trns23 * yb1 + TRNS33      * zb1;
            double xc1d                         = psA->trns11 * xc1 + psA->trns21 * yc1 + psA->trns31 * zc1;
            double yc1d                         = psA->trns12 * xc1 + psA->trns22 * yc1 + TRNS32      * zc1;
            double zc1d                         = psA->trns13 * xc1 + psA->trns23 * yc1 + TRNS33      * zc1;

            // Step2  A2_prime:
            double sinphi                       = za1d * cSim.ra_inv;
            double cosphi                       = sqrt(1.0 - sinphi * sinphi);
            double sinpsi                       = (zb1d - zc1d) / (cSim.rc2 * cosphi);
            double cospsi                       = sqrt(1.0 - sinpsi * sinpsi);
 
            double ya2d                         =  cSim.ra * cosphi;
            double xb2d                         = -cSim.rc * cospsi;
            double yb2d                         = -cSim.rb * cosphi - cSim.rc * sinpsi * sinphi;
            double yc2d                         = -cSim.rb * cosphi + cSim.rc * sinpsi * sinphi;
            xb2d                                = -0.5 * sqrt(cSim.hhhh - (yb2d-yc2d) * (yb2d - yc2d) - (zb1d - zc1d) * (zb1d - zc1d));

            // Step3  al,be,ga:
            double alpa                         = (xb2d * (xb0d-xc0d) + yb0d * yb2d + yc0d * yc2d);
            double beta                         = (xb2d * (yc0d-yb0d) + xb0d * yb2d + xc0d * yc2d);
            double gama                         = xb0d * yb1d - xb1d * yb0d + xc0d * yc1d - xc1d * yc0d;

            double al2be2                       =  alpa * alpa + beta * beta;
            double sinthe                       = (alpa * gama - beta * sqrt(al2be2 - gama * gama)) / al2be2;

            // Step4  A3_prime:
            double costhe                       =  sqrt(1.0 - sinthe * sinthe);
            double xa3d                         = -ya2d * sinthe;
            double ya3d                         =  ya2d * costhe;
            double za3d                         =  za1d;
            double xb3d                         =  xb2d * costhe - yb2d * sinthe;
            double yb3d                         =  xb2d * sinthe + yb2d * costhe;
            double zb3d                         =  zb1d;
            double xc3d                         = -xb2d * costhe - yc2d * sinthe;
            double yc3d                         = -xb2d * sinthe + yc2d * costhe;
            double zc3d                         =  zc1d;

            // Step5  A3:
            cSim.pImageX[shakeID.x]             = psA->xcom + psA->trns11 * xa3d + psA->trns12 * ya3d + psA->trns13 * za3d;
            cSim.pImageY[shakeID.x]             = psA->ycom + psA->trns21 * xa3d + psA->trns22 * ya3d + psA->trns23 * za3d;
            cSim.pImageZ[shakeID.x]             = psA->zcom + psA->trns31 * xa3d + TRNS32      * ya3d + TRNS33      * za3d;
            cSim.pImageX[shakeID.y]             = psA->xcom + psA->trns11 * xb3d + psA->trns12 * yb3d + psA->trns13 * zb3d;
            cSim.pImageY[shakeID.y]             = psA->ycom + psA->trns21 * xb3d + psA->trns22 * yb3d + psA->trns23 * zb3d;
            cSim.pImageZ[shakeID.y]             = psA->zcom + psA->trns31 * xb3d + TRNS32      * yb3d + TRNS33      * zb3d;
            cSim.pImageX[shakeID.z]             = psA->xcom + psA->trns11 * xc3d + psA->trns12 * yc3d + psA->trns13 * zc3d;
            cSim.pImageY[shakeID.z]             = psA->ycom + psA->trns21 * xc3d + psA->trns22 * yc3d + psA->trns23 * zc3d;
            cSim.pImageZ[shakeID.z]             = psA->zcom + psA->trns31 * xc3d + TRNS32      * yc3d + TRNS33      * zc3d;
            pos                                += gridDim.x * blockDim.x;                                     
        }
    }
#undef TRNS32
#undef TRNS33

    if (cSim.slowShakeConstraints > 0)
    {    
        while (pos < cSim.fastShakeOffset)
        {
            pos                                += gridDim.x * blockDim.x;
        }
        pos                                    -= cSim.fastShakeOffset;

        while (pos < cSim.slowShakeConstraints)
        {

            int shakeID1;
#if (__CUDA_ARCH__ >= 200)
#define TOLER psA->toler
#define SHAKEID2 psA->shakeID
#define SHAKEID2X psA->shakeID.x
#define SHAKEID2Y psA->shakeID.y
#define SHAKEID2Z psA->shakeID.z
#define SHAKEID2W psA->shakeID.w
#else
		    int4 shakeID2;
            double toler;
#define TOLER toler
#define SHAKEID2 shakeID2
#define SHAKEID2X shakeID2.x
#define SHAKEID2Y shakeID2.y
#define SHAKEID2Z shakeID2.z
#define SHAKEID2W shakeID2.w
#endif
            // Read SHAKE network data
            PMEAtom* psA                        = &sA[threadIdx.x];
            shakeID1                            = cSim.pImageSlowShakeID1[pos];
            SHAKEID2                            = cSim.pImageSlowShakeID2[pos];
            double2 shakeParm                   = cSim.pSlowShakeParm[pos];
        
            // Read SHAKE network components
#ifdef use_SPSP        
            double xi                           = tex1Dfetch(texref, shakeID1);
            double yi                           = tex1Dfetch(texref, shakeID1 + cSim.stride);
            double zi                           = tex1Dfetch(texref, shakeID1 + cSim.stride2);
            double xij                          = tex1Dfetch(texref, SHAKEID2X);
            double yij                          = tex1Dfetch(texref, SHAKEID2X + cSim.stride);
            double zij                          = tex1Dfetch(texref, SHAKEID2X + cSim.stride2);
            double xik                          = tex1Dfetch(texref, SHAKEID2Y);
            double yik                          = tex1Dfetch(texref, SHAKEID2Y + cSim.stride);
            double zik                          = tex1Dfetch(texref, SHAKEID2Y + cSim.stride2); 
            psA->xil                            = tex1Dfetch(texref, SHAKEID2Z);
            psA->yil                            = tex1Dfetch(texref, SHAKEID2Z + cSim.stride);
            psA->zil                            = tex1Dfetch(texref, SHAKEID2Z + cSim.stride2); 
            double xim                          = tex1Dfetch(texref, SHAKEID2W);
            double yim                          = tex1Dfetch(texref, SHAKEID2W + cSim.stride);
            double zim                          = tex1Dfetch(texref, SHAKEID2W + cSim.stride2); 
#elif defined(NODPTEXTURE) && (__CUDA_ARCH__ < 200)
            double xi                           = cSim.pForceX[shakeID1];
            double yi                           = cSim.pForceY[shakeID1];
            double zi                           = cSim.pForceZ[shakeID1];
            double xij                          = cSim.pForceX[SHAKEID2X];
            double yij                          = cSim.pForceY[SHAKEID2X];
            double zij                          = cSim.pForceZ[SHAKEID2X];
            double xik                          = cSim.pForceX[SHAKEID2Y];
            double yik                          = cSim.pForceY[SHAKEID2Y];
            double zik                          = cSim.pForceZ[SHAKEID2Y];  
            psA->xil                            = cSim.pForceX[SHAKEID2Z];
            psA->yil                            = cSim.pForceY[SHAKEID2Z];
            psA->zil                            = cSim.pForceZ[SHAKEID2Z]; 
            double xim                          = cSim.pForceX[SHAKEID2W];
            double yim                          = cSim.pForceY[SHAKEID2W];
            double zim                          = cSim.pForceZ[SHAKEID2W];       
#else        
            int2 ixi                            = tex1Dfetch(texref, shakeID1);
            int2 iyi                            = tex1Dfetch(texref, shakeID1 + cSim.stride);
            int2 izi                            = tex1Dfetch(texref, shakeID1 + cSim.stride2);
            int2 ixij                           = tex1Dfetch(texref, SHAKEID2X);
            int2 iyij                           = tex1Dfetch(texref, SHAKEID2X + cSim.stride);
            int2 izij                           = tex1Dfetch(texref, SHAKEID2X + cSim.stride2);
            int2 ixik                           = tex1Dfetch(texref, SHAKEID2Y);
            int2 iyik                           = tex1Dfetch(texref, SHAKEID2Y + cSim.stride);
            int2 izik                           = tex1Dfetch(texref, SHAKEID2Y + cSim.stride2);
            int2 ixil                           = tex1Dfetch(texref, SHAKEID2Z);
            int2 iyil                           = tex1Dfetch(texref, SHAKEID2Z + cSim.stride);
            int2 izil                           = tex1Dfetch(texref, SHAKEID2Z + cSim.stride2); 
            int2 ixim                           = tex1Dfetch(texref, SHAKEID2W);
            int2 iyim                           = tex1Dfetch(texref, SHAKEID2W + cSim.stride);
            int2 izim                           = tex1Dfetch(texref, SHAKEID2W + cSim.stride2); 
#endif 
            double xpi                          = cSim.pImageX[shakeID1];
            double ypi                          = cSim.pImageY[shakeID1];
            double zpi                          = cSim.pImageZ[shakeID1];
            double xpj                          = cSim.pImageX[SHAKEID2X];
            double ypj                          = cSim.pImageY[SHAKEID2X];
            double zpj                          = cSim.pImageZ[SHAKEID2X];
            double xpk                          = cSim.pImageX[SHAKEID2Y];
            double ypk                          = cSim.pImageY[SHAKEID2Y];
            double zpk                          = cSim.pImageZ[SHAKEID2Y];
            psA->xpl                            = cSim.pImageX[SHAKEID2Z];
            psA->ypl                            = cSim.pImageY[SHAKEID2Z];
            psA->zpl                            = cSim.pImageZ[SHAKEID2Z];
            double xpm                          = cSim.pImageX[SHAKEID2W];
            double ypm                          = cSim.pImageY[SHAKEID2W];
            double zpm                          = cSim.pImageZ[SHAKEID2W];
#if !defined(use_SPSP) && (!defined(NODPTEXTURE) || (__CUDA_ARCH__ >= 200))
            double xi                           = __hiloint2double(ixi.y, ixi.x);
            double yi                           = __hiloint2double(iyi.y, iyi.x);
            double zi                           = __hiloint2double(izi.y, izi.x);
            double xij                          = __hiloint2double(ixij.y, ixij.x);
            double yij                          = __hiloint2double(iyij.y, iyij.x);
            double zij                          = __hiloint2double(izij.y, izij.x);
            double xik                          = __hiloint2double(ixik.y, ixik.x);
            double yik                          = __hiloint2double(iyik.y, iyik.x);
            double zik                          = __hiloint2double(izik.y, izik.x);   
            psA->xil                            = __hiloint2double(ixil.y, ixil.x);
            psA->yil                            = __hiloint2double(iyil.y, iyil.x);
            psA->zil                            = __hiloint2double(izil.y, izil.x);   
            double xim                          = __hiloint2double(ixim.y, ixim.x);
            double yim                          = __hiloint2double(iyim.y, iyim.x);
            double zim                          = __hiloint2double(izim.y, izim.x);   
#endif                       
            psA->invMassI                       = shakeParm.x;
            TOLER                               = shakeParm.y;
            
            // Calculate unchanging quantities
            xij                                 = xi - xij;
            yij                                 = yi - yij;
            zij                                 = zi - zij;
            xik                                 = xi - xik;
            yik                                 = yi - yik;
            zik                                 = zi - zik; 
            psA->xil                            = xi - psA->xil;
            psA->yil                            = yi - psA->yil;
            psA->zil                            = zi - psA->zil;
            xim                                 = xi - xim;
            yim                                 = yi - yim;
            zim                                 = zi - zim;               
           
            bool done                           = false;
            for (int i = 0; i < 3000; i++)
            {
                done = true;
                
                // Calculate nominal distance squared
                double xpxx                     = xpi - xpj;
                double ypxx                     = ypi - ypj;
                double zpxx                     = zpi - zpj;
                double rpxx2                    = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
          
                // Apply correction to first hydrogen
                double diff                     = TOLER - rpxx2;
                if (abs(diff) >= TOLER * cSim.tol)
                {
                    done                        = false;
                   
                    // Shake resetting of coordinate is done here
                    double rrpr                 = xij * xpxx + yij * ypxx + zij * zpxx;     
                    if (rrpr >= TOLER * 1.0e-06)
                    {
                    
                        double acor             = diff / (rrpr * (double)2.0 * (psA->invMassI + cSim.invMassH));
                        double h                = xij * acor;
                        xpi                    += h * psA->invMassI;
                        xpj                    -= h * cSim.invMassH;
                        h                       = yij * acor;
                        ypi                    += h * psA->invMassI;
                        ypj                    -= h * cSim.invMassH;
                        h                       = zij * acor;
                        zpi                    += h * psA->invMassI;
                        zpj                    -= h * cSim.invMassH;             
                    }
                }
          
     
                xpxx                            = xpi - xpk;
                ypxx                            = ypi - ypk;
                zpxx                            = zpi - zpk;
                rpxx2                           = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
          
                // Apply correction to second hydrogen
                diff                            = TOLER - rpxx2;
                if (abs(diff) >= TOLER * cSim.tol)
                {
                    done                        = false;
               
                    // Shake resetting of coordinate is done here
                    double rrpr                 = xik * xpxx + yik * ypxx + zik * zpxx;     
                    if (rrpr >= TOLER * 1.0e-06)
                    {
                    
                        double acor             = diff / (rrpr * 2.0 * (psA->invMassI + cSim.invMassH));
                        double h                = xik * acor;
                        xpi                    += h * psA->invMassI;
                        xpk                    -= h * cSim.invMassH;
                        h                       = yik * acor;
                        ypi                    += h * psA->invMassI;
                        ypk                    -= h * cSim.invMassH;
                        h                       = zik * acor;
                        zpi                    += h * psA->invMassI;
                        zpk                    -= h * cSim.invMassH;             
                    }
                }
                
 
                xpxx                            = xpi - psA->xpl;
                ypxx                            = ypi - psA->ypl;
                zpxx                            = zpi - psA->zpl;        
                rpxx2                           = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
          
                // Apply correction to third hydrogen
                diff                            = TOLER - rpxx2;
                if (abs(diff) >= TOLER * cSim.tol)
                {
                    done                        = false;
                  
                    // Shake resetting of coordinate is done here
                    double rrpr                 = psA->xil * xpxx + psA->yil * ypxx + psA->zil * zpxx;     
                    if (rrpr >= TOLER * 1.0e-06)
                    {             
                        double acor             = diff / (rrpr * (double)2.0 * (psA->invMassI + cSim.invMassH));
                        double h                = psA->xil * acor;
                        xpi                    += h * psA->invMassI;
                        psA->xpl               -= h * cSim.invMassH;
                        h                       = psA->yil * acor;
                        ypi                    += h * psA->invMassI;
                        psA->ypl               -= h * cSim.invMassH;
                        h                       = psA->zil * acor;
                        zpi                    += h * psA->invMassI;
                        psA->zpl               -= h * cSim.invMassH;             
                    }
                }

                xpxx                            = xpi - xpm;
                ypxx                            = ypi - ypm;
                zpxx                            = zpi - zpm;        
                rpxx2                           = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;

                // Apply correction to third hydrogen
                diff                            = TOLER - rpxx2;
                if (abs(diff) >= TOLER * cSim.tol)
                {
                    done                        = false;
                  
                    // Shake resetting of coordinate is done here
                    double rrpr                 = xim * xpxx + yim * ypxx + zim * zpxx;     
                    if (rrpr >= TOLER * 1.0e-06)
                    {             
                        double acor             = diff / (rrpr * (double)2.0 * (psA->invMassI + cSim.invMassH));
                        double h                = xim * acor;
                        xpi                    += h * psA->invMassI;
                        xpm                    -= h * cSim.invMassH;
                        h                       = yim * acor;
                        ypi                    += h * psA->invMassI;
                        ypm                    -= h * cSim.invMassH;
                        h                       = zim * acor;
                        zpi                    += h * psA->invMassI;
                        zpm                    -= h * cSim.invMassH;             
                    }
                }

                
                
                // Check for convergence
                if (done)
                    break;
            }
          
            // Write out results if converged, but there's no really good
            // way to indicate failure so we'll let the simulation heading
            // off to Neptune do that for us.  Wish there were a better way,
            // but until the CPU needs something from the GPU, those are the
            // the breaks.  I guess, technically, we could just set a flag to NOP
            // the simulation from here and then carry that result through upon
            // the next ntpr, ntwc, or ntwx update, but I leave that up to you 
            // guys to implement that (or not). 
            if (done)
            {
                cSim.pImageX[shakeID1]          = xpi;
                cSim.pImageY[shakeID1]          = ypi;
                cSim.pImageZ[shakeID1]          = zpi;
                cSim.pImageX[SHAKEID2X]         = xpj;
                cSim.pImageY[SHAKEID2X]         = ypj;
                cSim.pImageZ[SHAKEID2X]         = zpj;
                cSim.pImageX[SHAKEID2Y]         = xpk;
                cSim.pImageY[SHAKEID2Y]         = ypk;
                cSim.pImageZ[SHAKEID2Y]         = zpk;
                cSim.pImageX[SHAKEID2Z]         = psA->xpl;
                cSim.pImageY[SHAKEID2Z]         = psA->ypl;
                cSim.pImageZ[SHAKEID2Z]         = psA->zpl;
                cSim.pImageX[SHAKEID2W]         = xpm;
                cSim.pImageY[SHAKEID2W]         = ypm;
                cSim.pImageZ[SHAKEID2W]         = zpm;
            }        

    
            pos                                += gridDim.x * blockDim.x;  
#undef TOLER
#undef SHAKEID2
#undef SHAKEID2X
#undef SHAKEID2Y
#undef SHAKEID2Z
#undef SHAKEID2W                                  
        }
    }

}

void kShake(gpuContext gpu)
{
    texref.normalized = 0;
    texref.filterMode = cudaFilterModePoint;
    texref.addressMode[0] = cudaAddressModeClamp;
    texref.channelDesc.x = 32;
#ifndef use_SPSP    
    texref.channelDesc.y = 32;
#else    
    texref.channelDesc.y = 0;
#endif
    texref.channelDesc.z = 0;
    texref.channelDesc.w = 0;
#ifndef use_SPSP
    cudaBindTexture(NULL, texref, (int2*)(gpu->sim.pForce), gpu->sim.stride3 * sizeof(int2));
#else
    cudaBindTexture(NULL, texref, gpu->sim.pForce, gpu->sim.stride3 * sizeof(float));
#endif
    if (gpu->bNeighborList)
    {
        kPMEShake_kernel<<<gpu->blocks, gpu->shakeThreadsPerBlock>>>();  
    }
    else
    {    
        kShake_kernel<<<gpu->blocks, gpu->shakeThreadsPerBlock>>>();  
    }
    LAUNCHERROR("kShake");
    cudaUnbindTexture(texref);
}
