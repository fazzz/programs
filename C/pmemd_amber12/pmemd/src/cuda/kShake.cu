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

#if (__CUDA_ARCH__ >= 200)
#define INVMASSI invMassI
#define XPL xpl
#define YPL ypl
#define ZPL zpl
#define XIL xil
#define YIL yil
#define ZIL zil
#define XCOM xcom
#define YCOM ycom
#define ZCOM zcom
#define TRNS11 trns11
#define TRNS12 trns12
#define TRNS13 trns13
#define TRNS21 trns21
#define TRNS22 trns22
#define TRNS23 trns23
#define TRNS31 trns31
#else
#define INVMASSI psA->invMassI
#define XPL psA->xpl
#define YPL psA->ypl
#define ZPL psA->zpl
#define XIL psA->xil
#define YIL psA->yil
#define ZIL psA->zil
#define XCOM psA->xcom
#define YCOM psA->ycom
#define ZCOM psA->zcom
#define TRNS11 psA->trns11
#define TRNS12 psA->trns12
#define TRNS13 psA->trns13
#define TRNS21 psA->trns21
#define TRNS22 psA->trns22
#define TRNS23 psA->trns23
#define TRNS31 psA->trns31
#endif

// Texture reference for double-precision coordinates (disguised as int2 to work around HW limitations)
texture<int2, 1, cudaReadModeElementType> texref;

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
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_SHAKE_THREADS_PER_BLOCK, SM_3X_SHAKE_BLOCKS)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_SHAKE_THREADS_PER_BLOCK, SM_2X_SHAKE_BLOCKS)
#else
__launch_bounds__(SM_13_SHAKE_THREADS_PER_BLOCK, SM_13_SHAKE_BLOCKS)
#endif
kShake_kernel()
{
#if (__CUDA_ARCH__ >= 200)
    double invMassI;
    double xpl;
    double ypl;
    double zpl;
    double xil;
    double yil;
    double zil;
#else
    __shared__ Atom sA[SM_13_SHAKE_THREADS_PER_BLOCK];
    Atom* psA                                   = &sA[threadIdx.x];
#endif
    unsigned int pos                            = blockIdx.x * blockDim.x + threadIdx.x;   


    if (pos < cSim.shakeOffset)
    {
        if (pos < cSim.shakeConstraints)
        {   
            // Read SHAKE network data
            int4 shakeID                        = cSim.pShakeID[pos];
            double2 shakeParm                   = cSim.pShakeParm[pos];
            
            // Read SHAKE network components

#if defined(NODPTEXTURE)
            double xi                           = cSim.pForceX[shakeID.x];
            double yi                           = cSim.pForceY[shakeID.x];
            double zi                           = cSim.pForceZ[shakeID.x];
            double xij                          = cSim.pForceX[shakeID.y];
            double yij                          = cSim.pForceY[shakeID.y];
            double zij                          = cSim.pForceZ[shakeID.y];
#else        
            int2 ixi                            = tex1Dfetch(texref, shakeID.x);
            int2 iyi                            = tex1Dfetch(texref, shakeID.x + cSim.stride);
            int2 izi                            = tex1Dfetch(texref, shakeID.x + cSim.stride2);
            int2 ixij                           = tex1Dfetch(texref, shakeID.y);
            int2 iyij                           = tex1Dfetch(texref, shakeID.y + cSim.stride);
            int2 izij                           = tex1Dfetch(texref, shakeID.y + cSim.stride2);
            double xi                           = __hiloint2double(ixi.y, ixi.x);
            double yi                           = __hiloint2double(iyi.y, iyi.x);
            double zi                           = __hiloint2double(izi.y, izi.x);
            double xij                          = __hiloint2double(ixij.y, ixij.x);
            double yij                          = __hiloint2double(iyij.y, iyij.x);
            double zij                          = __hiloint2double(izij.y, izij.x); 
#endif                  
            double xpi                          = cSim.pAtomX[shakeID.x];
            double ypi                          = cSim.pAtomY[shakeID.x];
            double zpi                          = cSim.pAtomZ[shakeID.x];
            double xpj                          = cSim.pAtomX[shakeID.y];
            double ypj                          = cSim.pAtomY[shakeID.y];
            double zpj                          = cSim.pAtomZ[shakeID.y];               
            INVMASSI                            = shakeParm.x;
            double toler                        = shakeParm.y;
        
        
            // Optionally read 2nd hydrogen
            double xpk, ypk, zpk, xik, yik, zik;
            if (shakeID.z != -1)
            {
#if defined(NODPTEXTURE)
                xik                             = cSim.pForceX[shakeID.z];
                yik                             = cSim.pForceY[shakeID.z];
                zik                             = cSim.pForceZ[shakeID.z];   
#else
                int2 ixik                       = tex1Dfetch(texref, shakeID.z);
                int2 iyik                       = tex1Dfetch(texref, shakeID.z + cSim.stride);
                int2 izik                       = tex1Dfetch(texref, shakeID.z + cSim.stride2);
                xik                             = __hiloint2double(ixik.y, ixik.x);
                yik                             = __hiloint2double(iyik.y, iyik.x);
                zik                             = __hiloint2double(izik.y, izik.x);  
#endif 
                xpk                             = cSim.pAtomX[shakeID.z];
                ypk                             = cSim.pAtomY[shakeID.z];
                zpk                             = cSim.pAtomZ[shakeID.z];               
            }
            
            // Optionally read 3rd hydrogen into shared memory
            if (shakeID.w != -1)
            {
#if defined(NODPTEXTURE) 
                XIL                             = cSim.pForceX[shakeID.w];
                YIL                             = cSim.pForceY[shakeID.w];
                ZIL                             = cSim.pForceZ[shakeID.w];       
#else            
                int2 ixil                       = tex1Dfetch(texref, shakeID.w);
                int2 iyil                       = tex1Dfetch(texref, shakeID.w + cSim.stride);
                int2 izil                       = tex1Dfetch(texref, shakeID.w + cSim.stride2); 
                XIL                             = __hiloint2double(ixil.y, ixil.x);
                YIL                             = __hiloint2double(iyil.y, iyil.x);
                ZIL                             = __hiloint2double(izil.y, izil.x);   
#endif               
                XPL                             = cSim.pAtomX[shakeID.w];
                YPL                             = cSim.pAtomY[shakeID.w];
                ZPL                             = cSim.pAtomZ[shakeID.w];           
            }
            
            // Calculate unchanging quantities
            xij                                 = xi - xij;
            yij                                 = yi - yij;
            zij                                 = zi - zij;
            
            if (shakeID.z != -1)
            {
                xik                             = xi - xik;
                yik                             = yi - yik;
                zik                             = zi - zik;
            }        
             
            if (shakeID.w != -1)
            {
                XIL                             = xi - XIL;
                YIL                             = yi - YIL;
                ZIL                             = zi - ZIL;
            }      
       
            bool done                           = false;
            for (int i = 0; i < 3000; i++)
            {
                done                            = true;
                
                // Calculate nominal distance squared
                double xpxx                     = xpi - xpj;
                double ypxx                     = ypi - ypj;
                double zpxx                     = zpi - zpj;
                double rpxx2                    = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
      
                // Apply correction
                double diff                     = toler - rpxx2;
                if (abs(diff) >= toler * cSim.tol)
                {
                    done                        = false;
                   
                    // Shake resetting of coordinate is done here
                    double rrpr                 = xij * xpxx + yij * ypxx + zij * zpxx;     
                    if (rrpr >= toler * 1.0e-06)
                    {
                    
                        double acor             = diff / (rrpr * 2.0 * (INVMASSI + cSim.invMassH));
                        double h                = xij * acor;
                        xpi                    += h * INVMASSI;
                        xpj                    -= h * cSim.invMassH;
                        h                       = yij * acor;
                        ypi                    += h * INVMASSI;
                        ypj                    -= h * cSim.invMassH;
                        h                       = zij * acor;
                        zpi                    += h * INVMASSI;
                        zpj                    -= h * cSim.invMassH;             
                    }
                }
      
                // Second bond if present
                if (shakeID.z != -1)
                {
                    xpxx                        = xpi - xpk;
                    ypxx                        = ypi - ypk;
                    zpxx                        = zpi - zpk;
                    rpxx2                       = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
          
                    // Apply correction
                    diff                        = toler - rpxx2;
                    if (abs(diff) >= toler * cSim.tol)
                    {
                        done                    = false;
                   
                        // Shake resetting of coordinate is done here
                        double rrpr             = xik * xpxx + yik * ypxx + zik * zpxx;     
                        if (rrpr >= toler * 1.0e-06)
                        {
                    
                            double acor         = diff / (rrpr * 2.0 * (INVMASSI + cSim.invMassH));
                            double h            = xik * acor;
                            xpi                += h * INVMASSI;
                            xpk                -= h * cSim.invMassH;
                            h                   = yik * acor;
                            ypi                += h * INVMASSI;
                            ypk                -= h * cSim.invMassH;
                            h                   = zik * acor;
                            zpi                += h * INVMASSI;
                            zpk                -= h * cSim.invMassH;             
                        }
                    }
                }
            
                // Third bond if present
                if (shakeID.w != -1)
                {
                    xpxx                        = xpi - XPL;
                    ypxx                        = ypi - YPL;
                    zpxx                        = zpi - ZPL;
                    rpxx2                       = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
          
                    // Apply correction
                    diff                        = toler - rpxx2;
                    if (abs(diff) >= toler * cSim.tol)
                    {
                        done                    = false;
                   
                        // Shake resetting of coordinate is done here
                        double rrpr             = XIL * xpxx + YIL * ypxx + ZIL * zpxx;     
                        if (rrpr >= toler * 1.0e-06)
                        {
                    
                            double acor         = diff / (rrpr * 2.0 * (INVMASSI + cSim.invMassH));
                            double h            = XIL * acor;
                            xpi                += h * INVMASSI;
                            XPL                -= h * cSim.invMassH;
                            h                   = YIL * acor;
                            ypi                += h * INVMASSI;
                            YPL                -= h * cSim.invMassH;
                            h                   = ZIL * acor;
                            zpi                += h * INVMASSI;
                            ZPL                -= h * cSim.invMassH;             
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
                cSim.pAtomX[shakeID.x]          = xpi;
                cSim.pAtomY[shakeID.x]          = ypi;
                cSim.pAtomZ[shakeID.x]          = zpi;
                PMEFloat2 xyi                   = {xpi, ypi};
                cSim.pAtomXYSP[shakeID.x]       = xyi;
                cSim.pAtomZSP[shakeID.x]        = zpi;
                
                cSim.pAtomX[shakeID.y]          = xpj;
                cSim.pAtomY[shakeID.y]          = ypj;
                cSim.pAtomZ[shakeID.y]          = zpj;
                PMEFloat2 xyj                   = {xpj, ypj};
                cSim.pAtomXYSP[shakeID.y]       = xyj;
                cSim.pAtomZSP[shakeID.y]        = zpj;

                if (shakeID.z != -1)
                {
                    cSim.pAtomX[shakeID.z]      = xpk;
                    cSim.pAtomY[shakeID.z]      = ypk;
                    cSim.pAtomZ[shakeID.z]      = zpk;
                    PMEFloat2 xyk               = {xpk, ypk};
                    cSim.pAtomXYSP[shakeID.z]   = xyk;
                    cSim.pAtomZSP[shakeID.z]    = zpk;
                }
            
                if (shakeID.w != -1)
                {
                    cSim.pAtomX[shakeID.w]      = XPL;
                    cSim.pAtomY[shakeID.w]      = YPL;
                    cSim.pAtomZ[shakeID.w]      = ZPL;
                    PMEFloat2 xyl               = {XPL, YPL};
                    cSim.pAtomXYSP[shakeID.w]   = xyl;
                    cSim.pAtomZSP[shakeID.w]    = ZPL;
                }
            
            }
        }
    } 
    else if (pos < cSim.fastShakeOffset)
    {    
        pos                                    -= cSim.shakeOffset;
        if (pos < cSim.fastShakeConstraints)
        {
            // Read atom data
            int4 shakeID                        = cSim.pFastShakeID[pos];
#if defined(NODPTEXTURE)
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
            double xp1                          = cSim.pAtomX[shakeID.x];
            double yp1                          = cSim.pAtomY[shakeID.x];
            double zp1                          = cSim.pAtomZ[shakeID.x];
            double xp2                          = cSim.pAtomX[shakeID.y];
            double yp2                          = cSim.pAtomY[shakeID.y];
            double zp2                          = cSim.pAtomZ[shakeID.y];
            double xp3                          = cSim.pAtomX[shakeID.z];
            double yp3                          = cSim.pAtomY[shakeID.z];
            double zp3                          = cSim.pAtomZ[shakeID.z];
    
            // Step1  A1_prime:
            double xb0                          = x2 - x1;
            double yb0                          = y2 - y1;
            double zb0                          = z2 - z1;
            double xc0                          = x3 - x1;
            double yc0                          = y3 - y1;
            double zc0                          = z3 - z1;

            XPL                                 = xp1 * cSim.wo_div_wohh + (xp2 + xp3) * cSim.wh_div_wohh;
            YPL                                 = yp1 * cSim.wo_div_wohh + (yp2 + yp3) * cSim.wh_div_wohh;
            ZPL                                 = zp1 * cSim.wo_div_wohh + (zp2 + zp3) * cSim.wh_div_wohh;

            double xa1                          = xp1 - XPL;
            double ya1                          = yp1 - YPL;
            double za1                          = zp1 - ZPL;
            double xb1                          = xp2 - XPL;
            double yb1                          = yp2 - YPL;
            double zb1                          = zp2 - ZPL;
            double xc1                          = xp3 - XPL;
            double yc1                          = yp3 - YPL;
            double zc1                          = zp3 - ZPL;

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
            cSim.pAtomX[shakeID.x]              = XPL + trns11 * xa3d + trns12 * ya3d + trns13 * za3d;
            cSim.pAtomY[shakeID.x]              = YPL + trns21 * xa3d + trns22 * ya3d + trns23 * za3d;
            cSim.pAtomZ[shakeID.x]              = ZPL + trns31 * xa3d + trns32 * ya3d + trns33 * za3d;
            cSim.pAtomX[shakeID.y]              = XPL + trns11 * xb3d + trns12 * yb3d + trns13 * zb3d;
            cSim.pAtomY[shakeID.y]              = YPL + trns21 * xb3d + trns22 * yb3d + trns23 * zb3d;
            cSim.pAtomZ[shakeID.y]              = ZPL + trns31 * xb3d + trns32 * yb3d + trns33 * zb3d;
            cSim.pAtomX[shakeID.z]              = XPL + trns11 * xc3d + trns12 * yc3d + trns13 * zc3d;
            cSim.pAtomY[shakeID.z]              = YPL + trns21 * xc3d + trns22 * yc3d + trns23 * zc3d;
            cSim.pAtomZ[shakeID.z]              = ZPL + trns31 * xc3d + trns32 * yc3d + trns33 * zc3d;                        
        }
    }
    else if ( pos < cSim.slowShakeOffset)
    {    
        pos                                    -= cSim.fastShakeOffset;

        if (pos < cSim.slowShakeConstraints)
        {
            int shakeID1;
		    int4 shakeID2;
            double toler;

            // Read SHAKE network data
            shakeID1                            = cSim.pSlowShakeID1[pos];
            shakeID2                            = cSim.pSlowShakeID2[pos];
            double2 shakeParm                   = cSim.pSlowShakeParm[pos];
        
            // Read SHAKE network components
#if defined(NODPTEXTURE)
            double xi                           = cSim.pForceX[shakeID1];
            double yi                           = cSim.pForceY[shakeID1];
            double zi                           = cSim.pForceZ[shakeID1];
            double xij                          = cSim.pForceX[shakeID2.x];
            double yij                          = cSim.pForceY[shakeID2.x];
            double zij                          = cSim.pForceZ[shakeID2.x];
            double xik                          = cSim.pForceX[shakeID2.y];
            double yik                          = cSim.pForceY[shakeID2.y];
            double zik                          = cSim.pForceZ[shakeID2.y];  
            XIL                                 = cSim.pForceX[shakeID2.z];
            YIL                                 = cSim.pForceY[shakeID2.z];
            ZIL                                 = cSim.pForceZ[shakeID2.z]; 
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
            double xi                           = __hiloint2double(ixi.y, ixi.x);
            double yi                           = __hiloint2double(iyi.y, iyi.x);
            double zi                           = __hiloint2double(izi.y, izi.x);
            double xij                          = __hiloint2double(ixij.y, ixij.x);
            double yij                          = __hiloint2double(iyij.y, iyij.x);
            double zij                          = __hiloint2double(izij.y, izij.x);
            double xik                          = __hiloint2double(ixik.y, ixik.x);
            double yik                          = __hiloint2double(iyik.y, iyik.x);
            double zik                          = __hiloint2double(izik.y, izik.x);   
            XIL                                 = __hiloint2double(ixil.y, ixil.x);
            YIL                                 = __hiloint2double(iyil.y, iyil.x);
            ZIL                                 = __hiloint2double(izil.y, izil.x);   
            double xim                          = __hiloint2double(ixim.y, ixim.x);
            double yim                          = __hiloint2double(iyim.y, iyim.x);
            double zim                          = __hiloint2double(izim.y, izim.x); 
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
            XPL                                 = cSim.pAtomX[shakeID2.z];
            YPL                                 = cSim.pAtomY[shakeID2.z];
            ZPL                                 = cSim.pAtomZ[shakeID2.z];
            double xpm                          = cSim.pAtomX[shakeID2.w];
            double ypm                          = cSim.pAtomY[shakeID2.w];
            double zpm                          = cSim.pAtomZ[shakeID2.w];                    
            INVMASSI                            = shakeParm.x;
            toler                               = shakeParm.y;
            
            // Calculate unchanging quantities
            xij                                 = xi - xij;
            yij                                 = yi - yij;
            zij                                 = zi - zij;
            xik                                 = xi - xik;
            yik                                 = yi - yik;
            zik                                 = zi - zik; 
            XIL                                 = xi - XIL;
            YIL                                 = yi - YIL;
            ZIL                                 = zi - ZIL;
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
                    
                        double acor             = diff / (rrpr * (double)2.0 * (INVMASSI + cSim.invMassH));
                        double h                = xij * acor;
                        xpi                    += h * INVMASSI;
                        xpj                    -= h * cSim.invMassH;
                        h                       = yij * acor;
                        ypi                    += h * INVMASSI;
                        ypj                    -= h * cSim.invMassH;
                        h                       = zij * acor;
                        zpi                    += h * INVMASSI;
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
                    
                        double acor             = diff / (rrpr * 2.0 * (INVMASSI + cSim.invMassH));
                        double h                = xik * acor;
                        xpi                    += h * INVMASSI;
                        xpk                    -= h * cSim.invMassH;
                        h                       = yik * acor;
                        ypi                    += h * INVMASSI;
                        ypk                    -= h * cSim.invMassH;
                        h                       = zik * acor;
                        zpi                    += h * INVMASSI;
                        zpk                    -= h * cSim.invMassH;             
                    }
                }
                
 
                xpxx                            = xpi - XPL;
                ypxx                            = ypi - YPL;
                zpxx                            = zpi - ZPL;        
                rpxx2                           = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
          
                // Apply correction to third hydrogen
                diff                            = toler - rpxx2;
                if (abs(diff) >= toler * cSim.tol)
                {
                    done                        = false;
                  
                    // Shake resetting of coordinate is done here
                    double rrpr                 = XIL * xpxx + YIL * ypxx + ZIL * zpxx;     
                    if (rrpr >= toler * 1.0e-06)
                    {             
                        double acor             = diff / (rrpr * (double)2.0 * (INVMASSI + cSim.invMassH));
                        double h                = XIL * acor;
                        xpi                    += h * INVMASSI;
                        XPL                    -= h * cSim.invMassH;
                        h                       = YIL * acor;
                        ypi                    += h * INVMASSI;
                        YPL                    -= h * cSim.invMassH;
                        h                       = ZIL * acor;
                        zpi                    += h * INVMASSI;
                        ZPL                    -= h * cSim.invMassH;             
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
                        double acor             = diff / (rrpr * (double)2.0 * (INVMASSI + cSim.invMassH));
                        double h                = xim * acor;
                        xpi                    += h * INVMASSI;
                        xpm                    -= h * cSim.invMassH;
                        h                       = yim * acor;
                        ypi                    += h * INVMASSI;
                        ypm                    -= h * cSim.invMassH;
                        h                       = zim * acor;
                        zpi                    += h * INVMASSI;
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
                cSim.pAtomX[shakeID2.z]         = XPL;
                cSim.pAtomY[shakeID2.z]         = YPL;
                cSim.pAtomZ[shakeID2.z]         = ZPL;
                PMEFloat2 xyl                   = {XPL, YPL};
                cSim.pAtomXYSP[shakeID2.z]      = xyl;
                cSim.pAtomZSP[shakeID2.z]       = ZPL;
                cSim.pAtomX[shakeID2.w]         = xpm;
                cSim.pAtomY[shakeID2.w]         = ypm;
                cSim.pAtomZ[shakeID2.w]         = zpm;
                PMEFloat2 xym                   = {xpm, ypm};
                cSim.pAtomXYSP[shakeID2.w]      = xym;
                cSim.pAtomZSP[shakeID2.w]       = zpm;
            }                                      
        }
    }
}


#if (__CUDA_ARCH__ < 200)
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
#if (__CUDA_ARCH__ >= 300)
__launch_bounds__(SM_3X_SHAKE_THREADS_PER_BLOCK, SM_3X_SHAKE_BLOCKS)
#elif (__CUDA_ARCH__ >= 200)
__launch_bounds__(SM_2X_SHAKE_THREADS_PER_BLOCK, SM_2X_SHAKE_BLOCKS)
#else
__launch_bounds__(SM_13_SHAKE_THREADS_PER_BLOCK, SM_13_SHAKE_BLOCKS)
#endif
kPMEShake_kernel()
{
#if (__CUDA_ARCH__ < 200)
__shared__ PMEAtom sA[SM_13_SHAKE_THREADS_PER_BLOCK];
#endif

    unsigned int pos                            = blockIdx.x * blockDim.x + threadIdx.x;  
    if (pos < cSim.shakeOffset)
    { 
        if (pos < cSim.shakeConstraints)
        {
#if (__CUDA_ARCH__ >= 200)
            double invMassI;
            double xpl;
            double ypl;
            double zpl;
            double xil;
            double yil;
            double zil;
#else
            PMEAtom* psA                        = &sA[threadIdx.x];
#endif
            // Read SHAKE network data
            int4 shakeID                    	= cSim.pImageShakeID[pos];
            double2 shakeParm                   = cSim.pShakeParm[pos];
            
            // Read SHAKE network components
#if defined(NODPTEXTURE)
            double xi                           = cSim.pForceX[shakeID.x];
            double yi                           = cSim.pForceY[shakeID.x];
            double zi                           = cSim.pForceZ[shakeID.x];
            double xij                          = cSim.pForceX[shakeID.y];
            double yij                          = cSim.pForceY[shakeID.y];
            double zij                          = cSim.pForceZ[shakeID.y];
#else        
            int2 ixi                            = tex1Dfetch(texref, shakeID.x);
            int2 iyi                            = tex1Dfetch(texref, shakeID.x + cSim.stride);
            int2 izi                            = tex1Dfetch(texref, shakeID.x + cSim.stride2);
            int2 ixij                           = tex1Dfetch(texref, shakeID.y);
            int2 iyij                           = tex1Dfetch(texref, shakeID.y + cSim.stride);
            int2 izij                           = tex1Dfetch(texref, shakeID.y + cSim.stride2);
            double xi                           = __hiloint2double(ixi.y, ixi.x);
            double yi                           = __hiloint2double(iyi.y, iyi.x);
            double zi                           = __hiloint2double(izi.y, izi.x);
            double xij                          = __hiloint2double(ixij.y, ixij.x);
            double yij                          = __hiloint2double(iyij.y, iyij.x);
            double zij                          = __hiloint2double(izij.y, izij.x); 
#endif 
            double xpi                          = cSim.pImageX[shakeID.x];
            double ypi                          = cSim.pImageY[shakeID.x];
            double zpi                          = cSim.pImageZ[shakeID.x];
            double xpj                          = cSim.pImageX[shakeID.y];
            double ypj                          = cSim.pImageY[shakeID.y];
            double zpj                          = cSim.pImageZ[shakeID.y];                      
            INVMASSI                            = shakeParm.x;
            double toler                        = shakeParm.y;
            
            
            // Optionally read 2nd hydrogen
            double xpk, ypk, zpk, xik, yik, zik;
            if (shakeID.z != -1)
            {
#if defined(NODPTEXTURE)
                xik                             = cSim.pForceX[shakeID.z];
                yik                             = cSim.pForceY[shakeID.z];
                zik                             = cSim.pForceZ[shakeID.z];    
#else
                int2 ixik                       = tex1Dfetch(texref, shakeID.z);
                int2 iyik                       = tex1Dfetch(texref, shakeID.z + cSim.stride);
                int2 izik                       = tex1Dfetch(texref, shakeID.z + cSim.stride2);
                xik                             = __hiloint2double(ixik.y, ixik.x);
                yik                             = __hiloint2double(iyik.y, iyik.x);
                zik                             = __hiloint2double(izik.y, izik.x);  
#endif 
                xpk                             = cSim.pImageX[shakeID.z];
                ypk                             = cSim.pImageY[shakeID.z];
                zpk                             = cSim.pImageZ[shakeID.z];                
            }
        
            // Optionally read 3rd hydrogen into shared memory
            if (shakeID.w != -1)
            {
#if defined(NODPTEXTURE)
                XIL                             = cSim.pForceX[shakeID.w];
                YIL                             = cSim.pForceY[shakeID.w];
                ZIL                             = cSim.pForceZ[shakeID.w];      
#else            
                int2 ixil                       = tex1Dfetch(texref, shakeID.w);
                int2 iyil                       = tex1Dfetch(texref, shakeID.w + cSim.stride);
                int2 izil                       = tex1Dfetch(texref, shakeID.w + cSim.stride2); 
                XIL                             = __hiloint2double(ixil.y, ixil.x);
                YIL                             = __hiloint2double(iyil.y, iyil.x);
                ZIL                             = __hiloint2double(izil.y, izil.x);
#endif          
                XPL                             = cSim.pImageX[shakeID.w];
                YPL                             = cSim.pImageY[shakeID.w];
                ZPL                             = cSim.pImageZ[shakeID.w];          
            }
        
            // Calculate unchanging quantities
            xij                                 = xi - xij;
            yij                                 = yi - yij;
            zij                                 = zi - zij;
            
            if (shakeID.z != -1)
            {
                xik                             = xi - xik;
                yik                             = yi - yik;
                zik                             = zi - zik; 
            }        
             
            if (shakeID.w != -1)
            {
                XIL                             = xi - XIL;
                YIL                             = yi - YIL;
                ZIL                             = zi - ZIL;
            }      
       
            bool done                           = false;
            for (int i = 0; i < 3000; i++)
            {
                done                            = true;
                
                // Calculate nominal distance squared
                double xpxx                     = xpi - xpj;
                double ypxx                     = ypi - ypj;
                double zpxx                     = zpi - zpj;
                double rpxx2                    = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
          
                // Apply correction
                double diff                     = toler - rpxx2;
                if (abs(diff) >= toler * cSim.tol)
                {
                    done                        = false;
                   
                    // Shake resetting of coordinate is done here
                    double rrpr                 = xij * xpxx + yij * ypxx + zij * zpxx;     
                    if (rrpr >= toler * 1.0e-06)
                    {
                    
                        double acor             = diff / (rrpr * (double)2.0 * (INVMASSI + cSim.invMassH));
                        double h                = xij * acor;
                        xpi                    += h * INVMASSI;
                        xpj                    -= h * cSim.invMassH;
                        h                       = yij * acor;
                        ypi                    += h * INVMASSI;
                        ypj                    -= h * cSim.invMassH;
                        h                       = zij * acor;
                        zpi                    += h * INVMASSI;
                        zpj                    -= h * cSim.invMassH;             
                    }
                }
      
                // Second bond if present
                if (shakeID.z != -1)
                {
                    xpxx                        = xpi - xpk;
                    ypxx                        = ypi - ypk;
                    zpxx                        = zpi - zpk;
                    rpxx2                       = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
          
                    // Apply correction
                    diff                        = toler - rpxx2;
                    if (abs(diff) >= toler * cSim.tol)
                    {
                        done                    = false;
                   
                        // Shake resetting of coordinate is done here
                        double rrpr             = xik * xpxx + yik * ypxx + zik * zpxx;     
                        if (rrpr >= toler * 1.0e-06)
                        {
                    
                            double acor         = diff / (rrpr * 2.0 * (INVMASSI + cSim.invMassH));
                            double h            = xik * acor;
                            xpi                += h * INVMASSI;
                            xpk                -= h * cSim.invMassH;
                            h                   = yik * acor;
                            ypi                += h * INVMASSI;
                            ypk                -= h * cSim.invMassH;
                            h                   = zik * acor;
                            zpi                += h * INVMASSI;
                            zpk                -= h * cSim.invMassH;             
                        }
                    }
                }
            
                // Third bond if present
                if (shakeID.w != -1)
                {
                    xpxx                        = xpi - XPL;
                    ypxx                        = ypi - YPL;
                    zpxx                        = zpi - ZPL;        
                    rpxx2                       = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
          
                    // Apply correction
                    diff                        = toler - rpxx2;
                    if (abs(diff) >= toler * cSim.tol)
                    {
                        done                    = false;
                   
                        // Shake resetting of coordinate is done here
                        double rrpr             = XIL * xpxx + YIL * ypxx + ZIL * zpxx;     
                        if (rrpr >= toler * 1.0e-06)
                        {
                    
                            double acor         = diff / (rrpr * (double)2.0 * (INVMASSI + cSim.invMassH));
                            double h            = XIL * acor;
                            xpi                += h * INVMASSI;
                            XPL                -= h * cSim.invMassH;
                            h                   = YIL * acor;
                            ypi                += h * INVMASSI;
                            YPL                -= h * cSim.invMassH;
                            h                   = ZIL * acor;
                            zpi                += h * INVMASSI;
                            ZPL                -= h * cSim.invMassH;             
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
                cSim.pImageX[shakeID.x]          = xpi;
                cSim.pImageY[shakeID.x]          = ypi;
                cSim.pImageZ[shakeID.x]          = zpi;
                cSim.pImageX[shakeID.y]          = xpj;
                cSim.pImageY[shakeID.y]          = ypj;
                cSim.pImageZ[shakeID.y]          = zpj;

                if (shakeID.z != -1)
                {
                    cSim.pImageX[shakeID.z]      = xpk;
                    cSim.pImageY[shakeID.z]      = ypk;
                    cSim.pImageZ[shakeID.z]      = zpk;
                }

                if (shakeID.w != -1)
                {
                    cSim.pImageX[shakeID.w]      = XPL;
                    cSim.pImageY[shakeID.w]      = YPL;
                    cSim.pImageZ[shakeID.w]      = ZPL;
                }
            }        
        }
    }
    else if (pos < cSim.fastShakeOffset)
    {    
        pos                                     -= cSim.shakeOffset;

        if (pos < cSim.fastShakeConstraints)
        {
#if (__CUDA_ARCH__ >= 200)
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
#else
            PMEFastAtom* psA                    = (PMEFastAtom*)&sA[threadIdx.x];
#endif            

            // Read atom data
            int4 shakeID                        = cSim.pImageFastShakeID[pos];
#if defined(NODPTEXTURE)
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
            double xp1                          = cSim.pImageX[shakeID.x];
            double yp1                          = cSim.pImageY[shakeID.x];
            double zp1                          = cSim.pImageZ[shakeID.x];
            double xp2                          = cSim.pImageX[shakeID.y];
            double yp2                          = cSim.pImageY[shakeID.y];
            double zp2                          = cSim.pImageZ[shakeID.y];
            double xp3                          = cSim.pImageX[shakeID.z];
            double yp3                          = cSim.pImageY[shakeID.z];
            double zp3                          = cSim.pImageZ[shakeID.z];

            // Step1  A1_prime:
            double xb0                          = x2 - x1;
            double yb0                          = y2 - y1;
            double zb0                          = z2 - z1;
            double xc0                          = x3 - x1;
            double yc0                          = y3 - y1;
            double zc0                          = z3 - z1;
            XCOM                                = xp1 * cSim.wo_div_wohh + (xp2 + xp3) * cSim.wh_div_wohh;
            YCOM                                = yp1 * cSim.wo_div_wohh + (yp2 + yp3) * cSim.wh_div_wohh;
            ZCOM                                = zp1 * cSim.wo_div_wohh + (zp2 + zp3) * cSim.wh_div_wohh;

            double xa1                          = xp1 - XCOM;
            double ya1                          = yp1 - YCOM;
            double za1                          = zp1 - ZCOM;
            double xb1                          = xp2 - XCOM;
            double yb1                          = yp2 - YCOM;
            double zb1                          = zp2 - ZCOM;
            double xc1                          = xp3 - XCOM;
            double yc1                          = yp3 - YCOM;
            double zc1                          = zp3 - ZCOM;
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

            TRNS11                              = xaksxd * axlng_inv;
            TRNS21                              = yaksxd * axlng_inv;
            TRNS31                              = zaksxd * axlng_inv;
            TRNS12                              = xaksyd * aylng_inv;
            TRNS22                              = yaksyd * aylng_inv;
			double trns32                       = zaksyd * aylng_inv;
            TRNS13                              = xakszd * azlng_inv;
            TRNS23                              = yakszd * azlng_inv;
			double trns33                       = zakszd * azlng_inv;

            double xb0d                         = TRNS11 * xb0 + TRNS21 * yb0 + TRNS31 * zb0;
            double yb0d                         = TRNS12 * xb0 + TRNS22 * yb0 + trns32 * zb0;
            double xc0d                         = TRNS11 * xc0 + TRNS21 * yc0 + TRNS31 * zc0;
            double yc0d                         = TRNS12 * xc0 + TRNS22 * yc0 + trns32 * zc0;
            double za1d                         = TRNS13 * xa1 + TRNS23 * ya1 + trns33 * za1;
            double xb1d                         = TRNS11 * xb1 + TRNS21 * yb1 + TRNS31 * zb1;
            double yb1d                         = TRNS12 * xb1 + TRNS22 * yb1 + trns32 * zb1;
            double zb1d                         = TRNS13 * xb1 + TRNS23 * yb1 + trns33 * zb1;
            double xc1d                         = TRNS11 * xc1 + TRNS21 * yc1 + TRNS31 * zc1;
            double yc1d                         = TRNS12 * xc1 + TRNS22 * yc1 + trns32 * zc1;
            double zc1d                         = TRNS13 * xc1 + TRNS23 * yc1 + trns33 * zc1;

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
            cSim.pImageX[shakeID.x]             = XCOM + TRNS11 * xa3d + TRNS12 * ya3d + TRNS13 * za3d;
            cSim.pImageY[shakeID.x]             = YCOM + TRNS21 * xa3d + TRNS22 * ya3d + TRNS23 * za3d;
            cSim.pImageZ[shakeID.x]             = ZCOM + TRNS31 * xa3d + trns32 * ya3d + trns33 * za3d;
            cSim.pImageX[shakeID.y]             = XCOM + TRNS11 * xb3d + TRNS12 * yb3d + TRNS13 * zb3d;
            cSim.pImageY[shakeID.y]             = YCOM + TRNS21 * xb3d + TRNS22 * yb3d + TRNS23 * zb3d;
            cSim.pImageZ[shakeID.y]             = ZCOM + TRNS31 * xb3d + trns32 * yb3d + trns33 * zb3d;
            cSim.pImageX[shakeID.z]             = XCOM + TRNS11 * xc3d + TRNS12 * yc3d + TRNS13 * zc3d;
            cSim.pImageY[shakeID.z]             = YCOM + TRNS21 * xc3d + TRNS22 * yc3d + TRNS23 * zc3d;
            cSim.pImageZ[shakeID.z]             = ZCOM + TRNS31 * xc3d + trns32 * yc3d + trns33 * zc3d;                                 
        }
    }
    else if (pos < cSim.slowShakeOffset)
    {    
        pos                                    -= cSim.fastShakeOffset;
        if (pos < cSim.slowShakeConstraints)
        {
            // Read SHAKE network data
#if (__CUDA_ARCH__ >= 200)
            double invMassI;
            double xpl;
            double ypl;
            double zpl;
            double xil;
            double yil;
            double zil;
#else
            PMEAtom* psA                        = &sA[threadIdx.x];
#endif
            int shakeID1                        = cSim.pImageSlowShakeID1[pos];
            int4 shakeID2                       = cSim.pImageSlowShakeID2[pos];
            double2 shakeParm                   = cSim.pSlowShakeParm[pos];
        
            // Read SHAKE network components
#if defined(NODPTEXTURE)
            double xi                           = cSim.pForceX[shakeID1];
            double yi                           = cSim.pForceY[shakeID1];
            double zi                           = cSim.pForceZ[shakeID1];
            double xij                          = cSim.pForceX[shakeID2.x];
            double yij                          = cSim.pForceY[shakeID2.x];
            double zij                          = cSim.pForceZ[shakeID2.x];
            double xik                          = cSim.pForceX[shakeID2.y];
            double yik                          = cSim.pForceY[shakeID2.y];
            double zik                          = cSim.pForceZ[shakeID2.y];  
            XIL                                 = cSim.pForceX[shakeID2.z];
            YIL                                 = cSim.pForceY[shakeID2.z];
            ZIL                                 = cSim.pForceZ[shakeID2.z]; 
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
            double xi                           = __hiloint2double(ixi.y, ixi.x);
            double yi                           = __hiloint2double(iyi.y, iyi.x);
            double zi                           = __hiloint2double(izi.y, izi.x);
            double xij                          = __hiloint2double(ixij.y, ixij.x);
            double yij                          = __hiloint2double(iyij.y, iyij.x);
            double zij                          = __hiloint2double(izij.y, izij.x);
            double xik                          = __hiloint2double(ixik.y, ixik.x);
            double yik                          = __hiloint2double(iyik.y, iyik.x);
            double zik                          = __hiloint2double(izik.y, izik.x);   
            XIL                                 = __hiloint2double(ixil.y, ixil.x);
            YIL                                 = __hiloint2double(iyil.y, iyil.x);
            ZIL                                 = __hiloint2double(izil.y, izil.x);   
            double xim                          = __hiloint2double(ixim.y, ixim.x);
            double yim                          = __hiloint2double(iyim.y, iyim.x);
            double zim                          = __hiloint2double(izim.y, izim.x);   
#endif 
            double xpi                          = cSim.pImageX[shakeID1];
            double ypi                          = cSim.pImageY[shakeID1];
            double zpi                          = cSim.pImageZ[shakeID1];
            double xpj                          = cSim.pImageX[shakeID2.x];
            double ypj                          = cSim.pImageY[shakeID2.x];
            double zpj                          = cSim.pImageZ[shakeID2.x];
            double xpk                          = cSim.pImageX[shakeID2.y];
            double ypk                          = cSim.pImageY[shakeID2.y];
            double zpk                          = cSim.pImageZ[shakeID2.y];
            XPL                                 = cSim.pImageX[shakeID2.z];
            YPL                                 = cSim.pImageY[shakeID2.z];
            ZPL                                 = cSim.pImageZ[shakeID2.z];
            double xpm                          = cSim.pImageX[shakeID2.w];
            double ypm                          = cSim.pImageY[shakeID2.w];
            double zpm                          = cSim.pImageZ[shakeID2.w];                      
            INVMASSI                            = shakeParm.x;
            double toler                        = shakeParm.y;
            
            // Calculate unchanging quantities
            xij                                 = xi - xij;
            yij                                 = yi - yij;
            zij                                 = zi - zij;
            xik                                 = xi - xik;
            yik                                 = yi - yik;
            zik                                 = zi - zik; 
            XIL                                 = xi - XIL;
            YIL                                 = yi - YIL;
            ZIL                                 = zi - ZIL;
            xim                                 = xi - xim;
            yim                                 = yi - yim;
            zim                                 = zi - zim;               
           
            bool done                           = false;
            for (int i = 0; i < 3000; i++)
            {
                done                            = true;
                
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
                    
                        double acor             = diff / (rrpr * (double)2.0 * (INVMASSI + cSim.invMassH));
                        double h                = xij * acor;
                        xpi                    += h * INVMASSI;
                        xpj                    -= h * cSim.invMassH;
                        h                       = yij * acor;
                        ypi                    += h * INVMASSI;
                        ypj                    -= h * cSim.invMassH;
                        h                       = zij * acor;
                        zpi                    += h * INVMASSI;
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
                    
                        double acor             = diff / (rrpr * 2.0 * (INVMASSI + cSim.invMassH));
                        double h                = xik * acor;
                        xpi                    += h * INVMASSI;
                        xpk                    -= h * cSim.invMassH;
                        h                       = yik * acor;
                        ypi                    += h * INVMASSI;
                        ypk                    -= h * cSim.invMassH;
                        h                       = zik * acor;
                        zpi                    += h * INVMASSI;
                        zpk                    -= h * cSim.invMassH;             
                    }
                }
                
 
                xpxx                            = xpi - XPL;
                ypxx                            = ypi - YPL;
                zpxx                            = zpi - ZPL;        
                rpxx2                           = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;
          
                // Apply correction to third hydrogen
                diff                            = toler - rpxx2;
                if (abs(diff) >= toler * cSim.tol)
                {
                    done                        = false;
                  
                    // Shake resetting of coordinate is done here
                    double rrpr                 = XIL * xpxx + YIL * ypxx + ZIL * zpxx;     
                    if (rrpr >= toler * 1.0e-06)
                    {             
                        double acor             = diff / (rrpr * (double)2.0 * (INVMASSI + cSim.invMassH));
                        double h                = XIL * acor;
                        xpi                    += h * INVMASSI;
                        XPL                    -= h * cSim.invMassH;
                        h                       = YIL * acor;
                        ypi                    += h * INVMASSI;
                        YPL                    -= h * cSim.invMassH;
                        h                       = ZIL * acor;
                        zpi                    += h * INVMASSI;
                        ZPL                    -= h * cSim.invMassH;             
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
                        double acor             = diff / (rrpr * (double)2.0 * (INVMASSI + cSim.invMassH));
                        double h                = xim * acor;
                        xpi                    += h * INVMASSI;
                        xpm                    -= h * cSim.invMassH;
                        h                       = yim * acor;
                        ypi                    += h * INVMASSI;
                        ypm                    -= h * cSim.invMassH;
                        h                       = zim * acor;
                        zpi                    += h * INVMASSI;
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
                cSim.pImageX[shakeID2.x]        = xpj;
                cSim.pImageY[shakeID2.x]        = ypj;
                cSim.pImageZ[shakeID2.x]        = zpj;
                cSim.pImageX[shakeID2.y]        = xpk;
                cSim.pImageY[shakeID2.y]        = ypk;
                cSim.pImageZ[shakeID2.y]        = zpk;
                cSim.pImageX[shakeID2.z]        = XPL;
                cSim.pImageY[shakeID2.z]        = YPL;
                cSim.pImageZ[shakeID2.z]        = ZPL;
                cSim.pImageX[shakeID2.w]        = xpm;
                cSim.pImageY[shakeID2.w]        = ypm;
                cSim.pImageZ[shakeID2.w]        = zpm;
            }                                    
        }
    }
}

void kShakeInitKernels(gpuContext gpu)
{
    if (gpu->sm_version >= SM_3X)
    {
        cudaFuncSetCacheConfig(kShake_kernel, cudaFuncCachePreferL1);
        cudaFuncSetCacheConfig(kPMEShake_kernel, cudaFuncCachePreferL1);
    }
}

void kShake(gpuContext gpu)
{
    texref.normalized = 0;
    texref.filterMode = cudaFilterModePoint;
    texref.addressMode[0] = cudaAddressModeClamp;
    texref.channelDesc.x = 32;
    texref.channelDesc.y = 32;
    texref.channelDesc.z = 0;
    texref.channelDesc.w = 0;
    cudaBindTexture(NULL, texref, (int2*)(gpu->sim.pForce), gpu->sim.stride3 * sizeof(int2));
    unsigned int totalConstraints = gpu->sim.slowShakeOffset;
    unsigned int totalBlocks = (totalConstraints + 63) / 64;
    if (gpu->bNeighborList)
    {
        kPMEShake_kernel<<<totalBlocks, 64>>>();  
    }
    else
    {    
        kShake_kernel<<<totalBlocks, 64>>>();  
    }
    LAUNCHERROR("kShake");
    cudaUnbindTexture(texref);
}
