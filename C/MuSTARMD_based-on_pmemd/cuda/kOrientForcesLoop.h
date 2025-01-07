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

#ifdef EP_NEIGHBORLIST
#define PATOMX(p) cSim.pImageX[p]
#define PATOMY(p) cSim.pImageY[p]
#define PATOMZ(p) cSim.pImageZ[p]
#if defined(EP_ONEPOINT)
#define EP1 index  
#if defined(EP_TYPE1)
            int index                       = cSim.pImageExtraPoint11Index[pos]; 
            int4 frame                      = cSim.pImageExtraPoint11Frame[pos];        
#elif defined(EP_TYPE2)
            int index                       = cSim.pImageExtraPoint12Index[pos]; 
            int4 frame                      = cSim.pImageExtraPoint12Frame[pos];
#endif
#elif defined(EP_TWOPOINTS)
#define EP1 index.x
#if defined(EP_TYPE1) 
            int2 index                      = cSim.pImageExtraPoint21Index[pos]; 
            int4 frame                      = cSim.pImageExtraPoint21Frame[pos];          
#elif defined(EP_TYPE2)
            int2 index                      = cSim.pImageExtraPoint22Index[pos]; 
            int4 frame                      = cSim.pImageExtraPoint22Frame[pos];        
#endif
#endif  
#else
#define PATOMX(p) cSim.pAtomX[p]
#define PATOMY(p) cSim.pAtomY[p]
#define PATOMZ(p) cSim.pAtomZ[p]
#if defined(EP_ONEPOINT)
#define EP1 index  
#if defined(EP_TYPE1)
            int index                       = cSim.pExtraPoint11Index[pos]; 
            int4 frame                      = cSim.pExtraPoint11Frame[pos];        
#elif defined(EP_TYPE2)
            int index                       = cSim.pExtraPoint12Index[pos]; 
            int4 frame                      = cSim.pExtraPoint12Frame[pos];
#endif
#elif defined(EP_TWOPOINTS)
#define EP1 index.x
#if defined(EP_TYPE1) 
            int2 index                      = cSim.pExtraPoint21Index[pos]; 
            int4 frame                      = cSim.pExtraPoint21Frame[pos];          
#elif defined(EP_TYPE2)
            int2 index                      = cSim.pExtraPoint22Index[pos]; 
            int4 frame                      = cSim.pExtraPoint22Frame[pos];        
#endif
#endif  
#endif
            // Handle 1st or only EP                     
            PMEDouble x1                    = PATOMX(EP1);
            PMEDouble y1                    = PATOMY(EP1);
            PMEDouble z1                    = PATOMZ(EP1);
            PMEDouble xp                    = PATOMX(frame.x);
            PMEDouble yp                    = PATOMY(frame.x);
            PMEDouble zp                    = PATOMZ(frame.x);               
            PMEDouble forceX                = cSim.pForceX[EP1];
            PMEDouble forceY                = cSim.pForceY[EP1];
            PMEDouble forceZ                = cSim.pForceZ[EP1];
            PMEDouble rx                    = x1 - xp;
            PMEDouble ry                    = y1 - yp;
            PMEDouble rz                    = z1 - zp;
            cSim.pForceX[EP1]               = (PMEDouble)0.0;
            cSim.pForceY[EP1]               = (PMEDouble)0.0;
            cSim.pForceZ[EP1]               = (PMEDouble)0.0;       
#ifdef EP_VIRIAL 
            cSim.pNBForceX[EP1]             = (PMEDouble)0.0;
            cSim.pNBForceY[EP1]             = (PMEDouble)0.0;
            cSim.pNBForceZ[EP1]             = (PMEDouble)0.0; 
            v11                            += forceX * rx;
            v22                            += forceY * ry;
            v33                            += forceZ * rz;
#endif
            PMEDouble torqueX               = ry * forceZ - rz * forceY;            
            PMEDouble torqueY               = rz * forceX - rx * forceZ;            
            PMEDouble torqueZ               = rx * forceY - ry * forceX;             
#ifdef EP_TWOPOINTS
            // Handle 2nd EP
            PMEDouble x2                    = PATOMX(index.y);
            PMEDouble y2                    = PATOMY(index.y);
            PMEDouble z2                    = PATOMZ(index.y);         
            PMEDouble forceX2               = cSim.pForceX[index.y];
            PMEDouble forceY2               = cSim.pForceY[index.y];
            PMEDouble forceZ2               = cSim.pForceZ[index.y];   
            rx                              = x2 - xp;
            ry                              = y2 - yp;
            rz                              = z2 - zp;
            cSim.pForceX[index.y]           = (PMEDouble)0.0;
            cSim.pForceY[index.y]           = (PMEDouble)0.0;
            cSim.pForceZ[index.y]           = (PMEDouble)0.0;
#ifdef EP_VIRIAL
            cSim.pNBForceX[index.y]         = (PMEDouble)0.0;
            cSim.pNBForceY[index.y]         = (PMEDouble)0.0;
            cSim.pNBForceZ[index.y]         = (PMEDouble)0.0;
            v11                            += forceX2 * rx;
            v22                            += forceY2 * ry;
            v33                            += forceZ2 * rz;                                             
#endif
            forceX                         += forceX2;
            forceY                         += forceY2;
            forceZ                         += forceZ2;            
            torqueX                        += ry * forceZ2 - rz * forceY2;            
            torqueY                        += rz * forceX2 - rx * forceZ2;            
            torqueZ                        += rx * forceY2 - ry * forceX2;         
#endif   

#ifdef EP_TYPE1
            PMEDouble apx                   = PATOMX(frame.y);
            PMEDouble apy                   = PATOMY(frame.y);
            PMEDouble apz                   = PATOMZ(frame.y);
            PMEDouble bpx                   = PATOMX(frame.z);
            PMEDouble bpy                   = PATOMY(frame.z);
            PMEDouble bpz                   = PATOMZ(frame.z);
            PMEDouble cpx                   = PATOMX(frame.w);
            PMEDouble cpy                   = PATOMY(frame.w);
            PMEDouble cpz                   = PATOMZ(frame.w);
#elif defined(EP_TYPE2)
            PMEDouble apx                   = PATOMX(frame.y);
            PMEDouble apy                   = PATOMY(frame.y);
            PMEDouble apz                   = PATOMZ(frame.y);
            PMEDouble px                    = PATOMX(frame.z);
            PMEDouble py                    = PATOMY(frame.z);
            PMEDouble pz                    = PATOMZ(frame.z);
            PMEDouble cpx                   = PATOMX(frame.w);
            PMEDouble cpy                   = PATOMY(frame.w);
            PMEDouble cpz                   = PATOMZ(frame.w);
            PMEDouble bpx                   = PATOMX(frame.x);
            PMEDouble bpy                   = PATOMY(frame.x);
            PMEDouble bpz                   = PATOMZ(frame.x);
            apx                             = (PMEDouble)0.5 * (apx + px);
            apy                             = (PMEDouble)0.5 * (apy + py);
            apz                             = (PMEDouble)0.5 * (apz + pz);
            cpx                             = (PMEDouble)0.5 * (cpx + px);
            cpy                             = (PMEDouble)0.5 * (cpy + py);
            cpz                             = (PMEDouble)0.5 * (cpz + pz);           
#endif            
            PMEDouble ux                    = apx - bpx;
            PMEDouble uy                    = apy - bpy;
            PMEDouble uz                    = apz - bpz;
            PMEDouble usiz                  = rsqrt(ux * ux + uy * uy + uz * uz);
            PMEDouble vx                    = cpx - bpx;
            PMEDouble vy                    = cpy - bpy;
            PMEDouble vz                    = cpz - bpz;
            PMEDouble vsiz                  = rsqrt(vx * vx + vy * vy + vz * vz);
            PMEDouble wx                    = uy * vz - uz * vy;
            PMEDouble wy                    = uz * vx - ux * vz;
            PMEDouble wz                    = ux * vy - uy * vx;
            ux                             *= usiz;
            uy                             *= usiz;
            uz                             *= usiz;  
            vx                             *= vsiz;
            vy                             *= vsiz;
            vz                             *= vsiz;                      
            PMEDouble wsiz                  = rsqrt(wx * wx + wy * wy + wz * wz);
            wx                             *= wsiz;
            wy                             *= wsiz;
            wz                             *= wsiz;
            PMEDouble dx                    = vx - ux;
            PMEDouble dy                    = vy - uy;
            PMEDouble dz                    = vz - uz;
            PMEDouble dotdu                 = ux * dx + uy * dy + uz * dz;
            PMEDouble dotdv                 = vx * dx + vy * dy + vz * dz;
            PMEDouble upx                   = dx - dotdu * ux;
            PMEDouble upy                   = dy - dotdu * uy;
            PMEDouble upz                   = dz - dotdu * uz;
            PMEDouble vpx                   = dx - dotdv * vx;
            PMEDouble vpy                   = dy - dotdv * vy;
            PMEDouble vpz                   = dz - dotdv * vz;
            PMEDouble upsiz                 = rsqrt(upx * upx + upy * upy + upz * upz);
            upx                            *= upsiz;
            upy                            *= upsiz;
            upz                            *= upsiz;
            PMEDouble vpsiz                 = rsqrt(vpx * vpx + vpy * vpy + vpz * vpz);
            vpx                            *= vpsiz;
            vpy                            *= vpsiz;
            vpz                            *= vpsiz;              
            PMEDouble c                     = ux * vx + uy * vy + uz * vz;
            PMEDouble s                     = rsqrt((PMEDouble)1.0 - c * c);
            PMEDouble uvdis                 = usiz * s;
            PMEDouble vudis                 = vsiz * s;           
            PMEDouble dphidu                = -(torqueX * ux + torqueY * uy + torqueZ * uz);
            PMEDouble dphidv                = -(torqueX * vx + torqueY * vy + torqueZ * vz);
            PMEDouble dphidw                = -(torqueX * wx + torqueY * wy + torqueZ * wz);                      
            dphidv                         *= uvdis;
            dphidu                         *= vudis;
            usiz                           *= (PMEDouble)0.5 * dphidw;
            vsiz                           *= (PMEDouble)0.5 * dphidw;
            PMEDouble dux                   = -wx * dphidv + upx * usiz;
            PMEDouble duy                   = -wy * dphidv + upy * usiz;
            PMEDouble duz                   = -wz * dphidv + upz * usiz;
            PMEDouble dvx                   =  wx * dphidu + vpx * vsiz;
            PMEDouble dvy                   =  wy * dphidu + vpy * vsiz;
            PMEDouble dvz                   =  wz * dphidu + vpz * vsiz;  
        
#ifdef EP_VIRIAL
            // Get torque contribution to virial:
            v11                            += dux * (apx - bpx) + dvx * (cpx - bpx);
            v22                            += duy * (apy - bpy) + dvy * (cpy - bpy);
            v33                            += duz * (apz - bpz) + dvz * (cpz - bpz);
#endif
                                      
#ifdef EP_TYPE1    
            cSim.pForceX[frame.y]          -= dux;
            cSim.pForceY[frame.y]          -= duy;
            cSim.pForceZ[frame.y]          -= duz;
            cSim.pForceX[frame.w]          -= dvx;
            cSim.pForceY[frame.w]          -= dvy;
            cSim.pForceZ[frame.w]          -= dvz; 
            cSim.pForceX[frame.z]          += dvx + dux + forceX;
            cSim.pForceY[frame.z]          += dvy + duy + forceY;
            cSim.pForceZ[frame.z]          += dvz + duz + forceZ;
#ifdef EP_VIRIAL
            cSim.pNBForceX[frame.y]        -= dux;
            cSim.pNBForceY[frame.y]        -= duy;
            cSim.pNBForceZ[frame.y]        -= duz;
            cSim.pNBForceX[frame.w]        -= dvx;
            cSim.pNBForceY[frame.w]        -= dvy;
            cSim.pNBForceZ[frame.w]        -= dvz; 
            cSim.pNBForceX[frame.z]        += dvx + dux + forceX;
            cSim.pNBForceY[frame.z]        += dvy + duy + forceY;
            cSim.pNBForceZ[frame.z]        += dvz + duz + forceZ;
#endif                           
#elif defined(EP_TYPE2)
            cSim.pForceX[frame.x]          += dvx + dux + forceX;
            cSim.pForceY[frame.x]          += dvy + duy + forceY;
            cSim.pForceZ[frame.x]          += dvz + duz + forceZ;      
            cSim.pForceX[frame.y]          -= (PMEDouble)0.5 * dux;
            cSim.pForceY[frame.y]          -= (PMEDouble)0.5 * duy;
            cSim.pForceZ[frame.y]          -= (PMEDouble)0.5 * duz;
            cSim.pForceX[frame.w]          -= (PMEDouble)0.5 * dvx;
            cSim.pForceY[frame.w]          -= (PMEDouble)0.5 * dvy;
            cSim.pForceZ[frame.w]          -= (PMEDouble)0.5 * dvz; 
            cSim.pForceX[frame.z]          -= (PMEDouble)0.5 * (dvx + dux);
            cSim.pForceY[frame.z]          -= (PMEDouble)0.5 * (dvy + duy);
            cSim.pForceZ[frame.z]          -= (PMEDouble)0.5 * (dvz + duz);  
#ifdef EP_VIRIAL            
            cSim.pNBForceX[frame.x]        += dvx + dux + forceX;
            cSim.pNBForceY[frame.x]        += dvy + duy + forceY;
            cSim.pNBForceZ[frame.x]        += dvz + duz + forceZ;      
            cSim.pNBForceX[frame.y]        -= (PMEDouble)0.5 * dux;
            cSim.pNBForceY[frame.y]        -= (PMEDouble)0.5 * duy;
            cSim.pNBForceZ[frame.y]        -= (PMEDouble)0.5 * duz;
            cSim.pNBForceX[frame.w]        -= (PMEDouble)0.5 * dvx;
            cSim.pNBForceY[frame.w]        -= (PMEDouble)0.5 * dvy;
            cSim.pNBForceZ[frame.w]        -= (PMEDouble)0.5 * dvz; 
            cSim.pNBForceX[frame.z]        -= (PMEDouble)0.5 * (dvx + dux);
            cSim.pNBForceY[frame.z]        -= (PMEDouble)0.5 * (dvy + duy);
            cSim.pNBForceZ[frame.z]        -= (PMEDouble)0.5 * (dvz + duz); 
#endif                    
#endif
#undef EP1
#undef PATOMX
#undef PATOMY
#undef PATOMZ
