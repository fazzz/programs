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


/*
 * This special version of sincos is designed for |a| < 6*PI. On a GTX 285
 * it is about 25% faster than sincos from the CUDA math library. Also uses
 * 8 fewer registers than the CUDA math library's sincos. Maximum observed
 * error is 2 ulps across range stated above. Infinities and negative zero
 * are not handled according to C99 specifications. NaNs are handled fine.
 */
//__device__ void faster_sincos2(double a, double *sptr, double *cptr) 
{
  double t, u, s, c, j, a2;
  int i;

  i = __double2int_rn (a * 6.3661977236758138e-1);
  j = (double)i;
  a = __fma_rn (-j, 1.57079632679489660e+000, a); /* PIO2_HI */
  a = __fma_rn (-j, 6.12323399573676600e-017, a); /* PIO2_LO */
  a2 = a * a;
  u =                  -1.136788825395985E-011;   
  u = __fma_rn (u, a2,  2.087588480545065E-009);
  u = __fma_rn (u, a2, -2.755731555403950E-007);
  u = __fma_rn (u, a2,  2.480158729365970E-005);
  u = __fma_rn (u, a2, -1.388888888888074E-003);
  u = __fma_rn (u, a2,  4.166666666666664E-002);
  u = __fma_rn (u, a2, -5.000000000000000E-001);
  u = __fma_rn (u, a2,  1.000000000000000E+000);
  t =                   1.5896230157221844E-010;
  t = __fma_rn (t, a2, -2.5050747762850355E-008);
  t = __fma_rn (t, a2,  2.7557313621385676E-006);
  t = __fma_rn (t, a2, -1.9841269829589539E-004);
  t = __fma_rn (t, a2,  8.3333333333221182E-003);
  t = __fma_rn (t, a2, -1.6666666666666630E-001);
  t = t * a2;
  t = __fma_rn(t, a, a);
  if (i & 1) {
    s = u;
    c = t;
  } else {
    s = t;
    c = u;
  }
  if (i & 2) {
    s = -s;
  }
  i++;
  if (i & 2) {
    c = -c;
  }
  *sptr = s;
  *cptr = c;
}

