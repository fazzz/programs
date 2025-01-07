#ifndef PTXMACROS_H
#define PTXMACROS_H
__device__ inline unsigned long long int llitoulli(long long int l)
{
    unsigned long long int u;
    asm("mov.b64    %0, %1;" : "=l"(u) : "l"(l));
    return u;
}

__device__ inline long long int ullitolli(unsigned long long int u)
{
    long long int l;
    asm("mov.b64    %0, %1;" : "=l"(l) : "l"(u));
    return l;
}

__device__ inline long long int lliroundf(float f)
{
    long long int l;
    asm("cvt.rni.s64.f32 	%0, %1;" : "=l"(l) : "f"(f));
    return l;
}

__device__ inline long long int lliroundd(double d)
{
    long long int l;
    asm("cvt.rni.s64.f64 	%0, %1;" : "=l"(l) : "d"(d));
    return l;
}

__device__ inline double float2todouble(float2 f)
{
    double d;
    asm("mov.b64 	%0, {%1, %2};" : "=d"(d) : "f"(f.x), "f"(f.y));
    return d;
}

__device__ inline float2 doubletofloat2(double d)
{
    float2 f;
    asm("mov.b64    {%0, %1}, %2;" : "=f"(f.x), "=f"(f.y) : "d"(d));
    return f;
}

__device__ inline unsigned int reduce(unsigned int x)
{
    unsigned int s = x;
    asm("shfl.up %0, %1, 0x1, 0x0;\n\t"
     "   add.u32 %1, %0, %1;\n\t"
        "shfl.up %0, %1, 0x2, 0x0;\n\t"
     "   add.u32 %1, %0, %1;\n\t"
        "shfl.up %0, %1, 0x4, 0x0;\n\t"
     "   add.u32 %1, %0, %1;\n\t"      
        "shfl.up %0, %1, 0x8, 0x0;\n\t"
     "   add.u32 %1, %0, %1;\n\t"
        "shfl.up %0, %1, 0x16, 0x0;\n\t"
     "   add.u32 %1, %0, %1;" 
     : "=r"(s) : "r"(x));
    return s;
}

static __device__ __inline__ unsigned int __shfl(unsigned int var, unsigned int srcLane, int width=32) {
	unsigned int ret, c;
	c = ((32-width) << 8) | 0x1f;
	asm volatile ("shfl.idx.b32 %0, %1, %2, %3;" : "=r"(ret) : "r"(var), "r"(srcLane), "r"(c));
	return ret;
}

static __device__ __inline__ unsigned int __shfl(unsigned int var, int srcLane, int width=32) {
	unsigned int ret, c;
	c = ((32-width) << 8) | 0x1f;
	asm volatile ("shfl.idx.b32 %0, %1, %2, %3;" : "=r"(ret) : "r"(var), "r"(srcLane), "r"(c));
	return ret;
}

static __device__ __inline__ float __shfl(float var, unsigned int srcLane, int width=32) {
	float ret;
    unsigned int c;
	c = ((32-width) << 8) | 0x1f;
	asm volatile ("shfl.idx.b32 %0, %1, %2, %3;" : "=f"(ret) : "f"(var), "r"(srcLane), "r"(c));
	return ret;
}

static __device__ __inline__ double __shfl(double var, unsigned int srcLane, int width=32) {
	double ret;
    unsigned int c;
	c = ((32-width) << 8) | 0x1f;
	asm volatile ("{\n\t"
                  ".reg .u32 dlo;\n\t"
                  ".reg .u32 dhi;\n\t"
                  " mov.b64 {dlo, dhi}, %1;\n\t"
                  " shfl.idx.b32 dlo, dlo, %2, %3;\n\t"
                  " shfl.idx.b32 dhi, dhi, %2, %3;\n\t"
                  " mov.b64 %0, {dlo, dhi};\n\t"
                  "}"
                  : "=d"(ret) : "d"(var), "r"(srcLane), "r"(c)); 
	return ret;
}

static __device__ __inline__ double __shfl(double var, int srcLane, int width=32) {
	double ret;
    unsigned int c;
	c = ((32-width) << 8) | 0x1f;
	asm volatile ("{\n\t"
                  ".reg .u32 dlo;\n\t"
                  ".reg .u32 dhi;\n\t"
                  " mov.b64 {dlo, dhi}, %1;\n\t"
                  " shfl.idx.b32 dlo, dlo, %2, %3;\n\t"
                  " shfl.idx.b32 dhi, dhi, %2, %3;\n\t"
                  " mov.b64 %0, {dlo, dhi};\n\t"
                  "}"
                  : "=d"(ret) : "d"(var), "r"(srcLane), "r"(c)); 
	return ret;
}

#if 0
static __device__ __inline__ long long int __shfl(long long int var, int srcLane, int width=32) {
	long long int ret;
    unsigned int c;
	c = ((32-width) << 8) | 0x1f;
	asm volatile ("{\n\t"
                  ".reg .u32 dlo;\n\t"
                  ".reg .u32 dhi;\n\t"
                  " mov.b64 {dlo, dhi}, %1;\n\t"
                  " shfl.idx.b32 dlo, dlo, %2, %3;\n\t"
                  " shfl.idx.b32 dhi, dhi, %2, %3;\n\t"
                  " mov.b64 %0, {dlo, dhi};\n\t"
                  "}"
                  : "=l"(ret) : "l"(var), "r"(srcLane), "r"(c)); 
	return ret;
}
#endif

#endif
