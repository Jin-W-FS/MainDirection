#ifndef REALNUM_WRAPPER_H
#define REALNUM_WRAPPER_H

#if defined(REALNUM_USE_FLOAT)

#include <math.h>

typedef float real;

#define int2real(x)		((float)(x))
#define real2int(x)		((int)(x))
#define realdecimal(x)	((x) - (int)x)

#define float2real(x)	(x)
#define real2float(x)   (x)
#define real(x)			((float)(x))

#define fixed2real(x, N)	((float)(x) / (float)(1 << N))

#define realabs(x)		(fabsf(x))
#define realround(x)	(roundf(x))

#define realmul(x, y)	((x) * (y))
#define realdiv(x, y)	((x) / (y))

#define realsqrt(x)		(sqrtf(x))

#elif defined(REALNUM_USE_FIXED)

#include "fixedpoint.h"

typedef fixed real;

#define int2real(x)		int2fixed(x)
#define real2int(x)		fixed2int(x)
#define realdecimal(x)	fixdecimal(x)

#define float2real(x)	float2fixed(x)
#define real2float(x)   fixed2float(x)
#define real(x)			fixed(x)

static inline real fixed2real(long x, int N) {
	int d = FIXEDLEN - N;
	return (d >= 0) ? ((BASETYPE)x << d) : ((BASETYPE)x >> -d);
}

#define realabs(x)		fixabs(x)
#define realround(x)	fixround(x)

#define realmul(x, y)	fixmul(x, y)
#define realdiv(x, y)	fixdiv(x, y)

#define realsqrt(x)		fixsqrt(x)

#endif

#endif