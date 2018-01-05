#ifndef _NMC_UCU_COMPLEX_H_
#define _NMC_UCU_COMPLEX_H_

#include <cuComplex.h>

__host__ __device__ static __inline__ cuFloatComplex ucuCaddf (cuFloatComplex x,
                                                              float y)
{
    return make_cuFloatComplex (cuCrealf(x) + y, 
                                cuCimagf(x));
}

__host__ __device__ static __inline__ cuFloatComplex ucuCmulf(float c,
															  cuFloatComplex x)
{
	return make_cuFloatComplex(c*cuCrealf(x), c*cuCimagf(x));
}

__host__ __device__ static __inline__ cuFloatComplex ucuCmulf(cuFloatComplex x,
															  float c)
{
	return make_cuFloatComplex(c*cuCrealf(x), c*cuCimagf(x));
}

__host__ __device__ static __inline__ float ucuCargf(cuFloatComplex cm)
{
	if(cuCabsf(cm) == 0.0) return 0.0;
	return atan2f(cuCimagf(cm), cuCrealf(cm));
}

__host__ __device__ static __inline__ cuFloatComplex ucuCpowf(cuFloatComplex cm, int n)
{
	float t = ucuCargf( cm ) * n;
	float w = powf(cuCabsf(cm), n);
	return make_cuFloatComplex(w * cosf(t), w * sinf(t));
}

__host__ __device__ static __inline__ cuFloatComplex ucuClogf(cuFloatComplex x)
{
	return make_cuFloatComplex(logf(cuCabsf(x)), ucuCargf(x));
}

std::ostream& operator<<(std::ostream& os, const cuFloatComplex& x){
	return os << "(" << cuCrealf(x) << ")" << " + i(" << cuCimagf(x) << ")";
}

#endif //_NMC_UCU_COMPLEX_H_