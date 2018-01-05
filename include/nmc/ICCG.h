#ifndef _NMC_ICCG_H_
#define _NMC_ICCG_H_

#include "nmc/Matrix.h"
#include "nmc/Vector.h"
#include "nmc/Complex.h"

namespace nmc {

template <class VAL_TYPE>
class ICCG{
public:
	/*!
		˜A—§ˆêŽŸ•û’öŽ®Ax=b‚ð‰ð‚­
		@param[in] src
		@param[in] b
		@param[out] x
	*/
	static bool Solve(const Matrix<VAL_TYPE>& src, const Vector<VAL_TYPE>& b,
		Vector<VAL_TYPE>& x);
};

template <class VAL_TYPE>
bool ICCG::Solve(
	const Matrix<VAL_TYPE> &src,
	const Vector<VAL_TYPE> &b,
	Vector<VAL_TYPE> &x)
{
	assert(src.row() == src.col());
	assert(src.row() == b.dim());
	assert(src.row() == x.dim());

	Vector<VAL_TYPE> r0_vec;
	r0_vec = 
}

};

#endif //_NMC_ICCG_H_