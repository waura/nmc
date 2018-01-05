#ifndef _NMC_BICGSTAB_H_
#define _NMC_BICGSTAB_H_

#include "nmc/Matrix.h"
#include "nmc/Vector.h"
#include "nmc/Complex.h"

namespace nmc{

template <class VAL_TYPE>
class BiCGSTAB{
public:
	/*!
		˜A—§ˆêŽŸ•û’öŽ®Ax=b‚ð‰ð‚­
		@param[in,out] conv_ration
		@param[in,out] iteration number of max loop
		@param[in] src matrix
		@param[in] r_vec = b - Ax_{0}
		@param[out] x_vec
	*/
	static bool Solve(double& conv_ratio, unsigned int& iteration, 
		const Matrix<VAL_TYPE>& src, Vector<VAL_TYPE>& r_vec, Vector<VAL_TYPE>& x);
};

template <class VAL_TYPE>
bool BiCGSTAB<VAL_TYPE>::Solve(
	double& conv_ratio,
	unsigned int& iteration,
	const Matrix<VAL_TYPE> &src,
	Vector<VAL_TYPE> &r_vec,
	Vector<VAL_TYPE> &x_vec)
{
	assert(src.row() == src.col());
	assert(src.row() == r_vec.dim());
	assert(src.row() == x_vec.dim());

	int i,j;
	double conv_ratio_tol = conv_ratio;
	Vector<VAL_TYPE>  s_vec( r_vec.dim() );
	Vector<VAL_TYPE> As_vec( r_vec.dim() );
	Vector<VAL_TYPE>  p_vec( r_vec.dim() );
	Vector<VAL_TYPE> Ap_vec( r_vec.dim() );
	Vector<VAL_TYPE> r0_conjugate_vec( r_vec.dim() );
	
	VAL_TYPE sq_inv_norm_res;
	{
		VAL_TYPE tmp=0.0;
		for(j=0; j<r_vec.dim(); j++){
			tmp += r_vec[j] * r_vec[j];
		}
		if( tmp.Abs() < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;
		}
		sq_inv_norm_res = static_cast<VAL_TYPE>(1.0) / tmp;
	}

	//
	r0_conjugate_vec = r_vec;
	//p_{0} = r_{0}
	p_vec = r_vec;

	for(i=1; i<iteration; i++){
		//calc (r, r0*)
		VAL_TYPE r_r0conj=0.0;
		for(j=0; j<r_vec.dim(); j++){
			r_r0conj += r_vec[j] * r0_conjugate_vec[j];
		}

		//calc Ap
		Ap_vec = src * p_vec;

		//alpha = (r, r0*)/(Ap, ro*)
		VAL_TYPE alpha;
		{
			VAL_TYPE dnm = 0.0;
			for(j=0; j<Ap_vec.dim(); j++){
				dnm += Ap_vec[j] * r0_conjugate_vec[j];
			}
			alpha = r_r0conj / dnm;
		}

		//s = r - alpha*Ap
		s_vec = r_vec;
		for(j=0; j<s_vec.dim(); j++){
			s_vec[j] -= alpha * Ap_vec[j]; 
		}

		//calc As
		As_vec = src * s_vec;

		//omega = (s, As)/(As, As)
		VAL_TYPE omega;
		{
			VAL_TYPE numr = 0.0;
			for(j=0; j<s_vec.dim(); j++){
				numr += s_vec[j] * As_vec[j];
			}
			VAL_TYPE dnm = 0.0;
			for(j=0; j<As_vec.dim(); j++){
				dnm += As_vec[j] * As_vec[j];
			}
			omega = numr / dnm;
		}

		//new x = x + alpha * p + omega * s
		for(j=0; j<x_vec.dim(); j++){
			x_vec[j] += alpha * p_vec[j] + omega * s_vec[j];
		}

		//new r = s - omega * As
		for(j=0; j<r_vec.dim(); j++){
			r_vec[j] = s_vec[j] - omega * As_vec[j];
		}

		{
			VAL_TYPE sq_norm_res = 0.0;
			for(j=0; j<r_vec.dim(); j++){
				sq_norm_res += r_vec[j] * r_vec[j];
			}
			VAL_TYPE sq_conv_ratio = sq_norm_res * sq_inv_norm_res;
			if(sq_conv_ratio < conv_ratio_tol*conv_ratio_tol){
				conv_ratio = sqrt( sq_conv_ratio.Abs() );
				iteration = i;
				return true;
			}
		}

		//beta = {(new r, r0*)/(old r, r0*)}*{(alpha/omega)} 
		VAL_TYPE beta;
		{
			VAL_TYPE numr = 0.0;
			for(j=0; j<r_vec.dim(); j++){
				numr += r_vec[j] * r0_conjugate_vec[j];
			}
			beta = numr * alpha / (r_r0conj * omega);
		}

		//p = r + beta(p - omega * Ap)
		for(j=0; j<p_vec.dim(); j++){
			p_vec[j] *= beta;
			p_vec[j] += r_vec[j];
			p_vec[j] -= beta * omega * Ap_vec[j];
		}
	}

	return true;
}

typedef BiCGSTAB< Double > BiCGSTABd;
typedef BiCGSTAB< Complexd > BiCGSTABcmd;

};

#endif //_NMC_BICGSTAB_H_