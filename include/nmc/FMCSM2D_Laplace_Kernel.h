#ifndef _FMCSM2D_LAPLACE_KERNEL_H_
#define _FMCSM2D_LAPLACE_KERNEL_H_

#include <cmath>

#include "nmc/def.h"
#include "nmc/FMCSM2DKernel.h"

namespace nmc
{

class FMCSM2D_Laplace_Kernel : public FMCSM2DKernel<Complexd>
{
public:

	void Init(FMCSM2DTree<Complexd>* tree);
	inline void ElementMPCefficients(FMCSM2DTreeNode<Complexd>* node,
		const Vector<Complexd>& c_vec);
	inline void M2M_Translation(FMCSM2DTreeNode<Complexd>* node);
	inline void L2L_Translation(FMCSM2DTreeNode<Complexd>* node);
	inline void M2L_Translation(const FMCSM2DTreeNode<Complexd>* src,
		FMCSM2DTreeNode<Complexd>* dst);
	inline void EvalMultipole(const FMCSM2DTreeNode<Complexd>* node,
		Vector<Complexd>& dst);
	inline void EvalDirect(const FMMPos& from_pos, const FMMPos& to_pos,
		const Vector<Complexd>& cf, Vector<Complexd>& dst);
};


static inline double binomial_coefficient(int n, int m)
{
	assert(n >= 0);
	assert(m >= 0);

	if(n == m) return 1.0;
	if(m == 0) return 1.0;

	int i;
	double nn=1.0, nm=1.0, mm=1.0;
	for(i=1; i<=n; i++){
		nn *= i;
	}
	for(i=1; i<=(n-m); i++){
		nm *= i;
	}
	for(i=1; i<=m; i++){
		mm *= i;
	}

	return nn /(nm*mm);
}

inline void FMCSM2D_Laplace_Kernel::Init(
	FMCSM2DTree<Complexd>* tree)
{
}

inline void FMCSM2D_Laplace_Kernel::ElementMPCefficients(
	FMCSM2DTreeNode<Complexd> *node,
	const Vector<Complexd>& c_vec)
{
	assert(node->m_pSourcePos);

	int s,t;
	double x,y;
	double phi;
	double rho;
	double J_s;
	Complexd c;
	Complexd z;

	int p = node->m_Multipole.dim()-1;
	for(t=0; t<node->m_pSourcePos->dim(); t++){
		x = node->m_pSourcePos->GetVal(t).val[X_KEY] - node->m_CenterPos[X_KEY];
		y = node->m_pSourcePos->GetVal(t).val[Y_KEY] - node->m_CenterPos[Y_KEY];
		z.SetNum(x, y);

		c = c_vec[ node->m_pSourcePos->GetVal(t).index ];

		//M_{0} += c
		node->m_Multipole[0] += c;
		for(s=1; s<=p; s++){
			//M_{s} += -c * (z^s) / s		
			node->m_Multipole[s] += (-1.0/s) * c * pow(z, s);
		}

		/*
		std::cout << z << std::endl;
		std::cout << c << std::endl;
		std::cout << node->m_Multipole << std::endl;
		*/
	}
}

inline void FMCSM2D_Laplace_Kernel::M2M_Translation(
	FMCSM2DTreeNode<Complexd> *node)
{
	int s,t;
	int child_index;
	double x,y;
	Complexd z;

	int p = node->m_Multipole.dim()-1;
	for(child_index=0; child_index<4; child_index++){
		if(node->m_pNextTreeNodes[child_index]){
			x = node->m_pNextTreeNodes[child_index]->m_CenterPos[X_KEY] - node->m_CenterPos[X_KEY];
			y = node->m_pNextTreeNodes[child_index]->m_CenterPos[Y_KEY] - node->m_CenterPos[Y_KEY];
			z.SetNum(x, y);

			//M_{0] += M_0
			node->m_Multipole[0] += node->m_pNextTreeNodes[child_index]->m_Multipole[0];
			for(s=1; s<=p; s++){
				//M_{s} += -(M_0 * z^s) / s
				node->m_Multipole[s] += (-1.0/s) * node->m_pNextTreeNodes[child_index]->m_Multipole[0] * pow(z, s);
				for(t=1; t<=s; t++){
					//M_{s} += M_t * z^(s-t) * C(s-1, t-1)
					node->m_Multipole[s] += node->m_pNextTreeNodes[child_index]->m_Multipole[t] *
						pow(z, s-t) * binomial_coefficient(s-1, t-1);
				}
			}

		}
	}

	//
	//std::cout << node->m_Multipole << std::endl;
}

inline void FMCSM2D_Laplace_Kernel::L2L_Translation(
	FMCSM2DTreeNode<Complexd> *node)
{
	int s,t;
	int child_index;
	double x,y;
	Complexd z;

	//std::cout << node->m_Local << std::endl;

	int p = node->m_Local.dim()-1;
	for(child_index=0; child_index<4; child_index++){
		if(node->m_pNextTreeNodes[child_index] &&
			node->m_pNextTreeNodes[child_index]->m_isBundaryPosInNode)
		{
			x = node->m_pNextTreeNodes[child_index]->m_CenterPos[X_KEY] - node->m_CenterPos[X_KEY];
			y = node->m_pNextTreeNodes[child_index]->m_CenterPos[Y_KEY] - node->m_CenterPos[Y_KEY];
			z.SetNum(x, y);

			for(s=0; s<=p; s++){
				for(t=s; t<=p; t++){
					//L_{s} += L_{t} * C(t, s) * z^(t - s)
					node->m_pNextTreeNodes[child_index]->m_Local[s] += node->m_Local[t] * binomial_coefficient(t, s) * pow(z, t-s);
				}
			}

			//std::cout << node->m_pNextTreeNodes[child_index]->m_Local << std::endl;
		}
	}
}

inline void FMCSM2D_Laplace_Kernel::M2L_Translation(
	const FMCSM2DTreeNode<Complexd> *src,
	FMCSM2DTreeNode<Complexd> *dst)
{
	int s,t;
	double sign;
	Complexd z;

	double x = src->m_CenterPos[X_KEY] - dst->m_CenterPos[X_KEY];
	double y = src->m_CenterPos[Y_KEY] - dst->m_CenterPos[Y_KEY];
	z.SetNum(x, y);

	int p = dst->m_Local.dim()-1;

	//L_{0} += M_{0} * log( -z );
	dst->m_Local[0] += src->m_Multipole[0] * std::log( nmc::abs(-z) );
	for(t=1; t<=p; t++){
		sign = ((t % 2) == 0) ? 1.0 : -1.0;
		//L_{0} += (-1)^t * M_{t} / z^t
		dst->m_Local[0] += (sign * src->m_Multipole[t]) / pow(z, t);
	}

	for(s=1; s<=p; s++){
		//L_{s} += - M_{0} / (s * z^s)
		dst->m_Local[s] +=  -src->m_Multipole[0] / (s * pow(z, s));
		for(t=1; t<=p; t++){
			//L_{s} += (-1)^t * M_{t} * C(s+t-1, t-1) / z^(s+t)
			sign = ((t % 2) == 0) ? 1.0 : -1.0;
			dst->m_Local[s] += (sign * src->m_Multipole[t] * binomial_coefficient(s+t-1, t-1)) / pow(z, s+t);
		}
	}

	//std::cout << src->m_Multipole << std::endl;
	//std::cout << dst->m_Local << std::endl;
}

inline void FMCSM2D_Laplace_Kernel::EvalMultipole(
	const FMCSM2DTreeNode<Complexd> *node,
	Vector<Complexd> &dst)
{
	assert(node->m_pBundaryPos);

	int s,t;
	int index;
	double x,y;
	Complexd tmp;
	Complexd z;

	int p = node->m_Local.dim()-1;
	for(t=0; t<node->m_pBundaryPos->dim(); t++){
		index = node->m_pBundaryPos->GetVal(t).index;
		x = node->m_pBundaryPos->GetVal(t).val[X_KEY] - node->m_CenterPos[X_KEY];
		y = node->m_pBundaryPos->GetVal(t).val[Y_KEY] - node->m_CenterPos[Y_KEY];
		z.SetNum(x, y);

		tmp = 0.0;
		for(s=0; s<=p; s++){
			//L_{s} * z^s
			tmp += node->m_Local[s] * pow(z, s);
		}
		dst[index] += tmp.Re();
	}

	//std::cout << dst << std::endl;
	//std::cout << node->m_Local << std::endl;
}

inline void FMCSM2D_Laplace_Kernel::EvalDirect(
	const FMMPos& from_pos,
	const FMMPos& to_pos,
	const Vector<Complexd>& cf,
	Vector<Complexd> &dst)
{
	double x,y,r;
	x = to_pos.val[X_KEY] - from_pos.val[X_KEY];
	y = to_pos.val[Y_KEY] - from_pos.val[Y_KEY];
	Complexd z(x, y);

	//cf * log( sqrt(x*x + y*y) )
	dst[to_pos.index] += cf[from_pos.index] * 0.5 * std::log(sqrt(x*x + y*y));
	//dst[to_pos.index] += cf[from_pos.index] * nmc::log(z);
	//dst.dump(); printf("\n");
}

};

#endif //_FMCSM2D_LAPLACE_KERNEL_H_
