#ifndef _CSM_LAPLACE2D_KERNEL_H_
#define _CSM_LAPLACE2D_KERNEL_H_

#include "nmc/Vector2D.h"

class CSM_Laplace2D_KernelOp {
public:
	nmc::Double operator()(const nmc::Vector2Dd& spt, const nmc::Vector2Dd& cpt)
	{
		double r = sqrt((spt[0]-cpt[0])*(spt[0]-cpt[0]) + (spt[1]-cpt[1])*(spt[1]-cpt[1]));
		nmc::Double ret( log(r) );
		return ret;
	}
};

#endif //_CSM_LAPLACE2D_KERNEL_H_
