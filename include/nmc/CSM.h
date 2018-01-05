#ifndef _NMC_CSM_H_
#define _NMC_CSM_H_

#include "Vector2D.h"
#include "Matrix.h"

namespace nmc
{

template <class VAL_TYPE>
class CSM2D {
public:
	////////////////////////////////////////
	/*!
		@param[in] source_points csm source points
		@param[in] collocation_points csm collocation points
		@param[out] csm_mat
		@param[out] csm_f
		@return
	*/
	template<class CSMBoundaryOp, class CSMKernelOp>
	static void CreateMatrix(
		CSMBoundaryOp& csm_boundaryOp,
		CSMKernelOp& csm_kernelOp,
		const VectorPoint2Dd& source_points,
		const VectorPoint2Dd& collocation_points,
		Matrix<VAL_TYPE>& csm_mat,
		Vector<VAL_TYPE>& csm_f)
	{
		assert(source_points.dim() == collocation_points.dim());

		int i,j;
		int n = source_points.dim();
		csm_mat.Init(n, n);
		csm_f.Init(n);

		Vector2Dd spt,cpt;

		//boundary condition
		for(i=0; i<n; i++){
			cpt = collocation_points[i];
			csm_f[i] = csm_boundaryOp(cpt);
		}

		//
		for(j=0; j<n; j++){
			cpt = collocation_points[j];
			for(i=0; i<n; i++){
				spt = source_points[i];
				csm_mat.SetVal(j, i, csm_kernelOp(cpt, spt));
			}
		}
	}

	//////////////////////////////////////////////
	/*!
	*/
	template<class CSMKernelOp>
	static VAL_TYPE u_N(
		CSMKernelOp& csm_kernelOp,
		const Vector<VAL_TYPE>& Q,
		const VectorPoint2Dd& source_points,
		const Vector2Dd& point)
	{
		assert(Q.dim() == source_points.dim());
		int i;
		VAL_TYPE ret=0.0;
		Vector2Dd spt;
		for(i=0; i<Q.dim(); i++){
		  ret += Q[i] * csm_kernelOp(point, source_points[i]);
		}
		return ret;
	}	
};

typedef CSM2D< Double > CSM2Dd;
typedef CSM2D< Float > CSM2Df;
typedef CSM2D< Int > CSM2Di;
typedef CSM2D< Complexd > CSM2Dcmd;

}

#endif //_NMC_CSM_H_
