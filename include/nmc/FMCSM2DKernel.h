#ifndef _FMCSM2D_KERNEL_H_
#define _FMCSM2D_KERNEL_H_

#include "nmc/FMCSM2DTree.h"

namespace nmc
{
	template<class VAL_TYPE>
	class FMCSM2DKernel
	{
	public:
		///////////////////////////////////////////////
		/*!
			initialize kernel
		*/
		virtual void Init(FMCSM2DTree<VAL_TYPE>* tree) =0;

		///////////////////////////////////////////////
		/*!
			calc node multiopole cefficients
			param[in,out] node
			param[in] x cefficients
		*/
		virtual void ElementMPCefficients(FMCSM2DTreeNode<VAL_TYPE>* node, 
			const Vector<VAL_TYPE>& x) =0;

		///////////////////////////////////////////////
		/*!
			translation multipole coefficient from children node
			@param[in,out] node
		*/
		virtual void M2M_Translation(FMCSM2DTreeNode<VAL_TYPE>* node) =0;

		////////////////////////////////////////////////
		/*!
			translation local coefficient to children node
			@param[in,out] node 
		*/
		virtual void L2L_Translation(FMCSM2DTreeNode<VAL_TYPE>* node) =0;

		/////////////////////////////////////////////////
		/*!
			@param[in] src
			@param[out dst
		*/
		virtual void M2L_Translation(const FMCSM2DTreeNode<VAL_TYPE>* src, 
			FMCSM2DTreeNode<VAL_TYPE>* dst) =0;

		//////////////////////////////////////////////////
		/*!
			@param[in] node
			@param[out] dst matrix vector product
		*/
		virtual void EvalMultipole(const FMCSM2DTreeNode<VAL_TYPE>* node,
			Vector<VAL_TYPE>& dst) =0;

		//////////////////////////////////////////////////
		/*!
			@param[in] from_pos
			@param[in] to_pos
			@param[in] cf coefficient
			@param[out] dst matrix vector product
		*/
		virtual void EvalDirect(const FMMPos& from_pos, const FMMPos& to_pos,
			const Vector<VAL_TYPE>& cf, Vector<VAL_TYPE>& dst) =0;
	};
};

#endif //_FMCSM_KERNEL_H_