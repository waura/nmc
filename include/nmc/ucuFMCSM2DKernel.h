#ifndef _NMC_UCU_FMCSM2D_KERNEL_H_
#define _NMC_UCU_FMCSM2D_KERNEL_H_

#include "nmc/ucuFMCSM2DTree.h"

namespace nmc
{
	class ucuFMCSM2DKernel
	{
	public:
		///////////////////////////////////////////////
		/*!
			initialize kernel
		*/
		virtual void Init(ucuFMCSM2DTree* tree) =0;

		///////////////////////////////////////////////
		/*!
			calc node multiopole cefficients
			@param[in] tree
			@param[in,out] node
			@param[in] c_vec
		*/
		virtual void ElementMPCefficients(
			const ucuFMCSM2DTree* tree,
			ucuFMCSM2DTreeNode* node,
			const float2* c_vec) =0;

		///////////////////////////////////////////////
		/*!
			translation multipole coefficient from children node
			@param[in] term_num
			@param[in,out] node parent node
		*/
		virtual void M2M_Translation(
			int term_num,
			ucuFMCSM2DTreeNode* node) =0;

		////////////////////////////////////////////////
		/*!
			translation local coefficient to children node
			@param[in] term_num
			@param[in,out] node 
		*/
		virtual void L2L_Translation(
			int term_num,
			ucuFMCSM2DTreeNode* node) =0;

		/////////////////////////////////////////////////
		/*!
			translation multipole coefficient to local coefficient
			@param[in] term_num
			@param[in] node
		*/
		virtual void M2L_Translation(
			int term_num,
			ucuFMCSM2DTreeNode* node) =0;

		//////////////////////////////////////////////////
		/*!
			@param[in] tree
			@param[in] node
			@param[out] dst matrix vector product
		*/
		virtual void EvalMultipole(
			const ucuFMCSM2DTree* tree,
			const ucuFMCSM2DTreeNode* node,
			float2* dst) =0;

		//////////////////////////////////////////////////
		/*!
			@param[in] tree
			@param[in] node
			@param[in] cf coefficient
			@param[out] dst matrix vector product
		*/
		virtual void EvalDirect(
			const ucuFMCSM2DTree* tree,
			const ucuFMCSM2DTreeNode* node,
			const float2* cf,
			float2* dst) =0;
	};
};

#endif //_NMC_UCU_FMCSM2D_KERNEL_H_