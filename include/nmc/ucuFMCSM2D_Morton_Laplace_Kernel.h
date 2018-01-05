#ifndef _NMC_UCU_FMCSM2D_MORTON_LAPLACE_KERNEL_H_
#define _NMC_UCU_FMCSM2D_MORTON_LAPLACE_KERNEL_H_

#include "nmc/ucuFMCSM2DKernel.h"
#include "nmc/ucuComplex.h"
#include "nmc/ucuFunc.h"

#include <cutil_inline.h>

namespace nmc
{
	class ucuFMCSM2D_Morton_Laplace_Kernel :
		public ucuFMCSM2DKernel
	{
	public:
		///////////////////////////////////////////////
		/*!
			initialize kernel
		*/
		void Init(ucuFMCSM2DTree* tree);

		///////////////////////////////////////////////
		/*!
			calc node multiopole cefficients
			@param[in] tree
			@param[in,out] node
			@param[in] c_vec
		*/
		void ElementMPCefficients(
			const ucuFMCSM2DTree* tree,
			ucuFMCSM2DTreeNode* node,
			const float2* c_vec);

		///////////////////////////////////////////////
		/*!
			translation multipole coefficient from children node
			@param[in] term_num
			@param[in,out] node parent node
		*/
		void M2M_Translation(
			int term_num,
			ucuFMCSM2DTreeNode* node);

		////////////////////////////////////////////////
		/*!
			translation local coefficient to children node
			@param[in] term_num
			@param[in,out] node 
		*/
		void L2L_Translation(
			int term_num,
			ucuFMCSM2DTreeNode* node);

		/////////////////////////////////////////////////
		/*!
			translation multipole coefficient to local coefficient
			@param[in] term_num
			@param[in] node
		*/
		void M2L_Translation(
			int term_num,
			ucuFMCSM2DTreeNode* node);

		//////////////////////////////////////////////////
		/*!
			@param[in] tree
			@param[in] node
			@param[out] dst matrix vector product
		*/
		void EvalMultipole(
			const ucuFMCSM2DTree* tree,
			const ucuFMCSM2DTreeNode* node,
			float2* dst);

		//////////////////////////////////////////////////
		/*!
			@param[in] tree
			@param[in] node
			@param[in] cf coefficient
			@param[out] dst matrix vector product
		*/
		void EvalDirect(
			const ucuFMCSM2DTree* tree,
			const ucuFMCSM2DTreeNode* node,
			const float2* cf,
			float2* dst);

	private:
		int m_block_size;
	};

	void ucuFMCSM2D_Morton_Laplace_Kernel::Init(
		nmc::ucuFMCSM2DTree *tree)
	{
		m_block_size = 64;
	}

	__global__ void M_Laplace2D_P2M_Kernel(
		int term_num,
		float cell_size,
		float min_center_x,
		float max_center_y,
		const int* startIndexes,
		const int* endIndexes,
		const float* pointsX,
		const float* pointsY,
		const int* pointIndexes,
		const float2* c_vec,
		float2* MultipoleCoeff)
	{
		int s,t;
		extern __shared__ float M[];

		unsigned int num_of_cell = gridDim.x * blockDim.x;
		unsigned int block_size = blockDim.x;
		unsigned int elem_size = blockDim.x * (term_num + 1);
		unsigned int thread_idx = threadIdx.x;

		unsigned int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
		unsigned short i, j;
		GetMortonCellPos(cell_idx, &i, &j);
		
		float center_x = i * cell_size + min_center_x;
		float center_y = max_center_y - j * cell_size;

		int start_index = startIndexes[ cell_idx ];
		int end_index = endIndexes[ cell_idx ];

		cuComplex z;
		cuComplex c;
		cuComplex tmp;

		//init zero
		for(s=0; s<=term_num; s++){
			M[0 * elem_size + s * block_size + thread_idx] = 0.0;
			M[1 * elem_size + s * block_size + thread_idx] = 0.0;
		}

		for(t=start_index; t<end_index; t++){
			int source_index = pointIndexes[t];
			if(source_index >= 0)
				continue;

			z = make_cuFloatComplex(
				pointsX[t] - center_x,
				pointsY[t] - center_y);

			c = c_vec[ -source_index-1 ];

			//M_{0} += c
			tmp = make_cuFloatComplex(
				M[0 * elem_size + 0 * block_size + thread_idx],
				M[1 * elem_size + 0 * block_size + thread_idx]);
			tmp = cuCaddf(tmp, c);
			M[0 * elem_size + 0 * block_size + thread_idx] = cuCrealf(tmp);
			M[1 * elem_size + 0 * block_size + thread_idx] = cuCimagf(tmp);

			for(s=1; s<=term_num; s++){
				//M_{s} += -c * (z^s) / s
				tmp = make_cuFloatComplex(
					M[0 * elem_size + s * block_size + thread_idx],
					M[1 * elem_size + s * block_size + thread_idx]);
				tmp = cuCaddf(
						tmp,
						ucuCmulf(
							-1.0/s,
							cuCmulf(c, ucuCpowf(z, s))));
				M[0 * elem_size + s * block_size + thread_idx] = cuCrealf(tmp);
				M[1 * elem_size + s * block_size + thread_idx] = cuCimagf(tmp);
			}
		}

		//copy to device memory
		for(s=0; s<=term_num; s++){
			tmp = make_cuFloatComplex(
				M[0 * elem_size + s * block_size + thread_idx],
				M[1 * elem_size + s * block_size + thread_idx]);
			MultipoleCoeff[s * num_of_cell + cell_idx] = tmp;
		}
	}

	void ucuFMCSM2D_Morton_Laplace_Kernel::ElementMPCefficients(
		const ucuFMCSM2DTree* tree,
		ucuFMCSM2DTreeNode *node,
		const float2* c_vec)
	{
		assert(tree->m_pStartIndexes);
		assert(tree->m_pEndIndexes);

		//copy to device memory
		cutilSafeCall(
			cudaMemcpy(
				tree->m_pCfVec_dev,
				c_vec,
				sizeof(float2) * tree->m_SourceNum/2,
				cudaMemcpyHostToDevice) );

		int node_w, node_h;
		node_w = node_h = std::pow((float)2, node->m_NodeLevel);
		unsigned int num_of_cell = node_w * node_h;
		float cell_size = node->m_CellSize;
		float min_center_x = -cell_size * (node_w/2) + (cell_size/2.0);
		float max_center_y = cell_size * (node_h/2) - (cell_size/2.0);

		int block_size = m_block_size; //power of 2
		int grid_size = num_of_cell / block_size;
		if(num_of_cell < grid_size){
			block_size = num_of_cell;
			grid_size = num_of_cell / block_size;
		}
		dim3 Dg(grid_size, 1, 1);
		dim3 Db(block_size, 1, 1);
		M_Laplace2D_P2M_Kernel<<<Dg, Db, 2*(tree->m_TermNum+1)*block_size*sizeof(float)>>>(
			tree->m_TermNum,
			node->m_CellSize,
			min_center_x,
			max_center_y,
			tree->m_pStartIndexes_dev,
			tree->m_pEndIndexes_dev,
			tree->m_pPointsX_dev,
			tree->m_pPointsY_dev,
			tree->m_pPointIndexes_dev,
			tree->m_pCfVec_dev,
			node->m_MultipoleCoeff_dev);
		cutilCheckMsg("run P2M kernel");
		cudaThreadSynchronize();

		////debug print
		//size_t num_of_coeff = num_of_cell * (tree->m_TermNum + 1);
		//cutilSafeCall(
		//	cudaMemcpy(
		//		node->m_MultipoleCoeff,
		//		node->m_MultipoleCoeff_dev,
		//		num_of_coeff * sizeof(float2),
		//		cudaMemcpyDeviceToHost) );
		//std::cout << std::endl;
		//for(int s=0; s<= tree->m_TermNum; s++){
		//	std::cout << node->m_MultipoleCoeff[s * num_of_cell + 245] << std::endl;
		//}
		//std::cout << std::endl;
		//std::cout << std::endl;
		//for(int s=0; s<= tree->m_TermNum; s++){
		//	std::cout << node->m_MultipoleCoeff[s * num_of_cell + 249] << std::endl;
		//}
		//std::cout << std::endl;
	}

	__global__ void M_Laplace2D_M2M_Kernel1(
		int term_num,
		float child_cell_size,
		unsigned int num_of_cell,
		unsigned int num_of_ccell,
		const float2* sMultipoleCoeff,
		float2* dMultipoleCoeff)
	{
		int s,t;
		extern __shared__ float M[];

		unsigned int block_size = blockDim.y * blockDim.x;
		unsigned int elem_size = block_size * (term_num + 1);
		unsigned int thread_idx = blockDim.y * threadIdx.x + threadIdx.y;
		unsigned int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
		unsigned int child_cell_idx = GetMortonChildCellIndex(cell_idx, threadIdx.y);
		
		float half_ccell_size = child_cell_size/2;
		float cxs[] = {-half_ccell_size, half_ccell_size, -half_ccell_size,  half_ccell_size};
		float cys[] = { half_ccell_size, half_ccell_size, -half_ccell_size, -half_ccell_size};

		cuComplex tmp;
		cuComplex M_t;
		cuComplex z = make_cuFloatComplex(cxs[threadIdx.y], cys[threadIdx.y]);

		//copy child Multipole coeff to shared memory
		for(t=0; t<=term_num; t++){
			M_t = sMultipoleCoeff[t * num_of_ccell + child_cell_idx];
			M[0 * elem_size + block_size * t + thread_idx] = cuCrealf(M_t);
			M[1 * elem_size + block_size * t + thread_idx] = cuCimagf(M_t);
		}

		//M_{0] += M_0
		//tmp = sMultipoleCoeff[0 * num_of_ccell + child_cell_idx];
		//M[0 * elem_size + blockDim.x *(0 + threadIdx.x) + threadIdx.y] = tmp.x;
		//M[1 * elem_size + blockDim.x *(0 + threadIdx.x) + threadIdx.y] = tmp.y;

		for(s=term_num; s>=1; s--){
			tmp = make_cuFloatComplex(0.0, 0.0);

			//M_{s} += -(M_0 * z^s) / s
			M_t = make_cuFloatComplex(
				M[0 * elem_size + block_size * 0 + thread_idx],
				M[1 * elem_size + block_size * 0 + thread_idx]);
			tmp = cuCaddf(
				tmp,
				ucuCmulf(-1.0/s, cuCmulf(M_t, ucuCpowf(z,s))));

			for(t=1; t<=s; t++){
				//M_t * z^(s-t) * C(s-1, t-1)
				M_t = make_cuFloatComplex(
					M[0 * elem_size + block_size * t + thread_idx],
					M[1 * elem_size + block_size * t + thread_idx]);
				tmp = cuCaddf(
					tmp,
					ucuCmulf(binomial_coefficientf(s-1, t-1),
						cuCmulf(M_t, ucuCpowf(z, s-t))));
			}

			//
			M[0 * elem_size + block_size * s + thread_idx] = cuCrealf(tmp);
			M[1 * elem_size + block_size * s + thread_idx] = cuCimagf(tmp);
		}
		
		//
		__syncthreads();

		//write multipole coeff to device memory
		if(threadIdx.y == 0){
			for(t=0; t<=term_num; t++){
				tmp = make_cuFloatComplex(0.0, 0.0);
				for(s=0; s<4; s++){
					M_t = make_cuFloatComplex(
						M[0 * elem_size + block_size * t + thread_idx + s],
						M[1 * elem_size + block_size * t + thread_idx + s]);
					tmp = cuCaddf(tmp, M_t);
				}
				dMultipoleCoeff[t * num_of_cell + cell_idx] = tmp;
			}
		}
	}

	__global__ void M_Laplace2D_M2M_Kernel2(
		int term_num,
		float child_cell_size,
		unsigned int num_of_cell,
		unsigned int num_of_ccell,
		const float2* sMultipoleCoeff,
		float2* dMultipoleCoeff)
	{
		int s,t;
		extern __shared__ float M[];

		unsigned int block_size = blockDim.y * blockDim.x;
		unsigned int elem_size = block_size * (term_num + 1);
		unsigned int thread_idx = blockDim.x * threadIdx.y + threadIdx.x;
		unsigned int cell_idx = blockIdx.x * blockDim.y + threadIdx.y;
		unsigned int child_cell_idx = GetMortonChildCellIndex(cell_idx, threadIdx.x);
		
		float half_ccell_size = child_cell_size/2;
		float cxs[] = {-half_ccell_size, half_ccell_size, -half_ccell_size,  half_ccell_size};
		float cys[] = { half_ccell_size, half_ccell_size, -half_ccell_size, -half_ccell_size};

		cuComplex tmp;
		cuComplex M_t;
		cuComplex z = make_cuFloatComplex(cxs[threadIdx.x], cys[threadIdx.x]);

		//copy child Multipole coeff to shared memory
		for(t=0; t<=term_num; t++){
			M_t = sMultipoleCoeff[t * num_of_ccell + child_cell_idx];
			M[0 * elem_size + block_size * t + thread_idx] = cuCrealf(M_t);
			M[1 * elem_size + block_size * t + thread_idx] = cuCimagf(M_t);
		}

		//M_{0] += M_0
		//tmp = sMultipoleCoeff[0 * num_of_ccell + child_cell_idx];
		//M[0 * elem_size + blockDim.x *(0 + threadIdx.x) + threadIdx.y] = tmp.x;
		//M[1 * elem_size + blockDim.x *(0 + threadIdx.x) + threadIdx.y] = tmp.y;

		for(s=term_num; s>=1; s--){
			tmp = make_cuFloatComplex(0.0, 0.0);

			//M_{s} += -(M_0 * z^s) / s
			M_t = make_cuFloatComplex(
				M[0 * elem_size + block_size * 0 + thread_idx],
				M[1 * elem_size + block_size * 0 + thread_idx]);
			tmp = cuCaddf(
				tmp,
				ucuCmulf(-1.0/s, cuCmulf(M_t, ucuCpowf(z,s))));

			for(t=1; t<=s; t++){
				//M_t * z^(s-t) * C(s-1, t-1)
				M_t = make_cuFloatComplex(
					M[0 * elem_size + block_size * t + thread_idx],
					M[1 * elem_size + block_size * t + thread_idx]);
				tmp = cuCaddf(
					tmp,
					ucuCmulf(binomial_coefficientf(s-1, t-1),
						cuCmulf(M_t, ucuCpowf(z, s-t))));
			}

			//
			M[0 * elem_size + block_size * s + thread_idx] = cuCrealf(tmp);
			M[1 * elem_size + block_size * s + thread_idx] = cuCimagf(tmp);
		}
		
		//
		__syncthreads();

		//write multipole coeff to device memory
		if(threadIdx.x == 0){
			for(t=0; t<=term_num; t++){
				tmp = make_cuFloatComplex(0.0, 0.0);
				for(s=0; s<4; s++){
					M_t = make_cuFloatComplex(
						M[0 * elem_size + block_size * t + thread_idx + s],
						M[1 * elem_size + block_size * t + thread_idx + s]);
					tmp = cuCaddf(tmp, M_t);
				}
				dMultipoleCoeff[t * num_of_cell + cell_idx] = tmp;
			}
		}
	}

	void ucuFMCSM2D_Morton_Laplace_Kernel::M2M_Translation(
		int term_num,
		nmc::ucuFMCSM2DTreeNode *node)
	{
		int node_w, node_h;
		node_w = node_h = std::pow((float)2, node->m_NodeLevel);
		int num_of_cell = node_w * node_h;
		int num_of_ccell = 4 * num_of_cell;
		ucuFMCSM2DTreeNode* child_node = node->m_pNextTreeNode;

		int block_size = m_block_size/4; //power of 2
		int grid_size = num_of_cell / block_size;
		if(num_of_cell < block_size){
			block_size = num_of_cell;
			grid_size = num_of_cell / block_size;
		}
		dim3 Dg(grid_size, 1, 1);

		//
		dim3 Db(block_size, 4, 1);
		M_Laplace2D_M2M_Kernel1<<<Dg, Db, 4*2*(term_num+1)*block_size*sizeof(float)>>>(
			term_num,
			child_node->m_CellSize,
			num_of_cell,
			num_of_ccell,
			child_node->m_MultipoleCoeff_dev,
			node->m_MultipoleCoeff_dev);

		//
		//dim3 Db(4, block_size, 1);
		//M_Laplace2D_M2M_Kernel2<<<Dg, Db, 4*2*(term_num+1)*block_size*sizeof(float)>>>(
		//	term_num,
		//	child_node->m_CellSize,
		//	num_of_cell,
		//	num_of_ccell,
		//	child_node->m_MultipoleCoeff_dev,
		//	node->m_MultipoleCoeff_dev);
		
		cutilCheckMsg("run M2M kernel");
		cudaThreadSynchronize();

		////debug print
		//if(node->m_NodeLevel == 1){
		//	size_t num_of_coeff = num_of_cell * (term_num + 1);
		//	cutilSafeCall( cudaMemcpy(
		//		node->m_MultipoleCoeff,
		//		node->m_MultipoleCoeff_dev,
		//		num_of_coeff * sizeof(float2),
		//		cudaMemcpyDeviceToHost) );
		//	std::cout << std::endl;
		//	for(int s=0; s<=term_num; s++){
		//		std::cout << node->m_MultipoleCoeff[s * num_of_cell + 0] << std::endl;
		//	}
		//	std::cout << std::endl;
		//	std::cout << std::endl;
		//	for(int s=0; s<=term_num; s++){
		//		std::cout << node->m_MultipoleCoeff[s * num_of_cell + 2] << std::endl;
		//	}
		//	std::cout << std::endl;
		//}
	}

	__global__ void M_Laplace2D_L2L_Kernel1(
		int term_num,
		float child_cell_size,
		unsigned int num_of_cell,
		unsigned int num_of_ccell,
		const float2* sLocalCoeff,
		float2* dLocalCoeff)
	{
		int s,t;
		extern __shared__ float L[];

		unsigned int block_size = blockDim.y * blockDim.x;
		unsigned int elem_size = block_size * (term_num + 1);
		unsigned int thread_idx = blockDim.y * threadIdx.x + threadIdx.y;
		unsigned int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
		unsigned int child_cell_idx = GetMortonChildCellIndex(cell_idx, threadIdx.y);

		float half_ccell_size = child_cell_size/2;
		float cxs[] = {-half_ccell_size, half_ccell_size, -half_ccell_size,  half_ccell_size};
		float cys[] = { half_ccell_size, half_ccell_size, -half_ccell_size, -half_ccell_size};

		cuComplex L_s;
		cuComplex L_t;
		cuComplex z = make_cuFloatComplex(cxs[threadIdx.y], cys[threadIdx.y]);

		//copy parent Local coeff to shared memory
		for(t=0; t<=term_num; t++){
			L_t = sLocalCoeff[t * num_of_cell + cell_idx];
			L[0 * elem_size + block_size * t + thread_idx] = cuCrealf(L_t);
			L[1 * elem_size + block_size * t + thread_idx] = cuCimagf(L_t);
		}

		//
		__syncthreads();

		//
		for(s=0; s<=term_num; s++){
			L_s = dLocalCoeff[s * num_of_ccell + child_cell_idx];
			for(t=s; t<=term_num; t++){
				//L_{s} += L_{t} * C(t, s) * z^(t - s)
				L_t = make_cuFloatComplex(
					L[0 * elem_size + block_size * t + thread_idx],
					L[1 * elem_size + block_size * t + thread_idx]);

				L_s = cuCaddf(
					L_s,
					ucuCmulf(
						binomial_coefficientf(t, s),
						cuCmulf(L_t, ucuCpowf(z, t-s))));
			}

			//write local coeff to device memory
			dLocalCoeff[s * num_of_ccell + child_cell_idx] = L_s;
		}
	}

	__global__ void M_Laplace2D_L2L_Kernel2(
		int term_num,
		float child_cell_size,
		unsigned int num_of_cell,
		unsigned int num_of_ccell,
		const float2* sLocalCoeff,
		float2* dLocalCoeff)
	{
		int s,t;
		extern __shared__ float L[];

		unsigned int block_size = blockDim.y * blockDim.x;
		unsigned int elem_size = block_size * (term_num + 1);
		unsigned int thread_idx = blockDim.x * threadIdx.y + threadIdx.x;
		unsigned int cell_idx = blockIdx.x * blockDim.y + threadIdx.y;
		unsigned int child_cell_idx = GetMortonChildCellIndex(cell_idx, threadIdx.x);

		float half_ccell_size = child_cell_size/2;
		float cxs[] = {-half_ccell_size, half_ccell_size, -half_ccell_size,  half_ccell_size};
		float cys[] = { half_ccell_size, half_ccell_size, -half_ccell_size, -half_ccell_size};

		cuComplex L_s;
		cuComplex L_t;
		cuComplex z = make_cuFloatComplex(cxs[threadIdx.x], cys[threadIdx.x]);

		//copy parent Local coeff to shared memory
		for(t=0; t<=term_num; t++){
			L_t = sLocalCoeff[t * num_of_cell + cell_idx];
			L[0 * elem_size + block_size * t + thread_idx] = cuCrealf(L_t);
			L[1 * elem_size + block_size * t + thread_idx] = cuCimagf(L_t);
		}

		//
		__syncthreads();

		//
		for(s=0; s<=term_num; s++){
			L_s = dLocalCoeff[s * num_of_ccell + child_cell_idx];
			for(t=s; t<=term_num; t++){
				//L_{s} += L_{t} * C(t, s) * z^(t - s)
				L_t = make_cuFloatComplex(
					L[0 * elem_size + block_size * t + thread_idx],
					L[1 * elem_size + block_size * t + thread_idx]);

				L_s = cuCaddf(
					L_s,
					ucuCmulf(
						binomial_coefficientf(t, s),
						cuCmulf(L_t, ucuCpowf(z, t-s))));
			}

			//write local coeff to device memory
			dLocalCoeff[s * num_of_ccell + child_cell_idx] = L_s;
		}
	}

	void ucuFMCSM2D_Morton_Laplace_Kernel::L2L_Translation(
		int  term_num,
		nmc::ucuFMCSM2DTreeNode *node)
	{
		cuComplex z;
		int node_w, node_h;
		node_w = node_h = std::pow((float)2, node->m_NodeLevel);
		int num_of_cell = node_w * node_h;
		int num_of_ccell = 4 * node_w * node_h;
		ucuFMCSM2DTreeNode* child_node = node->m_pNextTreeNode;

		int block_size = m_block_size/4; //power of 2
		int grid_size = num_of_cell / block_size;
		if(num_of_cell < block_size){
			block_size = num_of_cell;
			grid_size = num_of_cell / block_size;
		}
		dim3 Dg(grid_size, 1, 1);

		//this kernel is bug
		//dim3 Db(block_size, 4, 1);
		//M_Laplace2D_L2L_Kernel1<<<Dg, Db, 4*2*(term_num+1)*block_size*sizeof(float)>>>(
		//	term_num,
		//	child_node->m_CellSize,
		//	num_of_cell,
		//	num_of_ccell,
		//	node->m_LocalCoeff_dev,
		//	child_node->m_LocalCoeff_dev);
		
		dim3 Db(4, block_size, 1);
		M_Laplace2D_L2L_Kernel2<<<Dg, Db, 4*2*(term_num+1)*block_size*sizeof(float)>>>(
			term_num,
			child_node->m_CellSize,
			num_of_cell,
			num_of_ccell,
			node->m_LocalCoeff_dev,
			child_node->m_LocalCoeff_dev);
		cutilCheckMsg("run L2L kernel");
		cudaThreadSynchronize();

		////debug print
		//if(node->m_NodeLevel == 4){
		//	size_t num_of_coeff = num_of_cell * (term_num + 1);
		//	cutilSafeCall( cudaMemcpy(
		//		node->m_LocalCoeff,
		//		node->m_LocalCoeff_dev,
		//		num_of_coeff * sizeof(float2),
		//		cudaMemcpyDeviceToHost) );
		//	std::cout << std::endl;
		//	for(int s=0; s<=term_num; s++){
		//		std::cout << node->m_LocalCoeff[s * num_of_cell + 0] << std::endl;
		//	}
		//	std::cout << std::endl;
		//	std::cout << std::endl;
		//	for(int s=0; s<=term_num; s++){
		//		std::cout << node->m_LocalCoeff[s * num_of_cell + 39] << std::endl;
		//	}
		//	std::cout << std::endl;
		//}
	}

	__global__ void M_Laplace2D_M2L_Kernel1(
		int term_num,
		float cell_size,
		unsigned int node_size,
		unsigned int num_of_cell,
		const float2* MultipoleCoeff,
		float2* LocalCoeff)
	{
		int dif_idx_x[4][27] = {
			{
				-2, -1, 0, 1, 2, 3,
				-2,           2, 3,
				-2,           2, 3,
				-2,           2, 3,
				-2, -1, 0, 1, 2, 3,
				-2, -1, 0, 1, 2, 3,
			},
			{
				-3, -2, -1, 0, 1, 2,
				-3, -2,           2,
				-3, -2,           2,
				-3, -2,           2,
				-3, -2, -1, 0, 1, 2,
				-3, -2, -1, 0, 1, 2,
			},
			{
				-2, -1, 0, 1, 2, 3,
				-2, -1, 0, 1, 2, 3,
				-2,           2, 3,
				-2,           2, 3,
				-2,           2, 3,
				-2, -1, 0, 1, 2, 3,
			},
			{
				-3, -2, -1, 0, 1, 2,
				-3, -2, -1, 0, 1, 2,
				-3, -2,           2,
				-3, -2,           2,
				-3, -2,           2,
				-3, -2, -1, 0, 1, 2,
			},
		};
		int dif_idx_y[4][27] = {
			{
				-2, -2, -2, -2, -2, -2,
				-1,             -1, -1,
				 0,              0,  0,
				 1,              1,  1,
				 2,  2,  2,  2,  2,  2,
				 3,  3,  3,  3,  3,  3,
			},
			{
				-2, -2, -2, -2, -2, -2,
				-1, -1,             -1,
				 0,  0,              0,
				 1,  1,              1,
				 2,  2,  2,  2,  2,  2,
				 3,  3,  3,  3,  3,  3,
			},
			{
				-3, -3, -3, -3, -3, -3,
				-2, -2, -2, -2, -2, -2,
				-1,             -1, -1,
				 0,              0,  0,
				 1,              1,  1,
				 2,  2,  2,  2,  2,  2,
			},
			{
				-3, -3, -3, -3, -3, -3,
				-2, -2, -2, -2, -2, -2,
				-1, -1,             -1,
				 0,  0,              0,
				 1,  1,              1,
				 2,  2,  2,  2,  2,  2,
			},
		};

		int s,t;
		int k;
		unsigned short i,j;
		unsigned int block_size = blockDim.y * blockDim.x;
		unsigned int elem_size = block_size * (term_num + 1);
		unsigned int thread_idx = blockDim.x * threadIdx.y + threadIdx.x;
		unsigned int dcell_idx = blockIdx.x * blockDim.x + threadIdx.x;

		extern __shared__ float buffer[];
		float* L = &(buffer[0]);
		float* M = &(buffer[2*elem_size]);

		cuComplex z;
		cuComplex L_s;
		cuComplex M_t;

		unsigned int block_idx = dcell_idx % 4;

		//init zero
		for(s=0; s<=term_num; s++){
			L[0 * elem_size + block_size * s + thread_idx] = 0.0;
			L[1 * elem_size + block_size * s + thread_idx] = 0.0;
		}

		GetMortonCellPos(dcell_idx, &i, &j);

		for(k=0; k<27; k++){
			int difx = dif_idx_x[block_idx][k];
			int dify = dif_idx_y[block_idx][k];
			int idx_x = i + difx;
			int idx_y = j + dify;

			if(idx_x < 0 || node_size <= idx_x)
				continue;
			if(idx_y < 0 || node_size <= idx_y)
				continue;

			unsigned int scell_idx = GetMortonCellIndex(idx_x, idx_y);

			z = make_cuFloatComplex(
				cell_size * difx,
				-cell_size * dify);

			//copy Multi coeff to shared memory
			for(t=0; t<=term_num; t++){
				M_t = MultipoleCoeff[t * num_of_cell + scell_idx];
				M[0 * elem_size + block_size * t + thread_idx] = cuCrealf(M_t);
				M[1 * elem_size + block_size * t + thread_idx] = cuCimagf(M_t);
			}

			//L_{0} += M_{0} * log( -z )
			L_s = make_cuFloatComplex(
				L[0 * elem_size + block_size * 0 + thread_idx],
				L[1 * elem_size + block_size * 0 + thread_idx]);
			M_t = make_cuFloatComplex(
				M[0 * elem_size + block_size * 0 + thread_idx],
				M[1 * elem_size + block_size * 0 + thread_idx]);
			L_s = cuCaddf(
				L_s,
				ucuCmulf(
					M_t,
					logf( cuCabsf(ucuCmulf(-1.0, z)) )));
					//ucuClogf( ucuCmulf(-1.0, z) )));
					
			for(t=1; t<=term_num; t++){
				//L_{0} += (-1)^t * M_{t} / z^t
				//float sign = ((t % 2) == 0) ? 1.0 : -1.0;
				float sign = ((t & 0x1) == 0) ? 1.0 : -1.0;
				M_t = make_cuFloatComplex(
					M[0 * elem_size + block_size * t + thread_idx],
					M[1 * elem_size + block_size * t + thread_idx]);
				
				L_s = cuCaddf(
					L_s,
					cuCdivf(
						ucuCmulf(sign, M_t),
						ucuCpowf(z, t)));
			}
			L[0 * elem_size + block_size * 0 + thread_idx] = cuCrealf(L_s);
			L[1 * elem_size + block_size * 0 + thread_idx] = cuCimagf(L_s);

			for(s=1; s<=term_num; s++){
				//L_{s} += -M_{0} / (s * z^s)
				L_s = make_cuFloatComplex(
					L[0 * elem_size + block_size * s + thread_idx],
					L[1 * elem_size + block_size * s + thread_idx]);
				M_t = make_cuFloatComplex(
					M[0 * elem_size + block_size * 0 + thread_idx],
					M[1 * elem_size + block_size * 0 + thread_idx]);

				L_s = cuCaddf(
					L_s,
					cuCdivf(
						ucuCmulf(-1.0, M_t),
						ucuCmulf(static_cast<float>(s), ucuCpowf(z, s))));

				for(t=1; t<=term_num; t++){
					//L_{s} += (-1)^t * M_{t} * C(s+t-1, t-1) / z^(s+t)
					//float sign = ((t % 2) == 0) ? 1.0 : -1.0;
					float sign = ((t & 0x1) == 0) ? 1.0 : -1.0;
					M_t = make_cuFloatComplex(
						M[0 * elem_size + block_size * t + thread_idx],
						M[1 * elem_size + block_size * t + thread_idx]);

					L_s = cuCaddf(
						L_s,
						cuCdivf(
							ucuCmulf(
								sign * binomial_coefficientf(s+t-1, t-1),
								M_t),
							ucuCpowf(z, s+t)));
				}
				L[0 * elem_size + block_size * s + thread_idx] = cuCrealf(L_s);
				L[1 * elem_size + block_size * s + thread_idx] = cuCimagf(L_s);
			}
		}

		//copy Local coeff to device memory
		for(s=0; s<=term_num; s++){
			L_s = make_cuFloatComplex(
				L[0 * elem_size + block_size * s + thread_idx],
				L[1 * elem_size + block_size * s + thread_idx]);
			LocalCoeff[s * num_of_cell + dcell_idx] = L_s;
		}
	}

	void ucuFMCSM2D_Morton_Laplace_Kernel::M2L_Translation(
		int term_num,
		nmc::ucuFMCSM2DTreeNode *node)
	{
		int i,j,k;
		int s,t;
		cuComplex z;
		int node_w, node_h;
		node_w = node_h = std::pow((float)2, node->m_NodeLevel);
		int num_of_cell = node_w * node_h;

		int block_size = m_block_size; //power of 2
		int grid_size = num_of_cell / block_size;
		if(num_of_cell < block_size){
			block_size = num_of_cell;
			grid_size = num_of_cell / block_size;
		}
		dim3 Dg(grid_size, 1, 1);
		dim3 Db(block_size, 1, 1);
		M_Laplace2D_M2L_Kernel1<<<Dg, Db, 2*2*(term_num+1)*block_size*sizeof(float)>>>(
			term_num,
			node->m_CellSize,
			node_w,
			num_of_cell,
			node->m_MultipoleCoeff_dev,
			node->m_LocalCoeff_dev);
		cutilCheckMsg("run M2L kernel");
		cudaThreadSynchronize();

		////debug print
		//if(node->m_NodeLevel == 2){
		//	size_t num_of_coeff = num_of_cell * (term_num + 1);
		//	cutilSafeCall( cudaMemcpy(
		//		node->m_LocalCoeff,
		//		node->m_LocalCoeff_dev,
		//		num_of_coeff * sizeof(float2),
		//		cudaMemcpyDeviceToHost) );
		//	std::cout << std::endl;
		//	for(int s=0; s<=term_num; s++){
		//		std::cout << node->m_LocalCoeff[s * num_of_cell + 0] << std::endl;
		//	}
		//	std::cout << std::endl;
		//	std::cout << std::endl;
		//	for(int s=0; s<=term_num; s++){
		//		std::cout << node->m_LocalCoeff[s * num_of_cell + 1] << std::endl;
		//	}
		//	std::cout << std::endl;
		//}
	}

	__global__ void M_Laplace2D_L2P_Kernel(
		int term_num,
		float cell_size,
		float min_center_x,
		float max_center_y,
		const int* startIndexes,
		const int* endIndexes,
		const float* pointsX,
		const float* pointsY,
		const int* pointIndexes,
		const float2* LocalCoeff,
		float2* dst)
	{
		int s,t;
		extern __shared__ float L[];

		unsigned int num_of_cell = gridDim.x * blockDim.x;
		unsigned int block_size = blockDim.x;
		unsigned int elem_size = blockDim.x * (term_num + 1);
		unsigned int thread_idx = threadIdx.x;

		unsigned int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
		unsigned short i,j;
		GetMortonCellPos(cell_idx, &i, &j);

		float center_x = i * cell_size + min_center_x;
		float center_y = max_center_y - j * cell_size;

		int start_index = startIndexes[ cell_idx ];
		int end_index = endIndexes[ cell_idx ];

		cuComplex z;
		cuComplex sum;
		cuComplex L_s;

		//copy Local coeff to shared memory
		for(s=0; s<=term_num; s++){
			L_s = LocalCoeff[s * num_of_cell + cell_idx];
			L[0 * elem_size + s * block_size + thread_idx] = cuCrealf(L_s);
			L[1 * elem_size + s * block_size + thread_idx] = cuCimagf(L_s);
		}

		//
		for(t=start_index; t<end_index; t++){
			int index = pointIndexes[t];
			if(index < 0)
				continue;

			cuComplex z = make_cuFloatComplex(
				pointsX[t] - center_x,
				pointsY[t] - center_y);

			cuComplex sum = make_cuFloatComplex(0.0, 0.0);
			for(s=0; s<=term_num; s++){
				//L_{s} * z^s
				L_s = make_cuFloatComplex(
					L[0 * elem_size + s * block_size + thread_idx],
					L[1 * elem_size + s * block_size + thread_idx]);
				sum = cuCaddf(
						sum,
						cuCmulf(L_s, ucuCpowf(z, s)));
			}

			//dst[index] must be zero
			//dst[index] = cuCaddf(dst[index],tmp);
			dst[index] = sum;
		}
	}

	void ucuFMCSM2D_Morton_Laplace_Kernel::EvalMultipole(
		const ucuFMCSM2DTree* tree,
		const ucuFMCSM2DTreeNode *node,
		float2* dst)
	{
		int node_w, node_h;
		node_w = node_h = std::pow((float)2, node->m_NodeLevel);
		unsigned int num_of_cell = node_w * node_h;
		float cell_size = node->m_CellSize;
		float min_center_x = -cell_size * (node_w/2) + (cell_size/2.0);
		float max_center_y = cell_size * (node_h/2) - (cell_size/2.0);

		//init zero
		cutilSafeCall(
			cudaMemset(tree->m_pDst_dev, 0, sizeof(float2) * tree->m_SourceNum/2) );

		int block_size = m_block_size; //power of 2
		int grid_size = num_of_cell / block_size;
		if(num_of_cell < block_size){
			block_size = num_of_cell;
			grid_size = num_of_cell / block_size;
		}
		dim3 Dg(grid_size, 1, 1);
		dim3 Db(block_size, 1, 1);
		M_Laplace2D_L2P_Kernel<<<Dg, Db, 2*(tree->m_TermNum+1)*block_size*sizeof(float)>>>(
			tree->m_TermNum,
			node->m_CellSize,
			min_center_x,
			max_center_y,
			tree->m_pStartIndexes_dev,
			tree->m_pEndIndexes_dev,
			tree->m_pPointsX_dev,
			tree->m_pPointsY_dev,
			tree->m_pPointIndexes_dev,
			node->m_LocalCoeff_dev,
			tree->m_pDst_dev);
		cutilCheckMsg("run L2P kernel");
		cudaThreadSynchronize();

		////debug print
		//cutilSafeCall(
		//	cudaMemcpy(
		//		dst,
		//		tree->m_pDst_dev,
		//		sizeof(float2) * tree->m_SourceNum/2,
		//		cudaMemcpyDeviceToHost) );
		//int n = tree->m_SourceNum/2;
		//std::cout << std::endl;
		//for(int s=0; s<n; s++){
		//	std::cout << dst[s] << std::endl;
		//}
		//std::cout << std::endl;
	}

	__global__ void M_Laplace2D_EvalDirect_Kernel(
		int term_num,
		unsigned int node_size,
		const int* startIndexes,
		const int* endIndexes,
		const float* pointsX,
		const float* pointsY,
		const int* pointIndexes,
		const float2* c_vec,
		float2* dst)
	{
		unsigned int num_of_cell = gridDim.x * blockDim.x;
		unsigned int block_size = blockDim.x;
		unsigned int thread_idx = threadIdx.x;

		__shared__ int dif_idx_x[9];
		dif_idx_x[0] = -1; dif_idx_x[1] = 0; dif_idx_x[2] = 1; 
		dif_idx_x[3] = -1; dif_idx_x[4] = 0; dif_idx_x[5] = 1; 
		dif_idx_x[6] = -1; dif_idx_x[7] = 0; dif_idx_x[8] = 1; 

		__shared__ int dif_idx_y[9];
		dif_idx_y[0] = -1; dif_idx_y[1] = -1; dif_idx_y[2] = -1;
		dif_idx_y[3] =  0; dif_idx_y[4] =  0; dif_idx_y[5] =  0;
		dif_idx_y[6] =  1; dif_idx_y[7] =  1; dif_idx_y[8] =  1;

		/*
		extern __shared__ int buf[];
		int* dif_idx_x = &(buf[0]);
		dif_idx_x[0 * block_size + thread_idx] = -1;
		dif_idx_x[1 * block_size + thread_idx] = 0;
		dif_idx_x[2 * block_size + thread_idx] = 1; 
		dif_idx_x[3 * block_size + thread_idx] = -1;
		dif_idx_x[4 * block_size + thread_idx] = 0;
		dif_idx_x[5 * block_size + thread_idx] = 1; 
		dif_idx_x[6 * block_size + thread_idx] = -1;
		dif_idx_x[7 * block_size + thread_idx] = 0;
		dif_idx_x[8 * block_size + thread_idx] = 1; 

		int* dif_idx_y = &(buf[block_size * 9]);
		dif_idx_y[0 * block_size + thread_idx] = -1;
		dif_idx_y[1 * block_size + thread_idx] = -1;
		dif_idx_y[2 * block_size + thread_idx] = -1;
		dif_idx_y[3 * block_size + thread_idx] =  0;
		dif_idx_y[4 * block_size + thread_idx] =  0;
		dif_idx_y[5 * block_size + thread_idx] =  0;
		dif_idx_y[6 * block_size + thread_idx] =  1;
		dif_idx_y[7 * block_size + thread_idx] =  1;
		dif_idx_y[8 * block_size + thread_idx] =  1;
		*/

		unsigned int dcell_idx = blockIdx.x * blockDim.x + threadIdx.x;
		unsigned short i, j;
		GetMortonCellPos(dcell_idx, &i, &j);

		int dstart_idx = startIndexes[ dcell_idx ];
		int dend_idx = endIndexes[ dcell_idx ];

		int s,t;
		int k;

		for(k=0; k<9; k++){
			int idx_x = i + dif_idx_x[k * block_size + thread_idx];
			int idx_y = j + dif_idx_y[k * block_size + thread_idx];

			if(idx_x < 0 || node_size <= idx_x)
				continue;
			if(idx_y < 0 || node_size <= idx_y)
				continue;

			unsigned int scell_idx = GetMortonCellIndex(idx_x, idx_y);
			int sstart_idx = startIndexes[ scell_idx ];
			int send_idx = endIndexes[ scell_idx ];

			for(s=dstart_idx; s<dend_idx; s++){
				int d_index = pointIndexes[s];
				if(d_index < 0)
					continue;

				cuComplex sum = make_cuFloatComplex(0.0, 0.0);
				for(t=sstart_idx; t<send_idx; t++){
					int s_index = pointIndexes[t];
					if(s_index >= 0)
						continue;

					cuComplex z = make_cuFloatComplex(
						pointsX[s] - pointsX[t],
						pointsY[s] - pointsY[t]);

					sum =
						cuCaddf(
							sum,
							cuCmulf(
								//ucuCmulf(0.5, cf[-s_index-1]),
								c_vec[-s_index-1],
								ucuClogf(z)));
								//logf( cuCabsf(z) )));
				}
				dst[d_index] = cuCaddf(dst[d_index], sum);
			}
		}
	}

	void ucuFMCSM2D_Morton_Laplace_Kernel::EvalDirect(
		const ucuFMCSM2DTree* tree,
		const ucuFMCSM2DTreeNode *node,
		const float2* cf,
		float2* dst)
	{
		int node_size = std::pow((float)2, node->m_NodeLevel);
		int num_of_cell = node_size * node_size;

		int block_size = m_block_size; //power of 2
		int grid_size = num_of_cell / block_size;
		if(num_of_cell < block_size){
			block_size = num_of_cell;
			grid_size = num_of_cell / block_size;
		}
		dim3 Dg(grid_size, 1, 1);
		dim3 Db(block_size, 1, 1);
		M_Laplace2D_EvalDirect_Kernel<<<Dg, Db>>>(
			tree->m_TermNum,
			node_size,
			tree->m_pStartIndexes_dev,
			tree->m_pEndIndexes_dev,
			tree->m_pPointsX_dev,
			tree->m_pPointsY_dev,
			tree->m_pPointIndexes_dev,
			tree->m_pCfVec_dev,
			tree->m_pDst_dev);
		cutilCheckMsg("run EvalDirect kernel");
		cudaThreadSynchronize();

		//copy to host memory
		cutilSafeCall(
			cudaMemcpy(
				dst,
				tree->m_pDst_dev,
				sizeof(float2) * tree->m_SourceNum/2,
				cudaMemcpyDeviceToHost) );
	}
};

#endif //_UCU_FMCSM2D_MORTON_LAPLACE_KERNEL_H_