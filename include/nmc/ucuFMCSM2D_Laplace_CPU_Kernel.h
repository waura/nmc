#ifndef _NMC_UCU_FMCSM2D_LAPLACE_CPU_KERNEL_H_
#define _NMC_UCU_FMCSM2D_LAPLACE_CPU_KERNEL_H_

#include "nmc/ucuFMCSM2DKernel.h"
#include "nmc/ucuComplex.h"
#include "nmc/ucuFunc.h"

namespace nmc
{
	class ucuFMCSM2D_Laplace_CPU_Kernel :
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
	};

	void ucuFMCSM2D_Laplace_CPU_Kernel::Init(
		nmc::ucuFMCSM2DTree *tree)
	{
	}

	void ucuFMCSM2D_Laplace_CPU_Kernel::ElementMPCefficients(
		const ucuFMCSM2DTree* tree,
		ucuFMCSM2DTreeNode *node,
		const float2* c_vec)
	{
		assert(tree->m_pStartIndexes);
		assert(tree->m_pEndIndexes);

		int s,t;
		int i,j,k;
		int term_num = tree->m_TermNum;
		int node_w, node_h;
		node_w = node_h = std::pow((float)2, node->m_NodeLevel);
		int num_of_cell = node_w * node_h;

		float cell_w = node->m_CellSize;
		float cell_h = node->m_CellSize;
		float min_center_x = -cell_w * (node_w/2) + (cell_w/2.0);
		float min_center_y = -cell_h * (node_h/2) + (cell_h/2.0);

		float x,y;
		float center_x;
		float center_y;

		cuComplex z;
		cuComplex c;
		cuComplex tmp;

		//init zero
		memset(node->m_MultipoleCoeff, 0, sizeof(float2) * num_of_cell * (term_num+1));

		//
		unsigned int cell_idx;
		for(j=0; j<node_h; j++){
			for(i=0; i<node_w; i++){
				cell_idx = GetCellIndex(i, j, node->m_NodeLevel);
				center_x = i * cell_w + min_center_x;
				center_y = (node_h - j - 1) * cell_h + min_center_y;

				int start_index = tree->m_pStartIndexes[ cell_idx ];
				int end_index = tree->m_pEndIndexes[ cell_idx ];

				for(k=start_index; k<end_index; k++){
					int source_index = tree->m_pPointIndexes[k];
					if(source_index >= 0)
						continue;

					z = make_cuFloatComplex(
						x = tree->m_pPointsX[k] - center_x,
						y = tree->m_pPointsY[k] - center_y);

					c = c_vec[ -source_index-1 ];

					//M_{0} += c
					tmp = node->m_MultipoleCoeff[0 * num_of_cell + cell_idx];
					tmp = cuCaddf(tmp, c);
					node->m_MultipoleCoeff[0 * num_of_cell + cell_idx] = tmp;

					for(s=1; s<=term_num; s++){
						//M_{s} += -c * (z^s) / s
						tmp = node->m_MultipoleCoeff[s * num_of_cell + cell_idx];
						tmp = cuCaddf(
								tmp,
								ucuCmulf(
									-1.0/s,
									cuCmulf(c, ucuCpowf(z, s))));
						node->m_MultipoleCoeff[s * num_of_cell + cell_idx] = tmp;
					}

					////debug print
					//std::cout << std::endl;
					//for(s=0; s<= term_num; s++){
					//	std::cout << node->m_MultipoleCoeff[s * num_of_cell + cell_idx] << std::endl;
					//}
					//std::cout << std::endl;					
				}
			}
		}
		
		////debug print
		//std::cout << std::endl;
		//for(int s=0; s<= tree->m_TermNum; s++){
		//	std::cout << node->m_MultipoleCoeff[s * num_of_cell + 399] << std::endl;
		//}
		//std::cout << std::endl;
		//std::cout << std::endl;
		//for(int s=0; s<= tree->m_TermNum; s++){
		//	std::cout << node->m_MultipoleCoeff[s * num_of_cell + 461] << std::endl;
		//}
	}

	void ucuFMCSM2D_Laplace_CPU_Kernel::M2M_Translation(
		int term_num,
		nmc::ucuFMCSM2DTreeNode *node)
	{
		int i,j,k;
		int s,t;
		cuComplex z;
		int node_w, node_h;
		node_w = node_h = std::pow((float)2, node->m_NodeLevel);
		int num_of_cell = node_w * node_h;
		int num_of_ccell = 4 * num_of_cell;
		ucuFMCSM2DTreeNode* child_node = node->m_pNextTreeNode;
		
		float x = child_node->m_CellSize/2;
		float y = child_node->m_CellSize/2;
		float cxs[] = {-x, x, -x,  x};
		float cys[] = { y, y, -y, -y};

		//init zero
		memset(node->m_MultipoleCoeff, 0, sizeof(float2) * num_of_cell * (term_num+1));
		
		//
		unsigned int cell_idx;
		unsigned int child_cell_idx;
		for(j=0; j<node_h; j++){
			for(i=0; i<node_w; i++){
				cell_idx = GetCellIndex(i, j, node->m_NodeLevel);
				assert(0 <= cell_idx && cell_idx < num_of_cell);

				for(k=0; k<4; k++){
					child_cell_idx = GetChildCellIndex(i, j, node->m_NodeLevel, k);
					assert(0 <= child_cell_idx && child_cell_idx < num_of_ccell);

					z = make_cuFloatComplex(cxs[k], cys[k]);

					//M_{0] += M_0
					node->m_MultipoleCoeff[0 * num_of_cell + cell_idx] =
						cuCaddf(
							node->m_MultipoleCoeff[0 * num_of_cell + cell_idx],
							child_node->m_MultipoleCoeff[0 * num_of_ccell + child_cell_idx]);
					for(s=1; s<=term_num; s++){
						//M_{s} += -(M_0 * z^s) / s
						node->m_MultipoleCoeff[s * num_of_cell + cell_idx] =
							cuCaddf(
								node->m_MultipoleCoeff[s * num_of_cell + cell_idx],
								ucuCmulf(
									-1.0/s,
									cuCmulf(
										child_node->m_MultipoleCoeff[0 * num_of_ccell + child_cell_idx],
										ucuCpowf(z, s))));
						for(t=1; t<=s; t++){
							//M_{s} += M_t * z^(s_t) * C(s-1, t-1)
							node->m_MultipoleCoeff[s * num_of_cell + cell_idx] =
								cuCaddf(
									node->m_MultipoleCoeff[s * num_of_cell + cell_idx],
									ucuCmulf(
										binomial_coefficientf(s-1, t-1),
										cuCmulf(
											child_node->m_MultipoleCoeff[t * num_of_ccell + child_cell_idx],
											ucuCpowf(z, s-t))));
						}
					}
				}

				////debug print
				//std::cout << std::endl;
				//for(s=0; s<= term_num; s++){
				//	std::cout << node->m_MultipoleCoeff[s * num_of_cell + cell_idx] << std::endl;
				//}
				//std::cout << std::endl;
			}
		}

		////debug print
		//if(node->m_NodeLevel == 1){
		//	std::cout << std::endl;
		//	for(int s=0; s<= term_num; s++){
		//		std::cout << node->m_MultipoleCoeff[s * num_of_cell + 0] << std::endl;
		//	}
		//	std::cout << std::endl;
		//	std::cout << std::endl;
		//	for(int s=0; s<= term_num; s++){
		//		std::cout << node->m_MultipoleCoeff[s * num_of_cell + 2] << std::endl;
		//	}
		//	std::cout << std::endl;
		//}
	}

	void ucuFMCSM2D_Laplace_CPU_Kernel::L2L_Translation(
		int  term_num,
		nmc::ucuFMCSM2DTreeNode *node)
	{
		int i,j,k;
		int s,t;
		cuComplex z;
		int node_w, node_h;
		node_w = node_h = std::pow((float)2, node->m_NodeLevel);
		int num_of_cell = node_w * node_h;
		int num_of_ccell = 4 * node_w * node_h;
		ucuFMCSM2DTreeNode* child_node = node->m_pNextTreeNode;

		float x = child_node->m_CellSize/2;
		float y = child_node->m_CellSize/2;
		float cxs[] = {-x, x, -x,  x};
		float cys[] = { y, y, -y, -y};

		unsigned int cell_idx;
		unsigned int child_cell_idx;
		for(j=0; j<node_h; j++){
			for(i=0; i<node_w; i++){
				cell_idx = GetCellIndex(i, j, node->m_NodeLevel);

				////debug print
				//std::cout << std::endl;
				//for(s=0; s<= term_num; s++){
				//	std::cout << node->m_LocalCoeff[s * num_of_cell + cell_idx] << std::endl;
				//}
				//std::cout << std::endl;

				for(k=0; k<4; k++){
					child_cell_idx = GetChildCellIndex(i, j, node->m_NodeLevel, k);
					z = make_cuFloatComplex(cxs[k], cys[k]);

					for(s=0; s<=term_num; s++){
						for(t=s; t<=term_num; t++){
							//L_{s} += L_{t} * C(t, s) * z^(t - s)
							child_node->m_LocalCoeff[s * num_of_ccell + child_cell_idx] =
								cuCaddf(
									child_node->m_LocalCoeff[s * num_of_ccell + child_cell_idx],
									ucuCmulf(
										binomial_coefficientf(t, s),
										cuCmulf(
											node->m_LocalCoeff[t * num_of_cell + cell_idx],
											ucuCpowf(z, t-s))));
						}
					}

					////debug print
					//std::cout << std::endl;
					//for(s=0; s<= term_num; s++){
					//	std::cout << child_node->m_LocalCoeff[s * num_of_ccell + child_cell_idx] << std::endl;
					//}
					//std::cout << std::endl;
				}
			}
		}
	}

	void ucuFMCSM2D_Laplace_CPU_Kernel::M2L_Translation(
		int term_num,
		nmc::ucuFMCSM2DTreeNode *node)
	{
		int i,j,k;
		int s,t;
		cuComplex z;
		int node_w, node_h;
		node_w = node_h = std::pow((float)2, node->m_NodeLevel);
		int num_of_cell = node_w * node_h;

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

		//init zero
		memset(node->m_LocalCoeff, 0, sizeof(float2) * num_of_cell * (term_num+1));

		//
		unsigned int dcell_idx;
		unsigned int scell_idx;
		unsigned int block_idx;
		for(j=0; j<node_h; j++){
			for(i=0; i<node_w; i++){
				dcell_idx = GetCellIndex(i, j, node->m_NodeLevel);
				block_idx = (i & 0x01) + ((j & 0x01) << 1);
				for(k=0; k<27; k++){
					int difx = dif_idx_x[block_idx][k];
					int dify = dif_idx_y[block_idx][k];
					int idx_x = i + difx;
					int idx_y = j + dify;

					if(idx_x < 0 || node_w <= idx_x)
						continue;
					if(idx_y < 0 || node_h <= idx_y)
						continue;

					scell_idx = GetCellIndex(idx_x, idx_y, node->m_NodeLevel);

					////debug print
					//std::cout << std::endl;
					//for(s=0; s<= term_num; s++){
					//	std::cout << node->m_MultipoleCoeff[s * num_of_cell + scell_idx] << std::endl;
					//}
					//std::cout << std::endl;

					cuComplex z = make_cuFloatComplex(
						node->m_CellSize * difx,
						-node->m_CellSize * dify);
					//L_{0} += M_{0} * log( -z )
					node->m_LocalCoeff[0 * num_of_cell + dcell_idx] =
						cuCaddf(
							node->m_LocalCoeff[0 * num_of_cell + dcell_idx],
							ucuCmulf(
								node->m_MultipoleCoeff[0 * num_of_cell + scell_idx],
								logf( cuCabsf( ucuCmulf(-1.0, z) ))));
								//ucuClogf( ucuCmulf(-1.0, z) )));
					for(t=1; t<=term_num; t++){
						float sign = ((t % 2) == 0) ? 1.0 : -1.0;
						//L_{0} += (-1)^t * M_{t} / z^t
						node->m_LocalCoeff[0 * num_of_cell + dcell_idx] = 
							cuCaddf(
								node->m_LocalCoeff[0 * num_of_cell + dcell_idx],
								cuCdivf(
									ucuCmulf(
										sign,
										node->m_MultipoleCoeff[t * num_of_cell + scell_idx]),
									ucuCpowf(z, t)));
					}

					for(s=1; s<=term_num; s++){
						//L_{s} += -M_{0} / (s * z^s)
						node->m_LocalCoeff[s * num_of_cell + dcell_idx] = 
							cuCaddf(
								node->m_LocalCoeff[s * num_of_cell + dcell_idx],
								cuCdivf(
									ucuCmulf(
										-1.0,
										node->m_MultipoleCoeff[0 * num_of_cell + scell_idx]),
									ucuCmulf(
										static_cast<float>(s),
										ucuCpowf(z, s))));
						for(t=1; t<=term_num; t++){
							//L_{s} += (-1)^t * M_{t} * C(s+t-1, t-1) / z^(s+t)
							float sign = ((t % 2) == 0) ? 1.0 : -1.0;
							node->m_LocalCoeff[s * num_of_cell + dcell_idx] =
								cuCaddf(
									node->m_LocalCoeff[s * num_of_cell + dcell_idx],
									cuCdivf(
										ucuCmulf(
											sign * binomial_coefficientf(s+t-1, t-1),
											node->m_MultipoleCoeff[t * num_of_cell + scell_idx]),
										ucuCpowf(z, s+t)));
						}
					}
				}

				////debug print
				//std::cout << std::endl;
				//for(s=0; s<= term_num; s++){
				//	std::cout << node->m_LocalCoeff[s * num_of_cell + dcell_idx] << std::endl;
				//}
				//std::cout << std::endl;
			}
		}
	}

	void ucuFMCSM2D_Laplace_CPU_Kernel::EvalMultipole(
		const ucuFMCSM2DTree* tree,
		const ucuFMCSM2DTreeNode *node,
		float2* dst)
	{
		int s,t;
		int i,j,k;
		int term_num = tree->m_TermNum;
		int node_w, node_h;
		node_w = node_h = std::pow((float)2, node->m_NodeLevel);
		int num_of_cell = node_w * node_h;

		float cell_w = node->m_CellSize;
		float cell_h = node->m_CellSize;
		float min_center_x = -cell_w * (node_w/2) + (cell_w/2.0);
		float min_center_y = -cell_h * (node_h/2) + (cell_h/2.0);

		float x,y;
		float center_x;
		float center_y;

		//init zero
		memset(dst, 0, sizeof(float2) * tree->m_SourceNum/2);

		//
		unsigned int cell_idx;
		for(j=0; j<node_h; j++){
			for(i=0; i<node_w; i++){
				cell_idx = GetCellIndex(i, j, node->m_NodeLevel);
				center_x = i * cell_w + min_center_x;
				center_y = (node_h - j - 1) * cell_h + min_center_y;

				int start_index = tree->m_pStartIndexes[ cell_idx ];
				int end_index = tree->m_pEndIndexes[ cell_idx ];

				for(k=start_index; k<end_index; k++){
					int index = tree->m_pPointIndexes[k];
					if(index < 0)
						continue;

					////debug print
					//std::cout << std::endl;
					//for(s=0; s<= term_num; s++){
					//	std::cout << node->m_LocalCoeff[s * num_of_cell + cell_idx] << std::endl;
					//}
					//std::cout << std::endl;

					cuComplex z = make_cuFloatComplex(
						x = tree->m_pPointsX[k] - center_x,
						y = tree->m_pPointsY[k] - center_y);

					cuComplex tmp = make_cuFloatComplex(0.0, 0.0);
					for(s=0; s<=term_num; s++){
						//L_{s} * z^s
						tmp = cuCaddf(
								tmp,
								cuCmulf(
									node->m_LocalCoeff[s * num_of_cell + cell_idx],
									ucuCpowf(z, s)));
					}
					dst[index] = cuCaddf(dst[index],tmp);
				}
			}
		}
	}

	void ucuFMCSM2D_Laplace_CPU_Kernel::EvalDirect(
		const ucuFMCSM2DTree* tree,
		const ucuFMCSM2DTreeNode *node,
		const float2* cf,
		float2* dst)
	{
		int s,t;
		int i,j,k;
		int term_num = tree->m_TermNum;
		int node_w, node_h;
		node_w = node_h = std::pow((float)2, node->m_NodeLevel);
		int num_of_cell = node_w * node_h;

		int dif_idx_x[9] = {
			-1, 0, 1,
			-1, 0, 1,
			-1, 0, 1,
		};
		int dif_idx_y[9] = {
			-1, -1, -1,
			 0,  0,  0,
			 1,  1,  1,
		};

		//
		unsigned int dcell_idx;
		unsigned int scell_idx;
		int dstart_idx, dend_idx;
		int sstart_idx, send_idx;
		for(j=0; j<node_h; j++){
			for(i=0; i<node_w; i++){
				dcell_idx = GetCellIndex(i, j, node->m_NodeLevel);
				dstart_idx = tree->m_pStartIndexes[ dcell_idx ];
				dend_idx = tree->m_pEndIndexes[ dcell_idx ];

				for(k=0; k<9; k++){
					int idx_x = i + dif_idx_x[k];
					int idx_y = j + dif_idx_y[k];

					if(idx_x < 0 || node_w <= idx_x)
						continue;
					if(idx_y < 0 || node_h <= idx_y)
						continue;

					scell_idx = GetCellIndex(idx_x, idx_y, node->m_NodeLevel);
					sstart_idx = tree->m_pStartIndexes[ scell_idx ];
					send_idx = tree->m_pEndIndexes[ scell_idx ];

					for(s=dstart_idx; s<dend_idx; s++){
						int d_index = tree->m_pPointIndexes[s];
						if(d_index < 0)
							continue;

						for(t=sstart_idx; t<send_idx; t++){
							int s_index = tree->m_pPointIndexes[t];
							if(s_index >= 0)
								continue;

							cuComplex z = make_cuFloatComplex(
								tree->m_pPointsX[s] - tree->m_pPointsX[t],
								tree->m_pPointsY[s] - tree->m_pPointsY[t]);

							dst[d_index] = 
								cuCaddf(
									dst[d_index],
									cuCmulf(
										//ucuCmulf(0.5, cf[-s_index-1]),
										cf[-s_index-1],
										ucuClogf(z)));
										//logf( cuCabsf(z) )));
						}
					}
				}
			}
		}
	}
};

#endif //_NMC_UCU_FMCSM2D_LAPLACE_CPU_KERNEL_H_