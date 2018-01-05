#ifndef _NMC_UCU_FMCSM2D_MORTON_LAPLACE_CPU_KERNEL_H_
#define _NMC_UCU_FMCSM2D_MORTON_LAPLACE_CPU_KERNEL_H_

#include "nmc/ucuFMCSM2DKernel.h"
#include "nmc/ucuComplex.h"

namespace nmc
{
	class ucuFMCSM2D_Morton_Laplace_CPU_Kernel :
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
		static inline double binomial_coefficientf(int n, int m)
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
	};

	void ucuFMCSM2D_Morton_Laplace_CPU_Kernel::Init(
		nmc::ucuFMCSM2DTree *tree)
	{
	}

	void ucuFMCSM2D_Morton_Laplace_CPU_Kernel::ElementMPCefficients(
		const ucuFMCSM2DTree* tree,
		ucuFMCSM2DTreeNode *node,
		const float2* c_vec)
	{
		assert(tree->m_pStartIndexes);
		assert(tree->m_pEndIndexes);

		int s,t, k;
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

		unsigned short i,j;

		//init zero
		memset(node->m_MultipoleCoeff, 0, sizeof(float2) * num_of_cell * (term_num+1));

		//
		unsigned int cell_idx;
		for(cell_idx = 0; cell_idx < num_of_cell; cell_idx++){
			GetMortonCellPos(cell_idx, &i, &j);
			center_x = i * cell_w + min_center_x;
			center_y = (node_h - j - 1) * cell_h + min_center_y;

			int start_index = tree->m_pStartIndexes[ cell_idx ];
			int end_index = tree->m_pEndIndexes[ cell_idx ];

			for(k=start_index; k<end_index; k++){
				int source_index = tree->m_pPointIndexes[k];
				if(source_index >= 0)
					continue;

				cuComplex z = make_cuFloatComplex(
					x = tree->m_pPointsX[k] - center_x,
					y = tree->m_pPointsY[k] - center_y);

				cuComplex c = c_vec[ -source_index-1 ];

				//M_{0} += c
				node->m_MultipoleCoeff[0 * num_of_cell + cell_idx] =
					cuCaddf(
						node->m_MultipoleCoeff[0 * num_of_cell + cell_idx],
						c);
				for(s=1; s<=term_num; s++){
					//M_{s} += -c * (z^s) / s
					node->m_MultipoleCoeff[s * num_of_cell + cell_idx] =
						cuCaddf(
							node->m_MultipoleCoeff[s * num_of_cell + cell_idx],
							ucuCmulf(
								-1.0/s,
								cuCmulf(
									c, ucuCpowf(z, s))));

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
		//std::cout << std::endl;
		//for(int s=0; s<= tree->m_TermNum; s++){
		//	std::cout << node->m_MultipoleCoeff[s * num_of_cell + 245] << std::endl;
		//}
		//std::cout << std::endl;
		//std::cout << std::endl;
		//for(int s=0; s<= tree->m_TermNum; s++){
		//	std::cout << node->m_MultipoleCoeff[s * num_of_cell + 249] << std::endl;
		//}
	}

	void ucuFMCSM2D_Morton_Laplace_CPU_Kernel::M2M_Translation(
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
		for(j=0; j<node_h; j++){
			for(i=0; i<node_w; i++){
				int cell_idx = i + j*node_w;
				assert(0 <= cell_idx && cell_idx < num_of_cell);

				for(k=0; k<4; k++){
					int child_cell_idx = GetMortonChildCellIndex(cell_idx, k);
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
	}

	void ucuFMCSM2D_Morton_Laplace_CPU_Kernel::L2L_Translation(
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

		for(j=0; j<node_h; j++){
			for(i=0; i<node_w; i++){
				int cell_idx = i+ j*node_w;

				////debug print
				//std::cout << std::endl;
				//for(s=0; s<= term_num; s++){
				//	std::cout << node->m_LocalCoeff[s * num_of_cell + cell_idx] << std::endl;
				//}
				//std::cout << std::endl;

				for(k=0; k<4; k++){
					int child_cell_idx = GetMortonChildCellIndex(cell_idx, k);
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

		////debug print
		//if(node->m_NodeLevel == 4){
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

	void ucuFMCSM2D_Morton_Laplace_CPU_Kernel::M2L_Translation(
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
		for(j=0; j<node_h; j++){
			for(i=0; i<node_w; i++){
				//unsigned int dcell_idx = i+ j*node_w;
				unsigned int dcell_idx = GetMortonCellIndex(i, j);
				unsigned int block_idx = dcell_idx % 4;
				for(k=0; k<27; k++){
					int difx = dif_idx_x[block_idx][k];
					int dify = dif_idx_y[block_idx][k];
					int idx_x = i + difx;
					int idx_y = j + dify;

					if(idx_x < 0 || node_w <= idx_x)
						continue;
					if(idx_y < 0 || node_h <= idx_y)
						continue;

					unsigned int scell_idx = GetMortonCellIndex(idx_x, idx_y);

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

	void ucuFMCSM2D_Morton_Laplace_CPU_Kernel::EvalMultipole(
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
		for(j=0; j<node_h; j++){
			for(i=0; i<node_w; i++){
				//int cell_idx = i + j*node_w;
				int cell_idx = GetMortonCellIndex(i, j);
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

	void ucuFMCSM2D_Morton_Laplace_CPU_Kernel::EvalDirect(
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
		for(j=0; j<node_h; j++){
			for(i=0; i<node_w; i++){
				unsigned int dcell_idx = GetMortonCellIndex(i, j);
				int dstart_idx = tree->m_pStartIndexes[ dcell_idx ];
				int dend_idx = tree->m_pEndIndexes[ dcell_idx ];

				for(k=0; k<9; k++){
					int idx_x = i + dif_idx_x[k];
					int idx_y = j + dif_idx_y[k];

					if(idx_x < 0 || node_w <= idx_x)
						continue;
					if(idx_y < 0 || node_h <= idx_y)
						continue;

					unsigned int scell_idx = GetMortonCellIndex(idx_x, idx_y);
					int sstart_idx = tree->m_pStartIndexes[ scell_idx ];
					int send_idx = tree->m_pEndIndexes[ scell_idx ];

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

#endif //_NMC_UCU_FMCSM2D_MORTON_LAPLACE_CPU_KERNEL_H_