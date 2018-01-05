#ifndef _NMC_UCU_FMCSM3D_H_
#define _NMC_UCU_FMCSM3D_H_

#include "nmc/Vector3D.h"
#include "nmc/ucuFMCSM3DTree.h"
#include "nmc/ucuFMCSM3DKernel.h"

#include <algorithm>
#include <assert.h>


namespace nmc
{
	class ucuFMCSM3D{
	public:
		//////////////////////////////////////////
		/*!
			Create FMMTree from source points and bundary points
			@param[in] term_num order of multipole expansion
			@param[in] max_level max level of fmm tree
			@param[in] source_points
			@param[in] bundary_points
			@return created FMMTree created FMMTree
		*/
		static ucuFMCSM3DTree* CreateTree(
			int term_num,
			int max_level,
			const Vector< Vector3Dd >& source_points,
			const Vector< Vector3Dd >& collocation_points);

		//////////////////////////////////////////
		/*!
			Destory FMMTree
			@param[in,out] tree
		*/
		static void DestoryTree(ucuFMCSM3DTree* tree);

		//////////////////////////////////////////////////
		/*!
			calc A*cf
			@param[in,out] tree fmm tree to be used
			@param[in] kernel 
			@param[in] cf 
			@param[in,out] Ax matrix vector product
		*/
		static void Evaluate(
			ucuFMCSM3DTree* tree,
			ucuFMCSM3DKernel& kernel,
			const Vector<float2>& cf,
			Vector<float2>& Ax);

	private:

		//////////////////////////////////////////
		/*!
			create tree node recursively
			@param[in] tree_node
		*/
		static void CreateTreeNode(
			ucuFMCSM3DTree* tree,
			ucuFMCSM3DTreeNode* node,
			int last_level,
			float center_x,
			float center_y,
			float center_z,
			int cell_index_x,
			int cell_index_y,
			int cell_index_z,
			int node_level,
			int start_index,
			int end_index);

		//////////////////////////////////////////
		/*!
		*/
		static void Sort(
			float* main_points,
			float* sub_points1,
			float* sub_points2,
			int* point_indexs,
			int start_index, int end_index);

		//////////////////////////////////////////
		/*!
		*/
		static int GetMidIndex(
			const float center_val,
			const float* points,
			int start_index,
			int end_index);

		//////////////////////////////////////////
		/*!
			set root box size to fmm_tree
			@param[in] tree
		*/
		static void SetRootBoxSize(ucuFMCSM3DTree* tree);

		static void UpwardTree(ucuFMCSM3DKernel& kernel,
			const ucuFMCSM3DTree* tree, const Vector<float2>& cf);

		static void DownwardTree(ucuFMCSM3DKernel& kernel,
			ucuFMCSM3DTree* tree, const Vector<float2>& cf,
			Vector<float2>& Ax);
	};

	ucuFMCSM3DTree* ucuFMCSM3D::CreateTree(
		int term_num,
		int max_level,
		const nmc::Vector<Vector3Dd> &source_points,
		const nmc::Vector<Vector3Dd> &collocation_points)
	{
		//create tree
		ucuFMCSM3DTree* tree = ucuFMCSM3DTree::CreateTree(
			term_num,
			max_level,
			source_points,
			collocation_points);

		//create root node
		tree->m_pRootTreeNode = new ucuFMCSM3DTreeNode;
		tree->m_pRootTreeNode->m_NodeLevel = 0;
		tree->m_pRootTreeNode->m_pParentNode = NULL;

		SetRootBoxSize( tree );

		//create tree node
		int i;
		ucuFMCSM3DTreeNode* tree_node = tree->m_pRootTreeNode;
		for(i = 0; i < max_level; i++){
			tree_node = ucuFMCSM3DTreeNode::CreateNextNode(tree_node, tree->m_TermNum);
		}
		tree->m_pLeafTreeNode = tree_node;

		//
		CreateTreeNode(
			tree,
			tree->m_pRootTreeNode,
			max_level,
			0.0, //center_x
			0.0, //center_y
			0.0, //center_z
			0, //cell_index_x
			0, //cell_index_y
			0, //cell_index_z
			0, //node_level
			0, //start index
			tree->m_SourceNum); //end index

		//copy coordinate to gpu memory
		tree->CopyCoordToDevice();
		return tree;
	}

	void ucuFMCSM3D::CreateTreeNode(
			ucuFMCSM3DTree* tree,
			ucuFMCSM3DTreeNode* node,
			int max_level,
			float center_x,
			float center_y,
			float center_z,
			int cell_index_x,
			int cell_index_y,
			int cell_index_z,
			int node_level,
			int start_index,
			int end_index)
	{
		if(max_level == node_level){ //reached leaf
			unsigned int cell_index = GetCellIndex(
				cell_index_x,
				cell_index_y,
				cell_index_z,
				node_level);
			tree->m_pStartIndexes[ cell_index ] = start_index;
			tree->m_pEndIndexes[ cell_index ] = end_index;

			return;
		}

		//sort by z
		Sort(
			tree->m_pPointsZ,
			tree->m_pPointsX,
			tree->m_pPointsY,
			tree->m_pPointIndexes,
			start_index,
			end_index);

		//get middle value index by z
		int mid_z_index = GetMidIndex(
			center_z,
			tree->m_pPointsZ,
			start_index,
			end_index);

		if(start_index < mid_z_index){
			//sort by y
			Sort(
				tree->m_pPointsY,
				tree->m_pPointsZ,
				tree->m_pPointsX,
				tree->m_pPointIndexes,
				start_index,
				mid_z_index);

			//get middle value index by y
			int mid_y_index = GetMidIndex(
				center_y,
				tree->m_pPointsY,
				start_index,
				mid_z_index);

			if(start_index < mid_y_index){
				//sort by x
				Sort(
					tree->m_pPointsX,
					tree->m_pPointsY,
					tree->m_pPointsZ,
					tree->m_pPointIndexes,
					start_index,
					mid_y_index);

				//get middle value index by x
				int mid_x_index = GetMidIndex(
					center_x,
					tree->m_pPointsX,
					start_index,
					mid_y_index);

				if(start_index < mid_x_index){
					CreateTreeNode(
						tree,
						node->m_pNextTreeNode,
						max_level,
						center_x - node->m_pNextTreeNode->m_CellSize/2.0,
						center_y - node->m_pNextTreeNode->m_CellSize/2.0,
						center_z - node->m_pNextTreeNode->m_CellSize/2.0,
						(cell_index_x << 1) + 0,
						(cell_index_y << 1) + 1,
						(cell_index_z << 1) + 1,
						node_level + 1,
						start_index,
						mid_x_index);
				}
				if(mid_x_index < mid_y_index){
					CreateTreeNode(
						tree,
						node->m_pNextTreeNode,
						max_level,
						center_x + node->m_pNextTreeNode->m_CellSize/2.0,
						center_y - node->m_pNextTreeNode->m_CellSize/2.0,
						center_z - node->m_pNextTreeNode->m_CellSize/2.0,
						(cell_index_x << 1) + 1,
						(cell_index_y << 1) + 1,
						(cell_index_z << 1) + 1,
						node_level + 1,
						mid_x_index,
						mid_y_index);
				}
			}
			if(mid_y_index < mid_z_index){
				//sort by x
				Sort(
					tree->m_pPointsX,
					tree->m_pPointsY,
					tree->m_pPointsZ,
					tree->m_pPointIndexes,
					mid_y_index,
					mid_z_index);

				//get middle value index by x
				int mid_x_index = GetMidIndex(
					center_x,
					tree->m_pPointsX,
					mid_y_index,
					mid_z_index);

				if(mid_y_index < mid_x_index){
					CreateTreeNode(
						tree,
						node->m_pNextTreeNode,
						max_level,
						center_x - node->m_pNextTreeNode->m_CellSize/2.0,
						center_y + node->m_pNextTreeNode->m_CellSize/2.0,
						center_z - node->m_pNextTreeNode->m_CellSize/2.0,
						(cell_index_x << 1) + 0,
						(cell_index_y << 1) + 0,
						(cell_index_z << 1) + 1,
						node_level + 1,
						mid_y_index,
						mid_x_index);
				}
				if(mid_x_index < mid_z_index){
					CreateTreeNode(
						tree,
						node->m_pNextTreeNode,
						max_level,
						center_x + node->m_pNextTreeNode->m_CellSize/2.0,
						center_y + node->m_pNextTreeNode->m_CellSize/2.0,
						center_z - node->m_pNextTreeNode->m_CellSize/2.0,
						(cell_index_x << 1) + 1,
						(cell_index_y << 1) + 0,
						(cell_index_z << 1) + 1,
						node_level + 1,
						mid_x_index,
						mid_z_index);
				}
			}
		}

		if(mid_z_index < end_index){
			//sort by y
			Sort(
				tree->m_pPointsY,
				tree->m_pPointsZ,
				tree->m_pPointsX,
				tree->m_pPointIndexes,
				mid_z_index,
				end_index);

			//get middle value index by y
			int mid_y_index = GetMidIndex(
				center_y,
				tree->m_pPointsY,
				mid_z_index,
				end_index);

			if(mid_z_index < mid_y_index){
				//sort by x
				Sort(
					tree->m_pPointsX,
					tree->m_pPointsY,
					tree->m_pPointsZ,
					tree->m_pPointIndexes,
					mid_z_index,
					mid_y_index);

				//get middle value index by x
				int mid_x_index = GetMidIndex(
					center_x,
					tree->m_pPointsX,
					mid_z_index,
					mid_y_index);

				if(mid_z_index < mid_x_index){
					CreateTreeNode(
						tree,
						node->m_pNextTreeNode,
						max_level,
						center_x - node->m_pNextTreeNode->m_CellSize/2.0,
						center_y - node->m_pNextTreeNode->m_CellSize/2.0,
						center_z + node->m_pNextTreeNode->m_CellSize/2.0,
						(cell_index_x << 1) + 0,
						(cell_index_y << 1) + 1,
						(cell_index_z << 1) + 0,
						node_level + 1,
						mid_z_index,
						mid_x_index);
				}
				if(mid_x_index < mid_y_index){
					CreateTreeNode(
						tree,
						node->m_pNextTreeNode,
						max_level,
						center_x + node->m_pNextTreeNode->m_CellSize/2.0,
						center_y - node->m_pNextTreeNode->m_CellSize/2.0,
						center_z + node->m_pNextTreeNode->m_CellSize/2.0,
						(cell_index_x << 1) + 1,
						(cell_index_y << 1) + 1,
						(cell_index_z << 1) + 0,
						node_level + 1,
						mid_x_index,
						mid_y_index);
				}
			}
			if(mid_y_index < end_index){
				//sort by x
				Sort(
					tree->m_pPointsX,
					tree->m_pPointsY,
					tree->m_pPointsZ,
					tree->m_pPointIndexes,
					mid_y_index,
					end_index);

				//get middle value index by x
				int mid_x_index = GetMidIndex(
					center_x,
					tree->m_pPointsX,
					mid_y_index,
					end_index);

				if(mid_y_index < mid_x_index){
					CreateTreeNode(
						tree,
						node->m_pNextTreeNode,
						max_level,
						center_x - node->m_pNextTreeNode->m_CellSize/2.0,
						center_y + node->m_pNextTreeNode->m_CellSize/2.0,
						center_z + node->m_pNextTreeNode->m_CellSize/2.0,
						(cell_index_x << 1) + 0,
						(cell_index_y << 1) + 0,
						(cell_index_z << 1) + 0,
						node_level + 1,
						mid_y_index,
						mid_x_index);
				}
				if(mid_x_index < end_index){
					CreateTreeNode(
						tree,
						node->m_pNextTreeNode,
						max_level,
						center_x + node->m_pNextTreeNode->m_CellSize/2.0,
						center_y + node->m_pNextTreeNode->m_CellSize/2.0,
						center_z + node->m_pNextTreeNode->m_CellSize/2.0,
						(cell_index_x << 1) + 1,
						(cell_index_y << 1) + 0,
						(cell_index_z << 1) + 0,
						node_level + 1,
						mid_x_index,
						end_index);
				}
			}
		}
	}

	void ucuFMCSM3D::Sort(
		float* main_points,
		float* sub_points1,
		float* sub_points2,
		int* point_indexs,
		int start_index,
		int end_index)
	{
		assert(start_index < end_index);
		if(end_index - start_index == 1) return;

		int i,j;
		float mid_val = main_points[(end_index + start_index)/2];

		i = start_index;
		j = end_index - 1;
		for(;;){
			while(main_points[i] < mid_val) i++;
			while(main_points[j] > mid_val) j--;

			if(j >= i)
				break;

			std::swap(main_points[i], main_points[j]);
			std::swap(sub_points1[i], sub_points1[j]);
			std::swap(sub_points2[i], sub_points2[j]);
			std::swap(point_indexs[i], point_indexs[j]);
			i++; j--;
		}

		if(start_index < i){
			Sort(
				main_points,
				sub_points1,
				sub_points2,
				point_indexs,
				start_index,
				i);
		}
		if(j+1 < end_index){
			Sort(
				main_points,
				sub_points1,
				sub_points2,
				point_indexs,
				j+1,
				end_index);
		}
	}

	int ucuFMCSM3D::GetMidIndex(
		float center_val,
		const float* points,
		int start_index,
		int end_index)
	{
		int i;
		float dife = points[end_index-1] - center_val;
		float difs = center_val - points[start_index];
		if(dife > difs){
			for(i = start_index; i < end_index; i++){
				if(points[i] >= center_val){
					return i;
				}
			}
			return end_index + 1;
		}
		else{
			for(i = end_index-1; i >= start_index; i--){
				if(points[i] < center_val){
					return i+1;
				}
			}
		}
		return -1;
	}

	void ucuFMCSM3D::SetRootBoxSize(
		ucuFMCSM3DTree* tree)
	{
		float x, y, z;
		float max_x, min_x;
		float max_y, min_y;
		float max_z, min_z;
		float dif_x;
		float dif_y;
		float dif_z;

		max_x = min_x = tree->m_pPointsX[0];
		max_y = min_y = tree->m_pPointsY[0];
		max_z = min_z = tree->m_pPointsZ[0];

		int i;
		int n = tree->m_SourceNum;
		for(i = 1; i < n; i++){
			x = tree->m_pPointsX[i];
			y = tree->m_pPointsY[i];
			z = tree->m_pPointsZ[i];

			if(max_x < x) max_x = x;
			if(min_x > x) min_x = x;

			if(max_y < y) max_y = y;
			if(min_y > y) min_y = y;

			if(max_z < z) max_z = z;
			if(min_z > z) min_z = z;
		}

		dif_x = max_x - min_x;
		dif_y = max_y - min_y;
		dif_z = max_z - min_z;
		tree->m_pRootTreeNode->m_CellSize = max(dif_x, max(dif_y, dif_z));
	}

	void ucuFMCSM3D::UpwardTree(
		ucuFMCSM3DKernel& kernel,
		const ucuFMCSM3DTree* tree,
		const Vector<float2>& cf)
	{
		ucuFMCSM3DTreeNode* pNode = tree->m_pLeafTreeNode;

#ifdef _UCU_FMCSM_MEASURE_TIME_
		QPCTimer timer;
#endif //_UCU_FMCSM_MEASURE_TIME_

		//calc P2M (level: max)
		kernel.ElementMPCefficients(tree, pNode, &(cf[0]));

#ifdef _UCU_FMCSM_MEASURE_TIME_
		g_fmcsm_p2m_time = timer.elapsed();
		timer.restart();		
#endif //_UCU_FMCSM_MEASURE_TIME_

		//calc M2M (level: max-1 ~ 1)
		pNode = pNode->m_pParentNode;
		while(pNode && pNode->m_NodeLevel > 0){
			kernel.M2M_Translation(tree->m_TermNum, pNode);
			pNode = pNode->m_pParentNode;
		}

#ifdef _UCU_FMCSM_MEASURE_TIME_
		g_fmcsm_m2m_time = timer.elapsed();
		timer.restart();		
#endif //_UCU_FMCSM_MEASURE_TIME_
	}

	void ucuFMCSM3D::DownwardTree(
		ucuFMCSM3DKernel& kernel,
		ucuFMCSM3DTree* tree,
		const Vector<float2>& cf,
		Vector<float2>& Ax)
	{
		ucuFMCSM3DTreeNode* pNode = NULL;

#ifdef _UCU_FMCSM_MEASURE_TIME_
		QPCTimer timer;
#endif //_UCU_FMCSM_MEASURE_TIME_

		//calc M2L (level: 1 ~ max)
		pNode = tree->m_pRootTreeNode->m_pNextTreeNode;
		while(pNode != NULL){
			kernel.M2L_Translation(tree->m_TermNum, pNode);
			pNode = pNode->m_pNextTreeNode;
		}

#ifdef _UCU_FMCSM_MEASURE_TIME_
		g_fmcsm_m2l_time = timer.elapsed();
		timer.restart();		
#endif //_UCU_FMCSM_MEASURE_TIME_

		//calc L2L (level: 1 ~ max-1)
		pNode = tree->m_pRootTreeNode->m_pNextTreeNode;
		while(pNode != tree->m_pLeafTreeNode){
			kernel.L2L_Translation(tree->m_TermNum, pNode);
			pNode = pNode->m_pNextTreeNode;
		}

#ifdef _UCU_FMCSM_MEASURE_TIME_
		g_fmcsm_l2l_time = timer.elapsed();
		timer.restart();		
#endif //_UCU_FMCSM_MEASURE_TIME_

		//calc L2P (level: max)
		kernel.EvalMultipole(tree, pNode, &(Ax[0]));

#ifdef _UCU_FMCSM_MEASURE_TIME_
		g_fmcsm_l2p_time = timer.elapsed();
		timer.restart();		
#endif //_UCU_FMCSM_MEASURE_TIME_

		//calc direct potential (level: max)
		kernel.EvalDirect(tree, pNode, &(cf[0]), &(Ax[0]));

#ifdef _UCU_FMCSM_MEASURE_TIME_
		g_fmcsm_direct_time = timer.elapsed();
		timer.restart();		
#endif //_UCU_FMCSM_MEASURE_TIME_
	}
};

#endif //_NMC_UCU_FMCSM3D_H_