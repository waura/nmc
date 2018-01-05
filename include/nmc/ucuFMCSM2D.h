#ifndef _NMC_UCU_FMCSM2D_H_
#define _NMC_UCU_FMCSM2D_H_

#define _UCU_FMCSM_MEASURE_TIME_ 1

#ifdef _UCU_FMCSM_MEASURE_TIME_
	#ifdef _WIN32
		#include "nmc/QPCTimer.h"
	#endif //_WIN32

	double g_fmcsm_p2m_time = 0.0;
	double g_fmcsm_m2m_time = 0.0;
	double g_fmcsm_m2l_time = 0.0;
	double g_fmcsm_l2l_time = 0.0;
	double g_fmcsm_l2p_time = 0.0;
	double g_fmcsm_direct_time = 0.0;
	double g_fmcsm_evaluate_time = 0.0;
#endif //_UCU_FMCSM_MEASURE_TIME_

#include "nmc/Vector2D.h"
#include "nmc/ucuFMCSM2DTree.h"
#include "nmc/ucuFMCSM2DKernel.h"

#include <algorithm>
#include <assert.h>


namespace nmc
{
	class ucuFMCSM2D{
	public:
		//////////////////////////////////////////
		/*!
			Create Morton FMMTree from source points and bundary points
			@param[in] term_num order of multipole expansion
			@param[in] max_level max level of fmm tree
			@param[in] source_points
			@param[in] bundary_points
			@return created FMMTree created FMMTree
		*/
		static ucuFMCSM2DTree* CreateMortonTree(
			int term_num,
			int max_level,
			const Vector< Vector2Dd >& source_points,
			const Vector< Vector2Dd >& bundary_points);

		//////////////////////////////////////////
		/*!
			Create FMMTree from source points and bundary points
			@param[in] term_num order of multipole expansion
			@param[in] max_level max level of fmm tree
			@param[in] source_points
			@param[in] bundary_points
			@return created FMMTree created FMMTree
		*/
		static ucuFMCSM2DTree* CreateTree(
			int term_num,
			int max_level,
			const Vector< Vector2Dd >& source_points,
			const Vector< Vector2Dd >& bundary_points);

		//////////////////////////////////////////
		/*!
			Destory FMMTree
			@param[in,out] tree
		*/
		static void DestoryTree(ucuFMCSM2DTree* tree);

		//////////////////////////////////////////////////
		/*!
			calc A*cf
			@param[in,out] tree fmm tree to be used
			@param[in] kernel 
			@param[in] cf 
			@param[in,out] Ax matrix vector product
		*/
		static void Evaluate(
			ucuFMCSM2DTree* tree,
			ucuFMCSM2DKernel& kernel,
			const Vector<float2>& cf,
			Vector<float2>& Ax);

	private:

		//////////////////////////////////////////
		/*!
			create morton tree node recursively
			@param[in] tree_node
		*/
		static void CreateMortonTreeNode(
			ucuFMCSM2DTree* fmm_tree,
			ucuFMCSM2DTreeNode* tree_node,
			int last_level,
			float center_x,
			float center_y,
			int cell_index,
			int start_index,
			int end_index);

		//////////////////////////////////////////
		/*!
			create tree node recursively
			@param[in] tree_node
		*/
		static void CreateTreeNode(
			ucuFMCSM2DTree* tree,
			ucuFMCSM2DTreeNode* node,
			int max_level,
			float center_x,
			float center_y,
			int cell_index_x,
			int cell_index_y,
			int node_level,
			int start_index,
			int end_index);

		//////////////////////////////////////////
		/*!
		*/
		static void Sort(
			float* main_points,
			float* sub_points, 
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
			@param[in] fmm_tree
		*/
		static void SetRootBoxSize(ucuFMCSM2DTree* tree);

		static void UpwardTree(ucuFMCSM2DKernel& kernel,
			const ucuFMCSM2DTree* tree, const Vector<float2>& cf);

		static void DownwardTree(ucuFMCSM2DKernel& kernel,
			ucuFMCSM2DTree* tree, const Vector<float2>& cf,
			Vector<float2>& Ax);

	};

	ucuFMCSM2DTree* ucuFMCSM2D::CreateMortonTree(
		int term_num,
		int max_level,
		const Vector< Vector2Dd > &source_points,
		const Vector< Vector2Dd > &collocation_points)
	{
		//create tree
		ucuFMCSM2DTree* tree = ucuFMCSM2DTree::CreateTree(
			term_num,
			max_level,
			source_points,
			collocation_points);

		//create root nodes
		tree->m_pRootTreeNode = new ucuFMCSM2DTreeNode;
		tree->m_pRootTreeNode->m_NodeLevel = 0;
		tree->m_pRootTreeNode->m_pParentNode = NULL;

		SetRootBoxSize( tree );

		//create tree node
		int i;
		ucuFMCSM2DTreeNode* tree_node = tree->m_pRootTreeNode;
		for(i = 0; i < max_level; i++){
			tree_node = ucuFMCSM2DTreeNode::CreateNextNode(tree_node, tree->m_TermNum);
		}
		tree->m_pLeafTreeNode = tree_node;

		//
		CreateMortonTreeNode(
			tree,
			tree->m_pRootTreeNode,
			max_level,
			0.0, //center_x
			0.0, //center_y
			0, //cell_index
			0, //start index
			tree->m_SourceNum); //end index

		//copy coordinate to gpu memory
		tree->CopyCoordToDevice();
		return tree;
 	}

	ucuFMCSM2DTree* ucuFMCSM2D::CreateTree(
		int term_num,
		int max_level,
		const Vector< Vector2Dd > &source_points,
		const Vector< Vector2Dd > &collocation_points)
	{
		//create tree
		ucuFMCSM2DTree* tree = ucuFMCSM2DTree::CreateTree(
			term_num,
			max_level,
			source_points,
			collocation_points);

		//create nodes
		tree->m_pRootTreeNode = new ucuFMCSM2DTreeNode;
		tree->m_pRootTreeNode->m_NodeLevel = 0;
		tree->m_pRootTreeNode->m_pParentNode = NULL;

		SetRootBoxSize( tree );

		//create tree node
		int i;
		ucuFMCSM2DTreeNode* tree_node = tree->m_pRootTreeNode;
		for(i = 0; i < max_level; i++){
			tree_node = ucuFMCSM2DTreeNode::CreateNextNode(tree_node, tree->m_TermNum);
		}
		tree->m_pLeafTreeNode = tree_node;

		//
		CreateTreeNode(
			tree,
			tree->m_pRootTreeNode,
			max_level,
			0.0, //center_x
			0.0, //center_y
			0, //cell_index_x
			0, //cell_index_y
			0, //node_level;
			0, //start index
			tree->m_SourceNum); //end index


		//copy coordinate to gpu memory
		tree->CopyCoordToDevice();
		return tree;
 	}

	void ucuFMCSM2D::DestoryTree(
		ucuFMCSM2DTree* tree)
	{
		ucuFMCSM2DTreeNode *pNode, *tmp;
		pNode = tree->m_pRootTreeNode;
		while(pNode != NULL){
			tmp = pNode;
			pNode = pNode->m_pNextTreeNode;

			delete tmp;
		}
		delete tree;
	}

	void ucuFMCSM2D::Evaluate(
		ucuFMCSM2DTree* tree,
		ucuFMCSM2DKernel& kernel,
		const Vector<float2>& cf,
		Vector<float2>& Ax)
	{
#ifdef _UCU_FMCSM_MEASURE_TIME_
		QPCTimer timer;
#endif //_UCU_FMCSM_MEASURE_TIME_

		//
		UpwardTree(kernel, tree, cf);
		//
		DownwardTree(kernel, tree, cf, Ax);

#ifdef _UCU_FMCSM_MEASURE_TIME_
		g_fmcsm_evaluate_time = timer.elapsed();
#endif //_UCU_FMCSM_MEASURE_TIME_
	}

	void ucuFMCSM2D::CreateMortonTreeNode(
		ucuFMCSM2DTree* tree,
		ucuFMCSM2DTreeNode* tree_node,
		int last_level,
		float center_x,
		float center_y,
		int cell_index,
		int start_index,
		int end_index)
	{
		if(last_level == 0){ //reached leaf
			//
			tree->m_pStartIndexes[ cell_index ] = start_index;
			tree->m_pEndIndexes[ cell_index ] = end_index;

			////debug print
			//std::cout << "cell_index = " << cell_index <<
			//	" s = " << tree->m_pPointIndexes[start_index] <<
			//	" e = " << tree->m_pPointIndexes[end_index-1] << 
			//	" x = " << tree->m_pPointsX[start_index] <<
			//	" y = " << tree->m_pPointsY[start_index] << std::endl;
			
			return;
		}

		//sort by y
		Sort(
			tree->m_pPointsY,
			tree->m_pPointsX,
			tree->m_pPointIndexes,
			start_index,
			end_index);

		//get middle value index by y
		int mid_y_index = GetMidIndex(
			center_y,
			tree->m_pPointsY,
			start_index,
			end_index);

		if(start_index < mid_y_index){
			//sort bottom points by x
			Sort(
				tree->m_pPointsX,
				tree->m_pPointsY,
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
				//go to left bottom cell
				CreateMortonTreeNode(
					tree,
					tree_node->m_pNextTreeNode,
					last_level - 1,
					center_x - tree_node->m_pNextTreeNode->m_CellSize/2.0,
					center_y - tree_node->m_pNextTreeNode->m_CellSize/2.0,
					GetMortonChildCellIndex(cell_index, 2),
					start_index,
					mid_x_index);
			}
			if(mid_x_index < mid_y_index){
				//go to right bottom
				CreateMortonTreeNode(
					tree,
					tree_node->m_pNextTreeNode,
					last_level - 1,
					center_x + tree_node->m_pNextTreeNode->m_CellSize/2.0,
					center_y - tree_node->m_pNextTreeNode->m_CellSize/2.0,
					GetMortonChildCellIndex(cell_index, 3),
					mid_x_index,
					mid_y_index);
			}
		}

		if(mid_y_index < end_index){
			//sort top points by x
			Sort(
				tree->m_pPointsX,
				tree->m_pPointsY,
				tree->m_pPointIndexes,
				mid_y_index,
				end_index);

			//get middle value index by y
			int mid_x_index = GetMidIndex(
				center_x,
				tree->m_pPointsX,
				mid_y_index,
				end_index);

			if(mid_y_index < mid_x_index){
				//go to left top cell
				CreateMortonTreeNode(
					tree,
					tree_node->m_pNextTreeNode,
					last_level - 1,
					center_x - tree_node->m_pNextTreeNode->m_CellSize/2.0,
					center_y + tree_node->m_pNextTreeNode->m_CellSize/2.0,
					GetMortonChildCellIndex(cell_index, 0),
					mid_y_index,
					mid_x_index);
			}
			if(mid_x_index < end_index){
				//go to right top cell
				CreateMortonTreeNode(
					tree,
					tree_node->m_pNextTreeNode,
					last_level - 1,
					center_x + tree_node->m_pNextTreeNode->m_CellSize/2.0,
					center_y + tree_node->m_pNextTreeNode->m_CellSize/2.0,
					GetMortonChildCellIndex(cell_index, 1),
					mid_x_index,
					end_index);
			}
		}
	}

	void ucuFMCSM2D::CreateTreeNode(
			ucuFMCSM2DTree* tree,
			ucuFMCSM2DTreeNode* tree_node,
			int max_level,
			float center_x,
			float center_y,
			int cell_index_x,
			int cell_index_y,
			int node_level,
			int start_index,
			int end_index)
	{
		if(max_level == node_level){ //reached leaf
			//
			unsigned int cell_index = GetCellIndex(
				cell_index_x,
				cell_index_y,
				node_level);
			tree->m_pStartIndexes[ cell_index ] = start_index;
			tree->m_pEndIndexes[ cell_index ] = end_index;

			////debug print
			//std::cout << "cell_index = " << cell_index <<
			//	" s = " << tree->m_pPointIndexes[start_index] <<
			//	" e = " << tree->m_pPointIndexes[end_index-1] << 
			//	" x = " << tree->m_pPointsX[start_index] <<
			//	" y = " << tree->m_pPointsY[start_index] << std::endl;
			
			return;
		}

		//sort by y
		Sort(
			tree->m_pPointsY,
			tree->m_pPointsX,
			tree->m_pPointIndexes,
			start_index,
			end_index);

		//get middle value index by y
		int mid_y_index = GetMidIndex(
			center_y,
			tree->m_pPointsY,
			start_index,
			end_index);

		if(start_index < mid_y_index){
			//sort bottom points by x
			Sort(
				tree->m_pPointsX,
				tree->m_pPointsY,
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
				//go to left bottom cell
				CreateTreeNode(
					tree,
					tree_node->m_pNextTreeNode,
					max_level,
					center_x - tree_node->m_pNextTreeNode->m_CellSize/2.0,
					center_y - tree_node->m_pNextTreeNode->m_CellSize/2.0,
					(cell_index_x << 1) + 0,
					(cell_index_y << 1) + 1,
					node_level + 1,
					start_index,
					mid_x_index);
			}
			if(mid_x_index < mid_y_index){
				//go to right bottom
				CreateTreeNode(
					tree,
					tree_node->m_pNextTreeNode,
					max_level,
					center_x + tree_node->m_pNextTreeNode->m_CellSize/2.0,
					center_y - tree_node->m_pNextTreeNode->m_CellSize/2.0,
					(cell_index_x << 1) + 1,
					(cell_index_y << 1) + 1,
					node_level + 1,
					mid_x_index,
					mid_y_index);
			}
		}

		if(mid_y_index < end_index){
			//sort top points by x
			Sort(
				tree->m_pPointsX,
				tree->m_pPointsY,
				tree->m_pPointIndexes,
				mid_y_index,
				end_index);

			//get middle value index by y
			int mid_x_index = GetMidIndex(
				center_x,
				tree->m_pPointsX,
				mid_y_index,
				end_index);

			if(mid_y_index < mid_x_index){
				//go to left top cell
				CreateTreeNode(
					tree,
					tree_node->m_pNextTreeNode,
					max_level,
					center_x - tree_node->m_pNextTreeNode->m_CellSize/2.0,
					center_y + tree_node->m_pNextTreeNode->m_CellSize/2.0,
					(cell_index_x << 1) + 0,
					(cell_index_y << 1) + 0,
					node_level + 1,
					mid_y_index,
					mid_x_index);
			}
			if(mid_x_index < end_index){
				//go to right top cell
				CreateTreeNode(
					tree,
					tree_node->m_pNextTreeNode,
					max_level,
					center_x + tree_node->m_pNextTreeNode->m_CellSize/2.0,
					center_y + tree_node->m_pNextTreeNode->m_CellSize/2.0,
					(cell_index_x << 1) + 1,
					(cell_index_y << 1) + 0,
					node_level + 1,
					mid_x_index,
					end_index);
			}
		}
	}

	void ucuFMCSM2D::Sort(
		float* main_points,
		float* sub_points,
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

				if(i >= j)
					break;

				std::swap(main_points[i], main_points[j]);
				std::swap(sub_points[i], sub_points[j]);
				std::swap(point_indexs[i], point_indexs[j]);
				i++; j--;
		}

		if(start_index < i){
			Sort(
				main_points,
				sub_points,
				point_indexs,
				start_index,
				i);
		}
		if(j+1 < end_index){
			Sort(
				main_points,
				sub_points,
				point_indexs,
				j+1,
				end_index);
		}
	}

	int ucuFMCSM2D::GetMidIndex(
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

	void ucuFMCSM2D::SetRootBoxSize(
		ucuFMCSM2DTree* tree)
	{
		float x, y;
		float max_x, min_x;
		float max_y, min_y;
		float dif_x;
		float dif_y;

		max_x = min_x = tree->m_pPointsX[0];
		max_y = min_y = tree->m_pPointsY[0];

		int i;
		int n = tree->m_SourceNum;
		for(i = 1; i < n; i++){
			x = tree->m_pPointsX[i];
			y = tree->m_pPointsY[i];

			if(max_x < x) max_x = x;
			if(min_x > x) min_x = x;

			if(max_y < y) max_y = y;
			if(min_y > y) min_y = y;
		}

		dif_x = max_x - min_x;
		dif_y = max_y - min_y;
		tree->m_pRootTreeNode->m_CellSize = max(dif_x, dif_y);
	}

	void ucuFMCSM2D::UpwardTree(
		ucuFMCSM2DKernel& kernel,
		const ucuFMCSM2DTree* tree,
		const Vector<float2>& cf)
	{
		ucuFMCSM2DTreeNode* pNode = tree->m_pLeafTreeNode;

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

	void ucuFMCSM2D::DownwardTree(
		ucuFMCSM2DKernel& kernel,
		ucuFMCSM2DTree* tree,
		const Vector<float2>& cf,
		Vector<float2>& Ax)
	{
		ucuFMCSM2DTreeNode* pNode = NULL;

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

#endif //_NMC_UCU_FMCSM2D_H_