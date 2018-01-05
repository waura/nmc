#ifndef _NMC_FMCSM2D_H_
#define _NMC_FMCSM2D_H_

#include "nmc/FMCSM2DTree.h"
#include "nmc/FMCSM2DKernel.h"

#include <algorithm>

#ifdef _WIN32
	#define _WINDEBUG 1
	#ifdef _WINDEBUG
		#include "QPCTimer.h"
	#endif //_WINDEBUG
#endif //_WIN32

namespace nmc
{

	const int FMCSM_SOLVER_BiCGSTAB = 0;

	template <class VAL_TYPE>
	class FMCSM2D{
	public:

		//////////////////////////////////////////////////
		/*!
			Create FMMTree from source points and bundary points
			@param[in] term_num order of multipole expansion
			@param[in] source_points
			@param[in] bundary_points
			@return created FMMTree created FMMTree
		*/
		static FMCSM2DTree<VAL_TYPE>* CreateTree(
			int term_num,
			const Vector< Vector2Dd >& source_points,
			const Vector< Vector2Dd >& bundary_points);

		//////////////////////////////////////////////////
		/*!
			Destory FMMTree
			@param[in,out] tree
		*/
		static void DestoryTree(FMCSM2DTree<VAL_TYPE>* tree);

		//////////////////////////////////////////////////
		/*!
			calc A*cf
			@param[in,out] tree fmm tree to be used
			@param[in] kernel 
			@param[in] cf 
			@param[in,out] Ax matrix vector product
		*/
		static void Evaluate(
			FMCSM2DTree<VAL_TYPE>* tree,
			FMCSM2DKernel<VAL_TYPE>& kernel,
			const Vector<VAL_TYPE>& cf,
			Vector<VAL_TYPE>& Ax);

		/////////////////////////////////////////////
		/*!
			solve equation
			@param[in] tree
			@param[in] solver_id
			@param[in] r_vec
			@param[out] x_vec
		*/
		static void Solve(
			FMCSM2DTree<VAL_TYPE>* tree,
			FMCSM2DKernel<VAL_TYPE>& kernel,
			int solver_id, double& conv_ratio,
			unsigned int& iteration, Vector<VAL_TYPE>& r_vec, Vector<VAL_TYPE>& x_vec);

		/////////////////////////////////////////////
		/*!
		*/
		static void PlotTreeBox(const FMCSM2DTree<VAL_TYPE>* tree);
		static void PlotSourcePoint(const FMCSM2DTree<VAL_TYPE>* tree);
		static void PlotBundaryPoint(const FMCSM2DTree<VAL_TYPE>* tree);

	private:

		/////////////////////////////////////////////
		/*!
			create tree node recursively
			@param[in] term_num
			@param[in] tree_node
			@param[in] fmm_source_points
			@param[in] start_index
			@param[in] end_index
			@return 
		*/
		static int CreateTreeNode(int term_num,
			FMCSM2DTreeNode<VAL_TYPE>* tree_node,
			Vector< FMMPos >& fmm_source_points, 
			int start_index, int end_index);

		//////////////////////////////////////////////
		/*!
		*/
		static void Sort(Vector< FMMPos >& fmm_source_points, 
			int sort_key, int start_index, int end_index);

		//////////////////////////////////////////////
		/*!
			find index of middle value from fmm_source_points
			@param[in] center_point
			@param[in] fmm_source_points
			@param[in] find_key
			@param[in] start_index
			@param[in] end_index
		*/
		static int GetMidIndex(
			const Vector2Dd& center_point,
			const Vector< FMMPos >& fmm_source_points,
			int find_key, int start_index, int end_index);

		static void SetRootBoxSize(
			FMCSM2DTree<VAL_TYPE>* fmm_tree,
			Vector< FMMPos >& fmm_source_points);

		static void DestoryTreeNode(FMCSM2DTreeNode<VAL_TYPE>* node);

		static void UpwardTree(FMCSM2DKernel<VAL_TYPE>& kernel, 
			FMCSM2DTreeNode<VAL_TYPE>* node, const Vector<VAL_TYPE>& cf);

		static void DownwardTree(FMCSM2DKernel<VAL_TYPE>& kernel, 
			FMCSM2DTreeNode<VAL_TYPE>* node, const Vector<VAL_TYPE>& cf,
			Vector<VAL_TYPE>& Ax);

		static void EvalDirect(FMCSM2DKernel<VAL_TYPE>& kernel,
			const FMCSM2DTreeNode<VAL_TYPE>* from_node,
			const FMCSM2DTreeNode<VAL_TYPE>* to_node,
			const Vector<VAL_TYPE>& cf,
			Vector<VAL_TYPE>& Ax);

		static void GetM2LNode(const FMCSM2DTreeNode<VAL_TYPE>* src_node,
			FMCSM2DTreeNode<VAL_TYPE>** far_parent,
			FMCSM2DTreeNode<VAL_TYPE>** far_node);

		static void GetDirectEvalNode(const FMCSM2DTreeNode<VAL_TYPE>* src_node,
			FMCSM2DTreeNode<VAL_TYPE>** dst_nodes);

		static void PlotTreeNodeBox(const FMCSM2DTreeNode<VAL_TYPE>* node);
		static void PlotNodeSourcePoint(const FMCSM2DTreeNode<VAL_TYPE>* node);
		static void PlotNodeBundaryPoint(const FMCSM2DTreeNode<VAL_TYPE>* node);
	};


	template<class VAL_TYPE>
	FMCSM2DTree<VAL_TYPE>* FMCSM2D<VAL_TYPE>::CreateTree(
		int term_num,
		const Vector< Vector2Dd > &source_points,
		const Vector< Vector2Dd > &bundary_points)
	{
		assert(source_points.dim() == bundary_points.dim());

		Vector<FMMPos> fmm_points;
		fmm_points.Init(source_points.dim()*2);

		//copy to fmm point  
		int i;
		FMMPos fmm_pos;
		for(i=0; i<source_points.dim(); i++){
			fmm_pos.val[X_KEY] = source_points[i][X_KEY];
			fmm_pos.val[Y_KEY] = source_points[i][Y_KEY];
			fmm_pos.index = -i-1;
			fmm_points[i] = fmm_pos;
		}
		for(i=0; i<bundary_points.dim(); i++){
			fmm_pos.val[X_KEY] = bundary_points[i][X_KEY];
			fmm_pos.val[Y_KEY] = bundary_points[i][Y_KEY];
			fmm_pos.index = i;
			fmm_points[i+source_points.dim()] = fmm_pos;
		}

		//create nodes
		FMCSM2DTree<VAL_TYPE>* fmm_tree = new FMCSM2DTree<VAL_TYPE>;
		fmm_tree->m_TermNum = 2*term_num + 1;
		fmm_tree->m_SourceNum = fmm_points.dim();
		fmm_tree->m_pRootTreeNode = new FMCSM2DTreeNode<VAL_TYPE>;
		fmm_tree->m_pRootTreeNode->m_NodeLevel = 0;
		fmm_tree->m_pRootTreeNode->m_NodeIndex = 0;
		fmm_tree->m_pRootTreeNode->m_pParentNode = NULL;

		SetRootBoxSize( fmm_tree, fmm_points );
		fmm_tree->m_MaxNodeLevel = CreateTreeNode(
			2*term_num+1,
			fmm_tree->m_pRootTreeNode,
			fmm_points,
			0,
			fmm_points.dim()-1 );

		return fmm_tree;
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::DestoryTree(FMCSM2DTree<VAL_TYPE>* tree)
	{
		DestoryTreeNode(tree->m_pRootTreeNode);
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::Evaluate(
		FMCSM2DTree<VAL_TYPE>* tree,
		FMCSM2DKernel<VAL_TYPE>& kernel,
		const Vector<VAL_TYPE>& cf,
		Vector<VAL_TYPE>& Ax)
	{
		int i,j;
	#ifdef _WINDEBUG
		char str_lap[256];
		QPCTimer timer;
	#endif
		//init zero clear
		for(i=0; i<Ax.dim(); i++){
			Ax[i] = 0.0;
		}
	#ifdef _WINDEBUG
		sprintf(str_lap, "init Ax zero clear %f msec\n", timer.elapsed());
		::OutputDebugString(str_lap);
		timer.restart();
	#endif
		for(i=0; i<4; i++){
			if(tree->m_pRootTreeNode->m_pNextTreeNodes[i]){
				UpwardTree(kernel, tree->m_pRootTreeNode->m_pNextTreeNodes[i], cf);
			}
		}
	#ifdef _WINDEBUG
		sprintf(str_lap, "UpwardTree %f msec\n", timer.elapsed());
		::OutputDebugString(str_lap);
		timer.restart();
	#endif
		for(i=0; i<4; i++){
			if(tree->m_pRootTreeNode->m_pNextTreeNodes[i] &&
				tree->m_pRootTreeNode->m_pNextTreeNodes[i]->m_isBundaryPosInNode){
				//init local zero clear
				for(j=0; j<tree->m_pRootTreeNode->m_pNextTreeNodes[i]->m_Local.dim(); j++){
					tree->m_pRootTreeNode->m_pNextTreeNodes[i]->m_Local[j] = 0.0;
				}
				DownwardTree(kernel, tree->m_pRootTreeNode->m_pNextTreeNodes[i], cf, Ax);
			}
		}
	#ifdef _WINDEBUG
		sprintf(str_lap, "DownwardTree %f msec\n", timer.elapsed());
		::OutputDebugString(str_lap);
	#endif
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::Solve(
		FMCSM2DTree<VAL_TYPE>* tree,
		FMCSM2DKernel<VAL_TYPE>& kernel,
		int solver_id,
		double& conv_ratio,
		unsigned int& iteration,
		Vector<VAL_TYPE>& r_vec,
		Vector<VAL_TYPE>& x_vec)
	{
		assert(r_vec.dim() == x_vec.dim());

		int i,j;
		double conv_ratio_tol = conv_ratio;
		Vector<VAL_TYPE>  s_vec( r_vec.dim() );
		Vector<VAL_TYPE> As_vec( r_vec.dim() );
		Vector<VAL_TYPE>  p_vec( r_vec.dim() );
		Vector<VAL_TYPE> Ap_vec( r_vec.dim() );
		Vector<VAL_TYPE> r0_conjugate_vec( r_vec.dim() );
		
		VAL_TYPE sq_inv_norm_res;
		{
			VAL_TYPE tmp=0.0;
			for(j=0; j<r_vec.dim(); j++){
				tmp += r_vec[j] * r_vec[j];
			}
			if( tmp.Abs() < 1.0e-30 ){
				conv_ratio = 0.0;
				iteration = 0;
				return;
			}
			sq_inv_norm_res = static_cast<VAL_TYPE>(1.0) / tmp;
		}

		//
		r0_conjugate_vec = r_vec;
		//p_{0} = r_{0}
		p_vec = r_vec;

		for(i=1; i<iteration; i++){
			//calc (r, r0*)
			VAL_TYPE r_r0conj=0.0;
			for(j=0; j<r_vec.dim(); j++){
				r_r0conj += r_vec[j] * r0_conjugate_vec[j];
			}

			//calc Ap
			FMCSM2D::Evaluate(tree, kernel, p_vec, Ap_vec);

			//alpha = (r, r0*)/(Ap, ro*)
			VAL_TYPE alpha;
			{
				VAL_TYPE dnm = 0.0;
				for(j=0; j<Ap_vec.dim(); j++){
					dnm += Ap_vec[j] * r0_conjugate_vec[j];
				}
				alpha = r_r0conj / dnm;
			}

			//s = r - alpha*Ap
			s_vec = r_vec;
			for(j=0; j<s_vec.dim(); j++){
				s_vec[j] -= alpha * Ap_vec[j]; 
			}

			//calc As
			FMCSM2D::Evaluate(tree, kernel, s_vec, As_vec);

			//omega = (s, As)/(As, As)
			VAL_TYPE omega;
			{
				VAL_TYPE numr = 0.0;
				for(j=0; j<s_vec.dim(); j++){
					numr += s_vec[j] * As_vec[j];
				}
				VAL_TYPE dnm = 0.0;
				for(j=0; j<As_vec.dim(); j++){
					dnm += As_vec[j] * As_vec[j];
				}
				omega = numr / dnm;
			}

			//new x = x + alpha * p + omega * s
			for(j=0; j<x_vec.dim(); j++){
				x_vec[j] += alpha * p_vec[j] + omega * s_vec[j];
			}

			//new r = s - omega * As
			for(j=0; j<r_vec.dim(); j++){
				r_vec[j] = s_vec[j] - omega * As_vec[j];
			}

			{
				VAL_TYPE sq_norm_res = 0.0;
				for(j=0; j<r_vec.dim(); j++){
					sq_norm_res += r_vec[j] * r_vec[j];
				}
				VAL_TYPE sq_conv_ratio = sq_norm_res * sq_inv_norm_res;
				if(sq_conv_ratio < conv_ratio_tol*conv_ratio_tol){
					conv_ratio = sqrt( sq_conv_ratio.Abs() );
					iteration = i;
					return;
				}
			}

			//beta = {(new r, r0*)/(old r, r0*)}*{(alpha/omega)} 
			VAL_TYPE beta;
			{
				VAL_TYPE numr = 0.0;
				for(j=0; j<r_vec.dim(); j++){
					numr += r_vec[j] * r0_conjugate_vec[j];
				}
				beta = numr * alpha / (r_r0conj * omega);
			}

			//p = r + beta(p - omega * Ap)
			for(j=0; j<p_vec.dim(); j++){
				p_vec[j] *= beta;
				p_vec[j] += r_vec[j];
				p_vec[j] -= beta * omega * Ap_vec[j];
			}
		}

		return;
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::PlotTreeBox(
		const FMCSM2DTree<VAL_TYPE>* tree)
	{
		PlotTreeNodeBox( tree->m_pRootTreeNode );
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::PlotSourcePoint(
		const FMCSM2DTree<VAL_TYPE>* tree)
	{
		PlotNodeSourcePoint( tree->m_pRootTreeNode );
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::PlotBundaryPoint(
		const FMCSM2DTree<VAL_TYPE>* tree)
	{
		PlotNodeBundaryPoint( tree->m_pRootTreeNode );
	}

	template<class VAL_TYPE>
	int FMCSM2D<VAL_TYPE>::CreateTreeNode(
		int term_num,
		FMCSM2DTreeNode<VAL_TYPE> *tree_node,
		Vector< FMMPos >& fmm_source_points,
		int start_index,
		int end_index)
	{
		int max_node_level = tree_node->m_NodeLevel;
		double next_box_length = tree_node->m_BoxLength/2.0;

		if(start_index == end_index){
			FMMPos pos = fmm_source_points[end_index];
			//printf("index = %d x = %f y = %f\n", end_index, pos.val[X_KEY], pos.val[Y_KEY]);

			if(pos.index < 0){
				pos.index = -pos.index - 1;
				tree_node->m_pSourcePos = new Vector2D<FMMPos>;
				tree_node->m_pSourcePos->Init(1);
				tree_node->m_pSourcePos->SetVal(0, pos);
				tree_node->m_isBundaryPosInNode = false;
				tree_node->m_isSourcePosInNode = true;
			}
			else{
				tree_node->m_pBundaryPos = new Vector2D<FMMPos>;
				tree_node->m_pBundaryPos->Init(1);
				tree_node->m_pBundaryPos->SetVal(0, pos);
				tree_node->m_isBundaryPosInNode = true;
				tree_node->m_isSourcePosInNode = false;
			}
			return max_node_level;
		}

		//sort by y
		Sort( fmm_source_points, Y_KEY, start_index, end_index );
		//get middle value index
		int mid_y_index = GetMidIndex( tree_node->m_CenterPos, fmm_source_points, Y_KEY, start_index, end_index);
		
		if(start_index <= (mid_y_index - 1)){
			//sort bottom points by x
			Sort( fmm_source_points, X_KEY, start_index, mid_y_index-1 );
			int mid_x_index = GetMidIndex(
								tree_node->m_CenterPos,
								fmm_source_points,
								X_KEY,
								start_index,
								mid_y_index-1);
			
			if(start_index <= (mid_x_index-1)){
				tree_node->m_pNextTreeNodes[2] = new FMCSM2DTreeNode<VAL_TYPE>;
				tree_node->m_pNextTreeNodes[2]->Init(
					2,
					tree_node->m_NodeLevel+1,
					term_num,
					next_box_length,
					tree_node->m_CenterPos[X_KEY] - next_box_length/2.0,
					tree_node->m_CenterPos[Y_KEY] - next_box_length/2.0,
					tree_node);

				//create next tree node
				int level;
				level = CreateTreeNode(
					term_num,
					tree_node->m_pNextTreeNodes[2],
					fmm_source_points,
					start_index,
					mid_x_index-1);
				max_node_level = (level > max_node_level) ? level : max_node_level;
			}
			if(mid_x_index <= (mid_y_index-1)){
				tree_node->m_pNextTreeNodes[3] = new FMCSM2DTreeNode<VAL_TYPE>;
				tree_node->m_pNextTreeNodes[3]->Init(
					3,
					tree_node->m_NodeLevel+1,
					term_num,
					next_box_length,
					tree_node->m_CenterPos[X_KEY] + next_box_length/2.0,
					tree_node->m_CenterPos[Y_KEY] - next_box_length/2.0,
					tree_node);

				//create next tree node
				int level;
				level = CreateTreeNode(
					term_num,
					tree_node->m_pNextTreeNodes[3],
					fmm_source_points,
					mid_x_index,
					mid_y_index-1);
				max_node_level = (level > max_node_level) ? level : max_node_level;
			}
		}
		
		if(mid_y_index <= end_index){
			//sort top points by x
			Sort( fmm_source_points, X_KEY, mid_y_index, end_index );
			int mid_x_index = GetMidIndex( tree_node->m_CenterPos, fmm_source_points, X_KEY, mid_y_index, end_index);

			if(mid_y_index <= (mid_x_index-1)){
				tree_node->m_pNextTreeNodes[0] = new FMCSM2DTreeNode<VAL_TYPE>;
				tree_node->m_pNextTreeNodes[0]->Init(
					0,
					tree_node->m_NodeLevel+1,
					term_num,
					next_box_length,
					tree_node->m_CenterPos[X_KEY] - next_box_length/2.0,
					tree_node->m_CenterPos[Y_KEY] + next_box_length/2.0,
					tree_node);

				//create next tree node
				int level;
				level = CreateTreeNode(
					term_num,
					tree_node->m_pNextTreeNodes[0],
					fmm_source_points,
					mid_y_index,
					mid_x_index-1);
				max_node_level = (level > max_node_level) ? level : max_node_level;
			}
			if(mid_x_index <= end_index){
				tree_node->m_pNextTreeNodes[1] = new FMCSM2DTreeNode<VAL_TYPE>;
				tree_node->m_pNextTreeNodes[1]->Init(
					1,
					tree_node->m_NodeLevel+1,
					term_num,
					 next_box_length,
					tree_node->m_CenterPos[X_KEY] + next_box_length/2.0,
					tree_node->m_CenterPos[Y_KEY] + next_box_length/2.0,
					tree_node);

				//create next tree node
				int level;
				level = CreateTreeNode(
					term_num,
					tree_node->m_pNextTreeNodes[1],
					fmm_source_points,
					mid_x_index,
					end_index);
				max_node_level = (level > max_node_level) ? level : max_node_level;
			}


		}
		//
		int i;
		tree_node->m_isBundaryPosInNode = false;
		for(i=0; i<4; i++){
			if(tree_node->m_pNextTreeNodes[i] &&
				tree_node->m_pNextTreeNodes[i]->m_isBundaryPosInNode){
					tree_node->m_isBundaryPosInNode = true;
					break;
			}
		}
		tree_node->m_isSourcePosInNode = false;
		for(i=0; i<4; i++){
			if(tree_node->m_pNextTreeNodes[i] &&
				tree_node->m_pNextTreeNodes[i]->m_isSourcePosInNode){
					tree_node->m_isSourcePosInNode = true;
			}
		}
		return max_node_level;
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::Sort(
		Vector< FMMPos >& fmm_source_points,
		int sort_key,
		int start_index,
		int end_index)
	{
		assert(start_index <= end_index);
		if(start_index == end_index) return;

		int i,j;
		double mid_val = fmm_source_points[(end_index + start_index)/2].val[sort_key];

		i=start_index;
		j=end_index;

		for(;;){
			while(fmm_source_points[i].val[sort_key] < mid_val) i++;
			while(fmm_source_points[j].val[sort_key] > mid_val) j--;

			if(i >= j)
				break;

			std::swap(fmm_source_points[i], fmm_source_points[j]);
			i++; j--;
		}

		if(start_index < i-1)
			Sort(fmm_source_points, sort_key, start_index, i-1);
		if(j+1 < end_index)
			Sort(fmm_source_points, sort_key, j+1, end_index);
	}

	template<class VAL_TYPE>
	int FMCSM2D<VAL_TYPE>::GetMidIndex(
		const Vector2Dd& center_point,
		const Vector<FMMPos> &fmm_source_points,
		int find_key,
		int start_index,
		int end_index)
	{
		double dife = fmm_source_points[end_index].val[find_key] - center_point[find_key];
		double difs = center_point[find_key] - fmm_source_points[start_index].val[find_key];
		if(dife > difs){
			int i;
			for(i=start_index; i<=end_index; i++){
				if(fmm_source_points[i].val[find_key] >= center_point[find_key]){
					return i;
				}
			}
			return end_index+1;
		}
		else{
			int i;
			for(i=end_index; i>=start_index; i--){
				if(fmm_source_points[i].val[find_key] < center_point[find_key]){
					return i+1;
				}
			}
		}
		return -1;
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::SetRootBoxSize(
		FMCSM2DTree<VAL_TYPE>* fmm_tree,
		Vector<FMMPos>& fmm_source_points)
	{
		int i;
		double x, y;
		double max_x, min_x;
		double max_y, min_y;
		double width;
		double height;

		max_x = min_x = fmm_source_points[0].val[X_KEY];
		max_y = min_y = fmm_source_points[0].val[Y_KEY];

		for(i=1; i<fmm_source_points.dim(); i++){
			x = fmm_source_points[i].val[X_KEY];
			y = fmm_source_points[i].val[Y_KEY];

			if(max_x < x) max_x = x;
			if(min_x > x) min_x = x;
			
			if(max_y < y) max_y = y;
			if(min_y > y) min_y = y;
		}

		width = max_x - min_x;
		height = max_y - min_y;

		fmm_tree->m_pRootTreeNode->m_BoxLength = (width > height) ? width : height;
		fmm_tree->m_pRootTreeNode->m_CenterPos.SetVal(0, min_x + (width/2.0));
		fmm_tree->m_pRootTreeNode->m_CenterPos.SetVal(1, min_y + (height/2.0));
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::DestoryTreeNode(FMCSM2DTreeNode<VAL_TYPE> *node)
	{
		if(node == NULL) return;

		int i;
		for(i=0; i<4; i++){
			//
			DestoryTreeNode(node->m_pNextTreeNodes[i]);

			//destory child node
			if(node->m_pNextTreeNodes[i]){
				delete node->m_pNextTreeNodes[i];
				node->m_pNextTreeNodes[i] = NULL;
			}
			if(node->m_pSourcePos){
				delete node->m_pSourcePos;
				node->m_pSourcePos = NULL;
			}
		}
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::UpwardTree(
		FMCSM2DKernel<VAL_TYPE>& kernel,
		FMCSM2DTreeNode<VAL_TYPE>* node,
		const Vector<VAL_TYPE>& cf)
	{
		assert(node->m_NodeLevel >= 1);

		int i;
		//init multipole
		for(i=0; i<node->m_Multipole.dim(); i++){
			node->m_Multipole[i] = 0.0;
		}

		//
		if(node->m_pSourcePos){//this node is leaf
			kernel.ElementMPCefficients(node, cf);
			return;
		}

		for(i=0; i<4; i++){
			if(node->m_pNextTreeNodes[i]){
				UpwardTree(kernel, node->m_pNextTreeNodes[i], cf);

			}
		}

		if(node->m_isSourcePosInNode){
			//calc this node m2m
			kernel.M2M_Translation(node);
		}
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::DownwardTree(
		FMCSM2DKernel<VAL_TYPE>& kernel,
		FMCSM2DTreeNode<VAL_TYPE> *node,
		const Vector<VAL_TYPE>& cf,
		Vector<VAL_TYPE>& Ax)
	{
		assert(node->m_NodeLevel >= 1);

		int i,j,k;

		//calc far node m2l
		FMCSM2DTreeNode<VAL_TYPE>* far_parent[5];
		FMCSM2DTreeNode<VAL_TYPE>* far_node[7];
		GetM2LNode(node, far_parent, far_node);

		for(i=0; i<5; i++){
			if(far_parent[i] == NULL) continue;
			for(j=0; j<4; j++){
				if(far_parent[i]->m_pNextTreeNodes[j] &&
					far_parent[i]->m_pNextTreeNodes[j]->m_isSourcePosInNode){
					kernel.M2L_Translation(far_parent[i]->m_pNextTreeNodes[j], node);
				}
			}
		}
		for(i=0; i<7; i++){
			if(far_node[i] && far_node[i]->m_isSourcePosInNode){
				kernel.M2L_Translation(far_node[i], node);
			}
		}
		//
		//std::cout << node->m_Local << std::endl;

		if(node->m_pBundaryPos){//this node is a leaf
			//multipole evaluation
			kernel.EvalMultipole(node, Ax);

			//direct evaluation from near node
			FMCSM2DTreeNode<VAL_TYPE>* near_node[8];
			GetDirectEvalNode(node, near_node);
			for(i=0; i<8; i++){
				if(near_node[i] == NULL) continue;
				if(near_node[i]->m_isSourcePosInNode == false) continue;
				EvalDirect(kernel, near_node[i], node, cf, Ax);
			}
			if(node->m_pSourcePos && node->m_pBundaryPos){
				for(i=0; node->m_pSourcePos->dim(); i++){
					for(j=0; node->m_pBundaryPos->dim(); j++){
						kernel.EvalDirect(
							node->m_pSourcePos->GetVal(i),
							node->m_pBundaryPos->GetVal(j),
							cf, Ax);
					}
				}
			}
		}
		else{
			//init local zero clear
			for(i=0; i<4; i++){
				if(node->m_pNextTreeNodes[i]){
					for(j=0; j<node->m_pNextTreeNodes[i]->m_Local.dim(); j++){
						node->m_pNextTreeNodes[i]->m_Local[j] = 0.0;
					}
				}
			}

			// calc l2l
			kernel.L2L_Translation(node);

			//go to next node
			for(i=0; i<4; i++){
				if(node->m_pNextTreeNodes[i] && node->m_isBundaryPosInNode){
					DownwardTree(kernel, node->m_pNextTreeNodes[i], cf, Ax);
				}
			}
		}
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::EvalDirect(
		FMCSM2DKernel<VAL_TYPE> &kernel,
		const FMCSM2DTreeNode<VAL_TYPE> *from_node,
		const FMCSM2DTreeNode<VAL_TYPE> *to_node,
		const Vector<VAL_TYPE>& cf,
		Vector<VAL_TYPE>& Ax)
	{
		int i,j;
		if(from_node->m_pSourcePos){
			assert(to_node->m_pBundaryPos);

			for(i=0; i<to_node->m_pBundaryPos->dim(); i++){
				for(j=0; j<from_node->m_pSourcePos->dim(); j++){
					kernel.EvalDirect(
						from_node->m_pSourcePos->GetVal(j),
						to_node->m_pBundaryPos->GetVal(i),
						cf, Ax);
				}
			}
			return;
		}

		for(i=0; i<4; i++){
			if(from_node->m_pNextTreeNodes[i] &&
				from_node->m_pNextTreeNodes[i]->m_isSourcePosInNode){
				EvalDirect(
					kernel,
					from_node->m_pNextTreeNodes[i],
					to_node,
					cf,
					Ax);
			}
		}
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::GetM2LNode(
		const FMCSM2DTreeNode<VAL_TYPE> *src_node,
		FMCSM2DTreeNode<VAL_TYPE> **far_parent,
		FMCSM2DTreeNode<VAL_TYPE> **far_node)
	{
		switch(src_node->m_NodeIndex){
			case 0:
				{				
					far_parent[0] = src_node->GetParentNode()->GetUpperRightNode();
					far_parent[1] = src_node->GetParentNode()->GetRightNode();
					far_parent[2] = src_node->GetParentNode()->GetLowerLeftNode();
					far_parent[3] = src_node->GetParentNode()->GetLowerNode();
					far_parent[4] = src_node->GetParentNode()->GetLowerRightNode();

					FMCSM2DTreeNode<VAL_TYPE>* tmp;
					tmp = src_node->GetParentNode()->GetUpperLeftNode();
					far_node[0] = (tmp) ? tmp->m_pNextTreeNodes[0] : NULL;
					far_node[1] = (tmp) ? tmp->m_pNextTreeNodes[1] : NULL;
					far_node[2] = (tmp) ? tmp->m_pNextTreeNodes[2] : NULL;
					tmp = src_node->GetParentNode()->GetUpperNode();
					far_node[3] = (tmp) ? tmp->m_pNextTreeNodes[0] : NULL;
					far_node[4] = (tmp) ? tmp->m_pNextTreeNodes[1] : NULL;
					tmp = src_node->GetParentNode()->GetLeftNode();
					far_node[5] = (tmp) ? tmp->m_pNextTreeNodes[0] : NULL;
					far_node[6] = (tmp) ? tmp->m_pNextTreeNodes[2] : NULL;
				}
				break;
			case 1:
				{
					far_parent[0] = src_node->GetParentNode()->GetUpperLeftNode();
					far_parent[1] = src_node->GetParentNode()->GetLeftNode();
					far_parent[2] = src_node->GetParentNode()->GetLowerLeftNode();
					far_parent[3] = src_node->GetParentNode()->GetLowerNode();
					far_parent[4] = src_node->GetParentNode()->GetLowerRightNode();

					FMCSM2DTreeNode<VAL_TYPE>* tmp;
					tmp = src_node->GetParentNode()->GetUpperNode();
					far_node[0] = (tmp) ? tmp->m_pNextTreeNodes[0] : NULL;
					far_node[1] = (tmp) ? tmp->m_pNextTreeNodes[1] : NULL;
					tmp = src_node->GetParentNode()->GetUpperRightNode();
					far_node[2] = (tmp) ? tmp->m_pNextTreeNodes[0] : NULL;
					far_node[3] = (tmp) ? tmp->m_pNextTreeNodes[1] : NULL;
					far_node[4] = (tmp) ? tmp->m_pNextTreeNodes[3] : NULL;
					tmp = src_node->GetParentNode()->GetRightNode();
					far_node[5] = (tmp) ? tmp->m_pNextTreeNodes[1] : NULL;
					far_node[6] = (tmp) ? tmp->m_pNextTreeNodes[3] : NULL;
				}
				break;
			case 2:
				{
					far_parent[0] = src_node->GetParentNode()->GetUpperLeftNode();
					far_parent[1] = src_node->GetParentNode()->GetUpperNode();
					far_parent[2] = src_node->GetParentNode()->GetUpperRightNode();
					far_parent[3] = src_node->GetParentNode()->GetRightNode();
					far_parent[4] = src_node->GetParentNode()->GetLowerRightNode();

					FMCSM2DTreeNode<VAL_TYPE>* tmp;
					tmp = src_node->GetParentNode()->GetLeftNode();
					far_node[0] = (tmp) ? tmp->m_pNextTreeNodes[0] : NULL;
					far_node[1] = (tmp) ? tmp->m_pNextTreeNodes[2] : NULL;
					tmp = src_node->GetParentNode()->GetLowerLeftNode();
					far_node[2] = (tmp) ? tmp->m_pNextTreeNodes[0] : NULL;
					far_node[3] = (tmp) ? tmp->m_pNextTreeNodes[2] : NULL;
					far_node[4] = (tmp) ? tmp->m_pNextTreeNodes[3] : NULL;
					tmp = src_node->GetParentNode()->GetLowerNode();
					far_node[5] = (tmp) ? tmp->m_pNextTreeNodes[2] : NULL;
					far_node[6] = (tmp) ? tmp->m_pNextTreeNodes[3] : NULL;
				}
				break;
			case 3:
				{
					far_parent[0] = src_node->GetParentNode()->GetUpperLeftNode();
					far_parent[1] = src_node->GetParentNode()->GetUpperNode();
					far_parent[2] = src_node->GetParentNode()->GetUpperRightNode();
					far_parent[3] = src_node->GetParentNode()->GetLeftNode();
					far_parent[4] = src_node->GetParentNode()->GetLowerLeftNode();

					FMCSM2DTreeNode<VAL_TYPE>* tmp;
					tmp = src_node->GetParentNode()->GetRightNode();
					far_node[0] = (tmp) ? tmp->m_pNextTreeNodes[1] : NULL;
					far_node[1] = (tmp) ? tmp->m_pNextTreeNodes[3] : NULL;
					tmp = src_node->GetParentNode()->GetLowerNode();
					far_node[2] = (tmp) ? tmp->m_pNextTreeNodes[2] : NULL;
					far_node[3] = (tmp) ? tmp->m_pNextTreeNodes[3] : NULL;
					tmp = src_node->GetParentNode()->GetLowerRightNode();
					far_node[4] = (tmp) ? tmp->m_pNextTreeNodes[1] : NULL;
					far_node[5] = (tmp) ? tmp->m_pNextTreeNodes[2] : NULL;
					far_node[6] = (tmp) ? tmp->m_pNextTreeNodes[3] : NULL;
				}
				break;
			default:
				assert(0);
				break;
		}
	}


	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::GetDirectEvalNode(
		const FMCSM2DTreeNode<VAL_TYPE> *src_node,
		FMCSM2DTreeNode<VAL_TYPE> **dst_nodes)
	{
		dst_nodes[0] = src_node->GetUpperLeftNode();
		dst_nodes[1] = src_node->GetUpperNode();
		dst_nodes[2] = src_node->GetUpperRightNode();
		dst_nodes[3] = src_node->GetLeftNode();
		dst_nodes[4] = src_node->GetRightNode();
		dst_nodes[5] = src_node->GetLowerLeftNode();
		dst_nodes[6] = src_node->GetLowerNode();
		dst_nodes[7] = src_node->GetLowerRightNode();
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::PlotTreeNodeBox(
		const FMCSM2DTreeNode<VAL_TYPE> *node)
	{
		printf("%f %f \n",
			static_cast<double>(node->m_CenterPos[X_KEY] - node->m_BoxLength/2.0),
			static_cast<double>(node->m_CenterPos[Y_KEY] + node->m_BoxLength/2.0));
		printf("%f %f \n",
			static_cast<double>(node->m_CenterPos[X_KEY] + node->m_BoxLength/2.0),
			static_cast<double>(node->m_CenterPos[Y_KEY] + node->m_BoxLength/2.0));
		printf("%f %f \n",
			static_cast<double>(node->m_CenterPos[X_KEY] - node->m_BoxLength/2.0),
			static_cast<double>(node->m_CenterPos[Y_KEY] - node->m_BoxLength/2.0));
		printf("%f %f \n",
			static_cast<double>(node->m_CenterPos[X_KEY] + node->m_BoxLength/2.0),
			static_cast<double>(node->m_CenterPos[Y_KEY] - node->m_BoxLength/2.0));

		int i;
		for(i=0; i<4; i++){
			if(node->m_pNextTreeNodes[i]){
				PlotTreeNodeBox(node->m_pNextTreeNodes[i]);
			}
		}
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::PlotNodeSourcePoint(
		const FMCSM2DTreeNode<VAL_TYPE> *node)
	{
		int i;
		if(node->m_pSourcePos){
			for(i=0; i<node->m_pSourcePos->dim(); i++){
				printf("%f %f \n",
					static_cast<double>(node->m_pSourcePos->GetVal(i).val[X_KEY]),
					static_cast<double>(node->m_pSourcePos->GetVal(i).val[Y_KEY]));
			}
		}
		for(i=0; i<4; i++){
			if(node->m_pNextTreeNodes[i]){
				PlotNodeSourcePoint(node->m_pNextTreeNodes[i]);
			}
		}
	}

	template<class VAL_TYPE>
	void FMCSM2D<VAL_TYPE>::PlotNodeBundaryPoint(
		const FMCSM2DTreeNode<VAL_TYPE> *node)
	{
		int i;
		if(node->m_pBundaryPos){
			for(i=0; i<node->m_pBundaryPos->dim(); i++){
				printf("%f %f \n",
					static_cast<double>(node->m_pBundaryPos->GetVal(i).val[X_KEY]),
					static_cast<double>(node->m_pBundaryPos->GetVal(i).val[Y_KEY]));
			}
		}
		for(i=0; i<4; i++){
			if(node->m_pNextTreeNodes[i]){
				PlotNodeBundaryPoint(node->m_pNextTreeNodes[i]);
			}
		}
	}

};

#endif //_NMC_FMCSM2D_H_