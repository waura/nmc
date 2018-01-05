#ifndef _FMCSM2D_TREE_H_
#define _FMCSM2D_TREE_H_

#include "nmc/Vector2D.h"

namespace nmc
{
	const int X_KEY = 0;
	const int Y_KEY = 1;

	struct FMMPos{
		double val[2];
		int index;
	};

	template <class VAL_TYPE>
	class FMCSM2DTreeNode{
	public:
		int m_NodeIndex; //0~3
		int m_NodeLevel;
		bool m_isBundaryPosInNode; //
		bool m_isSourcePosInNode;
		FMCSM2DTreeNode* m_pParentNode;
		FMCSM2DTreeNode* m_pNextTreeNodes[4];
		double m_BoxLength;
		Vector2Dd m_CenterPos;
		Vector<FMMPos>* m_pSourcePos;
		Vector<FMMPos>* m_pBundaryPos;
		Vector<VAL_TYPE> m_Multipole;
		Vector<VAL_TYPE> m_Local;

		FMCSM2DTreeNode(){
			int i;
			for(i=0; i<4; i++)
				m_pNextTreeNodes[i] = NULL;

			m_pSourcePos = NULL;
			m_pBundaryPos = NULL;
		}
		void Init(
			int node_index,
			int node_level,
			int term_num,
			double box_length,
			double center_x,
			double center_y,
			FMCSM2DTreeNode* parent){
			this->m_NodeIndex = node_index;
			this->m_NodeLevel = node_level;
			this->m_Multipole.Init(term_num+1);
			this->m_Local.Init(term_num+1);
			this->m_BoxLength = box_length;
			this->m_CenterPos[X_KEY] = center_x;
			this->m_CenterPos[Y_KEY] = center_y;
			this->m_pParentNode = parent;
		}

		inline FMCSM2DTreeNode* GetParentNode() const{
			return m_pParentNode;
		}
		inline FMCSM2DTreeNode* GetUpperLeftNode() const{
			switch(m_NodeIndex){
				case 0:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetUpperLeftNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[3];
						}
					}
					break;
				case 1:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetUpperNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[2];
						}
					}
					break;
				case 2:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetLeftNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[1];
						}
					}
					break;
				case 3:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						if(parent){
							return parent->m_pNextTreeNodes[0];
						}
					}
					break;
				default:
					assert(0);
					break;
			}
			return NULL;
		}
		inline FMCSM2DTreeNode* GetUpperNode() const{
			switch(m_NodeIndex){
				case 0:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetUpperNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[2];
						}
					}
					break;
				case 1:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetUpperNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[3];
						}
					}
					break;
				case 2:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						if(parent){
							return parent->m_pNextTreeNodes[0];
						}
					}
					break;
				case 3:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						if(parent){
							return parent->m_pNextTreeNodes[1];
						}
					}
					break;
				default:
					assert(0);
					break;
			}
			return NULL;
		}
		inline FMCSM2DTreeNode* GetUpperRightNode() const{
			switch(m_NodeIndex){
				case 0:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetUpperNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[3];
						}
					}
					break;
				case 1:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetUpperRightNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[2];
						}
					}
					break;
				case 2:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						if(parent){
							return parent->m_pNextTreeNodes[1];
						}
					}
					break;
				case 3:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetRightNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[0];
						}
					}
					break;
				default:
					assert(0);
					break;
			}
			return NULL;
		}
		inline FMCSM2DTreeNode* GetLeftNode() const{
			switch(m_NodeIndex){
				case 0:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetLeftNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[1];
						}
					}
					break;
				case 1:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						if(parent){
							return parent->m_pNextTreeNodes[0];
						}
					}
					break;
				case 2:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetLeftNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[3];
						}
					}
					break;
				case 3:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						if(parent){
							return parent->m_pNextTreeNodes[2];
						}
					}
					break;
				default:
					assert(0);
					break;
			}
			return NULL;
		}
		inline FMCSM2DTreeNode* GetRightNode() const{
			switch(m_NodeIndex){
				case 0:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						if(parent){
							return parent->m_pNextTreeNodes[1];
						}
						break;
					}
				case 1:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetRightNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[0];
						}
					}
					break;
				case 2:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						if(parent){
						return parent->m_pNextTreeNodes[3];
						}
						break;
					}
				case 3:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetRightNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[2];
						}
					}
					break;
				default:
					assert(0);
					break;
			}
			return NULL;
		}
		inline FMCSM2DTreeNode* GetLowerLeftNode() const{
			switch(m_NodeIndex){
				case 0:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetLeftNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[3];
						}
					}
					break;
				case 1:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						if(parent){
							return parent->m_pNextTreeNodes[2];
						}
					}
					break;
				case 2:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetLowerLeftNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[1];
						}
					}
					break;
				case 3:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetLowerNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[0];
						}
					}
					break;
				default:
					assert(0);
					break;
			}
			return NULL;
		}
		inline FMCSM2DTreeNode* GetLowerNode() const{
			switch(m_NodeIndex){
				case 0:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						if(parent){
							return parent->m_pNextTreeNodes[2];
						}
					}
					break;
				case 1:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						if(parent){
							return parent->m_pNextTreeNodes[3];
						}
					}
					break;
				case 2:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetLowerNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[0];
						}
					}
					break;
				case 3:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetLowerNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[1];
						}
					}
					break;
				default:
					assert(0);
					break;
			}
			return NULL;
		}
		inline FMCSM2DTreeNode* GetLowerRightNode() const{
			switch(m_NodeIndex){
				case 0:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						if(parent){
							return parent->m_pNextTreeNodes[3];
						}
					}
					break;
				case 1:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetRightNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[2];
						}
					}
					break;
				case 2:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetLowerNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[1];
						}
					}
					break;
				case 3:
					{
						FMCSM2DTreeNode* parent = GetParentNode();
						parent = (parent) ? parent->GetLowerRightNode() : NULL;
						if(parent){
							return parent->m_pNextTreeNodes[0];
						}
					}
					break;
				default:
					assert(0);
					break;
			}
			return NULL;
		}
	};

	template <class VAL_TYPE>
	class FMCSM2DTree{
	public:
		int m_TermNum;
		int m_SourceNum;
		int m_MaxNodeLevel;
		FMCSM2DTreeNode<VAL_TYPE>* m_pRootTreeNode;
	};

};

#endif //_FMCSM2D_TREE_H_