#ifndef _NMC_UCU_FMCSM3D_TREE_H_
#define _NMC_UCU_FMCSM3D_TREE_H_


namespace nmc{

	__host__ __device__ static __inline__
	unsigned int GetCellIndex(unsigned int x, unsigned int y, unsigned int z, unsigned int level)
	{
		return (z << (level + 1)) + (y << level) + x;
	}

	__host__ __device__ static __inline__
	void GetCellPos(
		unsigned int cell_index,
		unsigned int level,
		unsigned int* x,
		unsigned int* y,
		unsigned int* z)
	{
		(*z) = cell_index >> (level + 1);
		(*y) = (cell_index - ((*z) << (level + 1))) >> level;
		(*x) = cell_index - ((*z) << (level + 1)) - ((*y) << level);
	}

	__host__ __device__ static __inline__
	unsigned int GetChildCellIndex(
		unsigned int x,
		unsigned int y,
		unsigned int z,
		unsigned int level,
		unsigned int sub_child_index)
	{
		return (x << 1) + (sub_child_index & 0x1) +
			(((y << 1) + ((sub_child_index & 0x2) >> 1)) << (level + 1)) +
			(((z << 1) + ((sub_child_index & 0x4) >> 2)) << (level + 2));
	}

	class ucuFMCSM3DTreeNode{
	public:
		ucuFMCSM3DTreeNode()
			: m_NodeLevel(-1),
			m_CellSize(0.0),
			m_MultipoleCoeff(NULL),
			m_LocalCoeff(NULL),
			m_MultipoleCoeff_dev(NULL),
			m_LocalCoeff_dev(NULL),
			m_pParentNode(NULL),
			m_pNextTreeNode(NULL){
		}
		~ucuFMCSM3DTreeNode(){
			if(m_MultipoleCoeff){
				delete m_MultipoleCoeff;
			}
			if(m_LocalCoeff){
				delete m_LocalCoeff;
			}
			if(m_MultipoleCoeff_dev){
				cudaFree(m_MultipoleCoeff_dev);
			}
			if(m_LocalCoeff_dev){
				cudaFree(m_LocalCoeff_dev);
			}
		}

		//–¢ŽÀ‘•
		static ucuFMCSM3DTreeNode* CreateNextNode(ucuFMCSM3DTreeNode* node, int term_num){
			node->m_pNextTreeNode = new ucuFMCSM3DTreeNode;

			return node->m_pNextTreeNode;
		}

		int m_NodeLevel;
		float m_CellSize;

		float2* m_MultipoleCoeff;
		float2* m_LocalCoeff;

		float2* m_MultipoleCoeff_dev;
		float2* m_LocalCoeff_dev;

		ucuFMCSM3DTreeNode* m_pParentNode;
		ucuFMCSM3DTreeNode* m_pNextTreeNode;
	};

	class ucuFMCSM3DTree{
	public:
		ucuFMCSM3DTree()
			: m_TermNum(-1),
			m_SourceNum(-1),
			m_MaxNodeLevel(-1),
			m_pStartIndexes(NULL),
			m_pEndIndexes(NULL),
			m_pStartIndexes_dev(NULL),
			m_pEndIndexes_dev(NULL),
			m_pPointsX(NULL),
			m_pPointsY(NULL),
			m_pPointsZ(NULL),
			m_pPointIndexes(NULL),
			m_pPointsX_dev(NULL),
			m_pPointsY_dev(NULL),
			m_pPointsZ_dev(NULL),
			m_pPointIndexes_dev(NULL),
			m_pCfVec_dev(NULL),
			m_pDst_dev(NULL),
			m_pRootTreeNode(NULL),
			m_pLeafTreeNode(NULL)
		{
		}

		~ucuFMCSM3DTree(){
			if(m_pStartIndexes){
				delete m_pStartIndexes;
			}
			if(m_pEndIndexes){
				delete m_pEndIndexes;
			}
			if(m_pPointsX){
				delete m_pPointsX;
			}
			if(m_pPointsY){
				delete m_pPointsY;
			}
			if(m_pPointsZ){
				delete m_pPointsZ;
			}
			if(m_pPointIndexes){
				delete m_pPointIndexes;
			}

			if(m_pStartIndexes_dev){
				cudaFree(m_pStartIndexes_dev);
			}
			if(m_pEndIndexes_dev){
				cudaFree(m_pEndIndexes_dev);
			}
			if(m_pPointsX_dev){
				cudaFree(m_pPointsX_dev);
			}
			if(m_pPointsY_dev){
				cudaFree(m_pPointsY_dev);
			}
			if(m_pPointsZ_dev){
				cudaFree(m_pPointsZ_dev);
			}
			if(m_pPointIndexes_dev){
				cudaFree(m_pPointIndexes_dev);
			}
			if(m_pCfVec_dev){
				cudaFree(m_pCfVec_dev);
			}
			if(m_pDst_dev){
				cudaFree(m_pDst_dev);
			}
		}

		static ucuFMCSM3DTree* CreateTree(
			int term_num,
			int max_level,
			const Vector< Vector3Dd > &source_points,
			const Vector< Vector3Dd > &collocation_points)
		{
			assert(max_level > 0);
			assert(source_points.dim() == collocation_points.dim());

			ucuFMCSM3DTree* tree = new ucuFMCSM3DTree;
			tree->m_TermNum = 2*term_num + 1;
			tree->m_MaxNodeLevel = max_level;
			tree->m_SourceNum = 2*source_points.dim();

			//allocate points coordinate
			int n = source_points.dim();
			tree->m_pPointsX = new float[ 2*n ];
			tree->m_pPointsY = new float[ 2*n ];
			tree->m_pPointsZ = new float[ 2*n ];
			tree->m_pPointIndexes = new int[ 2*n ];

			cutilSafeCall( cudaMalloc(&(tree->m_pPointsX_dev), 2*n*sizeof(float)) );
			cutilSafeCall( cudaMalloc(&(tree->m_pPointsY_dev), 2*n*sizeof(float)) );
			cutilSafeCall( cudaMalloc(&(tree->m_pPointIndexes_dev), 2*n*sizeof(int)) );

			cutilSafeCall( cudaMalloc(&(tree->m_pCfVec_dev), n * sizeof(float2)) );
			cutilSafeCall( cudaMalloc(&(tree->m_pDst_dev), n * sizeof(float2)) );

			//copy source and boundary points
			int i;
			for(i = 0 ; i < n; i++){
				tree->m_pPointsX[i] = source_points[i][0];
				tree->m_pPointsY[i] = source_points[i][1];
				tree->m_pPointsZ[i] = source_points[i][2];
				tree->m_pPointIndexes[i] = -i-1;
			}
			for(i = 0; i < n; i++){
				tree->m_pPointsX[n + i] = collocation_points[i][0];
				tree->m_pPointsY[n + i] = collocation_points[i][1];
				tree->m_pPointsZ[n + i] = collocation_points[i][2];
				tree->m_pPointIndexes[n + i] = i;
			}

			//create coordinate indexes array
			int num_of_max_level_cell = std::pow(8.0, max_level);
			tree->m_pStartIndexes = new int[num_of_max_level_cell];
			memset(tree->m_pStartIndexes, 0, sizeof(int) * num_of_max_level_cell);
			cutilSafeCall(
				cudaMalloc(
					&(tree->m_pStartIndexes_dev),
					num_of_max_level_cell * sizeof(int)) );
			
			tree->m_pEndIndexes = new int[num_of_max_level_cell];
			memset(tree->m_pEndIndexes, 0, sizeof(int) * num_of_max_level_cell);
			cutilSafeCall(
				cudaMalloc(
					&(tree->m_pEndIndexes_dev),
					num_of_max_level_cell * sizeof(int)) );

			return tree;
		}

		void CopyCoordToDevice(){
			int num_of_max_level_cell = std::pow(8.0, m_MaxNodeLevel);

			//copy cordinate to device memory
			cutilSafeCall(
				cudaMemcpy(
					m_pPointsX_dev,
					m_pPointsX,
					m_SourceNum * sizeof(float),
					cudaMemcpyHostToDevice) );
			cutilSafeCall(
				cudaMemcpy(
					m_pPointsY_dev,
					m_pPointsY,
					m_SourceNum * sizeof(float),
					cudaMemcpyHostToDevice) );
			cutilSafeCall(
				cudaMemcpy(
					m_pPointsZ_dev,
					m_pPointsZ,
					m_SourceNum * sizeof(float),
					cudaMemcpyHostToDevice) );
			cutilSafeCall(
				cudaMemcpy(
					m_pPointIndexes_dev,
					m_pPointIndexes,
					m_SourceNum * sizeof(int),
					cudaMemcpyHostToDevice) );
			cutilSafeCall(
				cudaMemcpy(
					m_pStartIndexes_dev,
					m_pStartIndexes,
					num_of_max_level_cell * sizeof(int),
					cudaMemcpyHostToDevice) );
			cutilSafeCall(
				cudaMemcpy(
					m_pEndIndexes_dev,
					m_pEndIndexes,
					num_of_max_level_cell * sizeof(int),
					cudaMemcpyHostToDevice) );
		}

		int m_TermNum;
		int m_SourceNum;
		int m_MaxNodeLevel;

		int* m_pStartIndexes;
		int* m_pEndIndexes;
		
		int* m_pStartIndexes_dev; //device memory
		int* m_pEndIndexes_dev; //

		float* m_pPointsX;
		float* m_pPointsY;
		float* m_pPointsZ;
		int* m_pPointIndexes;

		float* m_pPointsX_dev;
		float* m_pPointsY_dev;
		float* m_pPointsZ_dev;
		int* m_pPointIndexes_dev;

		float2* m_pCfVec_dev;
		float2* m_pDst_dev;

		ucuFMCSM3DTreeNode* m_pRootTreeNode;
		ucuFMCSM3DTreeNode* m_pLeafTreeNode;
	};
};

#endif //_NMC_UCU_FMCSM3D_TREE_H_