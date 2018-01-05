#ifndef _NMC_UCU_FMCSM2D_TREE_H_
#define _NMC_UCU_FMCSM2D_TREE_H_


namespace nmc{

	__host__ __device__ static __inline__
	unsigned int BitSeparate32(unsigned int n)
	{
		n = (n | (n << 8)) & 0x00ff00ff;
		n = (n | (n << 4)) & 0x0f0f0f0f;
		n = (n | (n << 2)) & 0x33333333;
		n = (n | (n << 1)) & 0x55555555;
		return n;
	}

	__host__ __device__ static __inline__
	unsigned int BitConnect32(unsigned int n)
	{
		n = ((n >> 1) & 0xaaaaaaaa) | (n & 0x55555555);
		n = ((n >> 2) & 0xcccccccc) | (n & 0x33333333);
		n = ((n >> 4) & 0xf0f0f0f0) | (n & 0x0f0f0f0f);
		n = ((n >> 8) & 0xff00ff00) | (n & 0x00ff00ff);
		return n;
	}

	__host__ __device__ static __inline__
	unsigned int GetMortonCellIndex(unsigned short x, unsigned short y)
	{
		return (BitSeparate32(x) | (BitSeparate32(y) << 1));
	}

	__host__ __device__ static __inline__
	void GetMortonCellPos(unsigned int cell_index, unsigned short* x, unsigned short* y)
	{
		(*x) = BitConnect32(cell_index);
		(*y) = BitConnect32(cell_index >> 1);
	}

	__host__ __device__ static __inline__
	unsigned int GetMortonChildCellIndex(unsigned int cell_index, int sub_child_index)
	{
		return (cell_index << 2) + sub_child_index;
	}

	__host__ __device__ static __inline__
	unsigned int GetMortonParentCellIndex(unsigned int cell_index)
	{
		return (cell_index >> 2);
	}

	__host__ __device__ static __inline__
	unsigned int GetCellIndex(unsigned int x, unsigned int y, unsigned int level)
	{
		return (y << level) + x;
	}

	__host__ __device__ static __inline__
	void GetCellPos(
		unsigned int cell_index,
		unsigned int level,
		unsigned int* x,
		unsigned int* y)
	{
		(*y) = cell_index >> level;
		(*x) = cell_index - ((*y) << level);
	}

	__host__ __device__ static __inline__
	unsigned int GetChildCellIndex(
		unsigned int x,
		unsigned int y,
		unsigned int level,
		unsigned int sub_child_index)
	{
		return ((x << 1) + (sub_child_index & 0x1)) + 
			(((y << 1) + ((sub_child_index & 0x2) >> 1)) << (level + 1));
	}

	__host__ __device__ static __inline__
	unsigned int GetRealCellIndex(unsigned int x, unsigned int y, unsigned int level)
	{
		return 0;
	}

	class ucuFMCSM2DTreeNode{
	public:
		ucuFMCSM2DTreeNode()
			: m_NodeLevel(-1),
			m_CellSize(0.0),
			m_MultipoleCoeff(NULL),
			m_LocalCoeff(NULL),
			m_MultipoleCoeff_dev(NULL),
			m_LocalCoeff_dev(NULL),
			m_pParentNode(NULL),
			m_pNextTreeNode(NULL){
		}
		~ucuFMCSM2DTreeNode(){
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

		//
		static ucuFMCSM2DTreeNode* CreateNextNode(ucuFMCSM2DTreeNode* node, int term_num){
			node->m_pNextTreeNode = new ucuFMCSM2DTreeNode;
			node->m_pNextTreeNode->m_NodeLevel = node->m_NodeLevel + 1;
			node->m_pNextTreeNode->m_CellSize = node->m_CellSize/2.0;
			node->m_pNextTreeNode->m_pParentNode = node;

			int num_of_cell = std::pow(4.0, node->m_pNextTreeNode->m_NodeLevel);
			size_t num_of_coeff = num_of_cell * (term_num + 1);
			node->m_pNextTreeNode->m_MultipoleCoeff = new float2[num_of_coeff];
			node->m_pNextTreeNode->m_LocalCoeff = new float2[num_of_coeff];

			cutilSafeCall(
				cudaMalloc(
					&(node->m_pNextTreeNode->m_MultipoleCoeff_dev),
					num_of_coeff * sizeof(float2)) );
			cutilSafeCall(
				cudaMalloc(
					&(node->m_pNextTreeNode->m_LocalCoeff_dev),
					num_of_coeff * sizeof(float2)) );
			return node->m_pNextTreeNode;
		}

		int m_NodeLevel;
		float m_CellSize;

		float2* m_MultipoleCoeff;
		float2* m_LocalCoeff;

		float2* m_MultipoleCoeff_dev;
		float2* m_LocalCoeff_dev;

		ucuFMCSM2DTreeNode* m_pParentNode;
		ucuFMCSM2DTreeNode* m_pNextTreeNode;
	};

	class ucuFMCSM2DTree{
	public:
		ucuFMCSM2DTree()
			: m_TermNum(-1),
			m_SourceNum(-1),
			m_MaxNodeLevel(-1),
			m_pStartIndexes(NULL),
			m_pEndIndexes(NULL),
			m_pStartIndexes_dev(NULL),
			m_pEndIndexes_dev(NULL),
			m_pPointsX(NULL),
			m_pPointsY(NULL),
			m_pPointIndexes(NULL),
			m_pPointsX_dev(NULL),
			m_pPointsY_dev(NULL),
			m_pPointIndexes_dev(NULL),
			m_pCfVec_dev(NULL),
			m_pDst_dev(NULL),
			m_pRootTreeNode(NULL),
			m_pLeafTreeNode(NULL)
		{
		}

		~ucuFMCSM2DTree(){
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
			if(m_pPointIndexes_dev){
				cudaFree(m_pPointIndexes_dev);
			}
			if(m_pCfVec_dev){
				cudaFree(m_pCfVec_dev);
			}
			if(m_pDst_dev){
				cudaFree(m_pCfVec_dev);
			}
		}

		static ucuFMCSM2DTree* CreateTree(
			int term_num,
			int max_level,
			const Vector< Vector2Dd > &source_points,
			const Vector< Vector2Dd > &collocation_points)
		{
			assert(max_level > 0);
			assert(source_points.dim() == collocation_points.dim());

			ucuFMCSM2DTree* tree = new ucuFMCSM2DTree;
			tree->m_TermNum = 2*term_num + 1;
			tree->m_MaxNodeLevel = max_level;
			tree->m_SourceNum = 2*source_points.dim();

			//allocate points coordinate
			int n = source_points.dim();
			tree->m_pPointsX = new float[ 2*n ];
			tree->m_pPointsY = new float[ 2*n ];
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
				tree->m_pPointIndexes[i] = -i-1;
			}
			for(i = 0; i < n; i++){
				tree->m_pPointsX[n + i] = collocation_points[i][0];
				tree->m_pPointsY[n + i] = collocation_points[i][1];
				tree->m_pPointIndexes[n + i] = i;
			}

			//create coordinate indexes array
			int num_of_max_level_cell = std::pow(4.0, max_level);
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
			int num_of_max_level_cell = std::pow(4.0, m_MaxNodeLevel);

			//copy coordinate to device memory
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
		int* m_pPointIndexes;

		float* m_pPointsX_dev; //
		float* m_pPointsY_dev; //
		int* m_pPointIndexes_dev; //

		float2* m_pCfVec_dev;
		float2* m_pDst_dev;

		ucuFMCSM2DTreeNode* m_pRootTreeNode;
		ucuFMCSM2DTreeNode* m_pLeafTreeNode;
	};
};

#endif //_NMC_UCU_FMCSM2D_TREEE_H_