#ifndef _FMCSM2D_HELMHOLTZ_KERNEL_H_
#define _FMCSM2D_HELMHOLTZ_KERNEL_H_

#include "nmc/def.h"
#include "nmc/FMCSM2DKernel.h"

#include <gsl/gsl_sf_bessel.h>

namespace nmc
{

class FMCSM2D_Helmholtz_Kernel : public FMCSM2DKernel<Complexd>
{
public:
	void SetWaveNumber(double k){
		m_WaveNumber = k;
	}

	void Init(FMCSM2DTree<Complexd>* tree);
	inline void ElementMPCefficients(FMCSM2DTreeNode<Complexd>* node,
		const Vector<Complexd>& c_vec);
	inline void M2M_Translation(FMCSM2DTreeNode<Complexd>* node);
	inline void L2L_Translation(FMCSM2DTreeNode<Complexd>* node);
	inline void M2L_Translation(const FMCSM2DTreeNode<Complexd>* src,
		FMCSM2DTreeNode<Complexd>* dst);
	inline void EvalMultipole(const FMCSM2DTreeNode<Complexd>* node,
		Vector<Complexd>& dst);
	inline void EvalDirect(const FMMPos& from_pos, const FMMPos& to_pos,
		const Vector<Complexd>& cf, Vector<Complexd>& dst);
private:
	double m_WaveNumber;
	Vector< Vectord > m_JnTable;
	Vector< Vectorcmd > m_ExpTable;
};

void FMCSM2D_Helmholtz_Kernel::Init(
	FMCSM2DTree<Complexd>* tree)
{
	int i,j;

	//create bessel function table
	int p = (tree->m_TermNum-1)/2;
	double rho;
	double r = (tree->m_pRootTreeNode->m_BoxLength)/2.0;
	m_JnTable.Init( tree->m_MaxNodeLevel );
	for(i=0; i<m_JnTable.dim(); i++){
		r /= 2.0;
		rho = sqrt(2.0*r*r);
		m_JnTable[i].Init(2*2*p+1);

		for(j=-2*p; j<=2*p; j++){
			m_JnTable[i].SetVal(j + 2*p, gsl_sf_bessel_Jn(j, rho));
		}
	}

	//create exp table
	double theta;
	Complexd exp;
	m_ExpTable.Init( 4 );
	for(i=0; i<4; i++){
		m_ExpTable[i].Init( 2*7 + 1 );

		theta = (2*i+1)*NMC_PI/4.0;
		for(j=-7; j<=7; j++){
			exp.SetNum(cos( j * theta), sin( j * theta) );
			m_ExpTable[i].SetVal(j+7, exp);
		}
	}
}

inline void FMCSM2D_Helmholtz_Kernel::ElementMPCefficients(
	FMCSM2DTreeNode<Complexd> *node,
	const Vector<Complexd>& c_vec)
{
	assert(node->m_pSourcePos);

	int s,t;
	double x,y;
	double phi;
	double rho;
	double J_s;
	Complexd c;
	Complexd exp;

	int p = (node->m_Multipole.dim()-1)/2;
	for(t=0; t<node->m_pSourcePos->dim(); t++){
		x = node->m_pSourcePos->GetVal(t).val[X_KEY] - node->m_CenterPos[X_KEY];
		y = node->m_pSourcePos->GetVal(t).val[Y_KEY] - node->m_CenterPos[Y_KEY];


		rho = sqrt(x*x + y*y);
		phi = atan2(y, x);
		c = c_vec[ node->m_pSourcePos->GetVal(t).index ];
		for(s=-p; s<=p; s++){
			//MP += c * J_s(k*rho) * exp(-i s phi)
			J_s = gsl_sf_bessel_Jn(s, m_WaveNumber*rho);
			exp.SetNum(cos(-s*phi), sin(-s*phi));			
			node->m_Multipole[s+p] += c * J_s * exp;
		}
	}
}

inline void FMCSM2D_Helmholtz_Kernel::M2M_Translation(
	FMCSM2DTreeNode<Complexd> *node)
{
	int s,t;
	int child_index;
	double J_st;
	Complexd exp;

	//double r = node->m_BoxLength/4.0;
	//double rho = sqrt(2.0*r*r);
	int st_idx;
	int exp_idx[4] = {1, 0, 2, 3}; //{3.0*PI/4, PI/4, 5.0*PI/4, 7.0*PI/4}
	//double phi[4] = {3.0*NMC_PI/4.0, NMC_PI/4.0, 5.0*NMC_PI/4.0, 7.0*NMC_PI/4.0};

	int p = (node->m_Multipole.dim()-1)/2;
	for(s=-p; s<=p; s++){
		for(t=-p; t<=p; t++){
			st_idx = (-(s-t)) % 8 + 7;
			//J_st = gsl_sf_bessel_Jn(s-t, m_WaveNumber*rho);
			J_st = m_JnTable[ node->m_NodeLevel ][ (s-t) + 2*p ];
			for(child_index=0; child_index<4; child_index++){
				if(node->m_pNextTreeNodes[child_index]){
					//new M_{s} += M_t * J_{s-t}(k rho) * exp(-i (s-t) phi)
					//exp.SetNum(cos(-(s-t)*phi[child_index]), sin(-(s-t)*phi[child_index]));
					exp = m_ExpTable[ exp_idx[child_index] ][ st_idx ];
					node->m_Multipole[s+p] += node->m_pNextTreeNodes[child_index]->m_Multipole[t+p] * J_st * exp;
				}
			}
		}
	}
}

inline void FMCSM2D_Helmholtz_Kernel::L2L_Translation(
	FMCSM2DTreeNode<Complexd> *node)
{
	int s,t;
	int child_index;
	double c;
	double J_ts;
	Complexd exp;

	//double r = node->m_BoxLength/4.0;
	//double rho = sqrt(2.0*r*r);
	int ts_idx;
	int exp_idx[4] = {3, 2, 0, 1}; //{7.0*NMC_PI/4.0, 5.0*NMC_PI/4.0, NMC_PI/4.0, 3.0*NMC_PI/4.0};
	//double phi[4] = {7.0*NMC_PI/4.0, 5.0*NMC_PI/4.0, NMC_PI/4.0, 3.0*NMC_PI/4.0};

	int p = (node->m_Local.dim()-1)/2;
	for(s=-p; s<=p; s++){
		for(t=-p; t<=p; t++){
			c = (((t-s) % 2) == 0) ? 1.0 : -1.0;
			ts_idx = (-(t-s)) % 8 + 7;
			//J_ts = gsl_sf_bessel_Jn(t-s, m_WaveNumber*rho);
			J_ts = m_JnTable[ node->m_NodeLevel ][ (t-s) + 2*p ];
			for(child_index=0; child_index<4; child_index++){
				if(node->m_pNextTreeNodes[child_index] &&
					node->m_pNextTreeNodes[child_index]->m_isBundaryPosInNode){
					// new L_{s} += (-1)^{t-s} * L_{t} * J_{t-s}(k rho) * exp(-i (t-s) phi)
					//exp.SetNum(cos(-(t-s)*phi[child_index]), sin(-(t-s)*phi[child_index]));
					exp = m_ExpTable[ exp_idx[child_index] ][ ts_idx ];
					node->m_pNextTreeNodes[child_index]->m_Local[s+p] += c * node->m_Local[t+p] * J_ts * exp;
				}
			}
		}
	}
}

inline void FMCSM2D_Helmholtz_Kernel::M2L_Translation(
	const FMCSM2DTreeNode<Complexd> *src,
	FMCSM2DTreeNode<Complexd> *dst)
{
	int s,t;
	double c;
	Complexd exp;
	Complexd H_st;

	
	double x = src->m_CenterPos[X_KEY] - dst->m_CenterPos[X_KEY]; //‹t‚¶‚á‚Ë?
	double y = src->m_CenterPos[Y_KEY] - dst->m_CenterPos[Y_KEY];
	double rho = sqrt(x*x + y*y);
	double phi = atan2(y, x);

	int p = (dst->m_Local.dim()-1)/2;
	for(s=-p; s<=p; s++){
		for(t=-p; t<=p; t++){
			// L_{s} += (-1)^{t} * M_{t} * H_{s+t}(k rho) * exp(i(s+t) phi)
			c = ((t % 2) == 0) ? 1.0 : -1.0;
			exp.SetNum(cos((s+t) * phi), sin((s+t) * phi));
			H_st.SetNum(gsl_sf_bessel_Jn(s+t, m_WaveNumber*rho), gsl_sf_bessel_Yn(s+t, m_WaveNumber*rho));
			dst->m_Local[s+p] += c * src->m_Multipole[t+p] * H_st * exp;		
		}
	}
}

inline void FMCSM2D_Helmholtz_Kernel::EvalMultipole(
	const FMCSM2DTreeNode<Complexd> *node,
	Vector<Complexd> &dst)
{
	assert(node->m_pBundaryPos);

	int s,t;
	int index;
	double x,y;
	double r,phi;
	Complexd exp;

	int p = (node->m_Local.dim()-1)/2;
	for(t=0; t<node->m_pBundaryPos->dim(); t++){
		index = node->m_pBundaryPos->GetVal(t).index;
		x = node->m_pBundaryPos->GetVal(t).val[X_KEY] - node->m_CenterPos[X_KEY];
		y = node->m_pBundaryPos->GetVal(t).val[Y_KEY] - node->m_CenterPos[Y_KEY];

		r = sqrt(x*x + y*y);
		phi = atan2(y, x);
		for(s=-p; s<=p; s++){
			//L_n * J_n(k r) * exp(-i s phi)
			exp.SetNum(cos(-s*phi), sin(-s*phi));
			dst[index] += node->m_Local[s+p] * gsl_sf_bessel_Jn(s, m_WaveNumber*r) * exp;
		}
	}
}

inline void FMCSM2D_Helmholtz_Kernel::EvalDirect(
	const FMMPos& from_pos,
	const FMMPos& to_pos,
	const Vector<Complexd>& cf,
	Vector<Complexd> &dst)
{
	double x,y,r;
	Complexd Hs;

	x = to_pos.val[X_KEY] - from_pos.val[X_KEY];
	y = to_pos.val[Y_KEY] - from_pos.val[Y_KEY];
	r = sqrt(x*x + y*y);
	Hs.SetNum(gsl_sf_bessel_J0(m_WaveNumber*r), gsl_sf_bessel_Y0(m_WaveNumber*r));
	dst[to_pos.index] += cf[from_pos.index] * Hs;
}

};

#endif //_FMCSM2D_HELMHOLTZ_KERNEL_H_