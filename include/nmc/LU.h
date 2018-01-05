#ifndef _NMC_LU_H_
#define _NMC_LU_H_

#include "nmc/Matrix.h"
#include "nmc/Vector.h"
#include "nmc/Complex.h"

#include <stdlib.h>
#include <math.h>

namespace nmc {

template <class VAL_TYPE>
class LU {
public:
	/*!
		LU分解
		A = LU
		@param[in] A
		@param[out] L
		@param[out] U
		@param[out] P インデックスベクトル
		@return LU分解が成功すればtrue、そうでなければfalseを返す
	*/
	static bool Decomp(const Matrix<VAL_TYPE>& A, Matrix<VAL_TYPE>& L, Matrix<VAL_TYPE>& U, Vectori& P);

	/*!
		LU分解
		@param[in,out] A
		@param[out] P
		@return
	*/
	static bool Decomp(Matrix<VAL_TYPE>& A, Vectori& P);

	/*!
		連立1次方程式LUx=bを解く
		@param[in] L 左下行列
		@param[in] U 右上行列
		@param[in] P インデックスベクトル
		@param[in] b 
		@param[out] x 解の出力先
	*/
	static void Solve(const Matrix<VAL_TYPE>& L, const Matrix<VAL_TYPE>& U, const Vectori& P, const Vector<VAL_TYPE>& b, Vector<VAL_TYPE>& x);

	/*!
		連立1次方程式Ax=bを解く
		@param[in] src
		@param[in] P
		@param[in] b
		@param[out] x
	*/
	static void Solve(const Matrix<VAL_TYPE>& src, const Vectori& P, const Vector<VAL_TYPE>& b, Vector<VAL_TYPE>& x);

	/*!
		行列の条件数を返す (名取-塚本)
		@param[in] L
		@param[in] U
		@param[in] P インデックスベクトル
		@param[in] mat_norm 行列のノルム
		@return 条件数を返す、1に近いほど行列の性質がよい。
	*/
	static VAL_TYPE Cond(const Matrix<VAL_TYPE>& L, const Matrix<VAL_TYPE>& U, const Vectori& P, VAL_TYPE mat_norm);

	/*!
		行列の条件数を返す (名取-塚本)
		@param[in] A
		@param[in] P インデックスベクトル
		@param[in] mat_norm 行列のノルム
		@return 条件数を返す、1に近いほど行列の性質がよい。
	*/
	static VAL_TYPE Cond(const Matrix<VAL_TYPE>& A, const Vectori& P, VAL_TYPE mat_norm);
};

template <class VAL_TYPE>
bool LU<VAL_TYPE>::Decomp(
	const Matrix<VAL_TYPE> &A,
	Matrix<VAL_TYPE> &L,
	Matrix<VAL_TYPE> &U,
	Vectori &P)
{
	assert(A.row() == A.col());
	assert(A.isEqualSize( L ));
	assert(A.isEqualSize( U ));
	assert(A.row() == P.dim());

	int i,j,k;
	//init
	A.CopyTo( U );

	for(i=0; i<P.dim(); i++){
		P[i] = i;
	}

	for(k=0; k<A.col()-1; k++){
		//pivoting
		int p = k;
		VAL_TYPE max_p = U.GetVal(P[k], k).Abs();
		for(i=k+1; i<A.row(); i++){
			VAL_TYPE val = U.GetVal(P[i], k).Abs();
			if(max_p < val){
				max_p = val;
				p = i;
			}
		}
		
		//swap pivot
		int tmp_p = P[p];
		P[p] = P[k];
		P[k] = tmp_p;

		//
		L.SetVal(P[k], k, static_cast<VAL_TYPE>(1.0));
		
		for(i=k+1; i<A.row(); i++){
			VAL_TYPE m = U.GetVal(P[i], k) / U.GetVal(P[k], k);
			L.SetVal(P[i], k, m);
			U.SetVal(P[i], k, static_cast<VAL_TYPE>(0.0));
			for(j=k+1; j<A.col(); j++){
				VAL_TYPE a = U.GetVal(P[i], j) - m * U.GetVal(P[k], j);
				U.SetVal(P[i], j, a);
			}
		}
	}

	L.SetVal(P[k], k, static_cast<VAL_TYPE>(1.0));
	return true;
}

template <class VAL_TYPE>
bool LU<VAL_TYPE>::Decomp(
	Matrix<VAL_TYPE> &A,
	Vectori &P)
{
	if(A.row() != A.col())
		return false;
	if(A.row() != P.dim())
		return false;

	int i,j,k;

	//init
	for(i=0; i<P.dim(); i++){
		P[i] = i;
	}

	for(k=0; k<A.col()-1; k++){
		//pivoting
		int p = k;
		VAL_TYPE max_p = A.GetVal(P[k], k).Abs();
		for(i=k+1; i<A.row(); i++){
			VAL_TYPE val = A.GetVal(P[i], k).Abs();
			if(max_p < val){
				max_p = val;
				p = i;
			}
		}
		//swap pivot
		int tmp_p = P[p];
		P[p] = P[k];
		P[k] = tmp_p;

		for(i=k+1; i<A.row(); i++){
			VAL_TYPE m = A.GetVal(P[i], k) / A.GetVal(P[k], k);
			A.SetVal(P[i], k, m);
			for(j=k+1; j<A.col(); j++){
				VAL_TYPE a = A.GetVal(P[i], j) - m * A.GetVal(P[k], j);
				A.SetVal(P[i], j, a);
			}
		}
	}

	return true;
}

template <class VAL_TYPE>
void LU<VAL_TYPE>::Solve(
	const Matrix<VAL_TYPE> &L,
	const Matrix<VAL_TYPE> &U,
	const Vectori &P,
	const Vector<VAL_TYPE> &b,
	Vector<VAL_TYPE> &x)
{
	int i,j;
	int P_i;
	//forward substitution
	Vector<VAL_TYPE> y(b.dim());
	for(i=0; i<y.dim(); i++){
		P_i = P[i];
		VAL_TYPE sum = static_cast<VAL_TYPE>(0.0);
		for(j=0; j<i; j++){
			sum += L.GetVal(P_i, j) * y[j];
		}
		y[i] = b[P_i] - sum;
	}

	//backward substitution
	for(i=b.dim()-1; i>=0; i--){
		P_i = P[i];
		VAL_TYPE sum = static_cast<VAL_TYPE>(0.0);
		for(j=i+1; j<b.dim(); j++){
			sum += U.GetVal(P_i, j) * x[j];
		}
		x[i] = (y[i] - sum)/U.GetVal(P_i, i);
	}
}

template <class VAL_TYPE>
void LU<VAL_TYPE>::Solve(
	const Matrix<VAL_TYPE> &src,
	const Vectori &P,
	const Vector<VAL_TYPE> &b,
	Vector<VAL_TYPE> &x)
{
	int i,j;
	int P_i;
	//forward substitution
	Vector<VAL_TYPE> y(b.dim());
	for(i=0; i<y.dim(); i++){
		P_i = P[i];
		VAL_TYPE sum=0.0;
		for(j=0; j<i; j++){
			sum += src.GetVal(P_i, j) * y[j];
		}
		y[i] = b[P_i] - sum;
	}

	//backward substitution
	for(i=b.dim()-1; i>=0; i--){
		P_i = P[i];
		VAL_TYPE sum=0.0;
		for(j=i+1; j<b.dim(); j++){
			sum += src.GetVal(P_i, j) * x[j];
		}
		x[i] = (y[i] - sum)/src.GetVal(P_i, i);
	}
}

template <class VAL_TYPE>
VAL_TYPE LU<VAL_TYPE>::Cond(
	const Matrix<VAL_TYPE> &L,
	const Matrix<VAL_TYPE> &U,
	const Vectori &P,
	VAL_TYPE mat_norm)
{
	int i,k;
	VAL_TYPE t = static_cast<VAL_TYPE>(0.0);
	VAL_TYPE max_y = static_cast<VAL_TYPE>(0.0);
	Vector<VAL_TYPE> e(P.dim());
	Vector<VAL_TYPE> y(P.dim());
	Vector<VAL_TYPE> v(P.dim());

	//forward substitution
	VAL_TYPE zero = static_cast<VAL_TYPE>(0.0);
	for(k=0; k<e.dim(); k++){
		t = zero;
		for(i=0; i<k; i++){
			t += U.GetVal(P[i], k) * v[i];
		}

		if(t < zero){
			e[k] = static_cast<VAL_TYPE>(1.0);
		}
		else{
			e[k] = static_cast<VAL_TYPE>(-1.0);
		}

		v[k] = (e[k] - t)/U.GetVal(P[k],k);
	}

	//backword substitution
	for(k=v.dim()-1; k>=0; k--){
		t = v[k];
		for(i=v.dim()-1; i >= k+1; i--){
			t -= L.GetVal(P[i], k) * y[ P[i] ];
		}

		y[ P[k] ] = t;

		if(max_y < t.Abs()){
			max_y = t.Abs();
		}
	}

	return max_y * mat_norm;
}

template <class VAL_TYPE> VAL_TYPE LU<VAL_TYPE>::Cond(
	const Matrix<VAL_TYPE> &A,
	const Vectori &P,
	VAL_TYPE mat_norm)
{
	int i,k;
	VAL_TYPE t = static_cast<VAL_TYPE>(0.0);
	VAL_TYPE max_y = static_cast<VAL_TYPE>(0.0);
	Vector<VAL_TYPE> e(P.dim());
	Vector<VAL_TYPE> y(P.dim());
	Vector<VAL_TYPE> v(P.dim());

	//forward substitution
	VAL_TYPE zero = static_cast<VAL_TYPE>(0.0);
	for(k=0; k<e.dim(); k++){
		t = zero;
		for(i=0; i<k; i++){
			t += A.GetVal(P[i], k) * v[i];
		}

		if(t < zero){
			e[k] = static_cast<VAL_TYPE>(1.0);
		}
		else{
			e[k] = static_cast<VAL_TYPE>(-1.0);
		}

		v[k] = (e[k] - t)/A.GetVal(P[k],k);
	}

	//backword substitution
	for(k=v.dim()-1; k>=0; k--){
		t = v[k];
		for(i=v.dim()-1; i >= k+1; i--){
			t -= A.GetVal(P[i], k) * y[ P[i] ];
		}

		y[ P[k] ] = t;

		if(max_y < t.Abs()){
			max_y = t.Abs();
		}
	}

	return max_y * mat_norm;
}

typedef LU< Double > LUd;
typedef LU< Complexd > LUcmd;

};

#endif //_NMC_LU_H_