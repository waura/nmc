#ifndef _NMC_MATRIX_H_
#define _NMC_MATRIX_H_

#include <stdlib.h>
#include <assert.h>

#include "nmc/Complex.h"
#include "nmc/Vector.h"

namespace nmc
{

template <class VAL_TYPE>
class Matrix {
public:
	Matrix(int row, int col);
	Matrix();
	Matrix(const Matrix &mat);
	~Matrix();

	void Init(int row, int col);

	inline int col() const;
	inline int row() const;

	inline VAL_TYPE l1_norm() const;
	
	inline bool isEqualSize(const Matrix& mat) const;
	
	inline VAL_TYPE GetVal(int row, int col) const;
	inline void SetVal(int row, int col, const VAL_TYPE& val);
	
	/*!
		正方行列の場合、単位行列にする
	*/
	bool Identity();

	void CopyTo(Matrix& mat) const;

	inline void Multiply(const Matrix& mat);

	inline Matrix& operator=(const Matrix& mat);

	friend Vector<VAL_TYPE> operator*(Matrix mat, Vector<VAL_TYPE> vec){
		assert(mat.col() == vec.dim());
		Vector<VAL_TYPE> ret;
		ret.Init(mat.row());

		int i,j;
		for(i=0; i < mat.row(); i++){
			VAL_TYPE sum = static_cast<VAL_TYPE>(0);
			for(j=0; j < mat.col(); j++){
				sum += mat.GetVal(i, j) * vec[j];
			}
			ret[i] = sum;
		}
		return ret;
	}

	friend Matrix<VAL_TYPE> operator*(Matrix mat1, Matrix mat2){
		Matrix<VAL_TYPE> ret;
		ret.Init(mat1.row(), mat2.col());
		int i,j,k;
		for(i=0; i<mat1.row(); i++){
			for(j=0; j<mat2.col(); j++){
				VAL_TYPE sum=0.0;
				for(k=0; k<mat1.col(); k++){
					sum += mat1.GetVal(i, k) * mat2.GetVal(k, j);
				}
				ret.SetVal(i, j, sum);
			}
		}
		return ret;
	}

	void dump() const;
private:
	int _col;
	int _row;
	VAL_TYPE* _d;
};

template <class VAL_TYPE> Matrix<VAL_TYPE>::Matrix(int row, int col)
{
	_d = NULL;
	Init(row, col);
}

template <class VAL_TYPE> Matrix<VAL_TYPE>::Matrix()
{
	_col = 0;
	_row = 0;
	_d = NULL;
}

template <class VAL_TYPE> Matrix<VAL_TYPE>::Matrix(const Matrix<VAL_TYPE>& mat)
{
	_d = NULL;
	Init(mat.row(), mat.col());

	int i,j;
	for(j=0; j<mat.col(); j++){
		for(i=0; i<mat.row(); i++){
			SetVal(i,j, mat.GetVal(i, j));
		}
	}
}

template <class VAL_TYPE> Matrix<VAL_TYPE>::~Matrix()
{
	if(_d){
		delete[] _d;
	}
}

template <class VAL_TYPE> void Matrix<VAL_TYPE>::Init(int row, int col)
{
	if(_d){
		if((row != _row) || (col != _col)){
			delete[] _d;
		}
	}

	_col = col;
	_row = row;
	_d = new VAL_TYPE[col*row];

	//zero clear
	if(_d){
		int i;
		for(i=0; i<col*row; i++){
			_d[i] = static_cast<VAL_TYPE>(0.0);
		}
	}
}

template <class VAL_TYPE> inline int Matrix<VAL_TYPE>::col() const
{
	return _col;
}

template <class VAL_TYPE> inline int Matrix<VAL_TYPE>::row() const
{
	return _row;
}

template <class VAL_TYPE> inline VAL_TYPE Matrix<VAL_TYPE>::l1_norm() const
{
	int i,j;
	VAL_TYPE norm = static_cast<VAL_TYPE>(0.0);
	VAL_TYPE s,t;

	for(j=0; j<_col; j++){
		t = static_cast<VAL_TYPE>(0.0);
		for(i=0; i<_row; i++){
			s = GetVal(i,j);
			t += s.Abs();
		}

		if( t > norm ) norm = t;
	}

	return norm;
}

template <class VAL_TYPE> inline bool Matrix<VAL_TYPE>::isEqualSize(const Matrix& mat) const
{
	if(_col == mat.col() && _row == mat.row()){
		return true;
	}
	return false;
}

template <class VAL_TYPE> inline VAL_TYPE Matrix<VAL_TYPE>::GetVal(int row, int col) const
{
	assert((0 <= row && row <= _row-1) && (0 <= col && col <= _col-1));
	return _d[row + col * _col];
}

template <class VAL_TYPE> inline void Matrix<VAL_TYPE>::SetVal(int row, int col, const VAL_TYPE& val)
{
	assert((0 <= row && row <= _row-1) && (0 <= col && col <= _col-1));
	_d[row + col * _col] = val;
}

template <class VAL_TYPE> bool Matrix<VAL_TYPE>::Identity()
{
	if(_row != _col)
		return false;

	int i,j;
	for(i=0; i<_row; i++){
		for(j=0; j<_col; j++){
			if(i == j)
				SetVal(i, j, 1);
			else
				SetVal(i, j, 0);
		}
	}
	return true;
}

template <class VAL_TYPE> void Matrix<VAL_TYPE>::CopyTo(Matrix &mat) const
{
	mat.Init(this->row(), this->col());

	int i,j;
	for(j=0; j<this->col(); j++){
		for(i=0; i<this->row(); i++){
			mat.SetVal(i,j, this->GetVal(i, j));
		}
	}
}

template <class VAL_TYPE> inline void Matrix<VAL_TYPE>::Multiply(const Matrix& mat)
{
	Matrix tmp(this->row(), mat.col());

	int i,j,k;
	for(i=0; i<tmp.row(); i++){
		for(j=0; j<tmp.col(); j++){
			VAL_TYPE sum=0.0;
			for(k=0; k<this->col(); k++){
				sum += this->GetVal(i, k) * mat.GetVal(k, j);
			}
			tmp.SetVal(i, j, sum);
		}
	}

	tmp.CopyTo( *this );
}

template <class VAL_TYPE> inline Matrix<VAL_TYPE>& Matrix<VAL_TYPE>::operator =(const Matrix<VAL_TYPE>& mat)
{
	Init(mat.row(), mat.col());

	int i,j;
	for(j=0; j<mat.col(); j++){
		for(i=0; i<mat.row(); i++){
			SetVal(i,j, mat.GetVal(i, j));
		}
	}
	return *this;
}

template <class VAL_TYPE> void Matrix<VAL_TYPE>::dump() const
{
	int i,j;
	for(i=0; i<row(); i++){
		printf("| ");
		for(j=0; j<col(); j++){
			VAL_TYPE val = GetVal(i, j);
			val.dump();
			printf(" ,");
		}
		printf("|\n");
	}
}

typedef Matrix< Double > Matrixd;
typedef Matrix< Float > Matrixf;
typedef Matrix< Int > Matrixi;
typedef Matrix< Complexd > Matrixcmd;

};

#endif //_NMC_MATRIX_H_ 