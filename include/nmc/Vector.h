#ifndef _VECTER_H_
#define _VECTER_H_

#include "Complex.h"

#include <assert.h>

namespace nmc
{

	template <class VAL_TYPE>
	class Vector{
	public:
		Vector(int dimensional);
		Vector();
		Vector(const Vector &vec);
		virtual ~Vector();

		void Init(int dimensional);

		inline int dim() const;
		inline void SetVal(int i, const VAL_TYPE& val);
		inline VAL_TYPE& GetVal(int i) const;
		inline void GetVal(int i, VAL_TYPE& out) const;
		inline void Swap(int i, int j);
		void CopyTo(Vector& vec) const;

		inline Vector& operator=(const Vector& vec);
		inline Vector& operator+=(const Vector& vec);
		inline Vector& operator-=(const Vector& vec);

		inline VAL_TYPE& operator[](int index) const;

		inline friend Vector operator+(const Vector& vec1, const Vector& vec2){
			assert(vec1.dim() == vec2.dim());
			
			int dim = vec1.dim();
			Vector vec( dim );

			int i;
			VAL_TYPE* ptr = vec._d;
			VAL_TYPE* pvec1 = vec1._d;
			VAL_TYPE* pvec2 = vec2._d;
			for(i=0; i<dim; i++){
				(*ptr++) = (*pvec1++) + (*pvec2++);
			}
			return vec;
		}
		inline friend Vector operator-(const Vector& vec1, const Vector& vec2){
			assert(vec1.dim() == vec2.dim());
			
			int dim = vec1.dim();
			Vector vec( dim );

			int i;
			VAL_TYPE* ptr = vec._d;
			VAL_TYPE* pvec1 = vec1._d;
			VAL_TYPE* pvec2 = vec2._d;
			for(i=0; i<dim; i++){
				(*ptr++) = (*pvec1++) - (*pvec2++);
			}
			return vec;
		}

		friend std::ostream& operator<<(std::ostream& os, const Vector& vec){
			int i;
			os << "|" << std::endl;
			for(i=0; i<vec.dim(); i++){
				VAL_TYPE val = vec.GetVal(i);
				os << "  " << val << "," << std::endl;
			}
			return os << "|" << std::endl;
		}

		//friend Vector operator+(Vector vec1, Vector vec2){
		//	assert(vec1.dim() == vec2.dim());

		//	int dim = vec1.dim();
		//	Vector ret(dim);
		//	int i;
		//	for(i=0; i<dim; i++){
		//		ret.SetVal(i, vec1.GetVal(i) + vec2.GetVal(i));
		//	}
		//	return ret;
		//}

	protected:
		VAL_TYPE* _d;

	private:
		int _dim; //dimensional
	};

	template <class VAL_TYPE>
	Vector<VAL_TYPE>::Vector(int dimensional)
	{
		_d = NULL;

		Init(dimensional);
	}

	template <class VAL_TYPE>
	Vector<VAL_TYPE>::Vector()
	{
		_dim = 0;
		_d = NULL;
	}

	template <class VAL_TYPE>
	Vector<VAL_TYPE>::Vector(const Vector &vec)
	{
		_d = NULL;
		Init(vec.dim());

		int i;
		for(i=0; i<vec.dim(); i++){
			SetVal(i, vec.GetVal(i));
		}
	}

	template <class VAL_TYPE>
	Vector<VAL_TYPE>::~Vector()
	{
		if(_d){
			delete[] _d;
			_d = NULL;
		}
	}

	template <class VAL_TYPE>
	void Vector<VAL_TYPE>::Init(int dimensional)
	{
		if(_d){
			if(_dim != dimensional){
				delete[] _d;
			}
		}

		_dim = dimensional;
		_d = new VAL_TYPE[_dim];
	}

	template <class VAL_TYPE>
	inline int Vector<VAL_TYPE>::dim() const
	{
		return _dim;
	}

	template <class VAL_TYPE>
	inline void Vector<VAL_TYPE>::SetVal(int i, const VAL_TYPE& val)
	{
		_d[i] = val;
	}

	template <class VAL_TYPE>
	inline VAL_TYPE& Vector<VAL_TYPE>::GetVal(int i) const
	{
		assert((0 <= i) && (i < this->dim()));
		return _d[i];
	}

	template <class VAL_TYPE>
	inline void Vector<VAL_TYPE>::GetVal(int i, VAL_TYPE &out) const
	{
		assert((0 <= i) && (i < this->dim()));
		out = _d[i];
	}

	template <class VAL_TYPE>
	inline void Vector<VAL_TYPE>::Swap(int i, int j)
	{
		VAL_TYPE tmp;
		tmp = _d[i];
		_d[i] = _d[j];
		_d[j] = tmp;
	}

	template <class VAL_TYPE>
	void Vector<VAL_TYPE>::CopyTo(Vector<VAL_TYPE>& vec) const
	{
		vec.Init(this->dim());

		int i;
		for(i=0; i<this->dim(); i++){
			vec.SetVal(i, this->GetVal(i));
		}
	}

	template <class VAL_TYPE>
	inline Vector<VAL_TYPE>& Vector<VAL_TYPE>::operator =(const Vector<VAL_TYPE>& vec)
	{
		Init(vec.dim());

		int i;
		for(i=0; i<vec.dim(); i++){
			this->SetVal(i, vec.GetVal(i));
		}
		return *this;
	}

	template <class VAL_TYPE>
	inline VAL_TYPE& Vector<VAL_TYPE>::operator [](int index) const
	{
		assert((0 <= index) && (index < this->dim()));
		return _d[index];
	}

	typedef Vector< Double > Vectord;
	typedef Vector< Float > Vectorf;
	typedef Vector< Int > Vectori;
	typedef Vector< Complexd > Vectorcmd;

};

#endif //_SALVEC_H_
