#ifndef _NMC_VECTOR2D_H_
#define _NMC_VECTOR2D_H_

#include "nmc/Vector.h"

namespace nmc
{

	template <class VAL_TYPE>
	class Vector2D : public Vector<VAL_TYPE>{
	public:
		Vector2D(){
			this->Init(2);
		}
		virtual ~Vector2D(){}

		inline bool operator>(const Vector<VAL_TYPE>& vec) const;
		inline bool operator<(const Vector<VAL_TYPE>& vec) const;
		inline bool operator>=(const Vector<VAL_TYPE>& vec) const;
		inline bool operator<=(const Vector<VAL_TYPE>& vec) const;
	};

	template <class VAL_TYPE>
	inline bool Vector2D<VAL_TYPE>::operator >(const nmc::Vector<VAL_TYPE> &vec) const
	{
		return (this->_d[0]*this->_d[0] + this->_d[1]*this->_d[1]) > (vec._d[0]*vec._d[0] + vec._d[1]*vec._d[1]);
	}

	template <class VAL_TYPE>
	inline bool Vector2D<VAL_TYPE>::operator <(const nmc::Vector<VAL_TYPE> &vec) const
	{
		return (this->_d[0]*this->_d[0] + this->_d[1]*this->_d[1]) < (vec._d[0]*vec._d[0] + vec._d[1]*vec._d[1]);
	}

	template <class VAL_TYPE>
	inline bool Vector2D<VAL_TYPE>::operator >=(const nmc::Vector<VAL_TYPE> &vec) const
	{
		return (this->_d[0]*this->_d[0] + this->_d[1]*this->_d[1]) >= (vec._d[0]*vec._d[0] + vec._d[1]*vec._d[1]);
	}

	template <class VAL_TYPE>
	inline bool Vector2D<VAL_TYPE>::operator <=(const nmc::Vector<VAL_TYPE> &vec) const
	{
		return (this->_d[0]*this->_d[0] + this->_d[1]*this->_d[1]) <= (vec._d[0]*vec._d[0] + vec._d[1]*vec._d[1]);
	}

	typedef Vector2D< Double > Vector2Dd;
	typedef Vector2D< Float > Vector2Df;
	typedef Vector2D< Int > Vector2Di;
	typedef Vector2D< Complexd > Vector2Dcmd;

	typedef Vector< Vector2Dd > VectorPoint2Dd;

};

#endif //_NMC_VECTOR2D_H_
