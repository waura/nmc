#ifndef _NMC_VECTOR3D_H_
#define _NMC_VECTOR3D_H_

#include "nmc/Vector.h"

namespace nmc
{

	template <class VAL_TYPE>
	class Vector3D : public Vector<VAL_TYPE>{
	public:
		Vector3D(){
			this->Init(3);
		}
		virtual ~Vector3D(){}

		inline bool operator>(const Vector<VAL_TYPE>& vec);
		inline bool operator<(const Vector<VAL_TYPE>& vec);
		inline bool operator>=(const Vector<VAL_TYPE>& vec);
		inline bool operator<=(const Vector<VAL_TYPE>& vec);
	};

	template <class VAL_TYPE>
	inline bool Vector3D<VAL_TYPE>::operator >(const nmc::Vector<VAL_TYPE> &vec)
	{
		return (_d[0]*_d[0] + _d[1]*_d[1] + _d[2]*_d[2]) >
			(vec._d[0]*vec._d[0] + vec._d[1]*vec._d[1] + vec._d[2]*vec._d[2]);
	}

	template <class VAL_TYPE>
	inline bool Vector3D<VAL_TYPE>::operator <(const nmc::Vector<VAL_TYPE> &vec)
	{
		return (_d[0]*_d[0] + _d[1]*_d[1] + _d[2]*_d[2]) <
			(vec._d[0]*vec._d[0] + vec._d[1]*vec._d[1] + vec._d[2]*vec._d[2]);
	}

	template <class VAL_TYPE>
	inline bool Vector3D<VAL_TYPE>::operator >=(const nmc::Vector<VAL_TYPE> &vec)
	{
		return (_d[0]*_d[0] + _d[1]*_d[1] + _d[2]*_d[2]) >=
			(vec._d[0]*vec._d[0] + vec._d[1]*vec._d[1] + vec._d[2]*vec._d[2]);
	}

	template <class VAL_TYPE>
	inline bool Vector3D<VAL_TYPE>::operator <=(const nmc::Vector<VAL_TYPE> &vec)
	{
		return (_d[0]*_d[0] + _d[1]*_d[1] + _d[2]*_d[2]) <=
			(vec._d[0]*vec._d[0] + vec._d[1]*vec._d[1] + vec._d[2]*vec._d[2]);
	}

	typedef Vector3D< Double > Vector3Dd;
	typedef Vector3D< Float > Vector3Df;
	typedef Vector3D< Int > Vector3Di;
	typedef Vector3D< Complexd > Vector3Dcmd;

	typedef Vector< Vector3Dd > VectorVector3Dd;

};

#endif //_NMC_VECTOR3D_H_