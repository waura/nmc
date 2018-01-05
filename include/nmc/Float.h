#ifndef _NMC_FLOAT_H_
#define _NMC_FLOAT_H_

#include <stdio.h>

#include <math.h>
#include <iostream>

class Float;

Float sqrt(Float& f);

namespace nmc {

class Float {
public:
	Float(): _num(0){};
	Float(const float d){_num = d;}
	~Float(){};

	inline void SetNum(const float d){
		_num = d;
	}
	inline float Abs() const {
		return fabs(_num);
	}

	inline operator float() const {
		return _num;
	}
	
	inline Float& operator=(const float d);
	inline Float& operator+=(const Float& d);
	inline Float& operator+=(const float d);
	inline Float& operator-=(const Float& d);
	inline Float& operator-=(const float d);
	inline Float& operator*=(const Float& d);
	inline Float& operator*=(const float d);
	inline Float& operator/=(const Float& d);
	inline Float& operator/=(const float d);

	inline bool operator>(const Float& d);
	inline bool operator>(const float d);
	inline bool operator<(const Float& d);
	inline bool operator<(const float d);
	inline bool operator>=(const Float& d);
	inline bool operator>=(const float d);
	inline bool operator<=(const Float& d);
	inline bool operator<=(const float d);

	inline friend Float operator+(const Float& d1, const Float& d2){
		return Float(d1._num + d2._num);
	}
	inline friend Float operator+(const Float& d1, const float d2){
		return Float(d1._num + d2);
	}
	inline friend Float operator+(const float d1, const Float& d2){
		return Float(d1 + d2._num);
	}

	inline friend Float operator-(const Float& d1, const Float& d2){
		return Float(d1._num - d2._num);
	}
	inline friend Float operator-(const Float& d1, const float d2){
		return Float(d1._num - d2);
	}
	inline friend Float operator-(const float d1, const Float& d2){
		return Float(d1 - d2._num);
	}

	inline friend Float operator*(const Float& d1, const Float& d2){
		return Float(d1._num * d2._num);
	}
	inline friend Float operator*(const Float& d1, const float d2){
		return Float(d1._num * d2);
	}
	inline friend Float operator*(const float d1, const Float& d2){
		return Float(d1 * d2._num);
	}

	inline friend Float operator/(const Float& d1, const Float& d2){
		return Float(d1._num / d2._num);
	}
	inline friend Float operator/(const Float& d1, const float d2){
		return Float(d1._num / d2);
	}
	inline friend Float operator/(const float d1, const Float& d2){
		return Float(d1 / d2._num);
	}

	inline friend std::ostream& operator<<(std::ostream& os, const Float& d){
		return os << d._num;
	}

private:
	float _num;
};


inline Float& Float::operator =(const float d)
{
	_num = d;
	return *this;
}

inline Float& Float::operator +=(const Float& d)
{
	_num += d._num;
	return *this;
}

inline Float& Float::operator +=(const float d)
{
	_num += d;
	return *this;
}

inline Float& Float::operator -=(const Float& d)
{
	_num -= d._num;
	return *this;
}
inline Float& Float::operator -=(const float d)
{
	_num -= d;
	return *this;
}

inline Float& Float::operator *=(const Float& d)
{
	_num *= d._num;
	return *this;
}

inline Float& Float::operator *=(const float d)
{
	_num *= d;
	return *this;
}

inline Float& Float::operator /=(const Float& d)
{
	_num /= d._num;
	return *this;
}

inline Float& Float::operator /=(const float d)
{
	_num /= d;
	return *this;
}

inline bool Float::operator >(const Float& d)
{
	return (_num > d._num);
}

inline bool Float::operator >(const float d)
{
	return (_num > d);
}

inline bool Float::operator <(const Float& d)
{
	return (_num < d._num);
}

inline bool Float::operator <(const float d)
{
	return (_num < d);
}

inline bool Float::operator >=(const Float& d)
{
	return (_num >= d._num);
}

inline bool Float::operator >=(const float d)
{
	return (_num >= d);
}

inline bool Float::operator <=(const Float& d)
{
	return (_num <= d._num);
}

inline bool Float::operator <=(const float d)
{
	return (_num <= d);
}

};

#endif //_NMC_DOUBLE_H_