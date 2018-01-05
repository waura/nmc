#ifndef _NMC_DOUBLE_H_
#define _NMC_DOUBLE_H_

#include <stdio.h>

#include <cmath>
#include <iostream>

class Double;

namespace nmc
{

class Double {
public:
	Double(): _num(0){};
	Double(const double d){_num = d;}
	Double(const int d){_num = d;}
	~Double(){};

	inline void SetNum(const double d){
		_num = d;
	}
	inline double Abs() const {
		return fabs(_num);
	}

	inline operator double() const{
		return _num;
	}
	
	inline Double& operator=(const double d);
	inline Double& operator=(const int d);
	inline Double& operator+=(const Double& d);
	inline Double& operator+=(const double d);
	inline Double& operator-=(const Double& d);
	inline Double& operator-=(const double d);
	inline Double& operator*=(const Double& d);
	inline Double& operator*=(const double d);
	inline Double& operator/=(const Double& d);
	inline Double& operator/=(const double d);

	inline bool operator>(const Double& d);
	inline bool operator>(const double d);
	inline bool operator<(const Double& d);
	inline bool operator<(const double d);
	inline bool operator>=(const Double& d);
	inline bool operator>=(const double d);
	inline bool operator<=(const Double& d);
	inline bool operator<=(const double d);

	inline friend Double operator+(const Double& d1, const Double& d2){
		return Double(d1._num + d2._num);
	}
	inline friend Double operator+(const Double& d1, const double d2){
		return Double(d1._num + d2);
	}
	inline friend Double operator+(const double d1, const Double& d2){
		return Double(d1 + d2._num);
	}

	inline friend Double operator-(const Double& d1, const Double& d2){
		return Double(d1._num - d2._num);
	}
	inline friend Double operator-(const Double& d1, const double d2){
		return Double(d1._num - d2);
	}
	inline friend Double operator-(const double d1, const Double& d2){
		return Double(d1 - d2._num);
	}

	inline friend Double operator*(const Double& d1, const Double& d2){
		return Double(d1._num * d2._num);
	}
	inline friend Double operator*(const Double& d1, const double d2){
		return Double(d1._num * d2);
	}
	inline friend Double operator*(const double d1, const Double& d2){
		return Double(d1 * d2._num);
	}

	inline friend Double operator/(const Double& d1, const Double& d2){
		return Double(d1._num / d2._num);
	}
	inline friend Double operator/(const Double& d1, const double d2){
		return Double(d1._num / d2);
	}
	inline friend Double operator/(const double d1, const Double& d2){
		return Double(d1 / d2._num);
	}

	inline friend std::ostream& operator<<(std::ostream& os, const Double& d){
		return os << d._num;
	}

private:
	double _num;
};


inline Double& Double::operator =(const double d)
{
	_num = d;
	return *this;
}

inline Double& Double::operator =(const int d)
{
	_num = d;
	return *this;
}

inline Double& Double::operator +=(const Double& d)
{
	_num += d._num;
	return *this;
}

inline Double& Double::operator +=(const double d)
{
	_num += d;
	return *this;
}

inline Double& Double::operator -=(const Double& d)
{
	_num -= d._num;
	return *this;
}
inline Double& Double::operator -=(const double d)
{
	_num -= d;
	return *this;
}

inline Double& Double::operator *=(const Double& d)
{
	_num *= d._num;
	return *this;
}

inline Double& Double::operator *=(const double d)
{
	_num *= d;
	return *this;
}

inline Double& Double::operator /=(const Double& d)
{
	_num /= d._num;
	return *this;
}

inline Double& Double::operator /=(const double d)
{
	_num /= d;
	return *this;
}

inline bool Double::operator >(const Double& d)
{
	return (_num > d._num);
}

inline bool Double::operator >(const double d)
{
	return (_num > d);
}

inline bool Double::operator <(const Double& d)
{
	return (_num < d._num);
}

inline bool Double::operator <(const double d)
{
	return (_num < d);
}

inline bool Double::operator >=(const Double& d)
{
	return (_num >= d._num);
}

inline bool Double::operator >=(const double d)
{
	return (_num >= d);
}

inline bool Double::operator <=(const Double& d)
{
	return (_num <= d._num);
}

inline bool Double::operator <=(const double d)
{
	return (_num <= d);
}


inline Double sqrt(const Double& d){
	double dd = static_cast<double>(d);
	dd = std::sqrt(dd);
	return dd;
}

};

#endif //_NMC_DOUBLE_H_