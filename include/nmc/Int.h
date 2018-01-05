#ifndef _NMC_INT_H_
#define _NMC_INT_H_

#include <stdio.h>

#include <math.h>
#include <iostream>

class Int;

Int sqrt(Int& d);

namespace nmc {

class Int {
public:
	Int(): _num(0){};
	Int(int d){_num = d;}
	~Int(){};

	inline void SetNum(int d){
		_num = d;
	}
	inline int Abs() const {
		return (_num >= 0) ? _num : -_num;
	}

	inline operator int() const {
		return _num;
	}
	
	inline Int& operator=(int d);
	inline Int operator+=(Int d);
	inline Int operator+=(int d);
	inline Int operator-=(Int d);
	inline Int operator-=(int d);
	inline Int operator*=(Int d);
	inline Int operator*=(int d);
	inline Int operator/=(Int d);
	inline Int operator/=(int d);

	inline bool operator>(Int d);
	inline bool operator>(int d);
	inline bool operator<(Int d);
	inline bool operator<(int d);
	inline bool operator>=(Int d);
	inline bool operator>=(int d);
	inline bool operator<=(Int d);
	inline bool operator<=(int d);

	inline friend Int operator+(Int d1, Int d2){
		return Int(d1._num + d2._num);
	}
	inline friend Int operator+(Int d1, int d2){
		return Int(d1._num + d2);
	}
	inline friend Int operator+(int d1, Int d2){
		return Int(d1 + d2._num);
	}

	inline friend Int operator-(Int d1, Int d2){
		return Int(d1._num - d2._num);
	}
	inline friend Int operator-(Int d1, int d2){
		return Int(d1._num - d2);
	}
	inline friend Int operator-(int d1, Int d2){
		return Int(d1 - d2._num);
	}

	inline friend Int operator*(Int d1, Int d2){
		return Int(d1._num * d2._num);
	}
	inline friend Int operator*(Int d1, int d2){
		return Int(d1._num * d2);
	}
	inline friend Int operator*(int d1, Int d2){
		return Int(d1 * d2._num);
	}

	inline friend Int operator/(Int d1, Int d2){
		return Int(d1._num / d2._num);
	}
	inline friend Int operator/(Int d1, int d2){
		return Int(d1._num / d2);
	}
	inline friend Int operator/(int d1, Int d2){
		return Int(d1 / d2._num);
	}

	inline friend std::ostream& operator<<(std::ostream& os, const Int& d){
		return os << d._num;
	}

private:
	int _num;
};

inline Int& Int::operator =(int d)
{
	_num = d;
	return *this;
}

inline Int Int::operator +=(Int d)
{
	_num += d._num;
	return *this;
}

inline Int Int::operator +=(int d)
{
	_num += d;
	return *this;
}

inline Int Int::operator -=(Int d)
{
	_num -= d._num;
	return *this;
}

inline Int Int::operator -=(int d)
{
	_num -= d;
	return *this;
}

inline Int Int::operator *=(Int d)
{
	_num *= d._num;
	return *this;
}

inline Int Int::operator *=(int d)
{
	_num *= d;
	return *this;
}

inline Int Int::operator /=(Int d)
{
	_num /= d._num;
	return *this;
}

inline Int Int::operator /=(int d)
{
	_num /= d;
	return *this;
}

inline bool Int::operator >(Int d)
{
	return (_num > d._num);
}

inline bool Int::operator >(int d)
{
	return (_num > d);
}

inline bool Int::operator <(Int d)
{
	return (_num < d._num);
}

inline bool Int::operator <(int d)
{
	return (_num < d);
}

inline bool Int::operator >=(Int d)
{
	return (_num >= d._num);
}

inline bool Int::operator >=(int d)
{
	return (_num >= d);
}

inline bool Int::operator <=(Int d)
{
	return (_num <= d._num);
}

inline bool Int::operator <=(int d)
{
	return (_num <= d);
}

};

#endif //_NMC_INT_H_
