#ifndef _NMC_COMPLEX_H_
#define _NMC_COMPLEX_H_

#include "nmc/Double.h"
#include "nmc/Float.h"
#include "nmc/Int.h"

#include <cmath>
#include <iostream>

namespace nmc {

template <class VAL_TYPE>
class Complex {
public:
	Complex();
	Complex(int scalar);
	Complex(double scalar);
	Complex(const VAL_TYPE& scalar);
	Complex(const VAL_TYPE& re, const VAL_TYPE& im);
	~Complex();

	inline void SetNum(const VAL_TYPE& re, const VAL_TYPE& im);

	inline VAL_TYPE Re() const;
	inline VAL_TYPE Im() const;

	//! 複素数の絶対値を返す
	inline VAL_TYPE Abs() const;
	//! 共役複素数を返す
	inline Complex Conjugate() const;

	//
	inline void Add(const Complex& complex_num);
	inline void Subtraction(const Complex& complex_num);
	inline void Multiply(const Complex& complex_num);
	inline void Multiply(const VAL_TYPE& scalar);
	inline void Multiply(const double scalar);
	inline void Division(const Complex& complex_num);
	inline void Division(const VAL_TYPE& scalar);
	inline void Division(const double scalar);

	inline Complex& operator=(const Complex& cm_num);
	inline Complex& operator=(const VAL_TYPE& scalar);
	inline Complex& operator=(const double scalar);
	inline Complex& operator+=(const Complex& cm_num);
	inline Complex& operator+=(const VAL_TYPE& scalar);
	inline Complex& operator+=(const double scalar);
	inline Complex& operator-=(const Complex& cm_num);
	inline Complex& operator-=(const VAL_TYPE& scalar);
	inline Complex& operator-=(const double scalar);
	inline Complex& operator*=(const Complex& cm_num);
	inline Complex& operator*=(const VAL_TYPE& scalar);
	inline Complex& operator*=(const double scalar);
	inline Complex& operator/=(const Complex& cm_num);
	inline Complex& operator/=(const VAL_TYPE& scalar);
	inline Complex& operator/=(const double scalar);

	inline bool operator>(const Complex& cm_num);
	inline bool operator>(const VAL_TYPE& scalar);
	inline bool operator>(double scalar);
	inline bool operator<(const Complex& cm_num);
	inline bool operator<(const VAL_TYPE& scalar);
	inline bool operator<(double scalar);
	inline bool operator>=(const Complex& cm_num);
	inline bool operator>=(const VAL_TYPE& scalar);
	inline bool operator>=(double scalar);
	inline bool operator<=(const Complex& cm_num);
	inline bool operator<=(const VAL_TYPE& scalar);
	inline bool operator<=(double scalar);

	inline Complex operator+() const {
		return *this;
	}

	inline Complex operator-() const {
		return Complex(-_re, -_im);
	}

	inline friend Complex operator+(const Complex& cm_num1, const Complex& cm_num2){
		Complex temp;
		temp.SetNum(cm_num1._re, cm_num1._im);
		temp.Add(cm_num2);
		return temp;
	}
	inline friend Complex operator+(const Complex& cm_num, const VAL_TYPE& scalar){
		Complex temp;
		temp.SetNum(cm_num._re, cm_num._im);
		temp.Add(scalar);
		return temp;
	}
	inline friend Complex operator+(const VAL_TYPE& scalar, const Complex& cm_num){
		Complex temp;
		temp.SetNum(scalar, 0.0);
		temp.Add(cm_num);
		return temp;
	}

	inline friend Complex operator-(const Complex& cm_num1, const Complex& cm_num2){
		Complex temp;
		temp.SetNum(cm_num1._re, cm_num1._im);
		temp.Subtraction(cm_num2);
		return temp;
	}
	inline friend Complex operator-(const Complex& cm_num, const VAL_TYPE& scalar){
		Complex temp;
		temp.SetNum(cm_num._re, cm_num._im);
		temp.Subtraction(scalar);
		return temp;
	}
	inline friend Complex operator-(const VAL_TYPE& scalar, const Complex& cm_num){
		Complex temp;
		temp.SetNum(scalar, 0.0);
		temp.Subtraction(cm_num);
		return temp;
	}
	
	inline friend Complex operator*(const Complex& cm_num1, const Complex& cm_num2){
		Complex temp;
		temp.SetNum(cm_num1._re, cm_num1._im);
		temp.Multiply(cm_num2);
		return temp;
	}
	inline friend Complex operator*(const Complex& cm_num, const VAL_TYPE& scalar){
		Complex temp;
		temp.SetNum(cm_num._re, cm_num._im);
		temp.Multiply(scalar);
		return temp;
	}
	inline friend Complex operator*(const VAL_TYPE& scalar, const Complex& cm_num){
		Complex temp;
		temp.SetNum(scalar, 0.0);
		temp.Multiply(cm_num);
		return temp;
	}
	inline friend Complex operator*(const Complex& cm_num, const double scalar){
		Complex temp;
		temp.SetNum(cm_num._re, cm_num._im);
		temp.Multiply(scalar);
		return temp;
	}
	inline friend Complex operator*(const double scalar, const Complex& cm_num){
		Complex temp;
		temp.SetNum(cm_num._re, cm_num._im);
		temp.Multiply(scalar);
		return temp;
	}

	inline friend Complex operator/(const Complex& cm_num1, const Complex& cm_num2){
		Complex temp;
		temp.SetNum(cm_num1._re, cm_num1._im);
		temp.Division(cm_num2);
		return temp;
	}
	inline friend Complex operator/(const Complex& cm_num, const VAL_TYPE& scalar){
		Complex temp;
		temp.SetNum(cm_num._re, cm_num._im);
		temp.Division(scalar);
		return temp;
	}
	inline friend Complex operator/(const VAL_TYPE& scalar, const Complex& cm_num){
		Complex temp;
		temp.SetNum(scalar, 0.0);
		temp.Division(cm_num);
		return temp;
	}

	inline friend std::ostream& operator<<(std::ostream& os, const Complex& cm){
		return os << "(" << cm._re << ")" << " + i(" << cm._im << ")";
	}

private:
	VAL_TYPE _re;
	VAL_TYPE _im;
};

template <class VAL_TYPE>
Complex<VAL_TYPE>::Complex()
{
	_re = static_cast<VAL_TYPE>(0.0);
	_im = static_cast<VAL_TYPE>(0.0);
}

template <class VAL_TYPE>
Complex<VAL_TYPE>::Complex(const VAL_TYPE& re, const VAL_TYPE& im)
{
	_re = re;
	_im = im;
}

template <class VAL_TYPE>
Complex<VAL_TYPE>::Complex(int scalar)
{
	_re = static_cast<VAL_TYPE>(scalar);
	_im = static_cast<VAL_TYPE>(0.0);
}

template <class VAL_TYPE>
Complex<VAL_TYPE>::Complex(double scalar)
{
	_re = static_cast<VAL_TYPE>(scalar);
	_im = static_cast<VAL_TYPE>(0.0);
}

template <class VAL_TYPE>
Complex<VAL_TYPE>::Complex(const VAL_TYPE& scalar)
{
	_re = scalar;
	_im = static_cast<VAL_TYPE>(0.0);
}

template <class VAL_TYPE>
Complex<VAL_TYPE>::~Complex()
{
}

template <class VAL_TYPE>
inline void Complex<VAL_TYPE>::SetNum(const VAL_TYPE& re, const VAL_TYPE& im)
{
	_re = re;
	_im = im;
}

template <class VAL_TYPE>
inline VAL_TYPE Complex<VAL_TYPE>::Re() const
{
	return _re;
}

template <class VAL_TYPE>
inline VAL_TYPE Complex<VAL_TYPE>::Im() const
{
	return _im;
}

template <class VAL_TYPE>
inline VAL_TYPE Complex<VAL_TYPE>::Abs() const
{
	return sqrt(_re*_re + _im*_im);
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE> Complex<VAL_TYPE>::Conjugate() const
{
	return Complex(_re, (-1.0*_im));
}

template <class VAL_TYPE>
inline void Complex<VAL_TYPE>::Add(const Complex& complex_num)
{
	_re += complex_num._re;
	_im += complex_num._im;
}

template <class VAL_TYPE>
inline void Complex<VAL_TYPE>::Subtraction(const Complex& complex_num)
{
	_re -= complex_num._re;
	_im -= complex_num._im;
}

template <class VAL_TYPE>
inline void Complex<VAL_TYPE>::Multiply(const Complex& complex_num)
{
	VAL_TYPE tmp = _re * complex_num._re - _im * complex_num._im;
	_im = _re * complex_num._im + _im * complex_num._re;
	_re = tmp;
}

template <class VAL_TYPE>
inline void Complex<VAL_TYPE>::Multiply(const VAL_TYPE& scalar)
{
	_re *= scalar;
	_im *= scalar;
}

template <class VAL_TYPE>
inline void Complex<VAL_TYPE>::Multiply(double scalar)
{
	_re *= scalar;
	_im *= scalar;
}

template <class VAL_TYPE>
inline void Complex<VAL_TYPE>::Division(const Complex& complex_num)
{
	VAL_TYPE d = complex_num._re*complex_num._re + complex_num._im*complex_num._im;
	Complex cm = complex_num.Conjugate();
	this->Multiply( cm );
	this->Division( d );
}

template <class VAL_TYPE>
inline void Complex<VAL_TYPE>::Division(const VAL_TYPE& scalar)
{
	_re /= scalar;
	_im /= scalar;
}

template <class VAL_TYPE>
inline void Complex<VAL_TYPE>::Division(double scalar)
{
	_re /= scalar;
	_im /= scalar;
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE>& Complex<VAL_TYPE>::operator =(const Complex<VAL_TYPE>& cm_num)
{
	_re = cm_num._re;
	_im = cm_num._im;
	return *this;
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE>& Complex<VAL_TYPE>::operator =(const VAL_TYPE& scalar)
{
	_re = scalar;
	_im = 0.0;
	return *this;
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE>& Complex<VAL_TYPE>::operator =(const double scalar)
{
	_re = scalar;
	_im = 0.0;
	return *this;
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE>& Complex<VAL_TYPE>::operator +=(const Complex<VAL_TYPE>& cm_num)
{
	_re += cm_num._re;
	_im += cm_num._im;
	return *this;
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE>& Complex<VAL_TYPE>::operator +=(const VAL_TYPE& scalar)
{
	_re += scalar;
	return *this;
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE>& Complex<VAL_TYPE>::operator +=(double scalar)
{
	_re += scalar;
	return *this;
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE>& Complex<VAL_TYPE>::operator -=(const Complex<VAL_TYPE>& cm_num)
{
	_re -= cm_num._re;
	_im -= cm_num._im;
	return *this;
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE>& Complex<VAL_TYPE>::operator -=(const VAL_TYPE& scalar)
{
	_re -= scalar;
	return *this;
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE>& Complex<VAL_TYPE>::operator -=(double scalar)
{
	_re -= scalar;
	return *this;
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE>& Complex<VAL_TYPE>::operator *=(const Complex<VAL_TYPE>& cm_num)
{
	this->Multiply(cm_num);
	return *this;
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE>& Complex<VAL_TYPE>::operator *=(const VAL_TYPE& scalar)
{
	this->Multiply(scalar);
	return *this;
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE>& Complex<VAL_TYPE>::operator *=(double scalar)
{
	this->Multiply(scalar);
	return *this;
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE>& Complex<VAL_TYPE>::operator /=(const Complex<VAL_TYPE>& cm_num)
{
	this->Division(cm_num);
	return *this;
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE>& Complex<VAL_TYPE>::operator /=(const VAL_TYPE& scalar)
{
	this->Division(scalar);
	return *this;
}

template <class VAL_TYPE>
inline Complex<VAL_TYPE>& Complex<VAL_TYPE>::operator /=(double scalar)
{
	this->Division(scalar);
	return *this;
}

template <class VAL_TYPE>
inline bool Complex<VAL_TYPE>::operator >(const Complex<VAL_TYPE>& cm_num)
{
	return (this->Abs() > cm_num.Abs());
}

template <class VAL_TYPE>
inline bool Complex<VAL_TYPE>::operator >(const VAL_TYPE& scalar)
{
	return (this->Abs() > scalar);
}

template <class VAL_TYPE>
inline bool Complex<VAL_TYPE>::operator >(double scalar)
{
	return (this->Abs() > static_cast<VAL_TYPE>(scalar));
}

template <class VAL_TYPE>
inline bool Complex<VAL_TYPE>::operator <(const Complex<VAL_TYPE>& cm_num)
{
	return (this->Abs() < cm_num.Abs());
}

template <class VAL_TYPE>
inline bool Complex<VAL_TYPE>::operator <(const VAL_TYPE& scalar)
{
	return (this->Abs() < scalar);
}

template <class VAL_TYPE>
inline bool Complex<VAL_TYPE>::operator <(double scalar)
{
	return (this->Abs() < static_cast<VAL_TYPE>(scalar));
}

template <class VAL_TYPE>
inline bool Complex<VAL_TYPE>::operator >=(const Complex<VAL_TYPE>& cm_num)
{
	return (this->Abs() >= cm_num.Abs());
}

template <class VAL_TYPE>
inline bool Complex<VAL_TYPE>::operator >=(const VAL_TYPE& scalar)
{
	return (this->Abs() >= scalar);
}

template <class VAL_TYPE>
inline bool Complex<VAL_TYPE>::operator >=(double scalar)
{
	return (this->Abs() >= static_cast<VAL_TYPE>(scalar));
}

template <class VAL_TYPE>
inline bool Complex<VAL_TYPE>::operator <=(const Complex<VAL_TYPE>& cm_num)
{
	return (this->Abs() <= cm_num.Abs());
}

template <class VAL_TYPE>
inline bool Complex<VAL_TYPE>::operator <=(const VAL_TYPE& scalar)
{
	return (this->Abs() <= scalar);
}

template <class VAL_TYPE>
inline bool Complex<VAL_TYPE>::operator <=(double scalar)
{
	return (this->Abs() <= static_cast<VAL_TYPE>(scalar));
}

template <class VAL_TYPE>
VAL_TYPE abs(const Complex<VAL_TYPE>& cm)
{
        return nmc::sqrt(cm.Re()*cm.Re() + cm.Im()*cm.Im());
}

template <class VAL_TYPE>
double carg(const Complex<VAL_TYPE>& cm)
{
	if(abs(cm) == 0.0) return 0.0;
	return atan2(cm.Im(), cm.Re());
}

template <class VAL_TYPE>
Complex<VAL_TYPE> pow(const Complex<VAL_TYPE>& cm, int n)
{
	VAL_TYPE t = carg( cm ) * n;
	VAL_TYPE w = std::pow(abs(cm), n);
	return Complex<VAL_TYPE>(w * cos(t), w * sin(t));
}

template <class VAL_TYPE>
Complex<VAL_TYPE> pow(const Complex<VAL_TYPE>& cm, double n)
{
	VAL_TYPE t = carg( cm ) * n;
	VAL_TYPE w = pow(abs(cm), n);
	return Complex<VAL_TYPE>(w * cos(t), w * sin(t));
}

template <class VAL_TYPE>
Complex<VAL_TYPE> log(const Complex<VAL_TYPE>& cm)
{
	return Complex<VAL_TYPE>(std::log( abs(cm) ), carg( cm ));
}


typedef Complex< Double> Complexd;
typedef Complex<Float> Complexf;
typedef Complex<Int> Complexi;

};

#endif //_SALCOMPLEX_H_
