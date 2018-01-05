#ifndef _SALCOMPLEXEF_H_
#define _SALCOMPLEXEF_H_

#include "exfloat.h"

class Complexef {
public:
	Complexef();
	Complexef(double scalar);
	Complexef(exfloat scalar);
	Complexef(exfloat re, exfloat im);
	~Complexef();

	void SetNum(exfloat re, exfloat im);

	exfloat Re();
	exfloat Im();

	//! •¡‘f”‚Ìâ‘Î’l‚ð•Ô‚·
	exfloat Abs() const;
	//! ‹¤–ð•¡‘f”‚ð•Ô‚·
	Complexef Conjugate() const;

	//
	void Add(const Complexef& complex_num);
	void Subtraction(const Complexef& complex_num);
	void Multiply(const Complexef& complex_num);
	void Multiply(exfloat scalar);
	void Division(const Complexef& complex_num);
	void Division(exfloat scalar);

	Complexef& operator=(exfloat scalar);
	Complexef operator+=(Complexef cm_num);
	Complexef operator+=(exfloat scalar);
	Complexef operator-=(Complexef cm_num);
	Complexef operator-=(exfloat scalar);

	bool operator>(Complexef cm_num);
	bool operator>(exfloat scalar);
	bool operator<(Complexef cm_num);
	bool operator<(exfloat scalar);
	bool operator>=(Complexef cm_num);
	bool operator>=(exfloat scalar);
	bool operator<=(Complexef cm_num);
	bool operator<=(exfloat scalar);

	friend Complexef operator+(Complexef cm_num1, Complexef cm_num2){
		Complexef temp;
		temp.SetNum(cm_num1._re, cm_num1._im);
		temp.Add(cm_num2);
		return temp;
	}
	friend Complexef operator+(Complexef cm_num, exfloat scalar){
		Complexef temp;
		temp.SetNum(cm_num._re, cm_num._im);
		temp.Add(scalar);
		return temp;
	}
	friend Complexef operator+(exfloat scalar, Complexef cm_num){
		Complexef temp;
		temp.SetNum(scalar, static_cast<exfloat>(0.0));
		temp.Add(cm_num);
		return temp;
	}

	friend Complexef operator-(Complexef cm_num1, Complexef cm_num2){
		Complexef temp;
		temp.SetNum(cm_num1._re, cm_num1._im);
		temp.Subtraction(cm_num2);
		return temp;
	}
	friend Complexef operator-(Complexef cm_num, exfloat scalar){
		Complexef temp;
		temp.SetNum(cm_num._re, cm_num._im);
		temp.Subtraction(scalar);
		return temp;
	}
	friend Complexef operator-(exfloat scalar, Complexef cm_num){
		Complexef temp;
		temp.SetNum(scalar, static_cast<exfloat>(0.0));
		temp.Subtraction(cm_num);
		return temp;
	}
	
	friend Complexef operator*(Complexef cm_num1, Complexef cm_num2){
		Complexef temp;
		temp.SetNum(cm_num1._re, cm_num1._im);
		temp.Multiply(cm_num2);
		return temp;
	}
	friend Complexef operator*(Complexef cm_num, exfloat scalar){
		Complexef temp;
		temp.SetNum(cm_num._re, cm_num._im);
		temp.Multiply(scalar);
		return temp;
	}
	friend Complexef operator*(exfloat scalar, Complexef cm_num){
		Complexef temp;
		temp.SetNum(scalar, static_cast<exfloat>(0.0));
		temp.Multiply(cm_num);
		return temp;
	}

	friend Complexef operator/(Complexef cm_num1, Complexef cm_num2){
		Complexef temp;
		temp.SetNum(cm_num1._re, cm_num1._im);
		temp.Division(cm_num2);
		return temp;
	}
	friend Complexef operator/(Complexef cm_num, exfloat scalar){
		Complexef temp;
		temp.SetNum(cm_num._re, cm_num._im);
		temp.Division(scalar);
		return temp;
	}
	friend Complexef operator/(exfloat scalar, Complexef cm_num){
		Complexef temp;
		temp.SetNum(scalar, static_cast<exfloat>(0.0));
		temp.Division(cm_num);
		return temp;
	}

	//void dump();
private:
	exfloat _re;
	exfloat _im;
};

Complexef::Complexef()
{
	_re = static_cast<exfloat>(0.0);
	_im = static_cast<exfloat>(0.0);
}

Complexef::Complexef(exfloat re, exfloat im)
{
	_re = re;
	_im = im;
}

Complexef::Complexef(double scalar)
{
	_re = static_cast<exfloat>(scalar);
	_im = 0;
}

Complexef::Complexef(exfloat scalar)
{
	_re = scalar;
	_im = 0;
}

Complexef::~Complexef()
{
}

void Complexef::SetNum(exfloat re, exfloat im)
{
	_re = re;
	_im = im;
}

exfloat Complexef::Re()
{
	return _re;
}

exfloat Complexef::Im()
{
	return _im;
}

exfloat Complexef::Abs() const
{
	return sqrt(_re*_re + _im*_im);
}

Complexef Complexef::Conjugate() const
{
	return Complexef(_re, (-1*_im));
}

void Complexef::Add(const Complexef& complex_num)
{
	_re += complex_num._re;
	_im += complex_num._im;
}

void Complexef::Subtraction(const Complexef& complex_num)
{
	_re -= complex_num._re;
	_im -= complex_num._im;
}

void Complexef::Multiply(const Complexef& complex_num)
{
	exfloat tmp = _re * complex_num._re - _im * complex_num._im;
	_im = _re * complex_num._im + _im * complex_num._re;
	_re = tmp;
}

void Complexef::Multiply(exfloat scalar)
{
	_re *= scalar;
	_im *= scalar;
}

void Complexef::Division(const Complexef& complex_num)
{
	exfloat d = complex_num._re*complex_num._re + complex_num._im*complex_num._im;
	Complexef cm = complex_num.Conjugate();
	this->Multiply( cm );
	this->Division( d );
}

void Complexef::Division(exfloat scalar)
{
	_re /= scalar;
	_im /= scalar;
}

Complexef& Complexef::operator =(exfloat scalar)
{
	_re = scalar;
	_im = static_cast<exfloat>(0.0);
	return *this;
}

Complexef Complexef::operator +=(Complexef cm_num)
{
	_re += cm_num._re;
	_im += cm_num._im;
	return *this;
}

Complexef Complexef::operator +=(exfloat scalar)
{
	_re += scalar;
	return *this;
}

Complexef Complexef::operator -=(Complexef cm_num)
{
	_re -= cm_num._re;
	_im -= cm_num._im;
	return *this;
}

Complexef Complexef::operator -=(exfloat scalar)
{
	_re -= scalar;
	return *this;
}

bool Complexef::operator >(Complexef cm_num)
{
	return (this->Abs() > cm_num.Abs());
}

bool Complexef::operator >(exfloat scalar)
{
	return (this->Abs() > scalar);
}

bool Complexef::operator <(Complexef cm_num)
{
	return (this->Abs() < cm_num.Abs());
}

bool Complexef::operator <(exfloat scalar)
{
	return (this->Abs() < scalar);
}

bool Complexef::operator >=(Complexef cm_num)
{
	return (this->Abs() >= cm_num.Abs());
}

bool Complexef::operator >=(exfloat scalar)
{
	return (this->Abs() >= scalar);
}

bool Complexef::operator <=(Complexef cm_num)
{
	return (this->Abs() <= cm_num.Abs());
}

bool Complexef::operator <=(exfloat scalar)
{
	return (this->Abs() <= scalar);
}

#endif //_SALCOMPLEXEF_H_