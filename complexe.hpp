#pragma once
#include <iostream>
#include <complex>

#include "entete objets.hpp"
#include "swap_T.hpp"

template<class T, class U> std::complex<T> operator*(const U& scalaire, const std::complex<T>& temp) {
	return std::complex<T>(scalaire * temp.real, scalaire * temp.imag);
};

template<class T>
std::complex<T> derivee(const std::complex<T>& element) {
	return std::complex<T>(derivee(element.real), derivee(element.imag));
};

template<typename T> class complexe {
public:
	T x;
	T y;

	friend complexe<T> derivee(const complexe<T>& element) {
		return complexe<T>(derivee(element.x), derivee(element.y));
	};


	explicit complexe(T c_x,T c_y) : x(c_x),y(c_y) {};

	complexe() : x(), y() {};

	complexe(complexe<T> const& copie) : x(copie.x), y(copie.y) {	};

	complexe(complexe<T>&& temp) {
		swap(*this, temp);
		return;
	};
	
	complexe(std::complex<T> const& temp) : x(temp.real) , y(temp.imag) {	};

	explicit inline operator bool() const {
		return (((bool)x) || ((bool)y));
	};

	complexe<T>& operator=(complexe<T> const& temp);

	complexe<T>& operator=(complexe<T>&& temp) {
		if (this == &temp)
			return *this;
		swap(*this, temp);
		return *this;
	};

	complexe<T>& operator*=(complexe<T> const& temp) {
		return (*this = (*this * temp));
	};

	template<class U>
	complexe<T>& operator*=(const U& scalaire) {
		x *= scalaire;
		y *= scalaire;
		return *this;
	};

	complexe<T>& operator+=(complexe<T> const& autre) {
		x += autre.x;
		y += autre.y;
		return *this;
	};

	friend complexe<T> operator*(const complexe<T>& z1, const complexe<T>& z2) {
		return complexe<T>((z1.x * z1.x) - (z1.y * z2.y), (z1.x * z2.y) + (z1.y * z2.x));
	};

	template<class U> friend complexe<T> operator*(const U& scalaire, const complexe<T>& z) {
		complexe<T> result(z);
		result.x = scalaire * result.x;
		result.y = scalaire * result.y;
		return result;
	};

	friend complexe<T> operator+(const complexe<T>& z1, const complexe<T>& z2) {
		return complexe<T>(z1.x + z2.x, z1.y + z2.y);
	};

	friend complexe<T> operator-(const complexe<T>& z1, const complexe<T>& z2) {
		return complexe<T>(z1.x - z2.x, z1.y - z2.y);
	};

	friend complexe<T> operator-(const complexe<T>& z1) {
		return complexe<T>(-z1.x, -z1.y);
	};

	complexe<T> conjugue() {
		return complexe<T>(x, -y);
	};

	friend complexe<T> operator/(const complexe<T>& z1, const complexe<T>& z2) {
		T vrai_ = unite(z1,true);

		T inv_module = vrai_ / ((z2.x * z2.x) + (z2.y * z2.y));
		return complexe<T>(inv_module * (z1.x * z2.x + z1.y * z2.y), inv_module * ((z1.y * z2.x) - (z1.x * z2.y)));
	};

	T norme() const {
		return (x * x) + (y * y);
	}

	friend bool operator==(const complexe<T>& temp1, const complexe<T>& temp2) {
		return ((temp1.x == temp2.x) && (temp1.y == temp2.y));
	}

	friend std::ostream& operator<<(std::ostream& os, const complexe<T>& element) {
		os << "(" << element.x << " + i*" << element.y << ")";
		return os;
	};

	
	operator std::complex<T>() {
		std::complex<T> temp(x, y);
		return temp;
	};

	friend void swap(complexe<T>& gauche, complexe<T>& droit) {
		swap_F(gauche.x, droit.x);
		swap_F(gauche.y, droit.y);
		return;
	};

};

template<class T>
complexe<T>& complexe<T>::operator=(complexe<T> const& temp) {
	if (this == &temp)
		return *this;
	x = temp.x;
	y = temp.y;
	return *this;
};
