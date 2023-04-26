#pragma once
#include <iostream>
#include <complex>

#include "entete objets.hpp"

template<class T, class U> std::complex<T> operator*(const U& scalaire, const std::complex<T>& temp) {
	return std::complex<T>(scalaire * temp.real, scalaire * temp.imag);
}

template<typename T> class complexe {
public:


	friend complexe<T> derivee(const complexe<T>& element) {
		return complexe<T>(derivee(element.x), derivee(element.y));
	};


	explicit complexe(T c_x,T c_y) : x(c_x),y(c_y) {};

	complexe() : x(), y() {};

	complexe(complexe<T>& copie) {
		x = copie.x;
		y = copie.y;
	};
	
	complexe(std::complex<T> const& temp) {
		x = temp.real();
		y = temp.imag();
	};

	explicit inline operator bool() {
		return (((bool)x) || ((bool)y));
	};

	/*
	template<class U> operator complexe<U>() {
		complexe<U> result;
		result.x = (U)x;
		result.y = (U)y;
		return result;
	};*/

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

	/*
	complexe<T>& operator=(bool test) {
		x = test;
		y = false;

		return *this;
	};
	*/

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

	T x;
	T y;
};

