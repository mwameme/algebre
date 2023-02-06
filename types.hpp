#pragma once

#include <complex>
#include "polynome.hpp"
#include "rationnel.hpp"
#include "corps quotient.hpp"
#include "anneau quotient.hpp"
#include "InfInt.h"
#include "erreur.hpp"
#include "polynome_n.hpp"
#include "complexe.hpp"
#include "precision/fprecision.h"
#include "polynome_n.hpp"
#include "polynome_n.hpp"

#include "entete objets.hpp"

template<class T> class complex;

template<class T> class erreur;
template<class T> class erreur_l;
template<class T> class anneau_quotient;
template<typename T> class complexe;
template<class T> class corps_quotient;
template<class T, class enable1 = void, class enable2 = void> class matrice;
template<typename T> class polynome;
template<class T> class polynome_n;
template<class T> class polynome_n_iter;
template<class T, class enable = void> class rationnel;

/*
template<typename, typename = void> constexpr int type_algebre{};
template<typename T> constexpr int type_algebre < rationnel<T>> = 0;
template<typename T> constexpr int type_algebre < polynome<T>, type_algebre<T> == 0> = 1;
*/

template<class T> class type_algebre {
public:
	static constexpr int get_type() { return 2; };
	static constexpr int get_approx() { return 0; };

	static constexpr int type{ get_type() }; //0 : division exacte ; 1 : division avec reste ; 2 : sans divisio,
	static constexpr int approx{ get_approx() }; //0 : exact ; 1 : approché 

};

template<class T> class type_algebre<polynome<T>> {
public:
	static constexpr int get_type() {
		if (type_algebre<T>::type == 0)
			return 1;
		else return 2;
	};
	static constexpr int get_approx() {
		return type_algebre<T>::approx;
	};

	static constexpr int type{ get_type() };
	static constexpr int approx{ get_approx() };
};

template<class T> class type_algebre<polynome_n<T>> {
public:
	static constexpr int get_type() {
		return 2;
	};
	static constexpr int get_approx() {
		return type_algebre<T>::approx;
	};

	static constexpr int type{ get_type() };
	static constexpr int approx{ get_approx() };
};

template<class T> class type_algebre<polynome_n_iter<T>> {
public:
	static constexpr int get_type() {
		return 2;
	};
	static constexpr int get_approx() {
		return type_algebre<T>::approx;
	};

	static constexpr int type{ get_type() };
	static constexpr int approx{ get_approx() };
};



template<class T> class type_algebre<anneau_quotient<T>> {
public:
	static constexpr int get_type() {
		return type_algebre<T>::type;
	};
	static constexpr int get_approx() {
		return type_algebre<T>::approx;
	};

	static constexpr int type{ get_type() };
	static constexpr int approx{ get_approx() };
};

template<class T> class type_algebre<corps_quotient<T>> {
public:
	static constexpr int get_type() {
		return 0;
	};

	static constexpr int get_approx() {
		return type_algebre<T>::approx;
	};

	static constexpr int type{ get_type() };
	static constexpr int approx{ get_approx() };
};

template<class T> class type_algebre<rationnel<T>> {
public:
	static constexpr int get_type() {
		return 0;
	};

	static constexpr int get_approx() {
		return type_algebre<T>::approx;
	};

	static constexpr int type{ get_type() };
	static constexpr int approx{ get_approx() };
};

template<class T> class type_algebre<erreur<T>> {
public:
	static constexpr int get_type() {
		return 0;
	};

	static constexpr int get_approx() {
		return type_algebre<T>::approx;
	};

	static constexpr int type{ get_type() };
	static constexpr int approx{ get_approx() };
};

template<class T> class type_algebre<erreur_l<T>> {
public:
	static constexpr int get_type() {
		return 0;
	};

	static constexpr int get_approx() {
		return type_algebre<T>::approx;
	};

	static constexpr int type{ get_type() };
	static constexpr int approx{ get_approx() };
};


template<class T> class type_algebre<complexe<T>> {
public:
	static constexpr int get_type() {
		if (type_algebre<T>::type == 0)
			return 0;
		else
			return 2;
	};

	static constexpr int get_approx() {
		return type_algebre<T>::approx;
	};

	static constexpr int type{ get_type() };
	static constexpr int approx{ get_approx() };

};

template<class T> class type_algebre<complex<T>> {
public:
	static constexpr int get_type() {
		if (type_algebre<T>::type == 0)
			return 0;
		else
			return 2;
	};

	static constexpr int get_approx() {
		return type_algebre<T>::approx;
	};

	static constexpr int type{ get_type() };
	static constexpr int approx{ get_approx() };

};


// ==========================================================
//types de base
template<> class type_algebre<int> {
public:
	static constexpr int type = 1;
	static constexpr int approx = 0;

};

template<> class type_algebre<long> {
public:
	static constexpr int type = 1;
	static constexpr int approx = 0;
};

template<> class type_algebre<long long> {
public:
	static constexpr int type = 1;
	static constexpr int approx = 0;

};



template<> class type_algebre<int_precision> {
public:
	static constexpr int type = 1;
	static constexpr int approx = 0;
};

template<> class type_algebre<InfInt> {
public:
	static constexpr int type = 1;
	static constexpr int approx = 0;

};

template<> class type_algebre<float> {
public:
	static constexpr int type = 0;
	static constexpr int approx = 1;

};

template<> class type_algebre<double> {
public:
	static constexpr int type = 0;
	static constexpr int approx = 1;

};
template<> class type_algebre<float_precision> {
public:
	static constexpr int type = 0;
	static constexpr int approx = 1;
};




//template type_algebre<anneau_quotient<polynome<rationnel<int>>>>()



/*

//0 : corps. 1 : anneau avec division et reste. 2 : anneau sans division
constexpr int type_algebre(int element) {
	return 1;
};

constexpr int type_algebre(long element) {
	return 1;
};

constexpr int type_algebre(long long element) {
	return 1;
};


template<class T> constexpr int type_algebre(InfInt& element) {
	return 1;
};

template<class T> constexpr int type_algebre(int_precision& element) {
	return 1;
};

constexpr int type_algebre(float element) {
	return 0;
};

constexpr int type_algebre(double element) {
	return 0;
};

constexpr int type_algebre(float_precision& element) {
	return 0;
};



template<class T> constexpr int type_algebre(polynome<T>& element) {
	if (type_algebre(T()) == 0)
		return 1;
	else
		return 2;
};

template<class T> constexpr int type_algebre(rationnel<T>& element) {
	return 0;
};


template<class T> constexpr int type_algebre(anneau_quotient<T>& element) {
	return (type_algebre(T()));
};


template<class T> constexpr int type_algebre(corps_quotient<T>& element) {
	if (type_algebre(T()) == 0)
		return 0;
	else
		if (type_alebre(T()) == 1)
			return 0;
	return 2; //interdit
};


template<class T> constexpr int type_algebre(matrice<T>& element) {
	if (type_algebre(T()) == 0)
		return 0;
	else
		return 2;
};

template<class T> constexpr int type_algebre(polynome_n<T>& element) {
	return 2;
};

template<class T> constexpr int type_algebre(complexe<T>& element) {
	if (type_algebre(T()) == 0)
		return 0;
	else
		return 2;
};

template<class T> constexpr int type_algebre(erreur<T>& element) {
	return type_algebre(T());
};



//0 : exact. 1 : approché
constexpr int type_approx(int element) {
	return 0;
};

constexpr int type_approx(long element) {
	return 0;
};

constexpr int type_approx(long long element) {
	return 0;
};

constexpr int type_approx(InfInt& element) {
	return 0;
};

constexpr int type_approx(int_precision& element) {
	return 0;
};

constexpr int type_approx(float element) {
	return 1;
};

constexpr int type_approx(double element) {
	return 1;
};

template<class T> constexpr int type_approx(float_precision& element) {
	return 1;
};

// composés
template<class T> constexpr int type_approx(erreur<T>& element) {
	return type_approx(T());
};

template<class T> constexpr int type_approx(polynome<T>& element) {
	return type_approx(T());
};

template<class T> constexpr int type_approx(matrice<T>& element) {
	return type_approx(T());
};

template<class T> constexpr int type_approx(rationnel<T>& element) {
	return type_approx(T());
};


template<class T> constexpr int type_approx(anneau_quotient<T>& element) {
	return type_approx(T());
};



template<class T> constexpr int type_approx(corps_quotient<T>& element) {
	return type_approx(T());
};


template<class T> constexpr int type_approx(polynome_n<T>& element) {
	return type_approx(T());
};

template<class T> constexpr int type_approx(complexe<T>& element) {
	return type_approx(T());
};


*/