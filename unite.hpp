#pragma once

#include "InfInt.h"
#include "precision/fprecision.h"
#include "precision/iprecision.h"
#include "erreur.hpp"

#include "fonctions template.hpp" //abs
#include "types.hpp" //pour rationnel
#include "entete objets.hpp"

template<class T> class erreur;
template<class T> class anneau_quotient;
template<typename T> class complexe;
template<class T> class corps_quotient;
template<class T, class enable1 = void, class enable2 = void> class matrice;
template<typename T> class polynome;
template<class T> class polynome_n;
template<class T> class polynome_n_iter;
template<class T, class enable = void> class rationnel;



template<class T> T unite(T const& element,bool test) {
	T vrai_ = element;
	vrai_ = test;
	return vrai_;
};

int unite(int const& element,bool test) {
	return test ? 1:0;
}

long unite(long const& element,bool test) {
	return test ? 1 : 0;
}

long long unite(long long const& element,bool test) {
	return test ? 1 : 0;
}

float unite(float const& element,bool test) {
	return test ? 1 : 0;
}

double unite(double const& element,bool test) {
	return test ? 1 : 0;
}

InfInt unite(InfInt const& element, bool test) {
	return test ? InfInt(1):InfInt(0);
}

int_precision unite(int_precision const& element, bool test) {
	return test ? int_precision(1) : int_precision(0);
}

float_precision unite(float_precision const& element, bool test) {
	return float_precision(test ? 1 : 0, element.precision(), element.mode());
}

template<class T> erreur<T> unite(erreur<T> const& temp, bool test) {
	if (test)
		return erreur<T>(unite(temp.valeur, test), precision_relative(temp.valeur));
	else
		return erreur<T>(unite(temp.valeur, test), 0.);
}

template<class T> polynome<T> unite(polynome<T> const& poly, bool test) {
	polynome<T> temp(unite(poly.coeffs[0],test));
	return temp;
}

template<class T> polynome_n<T> unite(polynome_n<T> const& poly_n, bool test) {
	T temp = unite(poly_n.element,test);
	polynome_n<T> temp_(poly_n.n_var, temp, poly_n.noms_variables);
	return temp_;
}

template<class T> polynome_n_iter<T> unite(polynome_n_iter<T> const& poly_n, bool test) {
	polynome_n_iter<T> temp(poly_n.coeffs.puissance, unite(poly_n.coeffs.data[0],test), poly_n.noms);
	return temp;
}

template<class T> anneau_quotient<T> unite(anneau_quotient<T> const& temp, bool test) {
	return anneau_quotient<T>(unite(temp.element,test), temp.quotient);
}

template<class T> corps_quotient<T> unite(corps_quotient<T> const& temp, bool test) {
	return corps_quotient<T>(unite(temp.element,test), temp.quotient);
}

template<class T> complexe<T> unite(complexe<T> const& c, bool test) {
	return complexe<T>(unite(c.x,test), unite(c.x,false));
}

template<class T> rationnel<T> unite(rationnel<T> const& temp, bool test) {
	return rationnel<T>(unite(temp.numerateur,test), unite(temp.numerateur,true));
}

template<class T> matrice<T> unite(matrice<T> const& m, bool test) {
	matrice<T> result(m.taille_l, m.taille_c, unite(m.coeffs[0][0],false));
	if (test) {
		T vrai_ = unite(m.coeffs[0][0], true);
		for (int i(0); i < min(m.taille_l, m.taille_c); ++i)
			result.coeffs[i][i] = vrai_;
	}
	return result;
}

/*
// >>>>>>>>>>>>>>>>>>>>   FAUX   <<<<<<<<<<<<<<<<<<<<

template<class T> T faux(T const& temp) {
	T faux_ = temp;
	faux_ = false;
	return faux_;
};

int faux(int const& element) {
	return 0;
}

long faux(long const& element) {
	return 0;
}

long long faux(long long const& element) {
	return 0;
}

float faux(float const& element) {
	return 0.;
}

double faux(double const& element) {
	return 0;
}

InfInt faux(InfInt const& element) {
	return InfInt(0);
}

int_precision faux(int_precision const& element) {
	return int_precision(0);
}

float_precision faux(float_precision const& element) {
	return float_precision(0, element.precision(), element.mode());
}

template<class T> erreur<T> faux(erreur<T> const& temp) {
	return erreur<T>(faux(temp.valeur), precision_relative(temp.valeur));
}


template<class T> polynome<T> faux(polynome<T> const& poly) {
	polynome<T> faux_(faux(poly.coeffs[0]));
	return faux_;
}

template<class T> polynome_n<T> faux(polynome_n<T> const& poly_n) {
	T faux_T = faux(poly_n.element);
	polynome_n<T> faux_(poly_n.n_var, faux_T, poly_n.noms_variables);
	return faux_;
}

template<class T> polynome_n_iter<T> faux(polynome_n_iter<T> const& poly_n) {
	polynome_n_iter<T> faux_(poly_n.coeffs.puissance, faux(poly_n.coeffs.data[0]), poly_n.noms);
	return faux_;
}

template<class T> anneau_quotient<T> faux(anneau_quotient<T> const& temp) {
	return anneau_quotient<T>(faux(temp.element), temp.quotient);
}

template<class T> corps_quotient<T> faux(corps_quotient<T> const& temp) {
	return corps_quotient<T>(faux(temp.element), temp.quotient);
}

template<class T> complexe<T> faux(complexe<T> const& c) {
	return complexe<T>(faux(c.x), faux(c.x));
}

template<class T> rationnel<T> faux(rationnel<T> const& temp) {
	return rationnel<T>(faux(temp.numerateur), vrai(temp.numerateur));
}

template<class T> matrice<T> faux(matrice<T> const& temp) {
	matrice<T> result(temp.taille_l, temp.taille_c, faux(temp.coeffs[0][0]));
	return result;
}

*/