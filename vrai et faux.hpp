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



template<class T> T vrai(T const& element) {
	T vrai_ = element;
	vrai_ = true;
	return vrai_;
};

int vrai(int const& element) {
	return 1;
}

long vrai(long const& element) {
	return 1;
}

long long vrai(long long const& element) {
	return 1;
}

float vrai(float const& element) {
	return 1.;
}

double vrai(double const& element) {
	return 1;
}

InfInt vrai(InfInt const& element) {
	return InfInt(1);
}

int_precision vrai(int_precision const& element) {
	return int_precision(1);
}

float_precision vrai(float_precision const& element) {
	return float_precision(1, element.precision(), element.mode());
}

template<class T> erreur<T> vrai(erreur<T> const& temp) {
	return erreur<T>(vrai(temp.valeur), precision_relative(temp.valeur));
}

template<class T> polynome<T> vrai(polynome<T> const& poly) {
	polynome<T> vrai_(vrai(poly.coeffs[0]));
	return vrai_;
}

template<class T> polynome_n<T> vrai(polynome_n<T> const& poly_n) {
	T vrai_T = vrai(poly_n.element);
	polynome_n<T> vrai_(poly_n.n_var, vrai_T, poly_n.noms_variables);
	return vrai_;
}

template<class T> polynome_n_iter<T> vrai(polynome_n_iter<T> const& poly_n) {
	polynome_n_iter<T> vrai_(poly_n.coeffs.puissance, vrai(poly_n.coeffs.data[0]), poly_n.noms);
	return vrai_;
}

template<class T> anneau_quotient<T> vrai(anneau_quotient<T> const& temp) {
	return anneau_quotient<T>(vrai(temp.element), temp.quotient);
}

template<class T> corps_quotient<T> vrai(corps_quotient<T> const& temp) {
	return corps_quotient<T>(vrai(temp.element), temp.quotient);
}

template<class T> complexe<T> vrai(complexe<T> const& c) {
	return complexe<T>(vrai(c.x), faux(c.x));
}

template<class T> rationnel<T> vrai(rationnel<T> const& temp) {
	return rationnel<T>(vrai(temp.numerateur), vrai(temp.numerateur));
}

template<class T> matrice<T> vrai(matrice<T> const& temp) {
	matrice<T> result(temp.taille_l, temp.taille_c, faux(temp.coeffs[0][0]));
	T vrai_ = vrai(temp.coeffs[0][0]);
	for (int i(0); i < min(temp.taille_l, temp.taille_c); ++i)
		result.coeffs[i][i] = vrai_;

	return result;
}

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
