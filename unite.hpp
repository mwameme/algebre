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
template<class T> class scalaire_vecteur;


template<class T> inline T unite(T const& element,bool test) {
	T vrai_ = element;
	vrai_ = test;
	return vrai_;
};

inline int unite(int const& element,bool test) {
	return test ? 1:0;
}

inline long unite(long const& element,bool test) {
	return test ? 1 : 0;
}

inline long long unite(long long const& element,bool test) {
	return test ? 1 : 0;
}

inline float unite(float const& element,bool test) {
	return test ? 1 : 0;
}

inline double unite(double const& element,bool test) {
	return test ? 1 : 0;
}

inline InfInt unite(InfInt const& element, bool test) {
	return test ? InfInt(1):InfInt(0);
}

inline int_precision unite(int_precision const& element, bool test) {
	return test ? int_precision(1) : int_precision(0);
}

inline float_precision unite(float_precision const& element, bool test) {
	return float_precision(test ? 1 : 0, element.precision(), element.mode());
}

template<class T> inline erreur<T> unite(erreur<T> const& temp, bool test) {
	if (test)
		return erreur<T>(unite(temp.valeur, test), precision_relative(temp.valeur));
	else
		return erreur<T>(unite(temp.valeur, test), 0.);
}

template<class T> inline polynome<T> unite(polynome<T> const& poly, bool test) {
	polynome<T> temp(unite(poly.coeffs[0],test));
	return temp;
}

template<class T> inline polynome_n<T> unite(polynome_n<T> const& poly_n, bool test) {
	T temp = unite(poly_n.element,test);
	polynome_n<T> temp_(poly_n.n_var, temp, poly_n.noms_variables);
	return temp_;
}

template<class T> inline polynome_n_iter<T> unite(polynome_n_iter<T> const& poly_n, bool test) {
	polynome_n_iter<T> temp(poly_n.coeffs.puissance, unite(poly_n.coeffs.data[0],test), poly_n.noms);
	return temp;
}

template<class T> inline anneau_quotient<T> unite(anneau_quotient<T> const& temp, bool test) {
	return anneau_quotient<T>(unite(temp.element,test), temp.quotient);
}

template<class T> inline corps_quotient<T> unite(corps_quotient<T> const& temp, bool test) {
	return corps_quotient<T>(unite(temp.element,test), temp.quotient);
}

template<class T> inline complexe<T> unite(complexe<T> const& c, bool test) {
	return complexe<T>(unite(c.x,test), unite(c.x,false));
}

template<class T> inline rationnel<T> unite(rationnel<T> const& temp, bool test) {
	return rationnel<T>(unite(temp.numerateur,test), unite(temp.numerateur,true));
}

template<class T> inline matrice<T> unite(matrice<T> const& m, bool test) {
	matrice<T> result(m.taille_l, m.taille_c, unite(m.coeffs[0][0],false));
	if (test) {
		T vrai_ = unite(m.coeffs[0][0], true);
		for (int i(0); i < min(m.taille_l, m.taille_c); ++i)
			result.coeffs[i][i] = vrai_;
	}
	return result;
}

template<class T> inline scalaire_vecteur<T> unite(scalaire_vecteur<T> const& temp, bool test) {
	T element = unite( temp.scalaire,test);
	scalaire_vecteur<T> retour(temp.vecteur.size());
	retour = element;
	return retour;
}