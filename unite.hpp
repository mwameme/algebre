#pragma once


#include "objets.hpp"
#include "entete objets.hpp"
#include "InfInt.h"
#include "precision/iprecision.h"
#include "precision/fprecision.h"
#include "simplifier polynome_n.hpp"

#include "erreur.hpp"
#include <complex>

#include "noms.hpp"

template<class T> class complex;

template<class T> class erreur_b;
template<class T> class erreur_l;
template<class T> class anneau_quotient;
template<class T> class complexe;
template<class T> class corps_quotient;
template<class T, class enable1 = void, class enable2 = void> class matrice;
template<class T> class polynome;
template<class T> class polynome_n_rec;
template<class T> class polynome_n_iter;
template<class T> class rationnel;
template<class T, int n> class polynome_n_fixe;

class InfInt;
class int_precision;
class float_precision;

template<class T> class scalaire_vecteur;

template<class T> class polynome_n_sparse;
template<class T> class monome;


template<class T>  T unite(T const& element,bool test) {
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

template<class T>  erreur_b<T> unite(erreur_b<T> const& temp, bool test) {
	if (test)
		return erreur_b<T>(unite(temp.valeur, test), precision_relative(temp.valeur));
	else
		return erreur_b<T>(unite(temp.valeur, test), 0.);
}

template<class T>  erreur_l<T> unite(erreur_l<T> const& temp, bool test) {
	if (test)
		return erreur_l<T>(unite(temp.valeur, test), precision_relative_l(temp.valeur));
	else
		return erreur_l<T>(unite(temp.valeur, test), 0.);
}


template<class T>  polynome<T> unite(polynome<T> const& poly, bool test) {
	polynome<T> temp(unite(poly.coeffs[0],test));
	return temp;
}

template<class T>  polynome_n_rec<T> unite(polynome_n_rec<T> const& poly_n, bool test) {
	T temp = unite(poly_n.element,test);
	polynome_n_rec<T> temp_(poly_n.n_var, temp);
	return temp_;
}

template<class T>  polynome_n_iter<T> unite(polynome_n_iter<T> const& poly_n, bool test) {
	polynome_n_iter<T> temp(poly_n.coeffs.puissance, unite(poly_n.coeffs.data[0],test));
	return temp;
}

template<class T>  anneau_quotient<T> unite(anneau_quotient<T> const& temp, bool test) {
	return anneau_quotient<T>(unite(temp.element,test), temp.quotient);
}

template<class T>  corps_quotient<T> unite(corps_quotient<T> const& temp, bool test) {
	return corps_quotient<T>(unite(temp.element,test), temp.quotient);
}

template<class T>  complexe<T> unite(complexe<T> const& c, bool test) {
	return complexe<T>(unite(c.x,test), unite(c.x,false));
}

template<class T>  std::complex<T> unite(std::complex<T> const& c, bool test) {
	return std::complex<T>(unite(c.real(), test), unite(c.real(), false));
};


template<class T>  rationnel<T> unite(rationnel<T> const& temp, bool test) {
	return rationnel<T>(unite(temp.numerateur, test), unite(temp.numerateur, true));
};

template<class T>  matrice<T> unite(matrice<T> const& m, bool test) {
	matrice<T> result(m.taille_l, m.taille_c, unite(m.coeffs[0][0], false));
	if (test) {
		T vrai_ = unite(m.coeffs[0][0], true);
		for (int i(0); i < min(m.taille_l, m.taille_c); ++i)
			result.coeffs[i][i] = vrai_;
	}
	return result;
};

template<class T>  scalaire_vecteur<T> unite(scalaire_vecteur<T> const& temp, bool test) {
	T element = unite(temp.scalaire, test);
	scalaire_vecteur<T> retour(temp.vecteur.size());
	retour = element;
	return retour;
};

template<class T, int n> polynome_n_fixe<T, n> unite(polynome_n_fixe<T, n> const& temp, bool test) {
	T element = temp.get_T();
	element = unite(element, test);
	return polynome_n_fixe<T, n>(element);
};

template<class T, int n> polynome_n_fixe<T, 1> unite(polynome_n_fixe<T, 1> const& temp, bool test) {
	T element = temp.get_T();
	element = unite(element, test);
	return polynome_n_fixe<T, 1>(element);
};

template<class T, int n> polynome_n_fixe<T, 0> unite(polynome_n_fixe<T, 0> const& temp, bool test) {
	T element = temp.get_T();
	element = unite(element, test);
	return polynome_n_fixe<T, 0>(element);
};

template<class T>
monome<T> unite(monome<T> const& temp, bool test) {
	monome<T> local(temp);
	for (int i(0); i < local.degres.size(); ++i)
		local.degres[i] = 0;
	local.element = unite(local.element, test);
	return local;
};

template<class T> 
polynome_n_sparse<T> unite(polynome_n_sparse<T> const& temp, bool test) {
	monome<T> local = unite(temp.monomes[0], test);
	return polynome_n_sparse<T>(local);
};