#pragma once

#include "fonctions template.hpp"
#include "entete objets.hpp"
#include "polynome_n_rec.hpp"

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

//définition
template<class T, typename Enable=void> class norme_T {
public:
};

//classes de base :
template<> class norme_T<int> {
public:
	static int norme(int x) {
		return abs(x);
	};
};

template<> class norme_T<long> {
public:
	static long norme(long x) {
		return abs(x);
	};
};

template<> class norme_T<long long > {
public:
	static long long norme(long long x) {
		return abs(x);
	};
};

template<> class norme_T<int_precision> {
public:

	static int_precision norme(int_precision const& x) {
		return abs(x);
	};
};

template<> class norme_T<InfInt> {
public:
	static InfInt norme(InfInt const& x) {
		return abs(x);
	};
};

template<> class norme_T<float> {
public:
	static float norme(float x) {
		return abs(x);
	};
};

template<> class norme_T<double> {
public:
	static float norme(double x) {
		return abs((float)x);
	};
};

template<> class norme_T<float_precision> {
public:
	static float norme(float_precision const& x) {
		return abs((float) x);
	};
};

template<class T> class norme_T<erreur_b<T>> {
public:
	static float norme(erreur_b<T> const& temp) {
		if ((bool) temp)
			return ((float) temp.valeur);
		else
			return 0.;
	};
};

template<class T> class norme_T<erreur_l<T>> {
public:
	static float norme(erreur_l<T> const& temp) {
		if ((bool)temp)
			return (float)temp.valeur;
		else
			return 0.;
	};
};


//class composées
template<class T> class norme_T<complexe<T>> {
public:
	static decltype(norme_T<T>::norme(T())) norme(complexe<T> const& temp) {
		return norme_T<T>::norme(temp.x) + norme_T<T>::norme(temp.y);
	};
};

template<class T> class norme_T<std::complex<T>> {
public:
	static decltype(norme_T<T>::norme(T())) norme(std::complex<T> const& temp) {
		return norme_T<T>::norme(temp.real) + norme_T<T>::norme(temp.imag);
	};
};

template<class T> class norme_T< polynome<T> > {
public:
	static decltype(norme_T<T>::norme(T())) norme(polynome<T> const& poly) {
		decltype(norme_T<T>::norme(T())) x = norme_T<T>::norme(poly.coeffs[0]);
		for (int i(1); i < poly.coeffs.size(); ++i) {
			x = x + norme_T<T>::norme(poly.coeffs[i]);
		}
		return x;
	};
};

/*
template<class T> class norme_T< polynome_n_rec<T> > {
public:
	static decltype(norme_T<T>::norme(T())) norme(polynome_n_rec<T> poly) {
		if (poly.n_var >= 1) {
			decltype(norme_T<T>::norme(T())) x = norme_T<T>::norme(unite(poly.element,false));
//			x = unite(x,false);
			for (int i(0); i < poly.poly.coeffs.size(); ++i) {
				x = x + norme_T<polynome_n_rec<T>>::norme(poly.poly.coeffs[i]);
			}
			return x;
		};
		return norme_T<T>::norme(poly.element);
	};
};*/

template<class T> class norme_T< polynome_n_rec<T> > {
public:
	static decltype(norme_T<T>::norme(T())) norme(polynome_n_rec<T> const& poly) {
		if (poly.n_var >= 1) {
			decltype(norme_T<T>::norme(T())) x = norme_T<T>::norme(unite(poly.element, false));
			//			x = unite(x,false);
			for (auto it=poly.cbegin(); (bool) it; ++it) {
				x  += norme_T<T>::norme(*it);
			}
			return x;
		};
		return norme_T<T>::norme(poly.element);
	};
};


template<class T> class norme_T< polynome_n_iter<T> > {
public:
	static decltype(norme_T<T>::norme(T())) norme(polynome_n_iter<T> const& poly) {
		decltype(norme_T<T>::norme(T())) x = norme_T<T>::norme(poly.coeffs.data[0]);
		for (int i(1); i < poly.coeffs.data.size(); ++i) {
			x = x + norme_T<T>::norme(poly.coeffs.data[i]);
		}
		return x;
	};
};


template<class T, int n> class norme_T< polynome_n_fixe<T, n> > {
public:
	static decltype(norme_T<T>::norme(T())) norme(polynome_n_fixe<T, n> const& poly) {
		if constexpr (n >= 1) {
			decltype(norme_T<T>::norme(T())) x = norme_T<T>::norme(unite(poly.get_T(), false));
			//			x = unite(x,false);
			for (auto it = poly.cbegin(); (bool)it; ++it) {
				x += norme_T<T>::norme(*it);
			}
			return x;
		}
		else
			return norme_T<T>::norme(poly.element);
	};
};

template<class T> class norme_T< polynome_n_fixe<T, 0> > {
public:
	static decltype(norme_T<T>::norme(T())) norme(polynome_n_fixe<T, 0> const& poly) {
			decltype(norme_T<T>::norme(T())) x = norme_T<T>::norme(poly.element);
			return x;
		};
};




//deux cas : si norme(T) a une division exacte ou non.
template<class T> class norme_T < rationnel<T>, typename std::enable_if< type_algebre<decltype(norme_T<T>::norme(T()))>::type == 0 > ::type > {
public:
	static decltype(norme_T<T>::norme(T())) norme(rationnel<T> const& frac) {
//		decltype(norme_T<T>::norme()) x = norme_T<T>::norme(frac.numerateur) / norme_T<T>::norme(frac.denominateur);
		return norme_T<T>::norme(frac.numerateur) / norme_T<T>::norme(frac.denominateur);
	};

};

template<class T> class norme_T < rationnel<T>, typename std::enable_if< type_algebre<decltype(norme_T<T>::norme(T()))>::type != 0 > ::type > {
public:
	static rationnel<decltype(norme_T<T>::norme(T()))> norme(rationnel<T> const& frac) {
//		rationnel<decltype(norme_T<T>::norme())> x = rationnel<decltype(norme_T<T>::norme()) >> (norme_T<T>::norme(frac.numerateur), norme_T<T>::norme(frac.denominateur));
		auto x= rationnel<decltype(norme_T<T>::norme(T()))>(norme_T<T>::norme(frac.numerateur), norme_T<T>::norme(frac.denominateur));
		return x;
	};

};

template<class T> class norme_T < anneau_quotient<polynome<T>>> {
public:
	static decltype(norme_T<T>::norme(T())) norme(anneau_quotient<polynome<T>> const& x) {
		//on suppose x normalisé : coeff_dominant =1 
		polynome<T> poly = x.element;
		T coeff = poly.coeff_dominant();
		poly.normaliser();
		polynome<T> pgcd = PGCD(poly, x.quotient); //dans cet ordre ... résultant proportionnel au coeff dominant de x.element.
		if (pgcd.degre >= 1) {
			return norme_T<T>::norme(unite(pgcd.coeffs[0], false));
		}
		return norme_T<polynome<T>>::norme(coeff * pgcd);
	};
};

template<class T> class norme_T < corps_quotient<polynome<T>>> {
public:
	static decltype(norme_T<T>::norme(T())) norme(corps_quotient<polynome<T>> const& x) {
		polynome<T> poly = x.element;
		T coeff = poly.coeff_dominant();
		poly.normaliser();
		polynome<T> pgcd = PGCD(poly, x.quotient); //dans cet ordre ... résultant proportionnel au coeff dominant de x.element.
		if (pgcd.degre >= 1) {
			return norme_T<T>::norme(unite(pgcd.coeffs[0], false));
		}
		return norme_T<polynome<T>>::norme(coeff * pgcd);
	};
};


//pour aller plus vite : la fonction directement !
template<class T> inline decltype(norme_T<T>::norme(T())) norme(T const& x) {
	return norme_T<T>::norme(x);
};

