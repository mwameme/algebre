#pragma once


#include "InfInt.h"
#include <complex>
#include "precision/fprecision.h"
#include "precision/iprecision.h"

#include "fonctions template.hpp" //abs
#include "types.hpp" //pour rationnel
#include "entete objets.hpp"

template<class T> class complex;

template<class T> class erreur;
template<class T> class anneau_quotient;
template<typename T> class complexe;
template<class T> class corps_quotient;
template<class T, class enable1 = void, class enable2 = void> class matrice;
template<typename T> class polynome;
template<class T> class polynome_n;
template<class T> class polynome_n_iter;
template<class T, class enable = void> class rationnel;



template<class T> class type_algebre;


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

	static int_precision norme(int_precision x) {
		return abs(x);
	};
};

template<> class norme_T<InfInt> {
public:
	static InfInt norme(InfInt x) {
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
	static double norme(double x) {
		return abs(x);
	};
};

template<> class norme_T<float_precision> {
public:
	static float_precision norme(float_precision x) {
		return abs(x);
	};
};

template<class T> class norme_T<erreur<T>> {
public:
	static float norme(erreur<T> temp) {
		return norme_T<float>::norme((float) temp.valeur);
	};
};

//class composées
template<class T> class norme_T<complexe<T>> {
public:
	static decltype(norme_T<T>::norme(T())) norme(complexe<T> temp) {
		return norme_T<T>::norme(temp.x) + norme_T<T>::norme(temp.y);
	};
};

template<class T> class norme_T<complex<T>> {
public:
	static decltype(norme_T<T>::norme(T())) norme(complex<T> temp) {
		return norme_T<T>::norme(temp.real()) + norme_T<T>::norme(temp.imag());
	};
};

template<class T> class norme_T< polynome<T> > {
public:
	static decltype(norme_T<T>::norme(T())) norme(polynome<T> poly) {
		decltype(norme_T<T>::norme(T())) x = norme_T<T>::norme(poly.coeffs[0]);
		for (int i(1); i < poly.coeffs.size(); ++i) {
			x = x + norme_T<T>::norme(poly.coeffs[i]);
		}
		return x;
	};
};

template<class T> class norme_T< polynome_n<T> > {
public:
	static decltype(norme_T<T>::norme(T())) norme(polynome_n<T> poly) {
		if (poly.n_var >= 1) {
			decltype(norme_T<T>::norme(T())) x = unite(norme_T<T>::norme(poly.element),false);
//			x = unite(x,false);
			for (int i(0); i < poly.coeffs.size(); ++i) {
				x = x + norme_T<polynome_n<T>>::norme(* poly.coeffs[i]);
			}
			return x;
		};
		return norme_T<T>::norme(poly.element);
	};
};

template<class T> class norme_T< polynome_n_iter<T> > {
public:
	static decltype(norme_T<T>::norme(T())) norme(polynome_n_iter<T> poly) {
		decltype(norme_T<T>::norme(T())) x = norme_T<T>::norme(poly.coeffs.data[0]);
		for (int i(1); i < poly.coeffs.data.size(); ++i) {
			x = x + norme_T<T>::norme(poly.coeffs.data[i]);
		}
		return x;
	};
};



//deux cas : si norme(T) a une division exacte ou non.
template<class T> class norme_T < rationnel<T>, typename std::enable_if< type_algebre<decltype(norme_T<T>::norme(T()))>::type == 0 > ::type > {
public:
	static decltype(norme_T<T>::norme(T())) norme(rationnel<T> frac) {
//		decltype(norme_T<T>::norme()) x = norme_T<T>::norme(frac.numerateur) / norme_T<T>::norme(frac.denominateur);
		return norme_T<T>::norme(frac.numerateur) / norme_T<T>::norme(frac.denominateur);
	};

};

template<class T> class norme_T < rationnel<T>, typename std::enable_if< type_algebre<decltype(norme_T<T>::norme(T()))>::type != 0 > ::type > {
public:
	static rationnel<decltype(norme_T<T>::norme(T()))> norme(rationnel<T> frac) {
//		rationnel<decltype(norme_T<T>::norme())> x = rationnel<decltype(norme_T<T>::norme()) >> (norme_T<T>::norme(frac.numerateur), norme_T<T>::norme(frac.denominateur));
		auto x= rationnel<decltype(norme_T<T>::norme(T()))>(norme_T<T>::norme(frac.numerateur), norme_T<T>::norme(frac.denominateur));
		return x;
	};

};

template<class T> class norme_T < anneau_quotient<polynome<T>>> {
public:
	static decltype(norme_T<T>::norme(T())) norme(anneau_quotient<polynome<T>> x) {
		polynome<T> pgcd = PGCD(x.element, x.quotient); //dans cet ordre ... résultant proportionnel au coeff dominant de x.element.
		if (pgcd.degre >= 1) {
			T temp = unite(pgcd.coeffs[0],false);
			return norme_T<T>::norme(temp);
		}
		return norme_T<polynome<T>>::norme(pgcd);
	};
};


//pour aller plus vite : la fonction directement !
template<class T> inline decltype(norme_T<T>::norme(T())) norme(T x) {
	return norme_T<T>::norme(x);
};

