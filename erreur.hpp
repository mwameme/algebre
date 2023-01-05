#pragma once
#include "fonctions template.hpp"
#include "polynome.hpp"
#include "rationnel.hpp"
#include "anneau quotient.hpp"
#include <cmath>
#include <iostream>
#include "precision/fprecision.h"
#include "types.hpp"


inline float precision_relative(double x) {
	return carre(1.0e-12);
};

inline float precision_relative(float x) {
	return carre(1.0e-6);
};

/*
float precision_relative(float_precision const& x) {
	int power = x.precision();
	float res = puissance((float) .1, 2*power);
	return res;
};*/

//test de précision ... Essayer avec 1. et 10 000. (par exemple)
inline float precision_relative(float_precision const& x) {
	return carre((float) x.epsilon());
};




template<class T> class erreur {
public :
	static_assert((type_algebre<T>::type == 0) && (type_algebre<T>::approx == 1), "erreur<T> : T doit être approché, avec division exacte");

	erreur() :precision(0) {
		valeur = false;
		precision = precision_relative(T(1.));
	};


	explicit erreur(T const& valeur_) {
		valeur = valeur_;
		precision = carre((float)valeur) * precision_relative(valeur);
	};

	/*
	erreur(bool const& test) {
		valeur = test;
		precision = precision_relative(valeur);
	};
	*/


	erreur(T const& valeur_, float const& precision_) {
		valeur = valeur_;
		precision = precision_;
		//		maj();
	};

	erreur(erreur<T> const& temp) {
		valeur = temp.valeur;
		precision = temp.precision;
	};

	operator bool() const {
		if (carre((float)valeur) <= (50. * precision))
			return false;
		else
			return true;
	};



	friend erreur<T> operator*(erreur<T> const& temp1, erreur<T> const& temp2) {
		erreur<T> result(temp1.valeur * temp2.valeur, 0.);
		result.precision = carre((float)temp1.valeur) * temp2.precision  +  carre((float)temp2.valeur) * temp1.precision; // +(temp1.precision * temp2.precision);
		return result;
	};

	//template<class U> 
	friend erreur<T> operator*(long const& scalaire, erreur<T> const& temp) {
		erreur<T> result(scalaire * temp.valeur, ((float)carre(scalaire)) * temp.precision); //((T)scalaire)
		return result;
	};

	friend erreur<T> operator/(erreur<T> const& temp1, erreur<T> const& temp2) {
		/*
		erreur<T> result(temp1.valeur / temp2.valeur);
		result.precision = temp1.precision/abs(temp2.valeur) + ( temp2.precision *abs(temp1.valeur/(temp2.valeur*temp2.valeur)) ) + ( abs(result.valeur)* precision_relative(result.valeur));
		*/
		erreur<T> result(temp1.valeur / temp2.valeur, 0.);
		result.precision = temp1.precision / (carre((float)temp2.valeur))  +(temp2.precision * carre((float)temp1.valeur / ((float)temp2.valeur * (float)temp2.valeur)));
	//		result.maj();
		return result;
	};


	//A MODIFIER
	/*
	friend erreur<T> operator+(erreur<T> const& temp1, erreur<T> const& temp2) {
		erreur<T> result(temp1.valeur + temp2.valeur, 0.);
		if (abs(temp1.valeur) > abs(temp2.valeur)) {
			if (abs((float)temp2.valeur) < abs((float)temp1.valeur) * precision_relative(temp1.valeur))
				result.precision = temp1.precision + carre((float)temp2.valeur);
			else
				result.precision = temp1.precision + temp2.precision + carre((float)temp1.valeur) * precision_relative(temp1.valeur);
		}
		else {
			if (abs((float)temp1.valeur) < abs((float)temp2.valeur) * precision_relative(temp2.valeur))
				result.precision = temp2.precision + carre((float)temp1.valeur);
			else
				result.precision = temp1.precision + temp2.precision + carre((float)temp2.valeur) * precision_relative(temp2.valeur);
		}
		//		result.maj();
		return result;
	};*/

	//A MODIFIER
	friend erreur<T> operator+(erreur<T> const& temp1, erreur<T> const& temp2) {
		erreur<T> result(temp1.valeur + temp2.valeur, 0.);
		if (carre((float)temp2.valeur) < carre((float)temp1.valeur) * precision_relative(temp1.valeur)) {
			result.precision = temp1.precision + carre((float)temp2.valeur);
			return result;
		}
		if (carre((float)temp1.valeur) < carre((float)temp2.valeur) * precision_relative(temp2.valeur)) {
			result.precision = temp2.precision + carre((float)temp1.valeur);
			return result;
		}

		float erreur1 = carre((float)temp1.valeur) * precision_relative(temp1.valeur); //erreur sur le calcul de x+y liée à x. registre
		float erreur2 = carre((float)temp2.valeur) * precision_relative(temp2.valeur); //erreur sur le calcul de x+y liée à y. registre
		if (erreur1 > erreur2) {
			if ((float)temp2.precision < erreur1) {
				result.precision = temp1.precision + erreur1;
				return result;
			}
			result.precision = temp1.precision + temp2.precision;
			return result;
		}
		if ((float)temp1.precision < erreur2) {
			result.precision = temp2.precision + erreur2;
			return result;
		}
		result.precision = temp2.precision + temp1.precision;
		return result;
	};

	//A MODIFIER
	friend erreur<T> operator-(erreur<T> const& temp1, erreur<T> const& temp2) {
		erreur<T> result(temp1.valeur - temp2.valeur,0.);
		//copié de operator+
		if (carre((float)temp2.valeur) < carre((float)temp1.valeur) * precision_relative(temp1.valeur)) {
			result.precision = temp1.precision + carre((float)temp2.valeur);
			return result;
		}
		if (carre((float)temp1.valeur) < carre((float)temp2.valeur) * precision_relative(temp2.valeur)) {
			result.precision = temp2.precision + carre((float)temp1.valeur);
			return result;
		}

		float erreur1 = carre((float)temp1.valeur) * precision_relative(temp1.valeur); //erreur sur le calcul de x+y liée à x
		float erreur2 = carre((float)temp2.valeur) * precision_relative(temp2.valeur); //erreur sur le calcul de x+y liée à y
		if (erreur1 > erreur2) {
			//		float x = carre((float)temp1.valeur) * precision_relative(temp1.valeur);
			if ((float)temp2.precision < erreur1) {
				result.precision = temp1.precision + erreur1;
				return result;
			}
			result.precision = temp1.precision + temp2.precision;
			return result;
		}
		//	float x = carre((float)temp2.valeur) * precision_relative(temp2.valeur);
		if ((float)temp1.precision < erreur2) {
			result.precision = temp2.precision + erreur2;
			return result;
		}
		result.precision = temp2.precision + temp1.precision;
		return result;
	};

	friend erreur<T> operator-(erreur<T> const& temp) {
		erreur<T> result(-temp.valeur,temp.precision);
		return result;
	};


	erreur<T>& operator=(erreur<T> const& temp) {
		valeur = temp.valeur;
		precision = temp.precision;
		return *this;
	};

	erreur<T>& operator=(bool const& test) {
		precision = precision_relative(valeur);
		valeur = test;
		return *this;
	};


	friend bool operator==(erreur<T> temp1, erreur<T> temp2) {
		return !((bool)(temp1 - temp2));
	};

	friend bool operator< (erreur<T> temp1, erreur<T> temp2) {
		return (temp1.valeur < temp2.valeur);
	};

	friend bool operator> (erreur<T> temp1, erreur<T> temp2) {
		return (temp1.valeur > temp2.valeur);
	};

	friend std::ostream& operator<<(std::ostream& os, const erreur<T>& element) {
		if ((bool) element)
			os << (double)element.valeur;
		else
			os << "0";
		return os;
	};

	operator T() {
		return valeur;
	}


	T valeur;
	float precision;
};