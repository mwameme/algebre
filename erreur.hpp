#pragma once
#include "fonctions template.hpp"
#include "polynome.hpp"
#include "rationnel.hpp"
#include "anneau quotient.hpp"
#include <cmath>
#include <iostream>
#include "precision/fprecision.h"
#include "types.hpp"

#include "swap_T.hpp"

inline float precision_relative(double x) {
	return carre(1.0e-15);
};

inline float precision_relative(float x) {
	return carre(1.0e-6);
};


//test de précision ... Essayer avec 1. et 10 000. (par exemple)
inline float precision_relative(float_precision const& x) {
	return carre((float) x.epsilon());
};


template<class T> class erreur_b {
public :
	static_assert((type_algebre<T>::type == 0) && (type_algebre<T>::approx == 1), "erreur_b<T> : T doit être approché, avec division exacte");

	T valeur;
	float precision;


	erreur_b() :precision(0) {
		valeur = 0;
		precision = 0.;
	};


	explicit erreur_b(T const& valeur_) {
		valeur = valeur_;
		precision = carre((float)valeur) * precision_relative(valeur);
	};


	erreur_b(T const& valeur_, float const& precision_) { //faire attention au carré
		valeur = valeur_;
		precision = precision_;
		//		maj();
	};

	erreur_b(erreur_b<T> const& temp) {
		valeur = temp.valeur;
		precision = temp.precision;
	};

	erreur_b(erreur_b<T>&& temp) {
		swap(*this, temp);
		return;
	}


	erreur_b<T>& operator=(erreur_b<T> const& temp) {
		if (this == &temp)
			return *this;
		valeur = temp.valeur;
		precision = temp.precision;
		return *this;
	};

	erreur_b<T>& operator=(erreur_b<T>&& temp) {
		if (this == &temp)
			return *this;
		swap(*this, temp);
		return *this;
	}

	erreur_b<T>& operator+=(erreur_b<T> const& autre) {
		return (*this = (*this + autre));
	}

	erreur_b<T>& operator*=(erreur_b<T> const& autre) {
		return (*this = (*this * autre));
	}

	template<class U>
	erreur_b<T>& operator*=(U const& scalaire) {
		valeur *= scalaire;
		precision *= carre((float)scalaire);
		return *this;
	};


	explicit inline operator bool() const {
		if (carre((float)valeur) <= (20. * precision))
			return false;
		else
			return true;
	};



	friend erreur_b<T> operator*(erreur_b<T> const& temp1, erreur_b<T> const& temp2) {
		erreur_b<T> result(temp1.valeur * temp2.valeur, 0.);
		result.precision = carre((float)temp1.valeur) * temp2.precision  +  carre((float)temp2.valeur) * temp1.precision; // +(temp1.precision * temp2.precision);
		return result;
	};

	template<class U> 
	friend erreur_b<T> operator*(U const& scalaire, erreur_b<T> const& temp) {
		erreur_b<T> result(scalaire * temp.valeur, (carre((float)scalaire)) * temp.precision); //((T)scalaire)
		return result;
	};

	friend erreur_b<T> operator/(erreur_b<T> const& temp1, erreur_b<T> const& temp2) {
		/*
		erreur_b<T> result(temp1.valeur / temp2.valeur);
		result.precision = temp1.precision/abs(temp2.valeur) + ( temp2.precision *abs(temp1.valeur/(temp2.valeur*temp2.valeur)) ) + ( abs(result.valeur)* precision_relative(result.valeur));
		*/
		erreur_b<T> result(temp1.valeur / temp2.valeur, 0.);
		result.precision = temp1.precision / (carre((float)temp2.valeur))  +(temp2.precision * carre((float)temp1.valeur / carre((float)temp2.valeur)));
	//		result.maj();
		return result;
	};


	//A MODIFIER
	friend erreur_b<T> operator+(erreur_b<T> const& temp1, erreur_b<T> const& temp2) {
		erreur_b<T> result(temp1.valeur + temp2.valeur, 0.);
		if (carre((float)temp2.valeur) < carre((float)temp1.valeur) * precision_relative(temp1.valeur)) {
			result.precision = temp1.precision + carre((float)temp2.valeur);
			return result;
		}
		if (carre((float)temp1.valeur) < carre((float)temp2.valeur) * precision_relative(temp2.valeur)) {
			result.precision = temp2.precision + carre((float)temp1.valeur);
			return result;
		}

		float erreur_b1 = carre((float)temp1.valeur) * precision_relative(temp1.valeur); //erreur_b sur le calcul de x+y liée à x. registre
		float erreur_b2 = carre((float)temp2.valeur) * precision_relative(temp2.valeur); //erreur_b sur le calcul de x+y liée à y. registre
		if (erreur_b1 > erreur_b2) {
			if ((float)temp2.precision < erreur_b1) {
				result.precision = temp1.precision + erreur_b1;
				return result;
			}
			result.precision = temp1.precision + temp2.precision;
			return result;
		}
		if ((float)temp1.precision < erreur_b2) {
			result.precision = temp2.precision + erreur_b2;
			return result;
		}
		result.precision = temp2.precision + temp1.precision;
		return result;
	};

	//A MODIFIER
	friend erreur_b<T> operator-(erreur_b<T> const& temp1, erreur_b<T> const& temp2) {
		erreur_b<T> result(temp1.valeur - temp2.valeur,0.);
		//copié de operator+
		if (carre((float)temp2.valeur) < carre((float)temp1.valeur) * precision_relative(temp1.valeur)) {
			result.precision = temp1.precision + carre((float)temp2.valeur);
			return result;
		}
		if (carre((float)temp1.valeur) < carre((float)temp2.valeur) * precision_relative(temp2.valeur)) {
			result.precision = temp2.precision + carre((float)temp1.valeur);
			return result;
		}

		float erreur_b1 = carre((float)temp1.valeur) * precision_relative(temp1.valeur); //erreur_b sur le calcul de x+y liée à x
		float erreur_b2 = carre((float)temp2.valeur) * precision_relative(temp2.valeur); //erreur_b sur le calcul de x+y liée à y
		if (erreur_b1 > erreur_b2) {
			//		float x = carre((float)temp1.valeur) * precision_relative(temp1.valeur);
			if ((float)temp2.precision < erreur_b1) {
				result.precision = temp1.precision + erreur_b1;
				return result;
			}
			result.precision = temp1.precision + temp2.precision;
			return result;
		}
		//	float x = carre((float)temp2.valeur) * precision_relative(temp2.valeur);
		if ((float)temp1.precision < erreur_b2) {
			result.precision = temp2.precision + erreur_b2;
			return result;
		}
		result.precision = temp2.precision + temp1.precision;
		return result;
	};

	friend erreur_b<T> operator-(erreur_b<T> const& temp) {
		erreur_b<T> result(-temp.valeur,temp.precision);
		return result;
	};


	friend bool operator==(erreur_b<T> temp1, erreur_b<T> temp2) {
		return !((bool)(temp1 - temp2));
	};

	friend bool operator< (erreur_b<T> temp1, erreur_b<T> temp2) {
		return (temp1.valeur < temp2.valeur);
	};

	friend bool operator> (erreur_b<T> temp1, erreur_b<T> temp2) {
		return (temp1.valeur > temp2.valeur);
	};

	friend std::ostream& operator<<(std::ostream& os, const erreur_b<T>& element) {
		if ((bool) element)
			os << (double)element.valeur;
		else
			os << "0";
		return os;
	};

	operator T() {
		return valeur;
	};

	friend void swap(erreur_b<T>& gauche, erreur_b<T>& droit) {
		swap_F(gauche.valeur, droit.valeur);
		std::swap(gauche.precision, droit.precision);
		return;
	};
};


//LA MEME EN LINEAIRE

inline float precision_relative_l(double x) {
	return 1.0e-15;
};

inline float precision_relative_l(float x) {
	return 1.0e-6;
};


//test de précision ... Essayer avec 1. et 10 000. (par exemple)
inline float precision_relative_l(float_precision const& x) {
	return (float)x.epsilon();
};


template<class T> class erreur_l {
public:
	static_assert((type_algebre<T>::type == 0) && (type_algebre<T>::approx == 1), "erreur_b<T> : T doit être approché, avec division exacte");

	T valeur;
	float precision;

	erreur_l() :precision(0.), valeur(0.) {	};


	explicit erreur_l(T const& valeur_) {
		valeur = valeur_;
		precision = ((float)valeur) * precision_relative_l(valeur);
	};


	erreur_l(T const& valeur_, float const& precision_) { //faire attention au carré
		valeur = valeur_;
		precision = precision_;
		//		maj();
	};

	erreur_l(erreur_l<T> const& temp) {
		valeur = temp.valeur;
		precision = temp.precision;
	};

	erreur_l(erreur_l<T>&& temp) {
		swap(*this, temp);
		return;
	}

	explicit inline operator bool() const {
		if (abs((float) valeur) <= (2. * precision))
			return false;
		else
			return true;
	};

	erreur_l<T>& operator=(erreur_l<T> const& temp) {
		if (this == &temp)
			return *this;
		valeur = temp.valeur;
		precision = temp.precision;
		return *this;
	};

	erreur_l<T>& operator=(erreur_l<T>&& temp) {
		if (this == &temp)
			return *this;
		swap(*this, temp);
		return *this;
	}

	erreur_l<T>& operator+=(erreur_l<T> const& autre) {
		return (*this = (*this + autre));
	};

	erreur_l<T>& operator*=(erreur_l<T> const& autre) {
		return (*this = (*this * autre));
	};

	template<class U>
	erreur_l<T>& operator*=(U const& scalaire) {
		valeur *= scalaire;
		precision *= abs((float)scalaire);
		return *this;
	}

	friend erreur_l<T> operator*(erreur_l<T> const& temp1, erreur_l<T> const& temp2) {
		erreur_l<T> result(temp1.valeur * temp2.valeur, 0.);
		result.precision = abs((float)temp1.valeur) * temp2.precision + abs((float)temp2.valeur) * temp1.precision; // +(temp1.precision * temp2.precision);
		return result;
	};

	template<class U> 
	friend erreur_l<T> operator*(U const& scalaire, erreur_l<T> const& temp) {
		erreur_l<T> result(scalaire * temp.valeur, (abs((float)scalaire)) * temp.precision); //((T)scalaire)
		return result;
	};

	friend erreur_l<T> operator/(erreur_l<T> const& temp1, erreur_l<T> const& temp2) {
		/*
		erreur_b<T> result(temp1.valeur / temp2.valeur);
		result.precision = temp1.precision/abs(temp2.valeur) + ( temp2.precision *abs(temp1.valeur/(temp2.valeur*temp2.valeur)) ) + ( abs(result.valeur)* precision_relative(result.valeur));
		*/
		erreur_l<T> result(temp1.valeur / temp2.valeur, 0.);
		result.precision = temp1.precision / (abs((float)temp2.valeur)) + (temp2.precision * abs((float)temp1.valeur / carre((float)temp2.valeur)));
		//		result.maj();
		return result;
	};


	//A MODIFIER
	friend erreur_l<T> operator+(erreur_l<T> const& temp1, erreur_l<T> const& temp2) {
		erreur_l<T> result(temp1.valeur + temp2.valeur, 0.);
		if (abs((float)temp2.valeur) < abs((float)temp1.valeur) * precision_relative_l(temp1.valeur)) {
			result.precision = temp1.precision + abs((float)temp2.valeur);
			return result;
		}
		if (carre((float)temp1.valeur) < abs((float)temp2.valeur) * precision_relative_l(temp2.valeur)) {
			result.precision = temp2.precision + abs((float)temp1.valeur);
			return result;
		}

		float erreur_b1 = abs((float)temp1.valeur) * precision_relative_l(temp1.valeur); //erreur_b sur le calcul de x+y liée à x. registre
		float erreur_b2 = abs((float)temp2.valeur) * precision_relative_l(temp2.valeur); //erreur_b sur le calcul de x+y liée à y. registre
		if (erreur_b1 > erreur_b2) {
			if ((float)temp2.precision < erreur_b1) {
				result.precision = temp1.precision + erreur_b1;
				return result;
			}
			result.precision = temp1.precision + temp2.precision;
			return result;
		}
		if ((float)temp1.precision < erreur_b2) {
			result.precision = temp2.precision + erreur_b2;
			return result;
		}
		result.precision = temp2.precision + temp1.precision;
		return result;
	};

	//A MODIFIER
	friend erreur_l<T> operator-(erreur_l<T> const& temp1, erreur_l<T> const& temp2) {
		erreur_l<T> result(temp1.valeur - temp2.valeur, 0.);
		//copié de operator+
		if (abs((float)temp2.valeur) < abs((float)temp1.valeur) * precision_relative_l(temp1.valeur)) {
			result.precision = temp1.precision + abs((float)temp2.valeur);
			return result;
		}
		if (abs((float)temp1.valeur) < abs((float)temp2.valeur) * precision_relative_l(temp2.valeur)) {
			result.precision = temp2.precision + abs((float)temp1.valeur);
			return result;
		}

		float erreur_b1 = abs((float)temp1.valeur) * precision_relative_l(temp1.valeur); //erreur_b sur le calcul de x+y liée à x
		float erreur_b2 = abs((float)temp2.valeur) * precision_relative_l(temp2.valeur); //erreur_b sur le calcul de x+y liée à y
		if (erreur_b1 > erreur_b2) {
			//		float x = carre((float)temp1.valeur) * precision_relative(temp1.valeur);
			if ((float)temp2.precision < erreur_b1) {
				result.precision = temp1.precision + erreur_b1;
				return result;
			}
			result.precision = temp1.precision + temp2.precision;
			return result;
		}
		//	float x = carre((float)temp2.valeur) * precision_relative(temp2.valeur);
		if ((float)temp1.precision < erreur_b2) {
			result.precision = temp2.precision + erreur_b2;
			return result;
		}
		result.precision = temp2.precision + temp1.precision;
		return result;
	};

	friend erreur_l<T> operator-(erreur_l<T> const& temp) {
		erreur_l<T> result(-temp.valeur, temp.precision);
		return result;
	};

	friend bool operator==(erreur_l<T> temp1, erreur_l<T> temp2) {
		return !((bool)(temp1 - temp2));
	};

	friend bool operator< (erreur_l<T> temp1, erreur_l<T> temp2) {
		return (temp1.valeur < temp2.valeur);
	};

	friend bool operator> (erreur_l<T> temp1, erreur_l<T> temp2) {
		return (temp1.valeur > temp2.valeur);
	};

	friend std::ostream& operator<<(std::ostream& os, const erreur_l<T>& element) {
		if ((bool)element)
			os << (double)element.valeur;
		else
			os << "0";
		return os;
	};

	operator T() {
		return valeur;
	};

	friend void swap(erreur_l<T>& gauche, erreur_l<T>& droit) {
		swap_F(gauche.valeur, droit.valeur);
		std::swap(gauche.precision, droit.precision);
		return;
	};

};

