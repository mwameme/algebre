#pragma once
#include "types.hpp"
#include "fonctions template.hpp"
#include <exception>
#include <vector>

#include <complex>


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

template<class T> class type_algebre;

template<class T> T PGCD(T x, T y);

template<class T> T PPCM(T const& a, T const& b);

template<class T> class calcul_objet;

template<class U> class calcul_objet<std::vector<U>> {
public:

	static decltype(calcul_objet<U>::PPCM(U())) PPCM(std::vector<U> const& vec) {
		//		if (vec.size() == 0)
		//			return decltype(calcul_objet<U>::ppcm(U()))();
		if (vec.size() == 0)
			throw std::domain_error("PPCM_objet : liste vide");

		decltype(calcul_objet<U>::PPCM(U())) temp = calcul_objet<U>::PPCM(vec[0]);
		for (int i(1); i < vec.size(); ++i)
			if ((bool) vec[i])
				temp = PPCM(temp, calcul_objet<U>::PPCM(vec[i]));
		return temp;
	};

	static decltype(calcul_objet<U>::PGCD(U())) PGCD(std::vector<U> const& vec) {
		//		if (vec.size() == 0)
		//			return decltype(calcul_objet<U>::ppcm(U()))();
		if (vec.size() == 0)
			throw std::domain_error("PPCM_objet : liste vide");

		decltype(calcul_objet<U>::PGCD(U())) temp = calcul_objet<U>::PGCD(vec[0]);
		for (int i(1); i < vec.size(); ++i)
			if((bool) vec[i])
				temp = PGCD(temp, calcul_objet<U>::PPCM(vec[i]));
		return temp;
	};

};

template<class T> class calcul_objet<rationnel<T>> {
public:
	static T PPCM(rationnel<T> const& ratio) {
		return ratio.denominateur;
	};

	static T PGCD(rationnel<T> const& ratio) {
		return ratio.numerateur;
	};
};

template<class T> class calcul_objet<complexe<T>> {
public:
	static decltype(calcul_objet<T>::PPCM(T())) PPCM(complexe<T> const& z) {
		return PPCM(calcul_objet<T>::PPCM(z.x), calcul_objet<T>::PPCM(z.y));
	};

	static decltype(calcul_objet<T>::PGCD(T())) PGCD(complexe<T> const& z) {
		return PGCD(calcul_objet<T>::PGCD(z.x), calcul_objet<T>::PGCD(z.y));
	};

};

template<class T> class calcul_objet<std::complex<T>> {
public:
	static decltype(calcul_objet<T>::PPCM(T())) PPCM(std::complex<T> const& z) {
		return PPCM(calcul_objet<T>::PPCM(z.real), calcul_objet<T>::PPCM(z.imag));
	};

	static decltype(calcul_objet<T>::PGCD(T())) PGCD(std::complex<T> const& z) {
		return PGCD(calcul_objet<T>::PGCD(z.real), calcul_objet<T>::PGCD(z.imag));
	};
};

template<class T> class calcul_objet<polynome_n_rec<T>> {
public:
	static decltype(calcul_objet<T>::PPCM(T())) PPCM(polynome_n_rec<T> const& poly) {
		if (poly.nul == false) {
			T temp = unite(poly.element, false);
			return calcul_objet<T>::PPCM(temp);
		}
		if (poly.n_var == 0) {
			return calcul_objet<T>::PPCM(poly.element);
		}
		return calcul_objet<std::vector<polynome_n_rec<T>>>::PPCM(poly.poly.coeffs);
	};

	static decltype(calcul_objet<T>::PGCD(T())) PGCD(polynome_n_rec<T> const& poly) {
		if (poly.nul == false) {
			T temp = unite(poly.element, false);
			return calcul_objet<T>::PGCD(temp);
		}
		if (poly.n_var == 0) {
			return calcul_objet<T>::PGCD(poly.element);
		}
		return calcul_objet<std::vector<polynome_n_rec<T>>>::PGCD(poly.poly.coeffs);
	};
};

template<class T, int n> class calcul_objet<polynome_n_fixe<T, n>> {
public:
	static decltype(calcul_objet<T>::PPCM(T())) PPCM(polynome_n_fixe<T, n> const& poly) {
		if (poly.nul == false) {
			T temp = unite(poly.get_T(), false);
			return calcul_objet<T>::PPCM(temp);
		}
		return calcul_objet<std::vector<polynome_n_fixe<T, n - 1>>>::PPCM(poly.poly.coeffs);
	};

	static decltype(calcul_objet<T>::PGCD(T())) PGCD(polynome_n_fixe<T, n> const& poly) {
		if (poly.nul == false) {
			T temp = unite(poly.get_T(), false);
			return calcul_objet<T>::PGCD(temp);
		}
		return calcul_objet<std::vector<polynome_n_fixe<T, n - 1>>>::PGCD(poly.poly.coeffs);
	};
};

template<class T> class calcul_objet<polynome_n_fixe<T, 0>> {
public:
	static decltype(calcul_objet<T>::PPCM(T())) PPCM(polynome_n_fixe<T, 0> const& poly) {
		return calcul_objet<T>::PPCM(poly.element);
	};

	static decltype(calcul_objet<T>::PGCD(T())) PGCD(polynome_n_fixe<T, 0> const& poly) {
		return calcul_objet<T>::PGCD(poly.element);
	};
};

template<class T> class calcul_objet<polynome_n_iter<T>> {
public:
	static decltype(calcul_objet<T>::PPCM(T())) PPCM(polynome_n_iter<T> const& poly) {
		return calcul_objet<std::vector<T>>::PPCM(poly.coeffs.data);
	};

	static decltype(calcul_objet<T>::PGCD(T())) PGCD(polynome_n_iter<T> const& poly) {
		return calcul_objet<std::vector<T>>::PGCD(poly.coeffs.data);
	};
};

template<class T> class calcul_objet<polynome<T>> {
public:
	static decltype(calcul_objet<T>::PPCM(T())) PPCM(polynome<T> const& poly) {
		return calcul_objet<std::vector<T>>::PPCM(poly.coeffs);
	};

	static decltype(calcul_objet<T>::PGCD(T())) PGCD(polynome<T> const& poly) {
		return calcul_objet<std::vector<T>>::PGCD(poly.coeffs);
	};
};

template<class T> class descendre_simplification;

template<class T> class descendre_simplification<rationnel<T>> {
public:
	static void simplifier(rationnel<T>& ratio) {
		if constexpr (type_algebre<T>::is_objet && type_algebre<type_algebre<T>::corps>::is_corps_ratio) {
			if constexpr (type_algebre<T>::is_ratio) {
				descendre_simplification<T>::simplifier(ratio.numerateur);
				descendre_simplification<T>::simplifer(ratio.denominateur);
				return;
			}
			else {
				auto ppcm = PPCM(calcul_objet<T>::PPCM(ratio.numerateur), calcul_objet<T>::PPCM(ratio.denominateur));
				auto pgcd = PGCD(calcul_objet<T>::PGCD(ratio.numerateur), calcul_objet<T>::PGCD(ratio.denominateur));
				auto frac = rationnel<decltype(pgcd)>(ppcm, pgcd);
				ratio.numerateur = frac * ratio.numerateur;
				ratio.denominateur = frac * ratio.denominateur;
				return;
			}
		}
		else
			return;
	};
};

template<class T> class descendre_simplification<polynome<T>> {
public:
	static void simplifier(polynome<T>& poly) {
		if constexpr (type_algebre<T>::is_ratio && type_algebre<type_algebre<T>::corps>::is_corps_ratio) {
			for (int i(0); i < poly.coeffs.size();++i)
				if ((bool)poly.coeffs[i])
					descendre_simplification<T>::simplifier(poly.coeffs[i]);
		}
		else
			return;
	};
};

template<class T> class descendre_simplification<polynome_n_iter<T>> {
public:
	static void simplifier(polynome_n_iter<T>& poly) {
		if constexpr (type_algebre<T>::is_ratio && type_algebre<type_algebre<T>::corps>::is_corps_ratio) {
			for (int i(0); i < poly.coeffs.data.size();++i)
				if ((bool)poly.coeffs.data[i])
					descendre_simplification<T>::simplifier(poly.coeffs.data[i]);
		}
		else
			return;
	};
};

template<class T> class descendre_simplification<polynome_n_rec<T>> {
public:
	static void simplifier(polynome_n_rec<T>& poly) {
		if constexpr (type_algebre<T>::is_ratio && type_algebre<type_algebre<T>::corps>::is_corps_ratio) {
			if (poly.n_var > 0) {
				for (int i(0); i < poly.poly.coeffs.size(); ++i)
					if ((bool)poly.poly.coeffs[i])
						descendre_simplification<polynome_n_rec<T>>::simplifier(poly.poly.coeffs[i]);
			}
			else
				descendre_simplification<T>::simplifier(poly.element);
		}
		else
			return;
	};
};

template<class T, int n> class descendre_simplification<polynome_n_fixe<T, n>> {
public:
	static void simplifier(polynome_n_fixe<T, n>& poly) {
		if constexpr (type_algebre<T>::is_ratio && type_algebre<type_algebre<T>::corps>::is_corps_ratio) {
			for (int i(0); i < poly.poly.coeffs.size();++i)
				if ((bool)poly.poly.coeffs[i])
					descendre_simplification<polynome_n_fixe<T, n - 1>>::simplifier(poly.poly.coeffs[i]);
		}
		else
			return;
	};
};

template<class T> class descendre_simplification<polynome_n_fixe<T, 0>> {
public:
	static void simplifier(polynome_n_fixe<T, 0>& poly) {
		if constexpr (type_algebre<T>::is_ratio && type_algebre<type_algebre<T>::corps>::is_corps_ratio) {
			descendre_simplification<T>::simplifier(poly.element);
		}
		else
			return;
	};
};

template<class T> class descendre_simplification<matrice<T>> {
public:
	static void simplifier(matrice<T>& mat) {
		if constexpr (type_algebre<T>::is_ratio && type_algebre<type_algebre<T>::corps>::is_corps_ratio) {
			for (int i(0); i < mat.taille_l; ++i)
				for (int j(0); j < mat.taille_c; ++j)
					if ((bool) mat.coeffs[i][j])
						descendre_simplification<T>::simplifier(mat.coeffs[i][j]);
		}
		else
			return;
	};
};

//pour les noyaux ... multiplication par un scalaire !
template<class T> class descendre_simplification<std::vector<T>> {
public:
	static void simplifier(std::vector<T>& vec) {
		if constexpr (type_algebre<T>::is_ratio && type_algebre<type_algebre<T>::corps>::is_corps_ratio) {
			for (int i(0); i < vec.size(); ++i)
				if ((bool) vec[i])
					descendre_simplification<T>::simplifier(vec[i]);
		}
		else if constexpr (type_algebre<T>::is_objet && type_algebre<type_algebre<T>::corps>::is_corps_ratio) {
			auto ppcm = calcul_objet<std::vector<T>>::PPCM(vec);
			auto pgcd = calcul_objet<std::vector<T>>::PGCD(vec);
			auto frac = rationnel<decltype(pgcd)>(ppcm, pgcd);
			for (int i(0); i < vec.size(); ++i)
				vec[i] = frac * vec[i];
			return;
		}
	};
};

template<class T> inline void simplifier_frac(T& element) {
	descendre_simplification<T>::simplifier(element);
};