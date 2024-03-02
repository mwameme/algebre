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

	static decltype(calcul_objet<U>::PPCM_objet(U())) PPCM_objet(std::vector<U> const& vec) {
		//		if (vec.size() == 0)
		//			return decltype(calcul_objet<U>::PPCM_objet(U()))();
		if (vec.size() == 0)
			throw std::domain_error("PPCM_objet : liste vide");

		decltype(calcul_objet<U>::PPCM_objet(U())) temp = calcul_objet<U>::PPCM_objet(vec[0]);
		for (int i(1); i < vec.size(); ++i)
			if ((bool)vec[i])
				temp = PPCM(temp, calcul_objet<U>::PPCM_objet(vec[i]));
		return temp;
	};

	static decltype(calcul_objet<U>::PGCD_objet(U())) PGCD_objet(std::vector<U> const& vec) {
		//		if (vec.size() == 0)
		//			return decltype(calcul_objet<U>::PPCM_objet(U()))();
		if (vec.size() == 0)
			throw std::domain_error("PPCM_objet : liste vide");

		decltype(calcul_objet<U>::PGCD_objet(U())) temp = calcul_objet<U>::PGCD_objet(vec[0]);
		for (int i(1); i < vec.size(); ++i)
			if ((bool)vec[i])
				temp = PGCD(temp, calcul_objet<U>::PPCM_objet(vec[i]));
		return temp;
	};

};

template<class T> class calcul_objet<rationnel<T>> {
public:
	static T PPCM_objet(rationnel<T> const& ratio) {
		return ratio.denominateur;
	};

	static T PGCD_objet(rationnel<T> const& ratio) {
		return ratio.numerateur;
	};
};

template<class T> class calcul_objet<complexe<T>> {
public:
	static decltype(calcul_objet<T>::PPCM_objet(T())) PPCM_objet(complexe<T> const& z) {
		return PPCM(calcul_objet<T>::PPCM_objet(z.x), calcul_objet<T>::PPCM_objet(z.y));
	};

	static decltype(calcul_objet<T>::PGCD_objet(T())) PGCD_objet(complexe<T> const& z) {
		return PGCD(calcul_objet<T>::PGCD_objet(z.x), calcul_objet<T>::PGCD_objet(z.y));
	};

};

template<class T> class calcul_objet<std::complex<T>> {
public:
	static decltype(calcul_objet<T>::PPCM_objet(T())) PPCM_objet(std::complex<T> const& z) {
		return PPCM(calcul_objet<T>::PPCM_objet(z.real), calcul_objet<T>::PPCM_objet(z.imag));
	};

	static decltype(calcul_objet<T>::PGCD_objet(T())) PGCD_objet(std::complex<T> const& z) {
		return PGCD(calcul_objet<T>::PGCD_objet(z.real), calcul_objet<T>::PGCD_objet(z.imag));
	};
};

template<class T> class calcul_objet<polynome_n_rec<T>> {
public:
	static decltype(calcul_objet<T>::PPCM_objet(T())) PPCM_objet(polynome_n_rec<T> const& poly) {
		auto it = poly.cbegin();
		auto result = calcul_objet<T>::PPCM_objet(*it);
		++it;
		for (; (bool)it; ++it)
			result = PPCM(result, calcul_objet<T>::PPCM_objet(*it));
		return result;
	};

	static decltype(calcul_objet<T>::PGCD_objet(T())) PGCD_objet(polynome_n_rec<T> const& poly) {
		auto it = poly.cbegin();
		auto result = calcul_objet<T>::PGCD_objet(*it);
		++it;
		for (; (bool)it; ++it)
			result = PGCD(result, calcul_objet<T>::PGCD_objet(*it));
		return result;
	};
};

template<class T> class calcul_objet<polynome_n_sparse<T>> {
public:
	static decltype(calcul_objet<T>::PPCM_objet(T())) PPCM_objet(polynome_n_sparse<T> const& poly) {
		auto it = poly.cbegin();
		auto result = calcul_objet<T>::PPCM_objet(*it.element);
		++it;
		for (; (bool)it; ++it)
			result = PPCM(result, calcul_objet<T>::PPCM_objet(*it.element));
		return result;
	};

	static decltype(calcul_objet<T>::PGCD_objet(T())) PGCD_objet(polynome_n_sparse<T> const& poly) {
		auto it = poly.cbegin();
		auto result = calcul_objet<T>::PGCD_objet(*it.element);
		++it;
		for (; (bool)it; ++it)
			result = PGCD(result, calcul_objet<T>::PGCD_objet(*it.element));
		return result;
	};
};


template<class T, int n> class calcul_objet<polynome_n_fixe<T, n>> {
public:
	static decltype(calcul_objet<T>::PPCM_objet(T())) PPCM_objet(polynome_n_fixe<T, n> const& poly) {
		auto it = poly.cbegin();
		auto result = calcul_objet<T>::PPCM_objet(*it);
		++it;
		for (; (bool)it; ++it)
			result = PPCM(result, calcul_objet<T>::PPCM_objet(*it));
		return result;
	};

	static decltype(calcul_objet<T>::PGCD_objet(T())) PGCD_objet(polynome_n_fixe<T, n> const& poly) {
		auto it = poly.cbegin();
		auto result = calcul_objet<T>::PGCD_objet(*it);
		++it;
		for (; (bool)it; ++it)
			result = PGCD(result, calcul_objet<T>::PGCD_objet(*it));
		return result;
	};
};

template<class T> class calcul_objet<polynome_n_fixe<T, 0>> {
public:
	static decltype(calcul_objet<T>::PPCM_objet(T())) PPCM_objet(polynome_n_fixe<T, 0> const& poly) {
		return calcul_objet<T>::PPCM_objet(poly.element);
	};

	static decltype(calcul_objet<T>::PGCD_objet(T())) PGCD_objet(polynome_n_fixe<T, 0> const& poly) {
		return calcul_objet<T>::PGCD_objet(poly.element);
	};
};

template<class T> class calcul_objet<polynome_n_iter<T>> {
public:
	static decltype(calcul_objet<T>::PPCM_objet(T())) PPCM_objet(polynome_n_iter<T> const& poly) {
		return calcul_objet<std::vector<T>>::PPCM_objet(poly.coeffs.data);
	};

	static decltype(calcul_objet<T>::PGCD_objet(T())) PGCD_objet(polynome_n_iter<T> const& poly) {
		return calcul_objet<std::vector<T>>::PGCD_objet(poly.coeffs.data);
	};
};

template<class T> class calcul_objet<polynome<T>> {
public:
	static decltype(calcul_objet<T>::PPCM_objet(T())) PPCM_objet(polynome<T> const& poly) {
		return calcul_objet<std::vector<T>>::PPCM_objet(poly.coeffs);
	};

	static decltype(calcul_objet<T>::PGCD_objet(T())) PGCD_objet(polynome<T> const& poly) {
		return calcul_objet<std::vector<T>>::PGCD_objet(poly.coeffs);
	};
};

template<class T> class descendre_simplification;

template<class T> class descendre_simplification<rationnel<T>> {
public:
	static void simplifier(rationnel<T>& ratio) {
		if constexpr (type_algebre<T>::is_objet && type_algebre<T>::is_corps_ratio) {
			if constexpr (type_algebre<T>::is_ratio) {
				descendre_simplification<T>::simplifier(ratio.numerateur);
				descendre_simplification<T>::simplifer(ratio.denominateur);
				return;
			}
			else {
				auto ppcm = PPCM(calcul_objet<T>::PPCM_objet(ratio.numerateur), calcul_objet<T>::PPCM_objet(ratio.denominateur));
				auto pgcd = PGCD(calcul_objet<T>::PGCD_objet(ratio.numerateur), calcul_objet<T>::PGCD_objet(ratio.denominateur));
				if (!(bool)pgcd)
					return;
				auto frac = rationnel<decltype(pgcd)>(ppcm, pgcd);
				ratio.numerateur *= frac;
				ratio.denominateur *= frac;
				frac = rationnel<decltype(pgcd)>(unite(pgcd, true), PGCD(ratio.numerateur) * PGCD(ratio.denominateur));
				if (!(bool)frac.denominateur)
					return;
				ratio.numerateur *= frac;
				ratio.denominateur *= frac;
#ifdef _DEBUG
				if (!(bool)denominateur)
					throw std::domain_error("un 0 au denominateur d'une fraction");
#endif
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
		if constexpr (type_algebre<T>::is_ratio && type_algebre<T>::is_corps_ratio) {
			for (int i(0); i < poly.coeffs.size(); ++i)
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
		if constexpr (type_algebre<T>::is_ratio && type_algebre<T>::is_corps_ratio) {
			for (int i(0); i < poly.coeffs.data.size(); ++i)
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
		if constexpr (type_algebre<T>::is_ratio && type_algebre<T>::is_corps_ratio) {
			for(auto it = poly.begin();(bool) it;++it)
				descendre_simplification<T>::simplifier(*it);
		}
		else
			return;
	};
};

template<class T, int n> class descendre_simplification<polynome_n_fixe<T, n>> {
public:
	static void simplifier(polynome_n_fixe<T, n>& poly) {
		if constexpr (type_algebre<T>::is_ratio && type_algebre<T>::is_corps_ratio) {
			for (auto it = poly.begin(); (bool)it; ++it)
				descendre_simplification<T>::simplifier(*it);

		}
		else
			return;
	};
};

template<class T> class descendre_simplification<polynome_n_fixe<T, 0>> {
public:
	static void simplifier(polynome_n_fixe<T, 0>& poly) {
		if constexpr (type_algebre<T>::is_ratio && type_algebre<T>::is_corps_ratio) {
			descendre_simplification<T>::simplifier(poly.element);
		}
		else
			return;
	};
};

template<class T> class descendre_simplification<matrice<T>> {
public:
	static void simplifier(matrice<T>& mat) {
		if constexpr (type_algebre<T>::is_ratio && type_algebre<T>::is_corps_ratio) {
			for (int i(0); i < mat.taille_l; ++i)
				for (int j(0); j < mat.taille_c; ++j)
					if ((bool)mat.coeffs[i][j])
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
		if constexpr (type_algebre<T>::is_ratio && type_algebre<T>::is_corps_ratio) {
			for (int i(0); i < vec.size(); ++i)
				if ((bool)vec[i])
					descendre_simplification<T>::simplifier(vec[i]);
		}
		else if constexpr (type_algebre<T>::is_objet && type_algebre<T>::is_corps_ratio) {
			auto ppcm = calcul_objet<std::vector<T>>::PPCM_objet(vec);
			auto pgcd = calcul_objet<std::vector<T>>::PGCD_objet(vec);
			auto frac = rationnel<decltype(pgcd)>(ppcm, pgcd);
#ifdef _DEBUG
			if (!(bool)frac.denominateur)
				throw std::domain_error("un 0 au denominateur d'une fraction");
#endif
			for (int i(0); i < vec.size(); ++i)
				vec[i] = frac * vec[i];
			return;
		}
	};
};

template<class T> inline void simplifier_frac(T& element) {
	descendre_simplification<T>::simplifier(element);
};