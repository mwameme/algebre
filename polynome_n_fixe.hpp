#pragma once
#include <vector>
#include <string>
#include <exception>

#include "polynome.hpp"
#include "polynome_n_rec.hpp"
#include "entete objets.hpp"

#include "swap_T.hpp"
#include "polynome_n_sparse.hpp"

template<class T, int n> class polynome_n_fixe;


template<class T> class polynome_n_rec;

template<class T, int n> class pointeur_n;

/*
template <auto Start, auto End, auto Inc, class F>
constexpr void constexpr_for(F&& f)
{
	if constexpr (Inc > 0) {
		if constexpr (Start <= End)
		{
			f(std::integral_constant<decltype(Start), Start>());
			constexpr_for<Start + Inc, End, Inc>(f);
		}
	}
	else 
		if constexpr (Start >= End)
		{
			f(std::integral_constant<decltype(Start), Start>());
			constexpr_for<Start + Inc, End, Inc>(f);
		}
}
*/

template <int start, int end, class F>
constexpr inline void constexpr_for(F&& f) //utilisé pour g, croissant
{
	if constexpr (start <= end)
	{
		(f(std::integral_constant<int, start>()));
		constexpr_for<start + 1, end>(f);
	}
	else return;

};

template <int start, int end, int end2, class F, class G>
constexpr inline void constexpr_for_2(F&& f, G&& g)
{
	if constexpr (start >= end)
	{
		if (f(std::integral_constant<int, start>()))
			constexpr_for_2<start - 1, end, end2>(f, g);
		else
			constexpr_for<start + 1, end2>(g);
	}
	else return;
};



template<class T, int n> class polynome_n_fixe {
public:

	polynome<polynome_n_fixe<T, n - 1>> poly;
	bool nul;

	polynome_n_fixe() {};
	polynome_n_fixe(T element) : poly(polynome_n_fixe<T, n - 1>(element)) ,nul((bool) element) {};
	polynome_n_fixe(monome<T> monome_) : polynome_n_fixe(monome_.degres, monome_.element) {};
	polynome_n_fixe(std::vector<int> vec, T element); //monome. vec : degrés
	polynome_n_fixe(int* vec, T element, T faux);  //monome interne
	polynome_n_fixe(const polynome_n_fixe<T, n>& copie) : poly(copie.poly), nul(copie.nul) {	};

	polynome_n_fixe(polynome_n_fixe<T, n>&& temp) {
		nul = temp.nul;
		swap(poly, temp.poly);
		return;
	};

	polynome_n_fixe<T, n>& operator=(const polynome_n_fixe<T, n>& copie) {
		if (this == &copie)
			return *this;
		nul = copie.nul;
		poly = copie.poly;
		return *this;
	};

	polynome_n_fixe<T, n>& operator=(polynome_n_fixe<T, n>&& temp) {
		if (this == &temp)
			return *this;
		nul = temp.nul;
		swap(poly, temp.poly);
		return *this;
	};

	polynome_n_fixe<T, n>& operator*=(const polynome_n_fixe<T, n>& temp) {
		return (*this = (*this * temp));
	};

	template<class U>
	polynome_n_fixe<T, n>& operator*=(U const& scalaire) {
		poly *= scalaire;
		nul = (bool)poly;
		return *this;
	};

	polynome_n_fixe<T, n>& operator+=(polynome_n_fixe<T, n> const& autre) {
		poly += autre.poly;
		nul = (bool)poly;
		return *this;
	};

	friend polynome_n_fixe<T, n> operator*(const polynome_n_fixe<T, n>& gauche, const polynome_n_fixe<T, n>& droit) {
		if (!gauche.nul)
			return gauche;
		if (!droit.nul)
			return droit;

		polynome_n_fixe<T, n> m_poly;
		m_poly.poly = gauche.poly * droit.poly;
		m_poly.nul = (bool)m_poly.poly;

		return m_poly;
	};

	friend polynome_n_fixe<T, n> operator+(const polynome_n_fixe<T, n>& gauche, const polynome_n_fixe<T, n>& droit) {
		if (!gauche.nul)
			return droit;
		if (!droit.nul)
			return gauche;
		polynome_n_fixe<T, n> m_poly;
		m_poly.poly = gauche.poly + droit.poly;
		m_poly.nul = (bool)m_poly.poly;
		return m_poly;
	};

	friend polynome_n_fixe<T, n> operator-(const polynome_n_fixe<T, n>& gauche, const polynome_n_fixe<T, n>& droit) {
		polynome_n_fixe<T, n> m_poly;
		m_poly.poly = gauche.poly - droit.poly;
		m_poly.nul = (bool)m_poly.poly;
		return m_poly;
	};

	friend polynome_n_fixe<T, n> operator-(const polynome_n_fixe<T, n>& element) {
		polynome_n_fixe<T, n> m_poly(element);
		m_poly.poly = -m_poly.poly;
		m_poly.nul = (bool)m_poly.poly;
		return m_poly;
	};

	template<class U> friend polynome_n_fixe<T, n> operator*(const U& scalaire, const polynome_n_fixe<T, n>& polynome_) {
		polynome_n_fixe<T, n> m_poly(polynome_);
		m_poly.poly = scalaire * m_poly.poly;
		m_poly.nul = (bool)m_poly.poly;

		return m_poly;
	};


	friend bool operator==(polynome_n_fixe<T, n> const& gauche, polynome_n_fixe<T, n> const& droit) {
		return gauche.poly == droit.poly;
	};

	T get_T() const {
		return poly.coeffs[0].get_T();
	};

	explicit operator bool() const {
		return nul;
	};

	friend void swap(polynome_n_fixe<T, n>& gauche, polynome_n_fixe<T, n>& droit) {
		std::swap(gauche.nul, droit.nul);
		swap(gauche.poly, droit.poly);
		return;
	};

	int lenght() {
		int somme = 0;
		for (int i(0); i < poly.coeffs.size(); ++i)
			somme += poly.coeffs[i].length();
		return somme+1;
	};

	template<class I>
	class iterator_polynome_n_fixe;

	using iterator = iterator_polynome_n_fixe<T>;
	using const_iterator = iterator_polynome_n_fixe<const T>;


	iterator begin() {
		return iterator(*this);
	};

	iterator end() {
		iterator it(*this);
		it.termine = false;
		return it;
	};


	const_iterator cbegin() const {
		return const_iterator(*this);
	};

	const_iterator cend() const {
		const_iterator it(*this);
		it.termine = false;
		return it;
	};

	template<class U>
	operator polynome_n_fixe<U, n>() {
		polynome_n_fixe<U,n> result();
		result.poly = (polynome<polynome_n_fixe<U, n - 1>>) poly;
		result.nul = (bool) result.poly;

		return result;
	}

	operator polynome_n_rec<T>() {
		polynome_n_rec<T> result();
		result.n_var = n;
		result.element = unite(get_T(), false);
		result.nul = nul;
		result.poly = (polynome<polynome_n_rec<T>>) poly;

		return result;
	};

	int max_degre(int i) const {
		if (i == 0)
			return poly.degre;
		int result = -1;
		for (int j(0); j < poly.coeffs.size(); ++j)
			result = max(result, poly.coeffs[i].max_degre(i-1));
		return result;
	}

	std::vector<int> max_degre() const {
		std::vector<int> result(n, -1);
		max_degre(result.data());
		return result;
	};

	void max_degre(int* p) const{
		if (poly.degre > p[0])
			p[0] = poly.degre;
		for (int i(0); i < poly.coeffs.size(); ++i)
			if((bool) poly.coeffs[i])
				poly.coeffs[i].max_degre(p + 1);
		return;
	};
};

template<class T, int n> class pointeur_n;

template<class T> class pointeur_n<T, 1> {
public:
	polynome_n_fixe<T, 1>* pointeur;

	template<int m>
	polynome_n_fixe<T, 1>*& get() {
		return pointeur;
	};

	pointeur_n(polynome_n_fixe<T, 1>& poly) {
		pointeur = &poly;
	};

};

template<class T, int n> class pointeur_n {
public:
	polynome_n_fixe<T, n>* pointeur;
	pointeur_n<T, n - 1> pointeurs;

	template<int m>
	auto*& get() {
		if constexpr (m == 0)
			return pointeur;
		else
			return pointeurs.get<m - 1>();
	};

	pointeur_n(polynome_n_fixe<T, n>& poly) : pointeur(&poly), pointeurs(poly.poly.coeffs[0]) {	};
};


template<class T, int n> class pointeur_n_const;

template<class T> class pointeur_n_const<T, 1> {
public:
	const polynome_n_fixe<T, 1>* pointeur;

	template<int m>
	const polynome_n_fixe<T, 1>*& get() {
		//		static_assert(m == 0);
		return pointeur;
	};

	pointeur_n_const(const polynome_n_fixe<T, 1>& poly) {
		pointeur = &poly;
	};

};

template<class T, int n> class pointeur_n_const {
public:
	const polynome_n_fixe<T, n>* pointeur;
	pointeur_n_const<T, n - 1> pointeurs;

	template<int m>
	auto*& get() {
		if constexpr (m == 0)
			return pointeur;
		else
			return pointeurs.get<m - 1>();
	};

	pointeur_n_const(const polynome_n_fixe<T, n>& poly) : pointeur(&poly), pointeurs(poly.poly.coeffs[0]) {	};
};

template<class T, int n> template<class I>
class polynome_n_fixe<T, n>::iterator_polynome_n_fixe {
public:
	using value_type = std::remove_const_t<I>;
	using poly_type = std::conditional_t< std::is_const_v<I>, const polynome_n_fixe<value_type, n>, polynome_n_fixe<value_type, n>	>;

	std::conditional_t < std::is_const_v<I>, pointeur_n_const<T, n>, pointeur_n<T, n> > pointeurs;
	std::vector<int> positions;
	bool termine;

	iterator_polynome_n_fixe& operator++() {
		constexpr_for_2<n - 1, 0, n - 1>([&](auto i) { //premiere fonction : mettre a jour les positions.
			++positions[i];
			if (positions[i] >= pointeurs.get<i>()->poly.coeffs.size()) {
				positions[i] = 0;
				if ((int) i == 0) {
					termine = false;
					return false;
				}
				return true;
			}
			else
				return false; //lors de return false, position i est bonne. recommencer à i+1

			}, [&](auto j) { //deuxieme fonction : mettre a jour les pointeurs.
				pointeurs.get<(int)j>() = &pointeurs.get<(int)j - 1>()->poly.coeffs[positions[(int)j - 1]];
				return;
			});

		return *this;
	};

	void go_position(std::vector<int> new_positions) {
		positions = new_positions;

		constexpr_for<1, n - 1>([&](auto j) { //deuxieme fonction : mettre a jour les pointeurs.
#ifdef ALGEBRA_USE_EXCEPTION
			if (positions[j - 1] >= pointeurs.get<j - 1>()->poly.coeffs.size())
				throw std::domain_error("poolynome_n_fixe::iterator, position hors domaine");
#endif
			pointeurs.get<j>() = &pointeurs.get<j - 1>()->poly.coeffs[positions[j - 1]];
			return;
			});
		return;
	}

	iterator_polynome_n_fixe(poly_type& temp) : pointeurs(temp), positions(n, 0), termine(true) {	};

	operator bool() const {
		return termine;
	};

	I& operator*() {
		return pointeurs.get<n>()->poly.coeffs[positions[n - 1]].element;
	};

	friend bool operator==(iterator_polynome_n_fixe const& it1, iterator_polynome_n_fixe const& it2) {
		if ((it1.termine == it2.termine) && (it1.pointeurs.get<0>() == it2.pointeurs.get<0>()) && (it1.positions == it1.positions))
			return true;
		else
			return false;
	};

	friend bool operator!=(iterator_polynome_n_fixe const& it1, iterator_polynome_n_fixe const& it2) {
		if (it1.termine != it2.termine)
			return true;
		if (it1.positions != it2.termine)
			return true;
		if (it1.pointeurs.get<0>() != it2.pointeurs.get<0>())
			return true;
		return false;
	};

};


template<class T> class polynome_n_fixe<T, 1> {
public:
	polynome<polynome_n_fixe<T, 0>> poly;
	bool nul;

	polynome_n_fixe() {};
	polynome_n_fixe(T element) : poly(polynome_n_fixe<T,0>(element)), nul((bool)element) { };
	polynome_n_fixe(std::vector<int> vec, T element); //monome
	polynome_n_fixe(int* vec, T element, T faux); //monome interne
	polynome_n_fixe(const polynome_n_fixe<T, 1>& copie) : poly(copie.poly), nul(copie.nul) { };

	polynome_n_fixe(polynome_n_fixe<T, 1>&& temp) {
		nul = temp.nul;
		swap(poly, temp.poly);
		return;
	}


	polynome_n_fixe<T, 1>& operator=(const polynome_n_fixe<T, 1>& copie) {
		if (this == &copie)
			return *this;
		nul = copie.nul;
		poly = copie.poly;
		return *this;
	};


	polynome_n_fixe<T, 1>& operator=(polynome_n_fixe<T, 1>&& temp) {
		nul = temp.nul;
		swap(poly, temp.poly);
		return *this;
	};

	polynome_n_fixe<T, 1>& operator*=(const polynome_n_fixe<T, 1>& temp) {
		return (*this = (*this * temp));
	};

	template<class U>
	polynome_n_fixe<T, 1>& operator*=(U const& scalaire) {
		poly *= scalaire;
		nul = (bool)poly;
		return *this;
	};

	polynome_n_fixe<T, 1>& operator+=(polynome_n_fixe<T, 1> const& autre) {
		poly += autre.poly;
		nul = (bool)poly;
		return *this;
	};

	friend polynome_n_fixe<T, 1> operator*(const polynome_n_fixe<T, 1>& gauche, const polynome_n_fixe<T, 1>& droit) {
		polynome_n_fixe<T, 1> m_poly;
		m_poly.poly = gauche.poly * droit.poly;
		m_poly.nul = (bool)m_poly.poly;

		return m_poly;

	};

	friend polynome_n_fixe<T, 1> operator+(const polynome_n_fixe<T, 1>& gauche, const polynome_n_fixe<T, 1>& droit) {
		polynome_n_fixe<T, 1> m_poly;
		m_poly.poly = gauche.poly + droit.poly;
		m_poly.nul = (bool)m_poly.poly;
		return m_poly;
	};

	friend polynome_n_fixe<T, 1> operator-(const polynome_n_fixe<T, 1>& gauche, const polynome_n_fixe<T, 1>& droit) {
		polynome_n_fixe<T, 1> m_poly;
		m_poly.poly = gauche.poly * droit.poly;
		m_poly.nul = (bool)m_poly.poly;
		return m_poly;
	};

	friend polynome_n_fixe<T, 1> operator-(const polynome_n_fixe<T, 1>& element) {
		polynome_n_fixe<T, 1> m_poly(element);
		m_poly.poly = -m_poly.poly;
		return m_poly;
	};

	friend polynome_n_fixe<T, 1> operator/(const polynome_n_fixe<T, 1>& gauche, const polynome_n_fixe<T, 1>& droit) {
		static_assert(type_algebre<T>::type == 0);

		polynome_n_fixe<T, 1> m_poly;
		m_poly.poly = gauche.poly / droit.poly;
		m_poly.nul = (bool)m_poly.poly;

		return m_poly;
	};

	friend polynome_n_fixe<T, 1> operator%(const polynome_n_fixe<T, 1>& gauche, const polynome_n_fixe<T, 1>& droit) {
		static_assert(type_algebre<T>::type == 0);

		polynome_n_fixe<T, 1> m_poly;
		m_poly.poly = gauche.poly % droit.poly;
		m_poly.nul = (bool)m_poly.poly;

		return m_poly;
	};

	template<class U> friend polynome_n_fixe<T, 1> operator*(const U& scalaire, const polynome_n_fixe<T, 1>& polynome_) {
		polynome_n_fixe<T, 1> m_poly(polynome_);
		m_poly.poly = scalaire * m_poly.poly;
		m_poly.nul = (bool)m_poly.poly;

		return m_poly;
	};

	friend bool operator==(polynome_n_fixe<T, 1> const& gauche, polynome_n_fixe<T, 1> const& droit) {
		return gauche.poly == droit.poly;
	};


	T get_T() const {
		return poly.coeffs[0].get_T();
	};

	explicit operator bool() const {
		return nul;
	};

	friend void swap(polynome_n_fixe<T, 1>& gauche, polynome_n_fixe<T, 1>& droit) {
		std::swap(gauche.nul, droit.nul);
		swap(gauche.poly, droit.poly);
		return;
	};

	int length() {
		return poly.coeffs.size();
	};


	template<class I>
	class iterator_polynome_n_fixe;

	using iterator = iterator_polynome_n_fixe<T>;
	using const_iterator = iterator_polynome_n_fixe<const T>;

	iterator begin() {
		return iterator(*this);
	};

	iterator end() {
		iterator it(*this);
		it.termine = false;
		return it;
	};


	const_iterator cbegin() const {
		return const_iterator(*this);
	};

	const_iterator cend() const {
		const_iterator it(*this);
		it.termine = false;
		return it;
	};


	template<class U>
	operator polynome_n_fixe<U, 1>() {
		polynome_n_fixe<U, 1> result();
		result.poly = (polynome<polynome_n_fixe<U, 0>>) poly;
		result.nul = (bool)result.poly;

		return result;
	};

	operator polynome_n_rec<T>() {
		polynome_n_rec<T> result();
		result.n_var = 1;
		result.element = unite(get_T(), false);
		result.nul = nul;
		result.poly = (polynome<polynome_n_rec<T>>) poly;

		return result;
	};

	int max_degre(int i) const {
#ifdef ALGEBRA_USE_EXCEPTION
		if (i != 0)
			throw std::domain_error("polynome_n_fixe : get_degre : i >= n");
#endif
		return poly.degre;
	}

	std::vector<int> max_degre() const {
		std::vector<int> result{ poly.degre };
		return result;
	};

	void max_degre(int* p) const {
		if (poly.degre > p[0])
			p[0] = poly.degre;
		return;
	};

	int lenght() {
		return poly.coeffs.size() + 1;
	}
};


template<class T, int n> polynome_n_fixe<T, n>::polynome_n_fixe(std::vector<int> vec, T element) {
#ifdef ALGEBRA_USE_EXCEPTION
	if (vec.size() != n)
		throw std::domain_error("constructeur de polynome_n_fixe : n ne correspond pas");
	for (int i(0); i < vec.size(); ++i)
		if (vec[i] < 0) {
//			poly = polynome<polynome_n_fixe<T, n - 1>>(polynome_n_fixe<T, n - 1>(unite(element, false)));
//			nul = false;
//			return;
			throw std::domain_error("constructeur de polynome_n_fixe : un degre est negatif");
		}
#endif
	if (!(bool)element) {
		poly = polynome<polynome_n_fixe<T, n - 1>>(polynome_n_fixe<T, n - 1>(element));
		nul = false;
		return;
	}
	nul = true;
	int* data = vec.data() + 1;
	int p = vec[0];
	T faux = unite(element, false);

	polynome_n_fixe<T, n - 1> sous_poly(faux);
	std::vector< polynome_n_fixe<T, n - 1>> vec_poly(p + 1, sous_poly);
	vec_poly[p] = polynome_n_fixe<T, n - 1>(data, element, faux);
	poly = polynome<polynome_n_fixe<T, n - 1>>(vec_poly);
};

template<class T, int n> polynome_n_fixe<T, n>::polynome_n_fixe(int* vec, T element, T faux) {  //monome interne
	int p = *vec;
	nul = true;

	polynome_n_fixe<T, n - 1> sous_poly(faux);
	std::vector< polynome_n_fixe<T, n - 1>> vec_poly(p + 1, sous_poly);
	vec_poly[p] = polynome_n_fixe<T, n - 1>(vec + 1, element, faux);
	poly = polynome<polynome_n_fixe<T, n - 1>>(vec_poly);
};

template<class T> polynome_n_fixe<T, 1>::polynome_n_fixe(std::vector<int> vec, T element) {
#ifdef ALGEBRA_USE_EXCEPTION
	if (vec.size() != 1)
		throw std::domain_error("constructeur de polynome_n_fixe : n ne correspond pas");
	if (vec[0] < 0) {
//			poly = polynome<polynome_n_fixe<T, 0>>(polynome_n_fixe<T, 0>(unite(element, false)));
//			nul = false;
//			return;
			throw std::domain_error("constructeur de polynome_n_fixe : un degre est negatif");
		}
#endif
	if (!(bool)element) {
		poly = polynome<polynome_n_fixe<T, 0>>(polynome_n_fixe<T, 0>(element));
		nul = false;
		return;
	}
	nul = true;
	int p = vec[0];
	T faux = unite(element, false);

	polynome_n_fixe<T, 0> sous_poly(faux);
	std::vector< polynome_n_fixe<T, 0>> vec_poly(p + 1, sous_poly);
	vec_poly[p] = polynome_n_fixe<T, 0>(element);
	poly = polynome<polynome_n_fixe<T, 0>>(vec_poly);
};


template<class T>
polynome_n_fixe<T, 1>::polynome_n_fixe(int* vec, T element, T faux) {  //monome interne
	int p = *vec;
	nul = true;

	polynome_n_fixe<T, 0> sous_poly(faux);
	std::vector< polynome_n_fixe<T, 0>> vec_poly(p + 1, sous_poly);
	vec_poly[p] = polynome_n_fixe<T, 0>(element);
	poly = polynome<polynome_n_fixe<T, 0>>(vec_poly);
};


template<class T> template<class I>
class polynome_n_fixe<T, 1>::iterator_polynome_n_fixe {
public:

	using value_type = std::remove_const_t<I>;
	using poly_type = std::conditional_t< std::is_const_v<I>, const polynome_n_fixe<value_type, 1>, polynome_n_fixe<value_type, 1>	>;

	poly_type* pointeur;

	int p;
	bool termine;

	iterator_polynome_n_fixe(poly_type& temp) : pointeur(&temp), termine(true), p(0) { };

	iterator_polynome_n_fixe(iterator_polynome_n_fixe const& copie) {
		p=copie.p;
		termine = copie.termine;
	}

	iterator_polynome_n_fixe& operator=(iterator_polynome_n_fixe const& copie) {
		p=copie.p;
		termine = copie.termine;
		return *this;
	}

	iterator_polynome_n_fixe& operator++() {
		++p;
		if (p >= pointeur->poly.coeffs.size()) {
			p = 0;
			termine = false;
		};
		return *this;
	};

	operator bool() const {
		return termine;
	};

	friend bool operator==(iterator_polynome_n_fixe const& it1, iterator_polynome_n_fixe const& it2) {
		if ((it1.termine != it2.termine) || (it1.pointeur != it2.pointeur) || (it1.p != it2.p))
			return false;
		else
			return true;
	}

	friend bool operator!=(iterator_polynome_n_fixe const& it1, iterator_polynome_n_fixe const& it2) {
		if (it1.termine != it2.termine)
			return true;
		if (it1.pointeur != it2.pointeur)
			return true;
		if (it1.p != it2.p)
			return true;
		return false;
	};

	void go_position(int j) {
		p = j;
		return;
	};

	I& operator*() {
		return pointeur->poly.coeffs[p].element;
	};

};


template<class T> class polynome_n_fixe<T, 0> {
public:
	T element;
	bool nul;

	polynome_n_fixe() {};
	polynome_n_fixe(T element_) : element(element_), nul((bool)element_) {};
	polynome_n_fixe(std::vector<int> vec, T element_) : element(element_), nul((bool)element_) {};  //monome
	polynome_n_fixe(int* vec, T element_, T faux) : element(element_), nul((bool)element_) {};  //monome interne
	polynome_n_fixe(const polynome_n_fixe<T, 0>& copie) : element(copie.element), nul(copie.nul) {	};
	polynome_n_fixe(polynome_n_fixe<T, 0>&& temp) {
		nul = temp.nul;
		swap_F(element, temp.element);
		return;
	};

	polynome_n_fixe<T, 0>& operator=(const polynome_n_fixe<T, 0>& copie){
		if (this == &copie)
			return *this;
		element = copie.element;
		nul = copie.nul;
		return *this;
	};

	polynome_n_fixe<T, 0>& operator=(polynome_n_fixe<T, 0>&& temp) {
		nul = temp.nul;
		swap_F(element, temp.element);
		return *this;
	};

	polynome_n_fixe<T, 0>& operator*=(const polynome_n_fixe<T, 0>& temp) {
		element *= temp.element;
		nul = (bool)element;
		return *this;
	};

	template<class U>
	polynome_n_fixe<T, 0>& operator*=(U const& scalaire) {
		element *= scalaire;
		nul = (bool)element;
		return *this;
	};

	polynome_n_fixe<T, 0>& operator+=(polynome_n_fixe<T, 0> const& autre) {
		element += autre.element;
		nul = (bool)element;
		return *this;
	}

	friend polynome_n_fixe<T, 0> operator*(const polynome_n_fixe<T, 0>& gauche, const polynome_n_fixe<T, 0>& droit) {
		return polynome_n_fixe<T, 0>(gauche.element * droit.element);
	};

	friend polynome_n_fixe<T, 0> operator+(const polynome_n_fixe<T, 0>& gauche, const polynome_n_fixe<T, 0>& droit) {
		return polynome_n_fixe<T, 0>(gauche.element + droit.element);
	};

	friend polynome_n_fixe<T, 0> operator-(const polynome_n_fixe<T, 0>& gauche, const polynome_n_fixe<T, 0>& droit) {
		return polynome_n_fixe<T, 0>(gauche.element - droit.element);
	};

	friend polynome_n_fixe<T, 0> operator-(const polynome_n_fixe<T, 0>& elem) {
		return polynome_n_fixe<T, 0>(- elem.element);
	}

	friend polynome_n_fixe<T, 0> operator/(const polynome_n_fixe<T, 0>& gauche, const polynome_n_fixe<T, 0>& droit) {
		static_assert((type_algebre<T>::type == 1) || (type_algebre<T>::type == 0));
		return polynome_n_fixe<T, 0>(gauche.element / droit.element);
	};

	friend polynome_n_fixe<T, 0> operator%(const polynome_n_fixe<T, 0>& gauche, const polynome_n_fixe<T, 0>& droit) {
		static_assert(type_algebre<T>::type == 1);
		return polynome_n_fixe<T, 0>(gauche.element % droit.element);
	};

	template<class U> friend polynome_n_fixe<T, 0> operator*(const U& scalaire, const polynome_n_fixe<T, 0>& poly) {
		return polynome_n_fixe<T, 0>(scalaire * poly.element);
	};

	friend bool operator==(polynome_n_fixe<T, 0> const& gauche, polynome_n_fixe<T, 0> const& droit) {
		return gauche.element == droit.element;
	}


	T get_T() const {
		return element;
	};

	explicit operator bool() const {
		return nul;
	};


	friend void swap(polynome_n_fixe<T, 0>& gauche, polynome_n_fixe<T, 0>& droit) {
		std::swap(gauche.nul, droit.nul);
		swap(gauche.element, droit.element);
		return;
	};

	int length() {
		return 1;
	};

	template<class U>
	operator polynome_n_fixe<U, 0>() {
		polynome_n_fixe<U, 0> result((U) element);
		result.nul = (bool) result.poly;

		return result;
	}

	operator polynome_n_rec<T>() {
		polynome_n_rec<T> result();
		result.n_var = 0;
		result.element = element;
		result.nul = nul;

		return result;
	};


};

