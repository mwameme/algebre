#pragma once
#include <vector>
#include <string>
#include <exception>
#include "polynome.hpp"

#include "entete objets.hpp"

template<class T, int n> class polynome_n_fixe;

template<class T, int n, int m> bool operator==(polynome_n_fixe<T, n> const& gauche, polynome_n_fixe<T, m> const& droit) {
	return false;
};

template<class T, int n> class polynome_n_fixe {
public:

	polynome<polynome_n_fixe<T, n - 1>> poly;
	bool nul;

	polynome_n_fixe() {};
	polynome_n_fixe(T element) : poly(polynome_n_fixe<T, n - 1>(element)) ,nul((bool) element) {};
	polynome_n_fixe(std::vector<int> vec, T element); //monome. vec : degrés
	polynome_n_fixe(int* vec, T element, T faux);  //monome interne
	polynome_n_fixe(const polynome_n_fixe<T, n>& copie) : poly(copie.poly), nul(copie.nul) {	};

	polynome_n_fixe(polynome_n_fixe<T, n>&& temp) {
		nul = temp.nul;
		swap(poly, temp.poly);
		return;
	}

	polynome_n_fixe<T, n>& operator=(const polynome_n_fixe<T, n>& copie) {
		if (this == &copie)
			return *this;
		nul = copie.nul;
		poly = copie.poly;
		return *this;
	};

	polynome_n_fixe<T, n>& operator=(const polynome_n_fixe<T, n>&& temp) {
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
		nul = false;
		for (int i(poly.coeffs.size() - 1); i >= 0; --i)
			if (poly.coeffs[i].nul)
				nul = true;
		return *this;
	};

	friend polynome_n_fixe<T, n> operator*(const polynome_n_fixe<T, n>& gauche, const polynome_n_fixe<T, n>& droit) {
		polynome_n_fixe<T, n> m_poly;
		m_poly.poly = gauche.poly * droit.poly;
		m_poly.nul = false;
		for (int i(m_poly.poly.coeffs.size() - 1); i >= 0; --i)
			if (m_poly.poly.coeffs[i].nul)
				m_poly.nul = true;
		return m_poly;
	};

	friend polynome_n_fixe<T, n> operator+(const polynome_n_fixe<T, n>& gauche, const polynome_n_fixe<T, n>& droit) {
		polynome_n_fixe<T, n> m_poly;
		m_poly.poly = gauche.poly + droit.poly;
		m_poly.nul = false;
		for (int i(m_poly.poly.coeffs.size() - 1); i >= 0; --i)
			if (m_poly.poly.coeffs[i].nul)
				m_poly.nul = true;
		return m_poly;
	};

	friend polynome_n_fixe<T, n> operator-(const polynome_n_fixe<T, n>& gauche, const polynome_n_fixe<T, n>& droit) {
		polynome_n_fixe<T, n> m_poly;
		m_poly.poly = gauche.poly - droit.poly;
		m_poly.nul = false;
		for (int i(m_poly.poly.coeffs.size() - 1); i >= 0; --i)
			if (m_poly.poly.coeffs[i].nul)
				m_poly.nul = true;
		return m_poly;
	};

	friend polynome_n_fixe<T, n> operator-(const polynome_n_fixe<T, n>& element) {
		polynome_n_fixe<T, n> m_poly(element);
		m_poly.poly = -m_poly.poly;
		return m_poly;
	};

	template<class U> friend polynome_n_fixe<T, n> operator*(const U& scalaire, const polynome_n_fixe<T, n>& polynome_) {
		polynome_n_fixe<T, n> m_poly(polynome_);
		m_poly.poly = scalaire * m_poly.poly;
		m_poly.nul = false;
		for (int i(m_poly.poly.coeffs.size() - 1); i >= 0; --i)
			if (m_poly.poly.coeffs[i].nul)
				m_poly.nul = true;
		return m_poly;
	};


	friend bool operator==(polynome_n_fixe<T, n>& gauche, polynome_n_fixe<T, n>& droit) {
		return gauche.poly == droit.poly;
	};

	T get_T() const {
		return poly.coeffs[0].get_T();
	};

	explicit operator bool() {
		return nul;
	};

	friend void swap(polynome_n_fixe<T, n>& gauche, polynome_n_fixe<T, n>& droit) {
		std::swap(gauche.nul, droit.nul);
		swap(gauche.poly, droit.poly);
		return;
	};

};



template<class T> class polynome_n_fixe<T, 1> {
public:
	polynome<polynome_n_fixe<T, 0>> poly;
	bool nul;

	polynome_n_fixe() {};
	polynome_n_fixe(T element) : poly(polynome_n_fixe<T,0>(element)), nul((bool)element) {};
	polynome_n_fixe(std::vector<int> vec, T element); //monome
	polynome_n_fixe(int* vec, T element, T faux); //monome interne
	polynome_n_fixe(const polynome_n_fixe<T, 1>& copie) : poly(copie.poly), nul(copie.nul) {	};

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


	polynome_n_fixe<T, 1>& operator=(const polynome_n_fixe<T, 1>&& temp) {
		nul = temp.nul;
		swap(poly, temp.poly);
		return *this;
	};

	polynome_n_fixe<T, 1>& operator*=(const polynome_n_fixe<T, 1>& temp) {
		return (*this = (*this * temp));
	};

	template<class U>
	polynome_n_fixe<T, n>& operator*=(U const& scalaire) {
		poly *= scalaire;
		nul = false;
		for (int i(poly.coeffs.size() - 1); i >= 0; --i)
			if (poly.coeffs[i].nul)
				nul = true;
		return *this;
	};


	friend polynome_n_fixe<T, 1> operator*(const polynome_n_fixe<T, 1>& gauche, const polynome_n_fixe<T, 1>& droit) {
		polynome_n_fixe<T, 1> m_poly;
		m_poly.poly = gauche.poly * droit.poly;
		m_poly.nul = false;
		for (int i(m_poly.poly.coeffs.size() - 1); i >= 0; --i)
			if (m_poly.poly.coeffs[i].nul)
				m_poly.nul = true;
		return m_poly;

	};

	friend polynome_n_fixe<T, 1> operator+(const polynome_n_fixe<T, 1>& gauche, const polynome_n_fixe<T, 1>& droit) {
		polynome_n_fixe<T, 1> m_poly;
		m_poly.poly = gauche.poly + droit.poly;
		m_poly.nul = false;
		for (int i(m_poly.poly.coeffs.size() - 1); i >= 0; --i)
			if (m_poly.poly.coeffs[i].nul)
				m_poly.nul = true;
		return m_poly;
	};

	friend polynome_n_fixe<T, 1> operator-(const polynome_n_fixe<T, 1>& gauche, const polynome_n_fixe<T, 1>& droit) {
		polynome_n_fixe<T, 1> m_poly;
		m_poly.poly = gauche.poly * droit.poly;
		m_poly.nul = false;
		for (int i(m_poly.poly.coeffs.size() - 1); i >= 0; --i)
			if (m_poly.poly.coeffs[i].nul)
				m_poly.nul = true;
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
		m_poly.nul = false;
		for (int i(m_poly.poly.coeffs.size() - 1); i >= 0; --i)
			if (m_poly.poly.coeffs[i].nul)
				m_poly.nul = true;
		return m_poly;
	};

	friend polynome_n_fixe<T, 1> operator%(const polynome_n_fixe<T, 1>& gauche, const polynome_n_fixe<T, 1>& droit) {
		static_assert(type_algebre<T>::type == 0);
		polynome_n_fixe<T, 1> m_poly;
		m_poly.poly = gauche.poly % droit.poly;
		m_poly.nul = false;
		for (int i(m_poly.poly.coeffs.size() - 1); i >= 0; --i)
			if (m_poly.poly.coeffs[i].nul)
				m_poly.nul = true;
		return m_poly;
	};

	template<class U> friend polynome_n_fixe<T, 1> operator*(const U& scalaire, const polynome_n_fixe<T, 1>& polynome_) {
		polynome_n_fixe<T, 1> m_poly(polynome_);
		m_poly.poly = scalaire * m_poly.poly;
		m_poly.nul = false;
		for (int i(m_poly.poly.coeffs.size() - 1); i >= 0; --i)
			if (m_poly.poly.coeffs[i].nul)
				m_poly.nul = true;
		return m_poly;
	};

	friend bool operator==(polynome_n_fixe<T, 1>& gauche, polynome_n_fixe<T, 1>& droit) {
		return gauche.poly == droit.poly;
	};


	T get_T() const {
		return poly.coeffs[0].get_T();
	};

	explicit operator bool() {
		return nul;
	};

	friend void swap(polynome_n_fixe<T, 1>& gauche, polynome_n_fixe<T, 1>& droit) {
		std::swap(gauche.nul, droit.nul);
		swap(gauche.poly, droit.poly);
		return;
	};
};


template<class T, int n> polynome_n_fixe<T, n>::polynome_n_fixe(std::vector<int> vec, T element) {
	if (vec.size() != n)
		throw std::domain_error("constructeur de polynome_n_fixe : n ne correspond pas");
	for (int i(0); i < vec.size(); ++i)
		if (vec[i] < 0) {
			poly = polynome<polynome_n_fixe<T, n - 1>>(polynome_n_fixe<T, n - 1>(unite(element, false)));
			nul = false;
			return;
		}
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


template<class T> class polynome_n_fixe<T, 0> {
public:
	T element;
	bool nul;

	polynome_n_fixe() {};
	polynome_n_fixe(T element_) : element(element_), nul((bool)element_) {};
	polynome_n_fixe(std::vector<int> vec, T element_) : element(element_), nul((bool) element_);  //monome
	polynome_n_fixe(int* vec, T element_, T faux) : element(element_), nul((bool)element_);  //monome interne
	polynome_n_fixe(const polynome_n_fixe<T, 0>& copie) : element(copie.element), nul(copie.nul) {	};
	polynome_n_fixe(polynome_n_fixe<T, 0>&& temp) {
		nul = temp.nul;
		swap(element, temp.element);
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
		swap(element, temp.element);
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

	friend bool operator==(polynome_n_fixe<T, 0>& gauche, polynome_n_fixe<T, 0>& droit) {
		return gauche.element == droit.element;
	}


	T get_T() const {
		return element;
	};

	explicit operator bool() {
		return nul;
	};


	friend void swap(polynome_n_fixe<T, 0>& gauche, polynome_n_fixe<T, 0>& droit) {
		std::swap(gauche.nul, droit.nul);
		swap(gauche.poly, droit.poly);
		return;
	};
};

