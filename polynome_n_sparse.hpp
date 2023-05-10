#pragma once

#include <vector>
#include <algorithm>
#include <exception>
#include <string>

#include "swap_T.hpp"
#include "unite.hpp"
#include "fonctions template.hpp"


template<class T> class polynome_n_rec;
template<class T> class polynome_n_iter;
template<class T, int n> class polynome_n_fixe;

template<class T> class monome {
public:
	std::vector<int> degres;
	T element;

	monome() {};

	monome(std::vector<int> const& degres_, T const& element_) : degres(degres_), element(element_) {
#ifdef ALGEBRA_USE_EXCEPTION
		for (int i(0); i < degres.size(); ++i)
			if (degres[i] < 0)
				throw std::domain_error("constructeur de monome : un degre est negatif");
#endif
	}

	monome(monome<T> const& copie) {
		degres = copie.degres;
		element = copie.element;
		return;
	};

	monome(monome<T>&& copie) {
		swap(*this, copie);
	};

	monome<T>& operator=(monome<T> const& copie) {
		degres = copie.degres;
		element = copie.element;
		return *this;
	};

	monome<T>& operator=(monome<T>&& copie) {
		if (this == &copie)
			return *this;
		swap(*this, copie);
		return *this;
	};

	friend void swap(monome<T>& gauche, monome<T>& droit) {
		swap(gauche.degres, droit.degres);
		swap_F(gauche.element, droit.element);
		return;
	};

	template<class U> 
	monome<T>& operator*=(U const& scalaire) {
		element *= scalaire;
		return *this;
	};

	template<class U>
	friend monome<T> operator*(U const& scalaire, monome<T> const& monome_) {
		return monome<T>(monome_.degres, scalaire * monome_.element);
	}


	friend monome<T> operator*(monome<T> const& gauche, monome<T> const& droit) {
		T temp = gauche.element * droit.element;
		if ((bool)temp) {
			std::vector<int> degres_temp = gauche.degres;
			for (int i(0); i < degres_temp.size(); ++i)
				degres_temp[i] += droit.degres[i];

			return monome<T>(degres_temp, temp);
		}
		std::vector<int> degres_temp(gauche.degres.size(), 0);
		return monome<T>(degres_temp, temp);
	};

	monome<T> operator-() const {
		return monome<T>(degres, -element);
	};

	friend bool operator<(monome<T> const& gauche, monome<T> const& droit) {
		return gauche.degres < droit.degres;
	};

	friend bool operator>(monome<T> const& gauche, monome<T> const& droit) {
		return gauche.degres > droit.degres;
	};

	friend monome<T> operator-(monome<T> const& temp) {
		return monome<T>(temp.degres, -temp.element);
	};

	template<class U>
	operator monome<U>() const {
		return monome<U>(degres, (U)element);
	};
};

template<class T> class polynome_n_sparse {
public:

	int n_var;
	std::vector<monome<T>> monomes;
	bool est_trie; //si true, la liste est triee ...

	using iterator = typename std::vector<monome<T>>::iterator;
	using const_iterator = typename std::vector<monome<T>>::const_iterator;

	polynome_n_sparse(int n) : n_var(n), est_trie(true) {};
	polynome_n_sparse(monome<T> const& temp) : n_var(temp.degres.size()), est_trie(true), monomes({ temp }) {	};
	polynome_n_sparse(int n, T const& element) : n_var(n), est_trie(true), monomes({ monome<T>(std::vector<int>(n,0),element) }) {	};
	polynome_n_sparse(polynome_n_sparse<T> const& copie) : monomes(copie.monomes), n_var(copie.n_var), est_trie(copie.est_trie) {	};
	polynome_n_sparse(polynome_n_sparse<T> && copie) : n_var(copie.n_var), est_trie(copie.est_trie) {
		swap(monomes, copie.monomes);
	};

	polynome_n_sparse<T>& operator=(polynome_n_sparse<T> const& copie) {
		if (this == &copie)
			return *this;

		monomes = copie.monomes;
		n_var = copie.n_var;
		est_trie = copie.est_trie;
		return *this;
	};

	polynome_n_sparse<T>& operator=(polynome_n_sparse<T> && copie) {
		swap(monomes, copie.monomes);

		n_var = copie.n_var;
		est_trie = copie.est_trie;
		return *this;
	};

	polynome_n_sparse<T>& operator+=(polynome_n_sparse<T> const& autre) {
		return (*this = (*this + autre));
	};

	polynome_n_sparse<T>& operator*=(polynome_n_sparse<T> const& autre) {
		return (*this = (*this * autre));
	};

	template<class U>
	polynome_n_sparse<T>& operator*=(U const& scalaire) {
		if((bool) scalaire)
			for (int i(0); i < monomes.size(); ++i)
				monomes[i].element *= scalaire;
		else {
			monome<T> faux = unite(monomes[0], false);
			monomes = { faux };
			est_trie = true;
		}
		return *this;
	}


	friend void swap(polynome_n_sparse<T>& gauche, polynome_n_sparse<T>& droit) {
		std::swap(gauche.monomes, droit.monomes);
		std::swap(gauche.n_var, droit.n_var);
		std::swap(gauche.est_trie, droit.est_trie);
		return;
	};

	void ajouter(monome<T> const& temp) {
#ifdef ALGEBRA_USE_EXCEPTION
		if (temp.degres.size() != n_var)
			throw std::domain_error("polynome_n_sparse : ajout de monome : n_var ne correspond pas");
#endif
		monomes.push_back(temp);
		est_trie = false;
		return;
	};

	void ajouter(polynome_n_sparse<T> const& temp) {
		ajouter(temp.monomes);
	};

	void ajouter(std::vector<monome<T>> const& temp) {
#ifdef ALGEBRA_USE_EXCEPTION
		for(int i(0);i<temp.size();++i)
			if (temp[i].degres.size() != n_var)
				throw std::domain_error("polynome_n_sparse : ajout de monome : n_var ne correspond pas");
#endif
		monomes.insert(monomes.end(), temp.begin(), temp.end());
		est_trie = false;
		return;
	};

	void simplifier() {
#ifdef ALGEBRA_USE_EXCEPTION
		if (monomes.size() == 0)
			throw std::domain_error("polynome_n_sparse vide");
#endif
		if (est_trie)
			return;
		monome<T> faux = unite(monomes[0], false);

		std::vector<monome<T>*> vec(monomes.size(), NULL);
		for (int i(0); i < vec.size(); ++i)
			vec[i] = &monomes[i];
		std::sort(vec.begin(), vec.end(), [](monome<T>* g, monome<T>* d) {return (g->degres) < (d->degres); }); //ordre croissant

		std::vector<monome<T>> vec_2(0);
		vec_2.reserve(vec.size());
		int j = 0;
		vec_2.push_back(*vec[0]);
		for (int i(1); i < vec.size(); ++i)
			if (vec[i]->degres == vec_2[j].degres)
				vec_2[j].element += vec[i]->element;
			else {
				if (!(bool)vec_2[j].element) {
					vec_2[j] = * (vec[i]);
				}
				else {
					++j;
					vec_2.push_back(*(vec[i]));
				}
			}
		if (!(bool)vec_2[j].element)
			vec_2.pop_back();

		if (vec_2.size() == 0)
			vec_2.push_back(faux);

		est_trie = true; //trie
		std::swap(monomes, vec_2);
		return;
	};

	
	template<class U>
	friend polynome_n_sparse<T> operator*(U const& scalaire, polynome_n_sparse<T> const& poly) {
		if ((bool)scalaire) {
			polynome_n_sparse<T> result(poly);
			for (int i(0); i < result.monomes.size(); ++i)
				result.monomes[i].element *= scalaire;
			return result;
		}
		return polynome_n_sparse<T>(n_var, unite(scalaire, false));
	};
	
	friend polynome_n_sparse<T> operator*(monome<T> const& monome_, polynome_n_sparse<T> const& poly) {
		if (!(bool)monome_.element)
			return polynome_n_sparse<T>(n_var, monome_.element);
		if ((bool)monome_.element) {
			polynome_n_sparse<T> result(poly);
			for (int i(0); i < result.monomes.size(); ++i) {
				for (int j(0); j < result.n_var; ++j)
					result.monomes[i].degres[j] += monome_.degres[j];
				result.monomes[i].element *= monome_.element;
			}
			return result;
		}
		else
			return polynome_n_sparse(n_var, unite(monomes[0], false));
	};

	friend inline polynome_n_sparse<T> operator*(polynome_n_sparse<T> const& poly, monome<T> const& monome_) {
		return monome_* poly;
	};


	friend polynome_n_sparse<T> operator*(polynome_n_sparse<T> const& gauche, polynome_n_sparse<T> const& droit) {
#ifdef ALGEBRA_USE_EXCEPTION
		if ( (gauche.monomes.size() == 0) || (droit.monomes.size() == 0)
			throw std::domain_error("polynome_n_sparse vide");
#endif
		polynome_n_sparse<T> result(gauche.n_var,gauche.noms);
		result.monomes.reserve(gauche.monomes.size() * droit.monomes.size());
		for (int i(0); i < gauche.monomes.size(); ++i)
			for (int j(0); j < droit.monomes.size(); ++j)
				result.monomes.push_back(gauche.monomes[i] * droit.monomes[j]);
		result.est_trie = false;
		result.simplifier(); //n^2 log(n^2)
		return result;
	};

	polynome_n_sparse<T> operator-() const {
		polynome_n_sparse<T> result(*this);
		for (int i(0); i < result.monomes.size(); ++i)
			result.monomes[i].element = -result.monomes[i].element;
		result.est_trie = est_trie; //conserve l'ordre
		return result;
	};

	friend polynome_n_sparse<T> operator+(polynome_n_sparse<T> const& gauche, polynome_n_sparse<T> const& droit) {
		if ((!gauche.est_trie) || (!droit.est_trie)) {
			polynome_n_sparse<T> result(gauche);
			result.ajouter(droit.monomes);
			result.simplifier();
			return result;
		}
		if (gauche.n_var == 0) {
			polynome_n_sparse<T> result(monome<T>(std::vector<int>(0), gauche.monomes[0].element + droit.monomes[0].element));
			result.noms = gauche.noms;
			return result;
		}

		int ig = 0; //gauche droit et final
		int id = 0;
		polynome_n_sparse<T> result(gauche.n_var, gauche.noms);
		result.monomes.reserve(gauche.monomes.size() + droit.monomes.size());
		monome<T> faux = unite(monomes[0], false);

		while (true) {
			if (gauche.monomes[ig].degres == droit.monomes[id].degres) {
				T temp = gauche.monomes[ig].element + droit.monomes[id].element;
				if ((bool) temp)
					result.monomes.push_back(monome<T>(gauche.monomes[ig].degres, temp));
				++ig;
				++id;
			}
			else if (gauche.monomes[ig].degres < droit.monomes[id].degres) {
				result.monomes.push_back(gauche.monomes[ig]);
				++ig;
			}
			else {
				result.monomes.push_back(droit.monomes[id]);
				++id;
			}

			if (ig == gauche.monomes.size()) {
				result.monomes.insert(result.monomes.end(), droit.monomes.begin() + id, droit.monomes.end());
				break;
			};
			if (id == droit.monomes.size()) {
				result.monomes.insert(result.monomes.end(), gauche.monomes.begin() + ig, gauche.monomes.end());
				break;
			};
		};
		result.est_trie = true;
		if (result.monomes.size() == 0)
			result.monomes.push_back(faux);
		result.noms = gauche.noms;
		return result;
	};

	friend polynome_n_sparse<T> operator-(polynome_n_sparse<T> const& gauche, polynome_n_sparse<T> const& droit) {
		if ((!gauche.est_trie) || (!droit.est_trie)) {
			polynome_n_sparse<T> result(gauche);
			int i = result.monomes.size();
			result.ajouter(droit.monomes);
			for (; i < result.monomes.size(); ++i)
				result.monomes[i].element = -result.monomes[i].element;
			result.simplifier();
			return result;
		}
		if (gauche.n_var == 0) {
			polynome_n_sparse<T> result(monome<T>(std::vector<int>(0), gauche.monomes[0].element - droit.monomes[0].element));
			result.noms = gauche.noms;
			return result;
		}

		int ig = 0; //gauche droit et final
		int id = 0;
		polynome_n_sparse<T> result(gauche.n_var, gauche.noms);
		result.monomes.reserve(gauche.monomes.size() + droit.monomes.size());
		monome<T> faux = unite(monomes[0], false);

		while (true) {
			if (gauche.monomes[ig].degres == droit.monomes[id].degres) {
				T temp = gauche.monomes[ig].element - droit.monomes[id].element;
				if ((bool)temp)
					result.monomes.push_back(monome<T>(gauche.monomes[ig].degres, temp));
				++ig;
				++id;
			}
			else if (gauche.monomes[ig].degres < droit.monomes[id].degres) {
				result.monomes.push_back(gauche.monomes[ig]);
				++ig;
			}
			else {
				result.monomes.push_back(- droit.monomes[id]);
				++id;
			}

			if (ig == gauche.monomes.size()) { //on ajoute la fin de -droit ...
				result.monomes.reserve(result.monomes.size() + droit.monomes.size() - id);
				for (int i(id); i < droit.monomes.size(); ++i)
					result.monomes.push_back(-droit.monomes[i]);
//				result.monomes.insert(result.monomes.end(), droit.monomes.begin() + id, droit.monomes.end());
				break;
			};
			if (id == droit.monomes.size()) { //on ajoute la fin de gauche
				result.monomes.insert(result.monomes.end(), gauche.monomes.begin() + ig, gauche.monomes.end());
				break;
			};
		};
		result.est_trie = true;
		if (result.monomes.size() == 0)
			result.monomes.push_back(faux);

		result.noms = gauche.noms;
		return result;
	};

	int max_degre(int i) const {
#ifdef ALGEBRA_USE_EXCEPTION
		if ((i < 0) || (i >= n_var))
			throw std::domain_error("polynome_n_sparse : i > n_var");
#endif
		int result = 0;
		for (int j(0); j < monomes.size(); ++j)
			result = max(result, monomes[j].degres[i]);
		return result;
	};

	std::vector<int> max_degre() const {
		std::vector<int> result(n_var);
		for (int j(0); j < monomes.size(); ++j)
			for (int i(0); i < n_var; ++i)
				result[i] = max(result[i], monomes[j].degres[i]);
		return result;
	};

	polynome_n_sparse<T> operator()(int i, T const& Y) const {
		polynome_n_sparse<T> poly(monome<T>(std::vector<int>(n_var - 1, 0), Y));
		return *this(i, poly);
	};

	polynome_n_sparse<T> operator()(int i, polynome_n_sparse<T> const& poly) const {
		polynome_n_sparse<T> result = unite(poly, false);
		int m_pow = max_degre(i);
#ifdef ALGEBRA_USE_EXCEPTION
		if (poly.n_var != (n_var - 1))
			throw std::domain_error("evaluation de polynome n_sparse : les degr�s ne correspondent pas");
#endif
		std::vector< polynome_n_sparse<T>> vec_pow;
		vec_pow.reserve(m_pow + 1);
		vec_pow.push_back(unite(poly, true));
		for (int j(1); j < vec_pow.size(); ++j)
			vec_pow[j] = vec_pow[j - 1] * poly;

		result.monomes.reserve(poly.monomes.size() * vec_pow[vec_pow.size() - 1].monomes.size());

		for (int j(0); j < monomes.size(); ++j) {
			std::vector<int> temp_degres = monomes[j].degres;
			int pow = temp_degres[i];
			temp_degres.erase(temp_degres.begin() + i);

			//			polynome_n_sparse<T> temp_poly(monome<T>(temp_degres, monomes[j].element));

//			polynome_n_sparse<T> temp_poly(vec_pow[pow]); //On fait � la main et on simplifie seulement � la fin : + rapide
			polynome_n_sparse<T> temp_poly = monome<T>(temp_degres, monomes[j].element) * vec_pow[pow];
			result.ajouter(temp_poly);  //result += monome<T>(temp_degres,monomes.element) * vec_pow[pow] 
		};

		result.simplifier();
		return result;
		//marche avec r�sultat scalaire ? je crois
	};

	
	template<class U>
	U operator()(int i, U const& poly) const {
		U result = unite(poly, false);
		int m_pow = max_degre(i);
		std::vector< U> vec_pow;
		vec_pow.reserve(m_pow + 1);
		vec_pow.push_back(unite(poly, true));
		for (int j(1); j < vec_pow.size(); ++j)
			vec_pow[j] = vec_pow[j - 1] * poly;

//		result.monomes.reserve(poly.monomes.size() * vec_pow[vec_pow.size() - 1].monomes.size());

		for (int j(0); j < monomes.size(); ++j) {
			std::vector<int> temp_degres = monomes[j].degres;
			int pow = temp_degres[i];
			temp_degres.erase(temp_degres.begin() + i);

			//			polynome_n_sparse<T> temp_poly(monome<T>(temp_degres, monomes[j].element));

			U temp_poly(vec_pow[pow]); //On fait � la main et on simplifie seulement � la fin : + rapide
			temp_poly = monome<T>(temp_degres, monomes[j].element) * temp_poly;
			result += temp_poly;  //result += monome<T>(temp_degres,monomes.element) * vec_pow[pow] 
		};

		return result;
		//marche avec r�sultat scalaire ? je crois
	};
	

	T get_T() const {
		return monomes[0].element;
	}

	operator polynome_n_rec<T>() const {
		polynome_n_rec<T> poly(n_var, unite(get_T(), false));
		for (int i(0); i < monomes.size(); ++i)
			poly += polynome_n_rec<T>(monomes[i].degres, monomes[i].element);

		return poly;
	};

	operator polynome_n_iter<T>() {
		polynome_n_iter<T> result(max_degre(), unite(get_T(), false));
		if (!est_trie)
			simplifier();
		for (int i(0); i < monomes.size(); ++i) {
			result.coeffs.data[result.coeffs.position(monomes[i].degres)] = monomes[i].element;
		};
		return result;
	};

	template<int n>
	operator polynome_n_fixe<T,n>() const {
#ifdef ALGEBRA_USE_EXCEPTION
		if (n != n_var)
			throw std::domain_error("convertisseur de polynome_n_fixe : n ne correspond pas");
#endif
		polynome_n_fixe<T,n> poly(n_var, unite(get_T(), false));
		for (int i(0); i < monomes.size(); ++i)
			poly += polynome_n_fixe<T, n>(monomes[i].degres, monomes[i].element);

		return poly;
	};

	operator bool() const {
#ifdef ALGEBRA_USE_EXCEPTION
		if (!est_trie)
			throw std::domain_error("polynome_n_sparse, (bool) : n'est pas trie");
		if (monomes.size() == 0)
			throw std::domain_error("polynome_n_sparse, (bool) : vide");
#endif
		return (bool) monomes[0].element;
	};

	template<class U>
	operator polynome_n_sparse<U>() {
		polynome_n_sparse<U> poly(n_var);
		poly.monomes.reserve(monomes.size());
		for (int i(0); i < monomes.size(); ++i)
			poly.monomes.push_back((monome<U>) monomes[i]);
		poly.est_trie = est_trie; //faire attention : si des �l�ments deviennent nuls ...
		return poly;
	}

	typename iterator begin() {
		return monomes.begin();
	};

	typename iterator end() {
		return monomes.end();
	}

	typename const_iterator cbegin()  const {
		return monomes.cbegin();
	};

	typename const_iterator cend() const {
		return monomes.cend();
	};

	friend std::ostream& operator<<(std::ostream& os, polynome_n_sparse<T> const& poly) {
		if ((bool)poly) {
			for (int i(0); i < monomes.size(); ++i) {
				os << poly.monomes[i].element << " *";
				for (int j(0); j < n_var; ++j)
					os << "X" << std::to_string(j) << "^" << std::to_string(monomes[i].degres[j]) << " ";
				if (i < monomes.size() - 1)
					os << " + ";
			}
			return os;
		}
		os << " 0 ";
		return os;
	};

};

