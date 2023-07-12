#pragma once
#include <vector>
#include <string>
#include <exception>
#include "polynome_n_iter.hpp"
#include "polynome.hpp"

#include "entete objets.hpp"
#include "unite.hpp"
#include "swap_T.hpp"

#include "polynome_n_fixe.hpp"
#include "polynome_n_sparse.hpp"

template<class T> class polynome_n_sparse;
template<class T> class monome;

template<class T, int n> class polynome_n_fixe;

template<class T>  T unite(T const& element, bool test);

//polynome à n variables

//template<class T> void swap(T& gauche, T& droit);

template<class T> class polynome_n_rec {
public:
	int n_var;
	polynome<polynome_n_rec<T>> poly;
	T element;
	bool nul; //true si non-nul

	using sous_type = typename T;

	polynome_n_rec() : n_var(0) ,nul(false) {};

	polynome_n_rec(int n, T temp) { //constructeur de base. polynome vide.
		n_var = n;
		nul = (bool) temp;
		if (n > 0) {
			element = unite(temp,false);
			poly = polynome<polynome_n_rec<T>>(polynome_n_rec<T>(n - 1, temp) );
		}
		else {
			element = temp;
		}
	};

	polynome_n_rec(T element_) { //crée rapidement un polynome de base. Utilisé lors de la récurrence.
		element = element_;
		n_var = 0;
		nul = (bool)element;
	};

	polynome_n_rec(polynome_n_rec<T> const& copie) : nul(copie.nul), element(copie.element),n_var(copie.n_var),poly(copie.poly) {	};
	polynome_n_rec(polynome_n_rec<T>&& copie) : nul(copie.nul), n_var(copie.n_var) {
		swap_F(element, copie.element);
		swap(poly, copie.poly);
		return;
	};

	polynome_n_rec(int n_var_, T element_, std::vector<polynome_n_rec<T>> const& vec) : n_var(n_var_), poly(vec), element(element_) { //à partir d'un vecteur de pointeurs. Utilisé pour construire (* et +)
		nul = (bool)poly;
		return;
	};

	polynome_n_rec(int n_var_, T element_, std::vector<polynome_n_rec<T>> && vec) : n_var(n_var_),  poly(vec), element(element_) { //à partir d'un vecteur de pointeurs. Utilisé pour construire (* et +)
		nul = (bool)poly;
		return;
	};

	polynome_n_rec(monome<T> monome_) : polynome_n_rec(monome_.degres, monome_.element) {};

	polynome_n_rec(std::vector<int> liste, const T& element_) { //liste : degrés ... Monome.
		n_var = liste.size();

		if (n_var == 0) {
			nul = (bool)element_;
			element = element_;
			return;
		}
		element = unite(element_, false);

		polynome_n_rec<T> poly_faux(polynome_n_rec<T>(n_var - 1, element));
#ifdef ALGEBRA_USE_EXCEPTION
		for (int i(0); i < liste.size(); ++i)
			if (liste[i] < 0)
				throw std::domain_error("creation de polynome_n_rec : un degre est negatif");
#endif
		if (!(bool)element_) {
			nul = false;
			poly = { poly_faux };
			return;
		}

		std::vector<polynome_n_rec<T>> m_vec(liste[0] + 1, poly_faux);
		m_vec[liste[0] + 1] = polynome_n_rec<T>(liste.data() + 1, n_var - 1, element_, element);

		poly = polynome<polynome_n_rec<T>>(m_vec);
		nul = (bool) poly;

		return;
	};

	polynome_n_rec(int* liste, int n, T const& element_, T const& faux) { //monome interne
		n_var = n;
		if (n_var > 0) {
			nul = true; //vérifié précédemment
			element = faux;

			polynome_n_rec<T> poly_faux(polynome_n_rec<T>(n_var - 1, faux));

			std::vector<polynome_n_rec<T>> m_vec(liste[0] + 1, poly_faux);
			m_vec[liste[0] + 1] = polynome_n_rec<T>(liste + 1, n_var - 1, element_, faux);

			poly = polynome<polynome_n_rec<T>>(m_vec);
		}
		else {
			element = element_;
			nul = true;
		}
		return;
	};

	~polynome_n_rec() {	};

	/*
	void getString(std::string puissance,std::string& resultat) const {
		if (!nul) //n'est possible, que lors du premier appel. Car sinon on teste avant.
			resultat = "0";
		else {
			if (n_var > 0) { //récurrence
				for (int i(0); i < poly.coeffs.size(); ++i)
					if (poly.coeffs[i].nul) {
						std::string puissance_locale = puissance + "*" + noms_variables[0] + "^" + std::to_string(i) + " ";
						poly.coeffs[i].getString(puissance_locale, resultat);
					}
			}
			else { //on modifie la chaine "résultat"
				std::string temp;
				std::stringstream ss;
				ss << element;
				temp=ss.str();
//				std::cout << "element "<< element << std::endl;
				if (resultat.size()>0)
					resultat = resultat + "+" + temp + puissance; // P_0 + (element) * puissance
				else
					resultat = temp + puissance;
			}
		}

		return;
	}*/

	friend std::ostream& operator<<(std::ostream& os, polynome_n_rec<T> const& element) {
		bool trouve = false;
		for (auto it = element.cbegin(); (bool)it; ++it) {
			if (!(bool)*it)
				continue;
			if (!trouve)
				trouve = true;
			else
				os << " + ";
			os << *it << "* ";
			for (int i(0); i < element.n_var; ++i)
				os << "X" << i << "^" << it.positions[i] << " ";
		}
		if (!trouve)
			os << " 0 ";
		return os;
	};


	polynome_n_rec<T>& operator=(polynome_n_rec<T> const& temp) {
		if (this == &temp)
			return *this;

		nul = temp.nul;
		element = temp.element;
		n_var = temp.n_var;
		poly = temp.poly;

		return *this;
	};

	polynome_n_rec<T>& operator=(polynome_n_rec<T> && temp) {
		if (this == &temp)
			return *this;
		swap(*this, temp);
		return *this;
	};

	polynome_n_rec<T>& operator*=(polynome_n_rec<T> const& temp) {
		return (*this = (*this * temp));
	};

	template<class U>
	polynome_n_rec<T>& operator*=(U const& scalaire) {
		if (n_var == 0) {
			element *= scalaire;
			nul = (bool)element;
		}
		else {
			poly *= scalaire;
			nul = (bool)poly;
		}
		return *this;
	};

	polynome_n_rec<T>& operator+=(polynome_n_rec<T> const& autre) {
#ifdef ALGEBRA_USE_EXCEPTION
		if (n_var != autre.n_var)
			throw std::domain_error("addition de polynome_n_rec : n_var ne correspondent pas");
#endif
		if (n_var == 0) {
			element += autre.element;
			return *this;
		}
		poly += autre.poly;
		nul = (bool)poly;
		return *this;
	};


	friend polynome_n_rec<T> operator*(polynome_n_rec<T> const& temp1, polynome_n_rec<T> const& temp2) {
#ifdef ALGEBRA_USE_EXCEPTION
		if (temp1.n_var != temp2.n_var)
			throw std::domain_error("polynomes à n variables de tailles différentes");
#endif
		if (temp1.n_var == 0)
			return polynome_n_rec<T>(temp1.element * temp2.element);

		polynome_n_rec<T> result;
		result.poly = temp1.poly * temp2.poly;
		result.element = temp1.element;
		result.n_var = temp1.n_var;

		result.nul = (bool)result.poly;

		return result;
	}

	friend polynome_n_rec<T> operator+(polynome_n_rec<T> const& temp1, polynome_n_rec<T> const& temp2) {
#ifdef ALGEBRA_USE_EXCEPTION
		if (temp1.n_var != temp2.n_var)
			throw std::domain_error("polynomes à n variables de tailles différentes");
#endif
		if (temp1.n_var == 0)
			return polynome_n_rec<T>(temp1.element + temp2.element);

		polynome_n_rec<T> result;
		result.poly = temp1.poly + temp2.poly;
		result.element = temp1.element;
		result.n_var = temp1.n_var;

		result.nul = (bool)result.poly;

		return result;

	};


	friend polynome_n_rec<T> operator-(polynome_n_rec<T> const& temp1, polynome_n_rec<T> const& temp2) {
#ifdef ALGEBRA_USE_EXCEPTION
		if (temp1.n_var != temp2.n_var)
			throw std::domain_error("polynomes à n variables de tailles différentes");
#endif
		if (temp1.n_var == 0)
			return polynome_n_rec<T>(temp1.element - temp2.element);
		polynome_n_rec<T> result;
		result.poly = temp1.poly - temp2.poly;
		result.element = temp1.element;
		result.n_var = temp1.n_var;

		result.nul = (bool)result.poly;

		return result;
	};

	friend polynome_n_rec<T> operator-(polynome_n_rec<T> const& temp) { //nul ne change pas. Donc on garde tout tel quel, le seul changement est : element -> - element. Et recopiage.

		if (temp.n_var == 0)
			return polynome_n_rec<T>(-temp.element);

		polynome_n_rec<T> result;
		result.poly = - temp.poly;
		result.element = temp.element;
		result.n_var = temp.n_var;

		result.nul = temp.nul;

		return result;
	};

	friend bool operator==(polynome_n_rec<T> const& temp1, polynome_n_rec<T> const& temp2) {
		if (temp1.n_var != temp2.n_var)
			return false;
		if (temp1.n_var > 0) {
			return (temp1.poly == temp2.poly);
		}
		else
			return (temp1.element == temp2.element);
	};

	template<class U> friend polynome_n_rec<T> operator*(U const& scalaire, polynome_n_rec<T> const& temp) {
		if (!(bool)temp) {
			return temp;
		}
		if ((bool)scalaire) {
			if (temp.n_var > 0) {
				polynome_n_rec<T> result;
				result.poly = scalaire * temp.poly;
				result.element = temp.element;
				result.n_var = temp.n_var;

				result.nul = (bool) result.poly;

				return result;
			}
			else 
				return polynome_n_rec<T>(scalaire * temp.element);
		} 
		else  //retourne le polynome nul.
			return polynome_n_rec<T>(temp.n_var, unite(temp.element, false), temp.noms_variables);
	};

	explicit inline operator bool() const {
		return nul;
	};

	template<class U>
	operator polynome_n_rec<U>() {
		if (n_var == 0)
			return polynome_n_rec<U>((U)element);
		polynome_n_rec<U> result;
		result.element = (U)element;
		result.poly = (polynome<polynome_n_rec<U>>) poly;
		result.n_var = n_var;

		result.nul = (bool)result.poly;

		return result;
	};

	int max_degre(int i) const {
		if (i == 0)
			return poly.degre;
		int result = -1;
		for (int j(0); j < poly.coeffs.size(); ++j)
			result = max(result, poly.coeffs[i].max_degre(i - 1));
		return result;
	}


	std::vector<int> max_degre() const { //ATTENTION : degrés et non dimensions !
		std::vector<int> tab(n_var, -1);
		max_degre(tab.data());
		return tab;
	};

	void max_degre(int* degres) const {
		int i = 0;

		if (*degres < poly.degre)
			*degres = poly.degre;
		if (n_var == 1)
			return;

		for (int j(0); j <poly.coeffs.size() ; ++j)
			if (poly.coeffs[j].nul)
				poly.coeffs[j].max_degre(degres + 1);
		return;
	};

	operator polynome_n_iter<T>() const {
		if (n_var == 0)
			return polynome_n_iter<T>(0, element);
		if (!nul)
			return polynome_n_iter<T>(n_var, unite(element, false));

		std::vector<int> degres = max_degre();
		T faux = unite(element,false);

		polynome_n_iter<T> result(degres, faux); //degrés ... se transforme en dimensions (+1).
		if (result.coeffs.data.size() == 1) {
			if (! nul)
				result.coeffs.data[0] = unite(result.coeffs.data[0],true);
		}
		else
			parcourir_convert(this, result.coeffs.data.data(), result.coeffs.puissances.data());
		return result;
	};

	void simplifier() { //reparcourt tout, met à jour nul, réduit les polynomes.
		if (n_var == 0) {
			nul = (bool)element;
			return nul;
		}
		for (int i(0); i < poly.coeffs.size(); ++i)
			poly.coeffs[i].simplifier();

		poly.getDegre();

		nul = (bool) poly;
	};

	friend void swap(polynome_n_rec<T>& gauche, polynome_n_rec<T>& droit) {
		std::swap(gauche.n_var, droit.n_var);
		std::swap(gauche.nul, droit.nul);
		swap_F(gauche.element, droit.element);
		swap(gauche.poly, droit.poly);
		return;
	};

	template<class I>
	class iterator_polynome_n_rec;

	using iterator = iterator_polynome_n_rec<T>;
	using const_iterator = iterator_polynome_n_rec<const T>;


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

	const_iterator end() const {
		const_iterator it(*this);
		it.termine = false;
		return it;
	};


	int length() { //taille en mémoire
		if (n_var == 0)
			return 1;
		if (n_var == 1)
			return poly.coeffs.size()+1;
		else {
			int somme = 0;
			for (int i(0); i < poly.coeffs.size(); ++i)
				somme += poly.coeffs[i].length();
			return somme+1; //+ cet espace local
		};
	};

	operator polynome_n_sparse<T>() const {
		polynome_n_sparse<T> poly(monome<T>(std::vector<int>(n_var, 0), unite(element, false)));

		for (iterator it = cbegin(); (bool)it; ++it)
			if ((bool)*it)
				poly.ajouter(monome<T>(it.positions, *it));
		poly.simplifier();//non nécessaire ...
	}

	
	polynome_n_rec<T> operator()(int i, polynome_n_rec<T> const& poly) const {
		polynome_n_rec<T> result = unite(poly, false);
		int m_pow = max_degre(i);
#ifdef ALGEBRA_USE_EXCEPTION
		if (poly.n_var != (n_var - 1))
			throw std::domain_error("evaluation de polynome n_sparse : les degrés ne correspondent pas");
#endif
		std::vector< polynome_n_rec<T>> vec_pow;
		vec_pow.reserve(m_pow + 1);
		vec_pow.push_back(unite(poly, true));
		for (int j(1); j < vec_pow.size(); ++j)
			vec_pow[j] = vec_pow[j - 1] * poly;

		for (iterator it = cbegin(); (bool) it; ++it) {
			if (!(bool)*it)
				continue;
			std::vector<int> temp_degres = it.positions;
			int pow = temp_degres[i];
			temp_degres.erase(temp_degres.begin() + i);

			//			polynome_n_sparse<T> temp_poly(monome<T>(temp_degres, monomes[j].element));

			polynome_n_rec<T> temp_poly(vec_pow[pow]); //On fait à la main et on simplifie seulement à la fin : + rapide
			polynome_n_rec<T> monome_(n_var - 1, temp_degres, *it);
			result += monome_ * temp_poly;
		};

		result.noms_variables = poly.noms_variables;
		return result;
		//marche avec résultat scalaire ? je crois
	};


	template<class U>
	U operator()(int i, U const& poly) const {
		U result = unite(poly, false);
		int m_pow = max_degre(i);
#ifdef ALGEBRA_USE_EXCEPTION
		if (poly.n_var != (n_var - 1))
			throw std::domain_error("evaluation de polynome n_sparse : les degrés ne correspondent pas");
#endif
		std::vector< U> vec_pow;
		vec_pow.reserve(m_pow + 1);
		vec_pow.push_back(unite(poly, true));
		for (int j(1); j < vec_pow.size(); ++j)
			vec_pow[j] = vec_pow[j - 1] * poly;

		for (iterator it = cbegin(); (bool)it; ++it) {
			if (!(bool)*it)
				continue;
			std::vector<int> temp_degres = it.positions;
			int pow = temp_degres[i];
			temp_degres.erase(temp_degres.begin() + i);

			//			polynome_n_sparse<T> temp_poly(monome<T>(temp_degres, monomes[j].element));

//			polynome_n_rec<T> temp_poly(vec_pow[pow]); //On fait à la main et on simplifie seulement à la fin : + rapide
			polynome_n_rec<T> monome_(n_var - 1, temp_degres, *it);
			result += monome_ * vec_pow[pow];
		};

		return result;
		//marche avec résultat scalaire ? je crois
	};

	polynome_n_rec<T> operator()(int i, T const& Y) const {
		polynome_n_rec<T> poly(n_var-1,Y);
		return *this(i, poly);
	};

};

//vérifier nul=true/false

template<class T> void parcourir_convert(polynome_n_rec<T> const* objet, T* data, int* puissances) {
	if (objet->n_var == 0) {
		data[0] = objet->element;
		return;
	}
	for (int i(objet->poly.coeffs.size() - 1); i >= 0; --i) {
		if (objet->poly.coeffs[i].nul)
			parcourir_convert(&objet->poly.coeffs[i], data + (i * (puissances[0])), puissances + 1);
	}
	return;
};

template<class T> template<class I >
class polynome_n_rec<T>::iterator_polynome_n_rec{
public:
	int n_var;
	using value_type = std::remove_const_t<I>;
	using poly_type = std::conditional_t< std::is_const_v<I>, const polynome_n_rec<value_type>, polynome_n_rec<value_type>	>;

	std::vector<int> positions;
	std::vector<poly_type *> pointeurs;
	bool termine;

	std::vector<int> get_position() {
		return positions;
	};


	iterator_polynome_n_rec& operator++() {
		if (n_var == 0) {
			termine = false;
			return *this;
		};

		int i = n_var - 1;
		++positions[i];
		while (positions[i] >= pointeurs[i]->poly.coeffs.size()) {
			positions[i] = 0;
			if (i == 0) {
				termine = false;
				return *this;
			}
			++positions[i - 1];
			--i;
		}
		for (int j = i + 1; j < n_var; ++j)
			pointeurs[j] = &pointeurs[j - 1]->poly.coeffs[positions[j - 1]];
		return *this;
	};

	bool operator==(iterator const& autre) const {
		return (pointeurs[0] == autre.pointeurs[0]) && (positions == autre.positions) && (termine = autre.termine);
	};

	bool operator!=(iterator const& autre) const {
		if (termine != autre.termine)
			return true;
		if (pointeurs[0] != autre.pointeurs[0])
			return true;
		if (n_var != 0)
			return (positions != autre.positions);
		else
			return false;
	};

	I& operator*() const {
		if (n_var > 0)
			return pointeurs[n_var - 1]->poly.coeffs[positions[n_var - 1]].element;
		else
			return pointeurs[0]->element;
	};

	iterator_polynome_n_rec(poly_type& temp) : positions(temp.n_var, 0), n_var(temp.n_var), pointeurs(temp.n_var, NULL), termine(true) {
		if (n_var > 0) {
			pointeurs[0] = &temp;
			for (int i(1); i < n_var; ++i)
				pointeurs[i] = &pointeurs[i - 1]->poly.coeffs[0];
		}
		else {
			pointeurs = { &temp };
		};
	};

	operator bool() const {
		return termine;
	};
};
