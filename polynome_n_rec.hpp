#pragma once
#include <vector>
#include <string>
#include <exception>
#include "polynome_n_iter.hpp"
#include "polynome.hpp"

#include "entete objets.hpp"
#include "unite.hpp"
#include "swap_T.hpp"


template<class T>  T unite(T const& element, bool test);

//polynome à n variables

//template<class T> void swap(T& gauche, T& droit);

template<class T> class polynome_n_rec {
public:

	T element;
	int n_var;
	bool nul; //true si non-nul

	std::string* noms_variables;
	polynome<polynome_n_rec<T>> poly;


	polynome_n_rec() {};

	polynome_n_rec(int n, T temp, std::string* noms) { //constructeur de base. polynome vide.
		n_var = n;
		nul = (bool) temp;
		if (n > 0) {
			noms_variables = noms;
			element = unite(temp,false);
			poly = polynome<polynome_n_rec<T>>(polynome_n_rec<T>(n - 1, temp, noms+1) );
		}
		else {
			element = temp;
			noms_variables = NULL;
		}
	};

	polynome_n_rec(T element_) { //crée rapidement un polynome de base. Utilisé lors de la récurrence.
		element = element_;
		n_var = 0;
		noms_variables = NULL;
		nul = (bool)element;
	};

	polynome_n_rec(polynome_n_rec<T> const& copie) : nul(copie.nul), element(copie.element), noms_variables(copie.noms_variables),n_var(copie.n_var),poly(copie.poly) {	};
	polynome_n_rec(polynome_n_rec<T>&& copie) : nul(copie.nul), noms_variables(copie.noms_variables), n_var(copie.n_var) {
		swap_F(element, copie.element);
		swap_F(poly, copie.poly);
		return;
	};

	polynome_n_rec(int n_var_, std::string* noms_, T element_, std::vector<polynome_n_rec<T>> const& vec) : n_var(n_var_),noms_variables(noms_),poly(vec),element(element_){ //à partir d'un vecteur de pointeurs. Utilisé pour construire (* et +)
		nul = (bool) poly;
		return;
	};

	polynome_n_rec(int n_var_, std::string* noms_, T element_, std::vector<polynome_n_rec<T>> && vec) : n_var(n_var_), noms_variables(noms_), poly(vec), element(element_) { //à partir d'un vecteur de pointeurs. Utilisé pour construire (* et +)
		nul = (bool)poly;
		return;
	};


	polynome_n_rec(std::vector<int> liste, std::string* noms, const T& element_) { //liste : degrés ... Monome.
		noms_variables = noms;
		n_var = liste.size();

		if (n_var == 0) {
			nul = (bool)element_;
			element = element_;
			return;
		}
		element = unite(element_, false);

		polynome_n_rec<T> poly_faux(polynome_n_rec<T>(n_var - 1, element, noms + 1));

		for(int i(0);i<liste.size();++i)
			if (liste[i] < 0) {
				nul = false;
				poly = { poly_faux };
				return;
			}

		if (!(bool)element_) {
			nul = false;
			poly = { poly_faux };
			return;
		}

		std::vector<polynome_n_rec<T>> m_vec(liste[0] + 1, poly_faux);
		m_vec[liste[0] + 1] = polynome_n_rec<T>(liste.data() + 1, n_var - 1, noms_variables + 1, element_, element);

		poly = polynome<polynome_n_rec<T>>(m_vec);
		nul = (bool) poly;

		return;
	};

	polynome_n_rec(int* liste, int n, std::string* noms, T const& element_, T const& faux) {
		n_var = n;
		if (n_var > 0) {
			noms_variables = noms;
			nul = true; //vérifié précédemment
			element = faux;

			polynome_n_rec<T> poly_faux(polynome_n_rec<T>(n_var - 1, faux, noms + 1));

			std::vector<polynome_n_rec<T>> m_vec(liste[0] + 1, poly_faux);
			m_vec[liste[0] + 1] = polynome_n_rec<T>(liste + 1, n_var - 1, noms_variables + 1, element_, faux);

			poly = polynome<polynome_n_rec<T>>(m_vec);
		}
		else {
			element = element_;
			nul = true;
			noms_variables = NULL;
		}
		return;
	};

	~polynome_n_rec() {	};

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
	}

	friend std::ostream& operator<<(std::ostream& os, polynome_n_rec<T> const& element) {
		std::string result = "";
		element.getString("", result);
		os << result;
		return os;
	};


	polynome_n_rec<T>& operator=(polynome_n_rec<T> const& temp) {
		if (this == &temp)
			return *this;

		nul = temp.nul;
		noms_variables = temp.noms_variables;
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
		if (temp1.n_var ==0)
			return polynome_n_rec<T>(temp1.element * temp2.element);
#endif
		polynome_n_rec<T> result;
		result.poly = temp1.poly * temp2.poly;
		result.noms_variables = temp1.noms_variables;
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
		result.noms_variables = temp1.noms_variables;
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
		result.noms_variables = temp1.noms_variables;
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
		result.noms_variables = temp.noms_variables;
		result.element = temp.element;
		result.n_var = temp.n_var;

		result.nul = temp.nul;

		return result;
	};


	void changer_nom(std::string* noms) {
		noms_variables = noms;
		if (n_var > 1)
			for (int i(0); i < poly.coeffs.size(); ++i)
				poly.coeffs[i].changer_nom(noms + 1);
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
				result.noms_variables = temp.noms_variables;
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

	template<class U> operator polynome_n_rec<U>() {
		if (n_var == 0)
			return polynome_n_rec<U>((U)element);
		polynome_n_rec<U> result;
		result.element = (U)element;
		result.noms_variables = noms_variables;
		result.poly = (polynome<polynome_n_rec<U>>) poly;
		result.n_var = n_var;

		result.nul = (bool)result.poly;

		return result;
	};

	std::vector<int> getDegre_tab() const { //ATTENTION : degrés et non dimensions !
		std::vector<int> tab(n_var, -1);
		getDegre(tab.data());
		return tab;
	};

	void getDegre(int* degres) const {
		int i = 0;
		for (i = poly.coeffs.size() - 1; i >= 0; --i)
			if (poly.coeffs[i].nul)
				break;

		if (*degres < i)
			*degres = i;
		if (n_var == 1)
			return;

		for (int j(i); j >= 0; --j)
			if (poly.coeffs[j].nul)
				poly.coeffs[j].getDegre(degres + 1);
		return;
	};

	operator polynome_n_iter<T>() const {
		if (n_var == 0)
			return polynome_n_iter<T>(0, element, NULL);
		if (!nul)
			return polynome_n_iter<T>(n_var, unite(element, false), noms_variables);

		std::vector<int> degres = getDegre_tab();
		T faux = unite(element,false);

		polynome_n_iter<T> result(degres, faux, noms_variables); //degrés ... se transforme en dimensions (+1).
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
		std::swap(gauche.noms_variables, droit.noms_variables);
		std::swap(gauche.nul, droit.nul);
		swap_F(gauche.element, droit.element);
		swap_F(gauche.poly, droit.poly);
		return;
	};

	class iterator {
	public:
		int n_var;

		std::vector<int> positions;
		std::vector<polynome_n_rec<T>*> pointeurs;
		bool termine;


		iterator& operator++() {
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

		T operator*() const {
			if (n_var > 0)
				return pointeurs[n_var - 1]->poly.coeffs[positions[n_var - 1]].element;
			else
				return pointeurs[0]->element;
		};

		iterator(polynome_n_rec<T> const& temp) : positions(temp.n_var, 0), n_var(temp.n_var), pointeurs(temp.n_var) ,termine(true){
			if (n_var > 0) {
				pointeurs[0] = &temp;
				for (int i(1); i < n_var; ++i)
					pointeurs[i] = &pointeurs[i - 1]->poly.coeffs[0];
			}
			else {
				pointeurs = { &temp };
			};
		};
	};

	iterator begin() {
		return iterator(*this);
	};

	iterator end() {
		iterator it(*this);
		it.termine = false;
		return it;
	};

	int length() {
		if (n_var == 0)
			return 1;
		if (n_var == 1)
			return poly.coeffs.size();
		else {
			int somme = 0;
			for (int i(0); i < poly.coeffs.size(); ++i)
				somme += poly.coeffs[i].length;
			return somme;
		};
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
			parcourir_convert(& objet->poly.coeffs[i], data + (i * (puissances[0])), puissances + 1);
	}
	return;
}
