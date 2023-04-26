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
		poly.coeffs.resize(0);
		nul = (bool)element;
	};

	polynome_n_rec(polynome_n_rec<T> const& copie) : nul(copie.nul), element(copie.element), noms_variables(copie.noms_variables),n_var(copie.n_var),poly(copie.poly) {	};
	polynome_n_rec(polynome_n_rec<T>&& copie) : nul(copie.nul), noms_variables(copie.noms_variables), n_var(copie.n_var) {
		swap_F(element, copie.element);
		swap_F(poly, copie.poly);
		return;
	};

	polynome_n_rec(int n_var_, std::string* noms_, T element_, std::vector<polynome_n_rec<T>> const& vec) : n_var(n_var_),noms_variables(noms_),poly(vec),element(element_){ //à partir d'un vecteur de pointeurs. Utilisé pour construire (* et +)
		nul = false;
		for (int i(poly.coeffs.size() - 1); i >= 0; --i)
			if ((bool)poly.coeffs[i])
				nul = true;
		return;
	};

	polynome_n_rec(int n_var_, std::string* noms_, T element_, std::vector<polynome_n_rec<T>> && vec) : n_var(n_var_), noms_variables(noms_), poly(vec), element(element_) { //à partir d'un vecteur de pointeurs. Utilisé pour construire (* et +)
		nul = false;
		for (int i(poly.coeffs.size() - 1); i >= 0; --i)
			if ((bool)poly.coeffs[i])
				nul = true;
		return;
	};


	polynome_n_rec(std::vector<int> liste, std::string* noms, const T element_) { //liste : degrés ... Monome.
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
			}

		if (!(bool)element_) {
			nul = false;
			poly = { poly_faux };
		}

		std::vector<polynome_n_rec<T>> m_vec(liste[0] + 1, poly_faux);
		m_vec[liste[0] + 1] = polynome_n_rec<T>(liste.data() + 1, n_var - 1, noms_variables + 1, element_, element);

		poly = polynome<polynome_n_rec<T>>(m_vec);

//		polynome_n_rec(int* liste, int n, std::string * noms, T const& element_, T const& faux
		return;
	};

	polynome_n_rec(int* liste, int n, std::string* noms, T const& element_, T const& faux) {
		n_var = n;
		if (n_var > 0) {
			noms_variables = noms;
			nul = true; //vérifié précédemment
			element = faux;

			polynome_n_rec<T> poly_faux(polynome_n_rec<T>(n_var - 1, element, noms + 1));

			std::vector<polynome_n_rec<T>> m_vec(liste[0] + 1, poly_faux);
			m_vec[liste[0] + 1] = polynome_n_rec<T>(liste + 1, n_var - 1, noms_variables + 1, element_, element);

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


	friend polynome_n_rec<T> operator*(polynome_n_rec<T> const& temp1, polynome_n_rec<T> const& temp2) {
		if (temp1.n_var != temp2.n_var)
			throw std::domain_error("polynomes à n variables de tailles différentes");
		
		if (temp1.n_var ==0)
			return polynome_n_rec<T>(temp1.element * temp2.element);

		polynome_n_rec<T> result;
		result.poly = temp1.poly * temp2.poly;
		result.noms_variables = temp1.noms_variables;
		result.element = temp1.element;
		result.n_var = temp1.n_var;

		result.nul = false;
		for (int i(result.poly.coeffs.size() - 1); i >= 0; --i)
			if ((bool)result.poly.coeffs[i])
				result.nul = true;

		return result;
	}

	friend polynome_n_rec<T> operator+(polynome_n_rec<T> const& temp1, polynome_n_rec<T> const& temp2) {
		if (temp1.n_var != temp2.n_var)
			throw std::domain_error("polynomes à n variables de tailles différentes");

		if (temp1.n_var == 0)
			return polynome_n_rec<T>(temp1.element + temp2.element);

		polynome_n_rec<T> result;
		result.poly = temp1.poly + temp2.poly;
		result.noms_variables = temp1.noms_variables;
		result.element = temp1.element;
		result.n_var = temp1.n_var;

		result.nul = false;
		for (int i(result.poly.coeffs.size() - 1); i >= 0; --i)
			if ((bool)result.poly.coeffs[i])
				result.nul = true;

		return result;

	};


	friend polynome_n_rec<T> operator-(polynome_n_rec<T> const& temp1, polynome_n_rec<T> const& temp2) {
		if (temp1.n_var != temp2.n_var)
			throw std::domain_error("polynomes à n variables de tailles différentes");

		if (temp1.n_var == 0)
			return polynome_n_rec<T>(temp1.element - temp2.element);

		polynome_n_rec<T> result;
		result.poly = temp1.poly - temp2.poly;
		result.noms_variables = temp1.noms_variables;
		result.element = temp1.element;
		result.n_var = temp1.n_var;

		result.nul = false;
		for (int i(result.poly.coeffs.size() - 1); i >= 0; --i)
			if ((bool)result.poly.coeffs[i])
				result.nul = true;

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
				result.nul = false;
				for (int i(result.poly.coeffs.size() - 1); i >= 0; --i)
					if ((bool)result.poly.coeffs[i])
						result.nul = true;
				result.noms_variables = temp.noms_variables;
				result.element = temp.element;
				result.n_var;
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

		result.nul = false;
		for (int i(result.poly.coeffs.size() - 1); i >= 0; --i)
			if ((bool)result.poly.coeffs[i])
				result.nul = true;
		result.n_var = n_var;
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

	void simplifier() { //reparcourt tout, met à jour nul.
		if (n_var == 0) {
			nul = (bool)element;
			return nul;
		}
		for (int i(0); i < poly.coeffs.size(); ++i)
			poly.coeffs[i].simplifier();

		int trouve_int = 0;
		bool trouve = false;
		for (int i(poly.coeffs.size() - 1); i >= 0; --i) {
			if (poly.coeffs[i].nul) {
				trouve_int = i;
				trouve = true;
				break;
			}
		}

		if (trouve_int + 1 < poly.coeffs.size())
			poly.coeffs.erase(poly.coeffs.begin() + trouve_int + 1, poly.coeffs.end());

		nul = trouve;
	};

	friend void swap(polynome_n_rec<T>& gauche, polynome_n_rec<T>& droit) {
		std::swap(gauche.n_var, droit.n_var);
		std::swap(gauche.noms_variables, droit.noms_variables);
		std::swap(gauche.nul, droit.nul);
		swap_F(gauche.element, droit.element);
		swap_F(gauche.poly, droit.poly);
		return;
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
