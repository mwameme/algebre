#pragma once
#include <vector>
#include <string>
#include <exception>
#include "polynome_n_iter.hpp"

#include "entete objets.hpp"
#include "unite.hpp"

//polynome à n variables

template<class T> class polynome_n {
public:

	polynome_n() {};

	polynome_n(int n, T temp, std::string* noms) { //constructeur de base. polynome vide.
		n_var = n;
		nul = (bool) temp;
		if (n > 0) {
			noms_variables = noms;
			element = unite(temp,false);
			coeffs = { new polynome_n<T>(n - 1, temp, noms+1) };
		}
		else {
			element = temp;
			noms_variables = NULL;
			coeffs.resize(0);
		}
	};

	/*
	polynome_n(int n, T temp, std::string* noms,int n_puissance,int puissance) : coeffs(0){//constructeur de X^n, où X est la n_puissance variable. et n la puissance.
		if (!(bool)temp) { //si c'est nul ... on a le modele de base avec temp=false. suppose n>0.
			n_var = n;
			nul = false;
			noms_variables = noms;
			element = temp;
			element = false;
			coeffs = { new polynome_n<T>(n - 1, temp, noms + 1) };
		}

		n_var = n;
		if (n > 0) {
			element = temp;
			element = false;
			if (n != n_puissance) {
				nul = true;
				noms_variables = noms;
				coeffs = { new polynome_n<T>(n - 1, temp, noms+1,n_puissance,puissance) }; //récurrence
			}
			else { //à partir de ce moment, on utilise le premier constructeur
				noms_variables = noms;
				T faux = temp;
				faux = false;
				nul = true;
				coeffs.reserve(puissance + 1); //puissance
				for(int i(0);i<puissance;++i)
					coeffs.push_back(new polynome_n<T>(n - 1, faux, noms+1));
				coeffs.push_back(new polynome_n<T>(n - 1, temp, noms+1));
				coeffs.reserve(0);
			}
		}
		else {
			element = temp;
			nul = true;
			noms_variables = NULL;
			coeffs.resize(0);
		}
	};
	*/

	polynome_n(T element_) { //crée rapidement un polynome de base. Utilisé lors de la récurrence.
		element = element_;
		n_var = 0;
		noms_variables = NULL;
		coeffs.resize(0);
		nul = (bool)element;
	};

	polynome_n(polynome_n<T> const& poly) {
		if (poly.n_var == 0) {
			element = poly.element;
			n_var = 0;
			noms_variables = NULL; // = NULL
			coeffs.resize(0);
			nul = poly.nul;
		}
		else {
			coeffs.resize(poly.coeffs.size());
			for (int i(0); i < coeffs.size(); ++i)
				coeffs[i] = new polynome_n<T>(*poly.coeffs[i]); //récurrence de la recopie.

			element = poly.element;
			n_var = poly.n_var;
			noms_variables = poly.noms_variables;
			nul = poly.nul;
		}
	};

	polynome_n(int n_var_, std::string* nom_, T element_, bool nul_, std::vector<polynome_n<T>*> vec) { //à partir d'un vecteur de pointeurs. Utilisé pour construire (* et +)
		n_var = n_var_;
		noms_variables = nom_;
		coeffs = vec;
		element = element_; //false
		nul = nul_;
	};

	polynome_n(std::vector<int> liste, std::string* noms, const T element_) : coeffs(0){ //liste : degrés ... Monome.
		bool test = true;
		if (!(bool)element)
			test = false;
		for (int i(0); i < liste.size(); ++i)
			if (liste[i] < 0) {
				test = false;
				break;
			}

		if (liste.size() == 0) {
			element = element_;
			nul = test;
			noms_variables = NULL;
			n_var = 0;
		}

		if (!test) { //constructeur de polynome "vide".
			element = unite(element_,false);
			nul = false;
			n_var = liste.size();
			if (n_var > 0) {
				noms_variables = noms;
				coeffs.push_back(new polynome_n<T>(n_var - 1, element, noms + 1));
				return;
			}
			noms_variables = NULL;
			return;
		}

		
		element = unite(element_,false);
		nul = true;
		n_var = liste.size();
		noms_variables = noms;
		coeffs.resize(liste[0] + 1);
		for (int i(0); i < liste[0]; ++i) {
			coeffs[i] = new polynome_n<T>(n_var - 1, element, noms + 1);
		}
		coeffs[liste[0]] =new polynome_n<T>(liste.data() + 1, n_var - 1, noms + 1, element_);
		return;
	};

	polynome_n(int* liste, int n, std::string* noms, T const element_) : coeffs(0) {
		n_var = n;
		if (n_var > 0) {
			noms_variables = noms;
			nul = true; //vérifié précédemment
			element = unite(element_,false);

			coeffs.resize(liste[0] + 1);
			for (int i(0); i < liste[0]; ++i) {
				coeffs[i] =new polynome_n<T>(n - 1, element, noms + 1); //polynomes nuls
			}
			coeffs[liste[0]] = new polynome_n<T>(liste + 1, n - 1, noms + 1, element_); //récurrence
		}
		else {
			element = element_;
			nul = (bool) element;
			noms_variables = NULL;
			coeffs.resize(0);
		}
		return;
	};

	~polynome_n() {
		for (int i(0); i < coeffs.size(); ++i) {
			if (coeffs[i])
				delete coeffs[i];
		}
	};

	void getString(std::string puissance,std::string& resultat) const {
		if (!nul) //n'est possible, que lors du premier appel. Car sinon on teste avant.
			resultat = "0";
		else {
			if (n_var > 0) { //récurrence
				for (int i(0); i < coeffs.size(); ++i)
					if (coeffs[i]->nul) {
						std::string puissance_locale = puissance + "*" + noms_variables[0] + "^" + std::to_string(i) + " ";
						coeffs[i]->getString(puissance_locale, resultat);
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

	friend std::ostream& operator<<(std::ostream& os, polynome_n<T> const& element) {
		std::string result = "";
		element.getString("", result);
		os << result;
		return os;
	};

	polynome_n<T>& operator=(bool const& test) {
		if (n_var > 0) {
			for (int i(1); i < coeffs.size(); ++i)
				delete coeffs[i];
			coeffs.resize(1);
			*coeffs[0] = test;
			nul = test;
		}
		else {
			if (test)
				element = unite(element,true);
			else
				element = unite(element,true);
			nul = test;
		}

		return *this;
	};

	polynome_n<T>& operator=(polynome_n<T> const& temp) {
		if (&temp == this)
			return *this;

		if (n_var > 0) {
			for (int i(0); i < coeffs.size(); ++i)
				if (coeffs[i])
					delete coeffs[i];
		}
		n_var = temp.n_var;
		nul = temp.nul;
		element = temp.element;
		noms_variables = temp.noms_variables;
		if (n_var > 0) {
			coeffs.resize(temp.coeffs.size());
			for (int i(0); i < coeffs.size(); ++i)
				coeffs[i] = new polynome_n<T>(*temp.coeffs[i]); //utilise le constructeur de copie ... 
		}
		else {
			coeffs.resize(0);
		}

		return *this;
	};

	friend polynome_n<T> operator*(polynome_n<T> const& temp1, polynome_n<T> const& temp2) {
		if (temp1.n_var != temp2.n_var)
			throw std::domain_error("polynomes à n variables de tailles différentes");
		T faux_ = unite(temp1.element,false);

		polynome_n<T> faux_p(temp1.n_var, faux_, temp1.noms_variables); //polynome vide

		if (temp1.n_var > 0) {
			bool trouve = false;
			int degre = temp1.coeffs.size() - 1 + temp2.coeffs.size() - 1;
			int trouve_int = 0;
			std::vector<polynome_n<T>*> vec_solution(degre + 1, NULL);

			if ((!temp1.nul) || (!temp2.nul)) {
				return faux_p; 
			}

			for (int i(degre); i >= 0; --i) {
				polynome_n<T> result = * faux_p.coeffs[0];

				for (int j(0); j <= i; ++j) {
					if ((j < temp1.coeffs.size()) && ((i - j) < temp2.coeffs.size()))
						if ((temp1.coeffs[j]->nul) && (temp2.coeffs[i-j]->nul)) //non-nuls
							result = result + ((*temp1.coeffs[j]) * (*temp2.coeffs[i - j]));
				}

				if ((bool)result) {
					if (!trouve) {
						trouve = true;
						trouve_int = i;
					}
				}
				if (trouve)
					vec_solution[i] = new polynome_n<T>(result);
			}

			if (!trouve) {
				vec_solution.clear();
				vec_solution.push_back(new polynome_n<T>(temp1.n_var - 1, faux_, temp1.noms_variables));
				trouve_int = 0;
			}

			if (trouve_int + 1 < vec_solution.size())
				vec_solution.erase(vec_solution.begin() + trouve_int + 1, vec_solution.end()); //faire attention

			return polynome_n<T>(temp1.n_var, temp1.noms_variables, temp1.element, trouve, vec_solution);
		}

		T result = temp1.element * temp2.element;
		return polynome_n<T>(result);
	}

	friend polynome_n<T> operator+(polynome_n<T> const& temp1, polynome_n<T> const& temp2) {
		if (temp1.n_var != temp2.n_var)
			throw std::domain_error("polynomes à n variables de tailles différentes");

		if (temp1.n_var > 0) {
			std::vector<polynome_n<T>*> result(max(temp1.coeffs.size(), temp2.coeffs.size()), NULL);
			for (int i(0); i < min(temp1.coeffs.size(), temp2.coeffs.size()); ++i) {
				result[i] = new polynome_n<T>((*temp1.coeffs[i]) + (*temp2.coeffs[i]));
			}
			if (temp1.coeffs.size() > temp2.coeffs.size()) {
				for (int i(temp2.coeffs.size()); i < temp1.coeffs.size(); ++i)
					result[i] = new polynome_n<T>(*temp1.coeffs[i]);
			}
			if (temp2.coeffs.size() > temp1.coeffs.size()) {
				for (int i(temp1.coeffs.size()); i < temp2.coeffs.size(); ++i)
					result[i] = new polynome_n<T>(*temp2.coeffs[i]);
			}

			int trouve_int = 0;
			bool trouve = false;
			for (int i(result.size() - 1); i >= 0; --i) {
				if ((bool) *result[i]) {
					trouve_int = i;
					trouve = true;
					break;
				}
			}

			for (int i(trouve_int + 1); i < result.size(); ++i)
				delete result[i];

			if (trouve_int + 1 < result.size())
				result.erase(result.begin() + trouve_int + 1, result.end());

			return polynome_n<T>(temp1.n_var, temp1.noms_variables, temp1.element, trouve, result);
		}

		T element = temp1.element + temp2.element;
		return polynome_n<T>(element);
	};


	friend polynome_n<T> operator-(polynome_n<T> const& temp1, polynome_n<T> const& temp2) {
		if (temp1.n_var != temp2.n_var)
			throw std::domain_error("polynomes à n variables de tailles différentes");

		if (temp1.n_var > 0) {
			std::vector<polynome_n<T>*> result(max(temp1.coeffs.size(), temp2.coeffs.size()), NULL);
			for (int i(0); i < min(temp1.coeffs.size(), temp2.coeffs.size()); ++i) {
				result[i] = new polynome_n<T>((*temp1.coeffs[i]) - (*temp2.coeffs[i]));
			}
			if (temp1.coeffs.size() > temp2.coeffs.size()) {
				for (int i(temp2.coeffs.size()); i < temp1.coeffs.size(); ++i)
					result[i] = new polynome_n<T>(*temp1.coeffs[i]);
			}
			if (temp2.coeffs.size() > temp1.coeffs.size()) {
				for (int i(temp1.coeffs.size()); i < temp2.coeffs.size(); ++i)
					result[i] = new polynome_n<T>( - *temp2.coeffs[i]);
			}

			int trouve_int = 0;
			bool trouve = false;
			for (int i(result.size() - 1); i >= 0; --i) {
				if ((bool) *result[i]) {
					trouve_int = i;
					trouve = true;
					break;
				}
			}

			for (int i(trouve_int + 1); i < result.size(); ++i)
				delete result[i];

			if (trouve_int + 1 < result.size())
				result.erase(result.begin() + trouve_int + 1, result.end());

			return polynome_n<T>(temp1.n_var, temp1.noms_variables, temp1.element, trouve, result);
		}

		T element = temp1.element - temp2.element;
		return polynome_n<T>(element);
	};

	friend polynome_n<T> operator-(polynome_n<T> const& temp) { //nul ne change pas. Donc on garde tout tel quel, le seul changement est : element -> - element. Et recopiage.
		if (temp.n_var > 0) {
			std::vector<polynome_n<T>*> result(0);
			result.resize(temp.coeffs.size());
			for (int i(0); i < temp.coeffs.size(); ++i)
				result[i] = new polynome_n<T>(-(*temp.coeffs[i]));
			return polynome_n<T>(temp.n_var, temp.noms_variables, temp.element, temp.nul, result);
		}
		else {
			return polynome_n<T>(-temp.element);
		}
	};


	void changer_nom(std::string* noms) {
		noms_variables = noms;
		if (n_var > 1)
			for (int i(0); i < coeffs.size(); ++i)
				coeffs[i]->changer_nom(noms + 1);
	};

	friend bool operator==(polynome_n<T> const& temp1, polynome_n<T> const& temp2) {
		if (temp1.n_var != temp2.n_var)
			return false;
		if (temp1.n_var > 0) {
			if (temp1.coeffs.size() != temp2.coeffs.size()) //suppose simplifié. Normalement OK. Tout est fait en interne (sauf monome, géré aussi).
				return false;
			for (int i(0); i < temp1.coeffs.size(); ++i)
				if (!(*temp1.coeffs[i] == *temp2.coeffs[i]))
					return false;
			return true;
		}
		else
			return (temp1.element == temp2.element);
	};

	template<class U> friend polynome_n<T> operator*(U const& scalaire, polynome_n<T> const& temp) {
		if (!(bool)temp) {
			return temp;
		}
		if ((bool)scalaire) {
			if (temp.n_var > 0) {
				std::vector<polynome_n<T>*> result(0);
				result.resize(temp.coeffs.size());
				for (int i(0); i < temp.coeffs.size(); ++i)
					result[i] = new polynome_n<T>(scalaire * (*temp.coeffs[i]));

				//au cas où le scalaire annule ...
				int trouve_int = 0;
				bool trouve = false;
				for (int i(result.size() - 1); i >= 0; --i) {
					if ((bool) *result[i]) {
						trouve_int = i;
						trouve = true;
						break;
					}
				}

				for (int i(trouve_int + 1); i < result.size(); ++i)
					delete result[i];

				if (trouve_int + 1 < result.size())
					result.erase(result.begin() + trouve_int + 1, result.end());

				return polynome_n<T>(temp.n_var, temp.noms_variables, temp.element, temp.nul, result);
			}
			else {
				return polynome_n<T>(scalaire * temp.element);
			}

		} 
		else { //retourne le polynome nul.
			T faux_ = unite( temp.element,false);
			return polynome_n<T>(temp.n_var, faux_, temp.noms_variables);
		}
	};

	operator bool() const {
		return nul;
	};

	template<class U> operator polynome_n<U>() {
		if (n_var == 0) {
			return polynome_n<U>((U)element);
		}
		std::vector<polynome_n<U>*> result(coeffs.size(), NULL);
		for (int i(0); i < coeffs.size(); ++i)
			result[i] = new polynome_n<U>((polynome_n<U>) * coeffs[i]);

		//idem, on vérifie ... au pire pas long.
		bool trouve = false;
		int trouve_int = 0;
		for (int j(coeffs.size()-1); j >= 0; --j)
			if ((bool) *result[j]) {
				trouve = true;
				trouve_int = j;
				break;
			}

		for (int i(trouve_int + 1); i < result.size(); ++i)
			delete result[i];

		if (trouve_int + 1 < result.size())
			result.erase(result.begin() + trouve_int + 1, result.end());

		return polynome_n<U>(n_var, noms_variables, (U)element, trouve, result);

	};

	std::vector<int> getDegre_tab() const { //ATTENTION : degrés et non dimensions !
		std::vector<int> tab(n_var, -1);
		getDegre(tab.data());
		return tab;
	};

	void getDegre(int* degres) const {
		int i = 0;
		for (i = coeffs.size() - 1; i >= 0; --i)
			if (coeffs[i]->nul)
				break;

		if (*degres < i)
			*degres = i;
		if (n_var == 1)
			return;

		for (int j(i); j >= 0; --j)
			if (coeffs[j]->nul)
				coeffs[j]->getDegre(degres + 1);
		return;
	};

	operator polynome_n_iter<T>() const {
		if (n_var == 0)
			return polynome_n_iter<T>(0, element, NULL);

		std::vector<int> degres = getDegre_tab();
		T faux_ = unite(element,false);

		polynome_n_iter<T> result(degres, faux_, noms_variables); //degrés ... se transforme en dimensions (+1).
		if (result.coeffs.data.size() == 1) {
			if (nul)
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
		for (int i(0); i < coeffs.size(); ++i)
			coeffs[i]->simplifier;

		int trouve_int = 0;
		bool trouve = false;
		for (int i(coeffs.size() - 1); i >= 0; --i) {
			if ((bool) *coeffs[i]) {
				trouve_int = i;
				trouve = true;
				break;
			}
		}

		for (int i(trouve_int + 1); i < coeffs.size(); ++i)
			delete coeffs[i];

		if (trouve_int + 1 < coeffs.size())
			coeffs.erase(coeffs.begin() + trouve_int + 1, coeffs.end());

		nul = trouve;
	};


	T element;
	int n_var;
	bool nul; //true si non-nul

	std::string* noms_variables;
	std::vector<polynome_n<T> *> coeffs;
};

//vérifier nul=true/false

template<class T> void parcourir_convert(polynome_n<T> const* objet, T* data, int* puissances) {
	if (objet->n_var == 0) {
		data[0] = objet->element;
		return;
	}
	for (int i(objet->coeffs.size() - 1); i >= 0; --i) {
		if (objet->coeffs[i]->nul)
			parcourir_convert(objet->coeffs[i], data + (i * (puissances[0])), puissances + 1);
	}
	return;
}
