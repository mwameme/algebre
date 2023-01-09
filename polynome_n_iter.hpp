#pragma once
#include "vecteur_n.hpp"
#include <vector>
#include <string>
#include "polynome_n.hpp"
#include <exception>

#include "entete objets.hpp"


template<class T> class polynome_n;

template<class T> polynome_n<T>* poly_n_convert_rec(T* data, int* dimensions, int* puissances, int n, std::string* noms);

template<class T> class polynome_n_iter {
public:

	polynome_n_iter() {};

	polynome_n_iter(int n_var, T element,std::string* noms_) {//polynome de degre 0. utilisé pour vrai/faux
		if (n_var <= 0)
			throw std::domain_error("polynome_n_iter : n_var <= 0.");

		std::vector<int> dimensions(n_var, 1);
		coeffs(dimensions);
		coeffs.data[0] = element;
		noms = noms_;
	};

	polynome_n_iter(std::vector<int> degres, T element, std::string* noms_){ //monome. degres >0. Transforme degres en dimensions ...  ATTENTION créé avec degrés
		if ( degres.size() ==0)
			throw std::domain_error("polynome_n_iter : liste vide.");

		for(int i(0);i<degres.size();++i)
			if (degres[i] < 0) {
				degres = std::vector<int>(degres.size(), 0);
				element = unite(element,false);
				break;
			}

		for (int i(0); i < degres.size(); ++i)
			degres[i] += 1; //dimensions
		coeffs = vecteur_n<T>(degres);
		T faux_ = unite(element,false);

		for (int i(0); i < coeffs.data.size(); ++i)
			coeffs.data[i] = faux_;
		coeffs.data[coeffs.data.size() - 1] = element; //monome
		noms = noms_;
	};

	polynome_n_iter(vecteur_n<T> tableau, std::string* noms_) { //au cas où ?
		coeffs = tableau;
		noms = noms_;
	};

	polynome_n_iter(polynome_n_iter const& temp) : coeffs(temp.coeffs) , noms(temp.noms){

	};

	polynome_n_iter<T>& operator=(polynome_n_iter<T> temp) {
		coeffs = temp.coeffs;
		noms = temp.noms;
		return *this;
	};

	polynome_n_iter<T>& operator=(bool test) {
		T temp;
		if (test)
			temp = unite(coeffs.data[0],true);
		else
			temp = unite(coeffs.data[0],false);

		std::vector<int> dimensions(coeffs.puissance, 1);
		coeffs = vecteur_n<T>(dimensions);
		coeffs.data[0] = temp;
		return *this;
	};

	friend polynome_n_iter<T> operator+(polynome_n_iter<T> const& gauche_, polynome_n_iter<T> const& droite_) {
		if (gauche_.coeffs.puissance != droite_.coeffs.puissance)
			throw std::domain_error("polynome_n_iter, addition : le nombre de variables ne correspond pas.");

		polynome_n_iter<T> gauche = gauche_;
		polynome_n_iter<T> droite = droite_;
		int n = gauche.coeffs.puissance;
		std::vector<int> dim(n, 0);
		for (int i(0); i < n; ++i)
			dim[i] = max(gauche.coeffs.dimensions[i], droite.coeffs.dimensions[i]);
		gauche.coeffs.modifier_dimension(dim);
		droite.coeffs.modifier_dimension(dim);
		for (int i(0); i < gauche.coeffs.data.size(); ++i)
			gauche.coeffs.data[i] = gauche.coeffs.data[i] + droite.coeffs.data[i];
		return gauche;
	};

	friend polynome_n_iter<T> operator-(polynome_n_iter<T> const& gauche_, polynome_n_iter<T> const& droite_) {
		if (gauche_.coeffs.puissance != droite_.coeffs.puissance)
			throw std::domain_error("polynome_n_iter, soustraction : le nombre de variables ne correspond pas.");

		polynome_n_iter<T> gauche = gauche_;
		polynome_n_iter<T> droite = droite_;
		int n = gauche.coeffs.puissance;
		std::vector<int> dim(n, 0);
		for (int i(0); i < n; ++i)
			dim[i] = max(gauche.coeffs.dimensions[i], droite.coeffs.dimensions[i]);
		gauche.coeffs.modifier_dimension(dim);
		droite.coeffs.modifier_dimension(dim);
		for (int i(0); i < gauche.coeffs.data.size(); ++i)
			gauche.coeffs.data[i] = gauche.coeffs.data[i] - droite.coeffs.data[i];
		return gauche;
	};

	friend polynome_n_iter<T> operator-(polynome_n_iter<T> const& temp) {
		polynome_n_iter<T> result = temp;
		for (int i(0); i < result.coeffs.data.size(); ++i)
			result.coeffs.data[i] = - result.coeffs.data[i];
		return result;
	};

	friend polynome_n_iter<T> operator*(polynome_n_iter<T> const& gauche, polynome_n_iter<T> const& droite) {
		if (gauche.coeffs.puissance != droite.coeffs.puissance)
			throw std::domain_error("polynome_n_iter, multiplication : le nombre de variables ne correspond pas.");

		int n = gauche.coeffs.puissance;
		std::vector<int> degres(n);
		for (int i(0); i < n;++i) 
			degres[i] = gauche.coeffs.dimensions[i] + droite.coeffs.dimensions[i] - 2;

		T faux_ = unite(gauche.coeffs.data[0],false);

		polynome_n_iter<T> result(degres , faux_ , gauche.noms); //ATTENTION degrés et non dimensions.

		std::vector<int> dimensions_gauche = gauche.coeffs.dimensions;
		std::vector<int> dimensions_droite = droite.coeffs.dimensions;

		std::vector<int> positions_gauche(n, 0);
		std::vector<int> positions_droite(n, 0);
		int position_gauche = 0;
		int position_droite = 0;
		int position_finale = 0;

		std::vector<int> puissances_finale = result.coeffs.puissances;
		std::vector<int> puissances_gauche = gauche.coeffs.puissances;
		std::vector<int> puissances_droite = droite.coeffs.puissances;


		bool fin = true;
		while (fin) {
			//on opère la multiplication ... chaque couple (gauche , droite) apparait une et une seule fois.
			if (((bool)gauche.coeffs.data[position_gauche]) && ((bool)droite.coeffs.data[position_droite])) 
				result.coeffs.data[position_finale] = result.coeffs.data[position_finale] + (gauche.coeffs.data[position_gauche] * droite.coeffs.data[position_droite]);
			++position_gauche;
			++positions_gauche[n - 1];
			++position_finale;
			int i = n - 1;
			bool fin_gauche = true;

			while (positions_gauche[i] >= dimensions_gauche[i]) { //On réajuste gauche, et donc finale aussi.
				position_finale = position_finale - puissances_finale[i] * positions_gauche[i];
//				position_gauche = position_gauche - puissances_gauche[i] * positions_gauche[i];

				positions_gauche[i] = 0;
				if (i == 0) {
					fin_gauche = false;
					break;
				}
				++positions_gauche[i - 1];
//				position_gauche = position_gauche + puissances_gauche[i - 1];
				position_finale = position_finale + puissances_finale[i - 1];
				--i;
			}

			fin = true;
			if (!fin_gauche) { //gauche a dépassé, et vaut  0 partout ... on incrémente à droite.
				++position_droite;
				++positions_droite[n - 1];
				++position_finale;
				int i = n - 1;

				while (positions_droite[i] >= dimensions_droite[i]) {
//					position_finale = position_finale - puissances_finale[i] * positions_droite[i];
					position_droite = position_droite - puissances_droite[i] * positions_droite[i];

					positions_droite[i] = 0;
					if (i == 0) {
						fin = false;
						break;
					}
					++positions_droite[i - 1];
//					position_droite = position_droite + puissances_droite[i - 1];
					position_finale = position_finale + puissances_finale[i - 1];
					--i;
				}
			}
			//fin du while
		}//Le calcule marche, car le vecteur finale est la sommme du vecteur gauche et droite, et cela se calcule avec les puissances (finales). 

		return result;

	};

	explicit operator bool() const {
		for (int i(0); i < coeffs.data.size(); ++i)
			if ((bool)coeffs.data[i])
				return true;
		return false;
	};

	operator polynome_n<T>() const {
		return *poly_n_convert_rec(coeffs.data.data(), coeffs.dimensions.data(), coeffs.puissances.data(), coeffs.puissance, noms);
	};

	friend std::ostream& operator<<(std::ostream& os, polynome_n_iter<T> const& element) {
		bool premier = false;
		for (n_for iter(element.coeffs.dimensions); (bool)iter; ++iter)
			if ((bool)element.coeffs.data[iter.position]) {
				std::string puissance="";
				for (int i(0); i < coeffs.puissance; ++i)
					puissance = puissance + "*" + element.noms[i] + "^" + std::to_string(iter.positions[i]) + " ";
				if (premier) 
					os << "+" << coeffs.data[iter.position] << " " << puissance;
				else{
					premier = true;
					os << element.coeffs.data[iter.position] << " " << puissance;
				}

			}
		return os;
	}
	/*
	*/


	std::string* noms;
	vecteur_n<T> coeffs;
};

template<class T> polynome_n<T>* poly_n_convert_rec(T* data, int* dimensions, int* puissances, int n, std::string* noms) {
	if (n == 0)
		return new polynome_n<T>(*data);

	std::vector<polynome_n<T>*> tab(0);
	for (int i(0); i < *dimensions; ++i)
		tab.push_back(poly_n_convert_rec(data + (i * puissances[0]), dimensions + 1, puissances + 1, n - 1, noms + 1));

	int i = 0;
	for (i(*dimensions - 1); i >= 1; ++i)
		if (tab[i]->nul)
			break;
		else
			delete tab[i];

	if (i < *dimensions - 1)
		tab.erase(tab.begin() + i + 1, tab.end());

	bool nul = true;
	if ((i == 0) && !(tab[i]->nul))
		nul = false;

	T faux_ = unite(*data,false);
	return new polynome_n<T>(n, noms, faux_, nul, tab);
};




// Verifier les dimensions : polynome nul. OK
// convertir en polynome_n, et réciproquement. OK
// mettre a jour types et norme.
// simplifier le vecteur_n OK
// vérifier puissances OK
// const OK
// ostream