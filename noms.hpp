#pragma once

/*
#include <vector>
#include <string>
#include <exception>

#include "polynome_n_rec.hpp"
#include "polynome_n_iter.hpp"
#include "polynome_n_sparse.hpp"
#include "unite.hpp"

template<class T>  T unite(T const& element, bool test);


template<class T> class polynome_n_rec;
template<class T> class polynome_n_iter;
template<class T> class polynome_n_sparse;
template<class T> class monome;


std::vector<std::string> union_nom(std::vector<std::string > const& gauche, std::vector<std::string> const& droit);


template<class T>
polynome_n_sparse<T> modifier_nom(std::vector<std::string > ancien_nom, polynome_n_sparse<T> const& ancien_poly, std::vector<std::string > nouveau_nom) {
	std::vector<int> permutation(ancien_nom.size(),0);
	for (int i(0); i < ancien_nom.size(); ++i) {
		int j=0;
		for (; j < nouveau_nom.size(); ++j)
			if (ancien_nom[i] == nouveau_nom[j])
				break;
		if (j == nouveau_nom.size())
			throw std::domain_error("changement de nom : le nouvel ensemble ne contient pas l'ancien");
		permutation[i] = j;
	}

	polynome_n_sparse<T> nouveau_poly(nouveau_nom.size());
	nouveau_poly.monomes.reserve(ancien_poly.monomes.size());
	for (int i(0); i < ancien_poly.monomes.size(); ++i) {
		std::vector<int> nouveau_degres(nouveau_nom.size(),0);
		for (int j(0); j < permutation.size(); ++j)
			nouveau_degres[permutation[j]] = ancien_poly.monomes[i].degres[j];
		nouveau_poly.monomes.push_back(monome<T>(nouveau_degres, ancien_poly.monomes[i].element));
	}
	nouveau_poly.est_trie = false;
	nouveau_poly.simplifier();
	return nouveau_poly;
};

template<class T>
polynome_n_iter<T> modifier_nom(std::vector<std::string > ancien_nom, polynome_n_iter<T> const& ancien_poly, std::vector<std::string > nouveau_nom) {
	std::vector<int> permutation(ancien_nom.size(), 0);
	for (int i(0); i < ancien_nom.size(); ++i) {
		int j=0;
		for (; j < nouveau_nom.size(); ++j)
			if (ancien_nom[i] == nouveau_nom[j])
				break;
		if (j == nouveau_nom.size())
			throw std::domain_error("changement de nom : le nouvel ensemble ne contient pas l'ancien");
		permutation[i] = j;
	}

	std::vector<int> dimensions(nouveau_nom.size(), 1);
	for (int i(0); i < ancien_nom.size(); ++i)
		dimensions[permutation[i]] = ancien_poly.coeffs.dimensions[i];
	polynome_n_iter<T> nouveau_poly(dimensions,unite(ancien_poly.coeffs.data[0],false),false);

	//poly créé ... maintenant on remplit.
	for (int i(0); i < ancien_poly.coeffs.data.size(); ++i) {
		if (!(bool)ancien_poly.coeffs.data[i])
			continue;
		//peut être accéléré, en modifiant au fur et à mesure le vecteur de positions (finales et initiales), comme dans polynome_n_iter::operator*.
		//seulement cette fonction n'est pas souvent utilisée ... peu de perte de temps.
		std::vector<int> ancien_degres = ancien_poly.coeffs.positions(i);
		std::vector<int> nouveau_degres(nouveau_nom.size(), 0);
		for (int j(0); j < permutation.size(); ++j)
			nouveau_degres[permutation[j]] = ancien_degres[j];
		nouveau_poly.coeffs.data[nouveau_poly.coeffs.position(nouveau_degres)] = ancien_poly.coeffs.data[i];
	}
	return nouveau_poly;
};

template<class T>
polynome_n_rec<T> modifier_nom(std::vector<std::string > ancien_nom, polynome_n_rec<T> const& ancien_poly, std::vector<std::string > nouveau_nom) {
	std::vector<int> permutation(ancien_nom.size(), 0);
	for (int i(0); i < ancien_nom.size(); ++i) {
		int j=0;
		for (; j < nouveau_nom.size(); ++j)
			if (ancien_nom[i] == nouveau_nom[j])
				break;
		if (j == nouveau_nom.size())
			throw std::domain_error("changement de nom : le nouvel ensemble ne contient pas l'ancien");
		permutation[i] = j;
	}

	polynome_n_rec<T> nouveau_poly(nouveau_nom.size(), unite(ancien_poly.element, false));
	for (auto it = ancien_poly.cbegin(); (bool)it; ++it) {
		std::vector<int> ancien_degres = it.positions;
		std::vector<int> nouveau_degres(nouveau_nom.size(), 0);
		for (int j(0); j < permutation.size(); ++j)
			nouveau_degres[permutation[j]] = ancien_degres[j];
		nouveau_poly += polynome_n_rec<T>(nouveau_degres, *it);
	}
	return nouveau_poly;
}

*/