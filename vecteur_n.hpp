#pragma once
#include <vector>
#include "fonctions template.hpp"
#include "n_for.hpp"
#include <exception>

template<class T> class vecteur_n {
public:
	vecteur_n() {};

	vecteur_n(std::vector<int> const& liste_dim) {
		if (liste_dim.size() == 0)
			throw std::domain_error("vecteur_n : liste vide");

		dimensions = liste_dim;
		puissance = dimensions.size();
		for (int i(0); i < puissance; ++i)
			if (dimensions[i] < 0)
				throw std::domain_error("vecteur_n : une dimension < 0");

		puissances= std::vector<int>(puissance);
		puissances[puissance - 1] = 1;
		for (int i(puissance-2); i >= 0; --i)
			puissances[i]= puissances[i+1]* dimensions[i+1];
		int taille = puissances[0] * dimensions[0];
		data = std::vector<T>(taille);
	};

	vecteur_n(std::vector<int> const& liste_dim,T const& element) {
		if (liste_dim.size() == 0)
			throw std::domain_error("vecteur_n : liste vide");

		dimensions = liste_dim;
		puissance = dimensions.size();
		for (int i(0); i < puissance; ++i)
			if (dimensions[i] < 0)
				throw std::domain_error("vecteur_n : une dimension < 0");

		puissances = std::vector<int>(puissance);
		puissances[puissance - 1] = 1;
		for (int i(puissance - 2); i >= 0; --i)
			puissances[i] = puissances[i + 1] * dimensions[i + 1];
		int taille = puissances[0] * dimensions[0];
		data = std::vector<T>(taille,element);
	};

	vecteur_n(vecteur_n<T> const& temp) {
		puissance = temp.puissance;
		dimensions = temp.dimensions;
		data = temp.data;
		puissances = temp.puissances;
	};
	
	vecteur_n<T>& operator=(vecteur_n<T> const& temp) {
		if (this == &temp)
			return *this;
		puissance = temp.puissance;
		dimensions = temp.dimensions;
		data = temp.data;
		puissances = temp.puissances;
		return *this;
	};

	/*
	vecteur_n<T>& operator=(bool test) {
		for (int i(0); i < data.size(); ++i)
			data[i] = test;
		return *this;
	};
	*/

	int accesseur(std::vector<int> const& positions) const { //position doit être de la bonne taille (et du bon format)
		int expo = 1;
		int pos = 0;
		for (int i(puissance - 1); i >= 0; --i) {
			pos = pos + positions[i] * puissances[i];
		}
		return pos;
	};

	std::vector<int> position(int pos) const{
		std::vector<int> retour(puissance, 0);

		for (int i(puissance - 1); i >= 0; --i) {
			retour[i] = pos % dimensions[i];
			pos = pos / dimensions[i];
		}

		return retour;
	};

	bool accesseur_possible(std::vector<int> const& position) const {
		if (position.size() != puissance)
			return false;
		for (int i(0); i < puissance; ++i)
			if (position[i] >= dimensions[i])
				return false;
		return true;
	};

	inline T& operator[](std::vector<int> position)  { //
		if (accesseur_possible(position))
			return data[accesseur(position)];
	};

	inline T operator[](std::vector<int> position) const { //
		if (accesseur_possible(position))
			return data[accesseur(position)];
		T element = unite(data[0],false);
		return element;
	};

	inline T& operator[](int position) {
		return data[position];
	};

	inline T operator[](int position) const {
		return data[position];
	};

	inline void modifier_dimension(std::vector<int> const& nouvelles_dimensions) {
		T element = unite(data[0],false);
		modifier_dimension(nouvelles_dimensions, element);
	};

	void modifier_dimension(std::vector<int> const& nouvelles_dimensions, T element) { //même puissance
		bool test=true;
		for(int i(0);i<puissance;++i)
			if (nouvelles_dimensions[i] != dimensions[i]) {
				test = false;
				break;
			}
		if (test)
			return;

		//faire la boucle de dimension2, et recopier (si possible). boucle intelligente.
		std::vector<int> pos(puissance, 0);
		int pos_int = 0;
		int pos_int2 = 0;

		std::vector<int> puissances2(puissance);
		puissances2[puissance - 1] = 1;

		for (int i(puissance - 2); i >= 0; --i) {
			puissances2[i] = puissances2[i + 1] * nouvelles_dimensions[i + 1];
		}
		int taille = puissances2[0] * nouvelles_dimensions[0];
		std::vector<T> data2(taille, element); //on met des 0 quand on peut ... en gros.


		while (true) {
			data2[pos_int2] = data[pos_int];
			++pos_int;
			++pos_int2;
			++pos[puissance - 1];

			int i = puissance - 1;
			while (pos[i] >= min(nouvelles_dimensions[i], dimensions[i])) {
				if (i == 0){
					i = -1;
					break;
				};

				pos_int = pos_int - pos[i] * puissances[i] + puissances[i - 1];
				pos_int2 = pos_int2 - pos[i] * puissances2[i] + puissances2[i - 1];
				pos[i] = 0;
				++pos[i - 1];
				--i;
			};

			if (i == -1)
				break;
		};

		data = data2;
		dimensions = nouvelles_dimensions;
		puissances = puissances2;

		return;
	};

	void simplifier() {
		if ((bool)data[data.size() - 1])
			return;
		std::vector<int> dimensions_max(puissance, 0);

		for (n_for iter(dimensions); (bool) iter; ++iter)
			if ((bool)data[iter.position])
				for (int i(0); i < puissance; ++i)
					if (iter.positions[i] > dimensions_max[i])
						dimensions_max[i] = iter.positions[i];

		bool modifie = false; //on passe à la taille .... dernier element +1
		for (int i(0); i < puissance; ++i) {
			++dimensions_max[i];
			if (dimensions_max[i] != dimensions[i])
				modifie = true;
		}

		if (modifie)
			modifier_dimension(dimensions_max);
	};

	std::vector<int> puissances;
	int puissance; //nombre de variables ... data[1][2]...[puissance]
	std::vector<int> dimensions;
	std::vector<T> data;
};

//verifier >= dimension.
//simplifier ...