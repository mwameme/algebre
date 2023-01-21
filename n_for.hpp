#pragma once
#include <vector>
#include "fonctions template.hpp"

class n_for {
public:
	n_for() {};

	n_for(std::vector<int> dim,bool croissant_=false) {
		dimensions = dim;
		positions.resize(dimensions.size());
		croissant = croissant_;
		OK = true;
		position = 0;

		for (int i(0); i < positions.size(); ++i)
			positions[i] = 0;

		int taille = dimensions.size();
		for (int i(0); i < taille; ++i) {
			if (dimensions[i] <= 0) {
				OK = false;
				break;
			}
		}
	};

	n_for(n_for const& temp) {
		croissant = temp.croissant;
		OK = temp.OK;
		positions = temp.positions;
		position = temp.position;
		dimensions = temp.dimensions;
	};

	n_for& operator=(n_for const& temp) {
		croissant = temp.croissant;
		OK = temp.OK;
		positions = temp.positions;
		position = temp.position;
		dimensions = temp.dimensions;
		return *this;
	};

	n_for& operator++() {
		++position;
		int i = positions.size() - 1;
		++positions[i];
		while (positions[i] >= dimensions[i]) {
			if (i == 0) {
				OK = false;
				return *this;
			}
			positions[i] = 0;
			++positions[i - 1];
			--i;
		}

		if (!croissant)
			return *this;

		/*
		std::vector<int> puissances(positions.size(), 1);
		for (int j(positions.size() - 2); j >= i+1; ++j) 
			puissances[j] = dimensions[j + 1] * puissances[j + 1];
			*/

		for (int j(i + 1); j < positions.size(); ++j)
			if (positions[j] < positions[j - 1]) {
//				position = position + (positions[j - 1] - positions[j]) * puissances[j];
				positions[j] = positions[j - 1];
				if (positions[j] >= dimensions[j]) {
					OK = false;
					return *this;
				}
			}

		return *this;
	};

	int getFact() const {
		if (!croissant)
			return 1;

		std::vector<int> multi(0);
		int i = 0;
		int j;
		for (j=1; j < dimensions.size(); ++j) {
			if (positions[j] != positions[j - 1]) {
				multi.push_back(j - i);
				i = j;
			}
		}
		multi.push_back(j - i);

		int result = factorielle(positions.size());
		int produit = 1;
		for (int i(0); i < multi.size(); ++i)
			produit = produit * factorielle(multi[i]);
		return result/produit;
	};

	inline int size() const {
		return positions.size();
	};

	int getPosition(std::vector<int> const& position) const { //position doit être de la bonne taille (et du bon format)
		int expo = 1;
		int pos = 0;
		for (int i(dimensions.size() - 1); i >= 0; --i) {
			pos = pos + expo * position[i];
			expo = expo * dimensions[i];
		}
		return pos;
	};

	int getPosition() const {
		if (!croissant)
			return position;

		std::vector<int> puissances(positions.size());
		puissances[positions.size() - 1] = 1;
		for (int j(positions.size() - 2); j >= 0; --j)
			puissances[j] = dimensions[j + 1] * puissances[j + 1];

		int pos = 0;
		for (int j(0); j < positions.size(); ++j)
			pos = pos + positions[j] * puissances[j];

		return pos;
	};

	std::vector<int> getPositions(int pos) const {
		std::vector<int> retour(dimensions.size(), 0);

		for (int i(dimensions.size() - 1); i >= 0; --i) {
			retour[i] = pos % dimensions[i];
			pos = pos / dimensions[i];
		}

		return retour;
	};

	explicit inline operator bool() const {// si OK=false, la boucle s'arrête.
		return OK;
	};

	bool croissant;
	bool OK;
	int position;
	std::vector<int> positions;
	std::vector<int> dimensions;

};