#pragma once
#include <vector>
#include <iostream>

class fact_for {
public:
//	fact_for() {};

	fact_for(int n_, bool signature_=true) : positions(n_,0) , permutation(n_,0) { //: dimensions(dim), positions(dimensions.size(),0),position(0),OK(true)
		if (signature_)
			signature = 1;
		else
			signature = 0;

		n = n_;
		if (n > 0)
			OK = true;
		else
			OK = false;
		
		dimensions.resize(n);
		for (int i(0); i < n; ++i)
			dimensions[i] = n - i;

		places_dispo.resize(n);
		for (int i(0); i < n; ++i) {
			for (int j(0); j < dimensions[i]; ++j) {
				places_dispo[i].push_back(i + j);
			}
		}

		if (signature) {
			signatures.resize(n);
			for (int i(0); i < n; ++i)
				signatures[i] = 1;
		}
		else
			signatures.resize(0);

		for (int i(0); i < n; ++i)
			permutation[i] = places_dispo[i][0];

		position = 0;
	};

	fact_for(fact_for const& temp) {
		signature=temp.signature; //calculer la signature !
		OK=temp.OK; //tant qu'on n'a pas tout parcouru
		n=temp.n; //pour fact(n)
		positions=temp.positions; // positions[i] < dimensions[i]
		dimensions=temp.dimensions; //vecteur n, n-1 , ...
		places_dispo=temp.places_dispo; //vecteur des choix possibles
		signatures=temp.signatures; //vecteur des signatures partielles.
		permutation = temp.permutation;
		position = temp.position;
	};

	fact_for& operator=(fact_for const& temp) {
		signature = temp.signature; //calculer la signature !
		OK = temp.OK; //tant qu'on n'a pas tout parcouru
		n = temp.n; //pour fact(n)
		positions = temp.positions; // positions[i] < dimensions[i]
		dimensions = temp.dimensions; //vecteur n, n-1 , ...
		places_dispo = temp.places_dispo; //vecteur des choix possibles
		signatures = temp.signatures; //vecteur des signatures partielles.
		permutation = temp.permutation;
		position = temp.position;

		return *this;
	};

	fact_for& operator++() {
		int i = positions.size() - 1;
		++position;
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

		for (int j(i + 1); j < n; ++j) {
			places_dispo[j] = places_dispo[j - 1];
			places_dispo[j].erase(places_dispo[j].begin() + positions[j - 1]);
		}

		for (int j(i); j < n; ++j) {
			permutation[j] = places_dispo[j][positions[j]];
		}

		if (signature) {
			for (int j(i); j < n; ++j) {
				if (j == 0)
					continue;
				signatures[j] = signatures[j - 1];
				for (int k(0); k < j; ++k)
					if (permutation[k] > permutation[j])
						signatures[j] = -signatures[j];
			}
			signature = signatures[n - 1];
		}

		return *this;
	};

	friend std::ostream& operator<<(std::ostream& os, fact_for const& element) {
		os << "(";
		for (int i(0); i < element.n; ++i) {
			os << element.permutation[i];
			if (i < element.n - 1)
				os << ",";
		}
		os << "), " << element.signature;
		return os;
	};

	explicit inline operator bool() const {// si OK=false, la boucle s'arrête.
		return OK;
	};

	int position;
	int signature; //calculer la signature !
	bool OK; //tant qu'on n'a pas tout parcouru
	int n; //pour fact(n)
	std::vector<int> positions; // positions[i] < dimensions[i]
	std::vector<int> dimensions; //vecteur n, n-1 , ...
	std::vector<std::vector<int>> places_dispo; //vecteur des choix possibles
	std::vector<int> signatures; //vecteur des signatures partielles.
	std::vector<int> permutation; //résultat final
};