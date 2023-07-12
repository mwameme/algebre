#pragma once

#include "polynome.hpp"
#include "rationnel.hpp"
#include "matrice.hpp"
#include <vector>
#include <set>

template<class T> class schema_global {
public:
	using gain = typename rationnel<polynome<T>>;

	int nbr_noeud;
	std::vector<noeud> noeuds;
	std::vector<diode> diodes;
	std::vector<AO> AOs;

	std::vector<int> entrees;

};