
#include "rationnel.hpp"
#include "polynome.hpp"
#include <complex>
#include "unite.hpp"
#include <vector>

class reseau {
public:
	noeud* entree;
	noeud* sortie;

	std::vector<noeud> noeuds;

};

class noeud {
public:
	std::vector<fil> liste_fil;
	std::vector<noeud> liste_noeud;
};

class impedance {
public:
	impedance(rationnel<polynome<std::complex<double>>> const& temp) : valeur(temp) {};
	rationnel<polynome<std::complex<double>>> valeur;

	impedance inverse() const;
	impedance ajouter(impedance const& autre) const;
	static impedance ajouter(std::vector<impedance> const& vec);
};

impedance impedance::inverse() const {
	rationnel<polynome<std::complex<double>>> vrai = unite(valeur, true);
	return impedance(vrai / valeur);
}

impedance impedance::ajouter(impedance const& autre) const {
	return impedance(valeur + autre.valeur);
}

impedance impedance::ajouter(std::vector<impedance> const& vec) {
	auto result = unite(vec[0].valeur,false);
	for (int i(0); i < vec.size(); ++i)
		result += vec[i].valeur;
	return impedance(result);
}

class fil {
public:
	noeud* gauche;
	noeud* droit;
	impedance resistance;
};

class AO {
	noeud* 
};