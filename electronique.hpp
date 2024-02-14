#pragma once

#include "polynome.hpp"
#include "rationnel.hpp"
#include "matrice.hpp"
#include "unite.hpp"

#include <vector>
#include <set>
#include <exception>

template<class T> class noeud;
template<class T> class AO;
template<class T> class diode;
template<class T> class impedance;

template<class T> class schema_global;
template<class T> class schema_diode;
template<class T> class schema_separe;
template<class T> class schema_AO;


template<class T> class schema_global {
public:
	int nbr_noeud;

	std::vector<noeud<T>> noeuds;
	std::vector<diode<T>> diodes;
	std::vector<AO<T>> AOs;
	std::vector<impedance<T>> impedances;

	std::vector<int> entrees; //numero des entrees
	noeud<T>* masse;
	std::set<schema_diode<T>>;

	void ajouter_impedance(rationnel<polynome<T>> valeur, int i, int j) {
		if (i == j)
			throw std::domain_error("ajouter impedance : i==j");
		if (max(i, j) >= noeuds.size()) {
			for (int k(noeuds.size()); k <= max(i, j); ++k)
				noeuds.push_back(noeud<T>(k));
			impedances.push_back(impedance<T>(&noeuds[i], &noeuds[j], valeur));
			impedances[impedances.size() - 1].numero = impedances.size() - 1;
			noeuds[i].impedances.push_back(&impedances[impedances.size() - 1]);
			noeuds[j].impedances.push_back(&impedances[impedances.size() - 1]);
			return;
		}

		noeud<T>* local_a = &noeuds[i];
		noeud<T>* local_b = &noeuds[j];

		bool trouve = false;
		for (int k(0); k < local_a->impedances.size(); ++k) {
			if ((local_a->impedances[k]->a == local_b) || (local_a->impedances[k]->b == local_b)){
				rationnel<polynome<T>> vrai = unite(local->impedances[k]->valeur, true);
				local_a->impedances[k]->valeur = vrai / ((vrai / local->impedances[k]->valeur) + (vrai / valeur));
				trouve = true;
				break;
			}
		}
		if (!trouve) {
			impedances.push_back(impedance<T>(local_a, local_b, valeur));
			impedances[impedances.size() - 1].numero = impedances.size() - 1;
			local_a->impedances.push_back(&impedances[impedances.size() - 1]);
			local_b->impedances.push_back(&impedances[impedances.size() - 1]);
		}
		return;
	}

	void ajouter_resistance(T valeur, int i, int j) {
		polynome<T> poly(valeur);
		rationnel<polynome<T>> temp( poly, unite(poly, true));
		ajouter_impedance(temp, i, j);
	}

	void ajouter_condensateur(T valeur, int i, int j) {
		polynome<T> poly(unite(valeur,false),valeur);
		rationnel<polynome<T>> temp(unite(poly, true), poly);
		ajouter_impedance(temp, i, j);
	}

	void ajouter_bobine(T valeur, int i, int j) {
		polynome<T> poly(unite(valeur, false), valeur);
		rationnel<polynome<T>> temp( poly , unite(poly, true));
		ajouter_impedance(temp, i, j);
	}

	void ajouter_diode(int i, int j) { // parfaite
		
	}


};

template<class T> class schema_diode {
public:
	std::vector<bool> diodes_ouvertes;
	int nombre_sous_schemas;

	std::set<schema_separe<T>> liste_schemas;
	std::vector< std::set<schema_separe<T>>::iterator> num_schema; // donne dans quel set est le noeud (ancien) que l'on recherche.
	//	std::vector<int> num_schema; //pour retrouver Va ...
	std::vector<int> place_dans_schema; //dans le set du schema, quel est le numero du noeud
};

template<class T> class schema_separe {
public:
	//on recopie tout
	int nbr_noeud;
	std::vector<noeud<T>> noeuds;
	std::vector<diode<T>> diodes;
	std::vector<AO<T>> AOs;
	std::vector<impedance<T>> impedances;

	std::vector<int> entrees; //numero des entrees
	noeud<T>* masse;
	std::set<schema_AO<T>>;
};

template<class T> class schema_AO {
public:
	std::vector<int> positions_AO;

	int nbr_noeud;
	std::vector<noeud<T>> noeuds;
	std::vector<diode<T>> diodes;
	std::vector<AO<T>> AOs; // basculé -1 0 1
	std::vector<impedance<T>> impedances;

	std::vector<int> entrees; //numero des entrees
	noeud<T>* masse;
};

template<class T> class noeud {
public:
	bool masse;
	bool fuite;
	int numero;

	std::vector<impedance<T>*> impedances;
	std::vector<diode<T>*> diodes;
	std::vector<AO<T>*> AOs;

	std::vector< rationnel<polynome<T>>> gains;
	polynome<T> transitoire; //racines

	// a la fin seulement
	std::vector<double> racines;
	std::vector<double> coeffs;
	std::vector<double> derivees;

	noeud(int num) {
		masse = false;
		fuite = false;
		numero = num;
	}

};

template<class T> class impedance {
public:
	rationnel<polynome<T>> valeur;
	noeud<T>* a;
	noeud<T>* b;

	int numero;

	impedance(noeud<T>* a_, noeud<T>* b_, rationnel<polynome<T>> val) {
		valeur = val;
		a = a_;
		b = b_;
		return;
	}
};

template<class T> class diode {
public:
	bool parfaite; // si true, est parfaite
	bool ouverte; //true si connectee
	noeud<T>* a;
	noeud<T>* b;

	rationnel<polynome<T>> resitance; // seulement si non-parfaite
	int numero; // pour la tension de seuil. Idem si non-parfaite

	//a la toute fin
	double Vsat;
};

template<class T> class AO {
public:
	bool parfait; //si true est parfait. Sinon, idéal.
	int position; //0 non-saturé, -1 = -Vsat, +1 = +Vsat

	noeud<T>* a;
	noeud<T>* b;
	noeud<T>* s;//sortie. Fuite.

	rationnel<polynome<T>> gain; //dans le cas non-parfait
	int numero; // si saturé

	double Vsat; //a la toute fin ...

	rationnel<polynome<T>> critere_stabilité; //(partie réelle doit être négative)
};