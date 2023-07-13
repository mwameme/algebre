#pragma once

#include "polynome.hpp"
#include "rationnel.hpp"
#include "matrice.hpp"
#include <vector>
#include <set>

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

	std::vector<impedance<T>*> fils;
	std::vector<diode<T>*> diodes;
	std::vector<AO<T>*> AOs;

	std::vector< rationnel<polynome<T>>> gains;
	polynome<T> transitoire; //racines

	// a la fin seulement
	std::vector<double> racines;
	std::vector<double> coeffs;
	std::vector<double> derivees;
};

template<class T> class impedance {
public:
	rationnel<polynome<T>> valeur;
	noeud<T>* a;
	noeud<T>* b;

	int numero;
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