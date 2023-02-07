#pragma once

#include <algorithm>
#include <vector>
#include <functional>
#include "TVD.hpp"

int compter_ordonne(std::vector<std::vector<int>> const& votes, std::vector<int> ordre);

std::vector<double> vote(std::vector<std::vector<int>> const& votes, double epsilon); //une seule fois 

int get_min(std::vector<std::vector<int>> votes, double epsilon, int d, bool asymetrique); //enleve le max, it�r�, jusqu'� trouver le min. epsilon=2.

std::vector<int> get_max_min(std::vector<std::vector<int>> votes, double epsilon, int d, bool asymetrique); //enleve le min (de get_min), it�r�, jusqu'� obtenir le max.
//appeler asym�trique=false, d=0, epsilon=2. Renvoit une liste ordonn�e (du plus grand au plus petit).

std::vector<int> ordre_alea(int n);
//pour g�n�rer des votes al�atoires ...

//ordre 
template<class T> class couple_index {
public:
	couple_index() {};

	int index;
	T element;
};

template<class T> std::vector<int> ordre(std::vector<T> liste_, std::function<bool(T, T)> const& f) {
	int n = liste_.size();
	std::vector<couple_index<T>> liste(n);
	for (int i(0); i < n; ++i) {
		liste[i].element = liste_[i];
		liste[i].index = i;
	}

	std::sort(liste.begin(), liste.end(), [&](couple_index<T> gauche, couple_index<T> droite) {return f(gauche.element, droite.element); });
	std::vector<int> resultat(n);
	for (int i(0); i < n; ++i)
		resultat[i] = liste[i].index;

	return resultat;
};

inline std::vector<int> ordre_double(std::vector<double> X) { //d�croissant
	return ordre<double>(X, [](double gauche, double droite) {return gauche >= droite; });
};
