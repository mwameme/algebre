#pragma once

#include <vector>
#include <algorithm>
#include <functional>

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

std::vector<int> ordre_double(std::vector<double> X); //décroissant


std::vector<int> ordre_alea(int n); // liste aléatoire ...
