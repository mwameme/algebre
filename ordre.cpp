#include "ordre.hpp"

std::vector<int> ordre_double(std::vector<double> X) { //décroissant
	return ordre<double>(X, [](double gauche, double droite) {return gauche >= droite; });
};

std::vector<int> ordre_alea(int n) {
	std::vector<double> liste(n, 0);
	for (int i(0); i < n; ++i)
		liste[i] = ((double)(rand()) / ((double)(RAND_MAX)));
	return ordre_double(liste);
};
