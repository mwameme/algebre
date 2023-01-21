#pragma once
#include "polynome.hpp"


template<class T> class bloc_interpolation {
public : 

	bloc_interpolation(int n_max, T* liste_coeffs, int debut_gauche, int largeur_) { //debut commence à 0 ... n_max taille de la liste. Largeur une puissance de 2
		largeur = largeur_;
		while (n_max <= (debut_gauche + (largeur / 2))) {
			largeur = largeur / 2;
		}
		if (largeur == 1) {
			T temp = unite(liste_coeffs[0],true);
			poly(std::vector<T>{-liste_coeffs[debut_gauche], temp});
			gauche = NULL;
			droite = NULL;
			debut_droit = debut_gauche;
		}
		else {
			gauche = new bloc_interpolation<T>(n_max, liste_coeffs, debut_gauche, largeur / 2);
			droite = new bloc_interpolation<T>(n_max, liste_coeffs, debut_gauche + (largeur/2), largeur / 2);
			poly = (gauche->poly * droite->poly);
			debut_droit = debut_gauche + largeur / 2;
		}
	}

	/*
	bloc_interpolation(int n_max, T* liste_coeffs) { //n_max taille de la liste
		int largeur = 1;
		while (largeur < n_max)
			largeur = largeur * 2;
		bloc_interpolation(n_max, liste_coeffs, 0, largeur);
	}*/

	~bloc_interpolation() {
		if (gauche != NULL)
			delete gauche;
		if (droite != NULL)
			delete droite;
	}

	polynome<T> recherche(int n) { //n de 0 à ....
		if (gauche == NULL) {
			if (n == debut_droit) {
				T temp = unite(poly.coeffs[0],true);
				return polynome<T>(temp);
			}
			else
				return poly;
		}
		else{
			if (n >= debut_droit) {
				return (gauche->poly) * (droite->recherche(n));
			}
			else {
				return (droite->poly) * (gauche->recherche(n));
			}
		}
	}


	bloc_interpolation<T>* gauche;
	bloc_interpolation<T>* droite;

	polynome<T> poly;

	int debut_droit; // commence à 0 ...
	int largeur; //2^n
};

template<class T> polynome<T> interpolation(T* liste_x, T* liste_y, int taille) {
	T faux = unite(liste_x[0],false);

	polynome<T> resultat(faux);

	int largeur = 1;
	while (largeur < taille)
		largeur = largeur * 2;
	bloc_interpolation<T> mon_interpolation(taille, liste_x, 0, largeur);

	for (int i(0); i < taille; ++i) {
		polynome<T> interpolation_locale = mon_interpolation.recherche(i);
		interpolation_locale = (liste_y[i] / interpolation_locale(liste_x[i])) * interpolation_locale;
		resultat = resultat + interpolation_locale;
	}

	return resultat;
};


template<class T> polynome<T> interpolation_newton(T* liste_x, T* liste_y, int taille) {

	int  n = taille - 1; // degré
	T* alphas = new T[taille];
	for (int i(0); i < taille; ++i)
		alphas[i] = liste_y[i];

	for (int i = 0; i < n; ++i)
		for (int j = n; j > i; --j)
			alphas[j] = (alphas[j] - alphas[j - 1]) / (liste_x[j] - liste_x[j - i - 1]);



//	polynome<T> resultat(liste_y[0]);
	T vrai = unite( liste_x[0],true);

	polynome<T> multiplication(vrai);
	polynome<T> resultat = unite(polynome<T>(vrai),false);

	for (int i(taille - 1); i >= 0; --i) {
		//		multiplication = multiplication * polynome<T>(-liste_x[i], vrai);
		//		resultat = resultat + (alphas[i] * multiplication) ;
		resultat = resultat * polynome<T>(-liste_x[i], vrai) + polynome<T>(alphas[i]);
	}

	delete alphas;

	return resultat;
}