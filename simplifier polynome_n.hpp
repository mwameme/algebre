#pragma once
#include "polynome_n_iter.hpp"
#include "matrice.hpp"
#include <vector>

template<class T> class scalaire_vecteur {
public:

	scalaire_vecteur() {};


	scalaire_vecteur(int n) {
		vecteur = std::vector<T>(n);
	};

	scalaire_vecteur<T>& operator=(int i) { //i commence à 0.
		T faux = unite(scalaire, false);
		T vrai = unite(scalaire, true);
		for (int j(0); j < vecteur.size(); ++j)
			vecteur[j] = faux;
		vecteur[i] = vrai;
		scalaire = faux;
		return *this;
	};

	scalaire_vecteur<T>& operator=(T element) {
		T faux = unite(scalaire, false);
		for (int i(0); i < vecteur.size(); ++i)
			vecteur[i] = faux;
		scalaire = element;
		return *this;
	};

	friend scalaire_vecteur<T> operator*(scalaire_vecteur<T> const& gauche, scalaire_vecteur<T> const& droit) {//gauche est scalaire, et droit est vecteur ou scalaire
		scalaire_vecteur<T> resultat(droit.vecteur.size());
		resultat.scalaire = gauche.scalaire * droit.scalaire;
		for (int i(0); i < resultat.vecteur.size(); ++i)
			resultat[i] = gauche.scalaire * droit.vecteur[i];
		return resultat;
	};

	friend scalaire_vecteur<T> operator+(scalaire_vecteur<T> const& gauche, scalaire_vecteur<T> const& droit) {//additionne membre à membre
		if (gauche.vecteur.size() != droit.vecteur.size())
			throw std::domain_error("addition de scalaire_vecteur : les dimensions ne correspondent pas");

		scalaire_vecteur<T> resultat(droit.vecteur.size());
		resultat.scalaire = gauche.scalaire + droit.scalaire;
		for (int i(0); i < resultat.vecteur.size(); ++i)
			resultat[i] = gauche.vecteur[i] + droit.vecteur[i];
		return resultat;
	};

	T scalaire;
	std::vector<T> vecteur;
};

template<class T> polynome_n_iter<T> simplifier(polynome_n_iter<T> const& num, polynome_n_iter<T> const& denom) {
	int n = num.coeffs.puissance;
	if (num.coeffs.puissance != denom.coeffs.puissance)
		throw std::domain_error("simplification de polynomes : n_var ne correspond pas");

	T faux = unite(num.coeffs.data[0], false);

	polynome_n_iter<T> vide(n, faux, num.noms); //renvoit le polynome nul si le calcul n'aboutit pas.

	std::vector<int> degres(n); //calcule les degres du polynome resultat
	for (int i(0); i < n; ++i)
		degres[i] = num.coeffs.dimensions[i] - denom.coeffs.dimensions[i];
	for (int i(0); i < n; ++i)
		if (degres[i] < 0)
			return vide;

	int puissance = 1; //nombre d'éléments dans le polynome retourné ... sert à scalaire_vecteur<T>
	for (int i(0); i < n; ++i)
		puissance *= (degres[i]+1);

	scalaire_vecteur<T> faux_T(0);
	faux_T = faux;

	vecteur_n<scalaire_vecteur<T>> denom_vec_T(denom.coeffs.dimensions, faux_T); //même dimension que denom. Mais de type différent ...
	polynome_n_iter<scalaire_vecteur<T>> denom_T();
	denom_T.coeffs = denom_vec_T;
	denom_T.noms = denom.noms;

	for (int i(0); i < denom.coeffs.data.size(); ++i)
		denom_T.coeffs.data[i] = denom.coeffs.data[i]; //scalaire ... on recopie membre à membre. Type T. (modifie le scalaire)
	
	faux_T.vecteur.resize(puissance); //on met à jour la taille du vecteur ... taille(faux_T) = taille(resultat_T)
	polynome_n_iter<scalaire_vecteur<T>> resultat_T(degres, faux_T, denom.noms); //on connait le degré du résultat à l'avance ...
	for (int i(0); i < resultat_T.coeffs.data.size(); ++i)
		resultat_T.coeffs.data[i] = i; //on génère les vecteurs de la "base"

	//on a le resultat, et le denom ... on les multiplie et extrait les équations !
	polynome_n_iter<scalaire_vecteur<T>> num_T = denom_T * resultat_T; //normalement les dimensions correspondent. Car pas de simplification (je parle de data ...)

	//construire la matrice. Il y a autant d'équations que le nombre(num_T) = nombre(num) : taille_l
	//le nombre de variables est nombre(resultat_T) : taille_c
	matrice<T> m_matrice(num_T.coeffs.data.size(), resultat_T.coeffs.data.size(), faux); //matrice remplie avec faux.
	std::vector<T> Y(m_matrice.taille_l); //pour resoudre ...
	for (int i(0); i < m_matrice.taille_l; ++i) {//pour chacune des lignes
		Y[i] = num.coeffs.data[i];
		for (int j(0); j < m_matrice.taille_c; ++j)
			m_matrice[i][j] = num_T.coeffs.data[i].vecteur[j];
	}

	std::vector<T> X = m_matrice.resoudre(Y);
	if (X.size() == 0) //non simplifiable
		return polynome_n_iter<T>(0,faux,NULL);

	polynome_n_iter<T> resultat(degres, faux, denom.noms);
	for (int i(0); i < resultat.coeffs.data.size(); ++i)
		resultat.coeffs.data[i] = X[i]; //on sait que e_i = X_i ... c'était le but.

	return resultat;
}