#pragma once
#include "polynome_n_iter.hpp"
#include "matrice.hpp"
#include "rationnel.hpp"
#include "unite.hpp"
#include <vector>
#include <iostream>


template<class T> class scalaire_vecteur {
public:

	scalaire_vecteur() : vecteur(0) { };


	scalaire_vecteur(int n) : vecteur(n) { };

	scalaire_vecteur(scalaire_vecteur<T> const& temp) {
		scalaire = temp.scalaire;
		vecteur = temp.vecteur;
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

	scalaire_vecteur<T>& operator=(T const& element) {
		T faux = unite(element, false);
		for (int i(0); i < vecteur.size(); ++i)
			vecteur[i] = faux;
		scalaire = element;
		return *this;
	};

	scalaire_vecteur<T>& operator=(scalaire_vecteur<T> const& temp) {
		if (this == &temp)
			return *this;
		scalaire = temp.scalaire;
		vecteur = temp.vecteur;
		return *this;
	};

	friend scalaire_vecteur<T> operator*(scalaire_vecteur<T> const& gauche, scalaire_vecteur<T> const& droit) {//gauche est scalaire, et droit est vecteur ou scalaire
		if ((!(bool)gauche.scalaire) && ((bool)droit.scalaire)) //car on calcule scalaire * sous_ev ...
			return droit * gauche;

		scalaire_vecteur<T> resultat(droit.vecteur.size());
		resultat.scalaire = gauche.scalaire * droit.scalaire;
		for (int i(0); i < resultat.vecteur.size(); ++i)
			resultat.vecteur[i] = gauche.scalaire * droit.vecteur[i];
		return resultat;
	};

	friend scalaire_vecteur<T> operator+(scalaire_vecteur<T> const& gauche, scalaire_vecteur<T> const& droit) {//additionne membre à membre
		if (gauche.vecteur.size() == 0)
			return droit;
		if (droit.vecteur.size() == 0)
			return gauche;

		if (gauche.vecteur.size() != droit.vecteur.size())
			throw std::domain_error("addition de scalaire_vecteur : les dimensions ne correspondent pas");

		scalaire_vecteur<T> resultat(droit.vecteur.size());
		resultat.scalaire = gauche.scalaire + droit.scalaire;
		for (int i(0); i < resultat.vecteur.size(); ++i)
			resultat.vecteur[i] = gauche.vecteur[i] + droit.vecteur[i];
		return resultat;
	};

	operator bool() const {
		if ((vecteur.size() == 0) && (!((bool)scalaire)))
			return false;
		return true;
	};

	T scalaire;
	std::vector<T> vecteur;
};

template<class T> polynome_n_iter<T> simplifier_poly(polynome_n_iter<T> const& num, polynome_n_iter<T> const& denom) {
	int n = num.coeffs.puissance;
	if (num.coeffs.puissance != denom.coeffs.puissance)
		throw std::domain_error("simplification de polynomes : n_var ne correspond pas");
	if (((num.scalaire) && (!denom.scalaire)) || ((denom.scalaire) && (!num.scalaire)))
		throw std::domain_error("multiplication de polynome_n_iter : scalaire * polynome");

	if (num.scalaire)
		return polynome_n_iter<T>(0, num.coeffs.data[0] / denom.coeffs.data[0], NULL);


	T faux = unite(num.coeffs.data[0], false);

	//	polynome_n_iter<T> vide(n, faux, num.noms); //renvoit le polynome nul si le calcul n'aboutit pas.

	std::vector<int> degres(n); //calcule les degres du polynome resultat
	for (int i(0); i < n; ++i)
		degres[i] = num.coeffs.dimensions[i] - denom.coeffs.dimensions[i];
	for (int i(0); i < n; ++i)
		if (degres[i] < 0)
			return polynome_n_iter<T>(0, faux, NULL);

	int puissance = 1; //nombre d'éléments dans le polynome retourné ... sert à scalaire_vecteur<T>
	for (int i(0); i < n; ++i)
		puissance *= (degres[i] + 1);

	scalaire_vecteur<T> faux_T(0);
	faux_T = faux;

	//denominateur sous forme scalaire_vecteur
	vecteur_n<scalaire_vecteur<T>> denom_vec_T(denom.coeffs.dimensions, faux_T); //même dimension que denom. Mais de type différent ...
	polynome_n_iter<scalaire_vecteur<T>> denom_T;
	denom_T.coeffs = denom_vec_T;
	denom_T.noms = denom.noms;
	denom_T.scalaire = false;

	for (int i(0); i < denom.coeffs.data.size(); ++i)
		denom_T.coeffs.data[i] = denom.coeffs.data[i]; //scalaire ... on recopie membre à membre. Type T. (modifie le scalaire)



	faux_T.vecteur.resize(puissance); //on met à jour la taille du vecteur ... taille(faux_T) = taille(resultat_T)
	polynome_n_iter<scalaire_vecteur<T>> resultat_T(degres, faux_T, denom.noms); //on connait le degré du résultat à l'avance ...
	for (int i(0); i < resultat_T.coeffs.data.size(); ++i)
		resultat_T.coeffs.data[i] = i; //on génère les vecteurs de la "base"

	//on a le resultat, et le denom ... on les multiplie et extrait les équations !
	polynome_n_iter<scalaire_vecteur<T>> num_T = denom_T * resultat_T; //normalement les dimensions correspondent. Car pas de simplification (je parle de data ...)

	//OK
//	std::cout << "tailles : " << denom_T.coeffs.data.size() << " ; " << denom.coeffs.data.size() << std::endl;
//	std::cout << faux_T.vecteur.size() << " ; " << resultat_T.coeffs.data.size() << " ; " << num_T.coeffs.data.size() << " ; " << num.coeffs.data.size() << std::endl;

	//construire la matrice. Il y a autant d'équations que le nombre(num_T) = nombre(num) : taille_l
	//le nombre de variables est nombre(resultat_T) : taille_c
	matrice<T> m_matrice(num_T.coeffs.data.size(), resultat_T.coeffs.data.size(), faux); //matrice remplie avec faux.
	std::vector<T> Y(m_matrice.taille_l); //pour resoudre ...
	for (int i(0); i < m_matrice.taille_l; ++i) {//pour chacune des lignes
		Y[i] = num.coeffs.data[i];
		for (int j(0); j < m_matrice.taille_c; ++j)
			m_matrice.coeffs[i][j] = num_T.coeffs.data[i].vecteur[j];
	}

	//	std::cout << m_matrice << std::endl;
		/*
		for (int i(0); i < m_matrice.taille_l; ++i)
			if ((bool)Y[i]) {
				for (int j(0); j < m_matrice.taille_c; ++j)
					std::cout << m_matrice.coeffs[i][j];
				std::cout << std::endl;
			}
		long question;
		std::cin >> question;
		*/

	std::vector<T> X = m_matrice.resoudre(Y);
	if (X.size() == 0) //non simplifiable
		return polynome_n_iter<T>(0, faux, NULL);

	polynome_n_iter<T> resultat(degres, faux, denom.noms);
	for (int i(0); i < resultat.coeffs.data.size(); ++i)
		resultat.coeffs.data[i] = X[i]; //on sait que e_i = X_i ... c'était le but.

	resultat.simplifier_2();

	return resultat;
};


template<class T> class Simplifier_T {
public:
	static void Simplifier(T & x) {
		return;
	};
};

template<class T> class Simplifier_T<rationnel<polynome_n_iter<T>>> {
public:
	static void Simplifier(rationnel<polynome_n_iter<T>> & x) {
		polynome_n_iter<T> simple = simplifier_poly<T>(x.numerateur, x.denominateur);
		if (simple.scalaire == true)
			return ;
		T vrai = unite(simple.coeffs.data[0], true);
		polynome_n_iter<T> vrai_T(simple.coeffs.puissance, vrai, simple.noms);
		x = rationnel<polynome_n_iter<T>>(simple, vrai_T);
		return;
	};
};

template<class T> inline void Simplifier_frac_poly(T & x) {
	Simplifier_T<T>::Simplifier(x);
};

template<class T> inline void simplifier_poly(polynome<T> &poly) {
	for (int i(0); i < poly.coeffs.size(); ++i)
		Simplifier_frac_poly(poly.coeffs[i]);
	poly.getDegre();
	return;
};