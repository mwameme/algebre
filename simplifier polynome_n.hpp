#pragma once
#include "polynome_n_iter.hpp"
#include "matrice.hpp"
#include "rationnel.hpp"
#include "unite.hpp"
#include <vector>
#include <iostream>
#include "n_for.hpp"
#include <algorithm>

#include "swap_T.hpp"

template<class T> class scalaire_vecteur {
public:

	T scalaire;
	std::vector<T> vecteur;

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

	scalaire_vecteur<T>& operator+=(scalaire_vecteur<T> const& temp) {
		return (*this = (*this + temp));
	}

	scalaire_vecteur<T>& operator*=(scalaire_vecteur<T> const& temp) {
		return (*this = (*this * temp));
	}

	template<class U>
	friend scalaire_vecteur<T> operator*(U const& scalaire, const scalaire_vecteur<T>& vec) {
		scalaire_vecteur<T> result(*this);
		result.scalaire *= scalaire;
		for (int i(0); i < result.vecteur.size(); ++i)
			result.vecteur[i] *= scalaire;
		return result;
	}

	template<class U>
	scalaire_vecteur<T>& operator*=(const U& scalaire_) {
		scalaire *= scalaire_;
		for (int i(0); i < vecteur.size(); ++i)
			vecteur[i] *= scalaire_;
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
#ifdef ALGEBRA_USE_EXCEPTION
		if (gauche.vecteur.size() != droit.vecteur.size())
			throw std::domain_error("addition de scalaire_vecteur : les dimensions ne correspondent pas");
#endif

		scalaire_vecteur<T> resultat(droit.vecteur.size());
		resultat.scalaire = gauche.scalaire + droit.scalaire;
		for (int i(0); i < resultat.vecteur.size(); ++i)
			resultat.vecteur[i] = gauche.vecteur[i] + droit.vecteur[i];
		return resultat;
	};

	operator bool() const {
		if ((bool)scalaire)
			return true;
		for (int i(0); i < vecteur.size(); ++i)
			if ((bool)vecteur[i])
				return true;
		return false;
	};

	friend void swap(scalaire_vecteur<T>& gauche, scalaire_vecteur<T>& droit) {
		std::swap(gauche.vecteur, droit.vecteur);
		swap_F(gauche.scalaire, droit.scalaire);
	};

};



template<class T> polynome_n_iter<T> simplifier_poly(polynome_n_iter<T> const& num, polynome_n_iter<T> const& denom) {
	int n = num.coeffs.puissance;
#ifdef ALGEBRA_USE_EXCEPTION
	if (num.coeffs.puissance != denom.coeffs.puissance)
		throw std::domain_error("simplification de polynome_n_iter : n_var ne correspond pas");
	if (((num.scalaire) && (!denom.scalaire)) || ((denom.scalaire) && (!num.scalaire)))
		throw std::domain_error("simplification de polynome_n_iter : scalaire / polynome");
#endif
	if (num.scalaire)
		return polynome_n_iter<T>(0, num.coeffs.data[0] / denom.coeffs.data[0]);


	T faux = unite(num.coeffs.data[0], false);

	//	polynome_n_iter<T> vide(n, faux, num.noms); //renvoit le polynome nul si le calcul n'aboutit pas.

	std::vector<int> degres(n); //calcule les degres du polynome resultat
	for (int i(0); i < n; ++i)
		degres[i] = num.coeffs.dimensions[i] - denom.coeffs.dimensions[i];
	for (int i(0); i < n; ++i)
		if (degres[i] < 0)
			return polynome_n_iter<T>(0, faux);

	int puissance = 1; //nombre d'éléments dans le polynome retourné ... sert à scalaire_vecteur<T>
	for (int i(0); i < n; ++i)
		puissance *= (degres[i] + 1);

	scalaire_vecteur<T> faux_T(0);
	faux_T = faux;

	//denominateur sous forme scalaire_vecteur
	vecteur_n<scalaire_vecteur<T>> denom_vec_T(denom.coeffs.dimensions, faux_T); //même dimension que denom. Mais de type différent ...
	polynome_n_iter<scalaire_vecteur<T>> denom_T;
	denom_T.coeffs = denom_vec_T;
	denom_T.scalaire = false;

	for (int i(0); i < denom.coeffs.data.size(); ++i)
		denom_T.coeffs.data[i] = denom.coeffs.data[i]; //scalaire ... on recopie membre à membre. Type T. (modifie le scalaire)

	//on fait une optimisation de coût n^2 (calcul total de coût n^3)
	//get positions_num
	std::vector<std::vector<int>> positions_num;
	positions_num.reserve(num.coeffs.data.size());
	for (n_for iter(num.coeffs.dimensions); (bool)iter; ++iter)
		if ((bool)num.coeffs.data[iter.position])
			positions_num.push_back(iter.positions);

	//get positions_denom
	std::vector<std::vector<int>> positions_denom;
	positions_denom.reserve(denom.coeffs.data.size());
	for (n_for iter(denom.coeffs.dimensions); (bool)iter; ++iter)
		if ((bool)denom.coeffs.data[iter.position])
			positions_denom.push_back(iter.positions);

	//get positions_R
	std::vector<std::vector<int>> positions_R;
	positions_R.reserve(puissance);
	for(int i(0);i<positions_num.size();++i)
		for (int j(0); j < positions_denom.size(); ++j) {
			std::vector<int> vec_temp(n, 0);
			bool test = true;
			for (int k(0); k < n; ++k) {
				int temp = positions_num[i][k] - positions_denom[j][k];
				if (temp < 0) {
					test = false;
					break;
				}
				vec_temp[k] = temp;
			}
			positions_R.push_back(vec_temp);
		}
	//trier
	std::sort(positions_R.begin(), positions_R.end());

	//unique (à la main)
	std::vector<std::vector<int>> positions_R2;
	positions_R2.reserve(positions_R.size());
	positions_R2.push_back(positions_R[0]);//vérifier au moins 1 ...
	std::vector<int> temp_vec = positions_R2[0];
	for (int i(0); i < positions_R.size(); ++i) {
		if (temp_vec == positions_R[i])
			continue;
		temp_vec = positions_R[i];
		positions_R2.push_back(temp_vec);
	}

	//unique a été fait !!! compter, et créer resultat_T

	faux_T.vecteur.resize(positions_R2.size()); //on met à jour la taille du vecteur ... taille(faux_T) = taille(resultat_T)
	polynome_n_iter<scalaire_vecteur<T>> resultat_T(degres, faux_T); //on connait le degré du résultat à l'avance ...
	for (int i(0); i < positions_R2.size(); ++i)
		resultat_T.coeffs.data[resultat_T.coeffs.position(positions_R2[i])] = i; //on génère les vecteurs de la "base"

	//on a le resultat, et le denom ... on les multiplie et extrait les équations !
	polynome_n_iter<scalaire_vecteur<T>> num_T = denom_T * resultat_T; //normalement les dimensions correspondent. Car pas de simplification (je parle de data ...)

	//construire la matrice. Il y a autant d'équations que le nombre(num_T) = nombre(num) : taille_l
	//le nombre de variables est nombre(resultat_T) : taille_c

	//on implémente num = num_T. On ne garde que : (bool) num.coeffs.data[i] == true
	matrice<T> m_matrice(positions_num.size(), positions_R2.size(), faux); //matrice remplie avec faux.
	std::vector<T> Y(m_matrice.taille_l); //pour resoudre ...
	int j = 0;
	for (int i(0); i < num.coeffs.data.size(); ++i) {//pour chacune des lignes
		if ((bool)num.coeffs.data[i]) {
			Y[j] = num.coeffs.data[i];
			for (int k(0); k < m_matrice.taille_c; ++k)
				m_matrice.coeffs[j][k] = num_T.coeffs.data[i].vecteur[k];
			++j;
		}
	}

	std::vector<T> X = m_matrice.resoudre(Y); //taille : positions_R2.size()
	if (X.size() == 0) //non simplifiable
		return polynome_n_iter<T>(0, faux); //exception ?

	polynome_n_iter<T> resultat(degres, faux);
	for (int i(0); i < positions_R2.size(); ++i)
		resultat.coeffs.data[resultat.coeffs.position(positions_R2[i])] = X[i]; //on sait que e_i = X_i ... c'était le but.

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