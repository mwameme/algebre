#pragma once
#include <vector>
#include <iostream>
#include <exception>
#include <initializer_list>

#include "entete objets.hpp"
#include "unite.hpp"
#include "fact_for.hpp"

//#include <typeinfo> 

template<class T> class sous_ev;


// matrice:
//cas type=1 ou 2 (static_assert)
//cas type=0, approx =1 ou 0
// 
//rationnel : cas type=1 ou 2
// 
//diagonalisation : cas approx=0 ou 1. Static_assert type=0
//

template<class T> class matrice<T,typename std::enable_if_t<type_algebre<T>::type==2>,void> { //sans division : determinant permutation
public:

	friend matrice<T> derivee(const matrice<T>& element) {
		matrice<T> result(element);
		for (int i(0); i < element.taille_l; ++i)
			for (int j(0); j < element.taille_c; ++j)
				result.coeffs[i][j] = derivee(result.coeffs[i][j]);
		return result;
	};

	explicit matrice(int m_taille_l,int m_taille_c, T element) {
		taille_l = m_taille_l;
		taille_c = m_taille_c;
		coeffs = std::vector<std::vector<T>>(m_taille_l, std::vector<T>(m_taille_c,element));
	};

	explicit matrice(int m_taille_l,int m_taille_c)  {
		taille_l = m_taille_l;
		taille_c = m_taille_c;

		coeffs = std::vector<std::vector<T>> (m_taille_l, std::vector<T>(m_taille_c));
	};

	matrice(matrice<T> const& copie) {
		taille_l = copie.taille_l;
		taille_c = copie.taille_c;
		coeffs = copie.coeffs;
	};

	explicit matrice(const std::vector<std::vector<T>>& vec) :  taille_l(vec.size()),taille_c(vec[0].size()) {
		coeffs = vec;
		for(int i(1);i<taille_l;++i)
			if (coeffs[i].size() != taille_c)
				throw std::domain_error("creation de matrice : non-rectangulaire");

	};

	explicit matrice() {};

	template<class U> explicit matrice(std::initializer_list< std::initializer_list<U>> liste) {
		int taille_l = liste.size();
		if (taille_l == 0)
			throw std::domain_error("initialisation de matrice : liste vide");
		int taille_c = liste[0].size();
		if (taille_c == 0)
			throw std::domain_error("initialisation de matrice : nombre de colonnes nul");

		for (int i(1); i < taille_l; ++i)
			if (liste[i].size() != taille_c)
				throw std::domain_error("initialisation de matrice : non-rectangulaire");

		coeffs(taille_l, std::vector<T>(taille_c));
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				coeffs[i][j] = T(liste[i][j]);

		return;
	}

	explicit operator bool() const {
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				if ((bool)coeffs[i][j])
					return true;
		return false;
	};

	template<class U> explicit operator matrice<U>() const {
		matrice<U> result(taille_l,taille_c);
		
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				result.coeffs[i][j] = (U)coeffs[i][j];

		return result;
	};


	matrice<T>& operator=(bool test) { //des 1 sur la diagonale, 0 partout ailleurs.
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				coeffs[i][j] = false;
		int taille = min(taille_l, taille_c);

		if (test)
			for (int i(0); i < taille; ++i)
				coeffs[i][i] = true;

		return *this;

	};

	template<class U> friend matrice<T> operator*(U scalaire, const matrice<T>& m_matrice) {
		matrice<T> result(m_matrice);
		for (int i(0); i < result.taille_l; ++i)
			for (int j(0); j < result.taille_c; ++j)
				result.coeffs[i][j] = scalaire * result.coeffs[i][j];
		return result;
	};

	matrice<T> operator*(const matrice<T>& autre) const {
		if (taille_c != autre.taille_l)
			throw std::domain_error("multiplication de matrices : les dimensions ne coïncident pas");

		matrice<T> result(taille_l,autre.taille_c);
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < autre.taille_c; ++j) {
				T somme = unite(coeffs[0][0],false);
				for (int k(0); k < taille_c; ++k)
					somme = somme + (coeffs[i][k] * autre.coeffs[k][j]);
				result.coeffs[i][j] = somme;
			}
		return result;
	};

	matrice<T> operator+(const matrice<T>& autre) const {
		if ((taille_l != autre.taille_l) || (taille_c != autre.taille_c))
			throw std::domain_error("addition de matrices : les dimensions ne coïncident pas");

		matrice<T> result(taille_l,taille_c);
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				result.coeffs[i][j] = coeffs[i][j] + autre.coeffs[i][j];
		return result;
	};

	friend matrice<T> operator-(matrice<T> const& gauche, matrice<T> const& droit) {
		if ((gauche.taille_l != droit.taille_l) || (gauche.taille_c != droit.taille_c))
			throw std::domain_error("soustraction de matrices : les dimensions ne coïncident pas");

		matrice<T> result(gauche.taille_l, gauche.taille_c);
		for (int i(0); i < gauche.taille_l; ++i)
			for (int j(0); j < gauche.taille_c; ++j)
				result.coeffs[i][j] = gauche.coeffs[i][j] - droit.coeffs[i][j];
		return result;
	};

	friend matrice<T> operator-(matrice<T> const& element) {
		matrice<T> result(element.taille_l,element.taille_c);
		for (int i(0); i < result.taille_l; ++i)
			for (int j(0); j < result.taille_c; ++j)
				result.coeffs[i][j] = - element.coeffs[i][j];
		return result;
	};
	
	T determinant_anneau() const {
		T resultat = unite(coeffs[0][0],false);
		T vrai_ = unite(resultat,true);

		if (taille_l != taille_c)
			throw std::domain_error("determinant : les dimensions ne coïncident pas");


		for (fact_for iter(taille_l); (bool)iter; ++iter) {
			T temp = vrai_;
			for (int i(0); i < taille_l; ++i)
				temp = temp * coeffs[i][iter.permutation[i]];
			temp = iter.signature * temp;
			resultat = resultat + temp;
		}

		return resultat;
	};
	
	inline T determinant() const {
		return determinant_anneau();
	};


	void echangerLigne(int i, int j) {
		if ((i < 0) || (i >= taille_l) || (j < 0) || (j >= taille_l))
			throw std::domain_error("echange de lignes : hors domaine");

		if (i == j)
			return;
		for (int k(0); k < taille_c; ++k) {
			T temp = coeffs[i][k];
			coeffs[i][k] = coeffs[j][k];
			coeffs[j][k] = temp;
		}
		return;
	};

	void ajouterLigne(int i, int j, T coefficient) { // ajouter ligne i * coeffs à la ligne j ...
		if ((i < 0) || (i >= taille_l) || (j < 0) || (j >= taille_l))
			throw std::domain_error("ajout de lignes : hors domaine");
		if (i == j)
			throw std::domain_error("ajout de lignes : même ligne");
		for (int k(0); k < taille_c; ++k) {
			coeffs[j][k] = coeffs[j][k] + (coefficient * coeffs[i][k]);
		}
		return;
	};

	void multiplierLigne(int i, T coefficient) {
		for (int j(0); j < taille_c; ++j)
			coeffs[i][j] = coefficient * coeffs[i][j];
	};

	polynome<T> polynomeCaracteristique() const {
		T mvrai = - unite(coeffs[0][0],true);

		if(taille_l != taille_c)
			throw std::domain_error("polynome caracteristique : les dimensions ne correspondent pas");
		int taille = taille_l;

//		polynome<T> m_poly = polynome<T>(vrai);

		matrice<polynome<T>> m_matrice(taille,taille);

		for (int i(0); i < taille; ++i)
			for (int j(0); j < taille; ++j)
				if (i != j)
					m_matrice.coeffs[i][j] = polynome<T>(coeffs[i][j]);
				else
					m_matrice.coeffs[i][j] = polynome<T>(std::vector<T>{ coeffs[i][j], mvrai });

		polynome<T> det = m_matrice.determinant();

		return det;
	};


	friend std::ostream& operator<<(std::ostream& os, const matrice<T>& element) {
		os << "{ ";
		for (int i(0); i < element.taille_l; ++i) {
			os << "{ ";
			for (int j(0); j < element.taille_c - 1; ++j)
				os << element.coeffs[i][j] << " , ";
			if (i < element.taille_l - 1)
				os << element.coeffs[i][element.taille - 1] << " } ,";
			else
				os << element.coeffs[i][element.taille - 1] << " }";

			if (i < element.taille_l - 1)
				os << std::endl;
		}
		os << "}" << std::endl;
		return os;
	};


	friend std::vector<T> operator*(std::vector<T> const& vec_ligne, matrice<T> const& m_matrice) {
		if ( vec_ligne.size() != m_matrice.taille_l)
			throw std::domain_error("multiplication vecteur*matrice : les dimensions ne correspondent pas");

		T faux_ = unite(m_matrice.coeffs[0][0],false);

		std::vector<T> resultat(m_matrice.taille_c);
		for (int i(0); i < m_matrice.taille_c; ++i) {
			T somme = faux_;
			for (int j(0); j < m_matrice.taille_l; ++j)
				somme = somme + (vec_ligne[j] * m_matrice.coeffs[j][i]);
			resultat[i] = somme;
		}

		return resultat;
	};

	friend std::vector<T> operator*(matrice<T> const& m_matrice, std::vector<T> const& vec_colonne) {
		if (vec_colonne.size() != m_matrice.taille_c)
			throw std::domain_error("multiplication matrice*vecteur : les dimensions ne correspondent pas");

		T faux_ = unite(m_matrice.coeffs[0][0],false);

		std::vector<T> resultat(m_matrice.taille_l, faux_);
		for (int i(0); i < m_matrice.taille_l; ++i) {
			T somme = faux_;
			for (int j(0); j < m_matrice.taille_c; ++j)
				somme = somme + (m_matrice.coeffs[i][j] * vec_colonne[j]);
			resultat[i] = somme;
		}
		return resultat;
	};

	int taille_l;
	int taille_c;
	std::vector < std::vector< T>> coeffs;
};


template<class T> class matrice<T, typename std::enable_if_t<type_algebre<T>::type == 1>, void> { //division avec reste
public:

	friend matrice<T> derivee(const matrice<T>& element) {
		matrice<T> result(element);
		for (int i(0); i < element.taille_l; ++i)
			for (int j(0); j < element.taille_c; ++j)
				result.coeffs[i][j] = derivee(result.coeffs[i][j]);
		return result;
	};

	explicit matrice(int m_taille_l, int m_taille_c, T element) {
		taille_l = m_taille_l;
		taille_c = m_taille_c;
		coeffs = std::vector<std::vector<T>>(m_taille_l, std::vector<T>(m_taille_c, element));
	};

	explicit matrice(int m_taille_l, int m_taille_c) { //déconseillé à l'utilisation, sauf si vous savez ce que vous faites.*
		taille_l = m_taille_l;
		taille_c = m_taille_c;

		coeffs = std::vector<std::vector<T>>(m_taille_l, std::vector<T>(m_taille_c));
	};

	matrice(matrice<T> const& copie) {
		taille_l = copie.taille_l;
		taille_c = copie.taille_c;
		coeffs = copie.coeffs;
	};

	explicit matrice(const std::vector<std::vector<T>>& vec) : coeffs(vec), taille_l(vec.size()), taille_c(vec[0].size()) {
		for (int i(1); i < taille_l; ++i)
			if (coeffs[i].size() != taille_c)
				throw std::domain_error("creation de matrice : non-rectangulaire");

	};

	explicit matrice() : taille_l(0), taille_c(0), coeffs(0, std::vector<T>(0)) {};

	explicit operator bool() const {
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				if ((bool)coeffs[i][j])
					return true;
		return false;
	};

	template<class U> explicit operator matrice<U>() const {
		matrice<U> result;
		result.taille_l = taille_l;
		result.taille_c = taille_c;
		result.coeffs.resize(taille_l);
		for (int i(0); i < taille_l; ++i)
			result.coeffs[i].resize(taille_c);


		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				result.coeffs[i][j] = (U)coeffs[i][j];

		return result;
	};


	matrice<T>& operator=(bool test) {
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				coeffs[i][j] = false;
		int taille = min(taille_l, taille_c);
		if (test)
			for (int i(0); i < taille; ++i)
				coeffs[i][i] = true;
		return *this;

	};

	template<class U> friend matrice<T> operator*(U scalaire, const matrice<T>& m_matrice) {
		matrice<T> result(m_matrice);
		for (int i(0); i < result.taille_l; ++i)
			for (int j(0); j < result.taille_c; ++j)
				result.coeffs[i][j] = scalaire * result.coeffs[i][j];
		return result;
	};

	matrice<T> operator*(const matrice<T>& autre) const {
		T faux_ = unite(coeffs[0][0],false);

		if (taille_c != autre.taille_l)
			throw std::domain_error("multiplication de matrices : les dimensions ne coïncident pas");

		matrice<T> result(taille_l, autre.taille_c);
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < autre.taille_c; ++j) {
				T somme(faux_);
				for (int k(0); k < taille_c; ++k)
					somme = somme + (coeffs[i][k] * autre.coeffs[k][j]);
				result.coeffs[i][j] = somme;
			}
		return result;
	};

	friend matrice<T> operator+(const matrice<T>& gauche, matrice<T> const& droit) {
		if ((gauche.taille_l != droit.taille_l) || (gauche.taille_c != droit.taille_c))
			throw std::domain_error("addition de matrices : les dimensions ne coïncident pas");

		matrice<T> result(gauche.taille_l, gauche.taille_c);
		for (int i(0); i < gauche.taille_l; ++i)
			for (int j(0); j < gauche.taille_c; ++j)
				result.coeffs[i][j] = gauche.coeffs[i][j] + droit.coeffs[i][j];
		return result;
	};

	matrice<T> operator-(const matrice<T>& autre) const {
		if ((taille_l != autre.taille_l) || (taille_c != autre.taille_c))
			throw std::domain_error("addition de matrices : les dimensions ne coïncident pas");

		matrice<T> result(taille_l, taille_c);
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				result.coeffs[i][j] = coeffs[i][j] - autre.coeffs[i][j];
		return result;
	};

	friend 	matrice<T> operator-(const matrice<T>& element) {
		matrice<T> result(element.taille_l, element.taille_c);
		for (int i(0); i < result.taille_l; ++i)
			for (int j(0); j < result.taille_c; ++j)
				result.coeffs[i][j] = -element.coeffs[i][j];
		return result;
	};

	T determinant_anneau() const {
		T resultat = unite(coeffs[0][0],false);
		T vrai_ = unite(resultat,true);
		if(taille_l != taille_c)
			throw std::domain_error("determinant matrice : les dimensions ne coïncident pas");
		int taille = taille_l;

		for (fact_for iter(taille); (bool)iter; ++iter) {
			T temp = vrai_;
			for (int i(0); i < taille; ++i)
				temp = temp * coeffs[i][iter.permutation[i]];
			temp = iter.signature * temp;
			resultat = resultat + temp;
		}

		return resultat;
	};

	
	T determinant() const {
		matrice<rationnel<T>> m_matrice(taille_l,taille_l);
		T vrai_ = unite(coeffs[0][0],true);

		if (taille_l != taille_c)
			throw std::domain_error("determinant matrice : les dimensions ne coïncident pas");
		int taille = taille_l;


		for (int i(0); i < taille; ++i)
			for (int j(0); j < taille; ++j)
				m_matrice.coeffs[i][j] = rationnel<T>(coeffs[i][j], vrai_); //passer par les rationnels. vrai_ optionnel

		rationnel<T> det = m_matrice.determinant();

		return det.numerateur / det.denominateur; //simplifier : division avec reste 
	};

	void echangerLigne(int i, int j) {
		if ((i < 0) || (i >= taille_l) || (j < 0) || (j >= taille_l))
			throw std::domain_error("echange de lignes : hors domaine");

		if (i == j)
			return;
		for (int k(0); k < taille_c; ++k) {
			T temp = coeffs[i][k];
			coeffs[i][k] = coeffs[j][k];
			coeffs[j][k] = temp;
		}
		return;
	};

	void ajouterLigne(int i, int j, T coefficient) { // ajouter ligne i * coeffs à la ligne j ...
		if ((i < 0) || (i >= taille_l) || (j < 0) || (j >= taille_l))
			throw std::domain_error("ajout de lignes : hors domaine");
		if (i == j)
			throw std::domain_error("ajout de lignes : même ligne");
		for (int k(0); k < taille_c; ++k) {
			coeffs[j][k] = coeffs[j][k] + (coefficient * coeffs[i][k]);
		}
		return;
	};

	void multiplierLigne(int i, T coefficient) {
		for (int j(0); j < taille_c; ++j)
			coeffs[i][j] = coefficient * coeffs[i][j];
	};

	polynome<T> polynomeCaracteristique() const {
		T mvrai = -unite(coeffs[0][0],true);

		if (taille_l != taille_c)
			throw std::domain_error("determinant matrice : les dimensions ne coïncident pas");
		int taille = taille_l;


		matrice<polynome<T>> m_matrice(taille, taille);
		for (int i(0); i < taille; ++i)
			for (int j(0); j < taille; ++j)
				if (i != j)
					m_matrice.coeffs[i][j] = polynome<T>(coeffs[i][j]);
				else
					m_matrice.coeffs[i][j] = polynome<T>(std::vector<T>{ coeffs[i][j], mvrai });

		polynome<T> det = m_matrice.determinant(); //passe par le determinant_anneau : sans utiliser la division.
		return det;
	};

	friend std::ostream& operator<<(std::ostream& os, const matrice<T>& element) {
		os << "{ ";
		for (int i(0); i < element.taille_l; ++i) {
			os << "{ ";
			for (int j(0); j < element.taille_c - 1; ++j)
				os << element.coeffs[i][j] << " , ";
			if (i < element.taille_l - 1)
				os << element.coeffs[i][element.taille - 1] << " } ,";
			else
				os << element.coeffs[i][element.taille - 1] << " }";

			if (i < element.taille_l - 1)
				os << std::endl;
		}
		os << "}" << std::endl;
		return os;
	};


	friend std::vector<T> operator*(std::vector<T> const& vec_ligne, matrice<T> const& m_matrice) {
		if (vec_ligne.size() != m_matrice.taille_l)
			throw std::domain_error("multiplication vecteur*matrice : les dimensions ne correspondent pas");

		T faux_ = unite(m_matrice.coeffs[0][0],false);

		std::vector<T> resultat(m_matrice.taille_c);
		for (int i(0); i < m_matrice.taille_c; ++i) {
			T somme = faux_;
			for (int j(0); j < m_matrice.taille_l; ++j)
				somme = somme + (vec_ligne[j] * m_matrice.coeffs[j][i]);
			resultat[i] = somme;
		}

		return resultat;
	};

	friend std::vector<T> operator*(matrice<T> const& m_matrice, std::vector<T> const& vec_colonne) {
		if (vec_colonne.size() != m_matrice.taille_c)
			throw std::domain_error("multiplication matrice*vecteur : les dimensions ne correspondent pas");

		T faux_ = unite(m_matrice.coeffs[0][0],false);

		std::vector<T> resultat(m_matrice.taille_l, faux_);
		for (int i(0); i < m_matrice.taille_l; ++i) {
			T somme = faux_;
			for (int j(0); j < m_matrice.taille_c; ++j)
				somme = somme + (m_matrice.coeffs[i][j] * vec_colonne[j]);
			resultat[i] = somme;
		}
		return resultat;
	};

	int taille_l;
	int taille_c;
	std::vector < std::vector< T>> coeffs;
};


template<class T> class matrice<T, typename std::enable_if_t<type_algebre<T>::type == 0>, typename std::enable_if_t<type_algebre<T>::approx==0>> { //division, exacte
public:

	friend matrice<T> derivee(const matrice<T>& element) {
		matrice<T> result(element);
		for (int i(0); i < element.taille_l; ++i)
			for (int j(0); j < element.taille_c; ++j)
				result.coeffs[i][j] = derivee(result.coeffs[i][j]);
		return result;
	};

	explicit matrice(int m_taille_l, int m_taille_c, T element) {
		taille_l = m_taille_l;
		taille_c = m_taille_c;
		coeffs = std::vector<std::vector<T>>(m_taille_l, std::vector<T>(m_taille_c, element));
	};

	explicit matrice(int m_taille_l, int m_taille_c) { //déconseillé à l'utilisation, sauf si vous savez ce que vous faites.*
		taille_l = m_taille_l;
		taille_c = m_taille_c;

		coeffs = std::vector<std::vector<T>>(m_taille_l, std::vector<T>(m_taille_c));
	};

	matrice(matrice<T> const& copie) {
		taille_l = copie.taille_l;
		taille_c = copie.taille_c;
		coeffs = copie.coeffs;
	};

	explicit matrice(const std::vector<std::vector<T>>& vec) : coeffs(vec), taille_l(vec.size()), taille_c(vec[0].size()) {
		for (int i(1); i < taille_l; ++i)
			if (coeffs[i].size() != taille_c)
				throw std::domain_error("creation de matrice : non-rectangulaire");

	};

	explicit matrice() : taille_l(0), taille_c(0), coeffs(0, std::vector<T>(0)) {}; //toujours remplire les matrices ...

	explicit operator bool() const {
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				if ((bool)coeffs[i][j])
					return true;
		return false;
	};

	template<class U> explicit operator matrice<U>() const {
		matrice<U> result;
		result.taille_l = taille_l;
		result.taille_c = taille_c;
		result.coeffs.resize(taille_l);
		for (int i(0); i < taille_l; ++i)
			result.coeffs[i].resize(taille_c);


		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				result.coeffs[i][j] = (U)coeffs[i][j];

		return result;
	};


	matrice<T>& operator=(bool test) {
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				coeffs[i][j] = false;
		int taille = min(taille_l, taille_c);
		if (test)
			for (int i(0); i < taille; ++i)
				coeffs[i][i] = true;
		return *this;

	};

	template<class U> friend matrice<T> operator*(U scalaire, const matrice<T>& m_matrice) {
		matrice<T> result(m_matrice);
		for (int i(0); i < result.taille_l; ++i)
			for (int j(0); j < result.taille_c; ++j)
				result.coeffs[i][j] = scalaire * result.coeffs[i][j];
		return result;
	};

	matrice<T> operator*(const matrice<T>& autre) const {
		T faux_ = unite(coeffs[0][0],false);

		if (taille_c != autre.taille_l)
			throw std::domain_error("multiplication de matrices : les dimensions ne coïncident pas");

		matrice<T> result(taille_l, autre.taille_c);
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < autre.taille_c; ++j) {
				T somme(faux_);
				for (int k(0); k < taille_c; ++k)
					somme = somme + (coeffs[i][k] * autre.coeffs[k][j]);
				result.coeffs[i][j] = somme;
			}
		return result;
	};

	matrice<T> operator+(const matrice<T>& autre) const {
		if ((taille_l != autre.taille_l) || (taille_c != autre.taille_c))
			throw std::domain_error("addition de matrices : les dimensions ne coïncident pas");

		matrice<T> result(taille_l, taille_c);
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				result.coeffs[i][j] = coeffs[i][j] + autre.coeffs[i][j];
		return result;
	};

	matrice<T> operator-(const matrice<T>& autre) const {
		if ((taille_l != autre.taille_l) || (taille_c != autre.taille_c))
			throw std::domain_error("addition de matrices : les dimensions ne coïncident pas");

		matrice<T> result(taille_l, taille_c);
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				result.coeffs[i][j] = coeffs[i][j] - autre.coeffs[i][j];
		return result;
	};

	friend 	matrice<T> operator-(const matrice<T>& element) {
		matrice<T> result(element.taille_l, element.taille_c);
		for (int i(0); i < result.taille_l; ++i)
			for (int j(0); j < result.taille_c; ++j)
				result.coeffs[i][j] = -element.coeffs[i][j];
		return result;
	};

	T determinant_anneau() {
		T resultat = unite( coeffs[0][0],false);
		T vrai_ = unite(resultat,true);

		if(taille_l != taille_c)
			throw std::domain_error("determinant de matrice : les dimensions ne coincident pas");
		int taille = taille_l;

		for (fact_for iter(taille); (bool) iter; ++iter) {
			T temp = vrai_;
			for (int i(0); i < taille; ++i)
				temp = temp * coeffs[i][iter.permutation[i]];
			temp = iter.signature * temp;
			resultat = resultat + temp;
		}

		return resultat;
	};

	T determinant() const {
		matrice<T> m_matrice(*this);
		T det = unite(coeffs[0][0],true);

		if (taille_l != taille_c)
			throw std::domain_error("determinant de matrice : les dimensions ne coincident pas");
		int taille = taille_l;

		for (int i(0); i < taille; ++i) {
			int j;
			for (j = i; j < taille; ++j) {
				if ((bool)m_matrice.coeffs[i][j])
					break;
			}
			if (j == taille)
				return (det = false);

			if (i != j) {
				m_matrice.echangerLigne(i, j);
				det = -det;
			}
			for (j = i + 1; j < taille; ++j) {
				m_matrice.ajouterLigne(i, j, ((-m_matrice.coeffs[j][i]) / m_matrice.coeffs[i][i]));
			}

		}

		for (int i(0); i < taille; ++i)
			det = (det * m_matrice.coeffs[i][i]);

		return det;
	};

	
	void echangerLigne(int i, int j) {
		if ((i < 0) || (i >= taille_l) || (j < 0) || (j >= taille_l))
			throw std::domain_error("echange de lignes : hors domaine");

		if (i == j)
			return;

		for (int k(0); k < taille_c; ++k) {
			T temp = coeffs[i][k];
			coeffs[i][k] = coeffs[j][k];
			coeffs[j][k] = temp;
		}
		return;
	};

	void ajouterLigne(int i, int j, T coefficient) { // ajouter ligne i * coeffs à la ligne j ...
		if ((i < 0) || (i >= taille_l) || (j < 0) || (j >= taille_l))
			throw std::domain_error("ajout de lignes : hors domaine");
		if (i == j)
			throw std::domain_error("ajout de lignes : même ligne. Utiliser la multiplication");
		for (int k(0); k < taille_c; ++k) {
			coeffs[j][k] = coeffs[j][k] + (coefficient * coeffs[i][k]);
		}
		return;
	};

	void multiplierLigne(int i, T coefficient) {
		if ((i < 0) || (i >= taille_l))
			throw std::domain_error("multiplier ligne : hors domaine");

		for (int j(0); j < taille_c; ++j)
			coeffs[i][j] = coefficient * coeffs[i][j];
	};

	matrice<T> inverse() const {
		T vrai_ = unite(coeffs[0][0],true);

		if( taille_l != taille_c)
			throw std::domain_error("inverse de matrice : dimensons ne coincident pas");
		int taille = taille_l;

		matrice<T> m_matrice(*this);
		matrice<T> resultat= unite(m_matrice,true);


		for (int i(0); i < taille; ++i) {
			int j(i);
			for (j = i; j < taille; ++j)
				if ((bool)coeffs[i][j])
					break;
			if (j == taille)
				throw std::domain_error("matrice non-inversible");

			if (!((bool)coeffs[i][j])) {
				return (resultat = false);
			}
			if (i != j) {
				m_matrice.echangerLigne(i, j);
				resultat.echangerLigne(i, j);
			}
			resultat.multiplierLigne(i, vrai_ / m_matrice.coeffs[i][i]);
			m_matrice.multiplierLigne(i, vrai_ / m_matrice.coeffs[i][i]);

			for (j = 0; j < taille; ++j) {
				if (j == i)
					continue;
				resultat.ajouterLigne(i, j, -m_matrice.coeffs[j][i]);
				m_matrice.ajouterLigne(i, j, -m_matrice.coeffs[j][i]);
			}
		}

		return resultat;
	};

	polynome<T> polynomeCaracteristique() const {
		T mvrai = -unite(coeffs[0][0],true);

		if (taille_l != taille_c)
			throw std::domain_error("polynome caracteristique : dimensons ne coincident pas");
		int taille = taille_l;

		matrice<polynome<T>> m_matrice(taille, taille);
		for (int i(0); i < taille; ++i)
			for (int j(0); j < taille; ++j)
				if (i != j)
					m_matrice.coeffs[i][j] = polynome<T>(coeffs[i][j]);
				else
					m_matrice.coeffs[i][j] = polynome<T>(coeffs[i][j], mvrai);

		auto det = m_matrice.determinant();

		return det;
	};
	//si T est un anneau, calculer pour matrice<rationnel<T>> et simplifier les fractions de polynome<rationnel<T>>

	std::vector<T> resoudre(std::vector<T> Y) {
		matrice<T> m_matrice(*this);
		if (m_matrice.taille_l != Y.size())
			throw std::domain_error("resolution equation lineaire : les dimensions ne correspondent pas");
		T vrai_ = unite(coeffs[0][0],true);

		// on parcourt les lignes. Le premier élément non-nul : annule toute la colonne
		for (int i(0); i < taille_l; ++i) {
			int j;
			for (j = 0; j < taille_c; ++j)
				if ((bool)m_matrice[i][j])
					break;
			if (j == taille_c)
				continue;
			T inv = vrai_ / m_matrice.coeffs[i][j]; //on est dans la ligne i, colonne j.
			for (int k(0); k < taille_l; ++k) { //ligne k != i.
				if (k == i)
					continue;
				Y[k] = Y[k] - inv * m_matrice.coeffs[k][j] * Y[i];
				m_matrice.ajouterLigne(i, k, -inv * m_matrice.coeffs[k][j]); //ajoute la ligne i à la ligne k. annnule [k][j]
			}
		}

		//on parcourt les lignes. Si elle est nulle, si Y[k] est non-nul, il y a une erreur
		for (int i(0); i < taille_l; ++i) {
			bool test = false;
			for (int j(0); j < taille_c; ++j)
				if ((bool)m_matrice[i][j]) {
					test = true;
					break;
				}
			if (test)
				continue;
			if ((bool)Y[i]) //ligne nulle, et Y non-nul. Renvoyer vide.
				return std::vector<T>(0);
		}

		std::vector<T> resultat(taille_c);

		//une solution est possible. La calculer. Pour ceci, tout d'abord enlever les "parametres libres". (virtuellement)
		//ensuite, pour chaque premier élément non-nul, = Y[k]
		for (int i(0); i < taille_l; ++i) {
			int j;
			for (j = 0; j < taille_c; ++j)
				if ((bool)m_matrice[i][j])
					break;
			if (j == taille_c) //ligne nulle. 
				continue;
			resultat[j] = Y[i] / m_matrice.coeffs[i][j];
		}

		return resultat;
	};

	sous_ev<T>* noyau() {
		sous_ev<T>* resultat = new sous_ev<T>();
		matrice<T> m_matrice(*this);

		T vrai_ = unite(coeffs[0][0],true);

		// on parcourt les lignes. Le premier élément non-nul : annule toute la colonne
		for (int i(0); i < taille_l; ++i) {
			int j;
			for (j = 0; j < taille_c; ++j)
				if ((bool)m_matrice[i][j])
					break;
			if (j == taille_c)
				continue;
			T inv = vrai_ / m_matrice.coeffs[i][j];
			for (int k(0); k < taille_l; ++k) {
				if (k == i)
					continue;
				m_matrice.ajouterLigne(i, k, -inv * m_matrice.coeffs[k][j]); //ajoute la ligne i à la ligne k. annnule [k][j]
			}
		}

		//on sépare les parametres libres. les colonnes nulles, et les éléments non-nuls apres un premier élément non-nul.
		std::vector<int> parametres_libres(0);
		std::vector<int> parametres_non_libres(0);// on s'en fout !!!
		// on commence : premier élément non-nul, et ceux qui suivent qui sont non-nuls sont des parametres libres.
		for (int i(0); i < taille_l; ++i) {
			int j;
			for (j = 0; j < taille_c; ++j)
				if ((bool)m_matrice[i][j])
					break;
			if (j == taille_c)
				continue;
			parametres_non_libres.push_back(j);
			for (int k(j + 1); k < taille_c; ++k)
				if ((bool)m_matrice[i][k])
					parametres_libres.push_back(k);
		}
		
		//puis, les colonnes nulles : parametres libres.
		for (int i(0); i < taille_c; ++i) {
			bool test = false;
			for(int j(0);j<taille_l;++j)
				if ((bool)m_matrice[j][i]) {
					test = true;
					break;
				}
			if (test)
				continue;
			parametres_libres.push_back(i); //la colonne i est nulle
		}

		parametres_libres.erase(std::unique(parametres_libres.begin(), parametres_libres.end()), parametres_libres.end()); //supprimer les doublons. A vérifier.
		for (int vec_k(0); vec_k < parametres_libres.size(); ++vec_k) { //parametres libre : on en prend un à la fois.
			T temp = unite(m_matrice.coeffs[0][0],false);

			std::vector<T> vec(taille_c, temp);
			temp = unite(temp,true);
			vec[parametres_libres[vec_k]] = temp; //faux partout sauf dans ce parametre libre.

			//on parcourt les premiers-éléments non-nuls. Et on calcule leurs valeurs.
			for (int i(0); i < taille_l; ++i) {
				int j;
				for (j = 0; j < taille_c; ++j)
					if ((bool)m_matrice[i][j])
						break;
				if (j == taille_c)
					continue;

				vec[j] = -m_matrice[i][parametres_libres[vec_k]] / m_matrice[i][j];
			}

			resultat->ajouter_vecteur(vec);
		}

		return resultat;
	}

	friend std::ostream& operator<<(std::ostream& os, const matrice<T>& element) {
		os << "{ ";
		for (int i(0); i < element.taille_l; ++i) {
			os << "{ ";
			for (int j(0); j < element.taille_c - 1; ++j)
				os << element.coeffs[i][j] << " , ";
			if (i < element.taille_l - 1)
				os << element.coeffs[i][element.taille_c - 1] << " } ,";
			else
				os << element.coeffs[i][element.taille_c - 1] << " }";

			if (i < element.taille_l - 1)
				os << std::endl;
		}
		os << "}" << std::endl;
		return os;
	};


	friend std::vector<T> operator*(std::vector<T> const& vec_ligne, matrice<T> const& m_matrice) {
		if (vec_ligne.size() != m_matrice.taille_l)
			throw std::domain_error("multiplication vecteur*matrice : les dimensions ne correspondent pas");

		T faux_ = unite(m_matrice.coeffs[0][0],false);

		std::vector<T> resultat(m_matrice.taille_c);
		for (int i(0); i < m_matrice.taille_c; ++i) {
			T somme = faux_;
			for (int j(0); j < m_matrice.taille_l; ++j)
				somme = somme + (vec_ligne[j] * m_matrice.coeffs[j][i]);
			resultat[i] = somme;
		}

		return resultat;
	};

	friend std::vector<T> operator*(matrice<T> const& m_matrice, std::vector<T> const& vec_colonne) {
		if (vec_colonne.size() != m_matrice.taille_c)
			throw std::domain_error("multiplication matrice*vecteur : les dimensions ne correspondent pas");

		T faux_ = unite(m_matrice.coeffs[0][0],false);

		std::vector<T> resultat(m_matrice.taille_l);
		for (int i(0); i < m_matrice.taille_l; ++i) {
			T somme = faux_;
			for (int j(0); j < m_matrice.taille_c; ++j)
				somme = somme + (m_matrice.coeffs[i][j] * vec_colonne[j]);
			resultat[i] = somme;
		}
		return resultat;
	};

	int taille_l;
	int taille_c;
	std::vector < std::vector< T>> coeffs;
};


template<class T> class matrice<T, typename std::enable_if_t<type_algebre<T>::type == 0>, typename std::enable_if_t<type_algebre<T>::approx == 1>> { //division, approx
public:

	friend matrice<T> derivee(const matrice<T>& element) {
		matrice<T> result(element);
		for (int i(0); i < element.taille_l; ++i)
			for (int j(0); j < element.taille_c; ++j)
				result.coeffs[i][j] = derivee(result.coeffs[i][j]);
		return result;
	};

	explicit matrice(int m_taille_l, int m_taille_c, T element) {
		taille_l = m_taille_l;
		taille_c = m_taille_c;
		coeffs = std::vector<std::vector<T>>(m_taille_l, std::vector<T>(m_taille_c, element));
	};

	explicit matrice(int m_taille_l, int m_taille_c) { //déconseillé à l'utilisation, sauf si vous savez ce que vous faites.*
		taille_l = m_taille_l;
		taille_c = m_taille_c;

		coeffs = std::vector<std::vector<T>>(m_taille_l, std::vector<T>(m_taille_c));
	};

	matrice(matrice<T> const& copie) {
		taille_l = copie.taille_l;
		taille_c = copie.taille_c;
		coeffs = copie.coeffs;
	};

	explicit matrice(const std::vector<std::vector<T>>& vec) : coeffs(vec), taille_l(vec.size()), taille_c(vec[0].size()) {
		for (int i(1); i < taille_l; ++i)
			if (coeffs[i].size() != taille_c)
				throw std::domain_error("creation de matrice : non-rectangulaire");
	};

	explicit matrice() : taille_l(0), taille_c(0), coeffs(0, std::vector<T>(0)) {};

	explicit operator bool() const {
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				if ((bool)coeffs[i][j])
					return true;
		return false;
	};

	template<class U> explicit operator matrice<U>() const {
		matrice<U> result(taille_l,taille_c);

		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				result.coeffs[i][j] = (U)coeffs[i][j];

		return result;
	};


	matrice<T>& operator=(bool test) {
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				coeffs[i][j] = false;
		
		int taille = min(taille_l, taille_c);
		if (test)
			for (int i(0); i < taille; ++i)
				coeffs[i][i] = true;
		return *this;
	};

	template<class U> friend matrice<T> operator*(U scalaire, const matrice<T>& m_matrice) {
		matrice<T> result(m_matrice);
		for (int i(0); i < result.taille_l; ++i)
			for (int j(0); j < result.taille_c; ++j)
				result.coeffs[i][j] = scalaire * result.coeffs[i][j];
		return result;
	};

	matrice<T> operator*(const matrice<T>& autre) const {
		T faux_ = unite(coeffs[0][0],false);

		if (taille_c != autre.taille_l)
			throw std::domain_error("multiplication de matrices : les dimensions ne coïncident pas");
		
		matrice<T> result(taille_l, autre.taille_c);
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < autre.taille_c; ++j) {
				T somme(faux_);
				for (int k(0); k < taille_c; ++k)
					somme = somme + (coeffs[i][k] * autre.coeffs[k][j]);
				result.coeffs[i][j] = somme;
			}
		return result;
	};

	matrice<T> operator+(const matrice<T>& autre) const {
		if ((taille_l != autre.taille_l) || (taille_c != autre.taille_c))
			throw std::domain_error("addition de matrices : les dimensions ne coïncident pas");

		matrice<T> result(taille_l, taille_c);
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				result.coeffs[i][j] = coeffs[i][j] + autre.coeffs[i][j];
		return result;
	};

	matrice<T> operator-(const matrice<T>& autre) const {
		if ((taille_l != autre.taille_l) || (taille_c != autre.taille_c))
			throw std::domain_error("addition de matrices : les dimensions ne coïncident pas");

		matrice<T> result(taille_l, taille_c);
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				result.coeffs[i][j] = coeffs[i][j] - autre.coeffs[i][j];
		return result;
	};

	friend 	matrice<T> operator-(const matrice<T>& element) {
		matrice<T> result(element.taille_l, element.taille_c);
		for (int i(0); i < result.taille_l; ++i)
			for (int j(0); j < result.taille_c; ++j)
				result.coeffs[i][j] = -element.coeffs[i][j];
		return result;
	};

	T determinant_anneau() {
		T resultat = unite(coeffs[0][0],false);
		T vrai_ = unite(resultat,true);

		if (taille_l != taille_c)
			throw std::domain_error("determinant de matrice : les dimensions ne coïncident pas");
		int taille = taille_l;

		for (fact_for iter(taille); (bool)iter; ++iter) {
			T temp = vrai_;
			for (int i(0); i < taille; ++i)
				temp = temp * coeffs[i][iter.permutation[i]];
			if (iter.signature < 0)
				temp = -temp;
			resultat = resultat + temp;
		}

		return resultat;
	};

	T determinant() const {
		matrice<T> m_matrice(*this);
		T det = unite(coeffs[0][0],true);
		T vrai_ = det;

		if (taille_l != taille_c)
			throw std::domain_error("determinant de matrice : les dimensions ne coïncident pas");
		int taille = taille_l;


		int taille_max = max(taille - 4, 0); //pour les 4 derniers : on utilise le determinant_anneau
		int taille_fin = taille - taille_max;

		for (int i(0); i < taille_max; ++i) {
			int j_max = -1;
			auto norme_max = norme(vrai_);
			for (int j = i; j < taille; ++j) {
				if ((bool) m_matrice.coeffs[i][j]) {
					auto norme_temp = norme(m_matrice.coeffs[i][j]);
					if (j_max == -1) {
						norme_max = norme_temp;
						j_max = j;
					}
					else
						if (norme_temp > norme_max) {
							norme_max = norme_temp;
							j_max = j;
						}
				}
			}
			if (j_max==-1) //0 sur une colonne
				return (det = false);

			if (i != j_max) {
				m_matrice.echangerLigne(i, j_max);
				det = -det;
			}
			for (int j = i + 1; j < taille; ++j) {
				m_matrice.ajouterLigne(i, j, ((-m_matrice.coeffs[j][i]) / m_matrice.coeffs[i][i]));
			}

		}

		for (int i(0); i < taille_max; ++i)
			det = (det * m_matrice.coeffs[i][i]);

		if (taille_fin > 0) {
			matrice<T> m_matrice_fin(taille_fin, taille_fin);
			for(int i(0);i<taille_fin;++i)
				for (int j(0); j < taille_fin; ++j) 
					m_matrice_fin.coeffs[i][j] = m_matrice.coeffs[taille_max + i][taille_max + j];
			det = det * m_matrice_fin.determinant_anneau();
		}

		return det;
	};


	void echangerLigne(int i, int j) {
		if ((i < 0) || (i >= taille_l) || (j < 0) || (j >= taille_l))
			throw std::domain_error("echange de lignes : hors domaine");

		if (i == j)
			return;
		for (int k(0); k < taille_c; ++k) {
			T temp = coeffs[i][k];
			coeffs[i][k] = coeffs[j][k];
			coeffs[j][k] = temp;
		}
		return;
	};

	void ajouterLigne(int i, int j, T coefficient) { // ajouter ligne i * coeffs à la ligne j ...
		if ((i < 0) || (i >= taille_l) || (j < 0) || (j >= taille_l))
			throw std::domain_error("ajout de lignes : hors domaine");
		if (i == j)
			throw std::domain_error("ajout de lignes : même ligne. Utiliser la multiplication.");

		for (int k(0); k < taille_c; ++k) {
			coeffs[j][k] = coeffs[j][k] + (coefficient * coeffs[i][k]);
		}
		return;
	};

	void multiplierLigne(int i, T coefficient) {
		for (int j(0); j < taille_c; ++j)
			coeffs[i][j] = coefficient * coeffs[i][j];
	};

	matrice<T> inverse() const {
		T vrai_ = unite(coeffs[0][0],true);

		if(taille_l != taille_c)
			throw std::domain_error("inverse de matrice : les dimensions ne coïncident pas.");
		int taille = taille_l;

		matrice<T> m_matrice(*this);
		matrice<T> resultat = unite(*this,true); //identité

		for (int i(0); i < taille; ++i) {
			int j_max = -1;
			auto norme_max = norme(vrai_);
//			auto norme_max = norme_T<T>::norme(vrai);
			for (int j = i; j < taille; ++j) {
				if ((bool)m_matrice.coeffs[i][j]) {
					auto norme_temp = norme(m_matrice.coeffs[i][j]);
					if (j_max == -1) {
						norme_max = norme_temp;
						j_max = j;
					}
					else
						if (norme_temp > norme_max) {
							norme_max = norme_temp;
							j_max = j;
						}
				}
			}
			if (j_max == -1)
				throw std::domain_error("matrice non-inversible");

			if (i != j_max) {
				m_matrice.echangerLigne(i, j_max);
				resultat.echangerLigne(i, j_max);
			}
			resultat.multiplierLigne(i, vrai_ / m_matrice.coeffs[i][i]);
			m_matrice.multiplierLigne(i, vrai_ / m_matrice.coeffs[i][i]);

			for (int j = 0; j < taille; ++j) {
				if (j == i)
					continue;
				resultat.ajouterLigne(i, j, -m_matrice.coeffs[j][i]);
				m_matrice.ajouterLigne(i, j, -m_matrice.coeffs[j][i]);
			}
		}

		return resultat;
	};

	polynome<T> polynomeCaracteristique() const {
		T mvrai = -unite(coeffs[0][0],true);

		if (taille_l != taille_c)
			throw std::domain_error("inverse de matrice : les dimensions ne coïncident pas.");
		int taille = taille_l;


		matrice<polynome<T>> m_matrice(taille,taille);

		for (int i(0); i < taille; ++i)
			for (int j(0); j < taille; ++j)
				if (i != j)
					m_matrice.coeffs[i][j] = polynome<T>(coeffs[i][j]);
				else
					m_matrice.coeffs[i][j] = polynome<T>(coeffs[i][j], mvrai);

		auto det = m_matrice.determinant();

		return det;
	};

	std::vector<T> resoudre(std::vector<T> Y) { //dans le cas approx.
		matrice<T> m_matrice(*this);
		if (m_matrice.taille_l != Y.size())
			throw std::domain_error("resolution equation lineaire : les dimensions ne correspondent pas");
		T vrai_ = unite(coeffs[0][0],true);

		std::vector<bool> ligne_faite(taille_l,true);

		// on parcourt les lignes. Le premier élément non-nul : annule toute la colonne
		for (int i(0); i < ligne_faite.size(); ++i) {
			if (!ligne_faite[i])
				continue;
			int j;
			for (j = 0; j < taille_c; ++j)
				if ((bool)m_matrice.coeffs[i][j])
					break;
			if (j == taille_c)
				continue;
			//le premier non-nul est la colonne j

			//on regarde la liste des lignes, qui ont la même colonne comme premier élément non-nul. Element de lignes_restantes.
			int i_max = i;
			auto norme_max = norme(m_matrice.coeffs[i][j]);

			for (int k(0); k < taille_l; ++k) {
				if (!ligne_faite[k])
					continue;
				if (k == i)
					continue;
				int l;
				for (l=0; l <= j; ++l)
					if ((bool)m_matrice.coeffs[k][l])
						break;
				if (l == j) {
					if (norme(m_matrice.coeffs[k][l]) > norme_max) {
						i_max = k;
						norme_max = norme(m_matrice.coeffs[k][l]);
					}
				}
			}
			
			//on a trouve notre ligne. c'est i_max. On enleve i_max de ligne_faite.
			ligne_faite[i_max] = false;

			T inv = vrai_ / m_matrice.coeffs[i_max][j];
			for (int k(0); k < taille_l; ++k) {
				if (k == i_max)
					continue;
				Y[k] = Y[k] - inv * m_matrice.coeffs[k][j] * Y[i_max]; // VERIFIER
				m_matrice.ajouterLigne(i_max, k, -inv * m_matrice.coeffs[k][j]); //ajoute la ligne i à la ligne k. annnule [k][j]
			}
			i = 0; //on reparcourt toutes les lignes non-faites. Ca va vite (test booléen)
		}

		//on parcourt les lignes. Si elle est nulle, si Y[k] est non-nul, il y a une erreur
		for (int i(0); i < taille_l; ++i) {
			bool test = false;
			for (int j(0); j < taille_c; ++j)
				if ((bool)m_matrice.coeffs[i][j]) {
					test = true;
					break;
				}
			if (test)
				continue;
			if ((bool)Y[i])
				return std::vector<T>(0); //alors on retourne le vide
		}

		T faux_ = unite(vrai_,false);
		std::vector<T> resultat(taille_c, faux_);

		//une solution est possible. La calculer. Pour ceci, tout d'abord enlever les "parametres libres". (virtuellement)
		//ensuite, pour chaque premier élément non-nul, = Y[k] / élément de matrice
		for (int i(0); i < taille_l; ++i) {
			int j;
			for (j = 0; j < taille_c; ++j)
				if ((bool)m_matrice.coeffs[i][j])
					break;
			if (j == taille_c)
				continue;
			resultat[j] = Y[i] / m_matrice.coeffs[i][j];
		}

		return resultat;
	}

	sous_ev<T>* noyau() {
		sous_ev<T>* resultat = new sous_ev<T>();
		matrice<T> m_matrice(*this);

		T vrai_ = unite(coeffs[0][0],true);

		std::vector<bool> ligne_faite(taille_l, true);

		// on parcourt les lignes. Le premier élément non-nul : annule toute la colonne
		for (int i(0); i < ligne_faite.size(); ++i) {
			if (!ligne_faite[i])
				continue;
			int j;
			for (j = 0; j < taille_c; ++j)
				if ((bool) m_matrice.coeffs[i][j])
					break;
			if (j == taille_c)
				continue;

			//on regarde la liste des lignes, qui ont la même colonne comme premier élément non-nul. Element de ligne_faite.
			int i_max = i;
			auto norme_max = norme(m_matrice.coeffs[i][j]);

			for (int k(0); k < taille_l; ++k) {
				if (!ligne_faite[k])
					continue;
				if (k == i)
					continue;
				int l;
				for (l=0; l < j + 1; ++l)
					if ((bool)m_matrice.coeffs[k][l])
						break;
				if (l == j) {
					if (norme(m_matrice.coeffs[k][l]) > norme_max) {
						i_max = k;
						norme_max = norme(m_matrice.coeffs[k][l]);
					}
				}
			}

			//on a trouve notre ligne. On enleve i_max de ligne_faite.
			ligne_faite[i_max] = false;

			T inv = vrai_ / m_matrice.coeffs[i_max][j];
			for (int k(0); k < taille_l; ++k) {
				if (k == i_max)
					continue;
				m_matrice.ajouterLigne(i_max, k, -inv * m_matrice.coeffs[k][j]); //ajoute la ligne i à la ligne k. annnule [k][j]
			}
			i = 0;
		}
		
		//on sépare les parametres libres. les colonnes nulles, et les éléments non-nuls apres un premier élément non-nul.
		std::vector<int> parametres_libres(0);
		std::vector<int> parametres_non_libres(0); //on s'en fiche en fait.
		// on commence : premier élément non-nul, et ceux qui suivent qui sont non-nuls
		for (int i(0); i < taille_l; ++i) {
			int j;
			for (j = 0; j < taille_c; ++j)
				if ((bool)m_matrice.coeffs[i][j])
					break;
			if (j == taille_c)
				continue;
			parametres_non_libres.push_back(j);
			for (int k(j + 1); k < taille_c; ++k)
				if ((bool)m_matrice.coeffs[i][k])
					parametres_libres.push_back(k);
		}

		//puis, les colonnes nulles : parametres libres.
		for (int i(0); i < taille_c; ++i) {
			bool test = false;
			for (int j(0); j < taille_l; ++j)
				if ((bool)m_matrice.coeffs[j][i]) {
					test = true;
					break;
				}
			if (test)
				continue;
			parametres_libres.push_back(i); //la colonne i est nulle
		}

		parametres_libres.erase(std::unique(parametres_libres.begin(), parametres_libres.end()), parametres_libres.end()); //supprimer les doublons. A vérifier.
		for (int vec_k(0); vec_k < parametres_libres.size(); ++vec_k) { //parametres libre : on en prend un à la fois.
			T temp = unite(m_matrice.coeffs[0][0],false);
			std::vector<T> vec(taille_c, temp);
			temp = unite(temp,true);
			vec[parametres_libres[vec_k]] = temp; //faux partout sauf dans ce parametre libre.

			//on parcourt les premiers-éléments non-nuls. Et on calcule leurs valeurs.
			for (int i(0); i < taille_l; ++i) {
				int j;
				for (j = 0; j < taille_c; ++j)
					if ((bool)m_matrice.coeffs[i][j])
						break;
				if (j == taille_c)
					continue;

				vec[j] = -m_matrice.coeffs[i][parametres_libres[vec_k]] / m_matrice.coeffs[i][j];
			}

			resultat->ajouter_vecteur(vec);
		}

		return resultat;
	}


	friend std::ostream& operator<<(std::ostream& os, const matrice<T>& element) {
		os << "{ ";
		for (int i(0); i < element.taille_l; ++i) {
			os << "{ ";
			for (int j(0); j < element.taille_c - 1; ++j)
				os << element.coeffs[i][j] << " , ";
			if (i < element.taille_l - 1)
				os << element.coeffs[i][element.taille_c - 1] << " } ,";
			else
				os << element.coeffs[i][element.taille_c - 1] << " }";

			if (i < element.taille_l - 1)
				os << std::endl;
		}
		os << "}" << std::endl;
		return os;
	};


	friend std::vector<T> operator*(std::vector<T> const& vec_ligne, matrice<T> const& m_matrice) {
		if (vec_ligne.size() != m_matrice.taille_l)
			throw std::domain_error("multiplication vecteur*matrice : les dimensions ne correspondent pas");

		T faux_ = unite(m_matrice.coeffs[0][0],false);

		std::vector<T> resultat(m_matrice.taille_c);
		for (int i(0); i < m_matrice.taille_c; ++i) {
			T somme = faux_;
			for (int j(0); j < m_matrice.taille_l; ++j)
				somme = somme + (vec_ligne[j] * m_matrice.coeffs[j][i]);
			resultat[i] = somme;
		}

		return resultat;
	};

	friend std::vector<T> operator*(matrice<T> const& m_matrice, std::vector<T> const& vec_colonne) {
		if (vec_colonne.size() != m_matrice.taille_c)
			throw std::domain_error("multiplication matrice*vecteur : les dimensions ne correspondent pas");

		T faux_ = unite(m_matrice.coeffs[0][0],false);

		std::vector<T> resultat(m_matrice.taille_l);
		for (int i(0); i < m_matrice.taille_l; ++i) {
			T somme = faux_;
			for (int j(0); j < m_matrice.taille_c; ++j)
				somme = somme + (m_matrice.coeffs[i][j] * vec_colonne[j]);
			resultat[i] = somme;
		}
		return resultat;
	};

	int taille_l;//nombre de lignes
	int taille_c; //nombre de colonnes
	std::vector < std::vector< T>> coeffs;
};

template<class T> class sous_ev {
public:
	sous_ev() : liste_vecteurs(0){
	}

	void ajouter_vecteur(std::vector<T> vecteur) {
		liste_vecteurs.push_back(vecteur);
	}

	std::vector<std::vector<T>> liste_vecteurs;
};

//template matrice<double>;

/*
template<class T> class matrice<T, typename std::enable_if<type_algebre<T>::type == 1>::type> : public matrice_base<T> {
public:
	using matrice_base<T>::coeffs;
	using matrice_base<T>::taille;

	explicit matrice() : taille(0), coeffs(0, std::vector<T>(0)) {
		std::cout << "copie";
	};

	explicit matrice(int m_taille, T element) : matrice<T>::coeffs(m_taille, std::vector<T>(m_taille, element)), matrice<T>::taille(m_taille) {
		std::cout << "copie";
	};

	virtual matrice<T>& operator=(matrice<T > const& copie) {
		taille = copie.taille;
		std::vector<std::vector<T>> coeffs = copie.coeffs;
		std::cout << "copie";
		return *this;
	};

	T determinant() const {
		matrice<rationnel<T>> m_matrice(taille);
		T vrai = coeffs[0][0];
		vrai = true;

		for (int i(0); i < taille; ++i)
			for (int j(0); j < taille; ++j)
				m_matrice.coeffs[i][j] = rationnel<T>(coeffs[i][j], vrai);

		rationnel<T> det = m_matrice.determinant();
		std::cout << "calcul determinant" << std::endl;

		return det.numerateur / det.denominateur;
	};

};

template<class T> class matrice<T, typename std::enable_if<type_algebre<T>::type == 0>::type> : public matrice_base<T> {
public:
	using matrice_base<T>::coeffs;
	using matrice_base<T>::taille;

};
*/

//template<class T> class matrice<T, typename std::enable_if<type_algebre<T>::type == 0, float>::type> : virtual public matrice<T, void> {
