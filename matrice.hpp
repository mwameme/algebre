#pragma once
#include <vector>
#include <iostream>
#include <exception>
#include <initializer_list>

#include "entete objets.hpp"
#include "fact_for.hpp"
#include "simplifier polynome_n.hpp"

#include <algorithm>

//#include <typeinfo> 

template<class T> class sous_ev;

template<class T, class enable1 = void, class enable2 = void> class matrice;


// matrice:
//cas type=1 ou 2 (static_assert)
//cas type=0, approx =1 ou 0
// 
//rationnel : cas type=1 ou 2
// 
//diagonalisation : cas approx=0 ou 1. Static_assert type=0
//

template<class T> class matrice<T,void,void> { //sans division : determinant permutation
public:

	friend matrice<T> derivee(const matrice<T>& element) {
		matrice<T> result(element);
		for (int i(0); i < element.taille_l; ++i)
			for (int j(0); j < element.taille_c; ++j)
				result.coeffs[i][j] = derivee(result.coeffs[i][j]);
		return result;
	};

	explicit matrice(int m_taille_l, int m_taille_c, T element) : taille_l(m_taille_l), taille_c(m_taille_c), coeffs(m_taille_l, std::vector<T>(m_taille_c, element)) { };

	explicit matrice(int m_taille_l,int m_taille_c) : taille_l(m_taille_l), taille_c(m_taille_c), coeffs(m_taille_l, std::vector<T>(m_taille_c)) {	};

	matrice(matrice<T> const& copie) :taille_l(copie.taille_l),taille_c(copie.taille_c),coeffs(copie.coeffs) {	};

	matrice(matrice<T>&& copie) {
		taille_l = copie.taille_l;
		taille_c = copie.taille_c;
		std::swap(coeffs, copie.coeffs);
		return;
	}

	explicit matrice(const std::vector<std::vector<T>>& vec) :  taille_l(vec.size()),taille_c(vec[0].size()) , coeffs(vec){
		for(int i(1);i<taille_l;++i)
			if (coeffs[i].size() != taille_c)
				throw std::domain_error("creation de matrice : non-rectangulaire");
	};

	explicit matrice(std::vector<std::vector<T>>&& vec) : taille_l(vec.size()), taille_c(vec[0].size()) {
		std::swap(coeffs, vec);
		for (int i(1); i < taille_l; ++i)
			if (coeffs[i].size() != taille_c)
				throw std::domain_error("creation de matrice : non-rectangulaire");
	};


	explicit matrice() {}; //A NE PAS UTILISER


	template<class U> explicit matrice(std::vector< std::vector<U>> liste) {
		taille_l = liste.size();
		if (taille_l == 0)
			throw std::domain_error("initialisation de matrice : liste vide");
		taille_c = liste[0].size();
		if (taille_c == 0)
			throw std::domain_error("initialisation de matrice : nombre de colonnes nul");

		for (int i(1); i < taille_l; ++i)
			if (liste[i].size() != taille_c)
				throw std::domain_error("initialisation de matrice : non-rectangulaire");

		coeffs = std::vector<std::vector<T>>(taille_l, std::vector<T>(taille_c));
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				coeffs[i][j] = T(liste[i][j]);

		return;
	};

	matrice<T>& operator=(const matrice<T>& temp);

	matrice<T>& operator=(matrice<T>&& temp) {
		swap(*this, temp);
		return *this;
	};

	matrice<T>& operator*=(matrice<T> const& temp) {
		return (*this = (*this * temp));
	};

	template<class U>
	matrice<T>& operator*=(U const& scalaire) {
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				coeffs[i][j] *= scalaire;
		return *this;
	};

	explicit inline operator bool() const {
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


	template<class U> friend matrice<T> operator*(U scalaire, const matrice<T>& m_matrice) {
		matrice<T> result(m_matrice);
		for (int i(0); i < result.taille_l; ++i)
			for (int j(0); j < result.taille_c; ++j)
				result.coeffs[i][j] = scalaire * result.coeffs[i][j];
		return result;
	};

	friend matrice<T> operator*(const matrice<T>& gauche,const matrice<T>& droit) {
		if (gauche.taille_c != droit.taille_l)
			throw std::domain_error("multiplication de matrices : les dimensions ne co�ncident pas");

		matrice<T> result(gauche.taille_l,droit.taille_c);
		T faux = unite(gauche.coeffs[0][0], false);
		for (int i(0); i < gauche.taille_l; ++i)
			for (int j(0); j < droit.taille_c; ++j) {
				T somme = faux;
				for (int k(0); k < gauche.taille_c; ++k)
					somme = somme + (gauche.coeffs[i][k] * droit.coeffs[k][j]);
				result.coeffs[i][j] = somme;
			}
		return result;
	};

	friend matrice<T> operator+(const matrice<T>& gauche, const matrice<T>& droit) {
		if ((gauche.taille_l != droit.taille_l) || (gauche.taille_c != droit.taille_c))
			throw std::domain_error("addition de matrices : les dimensions ne co�ncident pas");

		matrice<T> result(gauche.taille_l,gauche.taille_c);
		for (int i(0); i < gauche.taille_l; ++i)
			for (int j(0); j < gauche.taille_c; ++j)
				result.coeffs[i][j] = (gauche.coeffs[i][j] + droit.coeffs[i][j]);
		return result;
	};

	friend matrice<T> operator-(const matrice<T>& gauche, const matrice<T>& droit) {
		if ((gauche.taille_l != droit.taille_l) || (gauche.taille_c != droit.taille_c))
			throw std::domain_error("addition de matrices : les dimensions ne co�ncident pas");

		matrice<T> result(gauche.taille_l, gauche.taille_c);
		for (int i(0); i < gauche.taille_l; ++i)
			for (int j(0); j < gauche.taille_c; ++j)
				result.coeffs[i][j] = (gauche.coeffs[i][j] - droit.coeffs[i][j]);
		return result;
	};

	friend matrice<T> operator-(matrice<T> const& element) {
		matrice<T> result(element.taille_l,element.taille_c);
		for (int i(0); i < element.taille_l; ++i)
			for (int j(0); j < element.taille_c; ++j)
				result.coeffs[i][j] = - element.coeffs[i][j];
		return result;
	};

	friend matrice<T> operator/(matrice<T> const& gauche, matrice<T> const& droit) {
		return gauche * droit.inverse();
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

	void ajouterLigne(int i, int j, T coefficient) { // ajouter ligne i * coeffs � la ligne j ...
		if ((i < 0) || (i >= taille_l) || (j < 0) || (j >= taille_l))
			throw std::domain_error("ajout de lignes : hors domaine");
		if (i == j)
			throw std::domain_error("ajout de lignes : m�me ligne");
		for (int k(0); k < taille_c; ++k) {
			coeffs[j][k] = coeffs[j][k] + (coefficient * coeffs[i][k]);
		}
		return;
	};

	void multiplierLigne(int i, T coefficient) {
		for (int j(0); j < taille_c; ++j)
			coeffs[i][j] = coefficient * coeffs[i][j];
	};



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

	friend void swap(matrice<T>& gauche, matrice<T>& droit) {
		std::swap(gauche.taille_l, droit.taille_l);
		std::swap(gauche.taille_c, droit.taille_c);
		std::swap(gauche.coeffs, droit.coeffs);
		return;
	};

	friend bool operator==(matrice<T> const& gauche, matrice<T> const& droit) {
		if (gauche.taille_l != droit.taille_l)
			return false;
		if (gauche.taille_c != droit.taille_c)
			return false;
		for (int i(0); i < gauche.taille_l; ++i)
			for (int j(0); j < gauche.taille_c; ++j)
				if (!(gauche.coeffs[i][j] == droit.coeffs[i][j]))
					return false;
		return true;
	}


	T determinant_anneau() const {
		T resultat = unite(coeffs[0][0], false);
		T vrai_ = unite(resultat, true);

		if (taille_l != taille_c)
			throw std::domain_error("determinant : les dimensions ne co�ncident pas");


		for (fact_for iter(taille_l); (bool)iter; ++iter) {
			T temp = vrai_;
			for (int i(0); i < taille_l; ++i)
				temp = temp * coeffs[i][iter.permutation[i]];
			if (iter.signature == -1)
				temp = -temp;
			resultat = resultat + temp;
		}

		return resultat;
	};

	T determinant() const {
		if constexpr (type_algebre<T>::type == 2)
			return determinant_anneau();
		else if constexpr (type_algebre<T>::type == 1) {
			matrice<rationnel<T>> m_matrice(taille_l, taille_l);
			T vrai = unite(coeffs[0][0], true);

			if (taille_l != taille_c)
				throw std::domain_error("determinant matrice : les dimensions ne co�ncident pas");
			int taille = taille_l;


			for (int i(0); i < taille; ++i)
				for (int j(0); j < taille; ++j)
					m_matrice.coeffs[i][j] = rationnel<T>(coeffs[i][j], vrai); //passer par les rationnels. vrai_ optionnel

			rationnel<T> det = m_matrice.determinant();

			T result = det.numerateur / det.denominateur; //simplifier : division avec reste 
			Simplifier_frac_poly(result);
			return result;
		}
		else {
			if constexpr (type_algebre<T>::approx == 0) {
				matrice<T> m_matrice(*this);
				T det = unite(coeffs[0][0], true);
				T vrai = det;

				if (taille_l != taille_c)
					throw std::domain_error("determinant de matrice : les dimensions ne coincident pas");
				int taille = taille_l;

				for (int i(0); i < taille; ++i) {//colonne i
					int j;
					for (j = i; j < taille; ++j) {
						if ((bool)m_matrice.coeffs[j][i])
							break;
					}
					if (j == taille)
						return unite(det, false);

					if (i != j) {
						m_matrice.echangerLigne(i, j);
						det = -det;
					}

					T inv = vrai / m_matrice.coeffs[i][i];
					for (j = i + 1; j < taille; ++j)
						m_matrice.ajouterLigne(i, j, ((-m_matrice.coeffs[j][i]) * inv));  // ajouter ligne i * coeffs � la ligne j ...

				}

				for (int i(0); i < taille; ++i)
					det = (det * m_matrice.coeffs[i][i]);

				Simplifier_frac_poly(det);
				return det;
			}
			else {
				matrice<T> m_matrice(*this);
				T det = unite(coeffs[0][0], true);
				T vrai = det;

				if (taille_l != taille_c)
					throw std::domain_error("determinant de matrice : les dimensions ne co�ncident pas");
				int taille = taille_l;


				int taille_max = max(taille - 4, 0); //pour les 4 derniers : on utilise le determinant_anneau
				int taille_fin = taille - taille_max;

				for (int i(0); i < taille_max; ++i) {
					int j_max = -1;
					auto norme_max = norme(vrai);
					for (int j = i; j < taille; ++j) {
						if ((bool)m_matrice.coeffs[i][j]) {
							auto norme_temp = norme(m_matrice.coeffs[j][i]);
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
					if (j_max == -1) //0 sur une colonne
						return unite(det, false);

					if (i != j_max) {
						m_matrice.echangerLigne(i, j_max);
						det = - det;
					}
					T inv = vrai / m_matrice.coeffs[i][i];
					for (int j = i + 1; j < taille; ++j) {
						m_matrice.ajouterLigne(i, j, ((-m_matrice.coeffs[j][i]) * inv)); // ajouter ligne i * coeffs � la ligne j ...
					}

				}

				for (int i(0); i < taille_max; ++i)
					det = (det * m_matrice.coeffs[i][i]);

				if (taille_fin > 0) {
					matrice<T> m_matrice_fin(taille_fin, taille_fin);
					for (int i(0); i < taille_fin; ++i)
						for (int j(0); j < taille_fin; ++j)
							m_matrice_fin.coeffs[i][j] = m_matrice.coeffs[taille_max + i][taille_max + j];
					det = det * m_matrice_fin.determinant_anneau();
				}
				Simplifier_frac_poly(det);
				return det;
			}
		}

	};

	// type_algebre (2)
	polynome<T> polynomeCaracteristique() const {
		T mvrai = - unite(coeffs[0][0], true);

		if (taille_l != taille_c)
			throw std::domain_error("polynome caracteristique : les dimensions ne correspondent pas");
		int taille = taille_l;

		//		polynome<T> m_poly = polynome<T>(vrai);

		matrice<polynome<T>> m_matrice(taille, taille);

		for (int i(0); i < taille; ++i)
			for (int j(0); j < taille; ++j)
				if (i != j)
					m_matrice.coeffs[i][j] = polynome<T>(coeffs[i][j]);
				else
					m_matrice.coeffs[i][j] = polynome<T>(coeffs[i][j], mvrai);

		polynome<T> det = m_matrice.determinant();
		simplifier_poly(det);
		return det;
	};

	matrice<T> inverse() const {
		static_assert(type_algebre<T>::type == 0);
		if constexpr (type_algebre<T>::approx == 0) {
			T vrai = unite(coeffs[0][0], true);

			if (taille_l != taille_c)
				throw std::domain_error("inverse de matrice : dimensions ne coincident pas");
			int taille = taille_l;

			matrice<T> m_matrice(*this);
			matrice<T> resultat = unite(m_matrice, true);


			for (int i(0); i < taille; ++i) {
				int j;
				for (j = i; j < taille; ++j)
					if ((bool)coeffs[j][i])
						break;
				if (j == taille)
					throw std::domain_error("matrice non-inversible");

				if (i != j) {
					m_matrice.echangerLigne(i, j);
					resultat.echangerLigne(i, j);
				}
				T inv = vrai / m_matrice.coeffs[i][i];
				resultat.multiplierLigne(i, inv);
				m_matrice.multiplierLigne(i, inv);

				for (j = 0; j < taille; ++j) {
					if (j == i)
						continue;
					resultat.ajouterLigne(i, j, -m_matrice.coeffs[j][i]);  // ajouter ligne i * coeffs � la ligne j ...
					m_matrice.ajouterLigne(i, j, -m_matrice.coeffs[j][i]);
				}
			}

			return resultat;
		}
		else {
			T vrai = unite(coeffs[0][0], true);

			if (taille_l != taille_c)
				throw std::domain_error("inverse de matrice : les dimensions ne co�ncident pas.");
			int taille = taille_l;

			matrice<T> m_matrice(*this);
			matrice<T> resultat = unite(*this, true); //identit�

			for (int i(0); i < taille; ++i) {
				int j_max = -1;
				auto norme_max = norme(vrai);
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

				T inv = vrai / m_matrice.coeffs[i][i];
				resultat.multiplierLigne(i, inv);
				m_matrice.multiplierLigne(i, inv);

				for (int j = 0; j < taille; ++j) {
					if (j == i)
						continue;
					resultat.ajouterLigne(i, j, -m_matrice.coeffs[j][i]);
					m_matrice.ajouterLigne(i, j, -m_matrice.coeffs[j][i]);
				}
			}

			return resultat;
		}
	};

	std::vector<T> resoudre(std::vector<T> Y) const { //renvoit vecteur vide si erreur ...
		static_assert(type_algebre<T>::type == 0);
		if constexpr (type_algebre<T>::approx == 0) {
			matrice<T> m_matrice(*this);
			if (m_matrice.taille_l != Y.size())
				throw std::domain_error("resolution equation lineaire : les dimensions ne correspondent pas");
			T vrai = unite(m_matrice.coeffs[0][0], true);

			// on parcourt les lignes. Le premier �l�ment non-nul : annule toute la colonne
			for (int i(0); i < taille_l; ++i) {
				int j;
				for (j = 0; j < taille_c; ++j)
					if ((bool)m_matrice.coeffs[i][j])
						break;
				if (j == taille_c)
					continue;

				T inv = vrai / m_matrice.coeffs[i][j]; //on est dans la ligne i, colonne j. On annule la colonne j (sauf ligne i)
				for (int k(0); k < taille_l; ++k) { //ligne k != i.
					if (k == i)
						continue;
					if ((bool)m_matrice.coeffs[k][j]) {
						Y[k] = Y[k] - inv * m_matrice.coeffs[k][j] * Y[i];
						m_matrice.ajouterLigne(i, k, -inv * m_matrice.coeffs[k][j]); //ajoute la ligne i � la ligne k. annnule [k][j]
					}
				}
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
				if ((bool)Y[i]) //ligne nulle, et Y non-nul. Renvoyer vide.
					return std::vector<T>(0);
			}

			T faux = unite(vrai, false);
			std::vector<T> resultat(taille_c, faux);

			//une solution est possible. La calculer. Pour ceci, tout d'abord enlever les "parametres libres". (virtuellement)
			//ensuite, pour chaque premier �l�ment non-nul, = Y[k]
			for (int i(0); i < taille_l; ++i) {
				int j;
				for (j = 0; j < taille_c; ++j)
					if ((bool)m_matrice.coeffs[i][j])
						break;
				if (j == taille_c) //ligne nulle. 
					continue;
				resultat[j] = Y[i] / m_matrice.coeffs[i][j];
			}

			return resultat;
		}
		else {
			matrice<T> m_matrice(*this);
			if (m_matrice.taille_l != Y.size())
				throw std::domain_error("resolution equation lineaire : les dimensions ne correspondent pas");
			T vrai = unite(coeffs[0][0], true);
			T faux = unite(vrai, false);
			std::vector<bool> ligne_faite(taille_l, true);
			std::vector<int> prem_colonne(taille_l, 0);
			int n_restant = taille_l;

			for (int i(0); i < taille_l; ++i) {
				int j = 0;
				for (; j < taille_c; ++j)
					if ((bool)m_matrice.coeffs[i][j])
						break;
				prem_colonne[i] = j;
				if (j == taille_c) {
					ligne_faite[i] = false;
					--n_restant;
				}
			}

			while (n_restant > 0) {
				std::vector<decltype(norme(vrai))> normes(taille_c, norme(faux));
				std::vector<int> lignes(taille_c, -1);
				std::vector<int> colonnes(0);
				for (int i(0); i < taille_l; ++i) {
					if (!ligne_faite[i])
						continue;
					colonnes.push_back(prem_colonne[i]);
					if (lignes[prem_colonne[i]] == -1) {
						lignes[prem_colonne[i]] = i;
						normes[prem_colonne[i]] = norme(m_matrice.coeffs[i][prem_colonne[i]]);
					}
					else {
						auto norme_temp = norme(m_matrice.coeffs[i][prem_colonne[i]]);
						if (norme_temp > normes[prem_colonne[i]]) {
							normes[prem_colonne[i]] = norme_temp;
							lignes[prem_colonne[i]] = i;
						}
					}
				}

				std::sort(colonnes.begin(), colonnes.end());
				auto last = std::unique(colonnes.begin(), colonnes.end());
				colonnes.erase(last, colonnes.end());

				std::vector<decltype(norme(vrai))> normes_colonnes(taille_c, norme(faux));
				for (int k(0); k < colonnes.size(); ++k) {
					for (int i(0); i < taille_l; ++i) {
						if (i == lignes[colonnes[k]])
							continue;
						auto norme_temp = norme(m_matrice.coeffs[i][colonnes[k]]);
						if (norme_temp > normes_colonnes[colonnes[k]])
							normes_colonnes[colonnes[k]] = norme_temp;
					}
				}


				auto norme_min = norme(faux);
				int j_min = -1;
				for (int j(0); j < taille_c; ++j) {
					if (lignes[j] == -1) //VERIFIER
						continue;
					if (j_min == -1) {
						j_min = j;
						norme_min = normes_colonnes[j] / normes[j];
					}
					else {
						auto norme_temp = normes_colonnes[j] / normes[j];
						if (norme_temp < norme_min) {
							j_min = j;
							norme_min = norme_temp;
						}
					}
				}
				int i_min = lignes[j_min];

				T inv = vrai / m_matrice.coeffs[i_min][j_min];
				for (int k(0); k < taille_l; ++k) {
					if (k == i_min)
						continue;
					if (!(bool)m_matrice.coeffs[k][j_min])
						continue;
					Y[k] = Y[k] - inv * m_matrice.coeffs[k][j_min] * Y[i_min]; // VERIFIER
					m_matrice.ajouterLigne(i_min, k, -inv * m_matrice.coeffs[k][j_min]); //ajoute la ligne i � la ligne k. annnule [k][j_min]

					//ligne_faite et prem_colonne (taille_l)
					if (j_min == prem_colonne[k]) {
						int j = j_min;
						for (; j < taille_c; ++j)
							if ((bool)m_matrice.coeffs[k][j])
								break;
						prem_colonne[k] = j;
						if (j == taille_c) {
							ligne_faite[k] = false;
							--n_restant;
						}
					}
				}

				ligne_faite[i_min] = false;
				--n_restant;
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

			std::vector<T> resultat(taille_c, faux);

			//une solution est possible. La calculer. Pour ceci, tout d'abord enlever les "parametres libres". (virtuellement)
			//ensuite, pour chaque premier �l�ment non-nul, = Y[k] / �l�ment de matrice
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
	};

	sous_ev<T>* noyau() const {
		static_assert(type_algebre<T>::type == 0);
		if constexpr (type_algebre<T>::approx == 0) {

			sous_ev<T>* resultat = new sous_ev<T>();
			matrice<T> m_matrice(*this);

			T vrai = unite(coeffs[0][0], true);
			T faux = unite(vrai, false);

			// on parcourt les lignes. Le premier �l�ment non-nul : annule toute la colonne
			for (int i(0); i < taille_l; ++i) {
				int j;
				for (j = 0; j < taille_c; ++j)
					if ((bool)m_matrice.coeffs[i][j])
						break;
				if (j == taille_c)
					continue;

				T inv = vrai / m_matrice.coeffs[i][j];
				for (int k(0); k < taille_l; ++k) {
					if (k == i)
						continue;
					if ((bool)m_matrice.coeffs[k][j])
						m_matrice.ajouterLigne(i, k, -inv * m_matrice.coeffs[k][j]); //ajoute la ligne i � la ligne k. annnule [k][j]
				}
			}

			//on s�pare les parametres libres. les colonnes nulles, et les �l�ments non-nuls apres un premier �l�ment non-nul.
			std::vector<int> parametres_libres(0);
			std::vector<int> parametres_non_libres(0);// on s'en fout !!!
			// on commence : premier �l�ment non-nul, et ceux qui suivent qui sont non-nuls sont des parametres libres.
			for (int i(0); i < taille_l; ++i) {
				int j;
				for (j = 0; j < taille_c; ++j)
					if ((bool)m_matrice.coeffs[i][j])
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
				for (int j(0); j < taille_l; ++j)
					if ((bool)m_matrice.coeffs[j][i]) {
						test = true;
						break;
					}
				if (test)
					continue;
				parametres_libres.push_back(i); //la colonne i est nulle
			}

			std::sort(parametres_libres.begin(), parametres_libres.end());
			parametres_libres.erase(std::unique(parametres_libres.begin(), parametres_libres.end()), parametres_libres.end()); //supprimer les doublons. A v�rifier.
			for (int vec_k(0); vec_k < parametres_libres.size(); ++vec_k) { //parametres libre : on en prend un � la fois.
				T temp = unite(m_matrice.coeffs[0][0], false);

				std::vector<T> vec(taille_c, temp);
				temp = faux;
				vec[parametres_libres[vec_k]] = temp; //faux partout sauf dans ce parametre libre.

				//on parcourt les premiers-�l�ments non-nuls. Et on calcule leurs valeurs.
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
		else {
			sous_ev<T>* resultat = new sous_ev<T>();
			matrice<T> m_matrice(*this);

			T vrai = unite(coeffs[0][0], true);
			T faux = unite(vrai, false);
			std::vector<bool> ligne_faite(taille_l, true);
			std::vector<int> prem_colonne(taille_l, 0);
			int n_restant = taille_l;

			for (int i(0); i < taille_l; ++i) {
				int j = 0;
				for (; j < taille_c; ++j)
					if ((bool)m_matrice.coeffs[i][j])
						break;
				prem_colonne[i] = j;
				if (j == taille_c) {
					ligne_faite[i] = false;
					--n_restant;
				}
			}

			while (n_restant > 0) {
				std::vector<decltype(norme(vrai))> normes(taille_c, norme(faux));
				std::vector<int> lignes(taille_c, -1);
				std::vector<int> colonnes(0);
				for (int i(0); i < taille_l; ++i) {
					if (!ligne_faite[i])
						continue;
					colonnes.push_back(prem_colonne[i]);
					if (lignes[prem_colonne[i]] == -1) {
						lignes[prem_colonne[i]] = i;
						normes[prem_colonne[i]] = norme(m_matrice.coeffs[i][prem_colonne[i]]);
					}
					else {
						auto norme_temp = norme(m_matrice.coeffs[i][prem_colonne[i]]);
						if (norme_temp > normes[prem_colonne[i]]) {
							normes[prem_colonne[i]] = norme_temp;
							lignes[prem_colonne[i]] = i;
						}
					}
				}

				std::sort(colonnes.begin(), colonnes.end());
				auto last = std::unique(colonnes.begin(), colonnes.end());
				colonnes.erase(last, colonnes.end());

				std::vector<decltype(norme(vrai))> normes_colonnes(taille_c, norme(faux));
				for (int k(0); k < colonnes.size(); ++k) {
					for (int i(0); i < taille_l; ++i) {
						if (i == lignes[colonnes[k]])
							continue;
						auto norme_temp = norme(m_matrice.coeffs[i][colonnes[k]]);
						if (norme_temp > normes_colonnes[colonnes[k]])
							normes_colonnes[colonnes[k]] = norme_temp;
					}
				}


				auto norme_min = norme(faux);
				int j_min = -1;
				for (int j(0); j < taille_c; ++j) {
					if (lignes[j] == -1) //VERIFIER
						continue;
					if (j_min == -1) {
						j_min = j;
						norme_min = normes_colonnes[j] / normes[j];
					}
					else {
						auto norme_temp = normes_colonnes[j] / normes[j];
						if (norme_temp < norme_min) {
							j_min = j;
							norme_min = norme_temp;
						}
					}
				}
				int i_min = lignes[j_min];

				T inv = vrai / m_matrice.coeffs[i_min][j_min];
				for (int k(0); k < taille_l; ++k) {
					if (k == i_min)
						continue;
					if (!(bool)m_matrice.coeffs[k][j_min])
						continue;
					m_matrice.ajouterLigne(i_min, k, -inv * m_matrice.coeffs[k][j_min]); //ajoute la ligne i � la ligne k. annnule [k][j_min]

					//ligne_faite et prem_colonne (taille_l)
					if (j_min == prem_colonne[k]) {
						int j = j_min;
						for (; j < taille_c; ++j)
							if ((bool)m_matrice.coeffs[k][j])
								break;
						prem_colonne[k] = j;
						if (j == taille_c) {
							ligne_faite[k] = false;
							--n_restant;
						}
					}
				}

				ligne_faite[i_min] = false;
				--n_restant;
			}
			//MATRICE PREPAREE !!!

			//on s�pare les parametres libres. les colonnes nulles, et les �l�ments non-nuls apres un premier �l�ment non-nul.
			std::vector<int> parametres_libres(0);
			std::vector<int> parametres_non_libres(0); //on s'en fiche en fait.
			// on commence : premier �l�ment non-nul, et ceux qui suivent qui sont non-nuls
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

			std::sort(parametres_libres.begin(), parametres_libres.end());
			parametres_libres.erase(std::unique(parametres_libres.begin(), parametres_libres.end()), parametres_libres.end()); //supprimer les doublons. A v�rifier.
			for (int vec_k(0); vec_k < parametres_libres.size(); ++vec_k) { //parametres libre : on en prend un � la fois.
				T temp = unite(m_matrice.coeffs[0][0], false);
				std::vector<T> vec(taille_c, temp);
				temp = unite(temp, true);
				vec[parametres_libres[vec_k]] = temp; //faux partout sauf dans ce parametre libre.

				//on parcourt les premiers-�l�ments non-nuls. Et on calcule leurs valeurs.
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
	};



	template<class U>
	explicit operator matrice<U>() {
		matrice<U> result(taille_l, taille_c);
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				result.coeffs[i][j] = (U)coeffs[i][j];
		return result;
	};

	template<class U> matrice<U> operator()(U element) {
		matrice<U> result(taille_l, taille_c);
		for (int i(0); i < taille_l; ++i)
			for (int j(0); j < taille_c; ++j)
				result.coeffs[i][j] = coeffs[i][j](element);
	};


	int taille_l;
	int taille_c;
	std::vector < std::vector< T>> coeffs;
};

template<class T>
matrice<T>& matrice<T>::operator=(const matrice<T>& temp) {
	if (this == &temp)
		return *this;
	taille_l = temp.taille_l;
	taille_c = temp.taille_c;
	coeffs = temp.coeffs;
	return *this;
};

/*

	//type = 1. OK
	T determinant() const {
		matrice<rationnel<T>> m_matrice(taille_l,taille_l);
		T vrai_ = unite(coeffs[0][0],true);

		if (taille_l != taille_c)
			throw std::domain_error("determinant matrice : les dimensions ne co�ncident pas");
		int taille = taille_l;


		for (int i(0); i < taille; ++i)
			for (int j(0); j < taille; ++j)
				m_matrice.coeffs[i][j] = rationnel<T>(coeffs[i][j], vrai_); //passer par les rationnels. vrai_ optionnel

		rationnel<T> det = m_matrice.determinant();

		T result = det.numerateur / det.denominateur; //simplifier : division avec reste 
		Simplifier_frac_poly(result);
		return result;
	};


	//type = 0 approx = 0. OK
	T determinant() const {
		matrice<T> m_matrice(*this);
		T det = unite(coeffs[0][0],true);
		T vrai = det;

		if (taille_l != taille_c)
			throw std::domain_error("determinant de matrice : les dimensions ne coincident pas");
		int taille = taille_l;

		for (int i(0); i < taille; ++i) {//colonne i
			int j;
			for (j = i; j < taille; ++j) {
				if ((bool)m_matrice.coeffs[j][i])
					break;
			}
			if (j == taille)
				return unite(det,false);

			if (i != j) {
				m_matrice.echangerLigne(i, j);
				det = -det;
			}

			T inv = vrai / m_matrice.coeffs[i][i];
			for (j = i + 1; j < taille; ++j)
				m_matrice.ajouterLigne(i, j, ((-m_matrice.coeffs[j][i]) * inv));  // ajouter ligne i * coeffs � la ligne j ...

		}

		for (int i(0); i < taille; ++i)
			det = (det * m_matrice.coeffs[i][i]);

		Simplifier_frac_poly(det);
		return det;
	};

	//type = 0, approx = 1. OK
	T determinant() const {
		matrice<T> m_matrice(*this);
		T det = unite(coeffs[0][0], true);
		T vrai = det;

		if (taille_l != taille_c)
			throw std::domain_error("determinant de matrice : les dimensions ne co�ncident pas");
		int taille = taille_l;


		int taille_max = max(taille - 4, 0); //pour les 4 derniers : on utilise le determinant_anneau
		int taille_fin = taille - taille_max;

		for (int i(0); i < taille_max; ++i) {
			int j_max = -1;
			auto norme_max = norme(vrai);
			for (int j = i; j < taille; ++j) {
				if ((bool) m_matrice.coeffs[i][j]) {
					auto norme_temp = norme(m_matrice.coeffs[j][i]);
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
			if (j_max == -1) //0 sur une colonne
				return unite(det, false);

			if (i != j_max) {
				m_matrice.echangerLigne(i, j_max);
				det = -det;
			}
			T inv = vrai / m_matrice.coeffs[i][i];
			for (int j = i + 1; j < taille; ++j) {
				m_matrice.ajouterLigne(i, j, ((-m_matrice.coeffs[j][i]) * inv)); // ajouter ligne i * coeffs � la ligne j ...
			}

		}

		for (int i(0); i < taille_max; ++i)
			det = (det * m_matrice.coeffs[i][i]);

		if (taille_fin > 0) {
			matrice<T> m_matrice_fin(taille_fin, taille_fin);
			for (int i(0); i < taille_fin; ++i)
				for (int j(0); j < taille_fin; ++j)
					m_matrice_fin.coeffs[i][j] = m_matrice.coeffs[taille_max + i][taille_max + j];
			det = det * m_matrice_fin.determinant_anneau();
		}
		Simplifier_frac_poly(det);
		return det;
	};

	
	//type = 0  approx = 0
	matrice<T> inverse() const {
		static_assert(type_algebre<T>::type == 0);
		T vrai = unite(coeffs[0][0],true);

		if( taille_l != taille_c)
			throw std::domain_error("inverse de matrice : dimensions ne coincident pas");
		int taille = taille_l;

		matrice<T> m_matrice(*this);
		matrice<T> resultat= unite(m_matrice,true);


		for (int i(0); i < taille; ++i) {
			int j;
			for (j = i; j < taille; ++j)
				if ((bool)coeffs[j][i])
					break;
			if (j == taille)
				throw std::domain_error("matrice non-inversible");

			if (i != j) {
				m_matrice.echangerLigne(i, j);
				resultat.echangerLigne(i, j);
			}
			resultat.multiplierLigne(i, vrai / m_matrice.coeffs[i][i]);
			m_matrice.multiplierLigne(i, vrai / m_matrice.coeffs[i][i]);

			for (j = 0; j < taille; ++j) {
				if (j == i)
					continue;
				resultat.ajouterLigne(i, j, -m_matrice.coeffs[j][i]);  // ajouter ligne i * coeffs � la ligne j ...
				m_matrice.ajouterLigne(i, j, -m_matrice.coeffs[j][i]);
			}
		}

		return resultat;
	};

	//type_algebre = 0, approx = 1. OK
	matrice<T> inverse() const {
		T vrai = unite(coeffs[0][0], true);

		if (taille_l != taille_c)
			throw std::domain_error("inverse de matrice : les dimensions ne co�ncident pas.");
		int taille = taille_l;

		matrice<T> m_matrice(*this);
		matrice<T> resultat = unite(*this, true); //identit�

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
			resultat.multiplierLigne(i, vrai / m_matrice.coeffs[i][i]);
			m_matrice.multiplierLigne(i, vrai / m_matrice.coeffs[i][i]);

			for (int j = 0; j < taille; ++j) {
				if (j == i)
					continue;
				resultat.ajouterLigne(i, j, -m_matrice.coeffs[j][i]);
				m_matrice.ajouterLigne(i, j, -m_matrice.coeffs[j][i]);
			}
		}

		return resultat;
	};


	//type_algebre = 0, approx=0
	std::vector<T> resoudre(std::vector<T> Y) {
		matrice<T> m_matrice(*this);
		if (m_matrice.taille_l != Y.size())
			throw std::domain_error("resolution equation lineaire : les dimensions ne correspondent pas");
		T vrai = unite(m_matrice.coeffs[0][0],true);

		// on parcourt les lignes. Le premier �l�ment non-nul : annule toute la colonne
		for (int i(0); i < taille_l; ++i) {
			int j;
			for (j = 0; j < taille_c; ++j)
				if ((bool)m_matrice.coeffs[i][j])
					break;
			if (j == taille_c)
				continue;

			T inv = vrai / m_matrice.coeffs[i][j]; //on est dans la ligne i, colonne j. On annule la colonne j (sauf ligne i)
			for (int k(0); k < taille_l; ++k) { //ligne k != i.
				if (k == i)
					continue;
				if ((bool) m_matrice.coeffs[k][j]) {
					Y[k] = Y[k] - inv * m_matrice.coeffs[k][j] * Y[i];
					m_matrice.ajouterLigne(i, k, -inv * m_matrice.coeffs[k][j]); //ajoute la ligne i � la ligne k. annnule [k][j]
				}
			}
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
			if ((bool) Y[i]) //ligne nulle, et Y non-nul. Renvoyer vide.
				return std::vector<T>(0);
		}

		T faux = unite(vrai, false);
		std::vector<T> resultat(taille_c,faux);

		//une solution est possible. La calculer. Pour ceci, tout d'abord enlever les "parametres libres". (virtuellement)
		//ensuite, pour chaque premier �l�ment non-nul, = Y[k]
		for (int i(0); i < taille_l; ++i) {
			int j;
			for (j = 0; j < taille_c; ++j)
				if ((bool)m_matrice.coeffs[i][j])
					break;
			if (j == taille_c) //ligne nulle. 
				continue;
			resultat[j] = Y[i] / m_matrice.coeffs[i][j];
		}

		return resultat;
	};


	//type_algebre = 0, approx = 1. A MODIFIER
	std::vector<T> resoudre(std::vector<T> Y) { //dans le cas approx.
		matrice<T> m_matrice(*this);
		if (m_matrice.taille_l != Y.size())
			throw std::domain_error("resolution equation lineaire : les dimensions ne correspondent pas");
		T vrai = unite(coeffs[0][0], true);
		T faux = unite(vrai, false);
		std::vector<bool> ligne_faite(taille_l, true);
		std::vector<int> prem_colonne(taille_l, 0);
		int n_restant = taille_l;

		for (int i(0); i < taille_l; ++i) {
			int j = 0;
			for (; j < taille_c; ++j)
				if ((bool)m_matrice.coeffs[i][j])
					break;
			prem_colonne[i] = j;
			if (j == taille_c) {
				ligne_faite[i] = false;
				--n_restant;
			}
		}

		while (n_restant > 0) {
			std::vector<decltype(norme(vrai))> normes(taille_c, norme(faux));
			std::vector<int> lignes(taille_c, -1);
			std::vector<int> colonnes(0);
			for (int i(0); i < taille_l; ++i) {
				if (!ligne_faite[i])
					continue;
				colonnes.push_back(prem_colonne[i]);
				if (lignes[prem_colonne[i]] == -1) {
					lignes[prem_colonne[i]] = i;
					normes[prem_colonne[i]] = norme(m_matrice.coeffs[i][prem_colonne[i]]);
				}
				else {
					auto norme_temp = norme(m_matrice.coeffs[i][prem_colonne[i]]);
					if (norme_temp > normes[prem_colonne[i]]) {
						normes[prem_colonne[i]] = norme_temp;
						lignes[prem_colonne[i]] = i;
					}
				}
			}

			std::sort(colonnes.begin(), colonnes.end());
			auto last = std::unique(colonnes.begin(), colonnes.end());
			colonnes.erase(last, colonnes.end());

			std::vector<decltype(norme(vrai))> normes_colonnes(taille_c, norme(faux));
			for (int k(0); k < colonnes.size(); ++k) {
				for (int i(0); i < taille_l; ++i) {
					if (i == lignes[colonnes[k])
						continue;
					auto norme_temp = norme(m_matrice.coeffs[i][colonnes[k]]);
					if (norme_temp > normes_colonnes[colonnes[k]])
						normes_colonnes[colonnes[k]] = norme_temp;
				}
			}


			auto norme_min = norme(faux);
			int j_min = -1;
			for (int j(0); j < taille_c; ++j) {
				if (lignes[j] == -1) //VERIFIER
					continue;
				if (j_min == -1) {
					j_min = j;
					norme_min = normes_colonnes[j] / normes[j];
				}
				else {
					auto norme_temp = normes_colonnes[j] / normes[j];
					if (norme_temp < norme_min) {
						j_min = j;
						norme_min = norme_temp;
					}
				}
			}
			int i_min = lignes[j_min];

			T inv = vrai / m_matrice.coeffs[i_min][j_min];
			for (int k(0); k < taille_l; ++k) {
				if (k == i_min)
					continue;
				if (! (bool) m_matrice.coeffs[k][j_min])
					continue;
				Y[k] = Y[k] - inv * m_matrice.coeffs[k][j_min] * Y[i_min]; // VERIFIER
				m_matrice.ajouterLigne(i_min, k, -inv * m_matrice.coeffs[k][j_min]); //ajoute la ligne i � la ligne k. annnule [k][j_min]

				//ligne_faite et prem_colonne (taille_l)
				if (j_min == prem_colonne[k]) {
					int j = j_min;
					for (; j < taille_c; ++j)
						if ((bool)m_matrice.coeffs[k][j])
							break;
					prem_colonne[k] = j;
					if (j == taille_c) {
						ligne_faite[k] = false;
						--n_restant;
					}
				}
			}

			ligne_faite[i_min] = false;
			--n_restant;
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

		std::vector<T> resultat(taille_c, faux);

		//une solution est possible. La calculer. Pour ceci, tout d'abord enlever les "parametres libres". (virtuellement)
		//ensuite, pour chaque premier �l�ment non-nul, = Y[k] / �l�ment de matrice
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

	//A MODIFIER type_algebre = 0 ; approx = 0
	sous_ev<T>* noyau() {
		sous_ev<T>* resultat = new sous_ev<T>();
		matrice<T> m_matrice(*this);

		T vrai = unite(coeffs[0][0],true);
		T faux = unite(vrai, false);

		// on parcourt les lignes. Le premier �l�ment non-nul : annule toute la colonne
		for (int i(0); i < taille_l; ++i) {
			int j;
			for (j = 0; j < taille_c; ++j)
				if ((bool)m_matrice.coeffs[i][j])
					break;
			if (j == taille_c)
				continue;

			T inv = vrai / m_matrice.coeffs[i][j];
			for (int k(0); k < taille_l; ++k) {
				if (k == i)
					continue;
				if((bool) m_matrice.coeffs[k][j])
					m_matrice.ajouterLigne(i, k, -inv * m_matrice.coeffs[k][j]); //ajoute la ligne i � la ligne k. annnule [k][j]
			}
		}

		//on s�pare les parametres libres. les colonnes nulles, et les �l�ments non-nuls apres un premier �l�ment non-nul.
		std::vector<int> parametres_libres(0);
		std::vector<int> parametres_non_libres(0);// on s'en fout !!!
		// on commence : premier �l�ment non-nul, et ceux qui suivent qui sont non-nuls sont des parametres libres.
		for (int i(0); i < taille_l; ++i) {
			int j;
			for (j = 0; j < taille_c; ++j)
				if ((bool)m_matrice.coeffs[i][j])
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
				if ((bool)m_matrice.coeffs[j][i]) {
					test = true;
					break;
				}
			if (test)
				continue;
			parametres_libres.push_back(i); //la colonne i est nulle
		}

		std::sort(parametres_libres.begin(), parametres_libres.end());
		parametres_libres.erase(std::unique(parametres_libres.begin(), parametres_libres.end()), parametres_libres.end()); //supprimer les doublons. A v�rifier.
		for (int vec_k(0); vec_k < parametres_libres.size(); ++vec_k) { //parametres libre : on en prend un � la fois.
			T temp = unite(m_matrice.coeffs[0][0],false);

			std::vector<T> vec(taille_c, temp);
			temp = faux;
			vec[parametres_libres[vec_k]] = temp; //faux partout sauf dans ce parametre libre.

			//on parcourt les premiers-�l�ments non-nuls. Et on calcule leurs valeurs.
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



	//type_algebre = 0, approx = 1
	sous_ev<T>* noyau() {
		sous_ev<T>* resultat = new sous_ev<T>();
		matrice<T> m_matrice(*this);

		T vrai = unite(coeffs[0][0], true);
		T faux = unite(vrai, false);
		std::vector<bool> ligne_faite(taille_l, true);
		std::vector<int> prem_colonne(taille_l, 0);
		int n_restant = taille_l;

		for (int i(0); i < taille_l; ++i) {
			int j = 0;
			for (; j < taille_c; ++j)
				if ((bool)m_matrice.coeffs[i][j])
					break;
			prem_colonne[i] = j;
			if (j == taille_c) {
				ligne_faite[i] = false;
				--n_restant;
			}
		}

		while (n_restant > 0) {
			std::vector<decltype(norme(vrai))> normes(taille_c, norme(faux));
			std::vector<int> lignes(taille_c, -1);
			std::vector<int> colonnes(0);
			for (int i(0); i < taille_l; ++i) {
				if (!ligne_faite[i])
					continue;
				colonnes.push_back(prem_colonne[i]);
				if (lignes[prem_colonne[i]] == -1) {
					lignes[prem_colonne[i]] = i;
					normes[prem_colonne[i]] = norme(m_matrice.coeffs[i][prem_colonne[i]]);
				}
				else {
					auto norme_temp = norme(m_matrice.coeffs[i][prem_colonne[i]]);
					if (norme_temp > normes[prem_colonne[i]]) {
						normes[prem_colonne[i]] = norme_temp;
						lignes[prem_colonne[i]] = i;
					}
				}
			}

			std::sort(colonnes.begin(), colonnes.end());
			auto last = std::unique(colonnes.begin(), colonnes.end());
			colonnes.erase(last, colonnes.end());

			std::vector<decltype(norme(vrai))> normes_colonnes(taille_c, norme(faux));
			for (int k(0); k < colonnes.size(); ++k) {
				for (int i(0); i < taille_l; ++i) {
					if (i == lignes[colonnes[k])
						continue;
					auto norme_temp = norme(m_matrice.coeffs[i][colonnes[k]]);
					if (norme_temp > normes_colonnes[colonnes[k]])
						normes_colonnes[colonnes[k]] = norme_temp;
				}
			}


			auto norme_min = norme(faux);
			int j_min = -1;
			for (int j(0); j < taille_c; ++j) {
				if (lignes[j] == -1) //VERIFIER
					continue;
				if (j_min == -1) {
					j_min = j;
					norme_min = normes_colonnes[j] / normes[j];
				}
				else {
					auto norme_temp = normes_colonnes[j] / normes[j];
					if (norme_temp < norme_min) {
						j_min = j;
						norme_min = norme_temp;
					}
				}
			}
			int i_min = lignes[j_min];

			T inv = vrai / m_matrice.coeffs[i_min][j_min];
			for (int k(0); k < taille_l; ++k) {
				if (k == i_min)
					continue;
				if (!(bool)m_matrice.coeffs[k][j_min])
					continue;
				m_matrice.ajouterLigne(i_min, k, -inv * m_matrice.coeffs[k][j_min]); //ajoute la ligne i � la ligne k. annnule [k][j_min]

				//ligne_faite et prem_colonne (taille_l)
				if (j_min == prem_colonne[k]) {
					int j = j_min;
					for (; j < taille_c; ++j)
						if ((bool)m_matrice.coeffs[k][j])
							break;
					prem_colonne[k] = j;
					if (j == taille_c) {
						ligne_faite[k] = false;
						--n_restant;
					}
				}
			}

			ligne_faite[i_min] = false;
			--n_restant;
		}
		//MATRICE PREPAREE !!!

		//on s�pare les parametres libres. les colonnes nulles, et les �l�ments non-nuls apres un premier �l�ment non-nul.
		std::vector<int> parametres_libres(0);
		std::vector<int> parametres_non_libres(0); //on s'en fiche en fait.
		// on commence : premier �l�ment non-nul, et ceux qui suivent qui sont non-nuls
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

		std::sort(parametres_libres.begin(), parametres_libres.end());
		parametres_libres.erase(std::unique(parametres_libres.begin(), parametres_libres.end()), parametres_libres.end()); //supprimer les doublons. A v�rifier.
		for (int vec_k(0); vec_k < parametres_libres.size(); ++vec_k) { //parametres libre : on en prend un � la fois.
			T temp = unite(m_matrice.coeffs[0][0], false);
			std::vector<T> vec(taille_c, temp);
			temp = unite(temp, true);
			vec[parametres_libres[vec_k]] = temp; //faux partout sauf dans ce parametre libre.

			//on parcourt les premiers-�l�ments non-nuls. Et on calcule leurs valeurs.
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


	*/
	
template<class T> class sous_ev {
public:
	sous_ev() : liste_vecteurs(0){
	}

	void ajouter_vecteur(std::vector<T> vecteur) {
		liste_vecteurs.push_back(vecteur);
	}

	std::vector<std::vector<T>> liste_vecteurs;
};



// type_algebre<T>::type == 0;
// type_algebre<T>::approx == 0;
