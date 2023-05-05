#pragma once
#include "polynome.hpp"
#include "matrice.hpp"
#include "anneau quotient.hpp"
#include "fonctions template.hpp"
#include "types.hpp"
#include "norme.hpp"

#include <vector>
#include <algorithm>
#include <iostream>


template<class T> class sous_ev_condition; 
template<typename T, class enable = void> class calcul_condition;


template<class T> class calcul_condition<T,typename std::enable_if_t<type_algebre<T>::approx==0> > {
public:
	std::vector<polynome<T>> conditions; //polynomes conditions
	std::vector<matrice<anneau_quotient<polynome<T>>>> matrices; //là où on en est dans le calcul
	std::vector<int> ligne; // on se trouve à cette ligne ...
	std::vector<std::vector<int>> colonnes;

	sous_ev_condition<T>* m_sous_ev_condition;
	polynome<T> chi;


	static_assert(type_algebre<T>::type == 0);

	calcul_condition(matrice<T> m_matrice,bool vecteurs_propres = false) { //vecteurs propre : choisit de résoudre A.X=0 ou A^n.X=0
#ifdef ALGEBRA_USE_EXCEPTION
		if (m_matrice.taille_l != m_matrice.taille_c)
			throw std::domain_error("diagonalisation : dimensions ne coincident pas");
#endif
		polynome<T> chi = m_matrice.polynomeCaracteristique();
		polynome<T> chi_ = chi;

		/*
		polynome<T> R;
		while ((R = PGCD(chi_, 
		(chi_))).degre >= 1) {
			chi_ = chi_ / R;
		}//chi_ à racines simples */
		chi_ = chi_ / PGCD(chi_, derivee(chi_));

		int n = 1;
		polynome<T> chi_n = chi_;
		while ((bool) chi_n(m_matrice)) {
			++n;
			chi_n = chi_n * chi_;
		}
		if (vecteurs_propres)
			n = 1;

		int taille = m_matrice.taille_l;
		T e = unite(m_matrice.coeffs[0][0],true);
		anneau_quotient<polynome<T>> poly_true(polynome<T>(e), chi_);
		e = -e;
		matrice<anneau_quotient<polynome<T>>> matrice_poly(taille, taille,poly_true);
		for (int i(0); i < taille; ++i)
			for (int j(0); j < taille; ++j)
				matrice_poly.coeffs[i][j].element = polynome<T>(m_matrice.coeffs[i][j]);
		for (int i(0); i < taille; ++i)
			matrice_poly.coeffs[i][i].element = polynome<T>(m_matrice.coeffs[i][i], e);


		conditions = { chi_ };
		matrice_poly = puissance(matrice_poly, n);
		matrices = { matrice_poly };
		ligne = { 0 };
		colonnes = { std::vector<int>(m_matrice.taille_l,-1) };

		m_sous_ev_condition = NULL;
	};

	bool ajouter(polynome<T> const& Q, int i) { //ajoute une condition Q ... tous les polynomes condition divisent chi ! par récurrence (on commence avec chi). 
		//		if (aleatoire)
		//			return (bool) Q; // à vérifier ... je crois OK
		if (Q.degre <= 0)
			return (bool)Q;  // Q constant ... soit nul soit jamais nul.
		polynome<T> R = PGCD(Q, conditions[i]);

		if (R.degre <= 0) {
			return true; // l'élément testé n'est jamais nul ! pas de racine en commun avec la condition
		}
		else if (R.degre == conditions[i].degre) {
			return false;
		}
		else {
			matrices.push_back(matrices[i]);
			conditions.push_back(R);
			ligne.push_back(ligne[i]);
			colonnes.push_back(colonnes[i]);
			mettre_a_jour_quotient(matrices[matrices.size() - 1], R);

			conditions[i] = conditions[i] / R;
			conditions[i].normaliser();
			mettre_a_jour_quotient(matrices[i], conditions[i]);

			return true; //l'élément testé n'est jamais nul dans ce cas : la nouvelle condition n'a aucune racine en commun avec l'élément : on a séparé.
			//Mais la condition ajoutée annule l'élément testé :  dans l'espace ajouté, Q sera nul (la condition ajoutée est R, et Q = R*quelquechose)
		}
	};

	void calculer(int iter) {
//		matrice < anneau_quotient<polynome<T>>>& matrices[i] = matrices[i];
		int taille = matrices[0].taille_l;
		int i = ligne[iter];

		for (; i < taille; ++i) {
			ligne[iter] = i;
			int j = 0;
			for (; j < taille; ++j)
				if (ajouter(matrices[iter].coeffs[i][j].element, iter))
					break;
			if (j == taille)
				continue;
			colonnes[iter][i] = j;

			polynome<T> inv = inverse(matrices[iter].coeffs[i][j].element, conditions[iter]);
			for (int k(0); k < taille; ++k) {
				if (k == i)
					continue;
				if ((bool)matrices[iter].coeffs[k][j].element)
					matrices[iter].ajouterLigne(i, k, -inv * matrices[iter].coeffs[k][j]); //ajoute la ligne i à la ligne k. annnule [k][j]
			}
		}

		//on sépare les parametres libres. les colonnes nulles, et les éléments non-nuls apres un premier élément non-nul.
		std::vector<int> parametres_libres(0);
		parametres_libres.reserve(taille);
		// on commence : premier élément non-nul, et ceux qui suivent qui sont non-nuls sont des parametres libres.
		for (int i(0); i < taille; ++i) {
			for (int j(0); j < taille; ++j) {
				if (j == colonnes[iter][i])
					continue;
				if ((bool)matrices[iter].coeffs[i][j].element)
					parametres_libres.push_back(j);
			}
		}

		//puis, les colonnes nulles : parametres libres.
		for (int i(0); i < taille; ++i) {
			bool test = false;
			for (int j(0); j < taille; ++j)
				if ((bool)matrices[iter].coeffs[j][i].element) {
					test = true;
					break;
				}
			if (test)
				continue;
			parametres_libres.push_back(i); //la colonne i est nulle
		}

		std::sort(parametres_libres.begin(), parametres_libres.end());
		parametres_libres.erase(std::unique(parametres_libres.begin(), parametres_libres.end()), parametres_libres.end()); //supprimer les doublons. A vérifier.

		for (int vec_k(0); vec_k < parametres_libres.size(); ++vec_k) { //parametres libre : on en prend un à la fois.
			anneau_quotient<polynome<T>> faux = unite(matrices[iter].coeffs[0][0], false);
			anneau_quotient<polynome<T>> vrai = unite(faux, true);

			std::vector<anneau_quotient<polynome<T>>> vec(taille, faux);
			vec[parametres_libres[vec_k]] = vrai; //faux partout sauf dans ce parametre libre.

			//on parcourt les premiers-éléments non-nuls. Et on calcule leurs valeurs.
			for (int i(0); i < taille; ++i) {
				int j = colonnes[iter][i];
				if (j == -1)
					continue;

				vec[j] = -inverse(matrices[iter].coeffs[i][j].element, conditions[iter]) * matrices[iter].coeffs[i][parametres_libres[vec_k]];
			}

			std::vector<polynome<T>> vec_final(taille, faux.element);
			for (int j(0); j < taille; ++j)
				vec_final[j] = vec[j].element;

			m_sous_ev_condition->ajouter_vecteur(conditions[iter], vec_final);
		}



		return;

	};

	sous_ev_condition<T>* calculer() {
		if (m_sous_ev_condition) 
			delete m_sous_ev_condition;

		m_sous_ev_condition = new sous_ev_condition<T>;

		for(int i(0);i<conditions.size();++i)
			calculer(i);
		//on calcul les matrices reliées aux conditions qu'on ajoute au fur et à mesure !

		return m_sous_ev_condition;
	};


	void mettre_a_jour_quotient(matrice<anneau_quotient<polynome<T>>>& m_matrice, polynome<T> Q) {
		for (int i(0); i < m_matrice.taille_l; ++i)
			for (int j(0); j < m_matrice.taille_l; ++j)
				m_matrice.coeffs[i][j].changer_quotient(Q);
		return;
	}

};

template<class T> class calcul_condition<T, typename std::enable_if_t<type_algebre<T>::approx == 1>> {
public:
	std::vector<polynome<T>> conditions; //polynomes conditions
	std::vector<matrice<anneau_quotient<polynome<T>>>> matrices; //là où on en est dans le calcul
	std::vector<int> ligne; // on se trouve à cette ligne ...
	std::vector<std::vector<int>> colonnes;

	sous_ev_condition<T>* m_sous_ev_condition;
	polynome<T> chi;
	using norme_type = decltype(norme(T()));

	static_assert(type_algebre<T>::type == 0);

	calcul_condition(matrice<T> m_matrice, bool vecteurs_propres = false) { //vecteurs propre : choisit de résoudre A.X=0 ou A^n.X=0
#ifdef ALGEBRA_USE_EXCEPTION
		if (m_matrice.taille_l != m_matrice.taille_c)
			throw std::domain_error("diagonalisation : dimensions ne coincident pas");
#endif
		polynome<T> chi = m_matrice.polynomeCaracteristique();
		polynome<T> chi_ = chi;

		/*
		polynome<T> R;
		while ((R = PGCD(chi_,
		(chi_))).degre >= 1) {
			chi_ = chi_ / R;
		}//chi_ à racines simples */
		chi_ = chi_ / PGCD(chi_, derivee(chi_));

		int n = 1;
		polynome<T> chi_n = chi_;
		while ((bool)chi_n(m_matrice)) {
			++n;
			chi_n = chi_n * chi_;
		}
		if (vecteurs_propres)
			n = 1;

		int taille = m_matrice.taille_l;
		T e = unite(m_matrice.coeffs[0][0], true);
		anneau_quotient<polynome<T>> poly_true(polynome<T>(e), chi_);
		e = -e;
		matrice<anneau_quotient<polynome<T>>> matrice_poly(taille, taille, poly_true);
		for (int i(0); i < taille; ++i)
			for (int j(0); j < taille; ++j)
				matrice_poly.coeffs[i][j].element = polynome<T>(m_matrice.coeffs[i][j]);
		for (int i(0); i < taille; ++i)
			matrice_poly.coeffs[i][i].element = polynome<T>(m_matrice.coeffs[i][i], e);


		conditions = { chi_ };
		matrice_poly = puissance(matrice_poly, n);
		matrices = { matrice_poly };
		ligne = { 0 };
		colonnes = { std::vector<int>(m_matrice.taille_l,-1) };

		m_sous_ev_condition = NULL;
};

	norme_type ajouter(polynome<T> const& Q, int i) { //ajoute une condition Q ... tous les polynomes condition divisent chi ! par récurrence (on commence avec chi). 
		//		if (aleatoire)
		//			return (bool) Q; // à vérifier ... je crois OK
		if (Q.degre <= 0)
			return (bool)Q;  // Q constant ... soit nul soit jamais nul.
		polynome<T> R = PGCD(Q, conditions[i]);

		if (R.degre <= 0) {
			return norme(PGCD(inverse(Q,conditions[i]),conditions[i])); // l'élément testé n'est jamais nul ! pas de racine en commun avec la condition
		}
		else if (R.degre == conditions[i].degre)
			return 0.;
		else {
			matrices.push_back(matrices[i]);
			conditions.push_back(R);
			ligne.push_back(ligne[i]);
			colonnes.push_back(colonnes[i]);
			mettre_a_jour_quotient(matrices[matrices.size() - 1], R);

			conditions[i] = conditions[i] / R;
			conditions[i].normaliser();
			mettre_a_jour_quotient(matrices[i], conditions[i]);

			return norme(PGCD(inverse(Q, conditions[i]), conditions[i])); //l'élément testé n'est jamais nul dans ce cas : la nouvelle condition n'a aucune racine en commun avec l'élément : on a séparé.
			//Mais la condition ajoutée annule l'élément testé :  dans l'espace ajouté, Q sera nul (la condition ajoutée est R, et Q = R*quelquechose)
		};
	};

	void calculer(int iter) {
		//		matrice < anneau_quotient<polynome<T>>>& matrices[i] = matrices[i];
		int taille = matrices[0].taille_l;
		int i = ligne[iter];
		polynome<T> poly_faux = unite(matrices[0].coeffs[0][0].element, false);

		for (; i < taille; ++i) {
			ligne[iter] = i;
			int j_max = -1;
			norme_type norme_min = 0.;
			for (int j = 0; j < taille; ++j) {
				norme_type norme_temp = ajouter(matrices[iter].coeffs[i][j].element, iter);
				if (norme_temp == 0.)
					continue;
				if (j_max != -1) {
					if (norme_temp < norme_min) {
						norme_min = norme_temp;
						j_max = j;
					}
				}
				else {
					norme_min = norme_temp;
					j_max = j;
				}
			}
			if (j_max == -1)
				continue;
			colonnes[iter][i] = j_max;

			polynome<T> inv = inverse(matrices[iter].coeffs[i][j_max].element, conditions[iter]);
			for (int k(0); k < taille; ++k) {
				if (k == i)
					continue;
				if ((bool)matrices[iter].coeffs[k][j_max].element) {
					matrices[iter].ajouterLigne(i, k, -inv * matrices[iter].coeffs[k][j_max]); //ajoute la ligne i à la ligne k. annnule [k][j]
					matrices[iter].coeffs[k][j_max].element = poly_faux;
				}
			}
		}

		//on sépare les parametres libres. les colonnes nulles, et les éléments non-nuls apres un premier élément non-nul.
		std::vector<int> parametres_libres(0);
		parametres_libres.reserve(taille);
		// on commence : premier élément non-nul, et ceux qui suivent qui sont non-nuls sont des parametres libres.
		for (int i(0); i < taille; ++i) {
			for (int j(0); j < taille; ++j) {
				if (j == colonnes[iter][i])
					continue;
				if ((bool)matrices[iter].coeffs[i][j].element)
					parametres_libres.push_back(j);
			}
		}

		//puis, les colonnes nulles : parametres libres.
		for (int i(0); i < taille; ++i) {
			bool test = false;
			for (int j(0); j < taille; ++j)
				if ((bool)matrices[iter].coeffs[j][i].element) {
					test = true;
					break;
				}
			if (test)
				continue;
			parametres_libres.push_back(i); //la colonne i est nulle
		}

		std::sort(parametres_libres.begin(), parametres_libres.end());
		parametres_libres.erase(std::unique(parametres_libres.begin(), parametres_libres.end()), parametres_libres.end()); //supprimer les doublons. A vérifier.

		for (int vec_k(0); vec_k < parametres_libres.size(); ++vec_k) { //parametres libre : on en prend un à la fois.
			anneau_quotient<polynome<T>> faux = unite(matrices[iter].coeffs[0][0], false);
			anneau_quotient<polynome<T>> vrai = unite(faux, true);

			std::vector<anneau_quotient<polynome<T>>> vec(taille, faux);
			vec[parametres_libres[vec_k]] = vrai; //faux partout sauf dans ce parametre libre.

			//on parcourt les premiers-éléments non-nuls. Et on calcule leurs valeurs.
			for (int i(0); i < taille; ++i) {
				int j = colonnes[iter][i];
				if (j == -1)
					continue;

				vec[j] = -inverse(matrices[iter].coeffs[i][j].element, conditions[iter]) * matrices[iter].coeffs[i][parametres_libres[vec_k]];
			}

			std::vector<polynome<T>> vec_final(taille, faux.element);
			for (int j(0); j < taille; ++j)
				vec_final[j] = vec[j].element;

			m_sous_ev_condition->ajouter_vecteur(conditions[iter], vec_final);
		}

		return;

	};

	sous_ev_condition<T>* calculer() {
		if (m_sous_ev_condition)
			delete m_sous_ev_condition;

		m_sous_ev_condition = new sous_ev_condition<T>;

		for (int i(0); i < conditions.size(); ++i)
			calculer(i);
		//on calcul les matrices reliées aux conditions qu'on ajoute au fur et à mesure !

		return m_sous_ev_condition;
	};


	void mettre_a_jour_quotient(matrice<anneau_quotient<polynome<T>>>& m_matrice, polynome<T> Q) {
		for (int i(0); i < m_matrice.taille_l; ++i)
			for (int j(0); j < m_matrice.taille_l; ++j)
				m_matrice.coeffs[i][j].changer_quotient(Q);
		return;
	}

};


template<class T> class sous_ev_condition {
public :
	
	void ajouter_vecteur(polynome<T> condition, std::vector<polynome<T>> vecteur) {
		liste_vecteurs.push_back(vecteur);
		liste_conditions.push_back(condition);
	};

	sous_ev_condition() {
		liste_vecteurs.resize(0);
		liste_conditions.resize(0);
	};

	//fonction évaluer ... utilise racines d'un polynome.
	//est-ce que les vecteurs sont indépendants, pour chaque racine de la condition ?
	//ie, est-ce que les vecteurs sont bien libres étant donnée une condition ?

	int nombre_solutions() {
		int n = 0;
		for (int i(0); i < liste_conditions.size(); ++i)
			n += liste_conditions[i].degre();
		return n;
	};

	friend std::ostream& operator<<(std::ostream& os, const sous_ev_condition<T>& element) {
		for (int i(0); i < element.liste_vecteurs.size(); ++i) {
			os << "condition : " << element.liste_conditions[i] << " ; vecteur : (";
			for (int j(0); j < element.liste_vecteurs[i].size(); ++j) {
				os << element.liste_vecteurs[i][j];
				if (j < element.liste_vecteurs[i].size() - 1)
					os << ",";
			}
			os << ")" << std::endl;
		}
		return os;
	};

	std::vector<std::vector<polynome<T>>> liste_vecteurs;
	std::vector<polynome<T>> liste_conditions; //associées aux vecteurs

};

//epsilon : polynome.coeff et (bool) matrice.
//polynome<T> à n variables
//racines d'un polynome<double>

//diviseur de zero : friend ; rationnel  a==b
//dérivée : friend (n'a pas de sens pour _quotient)

//class double_erreur_b

//matrice<T> * vecteur<U>.
// => instaurer la multiplication par U à gauche comme à droite. pour T*T rien ne change.
//operator %=

