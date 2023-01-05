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

	calcul_condition(matrice<T> m_matrice,bool vecteurs_propres = false) { //vecteurs propre : choisit de résoudre A.X=0 ou A^n.X=0
		if (m_matrice.taille_l != m_matrice.taille_c)
			throw std::domain_error("diagonalisation : dimensions ne coincident pas");

		polynome<T> chi = m_matrice.polynomeCaracteristique();
		polynome<T> chi_ = chi;

		/*
		polynome<T> R;
		while ((R = PGCD(chi_, derivee(chi_))).degre >= 1) {
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
		T e = m_matrice.coeffs[0][0];
		e = true;
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
		iterateurs = { 0 };

		m_sous_ev_condition = NULL;
	};

	inline polynome<T> separer(polynome<T>& P, polynome<T>& Q) {
		polynome<T> R = PGCD(P, Q);
		if (R.degre <= 0)
			return R;
		P = P / R;
		Q = Q / R;
		polynome<T> R1;
		while ((R1 = PGCD(P, R)).degre >= 1)
			P = P / R1;

		return R;
	};

	bool ajouter(polynome<T> Q, int i) { //ajoute une condition Q ... tous les polynomes condition divisent chi ! par récurrence (on commence avec chi). 
//		if (aleatoire)
//			return (bool) Q; // à vérifier ... je crois OK
		if (Q.degre <= 0)
			return (bool) Q;  // Q constant ... soit nul soit jamais nul.
		
		polynome<T> poly_save = conditions[i];
		polynome<T> R = separer(Q, conditions[i]); //modifie la condition : référence

		if (R.degre <= 0) {
				return true; // l'élément testé n'est jamais nul ! pas de racine en commun avec la condition
		}
		else {
			if (conditions[i].degre <= 0) {
//				conditions[i] = R * conditions[i];
				conditions[i] = poly_save; //pour diminuer l'erreur. cas approx
				return false; //l'élément testé est toujours nul dans ce cas : Q = R * qqchose, et R est la nouvelle condition
			}
			else {
				matrices.push_back(matrices[i]);
				conditions.push_back(R);
				iterateurs.push_back(iterateurs[i]);

				mettre_a_jour_quotient(matrices[matrices.size() - 1], R);
				mettre_a_jour_quotient(matrices[i], conditions[i]);

				return true; //l'élément testé n'est jamais nul dans ce cas : la nouvelle condition n'a aucune racine en commun avec l'élément : on a séparé.
				//Mais la condition ajoutée annule l'élément testé :  dans l'espace ajouté, Q sera nul (la condition ajoutée est R, et Q = R*quelquechose)
			}
		}
	};

	void calculer(int i) {
//		matrice < anneau_quotient<polynome<T>>>& matrices[i] = matrices[i];
		polynome<T>& poly = conditions[i];
		int taille = matrices[i].taille_l;
		long question;

		for (int j = iterateurs[i]; j < taille; ++j) { //colonne j ...
			iterateurs[i] = j;
			int k = j;
			do {

				if (! ajouter(matrices[i].coeffs[k][j].element, i)) //si élément nul ...
					++k;
				else
					break;

			} while (k < taille);

			if (k == taille) //la colonne j nulle sous la diagonale
				continue;

			if (k != j)
				matrices[i].echangerLigne(j, k);


			polynome<T> inverse_ = inverse(matrices[i].coeffs[j][j].element, conditions[i]); //condition changée

			for (int k(0); k < taille; ++k) {
				if (k == j)
					continue;
				matrices[i].ajouterLigne(j, k, -(inverse_ * matrices[i].coeffs[k][j]));
			}
		}

		//il y a deux parties : une avec une diagonale non-nulle, une autre triangulaire supérieur : 
		//on peut vérifier que comme en-dessous de la diagonale c'est nul, cette partie reste triangulaire supérieur.
		// 
		//une partie a des 0 sur la diagonale ... 
		//on cree la liste des i qui ont 0 sur la diagonale
		std::vector<int> non_diagonale(0); 
		for (int j(0); j < taille; ++j)
			if (!(bool)matrices[i].coeffs[j][j]) //ces éléments ont été testés : toujours nuls ou jamais nuls ...
				non_diagonale.push_back(j); //si nul, on ajoute
		//non_diagonale est trié (ordre croissant), et triangulaire supérieur et diagonale nulle.


		//on cherche sur chaque ligne le premier élément non-nul ... on le teste bien-sûr
		for(int j(non_diagonale.size()-2) ; j>=0;--j)
			for(int k(j+1); k<non_diagonale.size();++k)
				if (ajouter(matrices[i].coeffs[non_diagonale[j]][non_diagonale[k]].element, i)) { //si l'élément n'est jamais nul ... On calcule son inverse etc.
					polynome<T> inverse_ = inverse(matrices[i].coeffs[non_diagonale[j]][non_diagonale[k]].element, conditions[i]); //condition changée
					for (int l(0); l < non_diagonale.size()-1; ++l) {
						if (l == j)
							continue;
						matrices[i].ajouterLigne(non_diagonale[j], non_diagonale[l], -(inverse_ * matrices[i].coeffs[non_diagonale[l]][non_diagonale[k]]));
						//si ce premier élément est non-nul, on supprime toute la colonne.
					}
					break;
				}

		//on ajoute les colonnes nulles exactement dans la liste des parametres libres. (dans le sous-ensemble non-diagonale)
		std::vector<int> parametres_libres(0);
		for (int j = 0; j < non_diagonale.size(); ++j) {
			bool test = true;
			for (int k(0); k < non_diagonale.size(); ++k)
				if ((bool)matrices[i].coeffs[non_diagonale[k]][non_diagonale[j]]){
					test = false;
					break;
				}
			if (test)
				parametres_libres.push_back(j);
		}

		//on ajoute les colonnes non-nulles (pas forcément exactement : sinon elles seraient dans la liste précédente)
		//à droite d'un premier nombre. (dans le sous-ensemble non-diagonale). Ce premier nombre est déjà testé : toujours non-nul.
		for (int j(0); j < non_diagonale.size(); ++j) {
			bool trouve = false;
			for (int k(j + 1); k < non_diagonale.size(); ++k)
				if (!trouve) { //trouvé le premier non-nul
					if ((bool) matrices[i].coeffs[non_diagonale[j]][non_diagonale[k]]) {
						trouve = true;
						continue;
					}
				}
				else {
					if ((bool)matrices[i].coeffs[non_diagonale[j]][non_diagonale[k]]) //on ajoute.
						parametres_libres.push_back(k);
				}
		}
		//en gros on a tout ajouté sauf le premier nombre non-nul de chaque ligne. dans le sous-espace non-diagonal.

		parametres_libres.erase(std::unique(parametres_libres.begin(), parametres_libres.end()), parametres_libres.end()); //supprimer les doublons. A vérifier.

		for (int vec_k(0); vec_k < parametres_libres.size(); ++vec_k) { //parametres libre c'est dans le sous-espace sans diagonale.
			anneau_quotient<polynome<T>> temp = matrices[i].coeffs[0][0];
			temp = false;

			std::vector<anneau_quotient<polynome<T>>> vec_non_diagonale(non_diagonale.size(),temp);
			temp = true;
			vec_non_diagonale[parametres_libres[vec_k]] = temp; //faux partout sauf dans ce parametre libre.

			//on regarde, pour chaque premier élément non-nul, sa valeur ... un seul parametre libre non-nul.
			for (int k(0); k < non_diagonale.size()-1; ++k)
				for(int j(k+1);j<non_diagonale.size();++j)
					if ((bool)matrices[i].coeffs[non_diagonale[k]][non_diagonale[j]]) { //un seul parametre libre vaut 1, les autres 0. Sur chaque ligne un seul parametre libre, donc ...
						vec_non_diagonale[j] = - (inverse(matrices[i].coeffs[non_diagonale[k]][non_diagonale[j]].element, conditions[i]) * matrices[i].coeffs[non_diagonale[k]][non_diagonale[parametres_libres[vec_k]]]);//condition changée
						break; // l'inverse existe : l'élément est testé et/donc toujours non-nul.
					}//dans les faits, ne compte que si l'élément parametre libre est à droite du premier élément non-nul.

			temp = false;
			//le vecteur entier ... on commence par recopier la partie déjà calculée.
			std::vector<anneau_quotient<polynome<T>>> vec(taille, temp);
			for (int j(0); j < non_diagonale.size(); ++j)
				vec[non_diagonale[j]] = vec_non_diagonale[j];

			for(int j(0);j<taille;++j)
				if ((bool)matrices[i].coeffs[j][j]) {// les éléments de diagonale non-nulle : on calcule leur valeur.
					anneau_quotient<polynome<T>> somme = matrices[i].coeffs[0][0];
					somme = false;
					for (int k(0); k < non_diagonale.size(); ++k) {
						somme = somme + (matrices[i].coeffs[j][non_diagonale[k]] * vec_non_diagonale[k]);
					}
					vec[j] = -(inverse(matrices[i].coeffs[j][j].element, conditions[i]) * somme);
				}
			
			//on en extrait seulement les polynomes.
			std::vector<polynome<T>> vec_final(taille, temp.element);
			for (int j(0); j < taille; ++j)
				vec_final[j] = vec[j].element;


			m_sous_ev_condition->ajouter_vecteur(conditions[i], vec_final);
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

	std::vector<polynome<T>> conditions; //polynomes conditions
	std::vector<matrice<anneau_quotient<polynome<T>>>> matrices; //là où on en est dans le calcul
	std::vector<int> iterateurs; // on se trouve à cette colonne ...

	sous_ev_condition<T>* m_sous_ev_condition;
	polynome<T> chi;
};

template<class T> class calcul_condition<T, typename std::enable_if_t<type_algebre<T>::approx == 1>> {
public:

	calcul_condition(matrice<T> m_matrice, bool vecteurs_propres = false) { //vecteurs propre : choisit de résoudre A.X=0 ou A^n.X=0
		if (m_matrice.taille_l != m_matrice.taille_c)
			throw std::domain_error("diagonalisation : dimensions ne coincident pas");

		polynome<T> chi = m_matrice.polynomeCaracteristique();
		polynome<T> chi_ = chi;

		/*
		polynome<T> R;
		while ((R = PGCD(chi_, derivee(chi_))).degre >= 1) {
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
		T e = m_matrice.coeffs[0][0];
		e = true;
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
		iterateurs = { 0 };

		m_sous_ev_condition = NULL;
	};

	inline polynome<T> separer(polynome<T>& P, polynome<T>& Q) {
		polynome<T> R = PGCD(P, Q);
		if (R.degre <= 0)
			return R;
		P = P / R;
		Q = Q / R;
		polynome<T> R1;
		while ((R1 = PGCD(P, R)).degre >= 1)
			P = P / R1;

		return R;
	};

	bool ajouter(polynome<T> Q, int i) { //ajoute une condition Q ... tous les polynomes condition divisent chi ! par récurrence (on commence avec chi). 
//		if (aleatoire)
//			return (bool)Q; // à vérifier ... je crois OK
		if (Q.degre <= 0)
			return (bool)Q;  // Q constant ... soit nul soit jamais nul.

		polynome<T> poly_save = conditions[i];
		polynome<T> R = separer(Q, conditions[i]); //modifie la condition : référence

		if (R.degre <= 0) {
			return true; // l'élément testé n'est jamais nul ! pas de racine en commun avec la condition
		}
		else {
			if (conditions[i].degre <= 0) {
				//				conditions[i] = R * conditions[i];
				conditions[i] = poly_save; //pour diminuer l'erreur. cas approx
				return false; //l'élément testé est toujours nul dans ce cas : Q = R * qqchose, et R est la nouvelle condition
			}
			else {
				matrices.push_back(matrices[i]);
				conditions.push_back(R);
				iterateurs.push_back(iterateurs[i]);

				mettre_a_jour_quotient(matrices[matrices.size() - 1], R);
				mettre_a_jour_quotient(matrices[i], conditions[i]);

				return true; //l'élément testé n'est jamais nul dans ce cas : la nouvelle condition n'a aucune racine en commun avec l'élément : on a séparé.
				//Mais la condition ajoutée annule l'élément testé :  dans l'espace ajouté, Q sera nul (la condition ajoutée est R, et Q = R*quelquechose)
			}
		}
	};

	void calculer(int i) {
		//		matrice < anneau_quotient<polynome<T>>>& matrices[i] = matrices[i];
		polynome<T>& poly = conditions[i];
		int taille = matrices[i].taille_l;
		long question;
		auto vrai = matrices[i].coeffs[0][0];
		vrai = true;

		for (int j = iterateurs[i]; j < taille; ++j) { //colonne j ...
			iterateurs[i] = j;
			int k = j;
			int j_max = -1;
			auto norme_max = norme(vrai);
			for(int k(j);k<taille;++k) {
				if (ajouter(matrices[i].coeffs[k][j].element, i)) { //si élément nul ...
					auto norme_temp = norme(matrices[i].coeffs[k][j]);
					if (j_max == -1) {
						j_max = k;
						norme_max = norme_temp;
					}
					else {
						if (norme_temp > norme_max) {
							j_max = k;
							norme_max = norme_temp;
						}
					}
				}
			}

			if (j_max==-1) //la colonne j nulle sous la diagonale
				continue;

			if (j_max != j)
				matrices[i].echangerLigne(j, j_max);


			polynome<T> inverse_ = inverse(matrices[i].coeffs[j][j].element, conditions[i]); //condition changée

			for (int k(0); k < taille; ++k) {
				if (k == j)
					continue;
				matrices[i].ajouterLigne(j, k, -(inverse_ * matrices[i].coeffs[k][j]));
			}
		}

		//il y a deux parties : une avec une diagonale non-nulle, une autre triangulaire supérieur : 
		//on peut vérifier que comme en-dessous de la diagonale c'est nul, cette partie reste triangulaire supérieur.
		// 
		//une partie a des 0 sur la diagonale ... 
		//on cree la liste des i qui ont 0 sur la diagonale
		std::vector<int> non_diagonale(0);
		for (int j(0); j < taille; ++j)
			if (!(bool)matrices[i].coeffs[j][j]) //ces éléments ont été testés : toujours nuls ou jamais nuls ...
				non_diagonale.push_back(j); //si nul, on ajoute
		//non_diagonale est trié (ordre croissant), et triangulaire supérieur et diagonale nulle.

		/*
		//on cherche sur chaque ligne le premier élément non-nul ... on le teste bien-sûr
		for (int j(non_diagonale.size() - 2); j >= 0; --j)
			for (int k(j + 1); k < non_diagonale.size(); ++k)
				if (ajouter(matrices[i].coeffs[non_diagonale[j]][non_diagonale[k]].element, i)) { //si l'élément n'est jamais nul ... On calcule son inverse etc.
					polynome<T> inverse_ = inverse(matrices[i].coeffs[non_diagonale[j]][non_diagonale[k]].element, conditions[i]); //condition changée
					for (int l(0); l < non_diagonale.size() - 1; ++l) {
						if (l == j)
							continue;
						matrices[i].ajouterLigne(non_diagonale[j], non_diagonale[l], -(inverse_ * matrices[i].coeffs[non_diagonale[l]][non_diagonale[k]]));
						//si ce premier élément est non-nul, on supprime toute la colonne.
					}
					break;
				}
				*/

		//on cherche sur chaque ligne le premier élément non-nul ... on le teste bien-sûr
		std::vector<bool> ligne_faite(non_diagonale.size() - 1, false);
		for (int j(non_diagonale.size() - 2); j >= 0; --j) { //ligne j
			if (ligne_faite[j])
				continue;
			for (int k(j + 1); k < non_diagonale.size(); ++k) {//colonne k
				if (ajouter(matrices[i].coeffs[non_diagonale[j]][non_diagonale[k]].element, i)) { //si l'élément n'est jamais nul ... On calcule son inverse etc.
					int ligne_max = -1;
					auto norme_max = norme(vrai);

					for (int ligne(0); ligne < non_diagonale.size() - 2; ++ligne) {//on re-parcourt les lignes
						for (int colonne(0); colonne <= k; ++colonne) { //qui ont le premier élément non-nul à k
							if ((colonne < k) && ajouter(matrices[i].coeffs[non_diagonale[ligne]][non_diagonale[colonne]].element, i))
								break;
							if ((colonne == k) && ajouter(matrices[i].coeffs[non_diagonale[ligne]][non_diagonale[k]].element, i)) {
								auto norme_temp = norme(matrices[i].coeffs[non_diagonale[ligne]][non_diagonale[k]]);//et on maximise la norme
								if (ligne_max == -1) {
									ligne_max = ligne;
									norme_max = norme_temp;
								}
								else {
									if (norme_temp > norme_max) {
										ligne_max = ligne;
										norme_max = norme_temp;
									}
								}
							}
						}
					}
					//la norme a été trouvée pour la ligne : ligne_max. colonne k.

					polynome<T> inverse_ = inverse(matrices[i].coeffs[non_diagonale[ligne_max]][non_diagonale[k]].element, conditions[i]); //condition changée
					for (int l(0); l < non_diagonale.size() - 1; ++l) {
						if (l == j)
							continue;
						matrices[i].ajouterLigne(non_diagonale[ligne_max], non_diagonale[l], -(inverse_ * matrices[i].coeffs[non_diagonale[l]][non_diagonale[k]]));
						//si ce premier élément est non-nul, on supprime toute la colonne.
					}
					ligne_faite[ligne_max] = true; //on ne re-visitera pas cette ligne (premier élément est fixé).
					j = non_diagonale.size() - 1; //on recommence la visite des lignes, pour trouver un premier élément non-nul.
					break;//on quitte la recherche de k.
				}
				if(k==non_diagonale.size()-1) //la ligne est entierement nulle : ne pas recommencer.
					ligne_faite[j] = true;
			}
		}//jusqu'à temps d'avoir visité toutes les lignes de premier-élément.
		
		 //les premiers éléments de chaque ligne sont testés non-nuls, et la colonne restante est annulée.


		//on ajoute les colonnes nulles exactement dans la liste des parametres libres. (dans le sous-ensemble non-diagonale)
		std::vector<int> parametres_libres(0);
		for (int j = 0; j < non_diagonale.size(); ++j) {
			bool test = true;
			for (int k(0); k < non_diagonale.size(); ++k)
				if ((bool)matrices[i].coeffs[non_diagonale[k]][non_diagonale[j]]) {
					test = false;
					break;
				}
			if (test)
				parametres_libres.push_back(j);
		}

		//on ajoute les colonnes non-nulles (pas forcément exactement : sinon elles seraient dans la liste précédente)
		//à droite d'un premier nombre. (dans le sous-ensemble non-diagonale). Ce premier nombre est déjà testé : toujours non-nul.
		for(int j(0); j < non_diagonale.size(); ++j) {
			bool trouve = false;
			for (int k(j + 1); k < non_diagonale.size(); ++k)
				if (!trouve) { //trouvé le premier non-nul
					if ((bool)matrices[i].coeffs[non_diagonale[j]][non_diagonale[k]]) {//on regarde si c'est un premier élément non-nul
						trouve = true;
						continue;
					}
				}
				else {
					if ((bool)matrices[i].coeffs[non_diagonale[j]][non_diagonale[k]]) //on ajoute. (apres avoir trouvé un premier non-nul).
						parametres_libres.push_back(k);
				}
		}
		//en gros on a tout ajouté sauf le premier nombre non-nul de chaque ligne. dans le sous-espace non-diagonal. Les parametres premier-non-nul sont fixés, les autres sont libres.

		parametres_libres.erase(std::unique(parametres_libres.begin(), parametres_libres.end()), parametres_libres.end()); //supprimer les doublons. A vérifier.

		for (int vec_k(0); vec_k < parametres_libres.size(); ++vec_k) { //parametres libre c'est dans le sous-espace sans diagonale.
			anneau_quotient<polynome<T>> temp = matrices[i].coeffs[0][0];
			temp = false;

			std::vector<anneau_quotient<polynome<T>>> vec_non_diagonale(non_diagonale.size(), temp);
			temp = true;
			vec_non_diagonale[parametres_libres[vec_k]] = temp; //faux partout sauf dans ce parametre libre.

			//on regarde, pour chaque premier élément non-nul, sa valeur ... un seul parametre libre non-nul.
			for (int k(0); k < non_diagonale.size() - 1; ++k)
				for (int j(k + 1); j < non_diagonale.size(); ++j)
					if ((bool)matrices[i].coeffs[non_diagonale[k]][non_diagonale[j]]) { //un seul parametre libre vaut 1, les autres 0. Sur chaque ligne un seul parametre libre, donc ...
						vec_non_diagonale[j] = -(inverse(matrices[i].coeffs[non_diagonale[k]][non_diagonale[j]].element, conditions[i]) * matrices[i].coeffs[non_diagonale[k]][non_diagonale[parametres_libres[vec_k]]]);//condition changée
						break; // l'inverse existe : l'élément est testé et/donc toujours non-nul.
					}//dans les faits, ne compte que si l'élément parametre libre est à droite du premier élément non-nul.

			temp = false;
			//le vecteur entier ... on commence par recopier la partie déjà calculée.
			std::vector<anneau_quotient<polynome<T>>> vec(taille, temp);
			for (int j(0); j < non_diagonale.size(); ++j)
				vec[non_diagonale[j]] = vec_non_diagonale[j];

			for (int j(0); j < taille; ++j)
				if ((bool)matrices[i].coeffs[j][j]) {// les éléments de diagonale non-nulle : on calcule leur valeur.
					anneau_quotient<polynome<T>> somme = matrices[i].coeffs[0][0];
					somme = false;
					for (int k(0); k < non_diagonale.size(); ++k) {
						somme = somme + (matrices[i].coeffs[j][non_diagonale[k]] * vec_non_diagonale[k]);
					}
					vec[j] = -(inverse(matrices[i].coeffs[j][j].element, conditions[i]) * somme);
				}

			//on en extrait seulement les polynomes.
			std::vector<polynome<T>> vec_final(taille, temp.element);
			for (int j(0); j < taille; ++j)
				vec_final[j] = vec[j].element;


			m_sous_ev_condition->ajouter_vecteur(conditions[i], vec_final);
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

	std::vector<polynome<T>> conditions; //polynomes conditions
	std::vector<matrice<anneau_quotient<polynome<T>>>> matrices; //là où on en est dans le calcul
	std::vector<int> iterateurs; // on se trouve à cette colonne ...

	sous_ev_condition<T>* m_sous_ev_condition;
	polynome<T> chi;
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

//class double_erreur

//matrice<T> * vecteur<U>.
// => instaurer la multiplication par U à gauche comme à droite. pour T*T rien ne change.
//operator %=

