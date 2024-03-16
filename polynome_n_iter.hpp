#pragma once
#include "vecteur_n.hpp"
#include <vector>
#include <string>
#include "polynome_n_rec.hpp"
#include <exception>
#include <initializer_list>

#include "entete objets.hpp"
#include <iostream>
#include "polynome_n_sparse.hpp"

template<class T> class vecteur_n;

template<class T> class polynome_n_rec;

template<class T> polynome_n_rec<T> poly_n_convert_rec(const T* data, const int* dimensions, const int* puissances, int n, std::string* noms);

template<class T> class polynome_n_sparse;
template<class T> class monome;


template<class T> class polynome_n_iter {
public:

	bool scalaire; //true si c'est juste un scalaire
	vecteur_n<T> coeffs;

	using sous_type = typename T;

	polynome_n_iter() {};

	polynome_n_iter(int n_var, T element) {//polynome de degre 0. utilisé pour unite
#ifdef _DEBUG
		if (n_var < 0)
			throw std::domain_error("polynome_n_iter : n_var < 0.");
#endif
		if (n_var == 0) {
			std::vector<int> dimensions(1, 1);
			coeffs = vecteur_n<T>(dimensions);
			coeffs.data[0] = element;
			scalaire = true;
		}
		else {
			std::vector<int> dimensions(n_var, 1);
			coeffs = vecteur_n<T>(dimensions);
			coeffs.data[0] = element;
			scalaire = false;
		}
	};
	polynome_n_iter(monome<T> monome_) : polynome_n_iter(monome_.degres, monome_.element, true) {};
	polynome_n_iter(std::vector<int> degres, T element, bool is_degre = true) {
		//monome. degres >0. Transforme degres en dimensions ...  
		//si true ATTENTION créé avec degrés (+1) ... Sinon avec dimensions ...

		if (degres.size() == 0) {
			degres = { 1 };
			coeffs = vecteur_n<T>(degres);
			coeffs.data[0] = element;
			scalaire = true;
		}
		else {
			if (is_degre) {
#ifdef _DEBUG
				for (int i(0); i < degres.size(); ++i)
					if (degres[i] < 0) {
						//						degres = std::vector<int>(degres.size(), 0);
						//						element = unite(element, false);
						//						break;
						throw std::domain_error("creation de polynome_n_iter : degre negatif");
					}
#endif
				for (int i(0); i < degres.size(); ++i)
					degres[i] += 1; //dimensions
			}
#ifdef _DEBUG
			for (int i(0); i < degres.size(); ++i)
				if (degres[i] <= 0)
					throw std::domain_error("declaration de polynome_n_iter avec un vecteur de dimensions : une dimension <= 0");
#endif
			coeffs = vecteur_n<T>(degres);
			T faux_ = unite(element, false);

			for (int i(0); i < coeffs.data.size(); ++i)
				coeffs.data[i] = faux_;
			coeffs.data[coeffs.data.size() - 1] = element; //monome
			scalaire = false;
		}
	};

	polynome_n_iter(const vecteur_n<T>& tableau) { //au cas où ?
		coeffs = tableau;
		scalaire = false;
	};

	polynome_n_iter(vecteur_n<T>&& tableau) { //au cas où ?
		swap(coeffs, tableau);
		scalaire = false;
	};


	polynome_n_iter(polynome_n_iter const& temp) : coeffs(temp.coeffs), scalaire(temp.scalaire) {
	};

	polynome_n_iter(polynome_n_iter&& temp) {
		swap(*this, temp);
	};

	polynome_n_iter<T>& operator=(polynome_n_iter<T> const& temp) {
		if (this == &temp)
			return *this;
		coeffs = temp.coeffs;
		scalaire = temp.scalaire;
		return *this;
	};

	polynome_n_iter<T>& operator=(polynome_n_iter<T>&& temp) {
		if (this == &temp)
			return *this;
		swap(*this, temp);
		return *this;
	};

	polynome_n_iter<T>& operator*=(polynome_n_iter<T> const& temp) {
		return (*this = (*this * temp));
	};

	template<class U>
	polynome_n_iter<T>& operator*=(U const& scalaire_) {
		if ((bool)scalaire_)
			for (int i(0); i < coeffs.data.size(); ++i)
				coeffs.data[i] *= scalaire_;
		else
			*this = polynome_n_iter(coeffs.puissance, unite(coeffs.data[0], false));

		return *this;
	};

	polynome_n_iter<T>& operator+=(polynome_n_iter<T> const& autre) {
#ifdef _DEBUG
		if (coeffs.puissance != autre.puissance)
			throw std::domain_error("addition de polynomes : puissances ne correspondent pas");
#endif
		std::vector<int> dim_finales(coeffs.puissance, 0);
		for (int i(0); i < coeffs.puissance; ++i)
			dim_finales[i] = max(coeffs.dimensions[i], autre.coeffs.dimensions[i]);

		coeffs.modifier_dimensions(dim_finales);

		int n = coeffs.puissance;
		int position_finale = 0;
		std::vector<int> puissances_finale = coeffs.puissances;
		std::vector<int> puissances = autre.coeffs.puissances;
		std::vector<int> dimensions = autre.coeffs.dimensions;
		bool fin = true;
		int position = 0;
		std::vector<int> positions(n, 0);

		while (fin) {
			coeffs.data[position_finale] += autre.coeffs.data[position];
			++position;
			++positions[n - 1];
			++position_finale;
			int i = n - 1;

			while (positions[i] >= dimensions[i]) { //On réajuste autre, et donc finale aussi.
				position_finale = position_finale - puissances_finale[i] * positions[i];
				positions[i] = 0;
				if (i == 0) {
					fin = false;
					break;
				}
				++positions[i - 1];
				position_finale = position_finale + puissances_finale[i - 1];
				--i;
			}
		}
		return *this;
	};

	friend polynome_n_iter<T> operator+(polynome_n_iter<T> const& gauche, polynome_n_iter<T> const& droite) {
#ifdef _DEBUG
		if (gauche.coeffs.puissance != droite.coeffs.puissance)
			throw std::domain_error("polynome_n_iter, addition : le nombre de variables ne correspond pas.");

		if (((gauche.scalaire) && (!droite.scalaire)) || ((droite.scalaire) && (!gauche.scalaire)))
			throw std::domain_error("addition de polynome_n_iter : scalaire + polynome");
#endif
		if (gauche.scalaire)
			return polynome_n_iter(0, gauche.coeffs.data[0] + droite.coeffs.data[0]);


		int n = gauche.coeffs.puissance;
		std::vector<int> deg(n, 0);
		for (int i(0); i < n; ++i)
			deg[i] = max(gauche.coeffs.dimensions[i], droite.coeffs.dimensions[i]) - 1;

		T faux = unite(gauche.coeffs.data[0], false);

		polynome_n_iter<T> result(deg, faux); //ATTENTION degrés et non dimensions.

		int position_finale = 0;
		std::vector<int> puissances_finale = result.coeffs.puissances;
		//gauche
		std::vector<int> puissances = gauche.coeffs.puissances;
		std::vector<int> dimensions = gauche.coeffs.dimensions;
		bool fin = true;
		int position = 0;
		std::vector<int> positions(n, 0);
		while (fin) {
			result.coeffs.data[position_finale] = gauche.coeffs.data[position];
			++position;
			++positions[n - 1];
			++position_finale;
			int i = n - 1;

			while (positions[i] >= dimensions[i]) { //On réajuste gauche, et donc finale aussi.
				position_finale = position_finale - puissances_finale[i] * positions[i];
				positions[i] = 0;
				if (i == 0) {
					fin = false;
					break;
				}
				++positions[i - 1];
				position_finale = position_finale + puissances_finale[i - 1];
				--i;
			}
		}
		//droite
		position = 0;
		positions = std::vector<int>(n, 0);
		puissances = droite.coeffs.puissances;
		fin = true;
		position_finale = 0;
		dimensions = droite.coeffs.dimensions;

		while (fin) {
			result.coeffs.data[position_finale] += droite.coeffs.data[position];
			++position;
			++positions[n - 1];
			++position_finale;
			int i = n - 1;

			while (positions[i] >= dimensions[i]) { //On réajuste droit, et donc finale aussi.
				position_finale = position_finale - puissances_finale[i] * positions[i];
				positions[i] = 0;
				if (i == 0) {
					fin = false;
					break;
				}
				++positions[i - 1];
				position_finale = position_finale + puissances_finale[i - 1];
				--i;
			}
		}

		return result;
	};

	friend polynome_n_iter<T> operator-(polynome_n_iter<T> const& gauche, polynome_n_iter<T> const& droite) {
#ifdef _DEBUG
		if (gauche.coeffs.puissance != droite.coeffs.puissance)
			throw std::domain_error("polynome_n_iter, addition : le nombre de variables ne correspond pas.");

		if (((gauche.scalaire) && (!droite.scalaire)) || ((droite.scalaire) && (!gauche.scalaire)))
			throw std::domain_error("addition de polynome_n_iter : scalaire + polynome");
#endif
		if (gauche.scalaire)
			return polynome_n_iter(0, gauche.coeffs.data[0] - droite.coeffs.data[0]);



		int n = gauche.coeffs.puissance;
		std::vector<int> deg(n, 0);
		for (int i(0); i < n; ++i)
			deg[i] = max(gauche.coeffs.dimensions[i], droite.coeffs.dimensions[i]) - 1;

		T faux = unite(gauche.coeffs.data[0], false);

		polynome_n_iter<T> result(deg, faux); //ATTENTION degrés et non dimensions.

		int position_finale = 0;
		std::vector<int> puissances_finale = result.coeffs.puissances;
		//gauche
		std::vector<int> puissances = gauche.coeffs.puissances;
		std::vector<int> dimensions = gauche.coeffs.dimensions;
		bool fin = true;
		int position = 0;
		std::vector<int> positions(n, 0);
		while (fin) {
			result.coeffs.data[position_finale] = gauche.coeffs.data[position];
			++position;
			++positions[n - 1];
			++position_finale;
			int i = n - 1;

			while (positions[i] >= dimensions[i]) { //On réajuste gauche, et donc finale aussi.
				position_finale = position_finale - puissances_finale[i] * positions[i];
				positions[i] = 0;
				if (i == 0) {
					fin = false;
					break;
				}
				++positions[i - 1];
				position_finale = position_finale + puissances_finale[i - 1];
				--i;
			}
		}
		//droite
		position = 0;
		positions = std::vector<int>(n, 0);
		puissances = droite.coeffs.puissances;
		fin = true;
		position_finale = 0;
		dimensions = droite.coeffs.dimensions;

		while (fin) {
			result.coeffs.data[position_finale] += - droite.coeffs.data[position];
			++position;
			++positions[n - 1];
			++position_finale;
			int i = n - 1;

			while (positions[i] >= dimensions[i]) { //On réajuste droit, et donc finale aussi.
				position_finale = position_finale - puissances_finale[i] * positions[i];
				positions[i] = 0;
				if (i == 0) {
					fin = false;
					break;
				}
				++positions[i - 1];
				position_finale = position_finale + puissances_finale[i - 1];
				--i;
			}
		}

		return result;
	}


	friend polynome_n_iter<T> operator-(polynome_n_iter<T> const& temp) {
		polynome_n_iter<T> result = temp;
		for (int i(0); i < result.coeffs.data.size(); ++i)
			result.coeffs.data[i] = -result.coeffs.data[i];
		return result;
	};

	friend polynome_n_iter<T> operator*(polynome_n_iter<T> const& gauche, polynome_n_iter<T> const& droite) { // mettre le polynome le plus sparse à droite ! optimisé.
#ifdef _DEBUG
		if (gauche.coeffs.puissance != droite.coeffs.puissance)
			throw std::domain_error("polynome_n_iter, multiplication : le nombre de variables ne correspond pas.");
		if (((gauche.scalaire) && (!droite.scalaire)) || ((droite.scalaire) && (!gauche.scalaire)))
			throw std::domain_error("multiplication de polynome_n_iter : scalaire * polynome");
#endif
		if (gauche.scalaire)
			return polynome_n_iter(0, gauche.coeffs.data[0] * droite.coeffs.data[0]);


		//mettre la plus sparse à droite
		int n_sparse = 0;
		for (int i(0); i < gauche.coeffs.data.size(); ++i)
			n_sparse += (bool)gauche.coeffs.data[i];

		int m_sparse = 0;
		for (int i(0); i < droite.coeffs.data.size(); ++i)
			m_sparse += (bool)droite.coeffs.data[i];

		if (n_sparse == 0)
			return polynome_n_iter<T>(gauche.coeffs.puissance, unite(gauche.coeffs.data[0], false));
		if (m_sparse == 0)
			return polynome_n_iter<T>(droite.coeffs.puissance, unite(droite.coeffs.data[0], false));


		if ((n_sparse * droite.coeffs.data.size()) < (m_sparse * gauche.coeffs.data.size()))
			return droite * gauche;

		int n = gauche.coeffs.puissance;
		std::vector<int> degres(n);
		for (int i(0); i < n; ++i)
			degres[i] = gauche.coeffs.dimensions[i] + droite.coeffs.dimensions[i] - 2;

		T faux_ = unite(gauche.coeffs.data[0], false);

		polynome_n_iter<T> result(degres, faux_); //ATTENTION degrés et non dimensions.

		std::vector<int> dimensions_gauche = gauche.coeffs.dimensions;
		std::vector<int> dimensions_droite = droite.coeffs.dimensions;

		std::vector<int> positions_gauche(n, 0);
		std::vector<int> positions_droite(n, 0);
		int position_gauche = 0;
		int position_droite = 0;
		int position_finale = 0;

		std::vector<int> puissances_finale = result.coeffs.puissances;
		std::vector<int> puissances_gauche = gauche.coeffs.puissances;
		std::vector<int> puissances_droite = droite.coeffs.puissances;

		bool fin = true;
		bool fin_gauche = true;
		int i;

		if (!(bool)droite.coeffs.data[0])
			goto increment_droite;

		while (fin) {
			//on opère la multiplication ... chaque couple (gauche , droite) apparait une et une seule fois.
			if (((bool)gauche.coeffs.data[position_gauche]) && ((bool)droite.coeffs.data[position_droite]))
				result.coeffs.data[position_finale] += (gauche.coeffs.data[position_gauche] * droite.coeffs.data[position_droite]);
			++position_gauche;
			++positions_gauche[n - 1];
			++position_finale;
			i = n - 1;
			fin_gauche = true;

			while (positions_gauche[i] >= dimensions_gauche[i]) { //On réajuste gauche, et donc finale aussi.
				position_finale = position_finale - puissances_finale[i] * positions_gauche[i];
				//				position_gauche = position_gauche - puissances_gauche[i] * positions_gauche[i];

				positions_gauche[i] = 0;
				if (i == 0) {
					fin_gauche = false;
					position_gauche = 0;
					break;
				}
				++positions_gauche[i - 1];
				//				position_gauche = position_gauche + puissances_gauche[i - 1];
				position_finale = position_finale + puissances_finale[i - 1];
				--i;
			}

			fin = true;
			if (!fin_gauche) { //gauche a dépassé, et vaut  0 partout ... on incrémente à droite.
			increment_droite:
				++position_droite;
				++positions_droite[n - 1];
				++position_finale;
				i = n - 1;

				while (positions_droite[i] >= dimensions_droite[i]) {
					//					position_finale = position_finale - puissances_finale[i] * positions_droite[i];
					//					position_droite = position_droite - puissances_droite[i] * positions_droite[i];
					position_finale = position_finale - puissances_finale[i] * positions_droite[i];

					positions_droite[i] = 0;
					if (i == 0) {
						fin = false;
						break;
					}
					++positions_droite[i - 1];
					//					position_droite = position_droite + puissances_droite[i - 1];
					position_finale = position_finale + puissances_finale[i - 1];
					--i;
				}
				if (!fin)
					break;

				if (!(bool)droite.coeffs.data[position_droite])
					goto increment_droite;
			}
			//fin du while
		}//Le calcule marche, car le vecteur finale est la sommme du vecteur gauche et droite, et cela se calcule avec les puissances (finales). 

		return result;

	};

	explicit inline operator bool() const {
		if (scalaire)
			return (bool)coeffs.data[0];

		if ((bool)coeffs.data[coeffs.data.size() - 1])
			return true;

		for (int i(0); i < coeffs.data.size(); ++i)
			if ((bool)coeffs.data[i])
				return true;
		//		std::vector<int> vec(coeffs.puissance, 1);
		//		coeffs.modifier_dimensions(vec);
		return false;
	};


	explicit inline operator bool() {
		if (scalaire)
			return (bool)coeffs.data[0];

		if ((bool)coeffs.data[coeffs.data.size() - 1])
			return true;

		for (int i(0); i < coeffs.data.size(); ++i)
			if ((bool)coeffs.data[i])
				return true;
		std::vector<int> vec(coeffs.puissance, 1);
		coeffs.modifier_dimensions(vec);
		return false;
	};


	operator polynome_n_rec<T>() const {
		if (scalaire)
			return polynome_n_rec<T>(0, coeffs.data[0]);
		polynome_n_rec<T> result = poly_n_convert_rec(coeffs.data.data(), coeffs.dimensions.data(), coeffs.puissances.data(), coeffs.puissance);
		return result;
	};

	friend std::ostream& operator<<(std::ostream& os, polynome_n_iter<T> const& element) {
		bool premier = false;
		for (n_for iter(element.coeffs.dimensions); (bool)iter; ++iter)
			if ((bool)element.coeffs.data[iter.position]) {
				std::string puissance = "";
				for (int i(0); i < element.coeffs.puissance; ++i)
					puissance = puissance + "* " + "X" + std::to_string(i) + "^" + std::to_string(iter.positions[i]) + " ";
				if (premier)
					os << "+" << element.coeffs.data[iter.position] << " " << puissance;
				else {
					premier = true;
					os << element.coeffs.data[iter.position] << " " << puissance;
				}

			}
		return os;
	};

	polynome_n_iter<T> operator() (int i, polynome_n_iter<T> const& Y) const {
#ifdef _DEBUG
		if ((i >= coeffs.puissance) || (i < 0))
			throw std::domain_error("evaluation de polynome_n_iter : i non-conforme");

		if (Y.coeffs.puissance != coeffs.puissance - 1)
			throw std::domain_error("evaluation de polynome : objet évalué non-conforme");

		if (Y.scalaire)
			throw std::domain_error("evaluation de polynome_n_iter : polynome=scalaire");
#endif
		if (coeffs.puissance == 1) {
			T result = unite(coeffs.data[0], false);
			T puissance = unite(result, true);
			T y = Y.coeffs.data[0];

			for (int i(0); i < coeffs.data.size(); ++i) {
				result = result + (coeffs.data[i] * puissance);
				puissance = puissance * y;
			}

			return polynome_n_iter(0, result);
		}

		int n = coeffs.puissance;
		int position = 0;
		bool fin = true;
		std::vector<int> positions = std::vector<int>(n, 0);

		//preparer les puissances de Y
		std::vector<polynome_n_iter<T>> pow_vector(coeffs.dimensions[i]);
		pow_vector[0] = unite(Y, true);
		for (int i(1); i < coeffs.dimensions[i]; ++i)
			pow_vector[i] = Y * pow_vector[i - 1];

		polynome_n_iter<T> result(n - 1, unite(coeffs.data[0], false), Y.noms); //polynome nul.
		while (fin) {

			int j = n - 1;
			if (j == i)
				--j;
			if (j < 0) //non-essentiel. Pour être sûr.
				break;
			do {
				++positions[j];
				position += coeffs.puissances[j];
				if (positions[j] < coeffs.dimensions[j])
					break;
				positions -= positions[j] * coeffs.puissances[j];
				--j;
				if (j == i)
					--j;
				if (j < 0)
					break;
			} while (true);
			if (j < 0)
				break;
			//on a la nouvelle position ... correcte avec 0 sur positions[i].

			std::vector<int> nouvelles_positions = positions;
			nouvelles_positions.erase(nouvelles_positions.begin() + i); //Pour le polynome à n-1 variables.
			polynome_n_iter<T> poly_fixe(nouvelles_positions, unite(coeffs.data[0], true)); //degrés et nouvelles_positions correspondent ... OK.

			//			polynome_n_iter<T> pow = unite(Y, true);
			for (int k(0); k < coeffs.dimensions[i]; ++k) {
				result = result + (coeffs.data[position + k * coeffs.puissances[i]] * (pow_vector[k] * poly_fixe));
				//				pow = pow * Y;
			}
		};

		return result;
	};

	template<class U>
	U operator() (int i, U const& Y) const { //Y doit être de type matrice de polynome_n_iter, ou rationnel de polynome_n_iter, pour que ça ait du sens ... 
		//Plus précisément polynome de taille n-1
		// Est utile pour évaluer une fraction de polynome_n_iter ... On évalue le numerateur puis le denominateur ...
#ifdef _DEBUG
		if ((i >= coeffs.puissance) || (i < 0))
			throw std::domain_error("evaluation de polynome_n_iter : i non-conforme");

		if (Y.coeffs.puissance != coeffs.puissance - 1)
			throw std::domain_error("evaluation de polynome : objet évalué non-conforme");

		if (Y.scalaire)
			throw std::domain_error("evaluation de polynome_n_iter : polynome=scalaire");
#endif
		if (coeffs.puissance == 1) {
			U result = unite(Y, false);
			U puissance = unite(Y, true);

			for (int i(0); i < coeffs.data.size(); ++i) {
				result = result + (coeffs.data[i] * puissance);
				puissance = puissance * Y;
			}

			return result;
		}

		int n = coeffs.puissance;
		int position = 0;
		bool fin = true;
		std::vector<int> positions = std::vector<int>(n, 0);

		//preparer les puissances de Y
		std::vector<U> pow_vector(coeffs.dimensions[i]);
		pow_vector[0] = unite(Y, true);
		for (int i(1); i < coeffs.dimensions[i]; ++i)
			pow_vector[i] = Y * pow_vector[i - 1];

		U result = unite(Y, false); //polynome nul.
		while (fin) {
			int j = n - 1;
			if (j == i)
				--j;
			if (j < 0) //non-essentiel. Pour être sûr.
				break;
			do {
				++positions[j];
				position += coeffs.puissances[j];
				if (positions[j] < coeffs.dimensions[j])
					break;
				position -= positions[j] * coeffs.puissances[j];
				--j;
				if (j == i)
					--j;
				if (j < 0)
					break;
			} while (true);
			if (j < 0)
				break;
			//on a la nouvelle position ... correcte avec 0 sur positions[i].

			std::vector<int> nouvelles_positions = positions;
			nouvelles_positions.erase(nouvelles_positions.begin() + i); //Pour le polynome à n-1 variables.
			polynome_n_iter<T> poly_fixe(nouvelles_positions, unite(coeffs.data[0], true), NULL); //degrés et nouvelles_positions correspondent ... OK.

			//			polynome_n_iter<T> pow = unite(Y, true);
			for (int k(0); k < coeffs.dimensions[i]; ++k) {
				result = result + (coeffs.data[position + k * coeffs.puissances[i]] * (poly_fixe * pow_vector[k]));
				//				pow = pow * Y;
			}
		};

		return result;
	};

	polynome_n_iter<T> operator() (int i, T element) const {
#ifdef _DEBUG
		if (coeffs.puissance <= 0)
			throw std::domain_error("evaluation de polynome_n_iter : puissance requise >0");
#endif
		polynome_n_iter<T> Y(coeffs.puissance - 1, element);
		return *this(i, Y);
	};

	template<class U>
	U operator() (std::vector<U> eval) const {
#ifdef _DEBUG
		if (eval.size() != coeffs.puissance)
			throw std::domain_error("evaluation de polynome_n_iter : les dimensions ne coincident pas");
#endif
		U neutre = unite(eval[0], true);
		std::vector<std::vector<U>> puissances(eval.size(), std::vector<U>(0));
		for (int i(0); i < eval.size(); ++i) {
			puissances[i] = std::vector<U>(coeffs.dimensions[i], neutre);
			for (int j(1); j < coeffs.dimensions[i]; ++j)
				puissances[i][j] = puissances[i][j - 1] * eval[i];
		}

		U resultat = unite(neutre, false);

		for (typename vecteur_n<T>::const_iterator it = coeffs.cbegin(); (bool) it; ++it) {
			if ((bool) *it)
			{
				U temp = neutre;
				for (int i(0); i < eval.size(); ++i) {
					temp *= puissances[i][it.mPositions[i]];
				}
				resultat += (*it) * temp;
			}
		}
		return resultat;
	}

	template<class U> friend polynome_n_iter<T> operator*(U const& scalaire, polynome_n_iter<T> const& poly) {
		polynome_n_iter<T> result(poly);
		for (int i(0); i < result.coeffs.data.size(); ++i)
			result.coeffs.data[i] *= result.coeffs.data[i];
		return result;
	};

	void simplifier() {
		if (!scalaire)
			coeffs.simplifier();
		return;
	};

	void simplifier_2() { //mieux
		*this = ((polynome_n_iter<T>) ((polynome_n_rec<T>) * this)); //deux conversions ... la premiere simplifie.
	};

	friend void swap(polynome_n_iter<T>& gauche, polynome_n_iter<T>& droit) {
		std::swap(gauche.scalaire, droit.scalaire);
		swap(gauche.coeffs, droit.coeffs);
		return;
	};

	friend bool operator==(polynome_n_iter<T> const& gauche, polynome_n_iter<T> const& droit) {
		if (gauche.scalaire != droit.scalaire)
			return false;
		return ((polynome_n_rec<T>) gauche) == ((polynome_n_rec<T>) droit); //plus simple ... optimise (en +)

	};

	using iterator = typename vecteur_n<T>::iterator;
	using const_iterator = typename vecteur_n<T>::const_iterator;

	iterator begin() {
		return coeffs.begin();
	};
	iterator end() {
		return coeffs.end();
	};

	const_iterator cbegin() const {
		return coeffs.cbegin();
	};
	const_iterator cend() const {
		return coeffs.cend();
	};


	operator polynome_n_sparse<T>() const {
		polynome_n_sparse<T> poly(monome<T>(std::vector<int>(coeffs.puissance, 0), unite(coeffs.data[0], false)));
		for (int i(0); i < coeffs.data.size(); ++i)
			if ((bool)coeffs.data[i])
				poly.ajouter(monome<T>(coeffs.positions(i), coeffs.data[i]));
		poly.simplifier();
		return poly;
	};

	template<class U>
	operator polynome_n_iter<U>() const { //degre element noms false
		if (scalaire)
			return polynome_n_iter<U>(std::vector<int>(0), (U)coeffs.data[0]);
		polynome_n_iter<U> poly(coeffs.dimensions, unite((U)coeffs.data[0], false), false); //false car ce sont les dimensions, et non les degrés ...
		for (int i(0); i < poly.coeffs.data.size(); ++i)
			poly.coeffs.data[i] = (U)coeffs.data[i];
		return poly;
	};

	std::vector<int> max_degre() const {
		std::vector<int> result(coeffs.puissance, -1);
		for (int i(0); i < coeffs.data.size(); ++i)
			if ((bool)coeffs.data[i]) {
				std::vector<int> positions = coeffs.positions(i);
				for (int j(0); j < coeffs.puissance; ++j)
					if (positions[j] > result[j])
						result[j] = positions[j];
			}
		return result;
	};

	int max_degre(int deg) const {
		int result = -1;
		int pow = coeffs.puissances[deg];
		int dim = coeffs.dimensions[deg];
		for (int i(0); i < coeffs.data.size(); ++i)
			if ((bool)coeffs.data[i]) {
				int pos = (i / pow) % dim;
				if (pos > result) {
					result = pos;
					if (result == (dim - 1))
						break;
				}
			}
		return result;
	};

	polynome_n_iter<T> derivee(int i) const {
		T nul = unite(get_T(), false);

		std::vector<int> nouvelles_dimensions = coeffs.dimensions;
		nouvelles_dimensions[i] -= 1;
		if (nouvelles_dimensions[i] <= 0)
			return polynome_n_iter<T>(coeffs.n_var, nul);

		polynome_n_iter<T> resultat = *this;

		for (const_iterator it = resultat.cbegin(); (bool)it; ++it) {
			int pos = it.mPosition;
			int pos2 = pos + coeffs.puissances[i];
			if (pos2 < coeffs.data.size())
				resultat.coeffs.data[pos] = (it.mPositions[i] +1) * coeffs.data[pos2];
			else
				resultat.coeffs.data[pos] = nul;
		}

		resultat.coeffs.modifier_dimension(nouvelles_dimensions);
		return resultat;
	}

	void deriver(int i) {
		int n = coeffs.puissance;
		int position = 0;
		bool fin = true;
		std::vector<int> positions = std::vector<int>(n, 0);
		T faux = unite(get_T(), false);

		while (fin) {
			int j = n - 1;
			if (j == i)
				--j;
			if (j < 0) //non-essentiel. Pour être sûr.
				break;
			do {
				++positions[j];
				position += coeffs.puissances[j];
				if (positions[j] < coeffs.dimensions[j])
					break;
				position -= positions[j] * coeffs.puissances[j];
				--j;
				if (j == i)
					--j;
				if (j < 0)
					break;
			} while (true);

			if (j < 0)
				break;
			//on a la nouvelle position ... correcte avec 0 sur positions[i].
			for (int k(1); k < coeffs.dimensions[i]; ++k) {
				int pos2 = position + coeffs.puissances[i] * k;
				if (pos2 >= coeffs.data.size()) {
					coeffs.data[position + (k - 1) * coeffs.puissances[i]] = faux;
					break;
				}
				else
					coeffs.data[position + (k - 1) * coeffs.puissances[i]] = k * coeffs.data[position + k * coeffs.puissances[i]];
			}
		}

		simplifier();
		return;

	}

	T get_T() {
		return coeffs.data[0];
	}

};

template<class T> polynome_n_rec<T> poly_n_convert_rec(const T* data, const int* dimensions, const int* puissances, int n) {
	if (n == 0)
		return polynome_n_rec<T>(*data);

	std::vector<polynome_n_rec<T>> tab(0);
	tab.reserve(*dimensions);
	for (int i(0); i < *dimensions; ++i)
		tab.push_back(poly_n_convert_rec(data + (i * puissances[0]), dimensions + 1, puissances + 1, n - 1));

	T faux = unite(*data, false);
	return polynome_n_rec<T>(n, faux, tab); //simplifie aussi.
};
