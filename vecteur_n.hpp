#pragma once
#include <vector>
#include "fonctions template.hpp"
#include "n_for.hpp"
#include <exception>

#include "unite.hpp"

template<class T>  T unite(T const& element, bool test);

template<class T> class vecteur_n {
public:

	std::vector<int> puissances;
	int puissance; //nombre de variables ... data[1][2]...[puissance]
	std::vector<int> dimensions;
	std::vector<T> data;

	vecteur_n() {};

	vecteur_n(std::vector<int> const& liste_dim) {
#ifdef ALGEBRA_USE_EXCEPTION
		if (liste_dim.size() == 0)
			throw std::domain_error("vecteur_n : liste vide");
#endif
		dimensions = liste_dim;
		puissance = dimensions.size();
#ifdef ALGEBRA_USE_EXCEPTION
		for (int i(0); i < puissance; ++i)
			if (dimensions[i] < 0)
				throw std::domain_error("vecteur_n : une dimension < 0");
#endif
		puissances= std::vector<int>(puissance);
		puissances[puissance - 1] = 1;
		for (int i(puissance-2); i >= 0; --i)
			puissances[i]= puissances[i+1]* dimensions[i+1];
		int taille = puissances[0] * dimensions[0];
		data = std::vector<T>(taille);
	};

	vecteur_n(std::vector<int> const& liste_dim,T const& element) {
#ifdef ALGEBRA_USE_EXCEPTION
		if (liste_dim.size() == 0)
			throw std::domain_error("vecteur_n : liste vide");
#endif
		dimensions = liste_dim;
		puissance = dimensions.size();
#ifdef ALGEBRA_USE_EXCEPTION
		for (int i(0); i < puissance; ++i)
			if (dimensions[i] < 0)
				throw std::domain_error("vecteur_n : une dimension < 0");
#endif
		puissances = std::vector<int>(puissance);
		puissances[puissance - 1] = 1;
		for (int i(puissance - 2); i >= 0; --i)
			puissances[i] = puissances[i + 1] * dimensions[i + 1];
		int taille = puissances[0] * dimensions[0];
		data = std::vector<T>(taille,element);
	};

	vecteur_n(vecteur_n<T> const& temp) {
		puissance = temp.puissance;
		dimensions = temp.dimensions;
		data = temp.data;
		puissances = temp.puissances;
	};

	vecteur_n(vecteur_n<T> && temp) {
		swap(*this, temp);
	};


	
	vecteur_n<T>& operator=(vecteur_n<T> const& temp) {
		if (this == &temp)
			return *this;
		puissance = temp.puissance;
		dimensions = temp.dimensions;
		data = temp.data;
		puissances = temp.puissances;
		return *this;
	};

	vecteur_n<T>& operator=(vecteur_n<T> && temp) {
		if (this == &temp)
			return *this;
		swap(*this, temp);
		return *this;
	};


	/*
	vecteur_n<T>& operator=(bool test) {
		for (int i(0); i < data.size(); ++i)
			data[i] = test;
		return *this;
	};
	*/

	int position(std::vector<int> const& positions_) const { //position doit être de la bonne taille (et du bon format)
		int expo = 1;
		int pos = 0;
		for (int i(puissance - 1); i >= 0; --i) {
			pos = pos + positions_[i] * puissances[i];
		}
		return pos;
	};

	std::vector<int> positions(int pos_) const{
		std::vector<int> retour(puissance, 0);

		for (int i(puissance - 1); i >= 0; --i) {
			retour[i] = pos_ % dimensions[i];
			pos_ = pos_ / dimensions[i];
		}

		return retour;
	};

	bool accesseur_possible(std::vector<int> const& positions_) const {
		if (positions_.size() != puissance)
			return false;
		for (int i(0); i < puissance; ++i)
			if (positions_[i] >= dimensions[i])
				return false;
		return true;
	};

	inline T& operator[](std::vector<int> positions_)  { //
		return data[position(positions_)];
	};

	inline T operator[](std::vector<int> positions_) const { //
		return data[position(positions_)];
	};

	inline T& operator[](int position) {
		return data[position];
	};

	inline T operator[](int position) const {
		return data[position];
	};

	inline void modifier_dimension(std::vector<int> const& nouvelles_dimensions) {
		T element = unite(data[0],false);
		modifier_dimension(nouvelles_dimensions, element);
	};

	void modifier_dimension(std::vector<int> const& nouvelles_dimensions, T element) { //même puissance
		bool test=true;
		for(int i(0);i<puissance;++i)
			if (nouvelles_dimensions[i] != dimensions[i]) {
				test = false;
				break;
			}
		if (test)
			return;

		//faire la boucle de dimension2, et recopier (si possible). boucle intelligente.
		std::vector<int> pos(puissance, 0);
		int pos_int = 0;
		int pos_int2 = 0;

		std::vector<int> puissances2(puissance);
		puissances2[puissance - 1] = 1;

		for (int i(puissance - 2); i >= 0; --i) {
			puissances2[i] = puissances2[i + 1] * nouvelles_dimensions[i + 1];
		}
		int taille = puissances2[0] * nouvelles_dimensions[0];
		std::vector<T> data2(taille, element); //on met des 0 quand on peut ... en gros.


		while (true) {
			data2[pos_int2] = data[pos_int];
			++pos_int;
			++pos_int2;
			++pos[puissance - 1];

			int i = puissance - 1;
			while (pos[i] >= min(nouvelles_dimensions[i], dimensions[i])) {
				if (i == 0){
					i = -1;
					break;
				};

				pos_int = pos_int - pos[i] * puissances[i] + puissances[i - 1];
				pos_int2 = pos_int2 - pos[i] * puissances2[i] + puissances2[i - 1];
				pos[i] = 0;
				++pos[i - 1];
				--i;
			};

			if (i == -1)
				break;
		};

		swap(data , data2);
		dimensions = nouvelles_dimensions;
		swap(puissances , puissances2);

		return;
	};

	void simplifier() {
		if ((bool)data[data.size() - 1])
			return;
		std::vector<int> dimensions_max(puissance, 0);

		for (n_for iter(dimensions); (bool) iter; ++iter)
			if ((bool)data[iter.position])
				for (int i(0); i < puissance; ++i)
					if (iter.positions[i] > dimensions_max[i])
						dimensions_max[i] = iter.positions[i];

		bool modifie = false; //on passe à la taille .... dernier element +1
		for (int i(0); i < puissance; ++i) {
			++dimensions_max[i];
			if (dimensions_max[i] != dimensions[i])
				modifie = true;
		}

		if (modifie)
			modifier_dimension(dimensions_max);
	};

	friend void swap(vecteur_n<T>& gauche, vecteur_n<T>& droit) {
		std::swap(gauche.puissance, droit.puissance);
		swap(gauche.data, droit.data);
		swap(gauche.puissances, droit.puissances);
		swap(gauche.dimensions, droit.dimensions);
		return;
	}

	template<class I>
	class iterator_vecteur_n;

	using iterator = iterator_vecteur_n<T>;
	using const_iterator = iterator_vecteur_n<const T>;

	iterator begin() {
		return iterator(*this);
	};

	iterator end() {
		iterator it(*this);
		it.go_position(data.size());
		return it;
	};

	const_iterator cbegin() const {
		return const_iterator(*this);
	};

	const_iterator cend() const {
		const_iterator it(*this);
		it.go_position(data.size());
		return it;
	};

};



template<class T> template<class I>
class vecteur_n<T>::iterator_vecteur_n {
public:
	using value_type = std::remove_const_t<I>;
	using vecteur_type = std::conditional_t< std::is_const_v<I>, const vecteur_n<value_type>, vecteur_n<value_type>	>;

	vecteur_type* pointeur;
	int mPosition;
	std::vector<int> mPositions;

	iterator_vecteur_n& operator++() {
		++mPosition;
		for (int i(pointeur->puissance - 1); i >= 0; --i) {
			++mPositions[i];
			if (mPositions[i] >= pointeur->dimensions[i])
				mPositions[i] = 0;
			else
				break;
		}
		return *this;
	};

	operator bool() const {
		return mPosition < pointeur->data.size();
	};

	I& operator*() {
		return pointeur->data[mPosition];
	};

	std::vector<int> positions() {
		return mPositions;
	};

	int position() {
		return mPosition;
	}

	void go_position(int i) {
		mPosition = i;
		mPositions = pointeur->positions(i);
	};

	iterator_vecteur_n<T>& operator+=(int i) {
		mPosition += i;
		mPositions = pointeur->positions(mPosition);
		return *this;
	}

	void go_position(std::vector<int> pos) {
		mPositions = pos;
		mPosition = pointeur->position(pos);
	};

	friend bool operator==(iterator_vecteur_n const& gauche, iterator_vecteur_n const& droit) {
		return ((gauche.pointeur == droit.pointeur) && (gauche.mPosition == droit.mPosition));
	};

	friend bool operator!=(iterator_vecteur_n const& gauche, iterator_vecteur_n const& droit) {
		return ((gauche.pointeur != droit.pointeur) || (gauche.mPosition != droit.mPosition));
	};

	iterator_vecteur_n(vecteur_type& vec) {
		pointeur = &vec;
		mPosition = 0;
		mPositions = std::vector<int>(puissance, 0);
	};

	iterator_vecteur_n& operator=(iterator_vecteur_n const& copie) {
		pointeur = copie.pointeur;
		mPosition = copie.mPosition;
		mPositions = copie.mPositions;
		return *this;
	}
};
//verifier >= dimension.
//simplifier ...