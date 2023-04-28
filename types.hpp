#pragma once
#include "entete objets.hpp"

#include <complex>


template<class T> class erreur;
template<class T> class erreur_l;
template<class T> class anneau_quotient;
template<class T> class complexe;
template<class T> class corps_quotient;
template<class T, class enable1 = void, class enable2 = void> class matrice;
template<class T> class polynome;
template<class T> class polynome_n_rec;
template<class T> class polynome_n_iter;
template<class T> class rationnel;
template<class T, int n> class polynome_n_fixe;

class InfInt;
class int_precision;
class float_precision;

template<class T> class type_algebre {
public:
	using corps = void; //le corps de base. En supposant que c'est bien défini (sinon ça bug).

	static constexpr int type=2; //0 : division exacte ; 1 : division avec reste ; 2 : sans divisio,
	static constexpr int approx=0; //0 : exact ; 1 : approché 
	static constexpr bool is_ratio = false; //dit si il y a un rationnel, en oubliant le rationnel<int> éventue de la fin
	static constexpr bool is_objet = true; //dit si il y a au moins un objet ... (polynomes)
	//	static constexpr bool is_corps_ratio = false; //dit si le corps de base est un rationnel d'entier. mal défini pour les autres types.
};

template<class T> class type_algebre<polynome<T>> {
public:
	using corps = typename type_algebre<T>::corps;

	static constexpr int type = (type_algebre<T>::type == 0? 1:2);
	static constexpr int approx= type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = true;
};

template<class T> class type_algebre<polynome_n_rec<T>> {
public:
	using corps = typename type_algebre<T>::corps;

	static constexpr int type = 2;
	static constexpr int approx = type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = true;
};

template<class T> class type_algebre<polynome_n_iter<T>> {
public:
	using corps = typename type_algebre<T>::corps;

	static constexpr int type = 2;
	static constexpr int approx = type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = true;
};

template<class T,int n> class type_algebre<polynome_n_fixe<T,n>> {
public:
	using corps = typename type_algebre<T>::corps;

	static constexpr int type = 2;
	static constexpr int approx = type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = true;
};



template<class T> class type_algebre<anneau_quotient<T>> {
public:
	using corps = typename type_algebre<T>::corps;


	static constexpr int type = type_algebre<T>::type;
	static constexpr int approx = type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = type_algebre<T>::is_objet;
};

template<class T> class type_algebre<corps_quotient<T>> {
public:
	using corps = typename std::conditional< std::is_same<typename type_algebre<T>::corps, void>::value, typename  corps_quotient<T>, typename type_algebre<T>::corps>::type;


	static constexpr int type=0;
	static constexpr int approx= type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = type_algebre<T>::is_objet;
	static constexpr bool is_corps_ratio = false;
};

/*
template<class T> class type_algebre<rationnel<T>> {
public:
	using corps = typename type_algebre<T>::corps;

	static constexpr int get_type() {
		return 0;
	};

	static constexpr int get_approx() {
		return type_algebre<T>::approx;
	};

	static constexpr int type{ get_type() };
	static constexpr int approx{ get_approx() };
};
*/

template<class T> class type_algebre<rationnel<T>> {
public:
	using corps = typename std::conditional< std::is_same<typename type_algebre<T>::corps,void>::value ,typename  rationnel<T> , typename type_algebre<T>::corps>::type;


	static constexpr bool get_ratio() {
		if constexpr (std::is_same<corps, rationnel<T>>::value)
			return false;
		return true;
	};


	static constexpr int type = 0;
	static constexpr int approx= type_algebre<T>::approx;
	static constexpr bool is_ratio = get_ratio();
	static constexpr bool is_objet = type_algebre<T>::is_objet;
	static constexpr bool is_corps_ratio = true;
};


template<class T> class type_algebre<erreur<T>> {
public:
	using corps = erreur<T>;

	static constexpr int type = 0;
	static constexpr int approx = 1;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
	static constexpr bool is_corps_ratio = false;
};

template<class T> class type_algebre<erreur_l<T>> {
public:
	using corps = erreur_l<T>;

	static constexpr int type = 0;
	static constexpr int approx = 1;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
	static constexpr bool is_corps_ratio = false;
};


template<class T> class type_algebre<complexe<T>> {
public:
	using corps = typename complexe<T>;

	static constexpr int type = 0;
	static constexpr int approx = type_algebre<T>::approx;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
	static constexpr bool is_corps_ratio = type_algebre<T>::is_corps_ratio;
};

template<class T> class type_algebre<std::complex<T>> {
public:
	using corps = typename std::complex<T>;

	static constexpr int type = 0;
	static constexpr int approx = type_algebre<T>::approx;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
	static constexpr bool is_corps_ratio = type_algebre<T>::is_corps_ratio;
};


// ==========================================================
//types de base
template<> class type_algebre<int> {
public:
	using corps = void;

	static constexpr int type = 1;
	static constexpr int approx = 0;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
};

template<> class type_algebre<long> {
public:
	using corps = void;

	static constexpr int type = 1;
	static constexpr int approx = 0;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
};

template<> class type_algebre<long long> {
public:
	using corps = void;

	static constexpr int type = 1;
	static constexpr int approx = 0;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
};



template<> class type_algebre<int_precision> {
public:
	using corps = void;

	static constexpr int type = 1;
	static constexpr int approx = 0;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
};

template<> class type_algebre<InfInt> {
public:
	using corps = void;

	static constexpr int type = 1;
	static constexpr int approx = 0;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
};

template<> class type_algebre<float> {
public:
	using corps = float;

	static constexpr int type = 0;
	static constexpr int approx = 1;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
	static constexpr bool is_corps_ratio = false;
};

template<> class type_algebre<double> {
public:
	using corps = double;

	static constexpr int type = 0;
	static constexpr int approx = 1;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
	static constexpr bool is_corps_ratio = false;
};
template<> class type_algebre<float_precision> {
public:
	using corps = float_precision;

	static constexpr int type = 0;
	static constexpr int approx = 1;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
	static constexpr bool is_corps_ratio = false;
};

/*
// ==========================================================
//redéfinit correctement corps dans ces cas ...
template<> class type_algebre<rationnel<int>> {
public:
	using corps = rationnel<int>;

	static constexpr int type = 0;
	static constexpr int approx = 0;
};



template<> class type_algebre<rationnel<long>> {
public:
	using corps = rationnel<long>;

	static constexpr int type = 0;
	static constexpr int approx = 0;
};


template<> class type_algebre<rationnel<long long>> {
public:
	using corps = rationnel<long long>;

	static constexpr int type = 0;
	static constexpr int approx = 0;
};


template<> class type_algebre<rationnel<int_precision>> {
public:
	using corps = rationnel<int_precision>;

	static constexpr int type = 0;
	static constexpr int approx = 0;
};


template<> class type_algebre<rationnel<InfInt>> {
public:
	using corps = rationnel<InfInt>;

	static constexpr int type = 0;
	static constexpr int approx = 0;
};
*/


/*

//0 : corps. 1 : anneau avec division et reste. 2 : anneau sans division
constexpr int type_algebre(int element) {
	return 1;
};

constexpr int type_algebre(long element) {
	return 1;
};

constexpr int type_algebre(long long element) {
	return 1;
};


template<class T> constexpr int type_algebre(InfInt& element) {
	return 1;
};

template<class T> constexpr int type_algebre(int_precision& element) {
	return 1;
};

constexpr int type_algebre(float element) {
	return 0;
};

constexpr int type_algebre(double element) {
	return 0;
};

constexpr int type_algebre(float_precision& element) {
	return 0;
};



template<class T> constexpr int type_algebre(polynome<T>& element) {
	if (type_algebre(T()) == 0)
		return 1;
	else
		return 2;
};

template<class T> constexpr int type_algebre(rationnel<T>& element) {
	return 0;
};


template<class T> constexpr int type_algebre(anneau_quotient<T>& element) {
	return (type_algebre(T()));
};


template<class T> constexpr int type_algebre(corps_quotient<T>& element) {
	if (type_algebre(T()) == 0)
		return 0;
	else
		if (type_alebre(T()) == 1)
			return 0;
	return 2; //interdit
};


template<class T> constexpr int type_algebre(matrice<T>& element) {
	if (type_algebre(T()) == 0)
		return 0;
	else
		return 2;
};

template<class T> constexpr int type_algebre(polynome_n_rec<T>& element) {
	return 2;
};

template<class T> constexpr int type_algebre(complexe<T>& element) {
	if (type_algebre(T()) == 0)
		return 0;
	else
		return 2;
};

template<class T> constexpr int type_algebre(erreur<T>& element) {
	return type_algebre(T());
};



//0 : exact. 1 : approché
constexpr int type_approx(int element) {
	return 0;
};

constexpr int type_approx(long element) {
	return 0;
};

constexpr int type_approx(long long element) {
	return 0;
};

constexpr int type_approx(InfInt& element) {
	return 0;
};

constexpr int type_approx(int_precision& element) {
	return 0;
};

constexpr int type_approx(float element) {
	return 1;
};

constexpr int type_approx(double element) {
	return 1;
};

template<class T> constexpr int type_approx(float_precision& element) {
	return 1;
};

// composés
template<class T> constexpr int type_approx(erreur<T>& element) {
	return type_approx(T());
};

template<class T> constexpr int type_approx(polynome<T>& element) {
	return type_approx(T());
};

template<class T> constexpr int type_approx(matrice<T>& element) {
	return type_approx(T());
};

template<class T> constexpr int type_approx(rationnel<T>& element) {
	return type_approx(T());
};


template<class T> constexpr int type_approx(anneau_quotient<T>& element) {
	return type_approx(T());
};



template<class T> constexpr int type_approx(corps_quotient<T>& element) {
	return type_approx(T());
};


template<class T> constexpr int type_approx(polynome_n_rec<T>& element) {
	return type_approx(T());
};

template<class T> constexpr int type_approx(complexe<T>& element) {
	return type_approx(T());
};


*/