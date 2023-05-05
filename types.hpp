#pragma once
#include "entete objets.hpp"

#include <complex>


template<class T> class erreur_b;
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
template<class T> class polynome_n_sparse;

class InfInt;
class int_precision;
class float_precision;

template<class T> class type_algebre {
public:
//	using corps =; //le corps de base. vaut void pour les entiers ...

//	static constexpr int type = ; //0 : division exacte ; 1 : division avec reste ; 2 : sans divisio,
//	static constexpr int approx = ; //0 : exact ; 1 : approché 
//	static constexpr bool is_ratio = //dit si il y a un rationnel, en oubliant le rationnel<int> éventue de la fin
//	static constexpr bool is_objet = //dit si il y a au moins un objet ... (polynomes)
//	static constexpr bool is_corps_ratio = false; //dit si le corps de base est un rationnel d'entier.
//	objet : désigne un des types de polynome ...
};

template<class T> class type_algebre<polynome<T>> {
public:
	using corps = typename type_algebre<T>::corps;

	static constexpr int type = (type_algebre<T>::type == 0? 1:2);
	static constexpr int approx= type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = true;
	static constexpr bool is_corps_ratio = type_algebre<corps>::is_corps_ratio;
};

template<class T> class type_algebre<polynome_n_sparse<T>> {
public:
	using corps = typename type_algebre<T>::corps;

	static constexpr int type = 2;
	static constexpr int approx = type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = true;
	static constexpr bool is_corps_ratio = type_algebre<corps>::is_corps_ratio;
};


template<class T> class type_algebre<polynome_n_rec<T>> {
public:
	using corps = typename type_algebre<T>::corps;

	static constexpr int type = 2;
	static constexpr int approx = type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = true;
	static constexpr bool is_corps_ratio = type_algebre<corps>::is_corps_ratio;
};

template<class T> class type_algebre<polynome_n_iter<T>> {
public:
	using corps = typename type_algebre<T>::corps;

	static constexpr int type = 2;
	static constexpr int approx = type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = true;
	static constexpr bool is_corps_ratio = type_algebre<corps>::is_corps_ratio;

};

template<class T,int n> class type_algebre<polynome_n_fixe<T,n>> {
public:
	using corps = typename type_algebre<T>::corps;

	static constexpr int type = 2;
	static constexpr int approx = type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = true;
	static constexpr bool is_corps_ratio = type_algebre<corps>::is_corps_ratio;
};

template<class T> class type_algebre<polynome_n_fixe<T, 1>> {
public:
	using corps = typename type_algebre<T>::corps;

	static constexpr int type = (type_algebre<T>::type == 0 ? 1 : 2);
	static constexpr int approx = type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = true;
	static constexpr bool is_corps_ratio = type_algebre<corps>::is_corps_ratio;
};

template<class T> class type_algebre<polynome_n_fixe<T, 0>> {
public:
	using corps = typename type_algebre<T>::corps;

	static constexpr int type = type_algebre<T>::type;
	static constexpr int approx = type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = type_algebre<T>::is_objet;
	static constexpr bool is_corps_ratio = type_algebre<T>::is_corps_ratio;
};



template<class T> class type_algebre<anneau_quotient<T>> {
public:
	using corps = typename type_algebre<T>::corps;


	static constexpr int type = type_algebre<T>::type;
	static constexpr int approx = type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = type_algebre<T>::is_objet;
	static constexpr bool is_corps_ratio = type_algebre<corps>::is_corps_ratio;
};

template<class T> class type_algebre<corps_quotient<T>> {
public:
	using corps = typename std::conditional< std::is_same<typename type_algebre<T>::corps, void>::value, typename  corps_quotient<T>, typename type_algebre<T>::corps>::type;

	static constexpr bool get_is_corps_ratio() {
		if constexpr (std::is_same< corps_quotient<T>, corps >::value)
			return false;
		else
			return type_algebre<T>::is_corps_ratio;
	}

	static constexpr int type=0;
	static constexpr int approx= type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = type_algebre<T>::is_objet;
	static constexpr bool is_corps_ratio = get_is_corps_ratio();
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
	static constexpr bool get_is_corps_ratio() {
		if constexpr (std::is_same< rationnel<T>, corps >::value)
			return true;
		else
			return type_algebre<T>::is_corps_ratio;
	}

	static constexpr int type = 0;
	static constexpr int approx= type_algebre<T>::approx;
	static constexpr bool is_ratio = get_ratio();
	static constexpr bool is_objet = type_algebre<T>::is_objet;
	static constexpr bool is_corps_ratio = get_is_corps_ratio();
};

template<class T> class type_algebre<matrice<T>> {
public:
	using corps = typename type_algebre<T>::corps;


	static constexpr int type = 1; //faire attention aux inverses
	static constexpr int approx = type_algebre<T>::approx;
	static constexpr bool is_ratio = type_algebre<T>::is_ratio;
	static constexpr bool is_objet = type_algebre<T>::is_objet;
	static constexpr bool is_corps_ratio = type_algebre<corps>::is_corps_ratio;
};



template<class T> class type_algebre<erreur_b<T>> {
public:
	using corps = erreur_b<T>;

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
//	static constexpr bool is_corps_ratio = false;

};

template<> class type_algebre<long> {
public:
	using corps = void;

	static constexpr int type = 1;
	static constexpr int approx = 0;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
//	static constexpr bool is_corps_ratio = false;

};

template<> class type_algebre<long long> {
public:
	using corps = void;

	static constexpr int type = 1;
	static constexpr int approx = 0;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
//	static constexpr bool is_corps_ratio = false;

};



template<> class type_algebre<int_precision> {
public:
	using corps = void;

	static constexpr int type = 1;
	static constexpr int approx = 0;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
//	static constexpr bool is_corps_ratio = false;

};

template<> class type_algebre<InfInt> {
public:
	using corps = void;

	static constexpr int type = 1;
	static constexpr int approx = 0;
	static constexpr bool is_ratio = false;
	static constexpr bool is_objet = false;
//	static constexpr bool is_corps_ratio = false;

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
