#pragma once
#include <initializer_list>
#include <vector> 


template<class T> class convertir_T {
public:
	static T convertir(const T element) {
		return element;
	};
};

template<class T> class convertir_T<std::initializer_list<T>> {
public:
	static std::vector< decltype(convertir_T<T>::convertir(T()))> convertir(const std::initializer_list<T> liste) {
		std::vector< decltype(convertir_T<T>::convertir(T()))> vec(0);
		for (auto it = liste.begin(); it != liste.end(); ++it)
			vec.push_back(convertir_T<T>::convertir(*it));

		return vec;
	};
};

template<class T> inline decltype(convertir_T<T>::convertir(T())) convertir(T x) {
	return convertir_T<T>::convertir(x);
};


//template<class T> inline decltype(convertir_T<T>::convertir(std::declval<T>())) convertir(const T x) {

/*
template<class T> std::vector< decltype(convertir_liste(T()))> convertir_liste(std::initializer_list<T> liste) {
	std::vector< decltype(convertir_liste(T()))> vec(0);
	for (auto it = liste.begin(); it != liste.end(); ++it)
		vec.push_back(convertir_liste(*it));

	return vec;
}
*/