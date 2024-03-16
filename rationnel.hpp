//#ifndef RATIONNEL_H
//#define RATIONNEL_H

#pragma once
#include <vector>
#include <iostream>
#include <exception>
#include <initializer_list>

#include "entete objets.hpp"
#include "swap_T.hpp"


template<class T> T unite(T const& element, bool test);

//template<class T> inline void swap_F(T& x, T& y);


template<class T> class rationnel {
public:
    T numerateur;
    T denominateur;

    rationnel<T> derivee(int i) const {
        rationnel<T> result;
        T x = denominateur.derivee(i);
        if ((bool)denominateur.derivee(i)) {
            result.numerateur = (numerateur.derivee(i) * denominateur) - (numerateur * x);
            result.denominateur = denominateur * denominateur;
            result.simplifier();
            return result;
        }
        else {
            result.numerateur = numerateur.derivee(i);
            result.denominateur = denominateur;
            return result;
        }
    }


    friend rationnel<T> derivee(rationnel<T> const& element) {
        rationnel<T> result;
        result.numerateur = (derivee(element.numerateur) * element.denominateur) - (element.numerateur * derivee(element.denominateur));
        result.denominateur = element.denominateur * element.denominateur;
        result.simplifier();
        return result;
    };

    explicit rationnel() : numerateur(), denominateur() {};

    explicit rationnel(T const& a, T const& b) : numerateur(a), denominateur(b) {
        simplifier();
    };

    explicit rationnel(T const& a) : numerateur(a), denominateur(unite(a, true)) {  };


    rationnel(const rationnel<T>& copie) : numerateur(copie.numerateur), denominateur(copie.denominateur) {    };

    rationnel(rationnel<T>&& temp) : numerateur(),denominateur()  {
        swap(*this, temp);
        return;
    };

    template<class U> explicit rationnel(std::vector<U> const& liste) {
        numerateur = T(liste[0]);
        if (liste.size() > 1)
            denominateur = T(liste[1]);
        else
            denominateur = unite(numerateur, true);
        simplifier();
    };

    rationnel<T>& operator=(const rationnel<T>& temp);

    rationnel<T>& operator=(rationnel<T>&& temp)  {
        if (this == &temp)
            return *this;
        swap(*this, temp);
        return *this;
    };

    rationnel<T>& operator*=(const rationnel<T>& temp) {
        numerateur *= temp.numerateur;
        denominateur *= temp.denominateur;
        simplifier();
        return *this;
    };

    template<class U>
    rationnel<T>& operator*=(U const& scalaire) {
        numerateur *= scalaire;
        simplifier();
        return *this;
    };

    rationnel<T>& operator+=(rationnel<T> const& temp) {
        numerateur = ((numerateur * temp.denominateur) + (denominateur * temp.numerateur));
        denominateur = denominateur * temp.denominateur;
        simplifier();
        return *this;
    };


    void simplifier();

    friend rationnel<T> operator* (const rationnel<T>& temp1, const rationnel<T>& temp2) {
        rationnel<T> result;
        result.numerateur = (temp1.numerateur * temp2.numerateur);
        result.denominateur = (temp1.denominateur * temp2.denominateur);
        result.simplifier();
        return result;
    };



    template<class U> friend rationnel<T> operator*(const U& scalaire, const rationnel<T>& temp) {
        //sert pour descendre les niveaux de construction. utilisé pour la dérivée : multiplier par un entier.
        rationnel<T> result;
        result.numerateur = scalaire * temp.numerateur;
        result.denominateur = temp.denominateur;
        result.simplifier();
        return result;
    };

    template<class U> U operator()(const U& element) const {
        U result = denominateur(element);
#ifdef _DEBUG
        if (!(bool)result)
            throw std::domain_error("division par zero lors d'une évaluation de fraction polynome");
#endif
        return (numerateur(element) / result);
    };

    friend rationnel<T> operator/(const rationnel<T>& temp1, const rationnel<T>& temp2) {
        rationnel<T> result;
        result.numerateur = (temp1.numerateur * temp2.denominateur);
        result.denominateur = (temp1.denominateur * temp2.numerateur);

        result.simplifier();

        return result;
    };

    friend rationnel<T> operator+(const rationnel<T>& temp1, const rationnel<T>& temp2) {
        rationnel<T> result;
        result.numerateur = ((temp1.numerateur * temp2.denominateur) + (temp2.numerateur * temp1.denominateur));
        result.denominateur = (temp1.denominateur * temp2.denominateur);

        result.simplifier();
        return result;
    };

    friend rationnel<T> operator-(const rationnel<T>& temp1, const rationnel<T>& temp2) {
        rationnel<T> result;
        result.numerateur = (temp1.numerateur * temp2.denominateur) - (temp2.numerateur * temp1.denominateur);
        result.denominateur = temp1.denominateur * temp2.denominateur;
        result.simplifier();
        return result;
    };

    friend rationnel<T> operator-(const rationnel<T>& temp) {
        rationnel<T> result;

        result.numerateur = -temp.numerateur;
        result.denominateur = temp.denominateur;
        return result;
    };

    friend bool operator==(rationnel<T> const& a, rationnel<T> const& b) {
        return (!((bool)((a.numerateur * b.denominateur) - (b.numerateur * a.denominateur))));
    };


    explicit inline operator bool() const {
        return ((bool)numerateur);
    };

    template<class U> explicit operator U() const {
        U result = (U)numerateur / (U)denominateur;
        return result;
    };


    template<class U> explicit operator rationnel<U>()const {
        rationnel<U> result((U)numerateur, (U)denominateur);
        return result;
    };



    friend std::ostream& operator<<(std::ostream& os, const rationnel<T>& element) {
        os << "(" << element.numerateur << "/" << element.denominateur << ")";
        return os;
    };

    friend void swap(rationnel<T>& gauche, rationnel<T>& droit) {
        swap_F(gauche.numerateur, droit.numerateur);
        swap_F(gauche.denominateur, droit.denominateur);
        return;
    };

    friend bool operator<(rationnel<T>& gauche, rationnel<T>& droit) { //seulement pour T = entier.
        static_assert(std::is_same<type_algebre<T>::corps, void>::value);

        bool test = true;
        if (gauche.denominateur < 0)
            test = !test;
        if (droit.denominateur < 0)
            test = !test;
        if (test)
            return gauche.numerateur * droit.denominateur < droit.numerateur * gauche.denominateur;
        else
            return !(gauche.numerateur * droit.denominateur < droit.numerateur * gauche.denominateur);
    };

    friend bool operator>(rationnel<T>& gauche, rationnel<T>& droit) { //seulement pour T = entier.
        static_assert(std::is_same<type_algebre<T>::corps, void>::value);

        bool test = true;
        if (gauche.denominateur < 0)
            test = !test;
        if (droit.denominateur < 0)
            test = !test;
        if (test)
            return gauche.numerateur * droit.denominateur > droit.numerateur * gauche.denominateur;
        else
            return !(gauche.numerateur * droit.denominateur > droit.numerateur * gauche.denominateur);
    };


};

template<class T>
rationnel<T>& rationnel<T>::operator=(const rationnel<T>& temp) {
    if (this == &temp)
        return *this;
    numerateur = temp.numerateur;
    denominateur = temp.denominateur;
    return *this;
};

template<class T>
void rationnel<T>::simplifier() {
    if (!(bool)numerateur) {
        denominateur = unite(denominateur, true);
        return;
    }

    if constexpr (type_algebre<T>::type == 1) {
        T pgcd = PGCD(numerateur, denominateur);

        if ((bool)pgcd) {
            numerateur = numerateur / pgcd;
            denominateur = denominateur / pgcd;
        }
        if (!(bool)numerateur)
            if ((bool)denominateur)
                denominateur = unite(denominateur, true);
        if (!(bool)denominateur)
            if ((bool)numerateur)
                numerateur = unite(numerateur, true);

        return;
    }
    return;
};