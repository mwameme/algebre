#pragma once

#include <exception>
#include <iostream>

#include "entete objets.hpp"
#include "types.hpp"
#include "swap_T.hpp"




template<class T> class anneau_quotient {
public:
    T element;
    T quotient;

    static_assert(type_algebre<T>::type == 1,"anneau_quotient<T> : T nécessite une division avec reste. Division exacte est futile");

    explicit anneau_quotient() :element(), quotient() {};


    explicit anneau_quotient(T element_, T quotient_) : element(element_), quotient(quotient_) {
        element %= quotient;
    };
    

    anneau_quotient(anneau_quotient<T> const& copie) : element(copie.element), quotient(copie.quotient){    };

    anneau_quotient(anneau_quotient<T>&& copie) {
        swap(*this, copie);
        return;
    };

    anneau_quotient<T>& operator=(anneau_quotient<T> const& temp);

    anneau_quotient<T>& operator=(anneau_quotient<T>&& temp) {
        if (this == &temp)
            return *this;
        swap(*this, temp);
        return *this;
    };

    anneau_quotient<T>& operator*=(anneau_quotient<T> const& droit) {
        element *= droit.element;
        element %= quotient;
        return *this;
    };

    template<class U>
    anneau_quotient<T>& operator*=(const U& scalaire) {
        element *= scalaire;
        element %= quotient;
        return *this;
    };

    anneau_quotient<T>& operator+=(anneau_quotient<T> const& droit) {
        element += droit.element;
        element %= quotient;
        return *this;
    };

    friend anneau_quotient<T> operator*(anneau_quotient<T> const& temp1, anneau_quotient<T> const& temp2) {
        anneau_quotient<T> result((temp1.element * temp2.element) % temp1.quotient, temp1.quotient);
        return result;
    };

    template<class U> friend anneau_quotient<T> operator*(U const& scalaire, anneau_quotient<T> const& temp) {
        anneau_quotient<T> result = temp;
        result.element = (scalaire * result.element);
        result.element %= result.quotient;
        return result;
    };

    friend anneau_quotient<T> operator+(anneau_quotient<T> const& temp1, anneau_quotient<T> const& temp2) {
        anneau_quotient<T> result((temp1.element + temp2.element) % temp1.quotient, temp1.quotient);
        return result;
    };

    friend anneau_quotient<T> operator-(anneau_quotient<T> const& temp1, anneau_quotient<T> const& temp2) {
        anneau_quotient<T> result((temp1.element - temp2.element) % temp1.quotient, temp1.quotient);
        return result;
    };

    friend anneau_quotient<T> operator%(anneau_quotient<T> const& temp1, anneau_quotient<T> const& temp2) {
        anneau_quotient<T> result((temp1.element % temp2.element) //% temp1.quotient
            , temp1.quotient);
        return result;
    };

    friend anneau_quotient<T> operator-(anneau_quotient<T> const& temp1) {
        anneau_quotient<T> result((-temp1.element) //% temp1.quotient
            , temp1.quotient);
        return result;
    };

    
    friend anneau_quotient<T> operator/(anneau_quotient<T> const& temp1, anneau_quotient<T> const& temp2) {
        anneau_quotient<T> result((temp1.element / temp2.element) // % temp1.quotient
            , temp1.quotient);
        return result;
    };

    explicit inline operator bool() const {
        return (bool) element;
    };

    template<class U> U operator()(const U& temp) const {
#ifdef ALGEBRA_USE_EXCEPTION
        if ((bool) quotient(temp))
            std::cerr << "ATTENTION : évaluation de polynome : quotient(element) != 0" << std::endl;
#endif
        return element(temp);
    };

    friend bool operator==(anneau_quotient<T> const& temp1, anneau_quotient<T> const& temp2) {
        if (temp1.quotient == temp2.quotient)
            return ! (bool) ( (temp1.element - temp2.element) % temp1.quotient);
        return (! (bool) ((temp1.element - temp2.element) % temp1.quotient)) && ( ! (bool) (temp1.quotient % temp2.quotient)) && (! (bool) (temp2.quotient % temp1.quotient)); 
        //dans ce dernier cas : vérifier que les 2 quotients sont équivalents ...
    };

    template<class U> explicit operator anneau_quotient<U>() const {
        anneau_quotient<U> result((U)element, (U)quotient);
        return result;
    };

    
    explicit operator T () const {
        return element;
    };


    void changer_quotient(T const& quotient_) {
        quotient = quotient_;
        element %= quotient;
    };

    friend std::ostream& operator<<(std::ostream& os, anneau_quotient<T> const& temp) {
        //        element.getDegre();
        os << "(" << temp.element << " % " << temp.quotient << ")";
        return os;
    };

    friend void swap(anneau_quotient<T>& gauche, anneau_quotient<T>& droit) {
        swap_F(gauche.element, droit.element);
        swap_F(gauche.quotient, droit.quotient);
        return;
    };

};

template<class T>
anneau_quotient<T>& anneau_quotient<T>::operator=(anneau_quotient<T> const& temp) {
    if (this == &temp)
        return *this;
    element = temp.element;
    quotient = temp.quotient;
    return *this;
};
