#pragma once

#include <exception>
#include <iostream>

#include "entete objets.hpp"
#include "types.hpp"


template<class T> class anneau_quotient {
public:

//    static_assert(type_algebre(T()) == 1,"anneau_quotient : sans division");
    static_assert(type_algebre<T>::type == 1,"anneau_quotient<T> : T nécessite une division avec reste. Division exacte est futile");

    explicit anneau_quotient() :element(), quotient() {};


    explicit anneau_quotient(T element_, T quotient_) : element(element_), quotient(quotient_) {
        element %= quotient;
    };
    

    anneau_quotient(anneau_quotient<T> const& copie) {
        element = copie.element;
        quotient = copie.quotient;
    };

    /*
    anneau_quotient<T>& operator=(bool test) {
        element = test;
        return *this;
    };
    */

    anneau_quotient<T>& operator=(anneau_quotient<T> const& temp) {
        if (this == &temp)
            return *this;
        element = temp.element;
        quotient = temp.quotient;
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

    /*
    friend std::ostream& operator<<(std::ostream& os, const anneau_quotient<T>& temp) {
        os << temp.element << " mod[" << temp.quotient << "] ";
        return os;
    };
    */

    template<class U> U operator()(const U& temp) const {
        if ((bool) quotient(temp))
            std::cerr << "ATTENTION : évaluation de polynome : quotient(element) != 0" << std::endl;
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
        result.element = result.element % result.quotient;
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

    T element;
    T quotient;
};

