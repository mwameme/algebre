#pragma once

#include <exception>
#include "types.hpp"

#include "entete objets.hpp"
#include "unite.hpp"
#include "swap_T.hpp"


template<class T> class corps_quotient {
public:
    static_assert(type_algebre<T>::type == 1, "corps_quotient<T> : T nécessite une division avec reste. Division exacte est futile");

    T element;
    T quotient;

    explicit corps_quotient() :element(),quotient(){};

    explicit corps_quotient(T element_, T quotient_) : element(element_), quotient(quotient_) {
        element %= quotient;
    };

    corps_quotient(corps_quotient<T> const& copie) : element(copie.element), quotient(copie.quotient) {    };

    corps_quotient(corps_quotient<T>&& copie)   {
        swap(*this, copie);
        return;
    };

    corps_quotient<T>& operator=(corps_quotient<T> const& temp);

    corps_quotient<T>& operator=(corps_quotient<T>&& temp)   {
        if (this == &temp)
            return *this;
        swap(*this, temp);
        return *this;
    };

    corps_quotient<T>& operator*=(corps_quotient<T> const& droit) {
        element *= droit.element;
        element %= quotient;
        return *this;
    };

    template<class U>
    corps_quotient<T>& operator*=(const U& scalaire) {
        element *= scalaire;
        element %= quotient;
        return *this;
    };

    corps_quotient<T>& operator+=(corps_quotient<T> const& autre) {
        element += autre.element;
        element %= quotient;
        return *this;
    };


    friend corps_quotient<T> operator*(corps_quotient<T> const& temp1, corps_quotient<T> const& temp2) {
        corps_quotient<T> result((temp1.element * temp2.element) % temp1.quotient, temp1.quotient);
        return result;
    };

    template<class U>
    friend corps_quotient<T> operator*(U const& scalaire, corps_quotient<T> const& temp) {
        corps_quotient<T> result = temp;
        result.element = scalaire * element;
        result.element %= result.quotient;
        return result;
    };

    friend corps_quotient<T> operator+(corps_quotient<T> const& temp1, corps_quotient<T> const& temp2) {
        corps_quotient<T> result((temp1.element + temp2.element) % temp1.quotient, temp1.quotient);
        return result;
    };

    friend corps_quotient<T> operator-(corps_quotient<T> const& temp1, corps_quotient<T> const& temp2) {
        corps_quotient<T> result((temp1.element - temp2.element) % temp1.quotient, temp1.quotient);
        return result;
    };

    friend corps_quotient<T> operator%(corps_quotient<T> const& temp1, corps_quotient<T> const& temp2) {
        corps_quotient<T> result((temp1.element % temp2.element) //% temp1.quotient
            , temp1.quotient);
        return result;
    };

    friend corps_quotient<T> operator-(corps_quotient<T> const& temp1) {
        corps_quotient<T> result((-temp1.element) //% temp1.quotient
            , temp1.quotient);
        return result;
    };

    friend corps_quotient<T> operator/(corps_quotient<T> const& temp1, corps_quotient<T> const& temp2) { //BEZOUT

        T p = temp1.quotient;
        T a, b, ub, vb, wb, c, q;

        a = p;
        b = temp2.element;
        if (!(bool)b)
            throw std::domain_error("corps_quotient : division par 0");
        vb = unite(b,true);
        ub = unite(b,false);

        c = vb;
        while ( (bool) c)
        {
            c = a % b;
            q = a / b;
//            wa = ua - q * va;
            wb = ub - q * vb;
            a = b;
            b = c;
//            ua = va;
            ub = vb;
//            va = wa;
            vb = wb;
        }
        ub = ub / a; //cas des polynomes.

        return corps_quotient<T>(temp1.element * ub, p);
    };


    explicit inline operator bool() const {
        return ((bool) element);
    };

    template<class U> U operator()(const U& temp) const {
#ifdef _DEBUG
        if ((bool)quotient(temp))
            throw std::domain_error("ATTENTION : évaluation de polynome : quotient(element) != 0");
#endif
        return element(temp);
    };

    friend bool operator==(corps_quotient<T> const& temp1, corps_quotient<T> const& temp2) {
        if (temp1.quotient == temp2.quotient)
            return ! (bool)((temp1.element - temp2.element) % temp1.quotient);
        return (! (bool)((temp1.element - temp2.element) % temp1.quotient)) && (! (bool)(temp1.quotient % temp2.quotient)) && !((bool)(temp2.quotient % temp1.quotient));
    };

    template<class U> explicit operator corps_quotient<U>() const {
        corps_quotient<U> result((U)element, (U)quotient);
        return result;
    };

    explicit operator T () const {
        return element;
    };

    void changer_quotient(T const& quotient_);

    friend std::ostream& operator<<(std::ostream& os, corps_quotient<T> const& temp) {
        //        element.getDegre();
        os << "(" << temp.element << " % " << temp.quotient << ")";
        return os;
    };

    friend void swap(corps_quotient<T>& gauche, corps_quotient<T>& droit) {
        swap_F(gauche.element, droit.element);
        swap_F(gauche.quotient, droit.quotient);
        return;
    };

};

template<class T>
corps_quotient<T>& corps_quotient<T>::operator=(corps_quotient<T> const& temp) {
    if (this == &temp)
        return *this;
    element = temp.element;
    quotient = temp.quotient;
    return *this;
};

template<class T>
void corps_quotient<T>::changer_quotient(T const& quotient_) {
    quotient = quotient_;
    element %= quotient;
};
