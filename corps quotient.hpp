#pragma once

#include <exception>
#include "types.hpp"

#include "entete objets.hpp"
#include "vrai et faux.hpp"


template<class T> class corps_quotient {
public:
    static_assert(type_algebre<T>::type == 1, "corps_quotient<T> : T nécessite une division avec reste. Division exacte est futile");


    explicit corps_quotient() :element(),quotient(){};

    explicit corps_quotient(T element_, T quotient_) : element(element_), quotient(quotient_) {
        element %= quotient;
    };

    corps_quotient(corps_quotient<T> const& copie) {
        element = copie.element;
        quotient = copie.quotient;
    };

    corps_quotient<T>& operator=(bool test) {
        element = test;
        return *this;
    };

    corps_quotient<T>& operator=(corps_quotient<T> const& temp) {
        element = temp.element;
        quotient = temp.quotient;
        return *this;
    };

    friend corps_quotient<T> operator*(corps_quotient<T> const& temp1, corps_quotient<T> const& temp2) {
        corps_quotient<T> result((temp1.element * temp2.element) % temp1.quotient, temp1.quotient);
        return result;
    };

    template<class U> friend corps_quotient<T> operator*(U const& scalaire, corps_quotient<T> const& temp) {
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

    friend corps_quotient<T> operator/(corps_quotient<T> const& temp1, corps_quotient<T> const& temp2) {

        T p = temp1.quotient;
        T a, b, //ua, va, wa, 
            ub, vb, wb, c, q;

        a = p;
        b = temp2.element;
//        ua = b; ua = true;
        vb = vrai(b);
        ub = faux(b);
//        va = ub;
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

        return corps_quotient<T>((temp1.element * ub) % p, p);
    };


    explicit operator bool() const {
        return ((bool) element);
    };

    /*
    friend std::ostream& operator<<(std::ostream& os, const corps_quotient<T>& temp) {
        os << temp.element << " mod[" << temp.quotient << "] ";
        return os;
    };
    */

    template<class U> U operator()(const U& temp) const {
        if ((bool)quotient(temp))
            std::cerr << "ATTENTION : évaluation de polynome : quotient(element) != 0" << std::endl;
        return element(temp);
    };

    friend bool operator==(corps_quotient<T> const& temp1, corps_quotient<T> const& temp2) {
        if (temp1.quotient == temp2.quotient)
            return ! (bool)((temp1.element - temp2.element) % temp1.quotient);
        return (! (bool)((temp1.element - temp2.element) % temp1.quotient)) && (! (bool)(temp1.quotient % temp2.quotient)) && !((bool)(temp2.quotient % temp1.quotient));
    };

    template<class U> explicit operator corps_quotient<U>() const {
        corps_quotient<U> result((U)element, (U)quotient);
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

    friend std::ostream& operator<<(std::ostream& os, corps_quotient<T> const& temp) {
        //        element.getDegre();
        os << "(" << temp.element << " % " << temp.quotient << ")";
        return os;
    };


    T element;
    T quotient;
};

