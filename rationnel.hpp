//#ifndef RATIONNEL_H
//#define RATIONNEL_H

#pragma once
#include <vector>
#include <iostream>
#include <exception>
#include <initializer_list>

#include "entete objets.hpp"
#include "norme.hpp"
#include "unite.hpp"
#include "types.hpp"



template<class T> T unite(T const& element,bool test);



template<class T, class enable = void> class rationnel;


template<class T> class rationnel<T,typename std::enable_if_t<type_algebre<T>::type==1>> {
public:

    friend rationnel<T> derivee(rationnel<T> const& element) {
        rationnel<T> result;
        result.numerateur = (derivee(element.numerateur) * element.denominateur) - (element.numerateur * derivee(element.denominateur));
        result.denominateur = element.denominateur * element.denominateur;
        result.simplifier();
        return result;
    };

    template<class U> explicit rationnel(std::initializer_list<U> liste) {
        if ((liste.size() == 0) || (liste.size() >2))
            throw std::domain_error("initialisation de rationnel : liste vide ou >2");
        if (liste.size() == 1) {
            numerateur = T(liste[0]);
            denominateur = unite(numerateur, true);
            simplifier();
        }
        if (liste.size() == 2) {
            numerateur =T(liste[0]);
            denominateur = T(liste[1]);
            simplifier();
        }
    }

    explicit rationnel() : numerateur(),denominateur(){};

    explicit rationnel(T const& a, T const& b) {
        numerateur = a;
        denominateur = b;
        simplifier();
    };

    explicit rationnel(T const& a) {
        numerateur = a;
        denominateur = unite(a,true);
    }


    rationnel(const rationnel<T>& copie) {
        numerateur = copie.numerateur;
        denominateur = copie.denominateur;
    };

    rationnel<T>& operator=(bool test) {
        numerateur = test;
        denominateur = true;
        return *this;

    };

    rationnel<T>& operator=(const rationnel<T>& temp) {
        numerateur = temp.numerateur;
        denominateur = temp.denominateur;

        return (*this);
    };

    /*
    friend rationnel<T> operator%(const rationnel<T>& temp1, const rationnel<T>& temp2) {
        rationnel<T> result = temp1;
        result = false;
        return result;
    };*/

    /*
    T PGCD() const {
        T x(numerateur);
        T y(denominateur);

        if (numerateur > denominateur) {
            x = numerateur;
            y = denominateur;
        }
        else {
            x = denominateur;
            y = numerateur;
        }

        T r(x);
        r = true;

        T faux(r);
        faux = false;

        if (! (bool) y)
            return faux;

        if (! (bool) x)
            return faux;

        while ((bool) r) {
            r = (x%y);
//            r = x - y * (x / y);

            x = y;
            y = r;
        };

        return x;
    };
    */

    void simplifier() {
        T pgcd = PGCD(numerateur, denominateur);

        if ((bool)pgcd) {
            numerateur = numerateur / pgcd;
            denominateur = denominateur / pgcd;
        }
        if (!(bool)numerateur)
            if ((bool)denominateur)
                denominateur = unite(denominateur,true);
        if (!(bool)denominateur)
            if ((bool)numerateur)
                numerateur = unite(numerateur,true);

        return;
    };
    //    template<typename U> 
    friend rationnel<T> operator* (const rationnel<T>& temp1, const rationnel<T>& temp2) {
        rationnel<T> result;
        result.numerateur = (temp1.numerateur * temp2.numerateur);
        result.denominateur = (temp1.denominateur * temp2.denominateur);
        result.simplifier();
        return result;
    };

    friend rationnel<T> operator*(const T& scalaire, const rationnel<T>& temp) {
        rationnel<T> result;
        result.numerateur = scalaire * temp.numerateur;
        result.denominateur = temp.denominateur;
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
        if (!(bool)result) {
//            std:cerr << "erreur : division par zero lors d'une évaluation de polynome" << std::endl;
            throw std::domain_error("division par zero lors d'une évaluation de fraction polynome");
        }
        return (numerateur(element) / result);
    };

    friend rationnel<T> operator/(const rationnel<T>& temp1, const rationnel<T>& temp2) {
        rationnel<T> result;
        result.numerateur = (temp1.numerateur * temp2.denominateur);
        result.denominateur = (temp1.denominateur * temp2.numerateur);

        result.simplifier();

        return result;
    };

    friend rationnel<T> operator+ (const rationnel<T>& temp1, const rationnel<T>& temp2) {
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
        return (!((bool) ((a.numerateur * b.denominateur) - (b.numerateur * a.denominateur))));
    };


    explicit inline operator bool() const {
        return ((bool) numerateur);
    };

//    explicit operator double() const;
    
    /*
    explicit operator double() const {
        double result = ((double)numerateur) / ((double)denominateur);
        return result;
    };*/
    
    template<class U> explicit operator U() const {
        U result = (U)numerateur / (U)denominateur;
        return result;
    };


    template<class U> explicit operator rationnel<U>()const {
        rationnel<U> result( (U) numerateur, (U) denominateur);
        return result;
    };

    

    friend std::ostream& operator<<(std::ostream& os, const rationnel<T>& element) {
        os << "(" << element.numerateur << "/" << element.denominateur << ")";
        return os;
    };

    T numerateur;
    T denominateur;

};


template<class T> class rationnel<T, typename std::enable_if_t<type_algebre<T>::type == 2>> {
public:

    friend rationnel<T> derivee(rationnel<T> const& element) {
        rationnel<T> result;
        result.numerateur = (derivee(element.numerateur) * element.denominateur) - (element.numerateur * derivee(element.denominateur));
        result.denominateur = element.denominateur * element.denominateur;
        result.simplifier();
        return result;
    };


    explicit rationnel() : numerateur(), denominateur() {};

    explicit rationnel(T const& a, T const& b) {
        numerateur = a;
        denominateur = b;
        simplifier();
    };

    explicit rationnel(T const& a) {
        numerateur = a;
        denominateur = unite(a);
    }


    rationnel(const rationnel<T>& copie) {
        numerateur = copie.numerateur;
        denominateur = copie.denominateur;
    };

    template<class U> explicit rationnel(std::initializer_list<U> liste) {
        if ((liste.size() == 0) || (liste.size() > 2))
            throw std::domain_error("initialisation de rationnel : liste vide ou >2");
        if (liste.size() == 1) {
            numerateur = T(liste[0]);
            denominateur = unite(numerateur, true);
            simplifier();
        }
        if (liste.size() == 2) {
            numerateur = T(liste[0]);
            denominateur = T(liste[1]);
            simplifier();
        }
    }

    rationnel<T>& operator=(bool test) {
        numerateur = test;
        denominateur = true;
        return *this;

    };

    rationnel<T>& operator=(const rationnel<T>& temp) {
        numerateur = temp.numerateur;
        denominateur = temp.denominateur;

        return (*this);
    };

    /*
    friend rationnel<T> operator%(const rationnel<T>& temp1, const rationnel<T>& temp2) {
        rationnel<T> result = temp1;
        result = false;
        return result;
    };*/

    /*
    T PGCD() const {
        T x(numerateur);
        T y(denominateur);

        if (numerateur > denominateur) {
            x = numerateur;
            y = denominateur;
        }
        else {
            x = denominateur;
            y = numerateur;
        }

        T r(x);
        r = true;

        T faux(r);
        faux = false;

        if (! (bool) y)
            return faux;

        if (! (bool) x)
            return faux;

        while ((bool) r) {
            r = (x%y);
//            r = x - y * (x / y);

            x = y;
            y = r;
        };

        return x;
    };
    */

    inline void simplifier() {
        /*
        T pgcd = PGCD(numerateur, denominateur);

        if ((bool)pgcd) {
            numerateur = numerateur / pgcd;
            denominateur = denominateur / pgcd;
        }
        if (!(bool)numerateur)
            if ((bool)denominateur)
                denominateur = true;
        if (!(bool)denominateur)
            if ((bool)numerateur)
                numerateur = true;
                */
        return;
    };
    //    template<typename U> 
    friend rationnel<T> operator* (const rationnel<T>& temp1, const rationnel<T>& temp2) {
        rationnel<T> result;
        result.numerateur = (temp1.numerateur * temp2.numerateur);
        result.denominateur = (temp1.denominateur * temp2.denominateur);
        result.simplifier();
        return result;
    };

    friend rationnel<T> operator*(const T& scalaire, const rationnel<T>& temp) {
        rationnel<T> result;
        result.numerateur = scalaire * temp.numerateur;
        result.denominateur = temp.denominateur;
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
        if (!(bool)result) {
            //            std:cerr << "erreur : division par zero lors d'une évaluation de polynome" << std::endl;
            throw std::domain_error("division par zero lors d'une évaluation de fraction de polynome");
        }
        return (numerateur(element) / result);
    };

    friend rationnel<T> operator/(const rationnel<T>& temp1, const rationnel<T>& temp2) {
        rationnel<T> result;
        result.numerateur = (temp1.numerateur * temp2.denominateur);
        result.denominateur = (temp1.denominateur * temp2.numerateur);

        result.simplifier();

        return result;
    };

    friend rationnel<T> operator+ (const rationnel<T>& temp1, const rationnel<T>& temp2) {
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
        return ((bool) numerateur);
    };

    //    explicit operator double() const;

        /*
        explicit operator double() const {
            double result = ((double)numerateur) / ((double)denominateur);
            return result;
        };*/

    template<class U> explicit operator U() const {
        U result = (U)numerateur / (U)denominateur;
        return result;
    };


    template<class U> explicit operator rationnel<U>()const {
        rationnel<U> result;
        result.numerateur = (U)numerateur;
        result.denominateur = (U)denominateur;
        return result;
    };



    friend std::ostream& operator<<(std::ostream& os, const rationnel<T>& element) {
        os << "(" << element.numerateur << "/" << element.denominateur << ")";
        return os;
    };

    T numerateur;
    T denominateur;

};


