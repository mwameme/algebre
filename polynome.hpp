//#ifndef POLYNOME_H
//#define POLYNOME_H

#pragma once 

#include <vector>
#include <iostream>
#include <exception>
#include <initializer_list>

#include "entete objets.hpp"
#include "types.hpp"
#include "unite.hpp"
#include "convertir_liste.hpp"

template<class T> class type_algebre;

template<class T, class enable1 = void, class enable2 = void> class matrice;

template<typename T> class polynome {
public:

    
    friend polynome<T> derivee(polynome<T> const& element) { // T = rationnel<InfInt>
        polynome<T> result(element);

        for (long i(0); i < result.coeffs.size() - 1; ++i) {
            result.coeffs[i] = ((i + 1)) * result.coeffs[i + 1];
        }
        result.coeffs[result.coeffs.size() - 1] = unite(result.coeffs[0],false);
        result.getDegre();
        return result;
    };
    

    explicit polynome(std::vector<T> tab) {
        coeffs = tab;
        if (coeffs.size() == 0) {
            throw std::domain_error("polynome : la taille du vecteur des coeffs doit ?tre > 0");
        }
        getDegre();
    };

    explicit polynome() { // : degre(-1), coeffs(1) {
    }; //probl?me si utilis? ... est juste l? pour les fonctions template.

    explicit polynome(const T& element) {
        coeffs.resize(1);
        coeffs[0] = element;
        getDegre();
    };

    polynome(const polynome<T>& copie) {
        coeffs = copie.coeffs;
        degre=copie.degre;
    };

    explicit polynome(T const& element1, T const& element2) {
        coeffs.resize(2);
        coeffs[0] = element1;
        coeffs[1] = element2;
        getDegre();
    };

    explicit polynome(T const& element1, T const& element2, T const& element3) {
        coeffs.resize(3);
        coeffs[0] = element1;
        coeffs[1] = element2;
        coeffs[2] = element3;
        getDegre();
    };

    /*
    template<class U> explicit polynome(std::initializer_list<U> liste) {
        std::vector<U> vec(liste);
        if (vec.size() == 0)
            throw std::domain_error("initialisation de polynome : liste vide");
        coeffs.resize(vec.size());
        for (int i(0); i < vec.size(); ++i)
            coeffs[i] = T(vec[i]);
        getDegre();
    }*/

    template<class U> explicit polynome(std::initializer_list<U> const & liste) {
        auto vec = convertir_T<std::initializer_list<U>>::convertir(liste);
        coeffs.resize(vec.size());
        for (int i(0); i < vec.size(); ++i)
            coeffs[i] = T(vec[i]);
        getDegre();
    };

    template<class U> explicit polynome(std::vector<U> const& vec) {
        coeffs.resize(vec.size());
        for (int i(0); i < vec.size(); ++i)
            coeffs[i] = T(vec[i]);
        getDegre();
    };



    polynome<T>& operator=(const polynome<T>& temp) {
        if (this == &temp)
            return *this;
        coeffs = temp.coeffs;
        degre = temp.degre;
        return (*this);
    };

    /*
    polynome<T>& operator=(bool test) {
        coeffs.resize(1);
        coeffs[0] = test;
        getDegre();
        return (*this);
    };
    */


    friend polynome<T> operator*(const polynome<T>& temp1, const polynome<T>& temp2) {

        polynome<T> result;
        T faux = unite(temp1.coeffs[0], false);

        if ((temp1.degre < 0) || (temp2.degre < 0))
            return unite(temp1,false);

        result.coeffs = std::vector<T>(temp1.degre + temp2.degre + 1,faux);
        /*
        for (int i(0); i < result.coeffs.size(); ++i) {
            T temp(faux);
            for (int j(0); j <= i; ++j) {
                if ((j < temp1.coeffs.size()) && ((i - j) < temp2.coeffs.size()))
                    temp = temp + (temp1.coeffs[j] * temp2.coeffs[i - j]);
            }
            result.coeffs[i] = temp;
        }*/

        for (int i(0); i <= temp1.degre; ++i)
            for (int j(0) ; j <= temp2.degre; ++j)
                result.coeffs[i + j] = result.coeffs[i + j] + (temp1.coeffs[i] * temp2.coeffs[j]);

        result.getDegre();

        return result;
    };

    friend polynome<T> operator*(const T& scalaire, const polynome<T>& temp) {
        polynome<T> result(temp);

        for (int i(0); i < temp.coeffs.size(); ++i)
            result.coeffs[i] = scalaire * result.coeffs[i];
        result.getDegre();
        return result;
    };

    template<class U> friend polynome<T> operator*(const U& scalaire, const polynome<T>& temp) { 
        //sert pour descendre les niveaux de construction. utilis? pour la d?riv?e : multiplier par un entier.
        polynome<T> result(temp);

        for (int i(0); i < temp.coeffs.size(); ++i)
            result.coeffs[i] = scalaire * result.coeffs[i];
        result.getDegre();
        return result;
    };

    friend polynome<T> operator%(const polynome<T>& temp1, const polynome<T>& temp2) {
        if (type_algebre<T>::type != 0)
            throw std::domain_error("modulo de polynomes : n?cessite une division exacte sur T");

        T faux = unite( temp1.coeffs[0],false);

        if (temp2.degre > temp1.degre)
            return temp1;
        if (temp2.degre < 0)
            return temp1;
        if (temp1.degre < 0)
            return temp1;

        polynome<T> A(temp1);
        while (A.degre >= temp2.degre) {
            T temp = A.coeffs[A.degre] / temp2.coeffs[temp2.degre];
            int degre = A.degre - temp2.degre;
            for (int i(0); i <= temp2.degre; ++i)
                A.coeffs[i + degre] = A.coeffs[i + degre] - temp * temp2.coeffs[i];
            
            A.coeffs[A.degre] = faux; // au cas o? ...
            A.getDegre();
        }

        return A;
    };

    polynome<T>& operator%=(polynome<T> Q) {
        if (Q.degre > degre)
            return *this;
        else
            return (*this = (*this % Q));
    }


    friend polynome<T> operator/( polynome<T> const& temp1, polynome<T> const& temp2) {

        if (type_algebre<T>::type != 0)
            throw std::domain_error("division de polynomes : n?cessite une division exacte sur T");


        if (temp2.degre < 0) 
            throw std::domain_error("division de polynomes : division par 0");


        if (temp1.degre < 0) {
            return unite(temp1, false);
        }

        if (temp2.degre > temp1.degre) {
            return unite(temp1, false);
        }


        polynome<T> A(temp1);
        polynome<T> B(temp2);
        A.getDegre();
        B.getDegre();

        T faux = unite(temp1.coeffs[0],false);

        polynome<T> result;
        result.coeffs = std::vector<T>(A.degre - B.degre + 1,faux);

        while (A.degre >= B.degre) {
            T temp = (A.coeffs[A.degre] / B.coeffs[B.degre]);
            result.coeffs[A.degre - B.degre] = temp;

            int degre = A.degre - B.degre;
            for (int i(0); i < B.degre; ++i) 
                A.coeffs[i + degre] = A.coeffs[i + degre] -  ( temp * B.coeffs[i]);

            A.coeffs[A.degre] = faux;
            A.getDegre();
        }

        result.getDegre();

        return result;
    };

    friend polynome<T> operator+(const polynome<T>& temp1, const polynome<T>& temp2) {
        polynome<T> result;

        if (temp1.degre < 0)
            return temp2;
        if (temp2.degre < 0)
            return temp1;

//        polynome<T> result(temp1.coeffs[0]);

        if (temp1.degre >= temp2.degre) {
            result.coeffs.resize(temp1.degre + 1);
            for (int i(0); i <= temp2.degre; ++i) {
                result.coeffs[i] = (temp1.coeffs[i] + temp2.coeffs[i]);
            }
            for (int i(temp2.degre + 1); i <= temp1.degre; ++i) {
                result.coeffs[i] = temp1.coeffs[i];
            }
        }
        else {
            result.coeffs.resize(temp2.degre + 1);
            for (int i(0); i <= temp1.degre; ++i) {
                result.coeffs[i] = (temp1.coeffs[i] + temp2.coeffs[i]);
            }
            for (int i(temp1.degre + 1); i <= temp2.degre; ++i) {
                result.coeffs[i] = temp2.coeffs[i];
            }

        }

        result.getDegre();
        return result;
    };

    friend polynome<T> operator-(const polynome<T>& temp1, const polynome<T>& temp2) {
        polynome<T> result;

        if (temp1.degre < 0)
            return (-temp2);
        if (temp2.degre < 0)
            return temp1;

//        polynome<T> result(temp1.coeffs[0]);

        if (temp1.degre >= temp2.degre) {
            result.coeffs.resize(temp1.degre + 1);
            for (int i(0); i <= temp2.degre; ++i) {
                result.coeffs[i] = (temp1.coeffs[i] - temp2.coeffs[i]);
            }
            for (int i(temp2.degre + 1); i <= temp1.degre; ++i) {
                result.coeffs[i] = temp1.coeffs[i];
            }
        }
        else {
            result.coeffs.resize(temp2.degre + 1);
            for (int i(0); i <= temp1.degre; ++i) {
                result.coeffs[i] = (temp1.coeffs[i] - temp2.coeffs[i]);
            }
            for (int i(temp1.degre + 1); i <= temp2.degre; ++i) {
                result.coeffs[i] = - temp2.coeffs[i];
            }

        }

        result.getDegre();
        return result;
    };



    friend polynome<T> operator-(const polynome<T>& temp) {
        polynome<T> result(temp);

        for (int i(0); i < result.coeffs.size(); ++i)
            result.coeffs[i] = - result.coeffs[i];

        result.getDegre(); // enlever ?
        return result;
    };

    template<class U> U operator()(const U& element) const{
        U result = unite(element,false);
        U puissance = unite(result,true);

        if (degre == -1)
            return result;

        for (int i(0); i <= degre; ++i) {
            result = result + (coeffs[i] * puissance);
            puissance = puissance * element;
        }

        return result;
    };


    friend bool operator==(const polynome<T>& a, const polynome<T>& b) {
        if (a.degre != b.degre)
            return false;
        for (int i(0); i <= a.degre; ++i)
            if (!(a.coeffs[i] == b.coeffs[i]))
                return false;
        return true;
    };

    friend bool operator>(const polynome<T>& a, const polynome<T>& b) { // ... ... ... 
        return (a.degre > b.degre);
    };



    explicit inline operator bool() const {
        if (degre < 0)
            return false;
        for (int i(coeffs.size() -1); i >= 0; --i) {
            if ((bool) coeffs[i]) {
                return true;
            }
        }
        return false;
    };

    template<class U> explicit operator polynome<U>() const {
        polynome<U> result();
        result.coeffs.resize(coeffs.size());
        for (int i(0); i < coeffs.size(); ++i) {
            result.coeffs[i] = (U)coeffs[i];
        }
        result.getDegre();
        return result;
    };


    void majListe() {
        int i(coeffs.size() - 1);

        while (! (bool) coeffs[i]) {
            --i;
            if (i == 0)
                break;
        };
        if( i < coeffs.size()-1)
            coeffs.erase(coeffs.begin() + i + 1, coeffs.end());
        return;
    };

    inline void getDegre() {
        majListe();

        if (coeffs.size() == 0)
            throw std::domain_error("la liste d'un polynome a 0 ?l?ments");

        int i;
        for (i = (coeffs.size() - 1) ; i >= 0; --i) 
            if ((bool)coeffs[i]) 
                break;
        
        degre = i;
    };

    friend std::ostream& operator<<(std::ostream& os, const polynome<T>& element) {
        //        element.getDegre();
        os << "(";
        if (element.degre == -1)
            os << "0 *X^0)";
        else {
            for (int i(0); i < element.degre; ++i) {
                os << element.coeffs[i] << " *X^" << i << " + ";
            }
            os << element.coeffs[element.degre] << " *X^" << element.degre << ")";
        }
        return os;
    };

    friend T resultant(const polynome<T>& poly1, const polynome<T>& poly2) {
        T faux_ = unite(poly1.coeffs[0],false);

        if ((poly1.degre < 0) || (poly2.degre < 0))
            return faux_;

        matrice<T> m_matrice(poly1.degre + poly2.degre, poly1.degre + poly2.degre, faux_);

        for (int i(0); i < poly2.degre; ++i)
            for (int j(0); j < poly1.coeffs.size(); ++j)
                m_matrice.coeffs[i][i + j] = poly1.coeffs[j];

        for (int i(0); i < poly1.degre; ++i)
            for (int j(0); j < poly2.coeffs.size(); ++j)
                m_matrice.coeffs[i + poly2.degre][i + j] = poly2.coeffs[j];

        return m_matrice.determinant();
    };

    int multiplicite_max() {
        if (degre < 0)
            return -1;
        if (degre == 0)
            return 0;

        /*
        polynome<T> R;
        polynome<T> P(*this);
        while ((R = PGCD(P, derivee(P))).degre >= 1) { // s'ex?cute une seule fois ?
            P = P / R; 
        }//P ? racines simples
        */

        polynome<T> P = *this / PGCD(*this, derivee(*this));

        polynome<T> P_n(P);
        int multiplicite = 1;
        while ((bool)(P_n % *this)) {
            P_n = P_n * P;
            ++multiplicite;
        }
        return multiplicite;
    };

    //protected:
    int degre;
    std::vector<T> coeffs;

};


