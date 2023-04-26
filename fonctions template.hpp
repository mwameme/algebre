#pragma once

#include "unite.hpp"

template<class T>  T unite(T const& element, bool test);


template<class T> T derivee(T const& element) {
    return unite(element,false);
};


template<class T> inline T carre(T temp) {
    return temp * temp;
};

template<class T> inline T abs(T temp) {
    if (temp > 0)
        return temp;
    else
        return -temp;
};

int factorielle(int n) {
	if (n <= 1)
		return 1;
	int fact = 1;
	for (int i(2); i <= n; ++i)
		fact = i * fact;
	return fact;
};


template<class T> T PGCD(T x, T y) {

    if (!(bool)y)
        return y;

    if (!(bool)x)
        return x;

    T r = x;
    do {
        r = (x % y);
        x = y;
        y = r;
    } while ((bool)r);

    return x;
};

template<class T> T PPCM(T const& a, T const& b) {
    return (a * b) / PGCD(a, b);
};

template<class T> T puissance(T element, int n) {
    T resultat = unite( element,true);

    T puissance_m = element;
    while (n > 0) {
        if ((n % 2) == 1)
            resultat = resultat * puissance_m;
        n = n / 2;
        if (n > 0)
            puissance_m = puissance_m * puissance_m;
    }
    return resultat;
};

template<class T> T  inverse(T b,T quotient) { // inverse de b dans T*q 
    T a, //b, //ua, va, wa, 
        ub, vb, wb, c, q;

    a = quotient;
//    b = x;
    //        ua = b; ua = true;
    vb = unite(b,true);
    ub = unite(b,false);
    //        va = ub;
    c = vb;
    while ((bool)c)
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
    return ub/a; //important pour les polynomes ... si le dernier reste est de degré 0, il peut ne pas être =1
};


template<class T> T max(T a, T b) {
    if (a > b) return a;
    else return b;
};


template<class T> T min(T a, T b) {
    if (a > b) return b;
    else return a;
};