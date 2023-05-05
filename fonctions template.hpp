#pragma once

#include "unite.hpp"
#include "types.hpp"
#include <tuple>

template<class T> class type_algebre;
//#include <pair>

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

template<class T>
T factorielle(T n) {
	if (n <= 1)
		return 1;
	int fact = 1;
	for (int i(2); i <= n; ++i)
		fact = i * fact;
	return fact;
};


template<class T> T PGCD(T x, T y) {
    static_assert(type_algebre<T>::type == 1);

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
    static_assert(type_algebre<T>::type == 1);

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
    static_assert(type_algebre<T>::type == 1);

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

//vérifier
template<class T> std::tuple<T, T,T> Bezout(T const& a, T const& b) {
    static_assert(type_algebre<T>::type == 1);

    T u, v, ua, va, wa,
        ub, vb, wb, r, q;

    ua = unite(a, true);;
    ub = unite(a, false);

    vb = unite(a, true);
    va = unite(a, false);
    r = vb;
    //u = ua*a + ub*b; 
    //v = va*a + vb*b;
    //étgalités conservées lors de l'algorithme
    while ((bool)r)
    {
        r = u % v; // r = u - q v
        q = u / v;
        wa = ua - q * va;
        wb = ub - q * vb;
        u = v;
        v = r;
        ua = va;
        ub = vb;
        va = wa;
        vb = wb;
    }
    //u = PGCD(a,b) = ua*a + ub*b.
    return std::make_pair(u, ua , ub ); // pour les polynomes : u peut être de degré 0, mais différent de 1 ...

};