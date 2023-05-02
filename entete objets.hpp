#pragma once

#include "objets.hpp"

template<class T> class complex;

template<class T> class erreur_b;
template<class T> class erreur_l;
template<class T> class anneau_quotient;
template<class T> class complexe;
template<class T> class corps_quotient;
template<class T, class enable1 = void, class enable2 = void> class matrice;
template<class T> class polynome;
template<class T> class polynome_n_rec;
template<class T> class polynome_n_iter;
template<class T> class rationnel;
template<class T, int n> class polynome_n_fixe;

class InfInt;
class int_precision;
class float_precision;

template<class T> class type_algebre;
template<class T, typename Enable = void> class norme_T;
template<class T>  T unite(T const& element, bool test);

template<class T> class scalaire_vecteur;

template<class T> class polynome_n_sparse;
template<class T> class monome;


#include "fonctions template.hpp"
#include "unite.hpp"

#include "types.hpp"
#include "norme.hpp"

#include "swap_T.hpp"

#include "objets.hpp"

//template<class T> class type_algebre;

