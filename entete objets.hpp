#pragma once
#include "anneau quotient.hpp"
#include "complexe.hpp"
#include "corps quotient.hpp"
#include "matrice.hpp"
#include "polynome.hpp"
#include "polynome_n_rec.hpp"
#include "polynome_n_iter.hpp"
#include "rationnel.hpp"
#include "norme.hpp"
#include "types.hpp"
#include "erreur.hpp"
#include "simplifier polynome_n.hpp"

#include "fonctions template.hpp"
#include <complex>


template<class T> class complex;

template<class T> class erreur;
template<class T> class anneau_quotient;
template<class T> class complexe;
template<class T> class corps_quotient;
template<class T, class enable1 = void, class enable2 = void> class matrice;
template<class T> class polynome;
template<class T> class polynome_n_rec;
template<class T> class polynome_n_iter;
template<class T> class rationnel;
template<class T, int n> class polynome_n_fixe;


//template<class T> class type_algebre;

