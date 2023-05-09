// algebre.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#include <iostream>
#include <vector>
#include "entete objets.hpp"


#include "interpolation.hpp"
#include "diagonalisation.hpp"

#include "InfInt.h"
#include "precision/fprecision.h"
#include "precision/iprecision.h"

#include "fact_for.hpp"
#include "simplifier polynome_n.hpp"

#include "convertir_liste.hpp"

#include "reduire_frac.hpp"

using namespace std;

#ifdef ALGEBRA_USE_EXCEPTION

#endif


int main()
{
    long question;
//    auto x = 4.5 % 4.1;

    if (true) {
        monome<rationnel<int>> monome_(std::vector<int>({ 1,2,3 }), rationnel<int>(3, 2));
        polynome_n_fixe<rationnel<int>, 3> poly(monome_);
        polynome_n_fixe<rationnel<int>, 3> poly2(std::vector<int>({ 1,2,3 }), rationnel<int>(3, 2));
        polynome_n_fixe<rationnel<int>, 3>::const_iterator it = poly.cbegin();
        polynome_n_fixe<rationnel<int>, 3>::const_iterator it2 = poly2.cbegin();

        while ((bool)it) {
            cout << *it << " , " << *it2 << " , " << (*it == *it2) << endl;
            ++it;
            ++it2;
        }

        cout << (poly == poly2) << endl;

        cin >> question;

    }

    if (false) {
        polynome<rationnel<int>> P(std::vector<int>({ 1,2,3,4 }));
        cout << type_algebre< polynome<rationnel<int>>>::type;

        cout << endl << "======FIN======" << endl;
        cin >> question;
    }

    if (false) {
        auto x = convertir_T<std::initializer_list<std::initializer_list<int>>>::convertir({ {1,2},{3,4}, {5,6} });
///        auto x = convertir({ {1,2},{3,4}, {5,6} });
        for (int i(0); i < 3; ++i)
            for (int j(0); j < 2; ++j)
                cout << x[i][j];
        cout << endl;

        cin >> question;

        polynome<rationnel<int>> poly(x);
//        polynome<int> poly({ 1,2,3 });
//        rationnel<int> ratio(std::vector<int>{1, 2});
//        cout << ratio << endl;
        cout << poly << endl;
        cin >> question;
    }

    //test simplifier et double conversion.
    if (false) {
        
        std::string noms[3] = { "X","Y","Z" };
        polynome_n_iter<rationnel<int>> poly1({ 3,1,1 }, rationnel<int> {3, 4}, noms);

        polynome_n_iter<rationnel<int>> poly2({ 1,3,1 }, rationnel<int> {5, 4}, noms);
        polynome_n_iter<rationnel<int>> poly3({ 1,1,3 }, rationnel<int> {2, 3}, noms);

        polynome_n_iter<rationnel<int>> poly = poly1 + poly2 + poly3;

        cout << "poly :\n" << poly << std::endl;

        polynome_n_iter<rationnel<int>> poly4 = poly * poly1;
        cout << "poly4 :\n" << poly4 << std::endl;
        poly4.simplifier_2();
        cout << "poly4 simplifie :\n" << poly4 << std::endl;
        polynome_n_rec<rationnel<int>> poly5 = (polynome_n_rec<rationnel<int>>) poly4;
        cout << poly5 << endl;

        cin >> question;
        polynome_n_iter<rationnel<int>> poly_simp = simplifier_poly(poly4, poly1);
        cout << "\n polynome simplifie fraction : \n" << poly_simp << endl;

        cout << endl << "======FIN======" << endl;
        cin >> question;
        
    }

    //constructeur de liste
    if (false) {
        
        matrice<int> mat (vector<vector<int>>{ {1,2,3},{4,5,6},{7,8,9} });

        cout << "matrice \n" << 2 * mat << endl;

        polynome<int> poly{ 1,2,3 };
        cout << poly << endl;

        rationnel<int> ratio{ 1,2 };
        cout << ratio << endl;

        polynome<rationnel<int>> poly_ratio{ rationnel<int>{1,2},rationnel<int>{3,4},rationnel<int>{5,6} };
        cout << poly_ratio << endl;

        cout << endl << "======FIN======" << endl;

        cin >> question;
        
    };

    //test polynome_n_iter
    if (false) {
        std::string noms[3] = { "X","Y","Z" };
        polynome_n_iter<int> poly1({ 3,1,1 }, 1, noms);

        polynome_n_iter<int> poly2({ 1,3,1 }, 2, noms);
        polynome_n_iter<int> poly3({ 1,1,3 }, 3, noms);

        polynome_n_iter<int> poly = poly1 + poly2 + poly3;

        cout << "poly :\n" << poly << std::endl;
        polynome_n_iter<int> poly4 = poly * poly1;

        cout << "poly4 :\n" << poly4 << std::endl;

        cout << endl << "======FIN======" << endl;
        cin >> question;
    }


    //test antoine
    if (false) {
        polynome<rationnel<int>> poly1(vector<rationnel<int>>{rationnel<int>(-1, 1), rationnel<int>(1, 1)});
        polynome<rationnel<int>> poly2(vector<rationnel<int>>{rationnel<int>(-1, 1), rationnel<int>(0, 1), rationnel<int>(1, 1)});
        cout << poly1 << endl;
        cout << poly2 << endl;

        rationnel<polynome<rationnel<int>>> ratio(poly2, poly1);
        ratio.simplifier();
        cout << ratio << endl;

        cout << "======FIN======" << endl;
        cin >> question;
    }


    //test determinant anneau
    if (false) {
        vector<vector<rationnel<InfInt>>> vec = { {rationnel<InfInt>(3,7) , rationnel<InfInt>(5,2)} ,
        { rationnel<InfInt>(0,1) , rationnel<InfInt>(8,9) } };

        matrice<rationnel<InfInt>> m_matrice(vec);

        cout << m_matrice.determinant() << endl;
        cout << m_matrice.determinant_anneau() << endl;

        cout << "======FIN======" << endl;
        cin >> question;
    }


    if (false) {
        int n;
        cout << "n=";
        cin >> n;

        for (fact_for iter(n); (bool)iter; ++iter) {
            cout << iter << endl;
        }

        cout << "======FIN======" << endl;
        cin >> question;
    }

    //tests polynome_n_rec<float>
    if (false) {
        float x = 4.3;
        string noms[3] = { "X","Y","Z" };
        vector<int> dim = { 1,2,3 };
        polynome_n_rec<float> poly1(dim,noms,x);
        cout << poly1 << endl;

        vector<int> dim2 = { 3,4,5 };
        polynome_n_rec<float> poly2(dim2, noms, 3.54);
        cout << poly2 << endl;

        polynome_n_rec<float> poly3 = poly1 + poly2;
        cout << poly3 << endl;

        cout << "n_var" << endl;
        cout << poly1.n_var << endl;
        cout << poly2.n_var << endl;
        cout << poly3.n_var << endl;

        poly1 = false;
        cout << poly1.n_var << endl;
        cout << poly1.poly.coeffs[0].n_var << endl;
        cin >> question;

        polynome_n_rec<float> poly4 = poly3 * poly2;
        cout << poly4 << endl;
        cout << poly4 + poly2 << endl;
//        auto f = norme(poly4);

        cout << norme(poly4 + poly2) << endl;

        cout << norme(polynome<float>(1.2, 1.3)) << endl;

        /*
        cout << poly1.nul << endl;
        cout << poly1.coeffs.size() << endl;
        cout << poly1.coeffs[1]->coeffs.size() << endl;

        cout << poly1.noms_variables[0] << endl;
        cout << poly1.coeffs[0]->noms_variables[0] << endl;
        cout << poly1.coeffs[1]->coeffs[2]->coeffs[3]->element << endl;
        
        cout << poly1.coeffs[1]->coeffs[2]->coeffs[2]->element << endl;
        cout << "tests nuls" << endl;
        cout << poly1.nul << endl;
        cout << poly1.coeffs[1]->nul << endl;
        cout << poly1.coeffs[1]->coeffs[2]->nul << endl;
        cout << poly1.coeffs[1]->coeffs[2]->coeffs[3]->nul << endl;
        */

        cout << "======FIN======" << endl;
        cin >> question;
    }
    if (false) {
        cout << derivee(5) << endl << derivee(1.5) << endl;

        cout << "======FIN======" << endl;
        cin >> question;
    }

    //test norme
    if (false) {
        /*
        vector<rationnel<InfInt>> mon_vecteur1 = { rationnel<InfInt>(8,1),rationnel<InfInt>(2,3), rationnel<InfInt>(1,1) }; //X^2 +2X +8
        vector<rationnel<InfInt>> mon_vecteur2 = { rationnel<InfInt>(7,1),rationnel<InfInt>(1,1) }; //X+3
        polynome<rationnel<InfInt>> mon_poly1(mon_vecteur1);
        polynome<rationnel<InfInt>> mon_poly2(mon_vecteur2);

        rationnel<InfInt> ratio1(8, 1);
        rationnel<InfInt> ratio2(1, 3);

        cout << norme_T<InfInt>::norme(InfInt(-2)) << endl;
        auto norme_T<rationnel<InfInt>>::norme(ratio1);
        */

        vector<rationnel<int_precision>> mon_vecteur1 = { rationnel<int_precision>(8,1),rationnel<int_precision>(2,3), rationnel<int_precision>(1,1) }; //X^2 +2X +8
        vector<rationnel<int_precision>> mon_vecteur2 = { rationnel<int_precision>(7,1),rationnel<int_precision>(1,1) }; //X+3
        polynome<rationnel<int_precision>> mon_poly1(mon_vecteur1);
        polynome<rationnel<int_precision>> mon_poly2(mon_vecteur2);

        rationnel<int_precision> ratio1(8, 1);
        rationnel<int_precision> ratio2(1, 3);

        cout << norme_T<int_precision>::norme(int_precision(-2)) << endl;

        polynome<int_precision> poly(std::vector<int_precision>{1, 3, 5, 7});
        cout << poly << endl;
        auto poly_norme =norme_T<polynome<int_precision>>::norme(poly);
        cout << "norme : " << poly_norme << endl;
        cout << typeid(decltype(norme_T<int_precision>::norme(int_precision()))).name() << endl;

        auto norme_ratio = norme_T<rationnel<int_precision>>::norme(ratio1);
        cout << norme_ratio << endl;

//        auto norme_T<rationnel<int_precision>>::norme(ratio1);
//        cout << typeid(ratio1).name() << endl;
//        cout << norme_T<InfInt>::norme(InfInt(-2)) << endl;
//        cout << norme_T<rationnel<InfInt>>::norme(ratio1) << endl;
//        cout << norme_T<rationnel<InfInt>>::norme(ratio2) << endl;

        auto m_norme = norme_T< polynome<rationnel<int_precision>>>::norme(mon_poly1);
        cout << typeid(m_norme).name() << endl;
        cout << m_norme << endl;

        rationnel<polynome<rationnel<int_precision>>> frac_poly(mon_poly1, mon_poly2);
        cout << norme_T<rationnel<polynome<rationnel<int_precision>>>>::norme(frac_poly) << endl;

        cout << "fonction norme" << endl;
        cout << norme(frac_poly) << endl;

        anneau_quotient<polynome<rationnel<int_precision>>> quotient(mon_poly1, mon_poly2);
        cout << "norme anneau_quotient" << endl;
        cout << norme(quotient) << endl;
        cout << "PGCD" << endl;
        cout << PGCD(quotient.element, quotient.quotient) << endl;
//        cout << norme_T< polynome<rationnel<InfInt>>>::norme(mon_poly2) << endl;

        cout << "======FIN======" << endl;

        cin >> question;

    }


    //test specialisation class template (futures matrices)
    if (false) {
        /*
        cout << "int" << endl;
        matrice_test<int> test1;

        cout << endl << "float" << endl;
        matrice_test<float> test2;

        cout << endl << "rationnel<int>" << endl;
        matrice_test<rationnel<int>> test3;
        cout << endl;

        cout << endl << "polynome<float>" << endl;
        matrice_test<polynome<float>> test4;
        cout << endl;
        */


        cout << "entier : type : " << type_algebre<int>::type << endl;
        cout << "entier : approx : " << type_algebre<int>::approx << endl << endl;

        cout << "rationnel<entier> : type : " << type_algebre<rationnel<int>>::type << endl;
        cout << "rationnel<entier> : approx : " << type_algebre<rationnel<int>>::approx << endl << endl;

        cout << "reel : type : " << type_algebre<float>::type << endl;
        cout << "reel : approx : " << type_algebre<float>::approx << endl;

        cout << "polynome<reel> : type : " << type_algebre<polynome<float>>::type << endl;
        cout << "polynome<reel> : approx : " << type_algebre<polynome<float>>::approx << endl;

        cout << "======FIN======" << endl;


        cin >> question;
    }

    if (false) {
        cout << type_algebre<int>::type << endl;
        cout << type_algebre<float>::type << endl;
        cout << type_algebre<rationnel<int>>::type << endl;
        cout << type_algebre<rationnel<polynome<rationnel<int>>>>::type << endl;
        cin >> question;

    }

    //PGCD et résultant
    if (false) {
        vector<rationnel<InfInt>> mon_vecteur1 = { rationnel<InfInt>(8,1),rationnel<InfInt>(2,3), rationnel<InfInt>(1,1) }; //X^2 +2X +8
        vector<rationnel<InfInt>> mon_vecteur2 = { rationnel<InfInt>(7,1),rationnel<InfInt>(1,1) }; //X+3
        polynome<rationnel<InfInt>> mon_poly1(mon_vecteur1);
        polynome<rationnel<InfInt>> mon_poly2(mon_vecteur2);
        auto mon_poly12 = 2*mon_poly1;
        auto mon_poly22 = 3*mon_poly2;
        
        cin >> question;
        polynome<rationnel<InfInt>> m_pgcd = PGCD(mon_poly1,mon_poly2);
        cout << "poly 1 :" << mon_poly1 << endl;
        cout << "poly 2 :" << mon_poly2 << endl;
        cout << "PGCD :" << m_pgcd << endl;
        cout << "resultant : " << resultant(mon_poly1, mon_poly2) << endl << endl;

        cout << "poly1 * 2, poly2 * 3" << endl;
        polynome<rationnel<InfInt>> m_pgcd2 = PGCD(mon_poly12, mon_poly22);
        cout << "poly 1 :" << mon_poly12 << endl;
        cout << "poly 2 :" << mon_poly22 << endl;
        cout << "PGCD :" << m_pgcd2 << endl;
        cout << "resultant : " << resultant(mon_poly12, mon_poly22) << endl;

        cout << "======FIN======" << endl;

        cin >> question;

    }


    //test précision_relative ( float_precision)
    if (false) {
        float_precision test1(100000.);
        float_precision test2(1.);
        cout << (float)test1.epsilon() << endl;
        cout << (float)test2.epsilon() << endl;
        int power = test1.precision();
        float res = puissance((float).1, power);
        cout << res << endl;

        cout << "======FIN======" << endl;

        cin >> question;
    };

    //test float_precision : diagonalisation
    if (false) {
        float_precision x11(1.23), x12(1.54), x21(0.27), x22(0.31);

        vector<vector<erreur_b<float_precision>>> vec2 = { {erreur_b<float_precision>(x11) , erreur_b<float_precision>(x12)},{erreur_b<float_precision>(x21),erreur_b<float_precision>(x22)} };

//        vector<vector<erreur_b<float_precision>>> vec2 = { {erreur_b<float_precision>(1.23) , erreur_b<float_precision>(1.54)},{erreur_b<float_precision>(0.27),erreur_b<float_precision>(0.31)} };
        matrice<erreur_b<float_precision>> m_matrice2(vec2);
        calcul_condition<erreur_b<float_precision>> calcul2(m_matrice2);
        cin >> question;

        cout << "matrice : " << m_matrice2 << endl;
        auto det2 = m_matrice2.determinant();
        cout << "determinant " << endl << det2 << endl << endl;

        auto inv2 = m_matrice2.inverse();
        cout << "inverse :" << endl << inv2 << endl;
        cout << "inverse * matrice : " << endl << inv2 * m_matrice2 << endl;

        polynome<erreur_b<float_precision>> polyCar2 = m_matrice2.polynomeCaracteristique();
        cout << "polynome Caracteristique" << endl << polyCar2 << endl << endl;

        cout << "multiplicte max " << polyCar2.multiplicite_max() << endl;

        cout << "polynome_car(matrice)" << endl << polyCar2(m_matrice2) << endl << endl;

        cin >> question;

        sous_ev_condition<erreur_b<float_precision>>* m_sous_ev2 = calcul2.calculer();
        cout << "afficher le résultat du calcul des conditions" << endl;
        cout << *m_sous_ev2 << endl;
        cout << "======FIN======" << endl;
        cin >> question;

    };

    //test erreur_b<double> : diagonalisation
    if (false) {
        vector<vector<erreur_b<double>>> vec2 = { {erreur_b<double>(1.23) , erreur_b<double>(1.54)},{erreur_b<double>(0.27),erreur_b<double>(0.31)} };
        matrice<erreur_b<double>> m_matrice2(vec2);
        calcul_condition<erreur_b<double>> calcul2(m_matrice2);
        cin >> question;

        cout << "matrice : " << m_matrice2 << endl;
        auto det2 = m_matrice2.determinant();
        cout << "determinant " << endl << det2 << endl << endl;

        auto inv2 = m_matrice2.inverse();
        cout << "inverse :" << endl << inv2 << endl;
        cout << "inverse * matrice : " << endl << inv2 * m_matrice2 << endl;

        polynome<erreur_b<double>> polyCar2 = m_matrice2.polynomeCaracteristique();
        cout << "polynome Caracteristique" << endl << polyCar2 << endl << endl;

        cout << "multiplicte max " << polyCar2.multiplicite_max() << endl;

        cout << "polynome_car(matrice)" << endl << polyCar2(m_matrice2) << endl << endl;

        cin >> question;

        sous_ev_condition<erreur_b<double>>* m_sous_ev2 = calcul2.calculer();
        cout << "afficher le résultat du calcul des conditions" << endl;
        cout << *m_sous_ev2 << endl;
        cout << "======FIN======" << endl;

        cin >> question;
    }
    
    //test tableau de string
    if (false) {
        string noms[3] = { "variable1","variable2","variable3" };
        cout << *noms << endl;
        cout << *(noms + 1) << endl;
        cout << *(noms + 2) << endl;
        cout << "======FIN======" << endl;

        cin >> question;
    }

    //test polynome et fraction rationnelle
    if (false) {
        //    polynome<rationnel<InfInt>> mon_polynome(mon_vecteur);
        vector<rationnel<InfInt>> mon_vecteur1 = { rationnel<InfInt>(6,1),rationnel<InfInt>(5,1), rationnel<InfInt>(1,1) }; //X^2 +5X +6
        vector<rationnel<InfInt>> mon_vecteur2 = { rationnel<InfInt>(3,1),rationnel<InfInt>(1,1) }; //X+3
        polynome<rationnel<InfInt>> mon_poly1(mon_vecteur1);
        polynome<rationnel<InfInt>> mon_poly2(mon_vecteur2);
        rationnel<polynome<rationnel<InfInt>>> ma_frac(mon_poly1, mon_poly2);

        cin >> question;
        polynome<rationnel<InfInt>> m_pgcd = PGCD(ma_frac.numerateur,ma_frac.denominateur);
        //polynome<rationnel<InfInt>> m_pgcd = (mon_poly1 % mon_poly2);

        cout << ma_frac << endl;
        cout << m_pgcd << endl;
        cin >> question;
        ma_frac.simplifier();
        cin >> question;
        cout << ma_frac << endl;
        cout << "fin simplifier" << endl;
        cout << ma_frac * ma_frac << endl;
        cout << "======FIN======" << endl;

        cin >> question;
    }


    //test corps_quotient et division
    if (false) {
        //tester corps_quotient<polynome> 
        polynome<rationnel<int>> poly(vector<rationnel<int>> {rationnel<int>(1, 1), rationnel<int>(2, 1), rationnel<int>(3, 1), rationnel<int>(1, 1)});
        polynome<rationnel<int>> poly2(vector<rationnel<int>> {rationnel<int>(1, 2), rationnel<int>(4, 1)});
        polynome<rationnel<int>> poly3(vector<rationnel<int>> {rationnel<int>(1, 1), rationnel<int>(1, 1), rationnel<int>(1, 1)});
        corps_quotient<polynome<rationnel<int>>> A(poly2, poly);
        corps_quotient<polynome<rationnel<int>>> B(poly3, poly);
        corps_quotient<polynome<rationnel<int>>> C = A / B;

        cout << C << endl;
        cout << C * B << endl;
        cout << "======FIN======" << endl;

        cin >> question;
    }


    //test diagonalisation : rationnel<InfInt>
    if (false) {
                                                
        vector<vector<rationnel<InfInt>>> vec = { {rationnel<InfInt>(1,1) , rationnel<InfInt>(1,1)} ,
                { rationnel<InfInt>(0,1) , rationnel<InfInt>(1,1) } };

        matrice<rationnel<InfInt>> m_matrice(vec);

        //    polynome<rationnel<InfInt>> m_poly = m_matrice.polynomeCaracteristique();
        cout << "matrice : " << m_matrice << endl;
        auto det = m_matrice.determinant();
        cout << "determinant " << endl << det << endl << endl;

        auto inv = m_matrice.inverse();
        cout << "inverse :" << endl << inv << endl;
        cout << "inverse * matrice : " << endl << inv * m_matrice << endl;

        polynome<rationnel<InfInt>> polyCar = m_matrice.polynomeCaracteristique();
        cout << "polynome Caracteristique" << endl << polyCar << endl << endl;

        cout << "multiplicte max " << polyCar.multiplicite_max() << endl;

        cout << "polynome_car(matrice)" << endl << polyCar(m_matrice) << endl << endl;

        calcul_condition<rationnel<InfInt>> calcul(m_matrice, true);
        cin >> question;

        sous_ev_condition<rationnel<InfInt>>* m_sous_ev = calcul.calculer();
        cout << "afficher le résultat du calcul des conditions" << endl;
        cout << *m_sous_ev << endl;

        cout << "======FIN======" << endl;

        cin >> question;

    }

    if (false)
    {
        vector<vector<rationnel<int_precision>>> vec = { {rationnel<int_precision>(5,4) , rationnel<int_precision>(2,11) , rationnel<int_precision>(0,1),rationnel<int_precision>(0,1)} ,
            { rationnel<int_precision>(4,1) , rationnel<int_precision>(5,1) , rationnel<int_precision>(0,1),rationnel<int_precision>(0,1)} ,
            {rationnel<int_precision>(0,1) , rationnel<int_precision>(0,1) , rationnel<int_precision>(5,4), rationnel<int_precision>(2,11)},
                {rationnel<int_precision>(0,1) , rationnel<int_precision>(0,1) , rationnel<int_precision>(4,1), rationnel<int_precision>(5,1)} };


        matrice<rationnel<int_precision>> m_matrice(vec);

        //    polynome<rationnel<InfInt>> m_poly = m_matrice.polynomeCaracteristique();
        cout << "matrice : " << m_matrice << endl;
        auto det = m_matrice.determinant();
        cout << "determinant " << endl << det << endl << endl;

        auto inv = m_matrice.inverse();
        cout << "inverse :" << endl << inv << endl;
        cout << "inverse * matrice : " << endl << inv * m_matrice << endl;
        cin >> question;

        polynome<rationnel<int_precision>> polyCar = m_matrice.polynomeCaracteristique();
        cout << "polynome Caracteristique" << endl << polyCar << endl << endl;
        cin >> question;

        cout << "multiplicte max " << polyCar.multiplicite_max() << endl;

        cout << "polynome_car(matrice)" << endl << polyCar(m_matrice) << endl << endl;
        cin >> question;

        calcul_condition<rationnel<int_precision>> calcul(m_matrice, true);
        cin >> question;

        sous_ev_condition<rationnel<int_precision>>* m_sous_ev = calcul.calculer();
        cout << "afficher le résultat du calcul des conditions" << endl;
        cout << *m_sous_ev << endl;

        cout << "======FIN======" << endl;

        cin >> question;
    }


    if (false)
    {
        vector<vector<rationnel<InfInt>>> vec = { {rationnel<InfInt>(5,4) , rationnel<InfInt>(2,11) , rationnel<InfInt>(0,1),rationnel<InfInt>(0,1)} ,
            { rationnel<InfInt>(4,1) , rationnel<InfInt>(5,1) , rationnel<InfInt>(0,1),rationnel<InfInt>(0,1)} ,
            {rationnel<InfInt>(0,1) , rationnel<InfInt>(0,1) , rationnel<InfInt>(5,4), rationnel<InfInt>(2,11)},
                {rationnel<InfInt>(0,1) , rationnel<InfInt>(0,1) , rationnel<InfInt>(4,1), rationnel<InfInt>(5,1)} };


        matrice<rationnel<InfInt>> m_matrice(vec);

        //    polynome<rationnel<InfInt>> m_poly = m_matrice.polynomeCaracteristique();
        cout << "matrice : " << m_matrice << endl;
        auto det = m_matrice.determinant();
        cout << "determinant " << endl << det << endl << endl;

        auto inv = m_matrice.inverse();
        cout << "inverse :" << endl << inv << endl;
        cout << "inverse * matrice : " << endl << inv * m_matrice << endl;
        cin >> question;

        polynome<rationnel<InfInt>> polyCar = m_matrice.polynomeCaracteristique();
        cin >> question;
        cout << "polynome Caracteristique" << endl << polyCar << endl << endl;
        cin >> question;

        cout << "multiplicte max " << polyCar.multiplicite_max() << endl;

        cout << "polynome_car(matrice)" << endl << polyCar(m_matrice) << endl << endl;
        cin >> question;

        calcul_condition<rationnel<InfInt>> calcul(m_matrice, true);
        cin >> question;

        sous_ev_condition<rationnel<InfInt>>* m_sous_ev = calcul.calculer();
        cout << "afficher le résultat du calcul des conditions" << endl;
        cout << *m_sous_ev << endl;

        cout << "======FIN======" << endl;

        cin >> question;
    }


    
    if (false) {
        vector<vector<corps_quotient<int>>> vec = { { corps_quotient<int>(1,11), corps_quotient<int>(2,11), corps_quotient<int>(3,11)},
            {corps_quotient<int>(5,11),corps_quotient<int>(7,11),corps_quotient<int>(8,11)}, 
            {corps_quotient<int>(4,11),corps_quotient<int>(2,11),corps_quotient<int>(1,11) } };
        matrice<corps_quotient<int>> m_matrice(vec);

        auto det = m_matrice.determinant();
        auto inv = m_matrice.inverse();
        auto polyCar = m_matrice.polynomeCaracteristique();

        cout << det << endl << endl;
        cout << inv << endl;
        cout << polyCar << endl;
        cout << m_matrice * inv << endl;

        cout << "======FIN======" << endl;

        cin >> question;
    }
    
    
}

// Exécuter le programme : Ctrl+F5 ou menu Déboguer > Exécuter sans débogage
// Déboguer le programme : F5 ou menu Déboguer > Démarrer le débogage

// Astuces pour bien démarrer : 
//   1. Utilisez la fenêtre Explorateur de solutions pour ajouter des fichiers et les gérer.
//   2. Utilisez la fenêtre Team Explorer pour vous connecter au contrôle de code source.
//   3. Utilisez la fenêtre Sortie pour voir la sortie de la génération et d'autres messages.
//   4. Utilisez la fenêtre Liste d'erreur_bs pour voir les erreur_bs.
//   5. Accédez à Projet > Ajouter un nouvel élément pour créer des fichiers de code, ou à Projet > Ajouter un élément existant pour ajouter des fichiers de code existants au projet.
//   6. Pour rouvrir ce projet plus tard, accédez à Fichier > Ouvrir > Projet et sélectionnez le fichier .sln.

