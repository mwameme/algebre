#include <algorithm>
#include <vector>
#include <functional>

#include "vote.hpp"
#include "TVD.hpp"
#include "ordre.hpp"

using namespace std;

vector<double> vote(vector<vector<int>> const& votes, double epsilon_){
	int taille = votes.size();//taille : le nombre de votants
	int n = votes[0].size(); //n : le nombre d'�lus / de vot�s
	vector<vector<int>> preferences(n, vector<int>(n, 0));

	for (int i(0); i < taille; ++i) //on parcourt la liste de votes
		for (int j(0); j < n - 1; ++j)
			for (int k(j + 1); k < n; ++k) //on regarde j pr�f�r� � k
				preferences[votes[i][j]][votes[i][k]] += 1;

	vector<vector<double>> P(n, vector<double>(n, 0)); // matrice de transition
	double inv_n = 1. / ((double)n);
	double epsilon = epsilon_ * inv_n;
	double n_2 = (double)n / 2.;

	for (int i(0); i < n; ++i) // on parcourt la matrice de transitions
		for (int j(0); j < n; ++j)
				P[i][j] = inv_n * (1. + epsilon * ((double)preferences[j][i] - n_2));
	
	for(int i(0);i<n;++i)
		P[i][i] = inv_n;

	double somme = 0; //calculer la stochasticit�
	for (int i(0); i < n; ++i) {
		double somme_temp = 0;
		for (int j(0); j < n; ++j)
			somme_temp += P[i][j];
		if (somme_temp > somme)
			somme = somme_temp;
	}

	double somme_inv = 1 / somme;
	for (int i(0); i < n; ++i)
		for (int j(0); j < n; ++j)
			P[i][j] = P[i][j] * somme_inv;

	//modifier les P_ii, toujours pour la stochasticit�
	for (int i(0); i < n; ++i) {
		somme = 0.;
		for (int j(0); j < n; ++j) {
			if (j == i)
				continue;
			somme += P[i][j];
		}
		P[i][i] = 1. - somme;
	}


	vector<double> X(n, inv_n); //calculer la proba stationnaire ...
	X = multiplication(X, puissance(P, 1024));

	somme = 0; //au cas o� X n'est pas strictement de masse 1.
	for (int i(0); i < n; +i)
		somme += X[i];
	somme_inv = 1. / somme;
	for (int i(0); i < n; ++i)
		X[i] *= somme_inv;

	/*
	//renvoyer le r�sultat ...
	vector<double> resultat(n, 0);
	double epsilon_inv = 1. / epsilon;
	for (int i(0); i < n; ++i)
		resultat[i] = (X[i] - inv_n) *epsilon_inv;
		*/
	return X;
}

vector<double> vote_iter(vector<vector<int>> const& votes, double epsilon,vector<double> const& poids) {
	//on oublie les votants, on ne garde que les vot�s :
	//on a donc n==taille
	int n = votes.size(); //nombre de vot�s 
							//NB : nombre de vot�s non-fix� (variable) ... tout parcourir !
	vector<vector<double>> preferences(n, vector<double>( n, 0.));

	for (int i(0); i < n; ++i) //on parcourt la liste de votes
		for (int j(0); j < n - 1; ++j)
			for (int k(j + 1); k < n; ++k) //on regarde j pr�f�r� � k
				preferences[votes[i][j]][votes[i][k]] += poids[i]; //j pr�f�r� � k : de poids i.

	vector<vector<double>> P(n, vector<double>(n, 0)); // matrice de transition
	double inv_n = 1. / ((double)n);
//	double epsilon = epsilon_ * inv_n;

	for (int i(0); i < n; ++i) // on parcourt la matrice de transitions
		for (int j(0); j < n; ++j)
			P[i][j] = inv_n * (1 + epsilon * (preferences[j][i] - .5));
	
	for(int i(0);i<n;++i)
		P[i][i] = inv_n;

	double somme = 0; //calculer la stochasticit�
	for (int i(0); i < n; ++i) {
		double somme_temp = 0;
		for (int j(0); j < n; ++j)
			somme_temp += P[i][j];
		if (somme_temp > somme)
			somme = somme_temp;
	}

	double somme_inv = 1 / somme;
	for (int i(0); i < n; ++i)
		for (int j(0); j < n; ++j)
			P[i][j] = P[i][j] * somme_inv;

	//modifier les P_ii, toujours pour la stochasticit�
	for (int i(0); i < n; ++i) {
		somme = 0.;
		for (int j(0); j < n; ++j) {
			if (j == i)
				continue;
			somme += P[i][j];
		}
		P[i][i] = 1. - somme;
	}


	vector<double> X(n, inv_n); //calculer la proba stationnaire ...
	X = multiplication(X, puissance(P, 1024));

	somme = 0; //au cas o� X n'est pas strictement de masse 1. r�duit l'al�a (lors des it�rations)
	for (int i(0); i < n; +i)
		somme += X[i];
	somme_inv = 1. / somme;
	for (int i(0); i < n; ++i)
		X[i] *= somme_inv;

	//renvoyer le r�sultat ...

	return X;
}

vector<double> vote_n_iter(vector<vector<int>> const& votes, double epsilon, int iter) {
	//rajouter un tableau de vote interne aux vot�s ...
	//comme �a deux �tapes : premier vote, puis vote_iter n fois ...

	int n = votes.size();
	vector<double> X(n, 1 / ((double)n));
	for (int i(0); i < iter; ++i)
		X = vote_iter(votes, epsilon, X);

	return X;
};

vector<double> vote_experience(vector<vector<int>> votes_total, vector<vector<int>> votes_elus, double epsilon, int iter) {
	vector<double> X = vote(votes_total, epsilon);

	for (int i(0); i < iter; ++i) 
		X = vote_iter(votes_elus, epsilon, X);
	
	return X;
}

/*
vector<int> ordre(vector<double> X) {
	int n = X.size();
	vector<couple> liste(n);
	for (int i(0); i < n; ++i) {
		liste[i].nombre = X[i];
		liste[i].index = i;
	}

	sort(liste.begin(), liste.end(), [](couple& x, couple& y) {return x.nombre > y.nombre; });

	vector<int> resultat(n, 0);
	for (int i(0); i < n; ++i)
		resultat[i] = liste[i].index;

	return resultat;
}*/


//ajouter un vote en deux �tapes : 
//d'abord les votant, ensuite les vot�s ...
// faire les tests !