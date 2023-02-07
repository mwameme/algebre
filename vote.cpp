#include <algorithm>
#include <vector>
#include <functional>


#include "vote.hpp"

using namespace std;

int compter_ordonne(vector<vector<int>> const& votes, vector<int> ordre) {
	//renvoit le nombre de fois que ce n'est pas dans l'ordre ... vis-à-vis des préférences.

	int taille = votes.size();//taille : le nombre de votants
	int n = votes[0].size(); //n : le nombre d'élus / de votés
	vector<vector<int>> preferences(n, vector<int>(n, 0));

	for (int i(0); i < taille; ++i) //on parcourt la liste de votes
		for (int j(0); j < n - 1; ++j)
			for (int k(j + 1); k < n; ++k) //on regarde j préféré à k
				preferences[votes[i][j]][votes[i][k]] += 1;

	int ordonne = 0;

	for (int i(0); i < n - 1; ++i)
		for (int j(i + 1); j < n; ++j)
			if (preferences[ordre[i]][ordre[j]] < preferences[ordre[j]][ordre[i]])
				++ordonne;

	return ordonne;
};

vector<double> vote(vector<vector<int>> const& votes, double epsilon_) {
	int taille = votes.size();//taille : le nombre de votants
	int n = votes[0].size(); //n : le nombre d'élus / de votés
	vector<vector<int>> preferences(n, vector<int>(n, 0));

	for (int i(0); i < taille; ++i) //on parcourt la liste de votes
		for (int j(0); j < n - 1; ++j)
			for (int k(j + 1); k < n; ++k) //on regarde j préféré à k
				preferences[votes[i][j]][votes[i][k]] += 1;

	vector<vector<double>> P(n, vector<double>(n, 0)); // matrice de transition
	double inv_n = 1. / ((double)n);
	double epsilon = epsilon_ * inv_n;
	double n_2 = (double)n / 2.;

	for (int i(0); i < n; ++i) // on parcourt la matrice de transitions
		for (int j(0); j < n; ++j)
			P[i][j] = inv_n * (1. + epsilon * ((double)preferences[j][i] - n_2));

	for (int i(0); i < n; ++i)
		P[i][i] = inv_n;

	double somme = 0; //calculer la stochasticité
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

	//modifier les P_ii, toujours pour la stochasticité
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

	somme = 0; //au cas où X n'est pas strictement de masse 1.
	for (int i(0); i < n; ++i)
		somme += X[i];
	somme_inv = 1. / somme;
	for (int i(0); i < n; ++i)
		X[i] *= somme_inv;

	return X;
}

int get_min(std::vector<std::vector<int>> votes, double epsilon, int d, bool asymetrique) {
	int n = votes[0].size();
	if (n == 1)
		return 0;

	std::vector<int> ordre;
	if (d == 0) {
		if (asymetrique == false)
			ordre = ordre_double(vote(votes, epsilon));
		else {
			ordre = ordre_double(vote(votes, epsilon));
			return ordre[n - 1];
		}
	}
	else
		ordre = get_max_min(votes, epsilon, d - 1, asymetrique);

	int i_max = ordre[0];

	for (int i(0); i < votes.size(); ++i)
		for (int j(0); j < n; ++j)
			if (votes[i][j] == i_max)
				votes[i].erase(votes[i].begin() + j);

	for (int i(0); i < votes.size(); ++i)
		for (int j(0); j < n; ++j)
			votes[i][j] = votes[i][j] > i_max ? votes[i][j] - 1 : votes[i][j];

	int i_min = get_min(votes, epsilon, d, asymetrique);
	i_min = i_min < i_max ? i_min : i_min + 1;

	return i_min;
}

std::vector<int> get_max_min(std::vector<std::vector<int>>  votes, double epsilon, int d, bool asymetrique) { //appeler asymétrique=false, d=0, epsilon=2. Renvoit une liste
	int n = votes[0].size();
	if (n == 1)
		return vector<int>(1, 0);

	int i_min = get_min(votes, epsilon, d, asymetrique);

	for (int i(0); i < votes.size(); ++i)
		for (int j(0); j < n; ++j)
			if (votes[i][j] == i_min)
				votes[i].erase(votes[i].begin() + j);

	for (int i(0); i < votes.size(); ++i)
		for (int j(0); j < n; ++j)
			votes[i][j] = votes[i][j] > i_min ? votes[i][j] - 1 : votes[i][j];

	vector<int> resultat = get_max_min(votes, epsilon, d, asymetrique);
	for (int i(0); i < resultat.size(); ++i)
		resultat[i] = resultat[i] < i_min ? resultat[i] : resultat[i] + 1;

	resultat.push_back(i_min);

	return resultat;
}

vector<int> ordre_alea(int n) {
	vector<double> liste(n, 0);
	for (int i(0); i < n; ++i)
		liste[i] = ((double)(rand()) / ((double)(RAND_MAX)));
	return ordre_double(liste);
}




//ajouter un vote en deux étapes : 
//d'abord les votant, ensuite les votés ...
// faire les tests !