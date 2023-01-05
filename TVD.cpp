//#include <Eigen/Dense>

#include "TVD.hpp"
#include <cmath>
//#include <Eigen/Geometry> 


using namespace std;
//using namespace Eigen;

int puissance2(int a) {
	int pow = 1;
	while (a > 0) {
		pow *= 2;
		--a;
	}
	return pow;
}

double distanceTVD(vector<double> const& pi, vector<double> const& X) {
	double somme = 0;
	int n = pi.size();
	for (int i(0); i < n; ++i)
		somme += abs(pi[i] - X[i]);
	somme *= 0.5;
	return somme;
}

double distanceChi2(vector<double> const& pi, vector<double> const& X) {
	double somme = 0;
	int n = pi.size();
	for (int i(0); i < n; ++i)
		somme += (pi[i] - X[i]) * (pi[i] - X[i]) / pi[i];
	return somme;
}

vector<double> multiplication(vector<double> const& X, vector<vector<double>> const& A) {
	int n = X.size();
	vector<double> Y(n, 0);
	for (int i(0); i < n; ++i) {
		double somme = 0;
		for (int j(0); j < n; ++j)
			somme += X[j] * A[j][i];
		Y[i] = somme;
	}

	return Y;
}

vector<vector<double>> multiplication(vector<vector<double>> const& A, vector<vector<double>> const& B) {
	int n = A.size();
	vector<vector<double>> C(n, vector<double>(n, 0.));
	for(int i(0);i<n;++i)
		for (int j(0); j < n; ++j) {
			double somme = 0;
			for (int k(0); k < n; ++k)
				somme += A[i][k] * B[k][j];
			C[i][j] = somme;
		}
	return C;
}


vector<vector<double>> puissance(vector<vector<double>> matrice, int puissance_n) {
	vector<vector<double>> resultat(matrice.size(), vector<double>(matrice.size(), 0));
	for (int i(0); i < matrice.size(); ++i)
		resultat[i][i] = 1;

	vector<vector<double>> puissance_m = matrice;
	while (puissance_n > 0) {
		if ((puissance_n % 2) == 1)
			resultat = multiplication(resultat, puissance_m);
		puissance_m = multiplication(puissance_m , puissance_m);
		puissance_n = puissance_n / 2;
	}
	return resultat;
}


int n_iter_distance(double epsilon, vector<vector<double>> matrice, vector<double> pi, function<double(vector<double>, vector<double>)> distance) {
	if (epsilon <= 0)
		return -1;
	int n = matrice.size();
	for (int i(0); i < n; ++i)
		if (matrice[i].size() != n)
			return -1;
	if (pi.size() != n)
		return -1;


	for (int i(0); i < n; ++i)
		for (int j(0); j < n; ++j)
			if ((matrice[i][j] < 0) || (matrice[i][j] > 1)) //stochastique
				return -1;

	for (int i(0); i < n; ++i) {
		double somme = 0;
		for (int j(0); j < n; ++j)
			somme += matrice[i][j];
		if ((somme < 1. - 0.0001) || (somme > 1. + 0.0001)) //stochastique
			return -1;
	}

	auto pi2 = multiplication(pi,matrice);
	double somme = distance(pi,pi2);

	if (abs(somme) > epsilon)
		return -1;

	vector<vector<vector<double>>> vec_M(0);
	vec_M.push_back(matrice);

	int n_iter = 0;
	bool continuer = true;
	auto matrice_temp = matrice;
	int puissance = 1;

	vector<vector<double>> vec_X(0);
	vec_X.reserve(n);
	for (int i(0); i < n; ++i) {
		vector<double> X (n,0);
		X [i] = 1;
		vec_X.push_back(X);
	}


	while (continuer) { //créer la liste des puissances
		++n_iter;
		if (n_iter > 25 )
			return -1;
		matrice_temp = multiplication(matrice_temp, matrice_temp);
		continuer = false;
		for (int i(0); i < vec_X.size(); ++i) {
			auto X = vec_X[i];
			X = multiplication(X,matrice_temp);
			somme = distance(pi,X);
//			for (int i(0); i < n; ++i)
//				somme += abs(X[i] - pi[i]);
			if (somme >= epsilon) {
				for (int j(0); j < i; ++j)
					vec_X.erase(vec_X.begin());
				continuer = true;
				break;
			}
		}
		if (continuer == false)
			break;
		puissance = puissance * 2;
		vec_M.push_back(matrice_temp);
	}

	if (puissance == 1)
		return 2;

	continuer = true;
	int n_max = vec_M.size() - 1;
	auto matrice_iter = vec_M[n_max];
	n_iter = n_max - 1;
	puissance = puissance2(n_max);
	while (continuer) {
		bool test = true;
		int n_test = n_iter;
		vector<vector<double>> matrice_locale;
		while (test) {
			bool test_local = true;
			matrice_locale = multiplication(matrice_iter,vec_M[n_test]);
			for (int i(0); i < vec_X.size(); ++i) {
				auto X = vec_X[i];
				X = multiplication(X,matrice_locale);
//				int j = 0;
				somme = distance(pi,X);
//				for (j = 0; i < n; ++j)
//					somme += abs(X[j] - pi[j]);
				if (somme > epsilon) {
					test_local = false;
					for (int j= 0; j < i; ++j)
						vec_X.erase(vec_X.begin());
					break;
				}
			}
			if (test_local == true) {
				n_test = n_test - 1;
				if (n_test < 0)
					break;
				test = true;
				continue;
			}
			else
				test = false;
		}
		if (n_test < 0)
			break;
		puissance = puissance + puissance2(n_test);
		matrice_iter = matrice_locale; // multiplication(matrice_iter, vec_M[n_test]);
		n_iter = n_test - 1;
		if (n_iter < 0)
			break;
	}


	puissance = puissance + 1;

	return puissance;
}