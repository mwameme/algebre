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

int n_iter_distance(double const epsilon, vector<vector<double>> matrice, vector<double> pi, function<double(vector<double>, vector<double>)> distance) {
	int n = matrice.size();


	vector<vector<vector<double>>> vec_M(0);
	vec_M.push_back(matrice);

	int n_iter = 0;
	bool continuer = true;
	auto matrice_temp = matrice;
	int puissance = 1;
	double somme;

	vector<vector<double>> vec_X(0);
	vec_X.reserve(n);
	for (int i(0); i < n; ++i) {
		vector<double> X(n, 0.);
		X[i] = 1.;
		vec_X.push_back(X);
	}


	while (continuer) { //créer la liste des puissances
		++n_iter;
		if (n_iter > 25)
			throw std::domain_error("n arrive pas à converger : n_iter_distance");
		continuer = false;
		for (int i(0); i < vec_X.size(); ++i) {
			auto X = vec_X[i];
			X = multiplication(X, matrice_temp);
			somme = distance(pi, X);
			//			for (int i(0); i < n; ++i)
			//				somme += abs(X[i] - pi[i]);
			if (somme >= epsilon) {
				if(i>0)
					vec_X.erase(vec_X.begin(),vec_X.begin() + (i-1));
				continuer = true;
				break;
			}
		}
		if (! continuer)
			break;
		puissance = puissance * 2;
		matrice_temp = multiplication(matrice_temp, matrice_temp);
		vec_M.push_back(matrice_temp);
	}

	if (puissance == 1) {
		return 1;
	}

	vec_M.pop_back();

	continuer = true;
	n_iter = vec_M.size() - 1;
	auto matrice_iter = vec_M[n_iter];
	puissance = puissance2(n_iter);


	while (true) {
		vector<vector<double>> matrice_locale(0);
		while (true) {
			bool test_local = true;
			matrice_locale = multiplication(matrice_iter, vec_M[n_iter]);
			for (int i(0); i < vec_X.size(); ++i) {
				vector<double> X = vec_X[i];
				X = multiplication(X, matrice_locale);
				somme = distance(pi, X);
				if (somme > epsilon) {
					test_local = false;
					if(i>0)
						vec_X.erase(vec_X.begin(),vec_X.begin() + (i-1));
					break;
				}
			}
			if (test_local) {
				n_iter -= 1;
				if (n_iter < 0)
					break;
				continue;
			}
			else
				break;
		}
		if (n_iter < 0)
			break;
		puissance = puissance + puissance2(n_iter);
		matrice_iter = matrice_locale;
	}


	puissance = puissance + 1;

	return puissance;
}
