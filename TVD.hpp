#pragma once
#include <vector>
#include <cmath>
#include <functional>
#include <exception>
#include <stdexcept>

//#include <iostream>


double distanceTVD(std::vector<double> const& , std::vector<double> const& );

double distanceChi2(std::vector<double> const&, std::vector<double> const&);

std::vector<double> multiplication(std::vector<double> const&, std::vector<std::vector<double>> const&);

std::vector<std::vector<double>> multiplication(std::vector<std::vector<double>> const&, std::vector<std::vector<double>> const&);

std::vector<std::vector<double>> puissance(std::vector<std::vector<double>> , int );


int n_iter_distance(double epsilon, std::vector<std::vector<double>> matrice, std::vector<double> pi,std::function<double(std::vector<double>,std::vector<double>)>);


int puissance2(int);
