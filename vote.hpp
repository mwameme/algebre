#pragma once
#include <algorithm>
#include <vector>

#include "TVD.hpp"
#include "ordre.hpp"

std::vector<double> vote(std::vector<std::vector<int>> const& votes, double epsilon);

std::vector<double> vote_iter(std::vector<std::vector<int>> const& votes, double epsilon,std::vector<double> const& poids);

std::vector<double> vote_n_iter(std::vector<std::vector<int>> const& votes, double epsilon,int iter);

std::vector<double> vote_experience(std::vector<std::vector<int>> votes_total, std::vector<std::vector<int>> votes_elus, double epsilon, int iter);


/*
class couple {
public:
	couple() {};

	double nombre;
	int index;
};

*/