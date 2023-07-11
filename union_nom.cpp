#include "union_nom.hpp"


std::vector<std::string> union_nom(std::vector<std::string > const& gauche, std::vector<std::string> const& droit) {
	std::vector<std::string> result = gauche;
	for (int i(0); i < droit.size(); ++i) {
		bool trouve = false;
		for (int j(0); j < gauche.size(); ++j)
			if (droit[i] == gauche[j]) {
				trouve = true;
				break;
			}
		if (trouve)
			continue;
		result.push_back(droit[i]);
	}
	return result;
};
