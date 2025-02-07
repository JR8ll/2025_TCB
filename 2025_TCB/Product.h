#pragma once

#include <vector>

class Product {
private:
	int id;

	std::vector<double> p;		// processing times
	std::vector<double> tcMax;

public:
	int getId() const;

	double getP(int stgIdx) const;
};

