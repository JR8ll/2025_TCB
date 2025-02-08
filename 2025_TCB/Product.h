#pragma once

#include <vector>

class Product {
private:
	int id;

	std::vector<int> route;
	std::vector<double> p;											// processing times
	std::vector<std::vector<std::pair<int, double>>> tcMaxBwd;		// time constraints from preceding steps [2nd stage]
	std::vector<std::vector<std::pair<int, double>>> tcMaxFwd;		// time cnostraints to successing steps [1st stage]

public:
	Product(int id, std::vector<int> route);

	int getId() const;
	int size() const;	// number of steps

	double getP(int stgIdx) const;

	void addTcMax(int first, int second, double duration);	// adds both backward and forward perspectives

};

