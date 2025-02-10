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

	const std::vector<std::pair<int, double>>& getTcMaxBwd(size_t stgIdx) const;
	const std::vector<std::pair<int, double>>& getTcMaxFwd(size_t stgIdx) const;

	int getWorkcenterId(size_t stgIdx) const;
	double getP(size_t stgIdx) const;

	void setProcessingTimes(std::vector<double> pTimes);	// add processing times for all stages
	void addTcMax(int first, int second, double duration);	// adds both backward and forward perspectives

};

