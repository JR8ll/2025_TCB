#pragma once

#include<string>
#include<vector>

#include "Job.h"
#include "Product.h"
#include "Schedule.h"

class Job;
class Product;
class Schedule;

class Problem {
private:
	std::string filename;

	std::vector<Product> products;
	std::vector<Job> jobs;

	unsigned seed;

	int n;		// number of jobs
	int stgs;	// number of stages/workcenters
	int F;		// number of products (families)

public:
	Problem();
	Problem(std::string filename);

	std::string getFilename();

	int getN() const;
	int getStgs() const;
	int getF() const;

	void loadFromDat(std::string filename);

};

