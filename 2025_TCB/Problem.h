#pragma once

#include<set>
#include<string>
#include<vector>

#include "Job.h"
#include "Product.h"
#include "Schedule.h"

class Job;
class Product;
class Schedule;

using pJob = std::unique_ptr<Job>;

class Problem {
private:
	std::string filename;

	std::vector<Product> products;
	std::vector<pJob> unscheduledJobs;

	unsigned seed;

	double omega;

	int n;		// number of jobs
	int stgs;	// number of stages/workcenters
	int F;		// number of products (families)

	std::vector<int> stages_1;	// stages (ids) with discrete processing
	std::vector<int> stages_b;	// stages (ids) with batching

	std::vector<int> m_o;		// number of machines at stage (workcenter)
	std::vector<int> m_B;		// batch capacity at stage (no discrimination of batching and discrete stages)

	std::vector<std::vector<double>> rm;	// release of machine [stg][mac]

	std::vector<std::vector<std::set<int>>> B_io;			// indices of preformed batches [product][stage]
	std::vector<std::vector<std::vector<int>>> B_iob;		// occupied capacity of preformed batches [product][stage][batch]
	std::vector<std::vector<std::vector<double>>> S_iob;	// start time of preformed batches [product][stage][batch]

	std::vector<double> jobs_d;							// due dates
	std::vector<double> jobs_r;							// release times
	std::vector<double> jobs_w;							// weights (priority)
	std::vector<int> jobs_f;							// products (families)
	std::vector<int> jobs_s;							// sizes

	std::vector<std::vector<int>> routes;				// indices of workstations in the routes [product]
	std::vector<std::vector<double>> pTimes;			// processing times [product][stage]

	std::vector<std::vector<std::vector<double>>> tc;	// time constraints [product][stage1][stage2]

public:
	Problem();
	Problem(std::string filename);

	std::string getFilename();

	int getN() const;
	int getStgs() const;
	int getF() const;

	void loadFromDat(std::string filename);

	std::pair<int, int> _tokenizeTupel(std::string tupel);

	std::unique_ptr<Schedule> getSchedule() const;	// machine environment + jobs to be scheduled

};

