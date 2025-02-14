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

struct ProbParams {
	unsigned seed;
	double omega;	// weighting factor for MILP objective function

	bool flowshop;	// true = flowshop, false = jobshop

	int n;			// number of jobs
	int stgs;		// number of stages (= number of workcenters in a job shop)
	int F;			// number of products
	int G;			// big integer

	double pReadyAtZero;	// percentage of jobs being ready at t=0
	double tcFlowFactor;	// time constraint flow factor >= 1.0

	std::pair<int, int> nTcInterval;	// lowest/highest number of time constraints per product
	int tcScenario;	// see Klemmt & Mönch

	std::pair<double, double> dueDateFF;	// lowest/highest factor multiplied with raw processing time to define due date (r_j + raw_processing_time * factor)

	std::pair<int, int> m_BIntervals;	// capacity per stage, dims: 1st Stage
	std::vector<int> m_BValues;			// capacity per stage, dims: 1st Stage (explicit alternative to m_BIntervals)
	std::pair<int, int> m_oIntervals;	// number of machines per stage, dims: 1st Stage

	std::vector<std::vector<int> > routes;	// indices of workstations, dims: 1st Product/Family

	// intervals first -> min, second -> max
	std::pair<double, double> pInterval;
	std::pair<double, double> wInterval;
	std::pair<double, double> rInterval;	// max will be multiplied with a factor u inside problem constructor
	std::pair<int, int> sInterval;
};

class Problem {
private:
	std::string filename;

	std::vector<Product> products;
	std::vector<pJob> unscheduledJobs;

	unsigned seed;

	double omega;

	bool flowshop;	// true = flowshop, false = jobshop

	int n;		// number of jobs
	int stgs;	// number of stages/workcenters
	int F;		// number of products (families)
	int G;		// big integer

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
	Problem(ProbParams& params);

	std::string getFilename();

	int getN() const;
	int getStgs() const;
	int getF() const;

	Product* getProduct(size_t productIdx);

	void loadFromDat(std::string filename);

	void _setG();	// set big integer

	std::pair<int, int> _tokenizeTupel(std::string tupel);

	std::unique_ptr<Schedule> getSchedule();	// machine environment + jobs to be scheduled

};

