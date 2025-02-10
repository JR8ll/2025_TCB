#pragma once

#include<vector>

class Schedule;

struct GA_params {
	int nPop;				// population size
	double pElt;			// elite fraction of population
	double pRpM;			// fraction of population to be replaced by mutants
	double rhoe;			// probability that offspring inherit allele from elite parent
	unsigned K;				// number of independent populations
	int maxThreads;			// number of parallel threads
};

class Solver_GA {
private:
	bool completed;
	std::vector<double> bestChromosome;

public:
	Solver_GA();
	~Solver_GA();

	double solveBRKGA_List_jobBased(Schedule& sched, int iTilimSeconds, GA_params& params);

	bool hasCompleted();
	std::vector<double> getBestChromosome();

	GA_params getDefaultParams();

};

