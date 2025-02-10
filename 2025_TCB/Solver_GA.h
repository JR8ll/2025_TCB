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
	GA_params* params;

public:
	Solver_GA(GA_params& params);
	~Solver_GA();

	double solveBRKGA_List_jobBased(Schedule& sched, int iTilimSeconds);

	bool hasCompleted();
	std::vector<double> getBestChromosome();

	static GA_params getDefaultParams();

};

