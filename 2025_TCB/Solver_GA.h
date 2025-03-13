#pragma once

#include<vector>
#include "Functions.h"

class Schedule;

struct GA_params {
	int nPop;				// population size
	double pElt;			// elite fraction of population
	double pRpM;			// fraction of population to be replaced by mutants
	double rhoe;			// probability that offspring inherit allele from elite parent
	unsigned K;				// number of independent populations
	int maxThreads;			// number of parallel threads
	int iterations;			// REPORTING: number of generations processed 
};

class Solver_GA {
private:
	bool completed;
	std::vector<double> bestChromosome;
	Sched_params* schedParams;
	GA_params* params;

public:
	Solver_GA(Sched_params& schedParams, GA_params& gaParams);
	~Solver_GA();

	double solveBRKGA_List_jobBased(Schedule& sched, int iTilimSeconds);

	bool hasCompleted();
	std::vector<double> getBestChromosome();

	static GA_params getDefaultParams();

};

