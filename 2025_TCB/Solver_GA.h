#pragma once


#include <vector>
#include "Functions.h"

class Schedule;

struct GA_params {
	int nPop;							// population size
	double pElt;						// elite fraction of population
	double pRpM;						// fraction of population to be replaced by mutants
	double rhoe;						// probability that offspring inherit allele from elite parent
	unsigned K;							// number of independent populations
	int maxThreads;						// number of parallel threads
	bool applyLocalSearchBestFit;		// if true, local search follows best fit, if false first fit strategy
	double localSearchFraction;			// percentage of the elite, where the best x percent of the non-elite are decoded with local search (elite fitness values are just copied, not decoded)
	int localSearchEveryNGenerations;	// local search is only applied every ith generation (i = localSearchEveryNGenerations)
	int iterations;						// REPORTING: number of generations processed 
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

