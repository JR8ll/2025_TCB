#pragma once

#include "Common_aliases.h"

class Schedule;

struct ILS_params {
	int nStarts;						// Multi-Start, number of new initializations
	int nPerturbationSteps;				// number of perturbations steps between two local search phases 
	bool applyBestFit;					// true: local search follows greedy/best fit strategy, false: local search follows first-fit strategy			
	bool randomizedLocalSearchSequence;
	int iterations;						// REPORTING: number of generations processed 
};

class Solver_ILS {
private:
	Sched_params* schedParams;
	ILS_params* params;
public:
	Solver_ILS(Sched_params& schedParams, ILS_params& ilsParams);
	double solveILS(Schedule& sched, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, double pWait = 0.0);

	static ILS_params getDefaultParams();

};

