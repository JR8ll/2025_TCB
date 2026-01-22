#pragma once

#include <Windows.h>
#include "Common_aliases.h"
#include<chrono>

class Schedule;

struct ILS_params {
	int nStarts;						// Multi-Start, number of new initializations
	int nPerturbationSteps;				// number of perturbations steps between two local search phases 
	bool applyBestFit;					// true: local search follows greedy/best fit strategy, false: local search follows first-fit strategy			
	bool randomizedLocalSearchSequence;
	int multiStartIterations;			// REPORTING: number of new initializations (MULTI-START)
	std::vector<size_t> ilsIterations;		// REPORTING: number of iterations (local search + perturbation)
	int bestAfterSeconds;				// REPORTING: best solution found after ... seconds
};

struct ILS_Thread {
	double bestTWT = DBL_MAX;
	std::unique_ptr<Schedule> bestSched;
	double bestAfterSeconds = 0.0;
	std::vector<size_t> ilsIterations;
	int multiStartIterations = 0;
};

class Solver_ILS {
private:
	Sched_params* schedParams;
	ILS_params* params;
public:
	Solver_ILS(Sched_params& schedParams, ILS_params& ilsParams);
	double solveILS(Schedule& sched, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, double pWait = 0.0);
	double solveILSparallelized(Schedule& sched, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, double pWait = 0.0);

	static void workerILS(DWORD coreIndex, std::unique_ptr<Schedule>& sched, Sched_params* schedParams, ILS_params* ilsParams,  initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, std::chrono::time_point<std::chrono::high_resolution_clock> start, ILS_Thread* localILS, double pWait = 0.0);
	static ILS_params getDefaultParams();

	static std::vector<DWORD> GetPCoreIndices();

};

