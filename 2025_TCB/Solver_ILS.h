#pragma once

#include <Windows.h>
#include "Common_aliases.h"
#include<chrono>

#include "GaDecoder.h"

class Schedule;

struct ILS_params {
	int nStarts;							// Multi-Start, number of new initializations
	int nPerturbationSteps;					// number of perturbations steps between two local search phases 
	bool applyBestFit;						// true: local search follows greedy/best fit strategy, false: local search follows first-fit strategy			
	bool randomizedLocalSearchSequence;
	double firstPhaseTimeLimitAllocation;	// percentage of time limit dedicated to a first phase (relevant for hybrid ILS)
	int multiStartIterations;				// REPORTING: number of new initializations (MULTI-START)
	std::vector<size_t> ilsIterations;		// REPORTING: number of iterations (local search + perturbation)
	int bestAfterSeconds;					// REPORTING: best solution found after ... seconds
};

struct ILS_Thread {
	double bestTWT = DBL_MAX;
	std::unique_ptr<Schedule> bestSched;
	double bestAfterSeconds = 0.0;
	std::vector<size_t> ilsIterations;
	int multiStartIterations = 0;
};

class Solver_ILS {
protected:
	Sched_params* schedParams;
	ILS_params* params;
public:
	Solver_ILS();
	Solver_ILS(Sched_params& schedParams, ILS_params& ilsParams);
	double solveILS(Schedule& sched, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, double pWait = 0.0);					// local search on schedule (insert job, swap job, merge batch)
	double solveILSparallelized(Schedule& sched, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, double pWait = 0.0);
	double solveILSonJobSequence(Schedule& sched, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, double pWait = 0.0);		// local search on sequence permutation (insert, swap)

	static void workerILS(DWORD coreIndex, std::unique_ptr<Schedule>& sched, Sched_params* schedParams, ILS_params* ilsParams,  initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, std::chrono::time_point<std::chrono::high_resolution_clock> start, ILS_Thread* localILS, double pWait = 0.0);
	static ILS_params getDefaultParams();

	static std::vector<DWORD> GetPCoreIndices();
};

class Solver_Sequence_ILS : public Solver_ILS {
private: 
	GaDecoderJobListSched decoder;	// means to decode a schedule from a given random key sequence
	std::vector<double> chr;		// chromosome = random key vector
	std::vector<double> bestChr;	
	double currentTWT;
	double bestTWT;

public: 
	Solver_Sequence_ILS(Schedule* schedule, Sched_params* schedParameters, ILS_params* ilsParameters);
	~Solver_Sequence_ILS();

	double solveILS(Schedule& sched, int iTilimSeconds);

	void insertJob(std::vector<double>& perm, size_t jobIdx, size_t posIdx);		// move random key from [jobIdx] to [posIdx]
	void swapJob(std::vector<double>& perm, size_t firstIdx, size_t secondIdx);		// swap keys at two indices
	double evaluateJobInsert(size_t jobIdx, size_t posIdx);
	double evaluateJobSwap(size_t firstIdx, size_t secondIdx);
	bool localSearchInsertJob(bool bestFit = true);
	bool localSearchSwapJob(bool bestFit = true);

	void perturbJobInsert();	// insert a randomly picked job to a randomly picked position
	void perturbJobSwap();		// swap two randomly chosen random keys
	void perturbRandomKey();	// reinitialize a random key at a randomly picked position

	void initRandomPermutation();

};

