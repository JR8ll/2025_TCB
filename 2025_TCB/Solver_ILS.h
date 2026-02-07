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
	bool randomizedLocalSearchSequence;		// NOT YET USED (controls different sorting order for local search steps)
	double firstPhaseTimeLimitAllocation;	// percentage of time limit dedicated to a first phase (relevant for hybrid ILS)
	double secondPhaseRandomizedFraction;	// percentage of threads in the 2nd phase to be started with a random solution
	double seqLSsearchDepth;				// this value multiplied with n defines the number of maximal local search steps evaluated during sequence based local search ILS (phase 1 in hybrid ILS)

	int multiStartIterations;				// REPORTING: number of new initializations (MULTI-START)
	std::vector<size_t> ilsIterations;		// REPORTING: number of iterations (local search + perturbation)
	int bestAfterSeconds;					// REPORTING: best solution found after ... seconds
	double bestTWTAfterPhase1;				// REPORTING: for hybrid ILS - best TWT after phase 1 (sequence based ILS)
};

struct ILS_Thread {
	double bestTWT = DBL_MAX;
	std::unique_ptr<Schedule> bestSched;
	double bestAfterSeconds = 0.0;
	std::vector<size_t> ilsIterations;
	int multiStartIterations = 0;
};

struct ThreadParams {
	Sched_params localSchedParams;
	ILS_params localIlsParams;
};

struct ILSseq_Thread {
	double bestTWT = DBL_MAX;
	std::vector<double> bestChr;
	double bestAfterSeconds = 0.0;
	std::vector<size_t> ilsIterations;
	int multiStartIterations = 0;
};

// ILS operating on a schedule-object
class Solver_ILS {
protected:
	Sched_params* schedParams;
	ILS_params* params;
public:
	Solver_ILS();
	Solver_ILS(Sched_params& schedParams, ILS_params& ilsParams);
	double solveILS(Schedule& sched, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, double pWait = 0.0);													// local search on schedule (insert job, swap job, merge batch)
	double solveILSparallelized(Schedule& sched, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, double pWait = 0.0);
	double solveILSparallelized(Schedule& sched, initializerRK<pJob> init, const std::vector<double>& randomKeys, int iTilimSeconds, double pWait = 0.0);					// random key init
	double solveILSparallelized(Schedule& sched, initializerRK<pJob> init, const std::vector<std::vector<double>>& randomKeySets, int iTilimSeconds, double pWait = 0.0);	// parallel random key inits
	
	static void workILS(DWORD coreIndex, std::unique_ptr<Schedule>& sched, Sched_params* schedParams, ILS_params* ilsParams,  initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, std::chrono::time_point<std::chrono::high_resolution_clock> start, uint64_t threadSeed, ILS_Thread* localBest, double pWait = 0.0);
	static void workILSrk(DWORD coreIndex, std::unique_ptr<Schedule>& sched, Sched_params* schedParams, ILS_params* ilsParams, initializerRK<pJob> init, const std::vector<double>& randomKeys, int iTilimSeconds, std::chrono::time_point<std::chrono::high_resolution_clock> start, uint64_t threadSeed, ILS_Thread* localBest, double pWait = 0.0);
	static ILS_params getDefaultParams();

	static std::vector<DWORD> GetPCoreIndices(bool showWarning = true);
};

// ILS operating on a sequence of jobs (represented by a random key permutation), list-scheduling is used for decoding 
class Solver_Sequence_ILS : public Solver_ILS {
private: 
	std::vector<double> currentChr;
	std::vector<double> bestChr;
	std::vector<std::vector<double>> bestChrParallel;		// for parallel processing, save best sequence from every thread
	double currentTWT;
	double bestTWT;
	Schedule* masterSched;

public: 
	Solver_Sequence_ILS(Schedule* schedule, Sched_params* schedParameters, ILS_params* ilsParameters, int nCores);
	~Solver_Sequence_ILS();

	double solveILSseq(Schedule& sched, int iTilimSeconds);
	double solveILSseqParallelized(Schedule& sched, int iTilimSeconds);

	void insertJob(size_t jobIdx, size_t posIdx);		// move random key from [jobIdx] to [posIdx]
	void swapJob(size_t firstIdx, size_t secondIdx);		// swap keys at two indices
	double evaluateJobInsert(size_t jobIdx, size_t posIdx);
	double evaluateJobSwap(size_t firstIdx, size_t secondIdx);
	bool localSearchInsertJob(std::chrono::time_point<std::chrono::high_resolution_clock> finishBy, bool bestFit = true);
	bool localSearchSwapJob(std::chrono::time_point<std::chrono::high_resolution_clock> finishBy, bool bestFit = true);

	void insertJob(std::vector<double>& perm, size_t jobIdx, size_t posIdx);		// move random key from [jobIdx] to [posIdx]
	void swapJob(std::vector<double>& perm, size_t firstIdx, size_t secondIdx);		// swap keys at two indices
	double evaluateJobInsert(std::vector<double>& chromosome, size_t jobIdx, size_t posIdx);
	double evaluateJobSwap(std::vector<double>& chromosome, size_t firstIdx, size_t secondIdx);
	bool localSearchInsertJob(std::vector<double>& chromosome, std::chrono::time_point<std::chrono::high_resolution_clock> finishBy, bool bestFit = true);
	bool localSearchSwapJob(std::vector<double>& chromosome, std::chrono::time_point<std::chrono::high_resolution_clock> finishBy, bool bestFit = true);

	std::vector<double> getBestChr();
	std::vector<std::vector<double>> getBestChrParallel();

	void perturbJobInsert();									// insert a randomly picked job to a randomly picked position
	void perturbJobSwap();										// swap two randomly chosen random keys
	void perturbRandomKey();									// reinitialize a random key at a randomly picked position
	void perturbJobInsert(std::vector<double>& chromosome);		// insert a randomly picked job to a randomly picked position
	void perturbJobSwap(std::vector<double>& chromosome);		// swap two randomly chosen random keys
	void perturbRandomKey(std::vector<double>& chromosome);		// reinitialize a random key at a randomly picked position

	void decode(Schedule* sched, std::vector<double>& chromosome);
	void formMasterSchedule(std::vector<double>& chromosome);
	double decodeAndGetTWT(std::vector<double>& chr);
	static double staticDecodeAndGetTWT(Schedule* sched, std::vector<double>& chr);

	void initRandomPermutation();
	static void initRandomPermutation(std::vector<double>& chromosome);

	static void workILSseq(DWORD coreIndex, std::unique_ptr<Schedule>& sched, Sched_params schedParams, ILS_params ilsParams, int iTilimSeconds, std::chrono::time_point<std::chrono::high_resolution_clock> start, uint64_t threadSeed, ILSseq_Thread* localILS);
};

class Solver_Hybrid_ILS {
private:
	Solver_Sequence_ILS phase1;
	Solver_ILS phase2;
	Sched_params* schedParams;
	ILS_params* ilsParams;
public:
	Solver_Hybrid_ILS(Schedule* schedule, Sched_params* schedParameters, ILS_params* ilsParameters, int nCores);
	~Solver_Hybrid_ILS();

	double solveILShybrid(Schedule& sched, int iTilimTotal);
};


