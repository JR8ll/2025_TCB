#include<chrono>
#include<thread>
#include<vector>

#include "MTRand.h"	// brkgaAPI
#include "BRKGA.h"	// brkgaAPI

#include "Solver_GA.h"
#include "Functions.h"
#include "GaDecoder.h"
#include "Schedule.h"

using namespace std;

Solver_GA::Solver_GA(Sched_params& schedParameters, GA_params& parameters) : completed(false) {
	completed = false;// : completed(false) {
	bestChromosome = vector<double>();
	schedParams = &schedParameters;
	params = &parameters;
}
Solver_GA::~Solver_GA() {}

double Solver_GA::solveBRKGA_List_jobBased(Schedule& sched, int iTilimSeconds) {
	auto start = chrono::high_resolution_clock::now();
	chrono::seconds usedTime;
	chrono::time_point<chrono::high_resolution_clock> stop;

	GaDecoderJobLs decoder(&sched, schedParams);
	const long unsigned rngSeed = TCB::seed;
	MTRand rng(rngSeed);

	int processor_count = std::thread::hardware_concurrency();
	const unsigned MAXT = max(1, processor_count);

	params->maxThreads = MAXT;

	BRKGA<GaDecoderJobLs, MTRand> algorithm((int)sched.getN(), params->nPop, params->pElt, params->pRpM, params->rhoe, decoder, rng, params->K, MAXT);

	int iterationCounter = 0;
	do {
		algorithm.evolve();
		++iterationCounter;
		stop = chrono::high_resolution_clock::now();
		usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
	} while (usedTime.count() < iTilimSeconds);

	TCB::logger.Log(Info, "Solver_GA::solveBRKGA_List_jobBased finished after " + to_string(iterationCounter) + " iterations.");
	decoder.formSchedule(algorithm.getBestChromosome());
	params->iterations = iterationCounter;
	completed = true;
	bestChromosome = algorithm.getBestChromosome();
	return sched.getTWT();
}

bool Solver_GA::hasCompleted() { return completed; }
std::vector<double> Solver_GA::getBestChromosome() { return bestChromosome; }

GA_params Solver_GA::getDefaultParams() {
	GA_params gaParams = GA_params();
	gaParams.nPop = 100;
	gaParams.pElt = 0.2;
	gaParams.pRpM = 0.1;
	gaParams.rhoe = 0.7;
	gaParams.K = 3;
	int processor_count = std::thread::hardware_concurrency();
	gaParams.maxThreads = max(1, processor_count);
	return gaParams;
}


