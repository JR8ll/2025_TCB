#include<vector>

#include "MTRand.h"	// brkgaAPI
#include "BRKGA.h"	// brkgaAPI

#include "Solver_GA.h"

using namespace std;

Solver_GA::Solver_GA() : completed(false) {
	bestChromosome = vector<double>();
}
Solver_GA::~Solver_GA() {}

double Solver_GA::solveBRKGA_List_jobBased(Schedule& sched, int iTilimSeconds)
{
	// TODO
	return 0.0;
}

bool Solver_GA::hasCompleted() { return completed; }
std::vector<double> Solver_GA::getBestChromosome() { return std::vector<double>(); }


