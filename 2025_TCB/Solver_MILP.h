#pragma once
#define NOMINMAX

#include "Common_aliases.h"

class IloIntTupleSet;
class IloIntArray;

class Job;
class Schedule;

class Solver_MILP {
public:
	Solver_MILP();
	~Solver_MILP();

	double solveDecompositionMILP(Schedule* schedule, int nDash = 10, int cplexTilim = 60);	// job based MILP by Cailloux & Mönch (MISTA 2019), iteratively considering max. nDash jobs, given sorting (static)
	//double solveDecompositionMILP(Schedule* schedule, prioRuleKappaT<Job> rule,int nDash = 10, int cplexTilim = 60);	// like above with dynamic sorting (ATC-like priority rule)


};

