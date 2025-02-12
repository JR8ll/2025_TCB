#pragma once
#define NOMINMAX

#include "Common_aliases.h"
#include "Functions.h"

class IloIntTupleSet;
class IloIntArray;

class Job;
class Schedule;

class Solver_MILP {
public:
	Solver_MILP();
	~Solver_MILP();

	double solveJobBasedMILP(Schedule* schedule, int nDash = 10, int cplexTilim = 60);	// job based MILP by Cailloux & Mönch (MISTA 2019), iteratively considering max. nDash jobs, given sorting (static)
	double solveDecompJobBasedMILP(Schedule* schedule, int nDash = 10, int cplexTilim = 60, prioRuleKappa<pJob> rule = nullptr, double kappa = 1.0);	// like above with dynamic sorting (ATC-like priority rule)
	double solveDecompJobBasedDynamicSortingMILP(Schedule* schedule, int nDash = 10, int cplexTilim = 60, prioRule<pJob> initRule = sortJobsByD, prioRuleKappa<pJob> dynRule = sortJobsByGATC, double kappa = 1.0);
	double solveDecompJobBasedDynamicSortingMILP(Schedule* schedule, int nDash = 10, int cplexTilim = 60, prioRule<pJob> initRule = sortJobsByD, prioRuleKappa<pJob> dynRule = sortJobsByGATC, std::vector<double> kappas = {1.0});
};

