#pragma once
#define NOMINMAX

#include "Common_aliases.h"
#include "Functions.h"

class IloIntTupleSet;
class IloIntArray;

class Job;
class Schedule;

static const int DECOMP_SOLVER_MILP = 1;		
static const int DECOMP_SOLVER_CP = 2;

static const int PRIORULE_EDD = 1;	// earliest due date
static const int PRIORULE_FCFS = 2;	// first come first served (release)

struct DECOMPMILP_params {
	int nDash;					// number of jobs being considered in each iteration
	int cplexTilim;				// time limit for MILP/CP solver
	int method;					// 1: MILP, 2: CP
	int initPrioRule;			// 1: EDD, 2: FCFS, ...
	double kappaLow;			// grid of kappa values considered for GATC-sorting
	double kappaHigh;
	double kappaStep;
};

class Solver_MILP {
private:
	DECOMPMILP_params* params;

public:
	Solver_MILP(DECOMPMILP_params& decompParams);
	~Solver_MILP();

	double solveJobBasedMILP(Schedule* schedule, int nDash = 10, int cplexTilim = 60);	// job based MILP by Cailloux & Mönch (MISTA 2019), iteratively considering max. nDash jobs, given sorting (static)
	double solveDecompJobBasedMILP(Schedule* schedule, int nDash = 10, int cplexTilim = 60, prioRuleKappa<pJob> rule = nullptr, double kappa = 1.0);	// like above with dynamic sorting (ATC-like priority rule)
	double solveDecompJobBasedDynamicSortingMILP(Schedule* schedule, int nDash = 10, int cplexTilim = 60, prioRule<pJob> initRule = sortJobsByD, prioRuleKappa<pJob> dynRule = sortJobsByGATC, double kappa = 1.0);
	double solveDecompJobBasedDynamicSortingGridMILP(Schedule* schedule, int nDash = 10, int cplexTilim = 60, prioRule<pJob> initRule = sortJobsByD, prioRuleKappa<pJob> dynRule = sortJobsByGATC);	// kappa grid from internal class member params

	static DECOMPMILP_params getDefaultParams();
};

