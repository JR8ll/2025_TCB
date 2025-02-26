#pragma once

#include "Common_aliases.h"

class Schedule;

class Solver_ILS {
private:
	Sched_params* schedParams;
public:
	Solver_ILS(Sched_params& schedParams);
	double solveILS(Schedule& sched, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, double pWait = 0.0);
};

