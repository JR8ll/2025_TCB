#include<map>

#include "Solver_CP.h"
#include "Functions.h"
#include "Problem.h"
#include "Schedule.h"

ILOSTLBEGIN

Solver_CP::Solver_CP() : sched(nullptr) {
	env = IloEnv();
}
Solver_CP::~Solver_CP() {}

double Solver_CP::solveCP(Schedule* schedule, int nDash, int tilim) {
	if(schedule->getProblem() == nullptr) throw ExcSched("solveDecompositionMILP missing problem reference");
	sched = schedule;
	initModel();

	vector<pJob> consideredJobs = vector<pJob>();
	for (size_t j = 0; j < nDash; ++j) {
		try {
			consideredJobs.push_back(move(schedule->get_pJob(0)));
		}
		catch (const out_of_range& ex) {
			break;	// all jobs considered
		}
	}

	
}

void Solver_CP::initModel() {
	int m = sched->size();
	int nJobs = sched->getN();
	for (int o = 0; o < m; ++o) {
		loadOfMachine.add(IloArray<IloCumulFunctionExpr>(env));
		stateOfMachine.add(IloArray<IloStateFunction>(env));
		for (int j = 0; j < nJobs; ++j) {
			loadOfMachine[o].add(IloCumulFunctionExpr(env));
			stateOfMachine[o].add(IloStateFunction(env));
		}
	}


}

void Solver_CP::setConstraints() {
}

void Solver_CP::updateModel() {
}
