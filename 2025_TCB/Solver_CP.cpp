#include<ilcp/cp.h>
#include<map>

#include "Solver_CP.h"
#include "Functions.h"
#include "Problem.h"
#include "Schedule.h"

ILOSTLBEGIN

using IloCumulFctArray2 = IloArray <IloArray<IloCumulFunctionExpr>>;
using IloStateFctArray2 = IloArray<IloArray<IloStateFunction>>;

Solver_CP::Solver_CP() {}
Solver_CP::~Solver_CP() {}

double Solver_CP::solveCP(Schedule* schedule, int nDash, int tilim) {
	if(schedule->getProblem() == nullptr) throw ExcSched("solveDecompositionMILP missing problem reference");
	vector<pJob> consideredJobs = vector<pJob>();
	for (size_t j = 0; j < nDash; ++j) {
		try {
			consideredJobs.push_back(move(schedule->get_pJob(0)));
		}
		catch (const out_of_range& ex) {
			break;	// all jobs considered
		}
	}

	IloEnv env;

	// PARAMETERS
	IloNum G = 9;	// TODO set omega
	IloInt m = schedule->size();

	int nProducts = schedule->getProblem()->getF();
	int nStages = schedule->size();
	int nJobs = consideredJobs.size();
	int nBatches = nJobs + schedule->getNumberOfScheduledJobs();

	IloIntTupleSet stageMachine(env, 2);
	for (int o = 0; o < nStages; ++o) {
		for (int l = 0; l < (*schedule)[o].size(); ++l) {
			stageMachine.add(IloIntArray(env, 2, o, l));
		}
	}

	map<pair<int, int>, int> mapStgMac = map<pair<int, int>, int>();	// for access to tuple "stageMachine"
	int index = 0;
	for (IloIntTupleSetIterator it(env, stageMachine); it.ok(); ++it) {
		mapStgMac.insert(make_pair(make_pair((*it)[0], (*it)[1]), index));
		index++;
	}
	const map<pair<int, int>, int> constMapStgMac = mapStgMac;

	IloIntArray B(env);	// capacity[stage]
	for (int o = 0; o < m; ++o) {
		int tempCap = schedule->getCapAtStageIdx(o);	// assumption: parallel identical machines
		if (tempCap > 1) {
			B.add(tempCap);
		}
	}

	IloCumulFctArray2 loadOfMachine;
	IloStateFctArray2 stateOfMachine;
	// CUMUL/STATE FUNCTIONS
	for (int o = 0; o < m; ++o) {
		loadOfMachine.add(IloArray<IloCumulFunctionExpr>(env));
		stateOfMachine.add(IloArray<IloStateFunction>(env));
		for (int j = 0; j < nJobs; ++j) {
			loadOfMachine[o].add(IloCumulFunctionExpr(env));
			stateOfMachine[o].add(IloStateFunction(env));
		}
	}

	// INTERVALS
	IloArray<IloArray<IloIntervalVarArray>> ops;					// interval [Products][Jobs][Stages]
	IloArray<IloArray<IloArray<IloIntervalVarArray>>> opsOnMachine;	// set if operation o of job i is scheduled on machine m
	for (int f = 0; f < nProducts; ++f) {
		ops.add(IloArray<IloIntervalVarArray>(env));
		opsOnMachine.add(IloArray<IloArray<IloIntervalVarArray>>(env));
		for (int j = 0; j < nJobs; ++j) {
			ops[f].add(IloIntervalVarArray(env));
			opsOnMachine[f].add(IloArray<IloIntervalVarArray>(env));
			for (int o = 0; o < nStages; ++o) {
				int processingTime = 0;
				if (consideredJobs[j]->getF() == (f + 1)) {
					processingTime = consideredJobs[j]->getP(o);
				}
				ops[f][j].add(IloIntervalVar(env, processingTime));
				opsOnMachine[f][j].add(IloIntervalVarArray(env));
				for (int k = 0; k < schedule[o].size(); ++k) {
					opsOnMachine[f][j][o].add(IloIntervalVar(env));
					opsOnMachine[f][j][o][k].setOptional();
					loadOfMachine[o][k] += IloPulse(opsOnMachine[f][j][o][k], 1);
				}
			}
		}
	}
}
