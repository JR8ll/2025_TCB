#pragma once

#include<ilcp/cp.h>
#include "Job.h"
#include <vector>
using std::vector;

struct BatchMock { // Find a better way
	int product;
	int stage;
	int machine;
	int start;
	int capacity;

	BatchMock(int _product, int _stage, int _machine, int _start, int _capacity) {
		product = _product;
		stage = _stage;
		machine = _machine;
		start = _start;
		capacity = _capacity;
	}
};

class Schedule;

using IloCumulFctArray2 = vector<vector<IloCumulFunctionExpr>>;
using IloStateFctArray2 = vector<vector<IloStateFunction>>;
using IloInterval2 = vector<vector<IloIntervalVarArray>>;
using IloInterval3 = vector<vector<vector<IloIntervalVarArray>>>;
using BatchVector = vector<vector<vector<BatchMock>>>;
using TimeVector = vector<vector<vector<double>>>;
using MachineVector = vector<vector<vector<int>>>;

class Solver_CP {
private:
	int decompStep = 0;		// track decompostion step
	int solveCount = 0;		// used for seed
	int solveLimit = 4;		// prevents infinite loops

	int timeLimit;
	bool modelCreated;
	Sched_params* schedParams;

	IloEnv env;
	IloModel model;
	IloCP cp;

	IloCumulFctArray2 loadOfMachine;
	IloStateFctArray2 stateOfMachine;
	IloInterval2 ops;
	IloInterval3 opsOnMachine;

	vector<pJob> consideredJobs;

	Schedule* sched;
	void initModel();															// parameters and decision variables (intervals)
	void setConstraints(vector<pJob>& consideredJobs, BatchVector& oldBatches = BatchVector());
	void updateModel(TimeVector& starts, TimeVector& ends, MachineVector& machines, BatchVector& oldBatches, double& obj);									// update constraints to partially generated schedule (cp solution)
	void staticUpdate(TimeVector& starts, TimeVector& ends, MachineVector& machines, BatchVector& oldBatches);												// update constraints to partially generated schedule (list scheduling)

public:
	Solver_CP(Sched_params& sched_params);
	~Solver_CP();

	bool solveCP(int tilim);																								 // 3. CP Solver does Iteration Step has a Solution?
	double solveStepCP(Schedule* schedule, BatchVector& oldBatches, int nDash = 5, int tilim = 60);							 // 2. Iteration Solver
	double solveDecompCP(Schedule* schedule, Sched_params* sched_params, int nDash = 5, int tilim = 60, int method = 1);	 // 1. Entry Point

	void solveLNS(Schedule* schedule);
	void setStartingPoint();
	void setListSchedWarmStart(Schedule& sched);	// JR
	void setFixPartialSolution(bool fixMachineAssignments, bool fixMachineSequences);
};