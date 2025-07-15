#pragma once

#include<ilcp/cp.h>

class Schedule;

using IloCumulFctArray2 = IloArray <IloArray<IloCumulFunctionExpr>>;
using IloStateFctArray2 = IloArray<IloArray<IloStateFunction>>;

class Solver_CP {
private:
	IloEnv env;
	IloCumulFctArray2 loadOfMachine;
	IloStateFctArray2 stateOfMachine;
	IloArray<IloArray<IloIntervalVarArray>> ops;
	IloArray<IloArray<IloArray<IloIntervalVarArray>>> opsOnMachine;

	Schedule* sched;
	void initModel();		// parameters and decision variables (intervals)
	void setConstraints();	
	void updateModel();		// update constraints to partially generated schedule

public:
	Solver_CP();
	~Solver_CP();

	double solveCP(Schedule* schedule, int nDash = 5, int tilim = 60);
};