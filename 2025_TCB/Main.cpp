#include <chrono>
#include <iostream>

#include "Functions.h"
#include "Problem.h"
#include "Solver_GA.h"
#include "Solver_MILP.h"

using namespace std;
using pSched = unique_ptr<Schedule>;

Problem* TCB::prob;
LogPriority Logger::verbosity = Info;
const char* Logger::filepath = "log.txt";
Logger TCB::logger;
unsigned TCB::seed = 123456789;
double TCB::precision = 0.1;
mt19937 TCB::rng = mt19937(123456789);

// argv[1] filename of problem instance to be solved
// argv[2] seed for pseudo random-number generator
// argv[3] int describing the solving method to be used: 1: Decomposition-MILP, 2: List Scheduling, 3: BRKGA, 4: BRKGA + final Decomposition-MILP with best sorting order (see constants in Functions.h)
// argv[4] time limit in seconds
// argv[5] console output on(=1)/off(=0)
// argv[6] filename of scheduling parameters
// argv[7] filename of ga parameters
// argv[8] filename of decompMILP parameters

int main(int argc, char* argv[]) {

	Problem::genInstancesTCB25_Feb25_exact();

	// PROCESS COMMAND LINE ARGUMENTS
	TCB::logger = Logger();
	Problem p = Problem();
	TCB::prob = &p;
	int iSolver = -1;
	int iTilimSeconds = 3600;
	bool bConsole = false;
	string solverName = "n/a";
	string objectiveName = "TWT";
	Sched_params schedParams = Sched_params();
	GA_params gaParams = GA_params();
	DECOMPMILP_params decompParams = DECOMPMILP_params();
	processCmd(argc, argv, iSolver, iTilimSeconds, bConsole, schedParams, gaParams, decompParams);

	// PREPARE 
	pSched sched = TCB::prob->getSchedule();

	// START TIME MEASUREMENT
	auto start = chrono::high_resolution_clock::now();
	chrono::seconds usedTime;
	chrono::time_point<chrono::high_resolution_clock> stop;

	// SOLVE
	switch (iSolver) {
	case ALG_ITERATEDMILP:
		solverName = "DecompMILP";
		{
			Solver_MILP cplex = Solver_MILP(schedParams, decompParams);
			cplex.solveDecompJobBasedDynamicSortingGridMILP(sched.get(), decompParams.nDash, iTilimSeconds, sortJobsByD, sortJobsByGATC);
		}
		break;
	case ALG_LISTSCHEDATC: 
		solverName = "ListSchedGATC";
		{
			sched->lSchedJobsWithSorting(sortJobsByGATC, schedParams);	
		}
		break;
	case ALG_BRKGALISTSCH:
		solverName = "BRKGA"; 
		{
			Solver_GA brkga = Solver_GA(gaParams);
			brkga.solveBRKGA_List_jobBased(*sched.get(), iTilimSeconds);
		}
		break;
	case ALG_BRKGALS2MILP:
		solverName = "BRKGA_MILP";
		break;
	default:
		TCB::logger.Log(Warning, "Program was executed with no valid algorithm key");
	} 

	// STOP TIME MEASUREMENT
	stop = chrono::high_resolution_clock::now();
	usedTime = chrono::duration_cast<chrono::seconds>(stop - start);

	// RESULT SUMMARY (FILE OUTPUT)
	writeSolutions(sched.get(), iSolver, solverName, objectiveName, iTilimSeconds, usedTime.count(), &schedParams, &gaParams, &decompParams);	// TODO measure time
	

	// CONSOLE OUTPUT
	if (bConsole) {
		cout << "Solved using " << solverName << " in " << " seconds with TWT = " << sched->getTWT() << "." << endl;
		cout << *sched;
	}

	return EXIT_SUCCESS;
}