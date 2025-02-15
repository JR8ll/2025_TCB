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
double TCB::precision = 0.001;
mt19937 TCB::rng = mt19937(123456789);

// argv[1] filename of problem instance to be solved
// argv[2] seed for pseudo random-number generator
// argv[3] int describing the solving method to be used
// argv[4] time limit in seconds
// argv[5] console output on(=1)/off(=0)
// argv[6] filename of ga parameters
// argv[7] filename of decompMILP parameters

int main(int argc, char* argv[]) {

	// PROCESS COMMAND LINE ARGUMENTS
	TCB::logger = Logger();
	Problem p = Problem();
	TCB::prob = &p;
	int iSolver = -1;
	int iTilimSeconds = 3600;
	bool bConsole = false;
	string solverName = "n/a";
	string objectiveName = "TWT";	
	GA_params gaParams = GA_params();
	DECOMPMILP_params decompParams = DECOMPMILP_params();
	processCmd(argc, argv, iSolver, iTilimSeconds, bConsole, gaParams, decompParams);

	// PREPARE 
	pSched sched = TCB::prob->getSchedule();

	// SOLVE
	switch (iSolver) {
	case ALG_ITERATEDMILP:
		solverName = "DecompMILP";
		{
			Solver_MILP cplex = Solver_MILP(decompParams);
			vector<double> kappas = getDoubleGrid(0.1, 2.5, 0.1);
			cplex.solveDecompJobBasedDynamicSortingGridMILP(sched.get(), 4, iTilimSeconds, sortJobsByD, sortJobsByGATC);
		}
		break;
	case ALG_LISTSCHEDATC: 
		solverName = "ListSchedGATC";
		{
			sched->lSchedJobsWithSorting(sortJobsByGATC, decompParams);	// TODO: make pWait a parameter
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

	// RESULT SUMMARY (FILE OUTPUT)
	writeSolutions(sched.get(), solverName, objectiveName, iTilimSeconds, 666, &gaParams, &decompParams);	// TODO measure time
	

	// CONSOLE OUTPUT
	if (bConsole) {
		cout << "Solved using " << solverName << " in " << " seconds with TWT = " << sched->getTWT() << "." << endl;
		cout << *sched;
	}

	return EXIT_SUCCESS;
}