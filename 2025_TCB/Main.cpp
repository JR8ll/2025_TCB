#include <iostream>

#include "Functions.h"
#include "Problem.h"
#include "Solver_GA.h"

using namespace std;
using pSched = unique_ptr<Schedule>;

Problem* TCB::prob;
LogPriority Logger::verbosity = Info;
const char* Logger::filepath = "log.txt";
Logger TCB::logger;
unsigned TCB::seed = 123456789;
double TCB::precision = 0.001;
mt19937 TCB::rng = mt19937(123456789);

static const int ALG_ITERATEDMILP = 1;		// iterated MILP solving
static const int ALG_LISTSCHEDATC = 2;		// simple List scheduling approach
static const int ALG_BRKGALISTSCH = 3;		// biased random-key ga with a list scheduling decoder
static const int ALG_BRKGALS2MILP = 4;		// get best sequence from brkga, then iteratively apply ops from this sequence to MILP

// argv[1] filename of problem instance to be solved
// argv[2] seed for pseudo random-number generator
// argv[3] int describing the solving method to be used
// argv[4] time limit in seconds
// argv[5] console output on(=1)/off(=0)
// argv[6] filename of ga parameters

int main(int argc, char* argv[]) {
	// PROCESS COMMAND LINE ARGUMENTS
	TCB::logger = Logger();
	Problem p = Problem();
	TCB::prob = &p;
	int iSolver = -1;
	int iTilimSeconds = 3600;
	bool bConsole = false;
	string solverName = "n/a";
	GA_params gaParams = GA_params();
	processCmd(argc, argv, iSolver, iTilimSeconds, bConsole, gaParams);

	// PREPARE 
	pSched sched = TCB::prob->getSchedule();

	// SOLVE
	switch (iSolver) {
	case ALG_ITERATEDMILP:
		solverName = "DecompMILP";
		break;
	case ALG_LISTSCHEDATC: 
		solverName = "ListSchedGATC";
		{
			vector<double> kappas = getDoubleGrid(0.1, 2.5, 0.1);
			sched->lSchedJobsWithSorting(sortJobsByGATC, kappas);
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
	double twt = sched->getTWT();

	// CONSOLE OUTPUT
	if (bConsole) {
		cout << "Solved using " << solverName << " in " << " seconds with TWT = " << twt << "." << endl;
		cout << *sched;
	}

	return EXIT_SUCCESS;
}