#include <iostream>

#include "Functions.h"
#include "Problem.h"

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
static const int ALG_LISTSCHEDWTG = 2;		// simple List scheduling approach
static const int ALG_BRKGALISTSCH = 3;		// biased random-key ga with a list scheduling decoder
static const int ALG_BRKGALS2MILP = 4;		// get best sequence from brkga, then iteratively apply ops from this sequence to MILP

int main(int argc, char* argv[]) {
	TCB::logger = Logger();

	Problem p = Problem();
	TCB::prob = &p;

	int iSolver = -1;
	int iTilimSeconds = 3600;

	// PROCESS COMMAND LINE ARGUMENTS
	if (argc > 1) {
		TCB::prob->loadFromDat(argv[1]);
	}
	else {
		TCB::prob->loadFromDat("debug_data.dat");
	}

	if (argc > 2) {
		TCB::seed = atoi(argv[2]);
		TCB::rng = mt19937(TCB::seed);	// global random engine (mersenne twister)
	}

	if (argc > 3) {
		iSolver = atoi(argv[3]);
	}
	if (argc > 4) {
		iTilimSeconds = atoi(argv[4]);
	}

	switch (iSolver) {
	case ALG_ITERATEDMILP:
		break;
	case ALG_LISTSCHEDWTG:
		break;
	case ALG_BRKGALISTSCH:
		break;
	case ALG_BRKGALS2MILP:
		break;
	default:
		TCB::logger.Log(Warning, "Program was executed with no valid algorithm key");
	} 


	// ************ DEBUGGING *****************
	pSched sched = TCB::prob->getSchedule();
	sched->lSchedJobsWithSorting(sortJobsByGATC, 1.5);
	cout << *sched;
	cout << "TWT = " << sched->getTWT();

	sched->reset();
	vector<double> kappas = getDoubleGrid(0.1, 2.5, 0.1);
	sched->lSchedJobsWithSorting(sortJobsByGATC, kappas);
	cout << *sched;
	cout << "TWT = " << sched->getTWT();
	
	sched->reset();
	
	sched->clearJobs();

	return EXIT_SUCCESS;
}