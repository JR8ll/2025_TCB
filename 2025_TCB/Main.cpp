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

	pSched sched = TCB::prob->getSchedule();


	return EXIT_SUCCESS;
}