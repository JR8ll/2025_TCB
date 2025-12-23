#include <chrono>
#include <iostream>

#include "Functions.h"
#include "Problem.h"
#include "Solver_ILS.h"
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

	//Problem::genInstancesTCB25_Jun25_exactMILPvsCP();
	//Problem::genInstancesEURO25_exact();
	//Problem::genInstancesEURO25();

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


	// DEBUGGING 2025-Apr-07
	// sort jobs 3 5 6 9 7 1 10 8 4 2
	//sched->sortUnscheduled(sortJobsDebugging);	// DEBUGGING change sortUnscheduled
	// solve with list-sched
	//sched->lSchedJobsInflated(0.0, 0.25, true, false);
	//sched->lSchedJobs(0.0);
	//cout << *sched;
	//sched->localSearchLeftShifting();	// Operation based leftShifting (too narrowly constrained to work properly)
	//cout << *sched;


	//DEBUGGING 


	sched->lSchedJobsWithSorting(sortJobsDebugging, 0);
	sched->saveJsonFactory("BEFORE_PERT");
	sched->perturbRandomJobSwap();
	sched->saveJsonFactory("AFTER_RND_SWAP");
	sched->perturbRandomRightShifting();
	sched->saveJsonFactory("AFTER_RND_RIGHTSHIFT");

	bool test = false;
	
	/*sched->saveJsonFactory("JOBWISE");

	sched->localSearchBatchConsolidation(false);
	sched->saveJsonFactory("AFTER_BATCHCONSOLIDATION");
	
	sched->localSearchConsolidateBatch(0, 0, 2, 0, 3, 0);

	sched->reset();
	sched->lSchedJobsStageWiseBackwardWithSorting(sortJobsDebugging, 0);
	sched->localSearchJobLeftShifting();
	sched->saveJsonFactory("AFTER_LEFTSHIFT");

	sched->reset();
	sched->lSchedJobsStageWiseWithSorting(sortJobsDebugging, 0);
	sched->saveJsonFactory("AFTER_STAGEWISE");


	sched->debugSetR(2, 40);
	sched->debugSetR(8, 0);
	sched->debugAddMachine(3);
	sched->saveJsonFactory("BEFORE");
	vector<vector<pair<double, double>>> testOption = vector<vector<pair<double, double>>>();

	Operation* testOp = sched->getScheduledJob(3)->getOpPtr(1);
	vector<pair<double, double>> rShiftOptions = sched->getRightShiftOptions(testOp, 25);
	pair<double, double> testR0 = sched->locSearchEvaluateJobRightShift(3, 1, 25, testOption);

	sched->localSearchJobLeftShifting(&sortJobsByD, false);
	sched->saveJsonFactory("AFTERLSHIFT");
	sched->localSearchJobSwapping();
	*/

	// SOLVE
	switch (iSolver) {
	case ALG_ITERATEDMILP:
		solverName = "DecompMILP";
		{
			Solver_MILP cplex = Solver_MILP(schedParams, decompParams);
			cplex.solveDecompJobBasedDynamicSortingGridMILP(sched.get(), decompParams.nDash, decompParams.cplexTilim, sortJobsByD, sortJobsByGATC);
		}
		break;
	case ALG_LISTSCHEDATC: 
		solverName = "ListSchedGATC";
		{
			sched->lSchedJobsWithSorting(sortJobsByGATC, schedParams);	
			// MISC REPORTING
			double twtBefore = sched->getTWT();
			sched->localSearchOpLeftShifting();	// try parameters sortJobsByC, sortJobsByStart, ...pWait
			double twtAfter = sched->getTWT();
			schedParams.leftShiftImprovement = (twtBefore - twtAfter) / twtBefore;
		}
		break;
	case ALG_BRKGALISTSCH:
		solverName = "BRKGA"; 
		{
			Solver_GA brkga = Solver_GA(schedParams, gaParams);
			brkga.solveBRKGA_List_jobBased(*sched.get(), iTilimSeconds);
		}
		break;
	case ALG_BRKGALS2MILP:
		solverName = "BRKGA_MILP";
		{
			pSched gaSched = sched->clone();
			Solver_GA brkga = Solver_GA(schedParams, gaParams);
			brkga.solveBRKGA_List_jobBased(*gaSched.get(), iTilimSeconds);
			if (brkga.hasCompleted()) {
				gaSched->saveJson("BRKGA");
				vector<double> bestChromosome = brkga.getBestChromosome();
				sched->sortUnscheduled(sortJobsByRK, bestChromosome);
				Solver_MILP cplex = Solver_MILP(schedParams, decompParams);
				cplex.solveDecompJobBasedMILP(sched.get(), decompParams.nDash, decompParams.cplexTilim);
				sched->saveJson("BRKGA2MILP");
			}
		}
		break;
	case ALG_ILS:
		solverName = "ILS";
		{
			Solver_ILS ils = Solver_ILS(schedParams);
			initializer<pJob> init = &Schedule::lSchedJobsWithSorting;
			ils.solveILS(*sched.get(), init, sortJobsRandomly, 1, iTilimSeconds);	// 4th parameter = 1 means no multistart
		}
		break;
	case ALG_ITMILPLSHIFT: 
		solverName = "MILP2LSHIFT";
		{
			Solver_MILP cplex = Solver_MILP(schedParams, decompParams);
			cplex.solveDecompJobBasedDynamicSortingGridMILP(sched.get(), decompParams.nDash, decompParams.cplexTilim, sortJobsByD, sortJobsByGATC);
			double twtBefore = sched->getTWT();
			sched->localSearchOpLeftShifting();
			double twtAfter = sched->getTWT();
			schedParams.leftShiftImprovement = (twtBefore - twtAfter) / twtBefore;
		}
	default:
		TCB::logger.Log(Error, "Program was executed with no valid algorithm key");
	} 

	// STOP TIME MEASUREMENT
	stop = chrono::high_resolution_clock::now();
	usedTime = chrono::duration_cast<chrono::seconds>(stop - start);

	// RESULT SUMMARY (FILE OUTPUT)
	writeSolutions(sched.get(), iSolver, solverName, objectiveName, iTilimSeconds, usedTime.count(), &schedParams, &gaParams, &decompParams);	
	sched->saveJsonFactory(solverName);

	// CONSOLE OUTPUT
	if (bConsole) {
		cout << "Solved using " << solverName << " in " << " seconds with TWT = " << sched->getTWT() << "." << endl;
		cout << *sched;
	}

	return EXIT_SUCCESS;
}