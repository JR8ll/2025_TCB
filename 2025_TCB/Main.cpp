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
// argv[3] int describing the solving method to be used: 1: Decomposition-MILP, 2: List Scheduling, 3: BRKGA, 4: BRKGA + final Decomposition-MILP with best sorting order (see constants in Functions.h), 5: ILS
// argv[4] time limit in seconds
// argv[5] console output on(=1)/off(=0)
// argv[6] filename of scheduling parameters
// argv[7] filename of ga parameters
// argv[8] filename of decompMILP parameters
// argv[9] filename of ILS parameters

int main(int argc, char* argv[]) {
	//Problem::genInstancesTCB26small_Testing();
	//Problem::genInstancesTCB26_Testing();
	//Problem::genInstancesTCB25_Jun25_exactMILPvsCP();
	//Problem::genInstancesEURO25_exact();
	//Problem::genInstancesEURO25_integer();
	Problem::genInstancesTCB26_BottleneckConfigs(100);
	Problem::genInstancesTCB26_BottleneckConfigs(200);

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
	ILS_params ilsParams = ILS_params();
	processCmd(argc, argv, iSolver, iTilimSeconds, bConsole, schedParams, gaParams, decompParams, ilsParams);

	if (!p.assertFeasibility()) {
		writeFailureReport();
		exit(EXIT_FAILURE);
	}

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
			// [JR-2026-Feb-23] 
			double omega = 0.9;
			if (decompParams.nDash >= p.getN()) {
				omega = 1.0;
			}

			Solver_MILP cplex = Solver_MILP(schedParams, decompParams, omega);
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
				sched->saveJsonFactory("BRKGA2MILP");
			}
		}
		break;
	case ALG_ILS:
		solverName = "ILS";																												// 5
		{
			Solver_ILS ils = Solver_ILS(schedParams, ilsParams);
			initializer<pJob> init = &Schedule::lSchedJobsWithSorting;
			ils.solveILS(*sched.get(), init, sortJobsRandomly, iTilimSeconds);
		}
		break;
	case ALG_ILS_PARALLELIZED:
		solverName = "ILS_PARALLEL";																									// 6
		{
			Solver_ILS ils = Solver_ILS(schedParams, ilsParams);
			initializer<pJob> init = &Schedule::lSchedJobsWithSorting;
			ils.solveILSparallelized(*sched.get(), init, sortJobsRandomly, iTilimSeconds);
		}
		break;
	case ALG_ILS_SEQUENCE: 
		solverName = "ILS_SEQ";																											// 7
		{
			Solver_Sequence_ILS ilsSeq = Solver_Sequence_ILS(sched.get(), &schedParams, &ilsParams, 1);									// [JR-2026-Feb-06] one core used
			ilsSeq.solveILSseq(*sched.get(), iTilimSeconds);
		}
		break;
	case ALG_ILS_SEQUENCE_PARALLELIZED:
		solverName = "ILS_SEQ_PARALLEL";																								// 8
		{
			int nCores = Solver_ILS::GetPCoreIndices(false).size();
			Solver_Sequence_ILS ilsSeq = Solver_Sequence_ILS(sched.get(), &schedParams, &ilsParams, nCores);
			ilsSeq.solveILSseqParallelized(*sched.get(), iTilimSeconds);
		}
		break;
	case ALG_ILS_HYBRID:
		solverName = "ILS_HYBRID";																										// 9
		{
			int nCores = Solver_ILS::GetPCoreIndices(false).size();
			Solver_Hybrid_ILS ilsHybrid = Solver_Hybrid_ILS(sched.get(), &schedParams, &ilsParams, nCores);
			ilsHybrid.solveILShybrid(*sched.get(), iTilimSeconds);
		}
		break;
	case ALG_ITMILPLSHIFT: 
		solverName = "MILP2LSHIFT";																										// 10
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
	writeSolutions(sched.get(), iSolver, solverName, objectiveName, iTilimSeconds, usedTime.count(), &schedParams, &gaParams, &decompParams, &ilsParams);	
	sched->saveJsonFactory(solverName);

	// CONSOLE OUTPUT
	cout << "Solved " << TCB::prob->filename << " using " << solverName << " in " << usedTime.count() << " seconds with TWT = " << sched->getTWT() << "." << endl;
	if (bConsole) {
		cout << *sched;
	}

	return EXIT_SUCCESS;
}