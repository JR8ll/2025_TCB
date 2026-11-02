#include <map>
#include <chrono>
#include <random>
#include <algorithm>

#include <random>
#include <limits>
#include <iostream>
#include <sstream>
#include <string>

#include "Solver_CP.h"
#include "Functions.h"
#include "Problem.h"
#include "Schedule.h"

ILOSTLBEGIN

Solver_CP::Solver_CP(Sched_params& sched_params) : sched(nullptr), schedParams(&sched_params) {
	env = IloEnv();
	modelCreated = false;
}
Solver_CP::~Solver_CP() {}

// Implemented by RK
bool Solver_CP::solveCP(int tilim) {
	bool solved = true;

	if (solveCount == solveLimit) {
		return false;
	}

	++solveCount;

	std::ofstream logfile("solver_output_cp.txt");

	cp.setOut(logfile);
	cp.setParameter(IloCP::Workers, 1);
	cp.setParameter(IloCP::LogVerbosity, IloCP::Normal);
	cp.setParameter(IloCP::TimeLimit, tilim);

	uniform_int_distribution<int> dist(0, INT_MAX);
	TCB::rng.seed(TCB::seed + solveCount);
	int randomSeed = dist(TCB::rng);
	cp.setParameter(IloCP::RandomSeed, randomSeed);

	// DEBUGGING (JR)
	//cp.exportModel("cpModel.cpo");
	if (!cp.solve()) {
		cout << "No solution found" << endl;
		return solved = false;
	}

	IloInt fails = cp.getInfo(IloCP::NumberOfFails);

	// CONSTRUCT SCHEDULE
	int nJobs = consideredJobs.size();
	int nStages = sched->size();
	int nProducts = sched->getProblem()->getF();

	for (int i = 0; i < nJobs; ++i) {
		int product = consideredJobs[i]->getF() - 1;
		int job = consideredJobs[i]->getId() - 1;

		for (int stage = 0; stage < nStages; ++stage) {
			double start = (double)cp.getStart(ops[product][job][stage]);

			int stageMaschines = (*sched)[stage].size();
			for (int machine = 0; machine < stageMaschines; ++machine) {
				if (cp.isPresent(opsOnMachine[product][job][stage][machine])) {
					bool inserted = false;
					auto& batches = (*sched)[stage][machine].getBatches();

					for (auto& batch : batches) { // Iteration über vorhandene Batches
						if (batch->getStart() == start && (batch->getF() - 1) == product) {
							batch->addOp(consideredJobs[i]->getOpPtr(stage));
							inserted = true;
							break;
						}
					}

					if (!inserted) { // Bildung eines neuen Batches
						int batchCapacity = sched->getCapAtStageIdx(stage);
						pBat newBatch = make_unique<Batch>(batchCapacity);
						newBatch->addOp(consideredJobs[i]->getOpPtr(stage));

						bool addBatch = (*sched)[stage][machine].addBatch(move(newBatch), start);
						if (!addBatch) {
							string errorMsg = "Problem in SolverCP::solve(): infeasible op assignment "
								+ to_string(job) + " at stage " + to_string(stage) + ": " + to_string(start);
							TCB::logger.Log(Error, errorMsg);
						}
					}
				}
			}
		}
	}

	return solved = true;
}

// Implemented by RK
double Solver_CP::solveStepCP(Schedule* schedule, BatchVector& oldBatches, int nDash, int tilim) {
	if (schedule->getProblem() == nullptr) throw ExcSched("solveDecompositionCP missing problem reference");

	// Consider the first nDash jobs
	if (consideredJobs.empty()) {
		for (size_t j = 0; j < nDash; ++j) {
			try {
				consideredJobs.push_back(move(schedule->get_pJob(0)));
			}
			catch (const out_of_range& ex) {
				break;
			}
		}

		setConstraints(consideredJobs, oldBatches);

		// JR WARMSTART
		setListSchedWarmStart(*schedule);
	}
}

// Implemented by RK
double Solver_CP::solveDecompCP(Schedule* schedule, Sched_params* schedParams, int nDash, int tilim, int method) {
	double obj = 0;
	sched = schedule;

	int nJobs = sched->getProblem()->getN();
	int nStages = sched->size();
	int nProducts = sched->getProblem()->getF();

	// Start- und Endzeiten
	TimeVector startTimes(nProducts, vector<vector<double>>(nJobs, vector<double>(nStages, 0)));
	TimeVector endTimes(nProducts, vector<vector<double>>(nJobs, vector<double>(nStages, 0)));

	// Batches und Verfügbarkeiten
	MachineVector machines(nProducts, vector<vector<int>>(nJobs, vector<int>(nStages, -1)));
	BatchVector oldBatches(nProducts, vector<vector<BatchMock>>(nStages));

	// Maschinenverfügbarkeiten
	vector<vector<double>> machineReleaseTimes(nStages);
	for (int stage = 0; stage < nStages; ++stage) {					// Iteration über Anzahl von Stages
		int stageMaschines = (*sched)[stage].size();
		machineReleaseTimes[stage].resize(stageMaschines, 0);		// Anzahl der Maschinen auf der Stage
	}

	double flowFactor = 1.6;
	double smoothingFactor = 0.7;
	schedule->initialWaitingTimeByFF(flowFactor);

	vector<double> kappas = getDoubleGrid(schedParams->kappaLow, schedParams->kappaHigh, schedParams->kappaStep);
	// OBTAIN INITIAL WAITING TIMES FROM A SCHEDULE WITH LIST SCHEDULING AND SORTING BY initRule (DEFAULT: EDD)
	unique_ptr<Schedule> initSched = schedule->clone();

	cout << "DEBUGGING - Solver_CP::solveDecompCP init scheduling commented. " << endl;
	//initSched->lSchedJobsWithSorting(sortJobsByD);
	//initSched->updateWaitingTimes();
	//schedule->mimicWaitingTimes(initSched.get());

	initModel();

	int iterationJR = 0;
	
	while (sched->getN() > 0 || !consideredJobs.empty()) {
		if (consideredJobs.empty()) {
			// DYNAMIC SORTING AND UPDATE WAITING TIMES

			if (method == 1) { // GATC + LTI
				unique_ptr<Schedule> tempSched = schedule->clone();
				unique_ptr<Schedule> bestSched = tempSched->lSchedJobsWithGATC(flowFactor, smoothingFactor, kappas);
				unordered_map<int, int> priorities = bestSched->getPriorities();
				schedule->mimicWaitingTimes(bestSched.get());
				schedule->sortUnscheduled(priorities);
			}

			if (method == 2) { // GATC
				unique_ptr<Schedule> tempSched = schedule->clone();
				double bestKappa = tempSched->lSchedJobsWithSorting(sortJobsByGATC, kappas);
				tempSched->updateWaitingTimes();
				schedule->mimicWaitingTimes(tempSched.get());
				schedule->sortUnscheduled(sortJobsByGATC, bestKappa);
			}

			if (method == 0) {	// DEBUGGING
				// no particular order
			}
		}

		obj += solveStepCP(sched, oldBatches, nDash, tilim);

		
		if (solveCP(tilim)) {
			updateModel(startTimes, endTimes, machines, oldBatches, obj);
			cp.exportModel("cpmodel.cpo");
			// DEBUGGING (JR)
			sched->saveJsonFactory("cpDebugging");
			// TODO DEBUGGING (JR) write model

			for (int i = 0; i < nDash; i++) {
				sched->markAsScheduled(move(consideredJobs[i]));
			}

			// DEBUGGING JR
			string debugScheduleName = "JR_DEBUG_STEP_" + to_string(iterationJR) + "___";
			sched->saveJsonFactory(debugScheduleName);
			iterationJR++;


			consideredJobs.clear();
			solveCount = 0;
			decompStep++;

		}
		else if (solveCount >= solveLimit) {
			staticUpdate(startTimes, endTimes, machines, oldBatches);
			solveCount = 0;
			decompStep++;
		}
	}

	return obj;
}

// Implemented by RK
void Solver_CP::initModel() {
	int nJobs = sched->getProblem()->getN();
	int nStages = sched->size();
	int nProducts = sched->getProblem()->getF();

	sched->sortUnscheduled(sortJobsById);	// JR

	// JR DEBUGGING
	ostringstream oss;

	// CUMUL & STATE FUNCTIONS
	for (int stage = 0; stage < nStages; ++stage) {
		int stageMaschines = (*sched)[stage].size();
		loadOfMachine.push_back(vector<IloCumulFunctionExpr>());	// load of machine on stage
		stateOfMachine.push_back(vector<IloStateFunction>());		// state of machine on stage
		for (int machine = 0; machine < stageMaschines; ++machine) {
			loadOfMachine[stage].push_back(IloCumulFunctionExpr(env));
			stateOfMachine[stage].push_back(IloStateFunction(env));	
			
			// JR DEBUGGING
			oss.str("");
			oss << "load_s" << stage + 1 << "_m" << machine + 1;
			loadOfMachine[stage][machine].setName(oss.str().c_str());
			oss.str("");
			oss << "family_s" << stage + 1 << "_m" << machine + 1;
			stateOfMachine[stage][machine].setName(oss.str().c_str());
			// JR DEBUGGING
		}
	}


	// INTERVALS
	for (int product = 0; product < nProducts; ++product) {
		ops.push_back(vector<IloIntervalVarArray>());
		opsOnMachine.push_back(vector<vector<IloIntervalVarArray>>());

		for (int job = 0; job < nJobs; ++job) {
			ops[product].push_back(IloIntervalVarArray(env));
			opsOnMachine[product].push_back(vector<IloIntervalVarArray>());

			for (int stage = 0; stage < nStages; ++stage) {
				int stageMaschines = (*sched)[stage].size();
				int processingTime = sched->getProblem()->getProduct(product)->getP(stage);

				if (job == 104) {
					int debugger = 999;
				}
				double releaseTime = sched->getJob(job).getR();

				ops[product][job].add(IloIntervalVar(env, processingTime));
				
				// JR DEBUGGING
				oss.str("");
				oss << "p" << product + 1 << "_j" << job + 1 << "_s" << stage + 1;
				ops[product][job][stage].setName(oss.str().c_str());	
				// JR DEBUGGING

				
						
				ops[product][job][0].setStartMin(releaseTime);		
				opsOnMachine[product][job].push_back(IloIntervalVarArray(env));
		


				for (int machine = 0; machine < stageMaschines; ++machine) {
					opsOnMachine[product][job][stage].add(IloIntervalVar(env));
					opsOnMachine[product][job][stage][machine].setOptional();
					loadOfMachine[stage][machine] += IloPulse(opsOnMachine[product][job][stage][machine], 1);

					// JR DEBUGGING
					oss.str("");
					oss << "opOnMac_p" << product + 1 << "_j" << job + 1 << "_s" << stage + 1 << "_m" << machine +1;
					opsOnMachine[product][job][stage][machine].setName(oss.str().c_str());

					oss.str("");
					oss << "load_s" << stage + 1 << "_m" << machine + 1;
					// JR DEBUGGING
				}
			}
		}
	}
}

// Implemented by RK
void Solver_CP::setConstraints(vector<pJob>& consideredJobs, BatchVector& oldBatches) {
	if (modelCreated) model.end();
	model = IloModel(env);
	int nStages = sched->size();
	int nProducts = sched->getProblem()->getF();

	cp = IloCP(model);
	cp.exportModel("EMPTY_CP_Modell.cpo");

	for (int i = 0; i < consideredJobs.size(); ++i) {
		int job = consideredJobs[i]->getId() - 1;
		int product = consideredJobs[i]->getF() - 1;
		for (int stage = 0; stage < nStages; ++stage) {
			int stageMaschines = (*sched)[stage].size();

			// Machines
			//model.add(IloAlternative(env, ops[0][0][0], opsOnMachine[0][0][0]));		// ???
			model.add(IloAlternative(env, ops[product][job][stage], opsOnMachine[product][job][stage]));

			// Liegezeitrestriktionen
			const vector<pair<int, double>>& tcMaxFwd = sched->getProblem()->getProduct(product)->getTcMaxFwd(stage);
			for (int successorStage = 1; successorStage < nStages; ++successorStage) {
				double tcDuration = DBL_MAX;
				for (auto tc = tcMaxFwd.begin(); tc != tcMaxFwd.end(); ++tc) {
					if (tc->first == successorStage) tcDuration = tc->second;
				} // TODO: find a better way to skip max tcs
				if (tcDuration < 999999) {
					model.add(IloStartBeforeStart(env, ops[product][job][successorStage], ops[product][job][stage], -tcDuration));
				}
			}

			// Kein Overlap beim Batching
			for (int machine = 0; machine < stageMaschines; ++machine) {
				model.add(IloAlwaysEqual(env, stateOfMachine[stage][machine], opsOnMachine[product][job][stage][machine], product, true, true));
			}
		}

		cp = IloCP(model);
		cp.exportModel("ONEJOB_CP_Modell.cpo");
	}

	// Vorrangsbeziehungen
	for (int i = 0; i < consideredJobs.size(); ++i) {
		int job = consideredJobs[i]->getId() - 1;
		int product = consideredJobs[i]->getF() - 1;
		for (int stage = 1; stage < nStages; ++stage) {
			model.add(IloEndBeforeStart(env, ops[product][job][stage - 1], ops[product][job][stage]));
		}
	}


	// Maximale Batching Kapazität
	for (int stage = 0; stage < nStages; ++stage) {
		int stageMaschines = (*sched)[stage].size();
		for (int machine = 0; machine < stageMaschines; ++machine) {
			int capacity = sched->getCapAtStageIdx(stage);
			model.add(loadOfMachine[stage][machine] <= capacity);
		}
	}


	// Alte Batches
	if (!oldBatches.empty()) {
		for (int stage = 0; stage < nStages; ++stage) {
			for (int product = 0; product < nProducts; ++product) {
				for (const auto& batch : oldBatches[product][stage]) {
					int capacityStage = sched->getCapAtStageIdx(stage);
					int processingTime = sched->getProblem()->getProduct(batch.product)->getP(batch.stage);
					int value = (batch.capacity == capacityStage) ? (product + 1) : batch.product;
					model.add(IloAlwaysEqual(env, stateOfMachine[stage][batch.machine], batch.start, batch.start + processingTime, value));
				}
			}
		}
	}

	// Zielfunktion
	IloExpr expr1(env);
	IloExpr expr2(env);
	for (int i = 0; i < consideredJobs.size(); ++i) {
		int job = consideredJobs[i]->getId() - 1;
		int product = consideredJobs[i]->getF() - 1;
		int dueDate = consideredJobs[i]->getD();
		double weight = consideredJobs[i]->getW();

		expr1 += IloMax(0, weight * (IloEndOf(ops[product][job][nStages - 1]) - dueDate));
		for (int stage = 0; stage < nStages; ++stage) expr2 += IloStartOf(ops[product][job][stage]);
	}


	static const int obj = 2; //TODO(RK): change to 2
	int omega = 9; //TODO: Change Omega
	if (obj == 0) model.add(IloMinimize(env, expr1));
	if (obj == 1) model.add(IloMinimize(env, IloStaticLex(env, expr1, expr2)));
	if (obj == 2) model.add(IloMinimize(env, omega * expr1 + expr2));

	expr1.end();
	expr2.end();

	if (modelCreated) cp.end();

	cp = IloCP(model);
	cp.exportModel("CP_Modell.cpo");

	modelCreated = true;

}

// Implemented by RK
void Solver_CP::updateModel(
	TimeVector& starts,							// Startzeiten
	TimeVector& ends,							// Endzeiten
	MachineVector& machines,					// Maschinenverfügbarkeit
	BatchVector& oldBatches,
	double& obj									// Zielfunktion
) {
	int nStages = sched->size();
	for (int i = 0; i < consideredJobs.size(); ++i) {
		int product = consideredJobs[i]->getF() - 1;
		int job = consideredJobs[i]->getId() - 1;

		for (int stage = 0; stage < nStages; ++stage) {
			starts[product][job][stage] = cp.getStart(ops[product][job][stage]);		// Startzeit
			ends[product][job][stage] = cp.getEnd(ops[product][job][stage]);			// Endzeit

			int machine = 0;
			while (cp.isAbsent(opsOnMachine[product][job][stage][machine])) {			// Prüfung ob unbelegtes Intervall
				++machine;
			}
			machines[product][job][stage] = machine;									// Welche Maschine ist frei

			loadOfMachine[stage][machine] += IloPulse(env, starts[product][job][stage], ends[product][job][stage], 1);

			int batchIndex = 0;
			auto& batches = oldBatches[product][stage];
			auto matches = [&](const BatchMock& b) {
				return  b.start == starts[product][job][stage] && b.machine == machines[product][job][stage];
				};

			while (batchIndex < batches.size() && !matches(batches[batchIndex])) {
				++batchIndex;
			}
			if (batchIndex == oldBatches[product][stage].size()) {
				BatchMock batch = BatchMock(product, stage, machine, starts[product][job][stage], 1);
				oldBatches[product][stage].push_back(batch);
			}
			else {
				++oldBatches[product][stage][batchIndex].capacity;
			}
		}
		int dueDate = consideredJobs[i]->getD();
		double weight = consideredJobs[i]->getW();
		if (ends[product][job][nStages - 1] > dueDate) {
			obj += weight * (ends[product][job][nStages - 1] - dueDate);
		}
	}
}

// Implemented by RK
void Solver_CP::staticUpdate(
	TimeVector& starts,							// Startzeiten
	TimeVector& ends,							// Endzeiten
	MachineVector& machines,					// Maschinenverfügbarkeit
	BatchVector& oldBatches
) {
	while (!consideredJobs.empty()) {
		Job* j = consideredJobs.front().get();

		sched->addJobFront(move(consideredJobs.front()));	//
		sched->lSchedFirstJob();							//

		for (const auto& wc : sched->getWorkcenters()) {
			int stage = wc->getId() - 1;
			int product = j->getF() - 1;
			int job = j->getId() - 1;
			auto* batch = j->getOpPtr(stage)->getBatch();
			int machine = batch->getMachine()->getIdx();

			double start = batch->getStart();
			double end = batch->getC();

			// Ich brauche Maschine, Start- und Endzeit
			starts[product][job][stage] = start;		// Startzeit
			ends[product][job][stage] = end;			// Endzeit
			machines[product][job][stage] = machine;	// Maschine

			loadOfMachine[stage][machine] += IloPulse(env, starts[product][job][stage], ends[product][job][stage], 1);

			int batchIndex = 0;
			auto& batches = oldBatches[product][stage];
			auto matches = [&](const BatchMock& b) {
				return  b.start == starts[product][job][stage] && b.machine == machines[product][job][stage];
				};

			while (batchIndex < batches.size() && !matches(batches[batchIndex])) {
				++batchIndex;
			}
			if (batchIndex == oldBatches[product][stage].size()) {
				BatchMock batch = BatchMock(product, stage, machine, starts[product][job][stage], 1);
				oldBatches[product][stage].push_back(batch);
			}
			else {
				++oldBatches[product][stage][batchIndex].capacity;
			}
		}

		consideredJobs.erase(consideredJobs.begin());

	}
}

void Solver_CP::setStartingPoint() {
	IloSolution startSol(env);

	int nJobs = sched->getNumberOfScheduledJobs();
	int nStages = sched->size();

	for (int i = 0; i < nJobs; ++i) {
		const Job* job = sched->getScheduledJob(i);
		int jobId = job->getId() - 1;
		int productId = job->getF() - 1;

		for (int stage = 0; stage < nStages; ++stage) {
			Operation* o = job->getOpPtr(stage);
			int stageId = o->getWorkcenterId() - 1;
			int machineId = o->getBatch()->getMachine()->getIdx();

			int start = o->getStart();
			int end = start + o->getP();

			IloIntervalVar op = opsOnMachine[productId][jobId][stageId][machineId];
			startSol.setPresent(op);
			startSol.setStart(op, start);
			startSol.setEnd(op, end);

			IloIntervalVar masterOp = ops[productId][jobId][stageId];
			startSol.setPresent(masterOp);
			startSol.setStart(masterOp, start);
			startSol.setEnd(masterOp, end);
		}

	}

	cp.setStartingPoint(startSol);
}

void Solver_CP::setListSchedWarmStart(Schedule& sched) {
	IloSolution warmStart(env);

	unique_ptr<Schedule> tempSched = sched.clone();
	tempSched->clearUnscheduledJobs();
	for (unsigned i = 0; i < consideredJobs.size(); ++i) {
		tempSched->addJob(move(consideredJobs[i]->clone()));		// CLONE??!
	}
		
	// list schedule considered jobs
	tempSched->lSchedJobsWithSorting(sortJobsDebugging);	// TODO sorting
	// retrieve solution and set as starting point for cpo

	tempSched->saveJsonFactory("WARMSTART");

	for (unsigned i = 0; i < consideredJobs.size(); ++i) {
		for (unsigned stageIdx = 0; stageIdx < tempSched->size(); ++stageIdx) {
			Operation* tempOp = tempSched->findInScheduledJobs(consideredJobs[i]->getOpPtr(stageIdx));
			double opStart = tempOp->getStart();
			double opCompletion = tempOp->getC();
			int opFamilyIdx = tempOp->getF() - 1;
			int opJobIdx = tempOp->getId() - 1;
			
			int machineIdx = tempOp->getBatch()->getMachine()->getIdx();
	

			IloIntervalVar op = opsOnMachine[opFamilyIdx][opJobIdx][stageIdx][machineIdx];
			warmStart.setPresent(op);
			warmStart.setStart(op, opStart);
			warmStart.setEnd(op, opCompletion);

			// optional interval: set !present = absent	JR DEBUGGING
			for (unsigned m = 0; m < sched.getWorkcenters()[stageIdx]->size(); ++m) {
				if (m != machineIdx) {
					opsOnMachine[opFamilyIdx][opJobIdx][stageIdx][m].setAbsent();
				}
			}

			IloIntervalVar masterOp = ops[opFamilyIdx][opJobIdx][stageIdx];
			warmStart.setPresent(masterOp);
			warmStart.setStart(masterOp, opStart);
			warmStart.setEnd(masterOp, opCompletion);
		}
	}
	cp.setStartingPoint(warmStart);
	cp.exportModel("AfterWarmStart.cpo");
}

void Solver_CP::setFixPartialSolution(bool fixMachineAssignments, bool fixMachineSequences) {
	int nJobs = sched->getNumberOfScheduledJobs();
	int nStages = sched->size();

	for (int i = 0; i < nJobs; ++i) {
		const Job* job = sched->getScheduledJob(i);
		int jobId = job->getId() - 1;
		int productId = job->getF() - 1;

		for (int stage = 0; stage < nStages; ++stage) {
			Operation* op = job->getOpPtr(stage);
			int stageId = op->getWorkcenterId() - 1;
			int machineId = op->getBatch()->getMachine()->getIdx();

			IloIntervalVar currOp = opsOnMachine[productId][jobId][stageId][machineId];

			if (fixMachineAssignments) model.add(IloPresenceOf(env, currOp) == 1);
			if (fixMachineSequences) {
				// Vorgänger-Operation auf derselben Maschine / im vorherigen Batch bestimmen
				Operation* pred = op->getPred();
				int prevProductId = pred->getF() - 1;
				int prevJobId = pred->getId() - 1;
				int prevStageId = pred->getStg() - 1;
				int prevMachineId = pred->getBatch()->getMachine()->getIdx();

				IloIntervalVar prevOp = opsOnMachine[prevProductId][prevJobId][prevStageId][prevMachineId];
				model.add(IloEndBeforeStart(env, prevOp, currOp));

				// Nachfolger-Operation auf derselben Maschine / im nächsten Batch bestimmen
				Operation* succ = op->getSucc(); // Nachfolger
				int nextProductId = succ->getF() - 1;
				int nextJobId = succ->getId() - 1;
				int nextStageId = succ->getStg() - 1;
				int nextMachineId = succ->getBatch()->getMachine()->getIdx();

				IloIntervalVar nextOp = opsOnMachine[nextProductId][nextJobId][nextStageId][nextMachineId];
				model.add(IloEndBeforeStart(env, currOp, nextOp));
			}
		}

	}
}

#pragma region LNS
void Solver_CP::solveLNS(Schedule* schedule) {
	sched = schedule;

	initModel();

	// Consider the first nDash jobs
	if (consideredJobs.empty()) {
		for (size_t j = 0; j < 50; ++j) {
			try {
				consideredJobs.push_back(move(schedule->get_pJob(0)));
			}
			catch (const out_of_range& ex) {
				break;
			}
		}
	}

	// DEBUGGER



	setConstraints(consideredJobs);

	sched->lSchedJobs(0.0);
	setStartingPoint();					// Warmstart
	setFixPartialSolution(true, false);	// Fixing

	std::ofstream logfile("solver_output_cp.txt");

	cp.setOut(logfile);
	cp.setParameter(IloCP::Workers, 1);
	cp.setParameter(IloCP::LogVerbosity, IloCP::Normal);
	cp.setParameter(IloCP::TimeLimit, 180);

	if (!cp.solve()) {
		cout << "No solution found" << endl;
	}
	else {
		cout << "Solution found " << cp.getObjValue();
	}

	cout << "";

}
#pragma endregion