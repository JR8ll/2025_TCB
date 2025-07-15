#include <algorithm>
#include <iostream>
#include <queue>
#include <map>
#include <Windows.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "Schedule.h"
#include "Functions.h"
#include "Problem.h"
#include "Solver_MILP.h"	// DECOMPMILP_params

namespace pt = boost::property_tree;
using namespace std;

Schedule::Schedule() {
	workcenters = vector<pWc>();
	unscheduledJobs = vector<pJob>();
	scheduledJobs = vector<pJob>();
	problem = nullptr;
}

ostream& operator<<(ostream& os, const Schedule& sched) {
	os << "---------SCHEDULE ----------" << endl;
	for (size_t wc = 0; wc < sched.size(); ++wc) {
		os << sched[wc] << endl;
	}
	return os;
}

Workcenter& Schedule::operator[](size_t idx) { return *workcenters[idx]; }
Workcenter& Schedule::operator[](size_t idx) const { return *workcenters[idx]; }

Job& Schedule::getJob(size_t idx) { return *unscheduledJobs[idx]; };
Job& Schedule::getJob(size_t idx) const { return *unscheduledJobs[idx]; }

pJob Schedule::get_pJob(size_t idx) {
	if (idx >= unscheduledJobs.size()) throw out_of_range("Schedule::get_pJob() out of range");
	pJob returnJob = move(unscheduledJobs[idx]);
	unscheduledJobs.erase(unscheduledJobs.begin() + idx);
	return returnJob;
}

std::unique_ptr<Schedule> Schedule::clone() const {
	auto newSchedule = make_unique<Schedule>();
	for (const auto& wc : workcenters) {
		newSchedule->addWorkcenter(wc->clone(newSchedule.get()));
	}
	for (const auto& job : unscheduledJobs) {
		newSchedule->addJob(move(job->clone()));
	}
	for (const auto& job : scheduledJobs) {
		newSchedule->markAsScheduled(move(job->clone()));
	}
	newSchedule->_reconstruct(this);
	return newSchedule;
}

void Schedule::_reconstruct(const Schedule* orig) {
	for (size_t wc = 0; wc < (*orig).size(); ++wc) {
		for (size_t m = 0; m < (*orig)[wc].size(); ++m) {
			for (size_t b = 0; b < (*orig)[wc][m].size(); ++b) {
				(*this)[wc][m].addBatch(move((*orig)[wc][m][b].clone()), (*orig)[wc][m][b].getStart());
				for (size_t j = 0; j < (*orig)[wc][m][b].size(); ++j) {
					Operation* remoteOp = &(*orig)[wc][m][b][j];
					Operation* localOp = findInScheduledJobs(remoteOp);
					if (localOp == nullptr) {
						localOp = findInUnscheduledJobs(remoteOp);
					}
					if (localOp == nullptr) throw ExcSched("Schedule::_reconstruct() operation not found");
					(*this)[wc][m][b].addOp(localOp);
				}
			}
		}
	}
	for (size_t wc = 0; wc < size(); ++wc) {
		for (size_t m = 0; m < workcenters[wc]->size(); ++m) {
			for (size_t b = 0; b < (*workcenters[wc])[m].size(); ++b) {
				(*this)[wc][m][b].setStart((*this)[wc][m][b].getStart());
			}
		}
	}
}

size_t Schedule::size() const { return workcenters.size();  }
size_t Schedule::getN() const { return unscheduledJobs.size(); }

bool Schedule::contains(Operation* op) const {
	for (size_t wc = 0; wc < size(); ++wc) {
		for (size_t m = 0; m < (*workcenters[wc]).size(); ++m) {
			for (size_t b = 0; b < (*workcenters[wc])[m].size(); ++b) {
				for (size_t o = 0; o < (*workcenters[wc])[m][b].size(); ++o) {
					if (op->getId() == (*workcenters[wc])[m][b][o].getId() && op->getStg() == (*workcenters[wc])[m][b][o].getStg()) {
						return true;
					}
				}
			}
		}
	}
	return false;
}

int Schedule::getCapAtStageIdx(size_t stgIdx) const {
	if (stgIdx >= size()) throw out_of_range("Schedule::getCapAtStageIdx() out of range");
	return (*workcenters[stgIdx])[0].getCap();
}

const vector<int> Schedule::getBatchingStages() const {
	vector<int> batchingStages = vector<int>();
	for (size_t o = 0; o < size(); ++o) {
		if ((*workcenters[o])[0].getCap() > 1) {	// assumption: parallel identical machines
			batchingStages.push_back(o + 1);
		}
	}
	return batchingStages;
}
const std::vector<int> Schedule::getDiscreteStages() const {
	vector<int> discreteStages = vector<int>();
	for (size_t o = 0; o < size(); ++o) {
		if ((*workcenters[o])[0].getCap() <= 1) {	// assumption: parallel identical machines
			discreteStages.push_back(o + 1);
		}
	}
	return discreteStages;
}

const std::vector<pWc>& Schedule::getWorkcenters() const {
	return workcenters;
}
void Schedule::addWorkcenter(pWc wc) {
	workcenters.push_back(move(wc));
}
void Schedule::addJob(pJob job) {
	unscheduledJobs.push_back(move(job));
}

void Schedule::schedOp(Operation* op, double pWait) {
	int wcIdx = op->getWorkcenterId() - 1;
 	workcenters[wcIdx]->schedOp(op, pWait);
}

Problem* Schedule::getProblem() const {
	return problem;
}

void Schedule::setProblemRef(Problem* prob) {
	problem = prob;
}

void Schedule::reset() {
	for (size_t wc = 0; wc < size(); ++wc) {
		for (size_t m = 0; m < workcenters[wc]->size(); ++m) {
			(*workcenters[wc])[m].removeAllBatches();
		}
	}
	while (!scheduledJobs.empty()) {
		shiftJobFromVecToVec(scheduledJobs, unscheduledJobs, 0);
	}

	for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
		for (size_t o = 0; o < (*unscheduledJobs[j]).size(); ++o) {
			(*unscheduledJobs[j])[o].setWait(0);
		}
	}
}
void Schedule::clearJobs() {
	for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
		for (size_t o = 0; o < (*unscheduledJobs[j]).size(); ++o) {	
			(*unscheduledJobs[j])[o].setPred(nullptr);
			(*unscheduledJobs[j])[o].setSucc(nullptr);
		}
	}

	for (size_t j = 0; j < scheduledJobs.size(); ++j) {
		for (size_t o = 0; o < (*scheduledJobs[j]).size(); ++o) {
			(*scheduledJobs[j])[o].setPred(nullptr);
			(*scheduledJobs[j])[o].setSucc(nullptr);
		}
	}
	unscheduledJobs.clear();
	scheduledJobs.clear();
}

void Schedule::sortUnscheduled(prioRule<pJob> rule) {
	rule(unscheduledJobs);
}

void Schedule::sortUnscheduled(prioRuleKappa<pJob> rule, double kappa) {
	double t = getMinMSP(0);
	rule(unscheduledJobs, t, kappa);
}

void Schedule::sortUnscheduled(prioRuleKeySet<pJob> rule, std::vector<double>& chr) {
	rule(unscheduledJobs, chr);
}

void Schedule::sortScheduled(prioRule<pJob> rule) {
	updateWaitingTimes();
	rule(scheduledJobs);
}

void Schedule::updateWaitingTimes() {
	for (size_t wc = 0; wc < size(); ++wc) {
		workcenters[wc]->updateWaitingTimes();
	}
}

void Schedule::mimicWaitingTimes(const Schedule* wtSched) {
	map<pair<int, int>, double> mapWait = map<pair<int, int>, double>();	// <job id, stage id>, waiting time 
	for (size_t wc = 0; wc < wtSched->size(); ++wc) {
		for (size_t m = 0; m < (*wtSched)[wc].size(); ++m) {
			for (size_t b = 0; b < (*wtSched)[wc][m].size(); ++b) {
				for (size_t op = 0; op < (*wtSched)[wc][m][b].size(); ++op) {
					mapWait.insert(make_pair(make_pair((*wtSched)[wc][m][b][op].getId(), (*wtSched)[wc][m][b][op].getStg()), (*wtSched)[wc][m][b][op].getWait()));
				}
			}
		}
	}

	for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
		for (size_t o = 0; o < unscheduledJobs[j]->size(); ++o) {
			(*unscheduledJobs[j])[o].setWait(mapWait[make_pair((*unscheduledJobs[j])[o].getId(), (*unscheduledJobs[j])[o].getStg())]);
		}
	}

	for (size_t j = 0; j < scheduledJobs.size(); ++j) {
		for (size_t o = 0; o < scheduledJobs[j]->size(); ++o) {
			(*scheduledJobs[j])[o].setWait(mapWait[make_pair((*scheduledJobs[j])[o].getId(), (*scheduledJobs[j])[o].getStg())]);
		}
	}
}

void Schedule::markAsScheduled(size_t jobIdx) {
	if (jobIdx >= unscheduledJobs.size()) throw out_of_range("Schedule::markAsScheduled() out of range");
	shiftJobFromVecToVec(unscheduledJobs, scheduledJobs, jobIdx);
}
void Schedule::markAsScheduled(pJob scheduledJob) {
	scheduledJobs.push_back(move(scheduledJob));
}

int Schedule::getNumberOfScheduledJobs() const {
	return scheduledJobs.size();
}

const Job* Schedule::getScheduledJob(size_t idx) const
{
	if (idx >= scheduledJobs.size()) throw out_of_range("Schedule::getScheduledJob() out of range");
	return scheduledJobs[idx].get();
}

Operation* Schedule::findInScheduledJobs(Operation* remoteOp) const {
	for (size_t j = 0; j < scheduledJobs.size(); ++j) {
		if (scheduledJobs[j]->getId() == remoteOp->getId()) {
			for (size_t o = 0; o < scheduledJobs[j]->size(); ++o) {
				if ((*scheduledJobs[j])[o].getStg() == remoteOp->getStg()) {
					return &(*scheduledJobs[j])[o];	// local op
				}
			}
		}
	}
	return nullptr;
}
Operation* Schedule::findInUnscheduledJobs(Operation* remoteOp) const {
	for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
		if (unscheduledJobs[j]->getId() == remoteOp->getId()) {
			for (size_t o = 0; o < unscheduledJobs[j]->size(); ++o) {
				if ((*unscheduledJobs[j])[o].getStg() == remoteOp->getStg()) {
					return &(*unscheduledJobs[j])[o];	// local op
				}
			}
		}
	}
	return nullptr;
}

void Schedule::lSchedFirstJob(double pWait) {
	for (size_t op = 0; op < (*unscheduledJobs.begin())->size(); ++op) {
		schedOp(&(**unscheduledJobs.begin())[op], pWait);
		cout << *this;
	}
	shiftJobFromVecToVec(unscheduledJobs, scheduledJobs, 0);
}
void Schedule::lSchedJobs(double pWait) {
	while(!unscheduledJobs.empty()) {
		lSchedFirstJob(pWait);
	}
}

void Schedule::lSchedJobs(vector<double> pWaitVec) {
	double bestTWT = DBL_MAX; // numeric_limits<double>::max();
	double bestWait = pWaitVec[0];
	for (size_t w = 0; w < pWaitVec.size(); ++w) {
		unique_ptr<Schedule> copySched = this->clone();
		copySched->lSchedJobs(pWaitVec[w]);
		double tempTWT = copySched->getTWT();
		if (tempTWT < bestTWT) {
			bestTWT = tempTWT;
			bestWait = pWaitVec[w];
		}
	}
	TCB::logger.Log(Info, "List Scheduling with best pWait " + to_string(bestWait));
	lSchedJobs(bestWait);
}

void Schedule::lSchedJobsWithSorting(prioRule<pJob> rule, double pWait) {
	rule(unscheduledJobs);
	lSchedJobs(pWait);
}
void Schedule::lSchedJobsWithSorting(prioRule<pJob> rule, Sched_params& sched_params) {
	vector<double> pWaitVec = getDoubleGrid(sched_params.pWaitLow, sched_params.pWaitHigh, sched_params.pWaitStep);
	if (pWaitVec.size() == 1) {
		lSchedJobsWithSorting(rule, pWaitVec[0]);
	} else {
		double bestTWT = DBL_MAX; // numeric_limits<double>::max();
		double bestWait = pWaitVec[0];
		for (size_t w = 0; w < pWaitVec.size(); ++w) {
			unique_ptr<Schedule> copySched = this->clone();
			copySched->lSchedJobsWithSorting(rule, pWaitVec[w]);
			double tempTWT = copySched->getTWT();
			if (tempTWT < bestTWT) {
				bestTWT = tempTWT;
				bestWait = pWaitVec[w];
			}
		}
		lSchedJobsWithSorting(rule, bestWait);
	}
}

void Schedule::lSchedJobsWithSorting(prioRuleKappa<pJob> rule, double kappa, double pWait) {
	double t = 0.0;	// dynamic computation of priority index (increase t)
	while (!unscheduledJobs.empty()) {
		t = getMinMSP(0);
		rule(unscheduledJobs, t, kappa);
		lSchedFirstJob(pWait);
	}
}

double Schedule::lSchedJobsWithSorting(prioRuleKappa<pJob> rule, const vector<double>& kappaGrid, double pWait, objectiveFunction objectiveFunction) {
	if (kappaGrid.size() == 1) {
		lSchedJobsWithSorting(rule, kappaGrid[0], pWait);
		return kappaGrid[0];
	}
	
	double bestObjectiveValue = DBL_MAX; //  numeric_limits<double>::max();
	double bestKappa = 0.0;
	for (size_t kappa = 0; kappa < kappaGrid.size(); ++kappa) {
		lSchedJobsWithSorting(rule, kappaGrid[kappa], pWait);
		double tempObjectiveValue = objectiveFunction(this);
		if (tempObjectiveValue < bestObjectiveValue) {
			bestObjectiveValue = tempObjectiveValue;
			bestKappa = kappaGrid[kappa];
		}
		reset();
	}
	lSchedJobsWithSorting(rule, bestKappa, pWait);
	TCB::logger.Log(Info, "Found a schedule with best kappa value = " + to_string(bestKappa));
	return bestKappa;
}

double Schedule::lSchedJobsWithSorting(prioRuleKappa<pJob> rule, Sched_params& sched_params, objectiveFunction objectiveFunction) {
	vector<double> pWaitVec = getDoubleGrid(sched_params.pWaitLow, sched_params.pWaitHigh, sched_params.pWaitStep);
	vector<double> kappas = getDoubleGrid(sched_params.kappaLow, sched_params.kappaHigh, sched_params.kappaStep);

	double bestTWT = DBL_MAX; //  numeric_limits<double>::max();
	double bestWait = pWaitVec[0];
	for (size_t w = 0; w < pWaitVec.size(); ++w) {
		unique_ptr<Schedule> copySched = this->clone();
		copySched->lSchedJobsWithSorting(rule, kappas, pWaitVec[w]);
		double tempTWT = copySched->getTWT();
		if (tempTWT < bestTWT) {
			bestTWT = tempTWT;
			bestWait = pWaitVec[w];
		}
	}
	return lSchedJobsWithSorting(rule, kappas, bestWait);
}

void Schedule::lSchedJobsWithRandomKeySorting(prioRuleKeySet<pJob> rule, const std::vector<double>& keys, double pWait) {
	rule(unscheduledJobs, keys);
	lSchedJobs(pWait);
}
void Schedule::lSchedJobsWithRandomKeySorting(prioRuleKeySet<pJob> rule, const std::vector<double>& keys, Sched_params& sched_params) {
	vector<double> pWaitVec = getDoubleGrid(sched_params.pWaitLow, sched_params.pWaitHigh, sched_params.pWaitStep);
	double bestTWT = DBL_MAX; // numeric_limits<double>::max();
	double bestWait = pWaitVec[0];
	for (size_t w = 0; w < pWaitVec.size(); ++w) {
		unique_ptr<Schedule> copySched = this->clone();
		copySched->lSchedJobsWithRandomKeySorting(rule, keys, pWaitVec[w]);
		double tempTWT = copySched->getTWT();
		if (tempTWT < bestTWT) {
			bestTWT = tempTWT;
			bestWait = pWaitVec[w];
		}
	}
	lSchedJobsWithRandomKeySorting(rule, keys, bestWait);
}

void Schedule::lSchedGifflerThompson(prioRule<pJob> rule, double pWait) {
	// 1) consider all "next" operations of jobs to be processed
	vector<queue<Operation*>> unscheduledOps = vector<queue<Operation*>>(unscheduledJobs.size());
	for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
		for (size_t o = 0; o < (*unscheduledJobs[j]).size(); ++o) {
			unscheduledOps[j].push(&(*unscheduledJobs[j])[o]);
		}
	}

	// 2) get earliest completion time at any suitable workcenter => op*: operation with min C, wc*: corresponding workcenter 
	bool allOpsScheduled = false;
	double earliestC = DBL_MAX; // numeric_limits<double>::max();
	while (!allOpsScheduled) {
		size_t nextOp = 0;
		allOpsScheduled = true;
		for (size_t j = 0; j < unscheduledOps.size(); ++j) {
			if (!unscheduledOps[j].empty()) {
				allOpsScheduled = false;

				size_t bestMacIdx = 0;
				size_t bestBatIdx = 0;
				double bestStart = 0.0;
				//double tempC getEarliestC(unscheduledOps[j].front());

			}
		}
	} 
	
	// 3) further consider all operations which can be started before the completion of op* at wc* (Overlapping set)

	// Original Giffler & Thompson algorithm: branch (construct all possible schedules with any of the ops of 3) to be processed next and loop back to 1))
	// List scheduling adaptation: assign op from set 3) by priority index, loop back to 1)


}

void Schedule::localSearchLeftShifting(prioRule<pJob> rule, double pWait) {
	bool bImproved = true;
	while (bImproved) {
		bImproved = false;
		rule(scheduledJobs);
		for (size_t j = 0; j < scheduledJobs.size(); ++j) {
			for (size_t o = 0; o < (*scheduledJobs[j]).size(); ++o) {
				size_t wcIdx = (*scheduledJobs[j])[o].getWorkcenterId() - 1;
				size_t mIdx = 0;
				size_t bIdx = 0;
				size_t jIdx = 0;
				if(workcenters[wcIdx]->locateOp(&(*scheduledJobs[j])[o], mIdx, bIdx, jIdx)) {
					if (workcenters[wcIdx]->leftShift(mIdx, bIdx, jIdx, pWait)) {
						bImproved = true;
					}
				}
			}
			/*for (int o = (*scheduledJobs[j]).size() - 1; o >= 0; --o) {
				size_t wcIdx = (*scheduledJobs[j])[o].getWorkcenterId() - 1;
				size_t mIdx = 0;
				size_t bIdx = 0;
				size_t jIdx = 0;
				if (workcenters[wcIdx]->locateOp(&(*scheduledJobs[j])[o], mIdx, bIdx, jIdx)) {
					if (workcenters[wcIdx]->leftShift(mIdx, bIdx, jIdx, pWait)) {
						bImproved = true;
					}
				}
			}*/
		}
	}
}

bool Schedule::isValid() const {
	// ALL OPERATIONS OF ALL JOBS ARE ASSIGNED
	for (size_t j = 0; j < problem->getN(); ++j) {
		for (size_t o = 0; o < (*problem)[j].size(); ++o) {
			if (!contains(&(*problem)[j][o])) {
				TCB::logger.Log(Error, "missing operation " + to_string((*problem)[j][o].getId()) + "." + to_string((*problem)[j][o].getStg()));
				cout << *this;
				return false;
			}
		}
	}

	// NO OVERLAPPING PROCESSING ON ANY MACHINE
	for (size_t wc = 0; wc < size(); ++wc) {
		for (size_t m = 0; m < (*workcenters[wc]).size(); ++m) {
			if ((*workcenters[wc])[m].hasOverlaps()) {
				TCB::logger.Log(Error, "overlapping processing of batches at machine " + to_string(workcenters[wc]->getId()) + "." + to_string((*workcenters[wc])[m].getId()));
				cout << *this;
				return false;
			}
		}
	}

	// NO OPERATION IS STARTED BEFORE ITS ROUTE PREDECESSOR IS COMPLETED + TIME CONSTRAINTS ARE MET
	for (size_t wc = 0; wc < size(); ++wc) {
		for (size_t m = 0; m < (*workcenters[wc]).size(); ++m) {
			for (size_t b = 0; b < (*workcenters[wc])[m].size(); ++b) {
				for (size_t o = 0; o < (*workcenters[wc])[m][b].size(); ++o) {
					if (!(*workcenters[wc])[m][b][o].checkProcessingOrder()) {
						TCB::logger.Log(Error, "Processing order violated");
						cout << *this;
						return false;
					}
					if (!(*workcenters[wc])[m][b][o].checkTimeConstraints()) {
						TCB::logger.Log(Error, "Time constraint violated");
						cout << *this;
						return false;
					}
				}
			}
		}
	}
	return true;
}

double Schedule::getTWT() const {
	double twt = 0;
	for (size_t wc = 0; wc < size(); ++wc) {
		twt += workcenters[wc]->getTWT();
	}
	return twt;
}

double Schedule::getMinMSP(size_t stgIdx) const {
	if (stgIdx >= size()) throw out_of_range("Schedule::getMSP() out of range");
	return workcenters[stgIdx]->getMinMSP();
}

void Schedule::saveJson(std::string solver) {
	pt::ptree treeFile;

	treeFile.put("Problem", TCB::prob->getFilename());
	treeFile.put("Solver", solver);
	stringstream ssTWT;
	ssTWT << setprecision(3) << fixed << getTWT();
	treeFile.put("TWT", ssTWT.str());

	pt::ptree treeSchedule;
	for (size_t wc = 0; wc < workcenters.size(); ++wc) {
		Workcenter* WC = &(*workcenters[wc]);
		pt::ptree treeWorkcenter;
		treeWorkcenter.put("stage", WC->getId());
		for (size_t m = 0; m < WC->size(); ++m) {
			Machine* MAC = &(*WC)[m];
			pt::ptree treeMachine;
			treeMachine.put("id", MAC->getId());
			treeMachine.put("capacity", MAC->getCap());
			for (size_t b = 0; b < MAC->size(); ++b) {
				Batch* BAT = &(*MAC)[b];
				pt::ptree treeBatch;
				stringstream ssStart;
				stringstream ssCompletion;
				ssStart << setprecision(3) << fixed << BAT->getStart();
				ssCompletion << setprecision(3) << fixed << BAT->getC();
				treeBatch.put("start", ssStart.str());
				treeBatch.put("completion", ssCompletion.str());
				for (size_t j = 0; j < BAT->size(); ++j) {
					Operation* OP = &(*BAT)[j];
					pt::ptree treeOp;
					treeOp.put("id", OP->getId());
					treeOp.put("stage", OP->getStg());
					treeBatch.add_child("operations.Operation", treeOp);
				}
				treeMachine.add_child("batches.Batch", treeBatch);
			}
			treeWorkcenter.add_child("machines.Machine", treeMachine);
		}
		treeSchedule.add_child("workcenters.Workcenter", treeWorkcenter);
	}
	treeFile.add_child("Schedule", treeSchedule);

	string probFileName = extractFileName(TCB::prob->getFilename());
	string probFileNameWithoutExtension = probFileName.substr(0, probFileName.find(".dat"));
	string solverName = solver;
	replaceWindowsSpecialCharsWithUnderscore(solverName);

	string sSeed = to_string(TCB::seed);
	string schedFileName = probFileNameWithoutExtension + "_" + solverName + sSeed + ".json";

	bool success = CreateDirectory(L".\\results", NULL);
	string pathAndFilename = string(".\\results\\").append(schedFileName);
	pt::write_json(pathAndFilename, treeFile);
}
void Schedule::saveJsonFactory(std::string solver) {
	pt::ptree treeFile;

	treeFile.put("Problem", TCB::prob->getFilename());
	treeFile.put("Solver", solver);
	stringstream ssTWT;
	ssTWT << setprecision(3) << fixed << getTWT();
	treeFile.put("TWT", ssTWT.str());

	pt::ptree treeFactory;
	pt::ptree arrayWorkareas;
	pt::ptree treeWorkarea;
	treeWorkarea.put("name", "Flow Shop");
	pt::ptree arrayWorkcenters;

	for (size_t wc = 0; wc < workcenters.size(); ++wc) {
		Workcenter* WC = &(*workcenters[wc]);
		pt::ptree treeWorkcenter;
		pt::ptree arrayResources;
		treeWorkcenter.put("name", "Stage " + to_string(WC->getId()));

		for (size_t m = 0; m < WC->size(); ++m) {
			Machine* MAC = &(*WC)[m];
			pt::ptree treeMachine;
			pt::ptree arrayLoad;
			treeMachine.put("name", "M" + to_string(MAC->getId()));
			
			for (size_t b = 0; b < MAC->size(); ++b) {
				Batch* BAT = &(*MAC)[b];
				pt::ptree treeBatch;
				pt::ptree arrayContent;
				stringstream ssStart;
				stringstream ssCompletion;
				stringstream ssP;
				ssStart << setprecision(3) << fixed << BAT->getStart();
				ssCompletion << setprecision(3) << fixed << BAT->getC();
				ssP << setprecision(3) << fixed << BAT->getP();
				treeBatch.put("id", "B" + to_string((b + 1)));
				treeBatch.put("type", "Batch");
				treeBatch.put("start", ssStart.str());
				treeBatch.put("C", ssCompletion.str());
				treeBatch.put("f", BAT->getF());
				treeBatch.put("p", ssP.str());
				treeBatch.put("capacity", BAT->getCap());
				for (size_t j = 0; j < BAT->size(); ++j) {
					Operation* OP = &(*BAT)[j];
					stringstream ssD;
					stringstream ssR;
					stringstream ssW;	
					ssD << setprecision(3) << fixed << OP->getD();
					ssR << setprecision(3) << fixed << OP->getR();
					ssW << setprecision(3) << fixed << OP->getW();
					pt::ptree treeOp;
					treeOp.put("id", to_string(OP->getId()) + "." + to_string(OP->getStg()));
					treeOp.put("type", "Operation");
					treeOp.put("d", ssD.str());
					treeOp.put("r", ssR.str());
					treeOp.put("p", to_string(OP->getP()));
					treeOp.put("f", to_string(OP->getF()));
					treeOp.put("s", to_string(OP->getS()));
					treeOp.put("w", ssW.str());
					arrayContent.push_back(make_pair("", treeOp));
				}
				treeBatch.add_child("content", arrayContent);
				arrayLoad.push_back(make_pair("", treeBatch));
			}
			treeMachine.add_child("load", arrayLoad);
			arrayResources.push_back(std::make_pair("", treeMachine));
		}
		treeWorkcenter.add_child("resources", arrayResources);
		arrayWorkcenters.push_back(make_pair("", treeWorkcenter));
		
	}
	treeWorkarea.add_child("workcenters", arrayWorkcenters);
	arrayWorkareas.push_back(make_pair("", treeWorkarea));
	treeFactory.add_child("workareas", arrayWorkareas);
	treeFile.add_child("Factory", treeFactory);

	string probFileName = extractFileName(TCB::prob->getFilename());
	string probFileNameWithoutExtension = probFileName.substr(0, probFileName.find(".dat"));
	string solverName = solver;
	replaceWindowsSpecialCharsWithUnderscore(solverName);

	string sSeed = to_string(TCB::seed);
	string schedFileName = probFileNameWithoutExtension + "_" + solverName + sSeed + ".json";

	bool success = CreateDirectory(L".\\results", NULL);
	string pathAndFilename = string(".\\results\\").append(schedFileName);
	pt::write_json(pathAndFilename, treeFile);
}
