#include <algorithm>
#include <iostream>
#include <queue>
#include <map>
#include <numeric>
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
void Schedule::schedOp(Operation* op, double pWait, double inflation, bool batchinStageInflationOnly, bool opsWithoutTcInflationOnly) {
	int wcIdx = op->getWorkcenterId() - 1;
	workcenters[wcIdx]->schedOp(op, pWait, inflation, batchinStageInflationOnly, opsWithoutTcInflationOnly);
}
void Schedule::schedOpDelayed(Operation* op, double startingAt) {
	int wcIdx = op->getWorkcenterId() - 1;
	workcenters[wcIdx]->schedOpDelayed(op, startingAt);
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
		//cout << *this;
	}
	shiftJobFromVecToVec(unscheduledJobs, scheduledJobs, 0);
}

void Schedule::lSchedFirstJobInflated(double pWait, double inflation, bool batchinStageInflationOnly, bool opsWithoutTcInflationOnly) {
	for (size_t op = 0; op < (*unscheduledJobs.begin())->size(); ++op) {
		schedOp(&(**unscheduledJobs.begin())[op], pWait, inflation, batchinStageInflationOnly, opsWithoutTcInflationOnly);
		cout << *this;
	}
	shiftJobFromVecToVec(unscheduledJobs, scheduledJobs, 0);
}
void Schedule::lSchedJobsStageWise(double pWait) {
	for (size_t i = 0; i < size(); ++i) {	
		
		for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
			Operation* op = &(*unscheduledJobs[j])[i];

			//DEBUG
			if (op->getId() == 5 && op->getStg() == 3) {
				int stop = 666;
			}


			int maxLookAhead = workcenters[i]->getCap() - op->getS();	// looking ahead, ops may be batched if their common start time thus is earlier than their latest start if not batched
			int lookAheadCapRqrmt = op->getS();
			
			if (!op->isScheduled()) {
				vector<Operation*> lookingAhead = vector<Operation*>();
				double earliest = op->getEarliestStart();
				double earliestC = earliest + op->getP();
				double latest = min(earliest + op->getP(), op->getLatestStartConsideringTc());
				double commonStart = earliest;

				if(maxLookAhead > 0) {
					// this operation should not wait if a full batch is ready to start at earliestC
					bool bFullBatchWaiting = false;
					/*int remainingCap = workcenters[i]->getCap();
					for (size_t k = j + 1; k < unscheduledJobs.size(); ++k) {
						Operation* lookAheadOp = &(*unscheduledJobs[k])[i];
						if (lookAheadOp->getF() == op->getF()) {
							if (!lookAheadOp->isScheduled()) {
								double tempEarliest = lookAheadOp->getEarliestStart();
								if (tempEarliest > earliest && tempEarliest <= earliestC) {
									remainingCap -= lookAheadOp->getS();
								}
							}
						}
						if (remainingCap <= 0) {
							bFullBatchWaiting = true;
							break;
						}
					}*/

					// look for operations worthy to wait for
					if (!bFullBatchWaiting) {
						for (size_t k = j + 1; k < unscheduledJobs.size(); ++k) {
							Operation* lookAheadOp = &(*unscheduledJobs[k])[i];
							if (lookAheadOp->getF() == op->getF()) {
								if (!lookAheadOp->isScheduled()) {
									double tempEarliest = lookAheadOp->getEarliestStart();
									if (tempEarliest <= latest) {
										if (workcenters[i]->getCap() >= lookAheadCapRqrmt + lookAheadOp->getS()) {
											// CANDIDATE FOUND
											lookAheadCapRqrmt += lookAheadOp->getS();
											// TODO: TRY MORE SOPHISTICATED DECISIONS (e.g. based on weights or anticipated wT)
											commonStart = max(earliest, tempEarliest);
											lookingAhead.push_back(lookAheadOp);
										}
									}
								}
							}
							if (lookingAhead.size() >= maxLookAhead || lookAheadCapRqrmt >= workcenters[i]->getCap()) break;
						}
					}
				}	
				schedOpDelayed(op, commonStart);
				for (size_t i = 0; i < lookingAhead.size(); ++i) {
					schedOpDelayed(lookingAhead[i], commonStart);
				}
				saveJsonFactory("DEBUGGING_STAGEWISE");
			}
		}
	}
	while (!unscheduledJobs.empty()) {
		shiftJobFromVecToVec(unscheduledJobs, scheduledJobs, 0);
	}
	if (!this->isValid()) {
		throw ExcSched("ERROR: invalid schedule after Schedule::lSchedStages.");
	}
}
void Schedule::lSchedJobsStageWiseWithSorting(prioRule<pJob> rule, double pWait) {
	rule(unscheduledJobs);
	lSchedJobsStageWise(pWait);
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
void Schedule::lSchedJobsInflated(double pWait, double inflation, bool batchinStageInflationOnly, bool opsWithoutTcInflationOnly) {
	while (!unscheduledJobs.empty()) {
		lSchedFirstJobInflated(pWait, inflation, batchinStageInflationOnly, opsWithoutTcInflationOnly);
	}
	cout << *this << endl;
	leftShiftBatches();	// deflate
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

void Schedule::leftShiftBatches() {
	bool bImproved = true;
	while (bImproved) {
		bImproved = false;
		for (size_t wc = 0; wc < size(); ++wc) {
			for (size_t i = 0; i < workcenters[wc]->size(); ++i) {
				for (size_t b = 0; b < (*workcenters[wc])[i].size(); ++b) {
					if (workcenters[wc]->leftShift(i, b)) {
						cout << *this << endl;
						bImproved = true;
					}
				}
			}
		}
		
	}
}

void Schedule::localSearchBatchLeftShifting() {
	cout << "Schedule::localSearchBatchLeftShifting() not yet implemented. " << endl;

}

void Schedule::localSearchOpLeftShifting(prioRule<pJob> rule, double pWait) {
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
void Schedule::localSearchJobSwapping(prioRule<pJob> rule, bool bestFit) {
	int debug = 0;
	this->saveJsonFactory("BEFORE");
	bool bImproved = true;
	int nJobs = scheduledJobs.size();
	while (bImproved) {

		debug++;

		bImproved = false;
		rule(scheduledJobs);
		bool swapFeasible = true;
		double bestImprovement = 0;
		size_t best1 = 0;
		size_t best2 = 0;
		for (size_t j = 0; j < nJobs; ++j) {
			for (int k = nJobs - 1; k >= 0; --k) {
				if (j != k) {
					double tempImprovement = locSearchEvaluateJobSwap(j, k, swapFeasible);
					if (swapFeasible && tempImprovement > bestImprovement) {
						if (!bestFit) { // FIRST FIT
							locSearchSwapJobs(j, k);

							// DEBUGGING
							TCB::prob->filename = TCB::prob->filename + "_STEP";
							this->saveJsonFactory("DEBUG");
							bImproved = true;
							break;
						} else {			// BEST FIT
							bestImprovement = tempImprovement;
							best1 = j;
							best2 = k;
						}
					}
				}
			}
			if (bImproved) {
				break; // continue with while loop
			}
		}
		if (bestFit && bestImprovement > 0) {
			locSearchSwapJobs(best1, best2);
			TCB::prob->filename = TCB::prob->filename + "STEP_";
			this->saveJsonFactory("DEBUG");
			bImproved = true;
		}
	}
}
void Schedule::localSearchJobLeftShifting(prioRule<pJob> rule, bool bestFit){
	bool bImproved = true;
	int nJobs = scheduledJobs.size();
	while (bImproved) {
		bImproved = false;
		rule(scheduledJobs);
		pair<double, double> bestImprovement = make_pair(0, 0);
		double bestTwtImprovement = 0.0;
		size_t best = 0;
		vector<vector<pair<double, double>>> bestPossibleLeftShifts = vector<vector<pair<double, double>>>(size());
		for (size_t j = 0; j < nJobs; ++j) {
			if (j == 7) {
				saveJsonFactory("preNineShift");
				int stop = 666;
			}
			vector<vector<pair<double, double>>> tempPossibleLeftShifts = vector<vector<pair<double, double>>>(size());
			pair<double, double> tempImprovement = locSearchEvaluateJobLeftShift(j, tempPossibleLeftShifts);
			double tempTwtImprovement = tempImprovement.first * scheduledJobs[j]->getW();
			if (tempTwtImprovement > bestTwtImprovement || (tempTwtImprovement == bestTwtImprovement && tempImprovement.second > bestImprovement.second)) {
				if (!bestFit) {
					locSearchLeftShiftJob(j, tempPossibleLeftShifts);
					bImproved = true;
					break;
				}
				else {
					bestImprovement = tempImprovement;
					bestTwtImprovement = tempTwtImprovement;
					bestPossibleLeftShifts = tempPossibleLeftShifts;
					best = j;
				}
			}
		}
		if (bestFit && (bestImprovement.first > 0 || bestImprovement.second > 0)) {
			locSearchLeftShiftJob(best, bestPossibleLeftShifts);
			bImproved = true;
		}
	}
}

double Schedule::locSearchEvaluateJobSwap(size_t idxFirst, size_t idxSecond, bool& feasible) {
	Job* job1 = scheduledJobs[idxFirst].get();
	Job* job2 = scheduledJobs[idxSecond].get();

	if (job1->getF() != job2->getF()) {
		feasible = false;
		return -1.0;
	}

	if (job1->getStart() < job2->getR()
		|| job2->getStart() < job1->getR()) {
		feasible = false;
		return -1.0;
	}

	feasible = true;
	
	double grossBetter = (max(job1->getC() - job1->getD(), 0) * job1->getW())
		+ (max(job2->getC() - job2->getD(), 0) * job2->getW());

	double grossWorse = (max(job1->getC() - job2->getD(), 0) * job2->getC())
		+ (max(job2->getC() - job1->getD(), 0) * job1->getW());

	
	return grossBetter - grossWorse;
}

pair<double, double> Schedule::locSearchEvaluateJobLeftShift(size_t idxFirst, vector<vector<pair<double, double>>>& possibleLeftShift) {
	Job* job = scheduledJobs[idxFirst].get();
	possibleLeftShift = vector<vector<pair<double, double>>>(size());		// [stage][shift option] true => continuous value (not depending on batching)

	// 1) at each stage check the potential for left shifting disregarding maximal time lags (time constraints) and previous stages
	double earliest = job->getR();
	for (size_t o = 0; o < size(); ++o) {
		possibleLeftShift[o] = this->getLeftShiftOptions(job->getOpPtr(o));
	}

	// 2) remove mutually exclusive options 	
	vector<double> currentLeeway = getLeeway(job);	// as of now, how much time between end and start of succeeding operations? => if job[i] is left shifted by t, then job[i-1] must be left shifted by at least max((t - leeway[i]), 0)
	vector<vector<pair<size_t, double>>> currentTcSlack = getTcSlack(job);
	bool possibleOverlaps = false;
	bool possibleTcViolations = false;
	do {
		possibleOverlaps = constrainLeftShiftOptionsFromOverlaps(possibleLeftShift, currentLeeway);
		possibleTcViolations = constrainLeftShiftOptionsFromTimeConstraints(possibleLeftShift, currentTcSlack);
	} while (possibleOverlaps || possibleTcViolations);

	// 3) evaluate best options
	pair<double, double> evaluation = make_pair(possibleLeftShift[possibleLeftShift.size()-1][0].first, 0);	
	for (size_t i = 0; i < (possibleLeftShift.size() - 1); ++i) {
		evaluation.second += possibleLeftShift[i][0].first;
	}

	return evaluation; // first = actual left shift of job (=last operation), second = secondary left shifts of intermediate operations
}
pair<double, double> Schedule::locSearchEvaluateJobRightShift(size_t idxJob, size_t idxStg, double time, std::vector<std::vector<std::pair<double, double>>>& possibleRightShifts)
{
	cout << "Schedule::locSearchEvaluateJobRightShift(...) not yet implemented." << endl;
	Job* job = scheduledJobs[idxJob].get();
	possibleRightShifts = vector<vector<pair<double, double>>>(size());
	// 1) at each stage check the possible right shifting disregarding maximal time lags (time constraints) and overlaps
	for (size_t o = 0; o < job->size(); ++o) {
		double delay = 0;
		if (o == idxStg) {
			delay = time;
		}
		possibleRightShifts[o] = getRightShiftOptions(&(*job)[o], delay);
	}

	// 2) remove mutually exclusive options
	vector<double> currentLeeway = getLeeway(job);
	vector<vector<pair<size_t, double>>> currentTcSlack = getTcSlack(job);
	bool possibleOverlaps = false;
	bool possibleTcViolations = false;

	// TODO eine Schleife statt zweier

	do {
		possibleOverlaps = constrainRightShiftOptionsFromOverlaps(possibleRightShifts, currentLeeway);
		//possibleTcViolations = constrainRightShiftOptionsFromTimeConstraints(possibleRightShifts, currentTcSlack);
	} while (possibleOverlaps);
		
	do {
		possibleTcViolations = constrainRightShiftOptionsFromTimeConstraints(possibleRightShifts, currentTcSlack);
	} while(possibleTcViolations);

	// 3) evaluate best option (== smallest right shift)
	
	return std::pair<double, double>();
}


pair<double, double> Schedule::locSearchEvaluateBatchLeftShift(Batch* batch, double time, bool& possible) {
	cout << "Schedule::locSearchEvaluateBatchLeftShift not yet implemented." << endl;
	pair<double, double> evaluation = make_pair(0, 0);

	for (size_t i = 0; i < batch->size(); ++i) {
		Operation* op = &(*batch)[i];
		size_t stgIdx = op->getStg() - 1;
		vector<vector<pair<size_t, double>>> tcSlack = getTcSlack(op->getJob());
		possible = true;
		for (size_t tc = 0; tc < tcSlack[stgIdx].size(); ++tc) {
			if (tcSlack[stgIdx][tc].second < time) {
				possible = false;
				break;
			}
		}
	}

	return evaluation;
}
pair<double, double> Schedule::locSearchEvaluateBatchRightShift(Batch* batch, double time, bool& possible) {
	cout << "Schedule::locSearchEvaluateBatchRightShift not yet implemented." << endl;
	pair<double, double> evaluation = make_pair(0, 0);
	for (size_t i = 0; i < batch->size(); ++i) {
		Operation* op = &(*batch)[i];
		
	
	}
	
	return evaluation;
}

double Schedule::locSearchEvaluateOpConsolidation(size_t idxJob, size_t idxStg, bool& feasible) {
	cout << "Schedule::locSearchEvaluateOpConsolidation() not yet implemented." << endl;
	double evaluation = 0.0;
	Job* job = scheduledJobs[idxJob].get();
	Operation* myOp = &(*job)[idxStg];
	Workcenter* myWc = workcenters[idxStg].get();
	if (myWc->getCap() > 2) {
		Batch* myBatch = myOp->getBatch();
		Machine* myMachine = myBatch->getMachine();
		size_t myBatchIdx = myBatch->getIdx();
		int requiredCap = myBatch->getCap() - myBatch->getAvailableCap();

		// consider consolidation of ops at the same machine
		if (myBatchIdx < myMachine->size() - 1) {
			Batch* succBatch = &(*myMachine)[myBatchIdx + 1];
			if (succBatch->getF() == myOp->getF() && succBatch->getAvailableCap() >= requiredCap) {
				bool bMovePossible = false;
				double succBatchR = succBatch->getR();
				double succBatchW = succBatch->getW();
				double delayOfMyOp = max(0, succBatchR - myOp->getStart());
				if (delayOfMyOp > 0) {
					// CONSOLIDATION POSSIBLE (LEAVE CASE delay <= 0 TO LOCAL SEARCH OPERATION LEFT SHIFT)
					double leftShift = succBatch->getStart() - succBatchR;
					pair<double, double> addedTime = locSearchEvaluateBatchRightShift(myBatch, delayOfMyOp, bMovePossible);
					if (bMovePossible) {
						pair<double, double> reducedTime = locSearchEvaluateBatchLeftShift(succBatch, leftShift, bMovePossible);
					}
	
				}
			}
		}

		// TODO consolidate ops at different machines
	}


	


	return evaluation;
}

bool Schedule::locSearchSwapJobs(size_t idxFirst, size_t idxSecond) {
	Job* job1 = scheduledJobs[idxFirst].get();
	Job* job2 = scheduledJobs[idxSecond].get();
	size_t mIdx1 = 0;
	size_t bIdx1 = 0;
	size_t jIdx1 = 0;
	size_t mIdx2 = 0;
	size_t bIdx2 = 0;
	size_t jIdx2 = 0;
	for (size_t o = 0; o < workcenters.size(); ++o) {
		Workcenter* wc = workcenters[o].get();
		Operation* op1 = job1->getOpPtr(o);
		Operation* op2 = job2->getOpPtr(o);
		if (!wc->locateOp(op1, mIdx1, bIdx1, jIdx1) || !wc->locateOp(op2, mIdx2, bIdx2, jIdx2)) {
			return false;
		}
		wc->swapOps(mIdx1, bIdx1, jIdx1, mIdx2, bIdx2, jIdx2);
	}
}

bool Schedule::locSearchLeftShiftJob(size_t jobIdx, vector<vector<pair<double, double>>>& options) {
	Job* job = scheduledJobs[jobIdx].get();
	for (size_t o = 0; o < job->size(); ++o) {
		executeLeftShiftOption(jobIdx, o, options[o][0]);	// first option must be largest shift
	}
	return true;
}

vector<pair<double, double>> Schedule::getLeftShiftOptions(Operation* op) {
	vector<pair<double, double>> options = vector<pair<double, double>>();
	options.push_back(make_pair(0.0, 0.0));
	
	double earliest = op->getRconsideringRawP();

	int stageIdx = op->getStg() - 1;
	int nMachines = this->getWorkcenters()[stageIdx]->size();
	for (size_t i = 0; i < nMachines; ++i) {
		int nBatches = (*this->getWorkcenters()[stageIdx])[i].size();
			// +++ options in existing batches (DISCRETE => option.first == option.second) +++
			for (size_t j = 0; j < nBatches; ++j) {
			double currentBatchStart = (*this->getWorkcenters()[stageIdx])[i][j].getStart();

			if (currentBatchStart < earliest) {
				continue;	// op cannot go into this batch, continue with next batches
			}

			if (currentBatchStart >= op->getStart()) {
				break;		// no favorable batches to insert the batch at this machine, continue with next machine
			}

			if((*this->getWorkcenters()[stageIdx])[i][j].getF() == op->getF() && (*this->getWorkcenters()[stageIdx])[i][j].getAvailableCap() >= op->getS()) {
				double tempOption = op->getStart() - (*this->getWorkcenters()[stageIdx])[i][j].getStart();
				options.push_back(make_pair(tempOption, tempOption));
			}
		}

		bool bOnlyOpInBatch = op->getBatch()->size() == 1;

		// +++ free time slots (CONTINUOUS => option.first <= option.second) +++
		if (nBatches == 0) {
			double tempOptionFrom = op->getStart() - earliest;
			double tempOptionTill = 0.0;
			options.push_back(make_pair(tempOptionFrom, tempOptionTill));
		}
		
		
		for (size_t j = 0; j < nBatches; ++j) {
			double currentBatchStart = (*this->getWorkcenters()[stageIdx])[i][j].getStart();
			if (j == 0) {
				// before 1st batch
				if (earliest + op->getP() <= currentBatchStart || (bOnlyOpInBatch && (*this->getWorkcenters()[stageIdx])[i][j].contains(op))) {
					double tempOptionFrom = op->getStart() - earliest;
					double tempOptionTill = op->getStart() - (currentBatchStart - op->getP());
					if ((bOnlyOpInBatch && (*this->getWorkcenters()[stageIdx])[i][j].contains(op))) {
						tempOptionTill = op->getStart() - currentBatchStart;
					}
					options.push_back(make_pair(tempOptionFrom, tempOptionTill));
				}
			}
			else {
				// inbetween batches
				double gapFrom = max((*this->getWorkcenters()[stageIdx])[i][j - 1].getC(), earliest);
				if (gapFrom + op->getP() <= currentBatchStart || (bOnlyOpInBatch && (*this->getWorkcenters()[stageIdx])[i][j].contains(op))) {
					double tempOptionFrom = op->getStart() - gapFrom;
					double tempOptionTill = op->getStart() - (currentBatchStart - op->getP());
					if ((bOnlyOpInBatch && (*this->getWorkcenters()[stageIdx])[i][j].contains(op))) {
						tempOptionTill = op->getStart() - currentBatchStart;
					}
					if (tempOptionFrom > 0.0) {
						options.push_back(make_pair(tempOptionFrom, tempOptionTill));
					}
				}
			}
		}
		// after last batch
		if (nBatches > 0) {
			if ((*this->getWorkcenters()[stageIdx])[i][nBatches - 1].getC() < op->getStart()) {
				double gapFrom = max(earliest, (*this->getWorkcenters()[stageIdx])[i][nBatches - 1].getC());
				double tempOptionFrom = op->getStart() - gapFrom;
				double tempOptionTill = 0.0;
				options.push_back(make_pair(tempOptionFrom, tempOptionTill));
			}
		}
	}

	// NON-INCREASING ORDER
	sort(options.begin(), options.end(), [&](pair<double, double> a, pair<double, double> b) {
		return a.first > b.first;
		});
	
	return options;
}
vector<pair<double, double>> Schedule::getLeftShiftOptions(Batch* batch) {
	vector<pair<double, double>> options = vector<pair<double, double>>();
	options.push_back(make_pair(0.0, 0.0));
	if (batch->size() == 0) {
		return options;
	}

	double earliest = batch->getRconsideringRawP();
	size_t myBatIdx = batch->getIdx();
	size_t myMacIdx = batch->getMachine()->getIdx();
	size_t myStgIdx = (*batch)[0].getStg() - 1;
	int nMachines = this->getWorkcenters()[myStgIdx]->size();
	for (size_t i = 0; i < nMachines; ++i) {
		int nBatches = (*workcenters[myStgIdx])[i].size();
		if (nBatches == 0) {
			double tempOptionFrom = batch->getStart() - earliest;
			double tempOptionTill = 0.0;
			options.push_back(make_pair(tempOptionFrom, tempOptionTill));
		}

		for (size_t j = 0; j < nBatches; ++j) {
			if (i != myMacIdx || j != myBatIdx) {
				double currentBatchStart = (*this->getWorkcenters()[myStgIdx])[i][j].getStart();
				if (j == 0) {
					// before 1st batch
					if (earliest + batch->getP() <= currentBatchStart) {
						double tempOptionFrom = batch->getStart() - earliest;
						double tempOptionTill = batch->getStart() - (currentBatchStart - batch->getP());
						options.push_back(make_pair(tempOptionFrom, tempOptionTill));
					}
				}
				else {
					// inbetween batches
					double gapFrom = max((*this->getWorkcenters()[myStgIdx])[i][j - 1].getC(), earliest);
					if (gapFrom + batch->getP() <= currentBatchStart) {
						double tempOptionFrom = batch->getStart() - gapFrom;
						double tempOptionTill = batch->getStart() - (currentBatchStart - batch->getP());
						if (tempOptionFrom > 0.0) {
							options.push_back(make_pair(tempOptionFrom, tempOptionTill));
						}
					}
				}
			}
		}

		// after last batch
		if (nBatches > 0) {
			if ((*this->getWorkcenters()[myStgIdx])[i][nBatches - 1].getC() < batch->getStart()) {
				double gapFrom = max(earliest, (*this->getWorkcenters()[myStgIdx])[i][nBatches - 1].getC());
				double tempOptionFrom = batch->getStart() - gapFrom;
				double tempOptionTill = 0.0;
				options.push_back(make_pair(tempOptionFrom, tempOptionTill));
			}
		}
	}
	return options;
}
vector<pair<double, double>> Schedule::getRightShiftOptions(Operation* op, double minDelay) {
	double lastStartAddendum = 10;
	vector<pair<double, double>> options = vector<pair<double, double>>();
	if (minDelay == 0) {
		options.push_back(make_pair(0.0, 0.0));
	}
	double newStart = op->getStart() + minDelay;
	int stageIdx = op->getStg() - 1;
	int nMachines = this->getWorkcenters()[stageIdx]->size();
	for (size_t i = 0; i < nMachines; ++i) {
		int nBatches = (*this->getWorkcenters()[stageIdx])[i].size();
		// +++ options in existing batches (DISCRETE => option.first == option.second) +++
		for (size_t j = 0; j < nBatches; ++j) {
			double currentBatchStart = (*this->getWorkcenters()[stageIdx])[i][j].getStart();

			if (currentBatchStart < newStart) {
				continue;	
			}

			if ((*this->getWorkcenters()[stageIdx])[i][j].getF() == op->getF() && (*this->getWorkcenters()[stageIdx])[i][j].getAvailableCap() >= op->getS()) {
				double tempOption = (*this->getWorkcenters()[stageIdx])[i][j].getStart() - op->getStart();
				options.push_back(make_pair(tempOption, tempOption));
			}
		}

		bool bOnlyOpInBatch = op->getBatch()->size() == 1;

		// +++ free time slots (CONTINUOUS => option.first <= option.second) +++
		if (nBatches == 0) {
			double tempOptionFrom = minDelay;
			double tempOptionTill = minDelay + (lastStartAddendum * op->getP());	// to infinity and beyond!
			options.push_back(make_pair(tempOptionFrom, tempOptionTill));
		}

		for (size_t j = 0; j < nBatches; ++j) {
			double currentBatchStart = (*this->getWorkcenters()[stageIdx])[i][j].getStart();
			if (j == 0) {
				// before 1st batch
				if (newStart + op->getP() <= currentBatchStart) {
					double tempOptionFrom = minDelay;
					double tempOptionTill = (currentBatchStart - op->getP()) - op->getStart();
					options.push_back(make_pair(tempOptionFrom, tempOptionTill));
				}
			} else {
				// inbetween batches
				double gapFrom = max((*this->getWorkcenters()[stageIdx])[i][j].getC(), newStart);
				if (gapFrom + op->getP() <= currentBatchStart) {
					double tempOptionFrom = gapFrom - op->getStart();
					double tempOptionTill = (currentBatchStart - op->getP()) - op->getStart();
					options.push_back(make_pair(tempOptionFrom, tempOptionTill));
				}
 			}
		}
		
		// after last batch
		if (nBatches > 0) {
			double tempOptionFrom = max((*this->getWorkcenters()[stageIdx])[i][nBatches-1].getC(), newStart) - op->getStart();
			double tempOptionTill = tempOptionFrom + (lastStartAddendum * op->getP());	// to infinity and beyond!
			options.push_back(make_pair(tempOptionFrom, tempOptionTill));
		}
	}
	
	// NON-DECREASING ORDER
	sort(options.begin(), options.end(), [&](pair<double, double> a, pair<double, double> b) {
		return a.first < b.first;
		});
	return options;
}
vector<double> Schedule::getLeeway(Job* job) {
	vector<double> leeway = vector<double>(job->size());
	leeway[0] = (*job)[0].getStart() - job->getR();
	for (size_t i = 1; i < job->size(); ++i) {
		leeway[i] = (*job)[i].getStart() - (*job)[i - 1].getC();
	}
	return leeway;
}
vector<double> Schedule::getRightSideLeeway(Job* job) {
	vector<double> rsLeeway = vector<double>(job->size());
	for (size_t i = 0; i < job->size() - 1; ++i) {
		rsLeeway[i] = (*job)[i + 1].getStart() - (*job)[i].getC();
	}
	rsLeeway[job->size() - 1] = 100 * (*job)[job->size() - 1].getP(); // infinity!
	return rsLeeway;
}
vector<vector<pair<size_t, double>>> Schedule::getTcSlack(Job* job) {
	vector<vector<pair<size_t, double>>> tcSlack = vector<vector<pair<size_t, double>>>(job->size());
	for (size_t i = 0; i < job->size(); ++i) {
		vector<pair<int, double>> tcFwd = job->getTcMaxFwd(i);
		tcSlack[i] = vector<pair<size_t, double>>();
		for (size_t j = 0; j < tcFwd.size(); ++j) {
			if (tcFwd[j].first > i && tcFwd[j].second < 999999) {
				double tempTcSlack = ((*job)[i].getStart() + tcFwd[j].second) - (*job)[tcFwd[j].first].getStart();
				tcSlack[i].push_back(make_pair(tcFwd[j].first, tempTcSlack));
			}
		}
	}
	return tcSlack;
}
bool Schedule::constrainLeftShiftOptionsFromOverlaps(std::vector<std::vector<std::pair<double, double>>>& options, std::vector<double>& leeway) {
	bool changeApplied = false;
	for (int o = options.size() - 1; o >= 1; --o) {
		for (int opt = options[o].size() - 1; opt >= 0; --opt) {
			double optFrom = options[o][opt].first;
			double optTill = options[o][opt].second;
						
			// 1) delete discrete options if infeasible
			if (optFrom == optTill) {
				bool feasible = false;
				for (size_t prevOpt = 0; prevOpt < options[o - 1].size(); ++prevOpt) {
					if (optFrom <= options[o - 1][prevOpt].first + leeway[o]) {
						feasible = true;
						break;
					}
				}
				if (!feasible) {
					changeApplied = true;
					options[o].erase(options[o].begin() + opt);
				}
			}
			
			// 2) reduce continuous options if necessary
			if (optFrom > optTill) {
				double maxLeftShift = 0.0;
				for (size_t prevOpt = 0; prevOpt < options[o - 1].size(); ++prevOpt) {
					maxLeftShift = max(maxLeftShift, options[o - 1][prevOpt].first + leeway[o]);
				}
				
				if (maxLeftShift < options[o][opt].first) {
					changeApplied = true;
					if (maxLeftShift < options[o][opt].second) {
						options[o].erase(options[o].begin() + opt);
					} else {
						options[o][opt].first = maxLeftShift;
					}
				}
			}
		}
	}
	return changeApplied;
}
bool Schedule::constrainLeftShiftOptionsFromTimeConstraints(std::vector<std::vector<std::pair<double, double>>>& options, vector<vector<pair<size_t, double>>>& tcSlack) {
	bool changeApplied = false;
	for (int o = 0; o < options.size() - 1; ++o) {
		for (int opt = 0; opt < options[o].size(); ++opt) {
			double optFrom = options[o][opt].first;
			double optTill = options[o][opt].second;

			// 1) reduce discrete options if infeasible
			if (optFrom == optTill) {
				bool allFeasible = true;
				bool tempFeasible = false;
				if (tcSlack[o].size() < 1) { 
					tempFeasible = true;
				}
				for (size_t tc = 0; tc < tcSlack[o].size(); ++tc) {
					tempFeasible = false;
					size_t succO = tcSlack[o][tc].first;	// successive stage connected by a time constraint
					for (size_t succOpt = 0; succOpt < options[succO].size(); ++succOpt) {
						if (options[o][opt].first <= options[succO][succOpt].first + tcSlack[o][tc].second) {
							tempFeasible = true;
							break;
						}
					}
					if (!tempFeasible) {
						allFeasible = false;
					} 
				}
				
				if (!allFeasible) {
					changeApplied = true;
					options[o].erase(options[o].begin() + opt);
				}
			}

			// 2) reduce continuous options if necessary
			if (optFrom > optTill) {
				double maxLeftShift = 0;
				for (size_t tc = 0; tc < tcSlack[o].size(); ++tc) {
					size_t succO = tcSlack[o][tc].first;
					for (size_t succOpt = 0; succOpt < options[succO].size(); ++succOpt) {
						maxLeftShift = max(maxLeftShift, options[succO][succOpt].first + tcSlack[o][tc].second);
					}
					if (maxLeftShift < options[o][opt].first) {
						changeApplied = true;
						if (maxLeftShift < options[o][opt].second) {
							options[o].erase(options[o].begin() + opt);
						} else {
							options[o][opt].first = maxLeftShift;
						}
					}
				}
			}

		}
	}
	return changeApplied;
}

bool Schedule::constrainRightShiftOptionsFromOverlaps(std::vector<std::vector<std::pair<double, double>>>& options, std::vector<double>& leeway) {
	bool changeApplied = false;
	for (int o = options.size() - 1; o >= 1; --o) {
		for (int opt = options[o].size() - 1; opt >= 0; --opt) {
			double optFrom = options[o][opt].first;
			double optTill = options[o][opt].second;

			// 1) delete discrete options if infeasible
			if (optFrom == optTill) {
				bool feasible = false;
				for (size_t prevOpt = 0; prevOpt < options[o - 1].size(); ++prevOpt) {
					if (optFrom >= options[o - 1][prevOpt].first + leeway[o]) {
						feasible = true;
						break;
					}
				}
				if (!feasible) {
					changeApplied = true;
					options[o].erase(options[o].begin() + opt);
				}
			}

			// 2) reduce continuous options if necessary
			if (optFrom < optTill) {
				double minRightShift = DBL_MAX;
				for (size_t prevOpt = 0; prevOpt < options[o - 1].size(); ++prevOpt) {
					minRightShift = min(minRightShift, options[o - 1][prevOpt].first + leeway[o]);
				}

				if (minRightShift > options[o][opt].first) {
					changeApplied = true;
					if (minRightShift > options[o][opt].second) {
						options[o].erase(options[o].begin() + opt);
					}
					else {
						options[o][opt].first = minRightShift;
					}
				}
			}
		}
	}
	return changeApplied;
}
bool Schedule::constrainRightShiftOptionsFromTimeConstraints(std::vector<std::vector<std::pair<double, double>>>& options, std::vector<std::vector<std::pair<size_t, double>>>& tcSlack) {
	bool changeApplied = false;
	for (int o = 0; o < options.size() - 1; ++o) {
		for (int opt = 0; opt < options[o].size(); ++opt) {
			double optFrom = options[o][opt].first;
			double optTill = options[o][opt].second;

			// 1) reduce discrete options if infeasible
			if (optFrom == optTill) {
				bool allFeasible = true;
				bool tempFeasible = false;
				if (tcSlack[o].size() < 1) {
					tempFeasible = true;
				}
				for (size_t tc = 0; tc < tcSlack[o].size(); ++tc) {
					tempFeasible = false;
					size_t succO = tcSlack[o][tc].first;	// successive stage connected by a time constraint
					for (size_t succOpt = 0; succOpt < options[succO].size(); ++succOpt) {
						if (options[o][opt].first >= options[succO][succOpt].first - tcSlack[o][tc].second) {
							tempFeasible = true;
							break;
						}
					}
					if (!tempFeasible) {
						allFeasible = false;
					}
				}

				if (!allFeasible) {
					changeApplied = true;
					options[o].erase(options[o].begin() + opt);
				}
			}

			// 2) reduce continuous options if necessary
			if (optFrom < optTill) {
				double minRightShift = DBL_MAX;
				for (size_t tc = 0; tc < tcSlack[o].size(); ++tc) {
					size_t succO = tcSlack[o][tc].first;
					for (size_t succOpt = 0; succOpt < options[succO].size(); ++succOpt) {
						minRightShift = min(minRightShift, options[succO][succOpt].first - tcSlack[o][tc].second);
					}
					if (minRightShift > options[o][opt].first) {
						changeApplied = true;
						if (minRightShift > options[o][opt].second) {
							options[o].erase(options[o].begin() + opt);
						}
						else {
							options[o][opt].first = minRightShift;
						}
					}
				}
			}

		}
	}
	return changeApplied;
}



void Schedule::executeLeftShiftOption(size_t jobIdx, size_t stgIdx, std::pair<double, double>& option) {
	if (option.first > 0) {
		Operation* operation = &(*scheduledJobs[jobIdx])[stgIdx];
		bool bIntoExistingBatch = option.first == option.second;
		double newStart = operation->getStart() - option.first;
		if (!workcenters[stgIdx]->moveOpDisregardingTc(operation, newStart, bIntoExistingBatch)) {
			throw(ExcSched("ERROR: Schedule::executeLeftShiftOption(...) invalid."));
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

void Schedule::debugSetR(size_t scheduledJobIdx, double newR) {
	cout << "+++ WARNING: INFEASIBLE OPERATION EXECUTED IN COURSE OF DEBUGGING +++" << endl;
	scheduledJobs[scheduledJobIdx]->setR(newR);
}
void Schedule::debugAddMachine(size_t stgIdx) {
	unique_ptr<Machine> newMac = make_unique<Machine>((*workcenters[stgIdx])[workcenters[stgIdx]->size()-1].getId() + 1, workcenters[stgIdx]->getCap(), &(*workcenters[stgIdx]));
	workcenters[stgIdx]->addMachine(move(newMac));
}

