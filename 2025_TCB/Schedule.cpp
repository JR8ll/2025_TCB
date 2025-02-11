#include <algorithm>
#include <iostream>

#include "Schedule.h"
#include "Functions.h"

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
	newSchedule->_reconstruct(this);
	return newSchedule;
}

void Schedule::_reconstruct(const Schedule* orig) {
	for (size_t wc = 0; wc < (*orig).size(); ++wc) {
		for (size_t m = 0; m < (*orig)[wc].size(); ++m) {
			for (size_t b = 0; b < (*orig)[wc][m].size(); ++b) {
				(*this)[wc][m].addBatch(move((*orig)[wc][m][b].clone()), (*orig)[wc][m][b].getStart());
				for (size_t j = 0; j < (*orig)[wc][m][b].size(); ++j) {
					Operation& op = (*orig)[wc][m][b][j];
					size_t jobIdx = op.getStg() - 1;
					(*this)[wc][m][b].addOp(&(*unscheduledJobs[wc])[jobIdx]);
				}
			}
		}
	}
}

size_t Schedule::size() const { return workcenters.size();  }
size_t Schedule::getN() const { return unscheduledJobs.size(); }

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

void Schedule::markAsScheduled(size_t jobIdx) {
	if (jobIdx >= unscheduledJobs.size()) throw out_of_range("Schedule::markAsScheduled() out of range");
	shiftJobFromVecToVec(unscheduledJobs, scheduledJobs, jobIdx);
}
void Schedule::markAsScheduled(pJob scheduledJob) {
	scheduledJobs.push_back(move(scheduledJob));
}

void Schedule::lSchedFirstJob(double pWait) {
	for (size_t op = 0; op < (*unscheduledJobs.begin())->size(); ++op) {
		schedOp(&(**unscheduledJobs.begin())[op], pWait);
	}
	shiftJobFromVecToVec(unscheduledJobs, scheduledJobs, 0);
}
void Schedule::lSchedJobs(double pWait) {
	while(!unscheduledJobs.empty()) {
		lSchedFirstJob(pWait);
	}
}

void Schedule::lSchedJobsWithSorting(prioRule<pJob> rule, double pWait) {
	rule(unscheduledJobs);
	lSchedJobs(pWait);
}

void Schedule::lSchedJobsWithSorting(prioRuleKappa<pJob> rule, double kappa, double pWait) {
	double t = 0.0;	// dynamic computation of priority index (increase t)
	while (!unscheduledJobs.empty()) {
		t = getMinMSP(0);
		rule(unscheduledJobs, t, kappa);
		lSchedFirstJob(pWait);
	}
}

void Schedule::lSchedJobsWithSorting(prioRuleKappa<pJob> rule, const std::vector<double>& kappaGrid, double pWait) {
	double bestTWT = numeric_limits<double>::max();
	double bestKappa = 0.0;
	for (size_t kappa = 0; kappa < kappaGrid.size(); ++kappa) {
		lSchedJobsWithSorting(rule, kappaGrid[kappa], pWait);
		double tempTWT = getTWT();
		if (tempTWT < bestTWT) {
			bestTWT = tempTWT;
			bestKappa = kappaGrid[kappa];
		}
		reset();
	}
	lSchedJobsWithSorting(rule, bestKappa, pWait);
	TCB::logger.Log(Info, "Found a schedule with best kappa value = " + to_string(bestKappa));
}

void Schedule::lSchedJobsWithRandomKeySorting(prioRuleKeySet<pJob> rule, const std::vector<double>& keys, double pWait) {
	rule(unscheduledJobs, keys);
	lSchedJobs(pWait);
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
