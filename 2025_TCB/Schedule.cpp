#include <iostream>

#include "Schedule.h"
#include "Functions.h"

using namespace std;

Schedule::Schedule() {
	workcenters = vector<pWc>();
	jobs = vector<pJob>();
	scheduledJobs = vector<pJob>();
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

Job& Schedule::getJob(size_t idx) { return *jobs[idx]; };
Job& Schedule::getJob(size_t idx) const { return *jobs[idx]; }

std::unique_ptr<Schedule> Schedule::clone() const {
	auto newSchedule = make_unique<Schedule>();
	for (const auto& wc : workcenters) {
		newSchedule->addWorkcenter(wc->clone(newSchedule.get()));
	}
	for (const auto& job : jobs) {
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
					(*this)[wc][m][b].addOp(&(*jobs[wc])[jobIdx]);
				}
			}
		}
	}
}

size_t Schedule::size() const { return workcenters.size();  }
size_t Schedule::getN() const { return jobs.size(); }

const std::vector<pWc>& Schedule::getWorkcenters() const {
	return workcenters;
}
void Schedule::addWorkcenter(pWc wc) {
	workcenters.push_back(move(wc));
}
void Schedule::addJob(pJob job) {
	jobs.push_back(move(job));
}

void Schedule::schedOp(Operation* op, double pWait) {
	int wcIdx = op->getWorkcenterId() - 1;
	workcenters[wcIdx]->schedOp(op, pWait);
}

void Schedule::lSchedJobs(std::vector<pJob>& unscheduled, double pWait) {
	while(!unscheduled.empty()) {
		for (size_t op = 0; op < (*unscheduled.begin())->size(); ++op) {
			schedOp(&(**unscheduled.begin())[op], pWait);
		}
		shiftJobFromVecToVec(unscheduled, scheduledJobs, 0);
	}
}

double Schedule::getTWT() const {
	double twt = 0;
	for (size_t wc = 0; wc < size(); ++wc) {
		twt += workcenters[wc]->getTWT();
	}
	return twt;
}