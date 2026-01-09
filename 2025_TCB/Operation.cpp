#include <algorithm>
#include <iostream>

#include "Operation.h"
#include "Batch.h"
#include "Functions.h"
#include "Job.h"
#include "Machine.h"
#include "Schedule.h"
#include "Workcenter.h"

using namespace std;

Operation::Operation(Job* j, int stg) : job(j), id(j->getId()), stg(stg), wait(0), batch(nullptr), pred(nullptr), succ(nullptr) {} 
Operation::~Operation() {
	resetLinks();
}

ostream& operator<<(ostream& os, const Operation& op) {
	os << op.getId() << "." << op.getStg();
	return os;
}

int Operation::getId() const { return id; }
size_t Operation::getIdxInBatch() const {
	return batch->findOp(this);
}
int Operation::getStg() const { return stg; }
int Operation::getWorkcenterId() const { return job->getWorkcenterId(stg-1); };

int Operation::getS() const { return job->getS(); }
int Operation::getF() const { return job->getF(); }
double Operation::getStart() const { 
	if (batch != nullptr) {
		if (batch->getMachine() != nullptr) {
			return batch->getStart();
		}
	}
	return numeric_limits<double>::max();	// operation has not yet been scheduled
}
double Operation::getC() const { 
	if (batch != nullptr) {
		if (batch->getMachine() != nullptr) {
			return batch->getC();
		}
	}
	return numeric_limits<double>::max();	// operation has not yet been scheduled
}
double Operation::getD() const { return job->getD(); }
double Operation::getP() const { 
	return job->getP(stg-1); 
}
double Operation::getR() const { return job->getR(); }
double Operation::getRconsideringRawP() const {
	double r = job->getR();
	for (size_t o = 0; o < (this->stg-1); ++o) {
		r += job->getOpPtr(o)->getP();
	}
	return r;
}
double Operation::getW() const { return job->getW(); }

bool Operation::isScheduled() const {
	if (batch != nullptr) {
		return batch->getMachine() != nullptr;
	}
	return false;
}
double Operation::getAvailability() const {
	if (pred == nullptr) {
		return getR();								// first op is available with job release => end recursion
	}
	
	if (pred->isScheduled()) {
		return pred->getC();						// op is available when its predecessor is completed => end recursion
	}
	return pred->getAvailability() + pred->getP();	// op is available when its predecessor could be completed => end recursion 
}
double Operation::getEarliestStart() const {
	
	double earliest = getAvailability();
	
	const vector<pair<int, double>>& tcFwd = getTcMaxFwd();
	for (size_t tc = 0; tc < tcFwd.size(); ++tc) {
		if (tcFwd[tc].second < 999999) {
			Operation* tcSucc = job->getOpPtr(tcFwd[tc].first);
			if (tcSucc->isScheduled()) {
				double constrainedStart = tcSucc->getStart() - tcFwd[tc].second;
				if (constrainedStart > earliest) {
					earliest = constrainedStart;
				}
			}
		}	
	}
	return earliest;
}
double Operation::getEarliestStartForBackwardScheduling() const {
	double earliest = getEarliestStart();
	//TODO take into account overlapping time constraints, i.e. time constraint between successor and predecessors
	Operation* succ = getSucc();
	while(succ != nullptr) {
		const vector<pair<int, double>>& tcBwd = succ->getTcMaxBwd();
		for (size_t tc = 0; tc < stg-1; ++tc) {
			if (tcBwd[tc].second < 999999) {
				double rawP = 0.0;
				for (size_t pred = tcBwd[tc].first; pred < stg-1; ++pred) {	
					rawP += (*job)[pred].getP();
				}
				earliest = max(earliest, succ->getStart() - tcBwd[tc].second + rawP);
			}
		}
		succ = succ->getSucc();
	}
	return earliest;
}
double Operation::getLatestStartConsideringBwdTc() const {
	double latest = DBL_MAX;
	const vector<pair<int, double>>* tcBwd = &getTcMaxBwd();
	for (size_t i = 0; i < tcBwd->size(); ++i) {
		if (!(*tcBwd)[i].second >= 999999) {
			latest = min(latest, (*job)[(*tcBwd)[i].first].getStart() + (*tcBwd)[i].second);
		}
	}
	return latest;
}
//double Operation::getLatestStartConsideringFwdTc() const {
//	double latest = DBL_MAX;
//	const vector<pair<int, double>>* tcFwd = &getTcMaxFwd();
//	for (size_t i = 0; i < tcFwd->size(); ++i) {
//		if (!(*tcFwd)[i].second >= 999999) {
//			latest = min(latest, (*job)[(*tcFwd)[i].first].getStart() - (*tcFwd)[i].second);
//		}
//	}
//	return latest;
//}

double Operation::getWait() const { return wait; }

double Operation::getGATC(double avgP, double t, double kappa) const {
	double slack = getD() - t - getP();
	Operation* next = succ;
	while (next != nullptr) {
		slack -= next->getWait() + next->getP();
		next = next->getSucc();
	}
	return (getW() / getP()) * exp(-1 * (max(slack, 0.0) / kappa * avgP));
}

Operation* Operation::getPred() const { return pred; }
Operation* Operation::getSucc() const { return succ; }
Operation* Operation::getSucc(size_t offsetStages) const {
	Operation* successor = succ;
	for (size_t i = 0; i < offsetStages; ++i) {
		if (successor == nullptr) break;
		successor = successor->getSucc();
	}
	return successor;
}

const std::vector<std::pair<int, double>>& Operation::getTcMaxBwd() const {
	return job->getTcMaxBwd(stg-1);
}
const std::vector<std::pair<int, double>>& Operation::getTcMaxFwd() const {
	return job->getTcMaxFwd(stg-1);
}

void Operation::setWait(double wt) { wait = wt; }
void Operation::computeWaitingTimeFromStart(double start) {
	double earliest = getEarliestStart();
	if (earliest - TCB::precision > start) {
		throw ExcSched("Negative waiting time");
	}
	wait = max(0.0, start - earliest);	// could be zero only by marginal value (precision)
}
void Operation::setPred(Operation* pre) { pred = pre; }
void Operation::setSucc(Operation* suc) { succ = suc; }

Job* Operation::getJob() const { return job; }
Batch* Operation::getBatch() const { return batch; }

bool Operation::checkProcessingOrder() const {
	if (pred != nullptr) {
		if (pred->getC() - TCB::precision > getC() - getP()) return false;
	}
	if (succ != nullptr) {
		if (succ->getStart() + TCB::precision < getC()) return false;
	}
	return true;
}

bool Operation::checkTimeConstraints() const {
	const vector<pair<int, double>>& tc = job->getTcMaxBwd(stg-1);
	for (size_t t = 0; t < tc.size(); ++t) {
		if (tc[t].second != 999999) {
			int steps = (stg-1) - tc[t].first;
			Operation* tcPred = pred;
			for (size_t step = 1; step < steps; ++step) {
				tcPred = pred->getPred();
			}
			if (tcPred->getStart() + tc[t].second + TCB::precision < getStart()) {
				return false;
			}
		}
	}
	return true;
}

void Operation::assignToBatch(Batch* newBatch) {
	batch = newBatch;
}
bool Operation::repairOverlaps() {
	bool bRepaired = false;
	if (batch != nullptr) {
		if (succ != nullptr) {
			if (succ->isScheduled()) {
				Batch* succBatch = succ->getBatch();
				if (succBatch != nullptr) {
					if (batch->getC() > succBatch->getStart() + TCB::precision) {
						Workcenter* wcSucc = succBatch->getMachine()->getWorkcenter();
						int mIdx = succBatch->getMachine()->getIdx();
						int bIdx = succBatch->getIdx();
						int jIdx = succ->getIdxInBatch();
						wcSucc->rightShiftOp(mIdx, bIdx, jIdx, batch->getC());
					}
				}
			}
		}
	}
	return bRepaired;	
}
bool Operation::repairTimeConstraints() {
	bool bRepaired = false;
	const vector<pair<int, double>>& tcMax = getTcMaxBwd();
	if (batch != nullptr) {
		for (size_t i = 0; i < tcMax.size(); ++i) {
			size_t predIdx = tcMax[i].first;
			Operation* predOp = &(*job)[predIdx];
			if (predOp != nullptr) {
				if (predOp->isScheduled()) {
					Batch* predBatch = predOp->getBatch();
					if (predBatch->getStart() < batch->getStart() - tcMax[i].second - TCB::precision) {
						// time constraint is violated
						double newStartForPred = batch->getStart() - tcMax[i].second;
						Workcenter* wc = predBatch->getMachine()->getWorkcenter();
						if (wc != nullptr) {
							int mIdx = predBatch->getMachine()->getIdx();
							int bIdx = predBatch->getIdx();
							int jIdx = predOp->getIdxInBatch();
							wc->rightShiftOp(mIdx, bIdx, jIdx, newStartForPred);
							bRepaired = true;
						}
						else {
							throw(ExcSched("Operation::repairTimeConstraint() missing workcenter reference"));
						}
					}
				}	
			}
		}
	}
	return bRepaired;
}

bool Operation::repairOverlapsFixedBatchFormation() {
	bool bRepaired = false;
	if (batch != nullptr) {
		if (succ != nullptr) {
			if (succ->isScheduled()) {
				Batch* succBatch = succ->getBatch();
				if (succBatch != nullptr) {
					if (batch->getC() > succBatch->getStart() + TCB::precision) {
						Workcenter* wcSucc = succBatch->getMachine()->getWorkcenter();
						int mIdx = succBatch->getMachine()->getIdx();
						int bIdx = succBatch->getIdx();
						wcSucc->rightShiftBatch(mIdx, bIdx, batch->getC(), true);	// if necessary push right succeding batches
						bRepaired = true;
					}
				}
			}
		}
	}
	return bRepaired;
}

bool Operation::repairTimeConstraintsFixedBatchFormation() {
	bool bRepaired = false;
	const vector<pair<int, double>>& tcMax = getTcMaxBwd();
	if (batch != nullptr) {
		for (size_t i = 0; i < tcMax.size(); ++i) {
			size_t predIdx = tcMax[i].first;
			Operation* predOp = &(*job)[predIdx];
			if (predOp != nullptr) {
				if (predOp->isScheduled()) {
					Batch* predBatch = predOp->getBatch();
					if (predBatch->getStart() < batch->getStart() - tcMax[i].second - TCB::precision) {
						// time constraint is violated

						double predAvailability = predBatch->getR();										// [JR-2026-Jan-09] pred->getAvailability changed to predBatch->getR()
						double newStartForPred = max(batch->getStart() - tcMax[i].second, predAvailability);	// [JR-2026-Jan-09] added max(.., getEarliestStart())
						Workcenter* wc = predBatch->getMachine()->getWorkcenter();
						if (wc != nullptr) {
							int mIdx = predBatch->getMachine()->getIdx();
							int bIdx = predBatch->getIdx();

							// DEBUGGING
							if (newStartForPred > 193.8 && newStartForPred < 193.9 && mIdx == 0 && bIdx == 7) {
								wc->getSchedule()->saveJsonFactory("DEBUGGING");
								int debugger = 666;
							}


							wc->rightShiftBatch(mIdx, bIdx, newStartForPred, true, false);		// [JR-2026-Jan-09] checkValidity == false
							bRepaired = true;
						}
						else {
							throw(ExcSched("Operation::repairTimeConstraint() missing workcenter reference"));
						}
					}
				}
			}
		}
	// TODO WORK IN PROGRESS
	}
	return bRepaired;
}

double Operation::getTWT() const {
	double twt = 0;
	if (succ == nullptr && batch != nullptr) {
		twt = max(0.0, batch->getC() - job->getD()) * job->getW();
	}
	return twt;
}

void Operation::resetLinks() {
	job = nullptr;
	batch = nullptr;
	pred = nullptr;
	succ = nullptr;
}
