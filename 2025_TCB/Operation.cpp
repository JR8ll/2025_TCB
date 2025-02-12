#include <algorithm>
#include <iostream>

#include "Operation.h"
#include "Batch.h"
#include "Functions.h"
#include "Job.h"
#include "Machine.h"
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
		Operation* tcSucc = job->getOpPtr(tcFwd[tc].first);
		if (tcSucc->isScheduled()) {
			double constrainedStart = tcSucc->getStart() - tcFwd[tc].second;
			if (constrainedStart > earliest) {
				earliest = constrainedStart;
			}
		}
	}
	return earliest;
}
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

const std::vector<std::pair<int, double>>& Operation::getTcMaxBwd() const {
	return job->getTcMaxBwd(stg-1);
}
const std::vector<std::pair<int, double>>& Operation::getTcMaxFwd() const {
	return job->getTcMaxFwd(stg-1);
}

void Operation::setWait(double wt) { wait = wt; }
void Operation::setPred(Operation* pre) { pred = pre; }
void Operation::setSucc(Operation* suc) { succ = suc; }



Job* Operation::getJob() const { return job; }
Batch* Operation::getBatch() const { return batch; }

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
						wcSucc->rightShift(mIdx, bIdx, jIdx, batch->getC());
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
							wc->rightShift(mIdx, bIdx, jIdx, newStartForPred);
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
