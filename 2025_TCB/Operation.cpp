#include <algorithm>
#include <iostream>

#include "Operation.h"
#include "Batch.h"
#include "Job.h"

using namespace std;

Operation::Operation(Job* j, int stg) : job(j), id(j->getId()), stg(stg), wait(0), batch(nullptr), pred(nullptr), succ(nullptr) {} 

ostream& operator<<(ostream& os, const Operation& op) {
	os << op.getId() << "." << op.getStg();
	return os;
}

int Operation::getId() const { return id; }
int Operation::getStg() const { return stg; }

double Operation::getP() const { return job->getP(stg); }
double Operation::getWait() const { return wait; }

void Operation::setWait(double wt) { wait = wt; }
void Operation::setPred(Operation* pre) { pred = pre; }
void Operation::setSucc(Operation* suc) { succ = suc; }

Job* Operation::getJob() const { return job; }
Batch* Operation::getBatch() const { return batch; }

void Operation::assignToBatch(Batch* batch) {
	batch = batch;
}

double Operation::getTWT() const {
	double twt = 0;
	if (succ != nullptr && batch != nullptr) {
		twt = max(0.0, batch->getC() - job->getD()) * job->getW();
	}
	return twt;
}
