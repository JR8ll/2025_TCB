#include "Operation.h"
#include "Job.h"

Operation::Operation(Job* j, int stg) : job(j), id(j->getId()), stg(stg), wait(0), batch(nullptr), pred(nullptr), succ(nullptr) {} 

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
