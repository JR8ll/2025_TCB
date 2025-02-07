#include "Operation.h"

Operation::Operation(Job* j) : job(j), batch(nullptr) {} 

Job* Operation::getJob() const { return job; }
Batch* Operation::getBatch() const { return batch; }

void Operation::assignToBatch(Batch* batch) {
	batch = batch;
}
