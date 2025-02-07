#include "Batch.h"
#include "Machine.h"
#include "Operation.h"

using namespace std;

const vector<Operation*>& Batch::getOps() const {
	return ops;
}

void Batch::addOp(Operation* op) {
	ops.push_back(op);
	op->assignToBatch(this);
}
void Batch::removeOp(Operation* op) {
	auto it = std::find(ops.begin(), ops.end(), op);
	if (it != ops.end()) {
		ops.erase(it);
		op->assignToBatch(nullptr);
	}
}