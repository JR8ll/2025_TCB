#include "Batch.h"
#include "Machine.h"
#include "Operation.h"

using namespace std;

unique_ptr<Batch> Batch::clone() const {
	auto newBatch = make_unique<Batch>();
	// shallow copy
	return newBatch;
}

Operation& Batch::operator[](size_t idx) { return *ops[idx]; }
Operation& Batch::operator[](size_t idx) const { return *ops[idx]; }

size_t Batch::size() const { return ops.size(); }

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