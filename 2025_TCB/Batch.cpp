#include <iostream>

#include "Batch.h"
#include "Machine.h"
#include "Operation.h"

using namespace std;

Batch::Batch() : machine(nullptr) {}

ostream& operator<<(ostream& os, const Batch& batch) {
	os << "S" << (int)batch.getStart() << "_C" << (int)batch.getC() << "[";
	for (size_t op = 0; op < batch.size() - 1; ++op) {
		os << batch[op] << ",";
	}
	if (!batch.isEmpty()) {
		os << batch[batch.size() - 1];
	}
	os << "]";
	return os;
}

unique_ptr<Batch> Batch::clone() const {
	auto newBatch = make_unique<Batch>();
	// shallow copy
	return newBatch;
}

Operation& Batch::operator[](size_t idx) { return *ops[idx]; }
Operation& Batch::operator[](size_t idx) const { return *ops[idx]; }

size_t Batch::size() const { return ops.size(); }
bool Batch::isEmpty() const { return ops.empty(); }
double Batch::getStart() const { return start; }
double Batch::getC() const { return c; }

int Batch::getF() const { return f; }
int Batch::getCap() const { return cap; }
int Batch::getAvailableCap() const { 
	int availableCap = cap;
	for (size_t op = 0; op < size(); ++op) {
		availableCap -= ops[op]->getS();
	}
	return availableCap;
}

const vector<Operation*>& Batch::getOps() const {
	return ops;
}

void Batch::assignToMachine(Machine* processor) { machine = processor; };

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

double Batch::getTWT() const {
	double twt = 0;
	for (size_t op = 0; op < size(); ++op) {
		twt += ops[op]->getTWT();
	}
	return twt;
}