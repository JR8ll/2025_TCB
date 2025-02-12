#include <iostream>

#include "Batch.h"
#include "Functions.h"
#include "Machine.h"
#include "Operation.h"

using namespace std;

Batch::Batch() : machine(nullptr), cap(0), start(0), f(0), c(0) {}
Batch::Batch(int cap) : machine(nullptr), cap(cap), start(0), f(0), c(0) {}
Batch::~Batch() { ops.clear(); }

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

size_t Batch::getIdx() const {
	return machine->findBatch(this);
}

size_t Batch::size() const { return ops.size(); }
size_t Batch::findOp(const Operation* op) const {
	for (size_t j = 0; j < size(); ++j) {
		if (ops[j] == op) return j;
	}
	throw ExcSched("Batch::findOp() Operation not found");
}
bool Batch::isEmpty() const { return ops.empty(); }
double Batch::getStart() const { return start; }
double Batch::getC() const { return c; }

double Batch::getP() const {
	double p = 0.0;
	for (size_t j = 0; j < size(); ++j) {
		double tempP = ops[j]->getP();
		if (tempP > p) p = tempP;
	}
	return p;
}

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
Machine* Batch::getMachine() const {
	return machine;
}

void Batch::setStart(double newStart, bool checkvalidity) {
	start = newStart;
	setC(newStart + getP(), checkvalidity);
}
void Batch::setC(double newC, bool checkvalidity) {
	if (checkvalidity) {
		for (size_t op = 0; op < size(); ++op) {
			double test = ops[op]->getEarliestStart();
			double test2 = ops[op]->getP();

			if (ops[op]->getEarliestStart() + ops[op]->getP() - TCB::precision > newC) throw ExcSched("Batch::setC() infeasible");
		}
	}
	c = newC;
}

void Batch::setCap(int newCap) {
	int requiredCap = cap - getAvailableCap();
	if (requiredCap > newCap) throw ExcSched("Batch::setCap() infeasible");
	cap = newCap;
}

void Batch::assignToMachine(Machine* processor) { machine = processor; };

bool Batch::addOp(Operation* op) {
	if (f == 0 || f == op->getF()) {
		if (getAvailableCap() >= op->getS()) {
			f = op->getF();
			ops.push_back(op);
			op->assignToBatch(this);
			return true;
			// TODO: for safety maybe consider operations being added with availability > batch.start
		}
	}
	return false;
}
void Batch::removeOp(Operation* op) {
	auto it = std::find(ops.begin(), ops.end(), op);
	if (it != ops.end()) {
		ops.erase(it);
		op->assignToBatch(nullptr);
	}
}

void Batch::removeAllOps() {
	for (size_t o = 0; o < size(); ++o) {
		ops[o]->assignToBatch(nullptr);
	}
	ops.clear();
}

double Batch::getTWT() const {
	double twt = 0;
	for (size_t op = 0; op < size(); ++op) {
		twt += ops[op]->getTWT();
	}
	return twt;
}