#include <iostream>
#include <string>

#include "Machine.h"
#include "Batch.h"
#include "Workcenter.h"

using namespace std;

Machine::Machine(int id, Workcenter* wc) : id(id), workcenter(wc) {}

ostream& operator<<(ostream& os, const Machine& machine) {
	os << "M" << to_string(machine.getId()) << ": ";
	for (size_t b = 0; b < machine.size(); ++b) {
		os << machine[b];
	}
	return os;
}

Batch& Machine::operator[] (size_t idx) { return *batches[idx]; }
Batch& Machine::operator[] (size_t idx) const { return *batches[idx]; }

size_t Machine::size() const { return batches.size(); }

int Machine::getId() const { return id; }

unique_ptr<Machine> Machine::clone(Workcenter* newWc) const {
	auto newMachine = make_unique<Machine>(id, newWc);
	for (const auto& bat : batches) {
		newMachine->addBatch(bat->clone(), bat->getStart());
	}
	return newMachine;
}

Workcenter* Machine::getWorkcenter() {
	return workcenter;
}
const std::vector<pBat>& Machine::getBatches() const {
	return batches;
}

void Machine::addBatch(pBat batch, double start){
	batches.push_back(move(batch));
	batch->assignToMachine(this);
}
pBat Machine::removeBatch(size_t idx) {
	if (idx >= batches.size()) {
		throw out_of_range("Machine::removeBatch out of range");
	}
	pBat removedBatch = move(batches[idx]);
	batches.erase(batches.begin() + idx);
	removedBatch->assignToMachine(nullptr);
	return removedBatch;
}

void Machine::assignToWorkcenter(Workcenter* wc) {
	workcenter = wc;
}

double Machine::getTWT() const {
	double twt = 0;
	for (size_t b = 0; b < size(); ++b) {
		twt += batches[b]->getTWT();
	}
	return twt;
}