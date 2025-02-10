#include <algorithm>
#include <iostream>
#include <string>

#include "Machine.h"
#include "Batch.h"
#include "Functions.h"
#include "Workcenter.h"

using namespace std;

Machine::Machine(int id, int cap, Workcenter* wc) : id(id), cap(cap), r(0), workcenter(wc) {}

ostream& operator<<(ostream& os, const Machine& machine) {
	os << "M" << to_string(machine.getId()) << ": ";
	for (size_t b = 0; b < machine.size(); ++b) {
		os << machine[b];
	}
	return os;
}

Batch& Machine::operator[] (size_t idx) { return *batches[idx]; }
Batch& Machine::operator[] (size_t idx) const { return *batches[idx]; }

pBat Machine::getBatch(size_t idx) {
	return move(batches[idx]);
}

size_t Machine::size() const { return batches.size(); }

size_t Machine::findBatch(const Batch* bat) const {
	for (size_t b = 0; b < size(); ++b) {
		if (batches[b].get() == bat) return b;
	}
	throw ExcSched("Machine::findBatch() Batch not found");
}

int Machine::getId() const { return id; }
size_t Machine::getIdx() const {
	return workcenter->findMachine(this);
}
int Machine::getCap() const { return cap; }

unique_ptr<Machine> Machine::clone(Workcenter* newWc) const {
	auto newMachine = make_unique<Machine>(id, cap, newWc);
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

double Machine::getEarliestSlot(double from, double duration) const {
	double slot = max(r, from);
	for (size_t b = 0; b < batches.size(); ++b) {
		if ((slot + duration) <= batches[b]->getStart()) {
			return slot;
		}
		slot = max(slot, batches[b]->getC());
	}
	return slot;
}

void Machine::addBatch(pBat batch, double start, bool checkValidity){

	// check capacity constraint
	if (cap != batch->getCap()) {
		batch->setCap(cap);
	}

	// check overlapping processing and insert in correct position
	double c = start + batch->getP();
	bool bInserted = false;
	for (auto it = batches.begin(); it != batches.end(); ++it) {
		if (c - TCB::precision <= it->get()->getStart()) {
			auto newBatch = batches.insert(it, move(batch));
			newBatch->get()->assignToMachine(this);
			newBatch->get()->setStart(start, checkValidity);	// automatically sets c	
			bInserted = true;
			break;
		}
		else if ((start + TCB::precision) < it->get()->getC()) {
			throw ExcSched("Machine::addBatch() infeasible due to overlap");
		}
	}

	if (!bInserted) {
		batches.push_back(move(batch));
		batches.back()->assignToMachine(this);
		batches.back()->setStart(start);
	}
}
void Machine::eraseNullptr(size_t batIdx) {
	if(batIdx >= batches.size()) throw out_of_range("Machine::eraseNullptr() out of range");
	if (batches[batIdx] == nullptr) {
		batches.erase(batches.begin() + batIdx);
	}
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

void Machine::removeAllBatches() {
	for (size_t b = 0; b < size(); ++b) {
		batches[b]->removeAllOps();
		batches[b]->assignToMachine(nullptr);
	}
	batches.clear();
}

void Machine::moveBatch(pBat batch, double newStart) {
	size_t fromIdx = batch->getIdx();
	size_t toIdx = 0;

	bool bGapFound = false;
	bool bToTheEnd = false;

	for (auto it = batches.begin(); it != batches.end(); ++it) {
		if ((newStart + batch->getP()) <= it->get()->getStart() || 
			(it->get() == batch.get() && newStart < it->get()->getStart() && newStart + batch->getP() > it->get()->getStart())) {
			if (it != batches.begin()) {
				if ((it - 1)->get()->getC() <= newStart + TCB::precision || (it - 1)->get() == batch.get()) {
					toIdx = it - batches.begin();
					bGapFound = true;
					break;
				}
			}
			else {
				toIdx = 0;
				bGapFound = true;
				break;
			}
		}
		else {
			if (it == (batches.end() - 1)) {
				if (it->get()->getC() <= newStart) {
					toIdx = it - batches.begin();
					bGapFound = true;
					bToTheEnd = true;
				}
			}
		}
	}

	if (!bGapFound) throw ExcSched("Machine::moveBatch(... newStart) no available time slot");

	if (fromIdx != toIdx) {
		pBat movingBatch = removeBatch(fromIdx);
		batches.erase(batches.begin() + fromIdx);
		if (fromIdx < toIdx && !bToTheEnd) {
			batches.insert(batches.begin() + toIdx - 1, move(movingBatch));
			batches[toIdx-1]->setStart(newStart);
		}
		else {
			batches.insert(batches.begin() + toIdx, move(movingBatch));
			batches[toIdx]->setStart(newStart);
		}
	} else {
		batches[toIdx]->setStart(newStart);
	}
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