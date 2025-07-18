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
	/*for (const auto& bat : batches) {
		newMachine->addBatch(bat->clone(), bat->getStart());
	}*/
	return newMachine;
}

Workcenter* Machine::getWorkcenter() {
	return workcenter;
}
const std::vector<pBat>& Machine::getBatches() const {
	return batches;
}

double Machine::getEarliestSlot(double from, const Operation& op) const {
	double slot = max(r, from);
	for (size_t b = 0; b < batches.size(); ++b) {
		if ((slot + op.getP()) <= batches[b]->getStart()
			|| batches[b]->size() == 1 && (*batches[b])[0].getId() == op.getId()) {	// [JR-2025-Jul-15] case: overlapping with self
			return slot;
		}
		if (batches[b]->size() != 1 || (*batches[b])[0].getId() != op.getId()) {	// [JR-2025-Jul-15] case: overlapping with self
			slot = max(slot, batches[b]->getC());
		}	
	}
	return slot;
}

bool Machine::addBatch(pBat batch, double start, bool checkValidity){

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
			return false;
		}
	}

	if (!bInserted) {
		batches.push_back(move(batch));
		batches.back()->assignToMachine(this);
		batches.back()->setStart(start);
	}
	return true;
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

void Machine::moveBatch(Batch* batch, double newStart) {
	size_t fromIdx = batch->getIdx();
	size_t toIdx = 0;

	bool bGapFound = false;
	bool bToTheEnd = false;

	for (auto it = batches.begin(); it != batches.end(); ++it) {
		if ((newStart + batch->getP()) <= it->get()->getStart() || 
			(it->get() == batch && newStart < it->get()->getStart() && newStart + batch->getP() > it->get()->getStart())) {
			if (it != batches.begin()) {
				if ((it - 1)->get()->getC() <= newStart + TCB::precision || (it - 1)->get() == batch) {
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
				if (it->get()->getC() <= newStart ||
					it->get() == batch) {	// [JR-2025-Jul-15] case: overlapping with self
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
		if (fromIdx < toIdx && !bToTheEnd) {
			batches.insert(batches.begin() + toIdx - 1, move(movingBatch));
			batches[toIdx-1]->assignToMachine(this);
			batches[toIdx-1]->setStart(newStart);
		}
		else {
			batches.insert(batches.begin() + toIdx, move(movingBatch));
			batches.back()->assignToMachine(this);
			batches[toIdx]->setStart(newStart);
		}
	} else {
		batches[toIdx]->setStart(newStart);
	}
}

void Machine::assignToWorkcenter(Workcenter* wc) {
	workcenter = wc;
}

bool Machine::hasOverlaps() const {
	for (size_t b1 = 0; b1 < size(); ++b1) {
		double start1 = batches[b1]->getStart(); 
		double completion1 = batches[b1]->getC();
		for (size_t b2 = 0; b2 < size(); ++b2) {
			if (b1 != b2) {
				double start2 = batches[b2]->getStart();
				double completion2 = batches[b2]->getC();
				if (start1 + TCB::precision < start2 && completion1 - TCB::precision > start2) {
					return true;
				}
				if (start1 - TCB::precision > start2 && completion1 + TCB::precision < completion2) {
					return true;
				}
			}
		}
	}
	return false;
}

void Machine::updateWaitingTimes() {
	for (size_t b = 0; b < size(); ++b) {
		batches[b]->updateWaitingTimes();
	}
}

double Machine::getTWT() const {
	double twt = 0;
	for (size_t b = 0; b < size(); ++b) {
		twt += batches[b]->getTWT();
	}
	return twt;
}

double Machine::getMSP() const {
	if (batches.empty()) return r;
	return batches.back()->getC();
}
