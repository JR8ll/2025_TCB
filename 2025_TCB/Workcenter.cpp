#include <algorithm>
#include <iostream>
#include <string>

#include "Workcenter.h"
#include "Functions.h"
#include "Operation.h"

using namespace std;
using pBat = unique_ptr<Batch>;

Workcenter::Workcenter(int id, Schedule* sched) : id(id), schedule(sched) {}

ostream& operator<<(ostream& os, const Workcenter& workcenter) {
	os << "WC" << to_string(workcenter.getId()) << ": " << endl;
	for (size_t m = 0; m < workcenter.size(); ++m) {
		os << workcenter[m] << endl;
	}
	return os;
}

Machine& Workcenter::operator[](size_t idx) { return *machines[idx]; }
Machine& Workcenter::operator[](size_t idx) const { return *machines[idx]; }

unique_ptr<Workcenter> Workcenter::clone(Schedule* newSchedule) const {
	auto newWc = make_unique<Workcenter>(id, newSchedule);
	for (const auto& mac : machines) {
		newWc->addMachine(mac->clone(newWc.get()));
	}
	return newWc;
}

size_t Workcenter::size() const { return machines.size(); }

size_t Workcenter::findMachine(const Machine* mac) const {
	for (size_t m = 0; m < size(); ++m) {
		if (machines[m].get() == mac) return m;
	}
	throw ExcSched("Workcenter::findMachine() Machine not found");
}

int Workcenter::getId() const { return id; }

int Workcenter::getCap() const {
	if (machines.empty()) return 0;
	return machines[0]->getCap();
}
 
Schedule* Workcenter::getSchedule() const {
	return schedule;
}
const std::vector<pMac>& Workcenter::getMachines() const {
	return machines;
}

void Workcenter::addMachine(pMac mac) {
	mac->assignToWorkcenter(this);
	machines.push_back(move(mac));
}


void Workcenter::schedOp(Operation* op, double pWait) {
	size_t bestMacIdx = -1;
	size_t bestBatIdx = -1;
	bool bNewBatch = true;
	double idealStart = op->getEarliestStart();
	double tempStart = numeric_limits<double>::max();

	findBestStart(op, bNewBatch, bestMacIdx, bestBatIdx, tempStart, pWait);

	// [JR-2025-02-25] findBestStart not yet tested
	// find best batch/time slot for operation TODO: call findBestStart instead
	//for (size_t m = 0; m < machines.size(); ++m) {
	//	Machine* mac = machines[m].get();
	//	// consider existing batches
	//	for (size_t b = 0; b < mac->size(); ++b) {
	//		Batch* bat = &(*mac)[b];
	//		if (bat->getStart() > tempStart) break; // earlier option already found
	//		if (bat->getStart() >= idealStart && bat->getF() == op->getF() && bat->getAvailableCap() >= op->getS()) {
	//			if (bat->getStart() < tempStart) {
	//				tempStart = bat->getStart();
	//				bestMacIdx = m;
	//				bestBatIdx = b;
	//				bNewBatch = false;
	//			}
	//			break;  // later batches on this machines are not considered
	//		}
	//	}

	//	// consider formation of a new batch
	//	double earliestSlot = mac->getEarliestSlot(idealStart, op->getP());
	//	if (earliestSlot + (op->getP() * pWait) < tempStart) {
	//		if (tempStart >= idealStart) {
	//			tempStart = earliestSlot;
	//			bestMacIdx = m;
	//			bNewBatch = true;
	//		}
	//	}
	//}

	// actually schedule operation
	Machine* bestMac = machines[bestMacIdx].get();
	if (!bNewBatch) {
		Batch* bestBat = &(*bestMac)[bestBatIdx];
		if(!bestBat->addOp(op)) throw ExcSched("Workcenter::schedOp -> op could not be added");
	} else {
		pBat newBatch = make_unique<Batch>(bestMac->getCap());
		newBatch->addOp(op);
		if(!bestMac->addBatch(move(newBatch), tempStart)) throw ExcSched("Workcenter::schedOp -> batch could not be added");
	}
	ensureValidity(op);
}

void Workcenter::ensureValidity(Operation* op) {
	bool bValid = false;
	bool bOverlaps = true;
	bool bTcViolations = true;
	while (!bValid) {
		bOverlaps = op->repairOverlaps();
		bTcViolations = op->repairTimeConstraints();
		bValid = !bOverlaps && !bTcViolations;
	}
}

bool Workcenter::leftShift(size_t mIdx, size_t bIdx, size_t jIdx, double pWait) {
	Operation* op = &(*machines[mIdx])[bIdx][jIdx];
	Batch* bat = op->getBatch();
	double currentStart = op->getStart();
	bool bNewBatch = false;
	size_t bestMacIdx = 0;
	size_t bestBatIdx = 0;
	double betterStart = currentStart;
	findBestStart(op, bNewBatch, bestMacIdx, bestBatIdx, betterStart, pWait);
	if (betterStart < currentStart) {
		// actually shift operation
		size_t currentBatIdx = op->getBatch()->getIdx();
		size_t currentMacIdx = op->getBatch()->getMachine()->getIdx();
		if (!bNewBatch) {
			auto tempOp = op;
			(*machines[currentMacIdx])[currentBatIdx].removeOp(op);
			(*machines[bestMacIdx])[bestBatIdx].addOp(tempOp);
			if ((*machines[currentMacIdx])[currentBatIdx].isEmpty()) {
				machines[currentMacIdx]->removeBatch(currentBatIdx);
			}
		}
		else {
			if ((*machines[currentMacIdx])[currentBatIdx].size() == 1) {
				moveBatch(bat, bestMacIdx, betterStart);

			}
			else {
				pBat newBatch = make_unique<Batch>(machines[bestMacIdx]->getCap());
				auto tempOp = op;
				bat->removeOp(op);
				newBatch->addOp(tempOp);
				machines[bestMacIdx]->addBatch(move(newBatch), betterStart);
			}
		}
		return true;
	}
	return false;
}
void Workcenter::rightShift(size_t mIdx, size_t bIdx, size_t jIdx, double from, double pWait) {
	// TODO outsource search for best scheduling option from schedOp and rightShift to new method + ensure validity
	Machine* mac = machines[mIdx].get();
	Batch* bat = &(*mac)[bIdx];
	Operation* op = &(*bat)[jIdx];

	double earliestStart = op->getEarliestStart();
	bool bOnlyOperation = bat->size() <= 1;
	double newStart = max(from, earliestStart);

	bool bNewBatch = true;
	size_t bestMacIdx = 0;
	size_t bestBatIdx = 0;
	double bestStart = numeric_limits<double>::max();

	findBestStart(op, bNewBatch, bestMacIdx, bestBatIdx, bestStart, pWait);

	//// [JR-2025-Feb-25] findBestStart not yet tested
	//// find best batch/time slot for operation TODO: call findBestStart instead
	//double tempStart = numeric_limits<double>::max();
	//for (size_t m = 0; m < size(); ++m) {
	//	Machine* tempMac = machines[m].get();
	//	for (size_t b = 0; b < tempMac->size(); ++b) {
	//		Batch* tempBat = &(*tempMac)[b];
	//		if (tempBat->getStart() > tempStart) {
	//			break;	// earlier option already found
	//		}
	//		// consider existing batches
	//		if (tempBat->getStart() >= newStart && tempBat->getF() == op->getF() && tempBat->getAvailableCap() >= op->getS()) {
	//			if (tempBat->getStart() < tempStart) {
	//				tempStart = bat->getStart();
	//				bestMacIdx = m;
	//				bestBatIdx = b;
	//				bNewBatch = false;
	//			}
	//			break;
	//		}
	//	}

	//	// consider formation of a new batch
	//	double earliestSlot = tempMac->getEarliestSlot(newStart, op->getP());
	//	if (earliestSlot + (op->getP() * pWait) < tempStart) {
	//		if (tempStart >= newStart) {
	//			tempStart = earliestSlot;
	//			bestMacIdx = m;
	//			bNewBatch = true;
	//		}
	//	}
	//}

	// actually shift operation
	if (!bNewBatch) {
		auto tempOp = op;
		bat->removeOp(op);
		(*machines[bestMacIdx])[bestBatIdx].addOp(tempOp);	// [JR-2025-Mar-03] replaced op with tempOp
		if (bOnlyOperation) {
			mac->removeBatch(bIdx);
		}
	}
	else {
		if (bOnlyOperation) {
			moveBatch(bat, bestMacIdx, bestStart);

		}
		else {
			pBat newBatch = make_unique<Batch>(machines[bestMacIdx]->getCap());
			auto tempOp = op;
			bat->removeOp(op);
			newBatch->addOp(tempOp);
			machines[bestMacIdx]->addBatch(move(newBatch), bestStart);
		}
	}

	ensureValidity(op);
}
void Workcenter::findBestStart(Operation* op, bool& bNewBatch, size_t& bestMacIdx, size_t& bestBatIdx, double& tempStart, double pWait) {
	double idealStart = op->getEarliestStart();
	tempStart = numeric_limits<double>::max();

	// find best batch/time slot for operation
	for (size_t m = 0; m < machines.size(); ++m) {
		Machine* mac = machines[m].get();
		// consider existing batches
		for (size_t b = 0; b < mac->size(); ++b) {
			Batch* bat = &(*mac)[b];
			if (bat->getStart() > tempStart) break; // earlier option already found
			if (bat->getStart() >= idealStart && bat->getF() == op->getF() && bat->getAvailableCap() >= op->getS()) {
				if (bat->getStart() < tempStart) {
					tempStart = bat->getStart();
					bestMacIdx = m;
					bestBatIdx = b;
					bNewBatch = false;
				}
				break;  // later batches on this machines are not considered
			}
		}

		// consider formation of a new batch
		double earliestSlot = mac->getEarliestSlot(idealStart, op->getP());
		if (earliestSlot + (op->getP() * pWait) < tempStart) {
			if (tempStart >= idealStart) {
				tempStart = earliestSlot;
				bestMacIdx = m;
				bNewBatch = true;
			}
		}
	}
}
bool Workcenter::locateOp(Operation* op, size_t& mIdx, size_t& batIdx, size_t& jIdx) {
	for (size_t m = 0; m < machines.size(); ++m) {
		for (size_t b = 0; b < (*machines[m]).size(); ++b) {
			for (size_t j = 0; j < (*machines[m])[b].size(); ++j) {
				Operation* candidate = &(*machines[m])[b][j];
				if (op->getId() == candidate->getId() && op->getStg() == candidate->getStg()) {
					mIdx = m;
					batIdx = b;
					jIdx = j;
					return true;
				}
			}
		}
	}
	return false;
}

bool Workcenter::localSearchLeftShift(double pWait) {
	// t´would be nice if we considered batches by order of non-decreasing start times
	size_t mIdx = 0;	
	size_t bIdx = 0;
	//double tempStart = DBL_MAX;
	//for (size_t m = 0; m < machines.size(); ++m) {
	//	if (machines[m]->size() > bIdx) {
	//		if ((*machines[m])[bIdx].getStart() < tempStart) {
	//			tempStart = (*machines[m])[bIdx].getStart();
	//			mIdx = m;
	//		}
	//	}
	//	bool bShifted = false;
	//	for (size_t j = 0; j < (*machines[mIdx])[bIdx].size(); ++j) {
	//		if (leftShift(mIdx, bIdx, j, pWait)) {
	//			bShifted = true;
	//		}
	//	}
	//}
	return false;
}

void Workcenter::moveBatch(Batch* batch, size_t tgtMac, double newStart) {
	size_t batIdx = batch->getIdx();
	size_t macIdx = batch->getMachine()->getIdx();
	if (tgtMac == macIdx) {
		// same machine
		return machines[tgtMac]->moveBatch(move(batch), newStart);
	}
	// different machine
	machines[tgtMac]->addBatch(move(machines[macIdx]->getBatch(batIdx)), newStart);
	machines[macIdx]->eraseNullptr(batIdx);
}

void Workcenter::updateWaitingTimes() {
	for (size_t m = 0; m < size(); ++m) {
		machines[m]->updateWaitingTimes();
	}
}

double Workcenter::getTWT() const {
	double twt = 0;
	for (size_t m = 0; m < size(); ++m) {
		twt += machines[m]->getTWT();
	}
	return twt;
}
double Workcenter::getMinMSP() const {
	double minMSP = numeric_limits<double>::max();
	for (size_t m = 0; m < size(); ++m) {
		double tempMSP = machines[m]->getMSP();
		if (tempMSP < minMSP) {
			minMSP = tempMSP;
		}
	}
	return minMSP;
}
