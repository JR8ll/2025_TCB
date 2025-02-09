#include <iostream>
#include <string>

#include "Workcenter.h"
#include "Operation.h"

using namespace std;

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

int Workcenter::getId() const { return id; }
 
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
		// TODO

	}
	
	// TODO
}

double Workcenter::getTWT() const {
	double twt = 0;
	for (size_t m = 0; m < size(); ++m) {
		twt += machines[m]->getTWT();
	}
	return twt;
}
