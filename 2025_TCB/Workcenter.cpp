#include <iostream>
#include <string>

#include "Workcenter.h"

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

double Workcenter::getTWT() const {
	double twt = 0;
	for (size_t m = 0; m < size(); ++m) {
		twt += machines[m]->getTWT();
	}
	return twt;
}
