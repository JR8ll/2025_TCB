#include "Workcenter.h"

using namespace std;

Workcenter::Workcenter(int id, Schedule* sched) : id(id), schedule(sched) {}

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
	machines.push_back(move(mac));
}
