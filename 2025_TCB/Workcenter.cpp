#include "Workcenter.h"

using namespace std;

Workcenter::Workcenter(int id, Schedule* sched) : id(id), schedule(sched) {}

Schedule* Workcenter::getSchedule() const {
	return schedule;
}
const std::vector<pMac>& Workcenter::getMachines() const {
	return machines;
}

void Workcenter::addMachine(pMac mac) {
	machines.push_back(move(mac));
}
