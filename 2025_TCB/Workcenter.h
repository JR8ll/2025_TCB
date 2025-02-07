#pragma once

#include<memory>
#include<vector>

class Machine;
class Schedule;

using pMac = std::unique_ptr<Machine>;

class Workcenter {
private:
	int id;
	std::vector<pMac> machines;
	Schedule* schedule;

public:
	Workcenter(int id, Schedule* sched);

	Schedule* getSchedule() const;
	const std::vector<pMac> getMachines() const;

	void addMachine(pMac mac);

};
