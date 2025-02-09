#pragma once

#include<memory>
#include<vector>

#include "Machine.h";

class Schedule;

using pMac = std::unique_ptr<Machine>;

class Workcenter {
private:
	int id;
	std::vector<pMac> machines;
	Schedule* schedule;

public:
	Workcenter(int id, Schedule* sched);

	friend std::ostream& operator<<(std::ostream& os, const Workcenter& workcenter);

	Machine& operator[](size_t idx);
	Machine& operator[](size_t idx) const;

	std::unique_ptr<Workcenter> clone(Schedule* newSchedule) const;

	size_t size() const;
	
	int getId() const;

	Schedule* getSchedule() const;
	const std::vector<pMac>& getMachines() const;

	void addMachine(pMac mac);

	void schedOp(Operation* op, double pWait);

	double getTWT() const;	// total weighted tardiness

};
