#pragma once

#include<memory>
#include<vector>

#include "Common_aliases.h"
#include "Machine.h";

class Schedule;

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
	size_t findMachine(const Machine* mac) const;	// returns index of mac
	
	int getId() const;
	int getCap() const;		// assumption: parallel identical machines

	Schedule* getSchedule() const;
	const std::vector<pMac>& getMachines() const;

	void addMachine(pMac mac);

	void schedOp(Operation* op, double pWait = 0.0);
	void ensureValidity(Operation* op);
	void rightShift(size_t mIdx, size_t bIdx, size_t jIdx, double from, double pWait = 0.0);	// indices identify op to be right-shifted, from is the new earliest starting time
	void findBestStart(Operation* op, bool& newBatch, size_t& bestMacIdx, size_t& bestBatIdx, double& bestStart, double pWait = 0.0);

	void moveBatch(Batch* batch, size_t tgtMac, double newStart);

	void updateWaitingTimes();

	double getTWT() const;		// total weighted tardiness
	double getMinMSP() const;	// smallest makespan (completion of last batch)

};
