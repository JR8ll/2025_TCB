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

	int size() const;
	size_t findMachine(const Machine* mac) const;	// returns index of mac
	
	int getId() const;
	int getCap() const;		// assumption: parallel identical machines

	Schedule* getSchedule() const;
	const std::vector<pMac>& getMachines() const;

	void addMachine(pMac mac);

	void schedOp(Operation* op, double pWait = 0.0);
	void schedOp(Operation* op, double pWait, double inflation, bool batchinStageInflationOnly = true, bool opsWithoutTcInflationOnly = true);
	void schedOpDelayed(Operation* op, double startingAt);
	void ensureValidity(Operation* op);
	void ensureValidityFixedBatchFormation(Operation* op);						// keeping batch formation unchanged
	bool leftShift(size_t mIdx, size_t bIdx, size_t jIdx, double pWait = 0.0);	// true if left-shift was performed
	bool leftShift(size_t mIdx, size_t bIdx, double pWait = 0.0);				// true if left-shift was performed
	void rightShiftOp(size_t mIdx, size_t bIdx, size_t jIdx, double from, double pWait = 0.0);	// indices identify op to be right-shifted, from is the new earliest starting time
	void rightShiftBatch(size_t mIdx, size_t bIdx, double from, bool pushingSuccessors = true, bool checkValidity = true);	// if pushing successors the batch may always be moved on its machine, if necessary pushing right its successors
	void findBestStart(Operation* op, bool& newBatch, size_t& bestMacIdx, size_t& bestBatIdx, double& bestStart, double pWait = 0.0);
	void findBestStart(Operation* op, bool& newBatch, size_t& bestMacIdx, size_t& bestBatIdx, double& bestStart, double pWait, double inflation);
	void findBestStartNotBefore(Operation* op, bool& newBatch, size_t& bestMacIdx, size_t& bestBatIdx, double& bestStart, double notBefore);
	void findBestStartNotBefore(Batch* batch, size_t& bestMacIdx, double& tempStart, double notBefore);
	double findLatestAvailableTimeSlotBefore(double latest, double duration);
	bool locateOp(Operation* op, size_t& mIdx, size_t& batIdx, size_t& jIdx);	// true if found
	bool swapOps(size_t mIdx1, size_t bIdx1, size_t jIdx1, size_t mIdx2, size_t bIdx2, size_t jIdx2);
	bool moveOpDisregardingTc(Operation* op, double newStart, bool& intoBatch);	// disregarding checks on tc violations

	bool localSearchLeftShift(double pWait = 0.0);	// true if at least one operation was shifted left

	void moveBatch(Batch* batch, size_t tgtMac, double newStart, bool checkValidity = true);

	void updateWaitingTimes();

	double getTWT() const;		// total weighted tardiness
	double getMinMSP() const;	// smallest makespan (completion of last batch)

	bool debugNoTwoBatchesStartingSimultaneously();
	
};
