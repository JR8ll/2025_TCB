#pragma once

#include <memory>
#include <vector>

#include "Batch.h"

class Workcenter;

using pBat = std::unique_ptr<Batch>;

class Machine {
private:
	int id;

	int cap;	// capacity
	double r;	// machine availability
	std::vector<pBat> batches;
	Workcenter* workcenter;

public:
	Machine(int id, int cap, Workcenter* wc);

	friend std::ostream& operator<<(std::ostream& os, const Machine& machine);

	Batch& operator[] (size_t idx);
	Batch& operator[] (size_t idx) const;

	pBat getBatch(size_t idx);

	size_t size() const;
	size_t findBatch(const Batch* bat) const;	// returns index of bat

	int getId() const;
	size_t getIdx() const;	// index in its workcenter
	int getCap() const;

	std::unique_ptr<Machine> clone(Workcenter* newWc) const;

	Workcenter* getWorkcenter();
	const std::vector<pBat>& getBatches() const;

	double getEarliestSlot(double from, double duration) const;

	void addBatch(pBat batch, double start, bool checkvalidity = true);
	void eraseNullptr(size_t batIdx);

	pBat removeBatch(size_t idx);
	void moveBatch(pBat batch, double newStart);

	void assignToWorkcenter(Workcenter* wc);

	double getTWT() const;	// total weighted tardiness
};
