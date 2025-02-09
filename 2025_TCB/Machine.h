#pragma once

#include <memory>
#include <vector>

#include "Batch.h"

class Workcenter;

using pBat = std::unique_ptr<Batch>;

class Machine {
private:
	int id;
	std::vector<pBat> batches;
	Workcenter* workcenter;

public:
	Machine(int id, Workcenter* wc);

	friend std::ostream& operator<<(std::ostream& os, const Machine& machine);

	Batch& operator[] (size_t idx);
	Batch& operator[] (size_t idx) const;

	size_t size() const;

	int getId() const;

	std::unique_ptr<Machine> clone(Workcenter* newWc) const;

	Workcenter* getWorkcenter();
	const std::vector<pBat>& getBatches() const;

	void addBatch(pBat batch, double start);
	pBat removeBatch(size_t idx);

	void assignToWorkcenter(Workcenter* wc);

	double getTWT() const;	// total weighted tardiness
};
