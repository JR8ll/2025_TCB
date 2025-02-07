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

	Batch& operator[] (size_t idx);
	Batch& operator[] (size_t idx) const;

	size_t size() const;

	std::unique_ptr<Machine> clone(Workcenter* newWc) const;

	Workcenter* getWorkcenter();
	const std::vector<pBat>& getBatches() const;

	void addBatch(pBat batch);

};
