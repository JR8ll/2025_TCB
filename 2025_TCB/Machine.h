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

	Workcenter* getWorkcenter();
	const std::vector<pBat>& getBatches() const;

	void addBatch(pBat batch);

};
