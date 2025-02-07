#pragma once

#include <vector>

class Operation;

class Batch {
private:
	std::vector<Operation*> ops;

public:
	const std::vector<Operation*>& getOps() const;

	void addOp(Operation* op);
	void removeOp(Operation* op);
};
