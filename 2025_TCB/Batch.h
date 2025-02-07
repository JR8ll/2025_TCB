#pragma once

#include <memory>
#include <vector>

class Operation;

class Batch {
private:
	std::vector<Operation*> ops;

public:
	std::unique_ptr<Batch> clone() const;	// shallow (no content)

	const std::vector<Operation*>& getOps() const;

	void addOp(Operation* op);
	void removeOp(Operation* op);
};
