#pragma once

#include <memory>
#include <vector>

class Operation;

class Batch {
private:
	std::vector<Operation*> ops;

public:
	std::unique_ptr<Batch> clone() const;	// shallow (no content)

	Operation& operator[](size_t idx);
	Operation& operator[](size_t idx) const;

	size_t size() const;

	const std::vector<Operation*>& getOps() const;

	void addOp(Operation* op);
	void removeOp(Operation* op);
};
