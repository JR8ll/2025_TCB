#pragma once

#include <memory>
#include <vector>

class Machine;
class Operation;

class Batch {
private:
	std::vector<Operation*> ops;
	Machine* machine;

	double start;	// start time
	double c;		// completion time

public:
	Batch();

	friend std::ostream& operator<<(std::ostream& os, const Batch& batch);

	std::unique_ptr<Batch> clone() const;	// shallow (no content)

	Operation& operator[](size_t idx);
	Operation& operator[](size_t idx) const;

	size_t size() const;
	bool isEmpty() const;

	double getStart() const;
	double getC() const;

	const std::vector<Operation*>& getOps() const;

	void assignToMachine(Machine* processor);

	void addOp(Operation* op);
	void removeOp(Operation* op);

	double getTWT() const;	// total weighted tardiness
};
