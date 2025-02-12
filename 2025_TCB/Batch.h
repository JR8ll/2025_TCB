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
	
	int f;			// product (family)
	int cap;		

public:
	Batch();
	Batch(int cap);
	~Batch();

	friend std::ostream& operator<<(std::ostream& os, const Batch& batch);

	std::unique_ptr<Batch> clone() const;	// shallow (no content)

	Operation& operator[](size_t idx);
	Operation& operator[](size_t idx) const;

	size_t getIdx() const;	// index in its machine
	size_t size() const;
	size_t findOp(const Operation* op) const;	// returns index of operation 
	bool isEmpty() const;

	double getStart() const;
	double getC() const;
	double getP() const;

	int getF() const;
	int getCap() const;
	int getAvailableCap() const;

	const std::vector<Operation*>& getOps() const;
	Machine* getMachine() const;

	void setStart(double newStart, bool checkValidity = true);
	void setC(double newC, bool checkvalidity = true);
	void setCap(int newCap);

	void assignToMachine(Machine* processor);

	bool addOp(Operation* op);
	void removeOp(Operation* op);
	void removeAllOps();

	void updateWaitingTimes();

	double getTWT() const;	// total weighted tardiness
};
