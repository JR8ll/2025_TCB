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

	bool operator==(const Batch& other) const;	// true, if content is identical
	
	Operation& operator[](size_t idx);
	Operation& operator[](size_t idx) const;


	size_t getIdx() const;	// index in its machine
	size_t size() const;
	size_t findOp(const Operation* op) const;	// returns index of operation 
	bool isEmpty() const;
	bool contains(const Operation* op) const;			

	double getStart() const;
	double getC() const;
	double getP() const;
	double getR() const;	// largest availability of operations at this stage
	double getRconsideringRawP() const;
	double getW() const;	// sum of job weights

	int getF() const;
	int getCap() const;
	int getAvailableCap() const;

	const std::vector<Operation*>& getOps() const;
	Machine* getMachine() const;

	void setStart(double newStart, bool checkValidity = true);
	void setC(double newC, bool checkvalidity = true);
	void setCap(int newCap);
	bool rightShiftFixedBatchFormationFixedMachineAssigment(double time);	// this move will keep all batch formations + machine assignments and enforce only backwards tc by right shifting predecessors in turn

	void assignToMachine(Machine* processor);

	bool addOp(Operation* op);
	void removeOp(Operation* op);
	void removeOp(size_t idx);
	void removeAllOps();

	void updateWaitingTimes();

	double getTWT() const;	// total weighted tardiness
};
