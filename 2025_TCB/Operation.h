#pragma once

#include <vector>

class Batch;
class Job;
class Product;

class Operation {
private:
	int id;		
	int stg;			// stage (workcenter)
	double wait;		// waiting time

	Job* job;
	Batch* batch;

	Operation* pred;	// route predecessor
	Operation* succ;	// route successor

public:
	Operation(Job* j, int stg); 
	~Operation();

	friend std::ostream& operator<<(std::ostream& os, const Operation& operation);

	int getId() const;
	size_t getIdxInBatch() const;
	int getStg() const;
	int getWorkcenterId() const;	

	int getS() const;			// size
	int getF() const;			// product id (family)
	double getStart() const;	// start of processing(infinite if op was not yet schedule)
	double getC() const;		// completion (infinite if op was not yet scheduled)
	double getD() const;		// external due date (of the job)
	double getP() const;
	double getR() const;		// external release date (of the job)
	double getW() const;
		
	bool isScheduled() const;
	double getAvailability() const;		// internal release (after predecessors completion)
	double getEarliestStart() const;	// constrained by job release, c of earlier steps, tc of later steps
	double getWait() const;

	double getGATC(double avgP, double t, double kappa) const;

	Operation* getPred() const;
	Operation* getSucc() const;

	const std::vector<std::pair<int, double>>& getTcMaxFwd() const;
	const std::vector<std::pair<int, double>>& getTcMaxBwd() const;

	void setWait(double wt);
	void setPred(Operation* pre);
	void setSucc(Operation* suc);

	Job* getJob() const;
	Batch* getBatch() const;

	void assignToBatch(Batch* batch);
	bool repairOverlaps();													// true, if an overlap was found and repaired
	bool repairTimeConstraints();											// true, if a tc violation was found and repaired

	double getTWT() const;


	void resetLinks();
};