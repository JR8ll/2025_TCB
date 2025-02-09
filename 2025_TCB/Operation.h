#pragma once

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

	friend std::ostream& operator<<(std::ostream& os, const Operation& operation);

	int getId() const;
	int getStg() const;

	int getS() const;
	double getD() const;	// external due date (of the job)
	double getP() const;
	double getW() const;
		
	
	double getWait() const;

	double getGATC(double avgP, double t, double kappa) const;

	Operation* getPred() const;
	Operation* getSucc() const;

	void setWait(double wt);
	void setPred(Operation* pre);
	void setSucc(Operation* suc);

	Job* getJob() const;
	Batch* getBatch() const;

	void assignToBatch(Batch* batch);

	double getTWT() const;

};