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

	int getId() const;
	int getStg() const;
		
	double getP() const;
	double getWait() const;

	void setWait(double wt);
	void setPred(Operation* pre);
	void setSucc(Operation* suc);
	

	Job* getJob() const;
	Batch* getBatch() const;

	void assignToBatch(Batch* batch);

};