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
	Product* product;

	Operation* pred;	// route predecessor
	Operation* succ;	// route successor

public:
	Operation(Job* j, Product* f, int stg); 

	int getId() const;
	int getStg() const;
	
	double getP() const;
	double getWait() const;
	

	Job* getJob() const;
	Batch* getBatch() const;

	void assignToBatch(Batch* batch);

};