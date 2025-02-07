#pragma once

class Batch;
class Job;

class Operation {
private:
	Job* job;
	Batch* batch;

public:
	Operation(Job* j); 

	Job* getJob() const;
	Batch* getBatch() const;

	void assignToBatch(Batch* batch);

};