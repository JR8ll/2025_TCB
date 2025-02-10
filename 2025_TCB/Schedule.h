#pragma once

#include<memory>
#include<vector>

#include "Job.h"
#include "Workcenter.h"

using pWc = std::unique_ptr<Workcenter>;
using pJob = std::unique_ptr<Job>;
using sharedJob = std::shared_ptr<Job>;
using sharedOp = std::shared_ptr<Operation>;

class Schedule {
private:
	std::vector<pWc> workcenters;
	std::vector<pJob> jobs;
	std::vector<pJob> scheduledJobs;

	std::vector<sharedOp> unscheduled;
	std::vector<sharedOp> scheduled;

public:
	Schedule();

	friend std::ostream& operator<<(std::ostream& os, const Schedule& schedule);

	Workcenter& operator[](size_t idx);			// access workcenter
	Workcenter& operator[](size_t idx) const;	// access workcenter

	Job& getJob(size_t idx);
	Job& getJob(size_t idx) const;

	std::unique_ptr<Schedule> clone() const;
	// Deep copy of Schedule, Workcenters and Jobs. 
	// Shallow copy of Machines (no batches are copied)
	// "Deep" copy is achieved by reconstruction of batch/op assignments -> _reconstruct(...)

	void _reconstruct(const Schedule* orig); 

	size_t size() const;	// number of workcenters
	size_t getN() const;	// number of jobs considered

	const std::vector<pWc>& getWorkcenters() const;
	void addWorkcenter(pWc wc);
	void addJob(pJob job);

	void schedOp(Operation* op, double pWait = 0.0);

	void lSchedJobs(double pWait = 0.0);	// List scheduling of jobs in member "jobs" in given order, pWait = accepted waiting time (ratio of processing time) if op can be added to exising batch

	double getTWT() const;	// total weighted tardiness

};