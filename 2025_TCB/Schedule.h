#pragma once

#include<memory>
#include<vector>

#include "Job.h"
#include "Workcenter.h"

using pWc = std::unique_ptr<Workcenter>;
using pJob = std::unique_ptr<Job>;
using sharedJob = std::shared_ptr<Job>;
using sharedOp = std::shared_ptr<Operation>;

template<typename T>
using prioRule = void(*)(std::vector<T>& vec);

template<typename T>
using prioRuleKappaT = void(*)(std::vector<T>& vec, double t, double kappa);

template<typename T>
using prioRuleKeySet = void(*)(std::vector<T>& vec, const std::vector<double>& keys);

class Schedule {
private:
	std::vector<pWc> workcenters;
	std::vector<pJob> unscheduledJobs;
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

	void reset();																										// clear all batches/machines and shift all jobs back to unscheduled
	void clearJobs();																									// clears unscheduled + scheduled jobs, operations and their references to products

	void lSchedFirstJob(double pWait = 0.0);
	void lSchedJobs(double pWait = 0.0);																				// List scheduling of jobs in member "jobs" in given order, pWait = accepted waiting time (ratio of processing time) if op can be added to exising batch
	void lSchedJobsWithSorting(prioRule<pJob> rule, double pWait = 0.0);												// non-parameter sorting (EDD, SPT, ...)
	void lSchedJobsWithSorting(prioRuleKappaT<pJob> rule, double kappa, double pWait = 0.0);							// Dynamic ATC-like sorting with parameters t and kappa
	void lSchedJobsWithSorting(prioRuleKappaT<pJob> rule, const std::vector<double>& kappaGrid, double pWait = 0.0);	// like above with a grid search of kappa values
	void lSchedJobsWithRandomKeySorting(prioRuleKeySet<pJob> rule, const std::vector<double>& keys, double pWait = 0.0);			// Sorting by given random keys

	double getTWT() const;											// total weighted tardiness
	double getMinMSP(size_t stgIdx) const;							// smallest makespan of the machines at stage (workcenter)

};