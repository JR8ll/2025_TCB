#pragma once

#include<memory>
#include<vector>

#include "Common_aliases.h"
#include "Job.h"
#include "Functions.h"
#include "Workcenter.h"

class Problem;

class Schedule {
private:
	std::vector<pWc> workcenters;
	std::vector<pJob> unscheduledJobs;
	std::vector<pJob> scheduledJobs;

	std::vector<sharedOp> unscheduled;
	std::vector<sharedOp> scheduled;
	Problem* problem;

public:
	Schedule();

	friend std::ostream& operator<<(std::ostream& os, const Schedule& schedule);

	Workcenter& operator[](size_t idx);			// access workcenter
	Workcenter& operator[](size_t idx) const;	// access workcenter

	Job& getJob(size_t idx);
	Job& getJob(size_t idx) const;

	pJob get_pJob(size_t idx);

	std::unique_ptr<Schedule> clone() const;
	// Deep copy of Schedule, Workcenters and Jobs. 
	// Shallow copy of Machines (no batches are copied)
	// "Deep" copy is achieved by reconstruction of batch/op assignments -> _reconstruct(...)

	void _reconstruct(const Schedule* orig); 

	size_t size() const;	// number of workcenters
	size_t getN() const;	// number of jobs considered

	bool contains(Operation* op) const;				// searches workcenters/machines/batches/operations (not scheduledJobs vector)

	int getCapAtStageIdx(size_t stgIdx) const;		// capacity (assumption: parallel identical machines)

	const std::vector<int> getBatchingStages() const;
	const std::vector<int> getDiscreteStages() const;

	const std::vector<pWc>& getWorkcenters() const;
	void addWorkcenter(pWc wc);
	void addJob(pJob job);

	void schedOp(Operation* op, double pWait = 0.0);

	Problem* getProblem() const;
	void setProblemRef(Problem* prob);

	void reset();																										// clear all batches/machines and shift all jobs back to unscheduled
	void clearJobs();																									// clears unscheduled + scheduled jobs, operations and their references to products

	void sortUnscheduled(prioRule<pJob> rule);
	void sortUnscheduled(prioRuleKappa<pJob> rule, double kappa);
	void sortUnscheduled(prioRuleKeySet<pJob> rule, std::vector<double>& chr);
	void sortScheduled(prioRule<pJob> rule);
	void updateWaitingTimes();
	void mimicWaitingTimes(const Schedule* wtSchedule);		// set waiting times of operations in this schedule according to the waiting times in the wtSchedule
	
	void markAsScheduled(size_t jobIdx);
	void markAsScheduled(pJob scheduledJob);																			// adds a job to the set of scheduled jobs
	int getNumberOfScheduledJobs() const;
	const Job* getScheduledJob(size_t idx) const;

	Operation* findInScheduledJobs(Operation* remoteOp) const;				// find operation in scheduled jobs by id and stage
	Operation* findInUnscheduledJobs(Operation* remoteOp) const;			// find operation in unscheduled jobs by id and stage

	// LIST SCHEDULING
	void lSchedFirstJob(double pWait = 0.0);
	void lSchedJobs(double pWait = 0.0);																				// List scheduling of jobs in member "jobs" in given order, pWait = accepted waiting time (ratio of processing time) if op can be added to exising batch
	void lSchedJobs(std::vector<double> pWaitVec = { 0.0 });
	void lSchedJobsWithSorting(prioRule<pJob> rule, double pWait = 0.0);												// non-parameter sorting (EDD, SPT, ...)
	void lSchedJobsWithSorting(prioRule<pJob> rule, Sched_params& sched_params);
	void lSchedJobsWithSorting(prioRuleKappa<pJob> rule, double kappa, double pWait = 0.0);							// Dynamic ATC-like sorting with parameters t and kappa
	double lSchedJobsWithSorting(prioRuleKappa<pJob> rule, const std::vector<double>& kappaGrid, double pWait = 0.0, objectiveFunction objectiveFunction = &getObjectiveTWT);	// Dynamic ATC-like sorting with the best kappa from a grid, returns best kappa
	double lSchedJobsWithSorting(prioRuleKappa<pJob> rule, Sched_params& sched_params, objectiveFunction objectiveFunction = &getObjectiveTWT);
	void lSchedJobsWithRandomKeySorting(prioRuleKeySet<pJob> rule, const std::vector<double>& keys, double pWait = 0.0);			// Sorting by given random keys
	void lSchedJobsWithRandomKeySorting(prioRuleKeySet<pJob> rule, const std::vector<double>& keys, Sched_params& sched_params);
	void lSchedGifflerThompson(prioRule<pJob> rule, double pWait = 0.0);

	void localSearchLeftShifting(prioRule<pJob> rule = &sortJobsByWaitingTimeDecr, double pWait = 0.0);

	bool isValid() const;

	double getTWT() const;											// total weighted tardiness
	double getMinMSP(size_t stgIdx) const;							// smallest makespan of the machines at stage (workcenter)

	void saveJson(std::string solver = "N/A");
	void saveJsonFactory(std::string solver = "N/A");				// format complying zui5_gantt viewer application

};