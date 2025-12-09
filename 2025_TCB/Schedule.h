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
	void schedOp(Operation* op, double pWait, double inflation, bool batchinStageInflationOnly = true, bool opsWithoutTcInflationOnly = true);

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
	void lSchedFirstJobInflated(double pWait, double inflation, bool batchinStageInflationOnly = true, bool opsWithoutTcInflationOnly = true);	// insert every operation (if it constitutes a new batch) delayed by (p * inflation)  
	void lSchedJobs(double pWait = 0.0);																				// List scheduling of jobs in member "jobs" in given order, pWait = accepted waiting time (ratio of processing time) if op can be added to exising batch
	void lSchedJobs(std::vector<double> pWaitVec = { 0.0 });
	void lSchedJobsInflated(double pWait, double inflation, bool batchinStageInflationOnly = true, bool opsWithoutTcInflationOnly = true);
	void lSchedJobsWithSorting(prioRule<pJob> rule, double pWait = 0.0);												// non-parameter sorting (EDD, SPT, ...)
	void lSchedJobsWithSorting(prioRule<pJob> rule, Sched_params& sched_params);
	void lSchedJobsWithSorting(prioRuleKappa<pJob> rule, double kappa, double pWait = 0.0);							// Dynamic ATC-like sorting with parameters t and kappa
	double lSchedJobsWithSorting(prioRuleKappa<pJob> rule, const std::vector<double>& kappaGrid, double pWait = 0.0, objectiveFunction objectiveFunction = &getObjectiveTWT);	// Dynamic ATC-like sorting with the best kappa from a grid, returns best kappa
	double lSchedJobsWithSorting(prioRuleKappa<pJob> rule, Sched_params& sched_params, objectiveFunction objectiveFunction = &getObjectiveTWT);
	void lSchedJobsWithRandomKeySorting(prioRuleKeySet<pJob> rule, const std::vector<double>& keys, double pWait = 0.0);			// Sorting by given random keys
	void lSchedJobsWithRandomKeySorting(prioRuleKeySet<pJob> rule, const std::vector<double>& keys, Sched_params& sched_params);
	void lSchedGifflerThompson(prioRule<pJob> rule, double pWait = 0.0);

	void leftShiftBatches();

	// LOCAL SEARCH
	// LOCAL SEARCH OPERATION BASED
	void localSearchOpLeftShifting(prioRule<pJob> rule = &sortJobsByWaitingTimeDecr, double pWait = 0.0);		
	
	// LOCAL SEARCH JOB BASED
	void localSearchJobSwapping(prioRule<pJob> rule = &sortJobsByD, bool bestFit = true);		// if not bestFit => firstFit
	void localSearchJobLeftShifting(prioRule<pJob> rule = &sortJobsByD, bool bestFit = true);	// if not bestFit => firstFit
	
	double locSearchEvaluateJobSwap(size_t idxFirst, size_t idxSecond, bool& feasible);			// positive return value => improvement
	std::pair<double, double> locSearchEvaluateJobLeftShift(size_t idxFirst, std::vector<std::vector<std::pair<double, double>>>& options);	// positive return value => improvement (first = job left shift (last op), second = sum of intermediate ops´ left shift) JOB BASED
	std::pair<double, double> locSearchEvaluateBatchLeftShift(Batch* batch, double time, bool& possible);
	std::pair<double, double> locSearchEvaluateBatchRightShift(Batch* batch, double time, bool& possible);
	//std::pair<double, double> locSearchEvaluateOpsLeftShift(size_t idxJob, size_t idxStg, std::vector<std::vector<std::pair<double, double>>>& options);
	//double locSearchEvaluateOpsRightShift(size_t idxJob, size_t idxStg, double delay);
	double locSearchEvaluateOpConsolidation(size_t idxJob, size_t idxStg, bool& feasible);						// positive return value => improvement
	


	// LOCAL SEARCH JOB BASED EXECUTION
	bool locSearchSwapJobs(size_t idxFirst, size_t idxSecond);
	bool locSearchLeftShiftJob(size_t jobIdx, std::vector<std::vector<std::pair<double, double>>>& options);

	// UTILITY FUNCTIONS FOR LOCAL SEARCH
	std::vector<std::pair<double, double>> getLeftShiftOptions(Operation* op);		   // return vector: [stage][possible improvement option from/to] Example 1: <0.0, 8.0> can be left shifted from 0 to 8 time units; Example 2
	std::vector<std::pair<double, double>> getRightShiftOptions(Operation* op, double minDelay);
	std::vector<double> getLeeway(Job* job);										   // gets leeway between a job´s operations´ processing (to define min left shifts for predecessors)
	std::vector<std::vector<std::pair<size_t, double>>> getTcSlack(Job* job);   // (to define max left shifts constrained by maximal time lags)
	
	bool constrainLeftShiftOptionsFromOverlaps(std::vector<std::vector<std::pair<double, double>>>& options, std::vector<double>& leeway);	// returns true if a change was applied to the options
	bool constrainLeftShiftOptionsFromTimeConstraints(std::vector<std::vector<std::pair<double, double>>>& options, std::vector<std::vector<std::pair<size_t, double>>>& tcSlack);
	void executeLeftShiftOption(size_t jobIdx, size_t stgIdx, std::pair<double, double>& option);									// left shift of single operation

	bool isValid() const;

	double getTWT() const;											// total weighted tardiness
	double getMinMSP(size_t stgIdx) const;							// smallest makespan of the machines at stage (workcenter)

	void saveJson(std::string solver = "N/A");
	void saveJsonFactory(std::string solver = "N/A");				// format complying zui5_gantt viewer application


	// DEBUGGING TOOLS (MAY LEAD TO INFEASIBLE SOLUTIONS)
	void debugSetR(size_t scheduledJobIdx, double newR);			// change r of a scheduled job
	void debugAddMachine(size_t stgIdx);							// add an empty machine
};