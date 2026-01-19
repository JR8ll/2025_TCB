#pragma once
#include <memory>
#include <vector>
#include "Functions.h"	// Sched_params

class Schedule;

class GaDecoderJobListSched {	// Random-key, job-based representation, List scheduling-based decoding
private:
	Schedule* masterSched;
	Sched_params* schedParams;
	GA_params* gaParams;

public:
	GaDecoderJobListSched(Schedule* schedule, Sched_params* schedParameters, GA_params* gaParameters);
	~GaDecoderJobListSched();
	double decode(const std::vector<double>& chr) const;
	void applyPersistentLocalSearch(std::unique_ptr<Schedule> sched, std::vector<double>& chr);					// it´s hard to tell if a swap or left shift in the chromosome is equivalent to a swap in the schedule because of integrated batching decisions
	double applyNonPersitentLocalSearch(std::unique_ptr<Schedule> sched) const;									// non-inheritabel local search (not encodable)
	void formSchedule(const std::vector<double>& chr);
};

