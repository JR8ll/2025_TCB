#pragma once
#include <chrono>
#include <memory>
#include <vector>
#include "Functions.h"	// Sched_params

class Schedule;

class GaDecoderJobListSched {	// Random-key, job-based representation, List scheduling-based decoding
private:
	Schedule* masterSched;
	Sched_params* schedParams;
	GA_params* gaParams;
	std::chrono::time_point<std::chrono::high_resolution_clock> finishBy;

public:
	GaDecoderJobListSched(Schedule* schedule, Sched_params* schedParameters, GA_params* gaParameters);
	GaDecoderJobListSched(Schedule* schedule, Sched_params* schedParameters, GA_params* gaParameters, std::chrono::time_point<std::chrono::high_resolution_clock> finishBy);
	~GaDecoderJobListSched();
	double decode(const std::vector<double>& chr) const;
	double decodeWithLocalSearch(const std::vector<double>& chr) const;
	double applyNonPersitentLocalSearch(Schedule* sched) const;			// non-inheritable local search (not encodable)
	void formSchedule(const std::vector<double>& chr);
	Schedule* returnSchedule();
};

