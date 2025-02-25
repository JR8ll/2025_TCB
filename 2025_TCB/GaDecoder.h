#pragma once

#include<vector>
#include "Functions.h"	// Sched_params

class Schedule;

class GaDecoderJobLs {	// Random-key, job-based representation, List scheduling-based decoding
private:
	Schedule* masterSched;
	Sched_params* schedParams;

public:
	GaDecoderJobLs(Schedule* schedule, Sched_params* schedParameters);
	~GaDecoderJobLs();
	double decode(const std::vector<double>& chr) const;
	void formSchedule(const std::vector<double>& chr);
};

