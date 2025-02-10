#pragma once

#include<vector>

class Schedule;

class GaDecoderJobLs {	// Random-key, job-based representation, List scheduling-based decoding
private:
	Schedule* masterSched;

public:
	GaDecoderJobLs(Schedule* schedule);
	~GaDecoderJobLs();
	double decode(const std::vector<double>& chr) const;
	void formSchedule(const std::vector<double>& chr);
};

