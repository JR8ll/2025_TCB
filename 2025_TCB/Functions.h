#pragma once

#include<random>

#include "Job.h"
#include "Log.h"

using pJob = std::unique_ptr<Job>;

class Problem;

namespace TCB {
	extern Logger logger;
	extern Problem* prob;
	extern unsigned seed;
	extern double precision;
	extern std::mt19937 rng;
}

class ExcSched {
private:
	std::string message;
public:
	ExcSched(const std::string& s) : message(s) {}
	const std::string& getMessage() const { return message; }
};

void sortJobsByD(std::vector<pJob>& jobs);
void sortJobsByR(std::vector<pJob>& jobs);
void sortJobsByGATC(std::vector<pJob>& jobs, double t, double kappa);

bool compJobsByD(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b);
bool compJobsByR(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b);

void shiftJobFromVecToVec(std::vector<pJob>& source, std::vector<pJob>& target, size_t sourceIdx);

double getAvgP(const std::vector<pJob>& jobs);

class CompJobsByGATC {
public:
	double avgP;
	double t;
	double kappa;
	CompJobsByGATC(double avgP, double t, double kappa) : avgP(avgP), t(t), kappa(kappa) {};
	bool operator() (const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b) const {
		double iA = a->getGATC(avgP, t, kappa);
		double iB = b->getGATC(avgP, t, kappa);
		if (iA == iB) {
			return a->getId() < b->getId();
		}
		return iA > iB;
	}
};