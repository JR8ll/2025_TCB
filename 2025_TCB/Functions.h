#pragma once

#include<random>

#include "Job.h"
#include "Common_aliases.h"
#include "Log.h"

static const int ALG_ITERATEDMILP = 1;		// iterated MILP solving
static const int ALG_LISTSCHEDATC = 2;		// simple List scheduling approach
static const int ALG_BRKGALISTSCH = 3;		// biased random-key ga with a list scheduling decoder
static const int ALG_BRKGALS2MILP = 4;		// get best sequence from brkga, then iteratively apply ops from this sequence to MILP
static const int ALG_ILS = 5;

struct GA_params;
struct DECOMPMILP_params;

class Problem;

namespace TCB {
	extern Logger logger;
	extern Problem* prob;
	extern unsigned seed;
	extern double precision;
	extern std::mt19937 rng;
}

struct Sched_params {	// solver indipendent generic parameters
	double pWaitLow;	// pWait is the fraction of an operation´s processing time, that is accepted as its waiting time if the operation can rather be inserted into some existing batch
	double pWaitHigh;
	double pWaitStep;

	double kappaLow;			// grid of kappa values considered for GATC-sorting
	double kappaHigh;
	double kappaStep;
};

class ExcSched {
private:
	std::string message;
public:
	ExcSched(const std::string& s) : message(s) {}
	const std::string& getMessage() const { return message; }
};

void processCmd(int argc, char* argv[], int& iSolver, int& iTilimSeconds, bool& bConsole, Sched_params& schedParams, GA_params& gaParams, DECOMPMILP_params& decompParams);
void writeSolutions(Schedule* solution, int solverType, std::string solverName, std::string objectiveName, int prescribedTime, int usedTime, Sched_params* schedParams, GA_params* gaParams, DECOMPMILP_params* decompParams);


void sortJobsByC(std::vector<pJob>& jobs);									// by completion (to be called on scheduled jobs)
void sortJobsByStart(std::vector<pJob>& jobs);								// by start	(to be called on scheduled jobs)
void sortJobsByWaitingTimeDecr(std::vector<pJob>& jobs);					// by waiting time (to be called on scheduled jobs)
void sortJobsByD(std::vector<pJob>& jobs);									// by due date
void sortJobsByR(std::vector<pJob>& jobs);									// by release time
void sortJobsByGATC(std::vector<pJob>& jobs, double t, double kappa);		// by global ATC
void sortJobsByRK(std::vector<pJob>& jobs, const std::vector<double>& chr);	// by random keys 

bool compJobsByC(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b);
bool compJobsByStart(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b);
bool compJobsByWaitingTimeDecr(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b);
bool compJobsByD(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b);
bool compJobsByR(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b);

void shiftJobFromVecToVec(std::vector<pJob>& source, std::vector<pJob>& target, size_t sourceIdx);

double getAvgP(const std::vector<pJob>& unscheduledJobs);

void loadSchedParams(Sched_params& schedParams, std::string filename);
void loadGaParams(GA_params& gaParams, std::string filename);
void loadDecompParams(DECOMPMILP_params& decompParams, std::string filename);
Sched_params getDefaultParams();

double getObjectiveTWT(const Schedule* sched);

std::vector<double> getDoubleGrid(double low, double high, double step);

std::string extractFileName(const std::string& fullPath);
void replaceWindowsSpecialCharsWithUnderscore(std::string& input);

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