#include <algorithm>
#include <string>
#include <fstream>

#include "Functions.h"
#include "Solver_GA.h"
#include "Problem.h"

using namespace std;

void processCmd(int argc, char* argv[], int& iSolver, int& iTilimSeconds, bool& bConsole, GA_params& gaParams) {
	// argv[1] filename of problem instance to be solved
	if (argc > 1) {
		TCB::prob->loadFromDat(argv[1]);
	}
	else {
		TCB::prob->loadFromDat("debug_data.dat");
	}
	// argv[2] seed for pseudo random-number generator
	if (argc > 2) {
		TCB::seed = atoi(argv[2]);
		TCB::rng = mt19937(TCB::seed);
	}
	
	// argv[3] int describing the solving method to be used
	if (argc > 3) {
		iSolver = atoi(argv[3]);
	}
	
	// argv[4] time limit in seconds
	if (argc > 4) {
		iTilimSeconds = atoi(argv[4]);
	}
	
	// argv[5] console output on(=1)/off(=0)
	if (argc > 5) {
		if (atoi(argv[5]) == 1) bConsole = true;
	}
	
	// argv[6] filename of ga parameters
	if (argc > 6) {
		try {
			loadGaParams(gaParams, argv[1]);
		}
		catch (...) {
			string warning = string("GA parameters file not found: ") + string(argv[1]) + string(" (default parameters initialized.)");
			gaParams = Solver_GA::getDefaultParams();
			TCB::logger.Log(Warning, warning);
		}
	} else {
		string warning = string("GA parameters filename not provided (default parameters initialized.)");
		gaParams = Solver_GA::getDefaultParams();
		TCB::logger.Log(Warning, warning);
	}
}

void sortJobsByD(vector<pJob>& unscheduledJobs) {
	sort(unscheduledJobs.begin(), unscheduledJobs.end(), compJobsByD);
}
void sortJobsByR(vector<pJob>& unscheduledJobs) {
	sort(unscheduledJobs.begin(), unscheduledJobs.end(), compJobsByR);
}
void sortJobsByGATC(vector<pJob>& unscheduledJobs, double t, double kappa) {
	double avgP = getAvgP(unscheduledJobs);
	sort(unscheduledJobs.begin(), unscheduledJobs.end(), CompJobsByGATC(avgP, t, kappa));
}

void sortJobsByRK(std::vector<pJob>& unscheduledJobs, const std::vector<double>& chr) {
	if (unscheduledJobs.size() != chr.size()) throw ExcSched("sortJobsByRK error different sizes of jobs and chromosome");

	vector<pJob> sortedJobs = vector<pJob>(chr.size());
	vector<pair<double, size_t>> ranking(chr.size());
	for (size_t j = 0; j < chr.size(); ++j) {
		ranking[j] = make_pair(chr[j], j);
	}
	sort(ranking.begin(), ranking.end());
	size_t pos = 0;
	for (auto rank = ranking.begin(); rank != ranking.end(); ++rank) {
		sortedJobs[pos] = move(unscheduledJobs[rank->second]);	// JENS assure that unscheduled size remains unchanged
		++pos;
	}

	for (size_t j = 0; j < chr.size(); ++j) {
		unscheduledJobs[j] = move(sortedJobs[j]);
	}
}

bool compJobsByD(const unique_ptr<Job>& a, const unique_ptr<Job>& b) {
	if (a->getD() == b->getD()) {
		return a->getId() < b->getId();
	}
	return a->getD() < b->getD();
}
bool compJobsByR(const unique_ptr<Job>& a, const unique_ptr<Job>& b) {
	if (a->getR() == b->getR()) {
		return a->getId() < b->getId();
	}
	return a->getR() < b->getR();
}

void shiftJobFromVecToVec(vector<pJob>& source, vector<pJob>& target, size_t sourceIdx) {
	if (sourceIdx >= source.size()) throw out_of_range("shiftJobFromVecToVec() out of range");
	target.push_back(move(source[sourceIdx]));
	source.erase(source.begin() + sourceIdx);
}

double getAvgP(const vector<pJob>& unscheduledJobs) {
	double totalP = 0.0;
	for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
		totalP += unscheduledJobs[j]->getTotalP();
	}
	return totalP / (double)unscheduledJobs.size();	// TODO check if to be divided by nJobs or nOps
}

void loadGaParams(GA_params& gaParams, std::string filename) {

}

vector<double> getDoubleGrid(double low, double high, double step) {
	vector<double> grid = vector<double>();
	while (low <= high + TCB::precision) {
		grid.push_back(low);
		low = low + step;
	}
	return grid;
}

