#include <algorithm>

#include "Functions.h"
#include "Problem.h"

using namespace std;

void processCmd(int argc, char* argv[], int& iSolver, int& iTilimSeconds, bool& bConsole) {
	if (argc > 1) {
		TCB::prob->loadFromDat(argv[1]);
	}
	else {
		TCB::prob->loadFromDat("debug_data.dat");
	}

	if (argc > 2) {
		TCB::seed = atoi(argv[2]);
		TCB::rng = mt19937(TCB::seed);
	}

	if (argc > 3) {
		iSolver = atoi(argv[3]);
	}
	if (argc > 4) {
		iTilimSeconds = atoi(argv[4]);
	}
	if (argc > 5) {
		if (atoi(argv[5]) == 1) bConsole = true;
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

vector<double> getDoubleGrid(double low, double high, double step) {
	vector<double> grid = vector<double>();
	while (low <= high + TCB::precision) {
		grid.push_back(low);
		low = low + step;
	}
	return grid;
}

