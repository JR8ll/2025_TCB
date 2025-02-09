#include <algorithm>

#include "Functions.h"

using namespace std;

void sortJobsByD(vector<pJob>& jobs) {
	sort(jobs.begin(), jobs.end(), compJobsByD);
}
void sortJobsByR(vector<pJob>& jobs) {
	sort(jobs.begin(), jobs.end(), compJobsByR);
}
void sortJobsByGATC(vector<pJob>& jobs, double t, double kappa) {
	double avgP = getAvgP(jobs);
	sort(jobs.begin(), jobs.end(), CompJobsByGATC(avgP, t, kappa));
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

double getAvgP(const vector<pJob>& jobs) {
	double totalP = 0.0;
	for (size_t j = 0; j < jobs.size(); ++j) {
		totalP += jobs[j]->getTotalP();
	}
	return totalP / (double)jobs.size();	// TODO check if to be divided by nJobs or nOps
}

