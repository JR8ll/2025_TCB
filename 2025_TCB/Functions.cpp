#include <algorithm>

#include "Functions.h"

using namespace std;

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

