#include "Schedule.h"

using namespace std;

Schedule::Schedule() {
	workcenters = vector<pWc>();
}

const std::vector<pWc>& Schedule::getWorkcenters() const {
	return workcenters;
}
void Schedule::addWorkcenter(pWc wc) {
	workcenters.push_back(move(wc));
}