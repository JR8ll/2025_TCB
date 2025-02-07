#include "Job.h"

using namespace std;

void Job::addOp(pOp op) {
	ops.push_back(move(op));
}

const std::vector<pOp>& Job::getOps() const {
	return ops;
}
