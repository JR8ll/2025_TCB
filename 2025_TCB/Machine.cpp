#include "Machine.h"

using namespace std;

Machine::Machine(int id, Workcenter* wc) : id(id), workcenter(wc) {}

Workcenter* Machine::getWorkcenter() {
	return workcenter;
}
const std::vector<pBat>& Machine::getBatches() const {
	return batches;
}

void Machine::addBatch(pBat batch){
	batches.push_back(move(batch));
}