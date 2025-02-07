#include "Machine.h"
#include "Batch.h"
#include "Workcenter.h"

using namespace std;

Machine::Machine(int id, Workcenter* wc) : id(id), workcenter(wc) {}

Batch& Machine::operator[] (size_t idx) { return *batches[idx]; }
Batch& Machine::operator[] (size_t idx) const { return *batches[idx]; }

size_t Machine::size() const { return batches.size(); }

unique_ptr<Machine> Machine::clone(Workcenter* newWc) const {
	auto newMachine = make_unique<Machine>(id, newWc);
	for (const auto& bat : batches) {
		newMachine->addBatch(bat->clone());
	}
	return newMachine;
}

Workcenter* Machine::getWorkcenter() {
	return workcenter;
}
const std::vector<pBat>& Machine::getBatches() const {
	return batches;
}

void Machine::addBatch(pBat batch){
	batches.push_back(move(batch));
}