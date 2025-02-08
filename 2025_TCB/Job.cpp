#include "Job.h"
#include "Operation.h"
#include "Product.h"

using namespace std;

Job::Job(int id, int s, Product* f, double r, double d, double w) : id(id), s(s), product(f), r(r), d(d), w(w) {
	ops = vector<pOp>();
}

unique_ptr<Job> Job::clone() const {
	auto newJob = make_unique<Job>(id, s, product.get(), r, d, w);
	for (const auto& op : ops) {
		auto newOp = make_unique<Operation>(newJob.get(), op->getStg());
		//newOp->assignToBatch(op->getBatch()); // copy not to be assigned to a batch
		newOp->setWait(op->getWait());
		newJob->addOp(std::move(newOp));
	}
	auto& newOps = newJob->getOps();
	for (size_t i = 0; i < newOps.size(); ++i) {
		if (i > 0) newOps[i]->setPred(newOps[i - 1].get());
		if (i < newOps.size() - 1) newOps[i]->setSucc(newOps[i + 1].get());
	}
	return newJob;
}

Operation& Job::operator[](size_t idx) { return *ops[idx]; }
Operation& Job::operator[](size_t idx) const { return *ops[idx]; }

size_t Job::size() const { return ops.size(); }

int Job::getId() const { return id; }
int Job::getS() const { return s; }
int Job::getF() const { return product->getId(); }

double Job::getR() const { return r; }
double Job::getD() const { return d; }
double Job::getW() const { return w; }
double Job::getP(int stgIdx) const { return product->getP(stgIdx); }

void Job::setD(double dueDate) { d = dueDate; }
void Job::setR(double release) { r = release; }
void Job::setW(double weight) { w = weight; }
void Job::setS(int size) { s = size; }

void Job::addOp(pOp op) {
	ops.push_back(move(op));
}

const std::vector<pOp>& Job::getOps() const {
	return ops;
}
