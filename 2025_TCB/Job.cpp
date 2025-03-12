#include "Job.h"
#include "Operation.h"
#include "Product.h"

using namespace std;

Job::Job(int id, int s, Product* f, double r, double d, double w) : id(id), s(s), product(f), r(r), d(d), w(w) {
	ops = vector<pOp>();
}
Job::~Job() {
	resetLinks();
}

unique_ptr<Job> Job::clone() const {
	auto newJob = make_unique<Job>(id, s, product, r, d, w);
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
int Job::getWorkcenterId(size_t stgIdx) const {
	return product->getWorkcenterId(stgIdx);
}

double Job::getR() const { return r; }
double Job::getD() const { return d; }
double Job::getW() const { return w; }
double Job::getP(size_t stgIdx) const { return product->getP(stgIdx); }
double Job::getTotalP() const {
	double totalP = 0.0;
	for (size_t o = 0; o < size(); ++o) {
		totalP += getP(o);
	}
	return totalP;
}

double Job::getC() const {
	return ops[ops.size()-1]->getC();
}
double Job::getStart() const {
	return ops[0]->getStart();
}
double Job::getWait() const {
	return ops[0]->getWait();
}


double Job::getGATC(double avgP, double t, double kappa) const {
	double gAtc = 0.0;
	for (size_t op = 0; op < size(); ++op) {
		gAtc += ops[op]->getGATC(avgP, t, kappa);
	}
	return gAtc;
}

const vector<std::pair<int, double>>& Job::getTcMaxFwd(size_t stgIdx) const {
	return product->getTcMaxFwd(stgIdx);
}
const vector<std::pair<int, double>>& Job::getTcMaxBwd(size_t stgIdx) const {
	return product->getTcMaxBwd(stgIdx);
}

void Job::setD(double dueDate) { d = dueDate; }
void Job::setR(double release) { r = release; }
void Job::setW(double weight) { w = weight; }
void Job::setS(int size) { s = size; }

void Job::setProduct(Product* prod) { 
	product = prod; 
}

Operation* Job::getOpPtr(size_t stgIdx) const {
	if (stgIdx >= ops.size()) throw out_of_range("Job::getOpPtr() out of range");
	return ops[stgIdx].get();
}

void Job::addOp(pOp op) {
	ops.push_back(move(op));
}

const std::vector<pOp>& Job::getOps() const {
	return ops;
}

void Job::resetLinks() {
	product = nullptr;
	ops.clear();
}
