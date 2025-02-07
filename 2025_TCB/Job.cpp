#include "Job.h"
#include "Operation.h"
#include "Product.h"

using namespace std;

Job::Job(int id, int s, Product* f, double r, double d, double w) : id(id), s(s), product(f), r(r), d(d), w(w) {
	ops = vector<pOp>();
}

int Job::getId() const { return id; }
int Job::getS() const { return s; }
int Job::getF() const { return product->getId(); }

double Job::getR() const { return r; }
double Job::getD() const { return d; }
double Job::getW() const { return w; }

double Job::getP(int stgIdx) const { return product->getP(stgIdx); }

void Job::addOp(pOp op) {
	ops.push_back(move(op));
}

const std::vector<pOp>& Job::getOps() const {
	return ops;
}
