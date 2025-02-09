#include "Product.h"

using namespace std;

Product::Product(int id, std::vector<int> route) : id(id), route(route) {
	tcMaxBwd = vector<vector<pair<int, double>>>(route.size());
	tcMaxFwd = vector<vector<pair<int, double>>>(route.size());
}

int Product::getId() const { return id; }
int Product::size() const { return route.size(); }

const vector<pair<int, double>>& Product::getTcMaxBwd(size_t stgIdx) const {
	if (stgIdx >= tcMaxBwd.size()) throw out_of_range("Product::getTcMaxBwd(...) out of range");
	return tcMaxBwd[stgIdx];
}
const vector<pair<int, double>>& Product::getTcMaxFwd(size_t stgIdx) const {
	if (stgIdx >= tcMaxFwd.size()) throw out_of_range("Product::getTcMaxBwd(...) out of range");
	return tcMaxFwd[stgIdx];
}

int Product::getWorkcenterId(size_t stgIdx) const {
	if (stgIdx >= route.size()) throw out_of_range("Product::getWorkcenterId out of range");
	return route[stgIdx];
}
double Product::getP(size_t stgIdx) const { return p[stgIdx]; }

void Product::addTcMax(int first, int second, double duration) {
	tcMaxBwd[second].push_back(make_pair(first, duration));
	tcMaxFwd[first].push_back(make_pair(second, duration));
}