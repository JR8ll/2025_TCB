#include "Product.h"

using namespace std;

Product::Product(int id, std::vector<int> route) : id(id), route(route) {
	tcMaxBwd = vector<vector<pair<int, double>>>(route.size());
	tcMaxFwd = vector<vector<pair<int, double>>>(route.size());
}

int Product::getId() const { return id; }
int Product::size() const { return route.size(); }

double Product::getP(int stgIdx) const { return p[stgIdx]; }

void Product::addTcMax(int first, int second, double duration) {
	tcMaxBwd[second].push_back(make_pair(first, duration));
	tcMaxFwd[first].push_back(make_pair(second, duration));
}