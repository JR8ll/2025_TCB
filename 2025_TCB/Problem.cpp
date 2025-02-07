#include "Problem.h"

#include<fstream>

#include "Functions.h"

using namespace std;

Problem::Problem() : filename("n/a"), seed(0), n(0), stgs(0), F(0) {}
Problem::Problem(string filename) : filename(filename), seed(0), n(0), stgs(0), F(0) {
	this->loadFromDat(filename);
}

string Problem::getFilename() { return filename; }

int Problem::getN() const { return n; }
int Problem::getStgs() const { return stgs; }
int Problem::getF() const { return F; }

void Problem::loadFromDat(string filename) {
	ifstream input(filename);
	if (!input) {
		TCB::logger.Log(Error, "Could not open " + filename + ".");
		throw ExcSched("Could not open " + filename + ".");
	}
}
