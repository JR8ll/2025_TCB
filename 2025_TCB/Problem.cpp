#include "Problem.h"

using namespace std;

Problem::Problem() {
	env = Schedule();
	products = vector<Product>();
	jobs = vector<Job>();
}