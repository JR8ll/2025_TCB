#pragma once

#include<vector>

#include "Job.h"
#include "Product.h"
#include "Schedule.h"

class Job;
class Product;
class Schedule;

class Problem {
private:
	Schedule env;
	std::vector<Product> products;
	std::vector<Job> jobs;

public:
	Problem();

};

