#pragma once

#include<vector>

class Schedule;

class Solver_GA {
private:
	bool completed;
	std::vector<double> bestChromosome;

public:
	Solver_GA();
	~Solver_GA();

	double solveBRKGA_List_jobBased(Schedule& sched, int iTilimSeconds);

	bool hasCompleted();
	std::vector<double> getBestChromosome();

};

