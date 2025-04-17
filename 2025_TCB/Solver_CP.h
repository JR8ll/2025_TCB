#pragma once

class Schedule;

class Solver_CP {
public:
	Solver_CP();
	~Solver_CP();

	double solveCP(Schedule* schedule, int nDash = 5, int tilim = 60);
};