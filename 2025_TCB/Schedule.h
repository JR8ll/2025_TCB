#pragma once

#include<memory>
#include<vector>

class Workcenter;

using pWc = std::unique_ptr<Workcenter>;

class Schedule {
private:
	int id;
	std::vector<pWc> workcenters;

public:
	Schedule(int id);

	const std::vector<pWc> getWorkcenters() const;
	void addWorkcenter(pWc wc);
};