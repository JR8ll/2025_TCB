#pragma once

#include<memory>
#include<vector>

#include "Workcenter.h"

using pWc = std::unique_ptr<Workcenter>;

class Schedule {
private:
	std::vector<pWc> workcenters;

public:
	Schedule();

	const std::vector<pWc>& getWorkcenters() const;
	void addWorkcenter(pWc wc);
};