#pragma once

#include<memory>
#include<vector>

class Operation;

using pOp = std::unique_ptr<Operation>;

class Job {
private:
	std::vector<pOp> ops;

public:
	void addOp(pOp op);
	const std::vector<pOp>& getOps() const;
};
