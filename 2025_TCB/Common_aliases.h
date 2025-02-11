#pragma once
#include<memory>
#include<vector>

class Batch;
class Machine;
class Job;
class Operation;
class Workcenter;

template<typename T>
using prioRule = void(*)(std::vector<T>& vec);

template<typename T>
using prioRuleKappa = void(*)(std::vector<T>& vec, double t, double kappa);

template<typename T>
using prioRuleKeySet = void(*)(std::vector<T>& vec, const std::vector<double>& keys);

using pOp = std::unique_ptr<Operation>;
using pBat = std::unique_ptr<Batch>;
using pMac = std::unique_ptr<Machine>;
using pWc = std::unique_ptr<Workcenter>;

using pJob = std::unique_ptr<Job>;
using sharedJob = std::shared_ptr<Job>;
using sharedOp = std::shared_ptr<Operation>;