#pragma once

#include<memory>
#include<vector>

#include "Operation.h";
#include "Common_aliases.h"

class Product;

class Job {
private:
	int id;		
	int s;		// size

	double r;	// release time
	double d;	// due date
	double w;	// weight (priority)

	Product* product;
	std::vector<pOp> ops;

public:
	Job(int id, int s, Product* f, double r, double d, double w);
	~Job();

	std::unique_ptr<Job> clone() const;	// deep copy

	Operation& operator[] (size_t idx);
	Operation& operator[] (size_t idx) const;

	size_t size() const;

	int getId() const;
	int getS() const;
	int getF() const;
	int getWorkcenterId(size_t stgIdx) const;

	double getR() const;
	double getD() const;
	double getW() const;	
	double getP(size_t stgIdx) const;
	double getTotalP() const;

	double getGATC(double avgP, double t, double kappa) const;

	const std::vector<std::pair<int, double>>& getTcMaxFwd(size_t stgIdx) const;
	const std::vector<std::pair<int, double>>& getTcMaxBwd(size_t stgIdx) const;

	void setD(double dueDate);
	void setR(double release);
	void setW(double weight);
	void setS(int size);

	void setProduct(Product* prod);

	Operation* getOpPtr(size_t stgIdx) const;

	void addOp(pOp op);
	const std::vector<pOp>& getOps() const;

	void resetLinks();
};
