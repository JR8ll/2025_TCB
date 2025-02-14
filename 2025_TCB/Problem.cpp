#include "Problem.h"

#include<fstream>

#include "Functions.h"
#include "Schedule.h"

using namespace std;

Problem::Problem() : filename("n/a"), seed(0), omega(0), n(0), stgs(0), F(0) {
	stages_1 = vector<int>();
	stages_b = vector<int>();
}
Problem::Problem(string filename) : filename(filename), seed(0), omega(0), n(0), stgs(0), F(0) {
	this->loadFromDat(filename);
}
Problem::Problem(ProbParams& params) {
	seed = params.seed;
	omega = params.omega;
	flowshop = params.flowshop;
	n = params.n;
	stgs = params.stgs;
	F = params.F;
	routes = params.routes;

	stages_1 = vector<int>();
	stages_b = vector<int>();

	m_o = vector<int>(stgs);
	m_B = vector<int>(stgs);

	bool bCapacityValues = params.m_BValues.size() == stgs;

	uniform_int_distribution<int> moDist(params.m_oIntervals.first, params.m_oIntervals.second);
	uniform_int_distribution<int> mbDist(params.m_BIntervals.first, params.m_BIntervals.second);
	for (int o = 0; o < stgs; ++o) {
		m_o[o] = moDist(TCB::rng);
		int myB = mbDist(TCB::rng);
		if (bCapacityValues) {
			myB = params.m_BValues[o];
		}
		m_B[o] = myB;
		if (myB <= 1) {
			stages_1.push_back(o + 1);
		}
		else {
			stages_b.push_back(o + 1);
		}
	}

	// all empty (no preformed batches)
	B_io = vector<vector<set<int> > >(F);
	B_iob = vector<vector<vector<int> > >(F);
	S_iob = vector<vector<vector<double> > >(F);
	for (int i = 0; i < F; ++i) {
		B_io[i] = vector<set<int> >(stgs);
		B_iob[i] = vector<vector<int> >(stgs);
		S_iob[i] = vector<vector<double> >(stgs);
		for (int o = 0; o < stgs; ++o) {
			B_io[i][o] = set<int>();
			B_iob[i][o] = vector<int>(n);
			S_iob[i][o] = vector<double>(n);
		}
	}
	rm = vector<vector<double> >(stgs);
	for (int o = 0; o < stgs; ++o) {
		rm[o] = vector<double>(m_o[o]);
		for (int m = 0; m < m_o[o]; ++m) {
			rm[o][m] = 0.0;
		}
	}

	// PROCESSING TIMES
	uniform_real_distribution<double> pDist(params.pInterval.first, params.pInterval.second);
	pTimes = vector<vector<double> >(F);
	for (int i = 0; i < F; ++i) {
		pTimes[i] = vector<double>(stgs);
		for (int o = 0; o < stgs; ++o) {
			pTimes[i][o] = pDist(TCB::rng);
		}
	}

	uniform_int_distribution<int> nTcDist(params.nTcInterval.first, params.nTcInterval.second);
	uniform_int_distribution<int> stepDist(0, (stgs - 1));
	
	// TIME CONSTRAINTS
	tc = vector<vector<vector <double> > >(F);
	for (int i = 0; i < F; ++i) {
		tc[i] = vector<vector<double> >(stgs);
		for (int o1 = 0; o1 < stgs; ++o1) {
			tc[i][o1] = vector<double>(stgs);
			for (int o2 = 0; o2 < stgs; ++o2) {
				tc[i][o1][o2] = 999999;	// dummy value (no constraint)
			}
		}
	}

	switch (params.tcScenario) {
	case 1:
		for (int i = 0; i < F; ++i) {
			set<pair<int, int> > constraints = set<pair<int, int> >();
			int nTc = nTcDist(TCB::rng);	// number of time constraints for this product
			for (int t = 0; t < nTc; ++t) {
				// select two stages arbitrarily
				int low = stepDist(TCB::rng);
				int high = stepDist(TCB::rng);
				if (high < low) {
					int temp = high;
					high = low;
					low = temp;
				}
				while (low == high || constraints.find(make_pair(low, high)) != constraints.end()) {
					high = stepDist(TCB::rng);
					low = stepDist(TCB::rng);
					if (high < low) {
						int temp = high;
						high = low;
						low = temp;
					}
				}
				constraints.insert(make_pair(low, high));

				// define tc length 
				double p = 0;
				for (int o = low; o < high; ++o) {
					p += pTimes[i][o];
				}
				double tcLength = p * params.tcFlowFactor;
				tc[i][low][high] = tcLength;
			}
		}
		break;
	default:
		TCB::logger.Log(Error, "missing or wrong parameter 'tcScenario' in Problem::Problem(...).");
	}

	double u = 0;
	for (int o = 0; o < stgs; ++o) {
		double uTemp = 0;
		for (int i = 0; i < F; ++i) {
			uTemp += pTimes[i][o] * ((double)n / (double)F);
		}
		uTemp = uTemp / (m_o[o] * F);
		if (uTemp > u) {
			u = uTemp;
		}
	}

	// CREATE PRODUCTS
	products = vector<Product>();
	for (size_t f = 0; f < F; ++f) {
		products.push_back(Product(f+1, routes[f]));
		products.back().setProcessingTimes(pTimes[f]);
		for (size_t o1 = 0; o1 < tc[f].size(); ++o1) {
			for (size_t o2 = 0; o2 < tc[f][o1].size(); ++o2) {
				products.back().addTcMax(o1, o2, tc[f][o1][o2]);
			}
		}
	}

	// SET JOB PARAMETERS
	uniform_real_distribution<double> probabilty(0.0, 1.0);
	uniform_real_distribution<double> rDist(params.rInterval.first, params.rInterval.second* u);
	uniform_real_distribution<double> wDist(params.wInterval.first, params.wInterval.second);
	uniform_real_distribution<double> ddFFDist(params.dueDateFF.first, params.dueDateFF.second);
	uniform_int_distribution<int> sDist(params.sInterval.first, params.sInterval.second);

	jobs_d = vector<double>(n);
	jobs_r = vector<double>(n);
	jobs_w = vector<double>(n);
	jobs_f = vector<int>(n);
	jobs_s = vector<int>(n);

	for (size_t j = 0; j < n; ++j) {
		double prob = probabilty(TCB::rng);
		if (prob < params.pReadyAtZero) {		// 25% of jobs are initially ready (Klemmt & Mönch)
			jobs_r[j] = 0.0;
		}
		else {
			jobs_r[j] = rDist(TCB::rng);
		}
		jobs_w[j] = wDist(TCB::rng);
		jobs_f[j] = ((j + F) % F) + 1;			// jobs are equally distributed accross products/families
		jobs_s[j] = sDist(TCB::rng);

		int nSteps = routes[jobs_f[j] - 1].size();
		double myD = jobs_r[j];
		double tempP = 0;
		for (size_t o = 0; o < nSteps; ++o) {
			tempP += pTimes[jobs_f[j] - 1][o];
		}
		double dueDateFF = ddFFDist(TCB::rng);
		myD += dueDateFF * tempP;				// Klemmt & Mönch: r_j + 2 x raw_processing_time
		jobs_d[j] = myD;
	}


	// CREATE JOBS
	unscheduledJobs = vector<pJob>();
	for (size_t j = 0; j < n; ++j) {
		pJob newJob = make_unique<Job>(j+1, jobs_s[j], &products[jobs_f[j]-1], jobs_r[j], jobs_d[j], jobs_w[j]);
		
		int nSteps = routes[jobs_f[j] - 1].size();
		for (size_t o = 0; o < nSteps; ++o) {
			pOp newOp = make_unique<Operation>(newJob.get(), o + 1);
			newJob->addOp(move(newOp));
		}
		for (size_t o = 0; o < newJob->size(); ++o) {
			if (o > 0) {
				(*newJob)[o].setPred(&(*newJob)[o - 1]);
			}
			if (o < (*newJob).size() - 1) {
				(*newJob)[o].setSucc(&(*newJob)[o + 1]);
			}
		}
		unscheduledJobs.push_back(move(newJob));
	}
	_setG();
}

string Problem::getFilename() { return filename; }

int Problem::getN() const { return n; }
int Problem::getStgs() const { return stgs; }
int Problem::getF() const { return F; }

Product* Problem::getProduct(size_t productIdx) {
	if (productIdx >= products.size()) throw out_of_range("Problem::getProduct() out of range");
	return &(products[productIdx]);

}

void Problem::loadFromDat(string filename) {
	this->filename = filename;
	ifstream input(filename);
	if (!input) {
		TCB::logger.Log(Error, "Could not open " + filename + ".");
		throw ExcSched("Could not open " + filename + ".");
	}
	filename = filename;
	string dummy;

	input >> dummy;	// "//seed="  
	input >> seed;	// <seed value>

	input >> dummy; // "n="
	input >> n;		// <n value>
	input >> dummy; // ";"

	input >> dummy;	// "m="
	input >> stgs;	// <m value>
	input >> dummy; // ";"

	input >> dummy;	// "F="
	input >> F;		// <F value>
	input >> dummy; // ";"

	input >> dummy;	// "omega="
	input >> omega;	// <omega value>
	input >> dummy; // ";"

	input >> dummy; // "Stages_1={"
	input >> dummy;	// <value>?

	while (dummy != "};") {
		stages_1.push_back(stoi(dummy));
		input >> dummy;
	}				// "};"

	input >> dummy;	// "Stages_b={"
	input >> dummy;	// <value>?

	int stages_bSize = 0;
	while (dummy != "};") {
		stages_b.push_back(stoi(dummy));
		stages_bSize++;
		input >> dummy;
	}				// "};"

	input >> dummy; // "m_o["
	for (size_t i = 0; i < stgs; ++i) {
		input >> dummy;
		m_o.push_back(stoi(dummy));
	}
	input >> dummy;	// "];"

	vector<int> m_B_temp = vector<int>(stages_bSize);
	m_B = vector<int>();

	rm = vector<vector<double> >(stgs);
	for (size_t i = 0; i < stgs; ++i) {
		rm[i] = vector<double>(m_o[i]);
	}

	input >> dummy;	// "rm=#["
	while (dummy != "]#;") {
		input >> dummy;
		if (dummy[0] == '<') {
			pair<int, int> stageMachine = Problem::_tokenizeTupel(dummy);
			input >> dummy;
			double ready = stod(dummy);
			rm[stageMachine.first - 1][stageMachine.second - 1] = ready;
		}
	}
	input >> dummy; // "]#;"

	B_io = vector<vector<set<int> > >(F);
	for (size_t i = 0; i < F; ++i) {
		B_io[i] = vector<set<int> >(stgs);
		for (size_t j = 0; j < stgs; ++j) {
			B_io[i][j] = set<int>();
			do {
				input >> dummy;
				int value;
				try {
					value = stoi(dummy);
					B_io[i][j].insert(value);
				}
				catch (invalid_argument const& exc) {
					break; // "B_io=[[{"
				}
			} while (dummy != "}{" && dummy != "}][{" && dummy != "}]];");
		}
	}

	input >> dummy; // "B_iob=[[["
	B_iob = vector<vector<vector<int> > >(F);
	for (size_t i = 0; i < F; ++i) {
		B_iob[i] = vector<vector<int> >(stgs);
		for (size_t j = 0; j < stgs; ++j) {
			B_iob[i][j] = vector<int>(n);
			for (size_t k = 0; k < n; ++k) {
				input >> dummy;
				B_iob[i][j][k] = stoi(dummy);
			}
			input >> dummy;		// "][" / "]][[" / "]]];
		}
	}

	input >> dummy; // "S_iob=["
	S_iob = vector<vector<vector<double> > >(F);
	for (size_t i = 0; i < F; ++i) {
		S_iob[i] = vector<vector<double> >(stgs);
		for (size_t j = 0; j < stgs; ++j) {
			S_iob[i][j] = vector<double>(n);
			for (size_t k = 0; k < n; ++k) {
				input >> dummy;
				S_iob[i][j][k] = stod(dummy);
			}
			input >> dummy;		// "][" / "]][[" / "]]];
		}
	}

	input >> dummy;	// "B=["
	for (size_t i = 0; i < stages_b.size(); ++i) {
		input >> dummy;
		m_B_temp[i] = stoi(dummy);
	}
	for (size_t i = 0; i < stgs; ++i) {
		bool bBatchingStage = false;
		for (size_t o = 0; o < stages_b.size(); ++o) {
			if (i + 1 == stages_b[o]) {
				bBatchingStage = true;	
				m_B.push_back(m_B_temp[o]);
			}	
		} 
		if (!bBatchingStage) {
			m_B.push_back(1);
		}
	}
	input >> dummy; // "];"

	input >> dummy;	// "d=["
	for (size_t i = 0; i < n; ++i) {
		input >> dummy;
		jobs_d.push_back(stod(dummy));
	}
	input >> dummy; // "];"

	input >> dummy;	// "r=["
	for (size_t i = 0; i < n; ++i) {
		input >> dummy;
		jobs_r.push_back(stod(dummy));
	}
	input >> dummy; // "];"

	input >> dummy;	// "w=["
	for (size_t i = 0; i < n; ++i) {
		input >> dummy;
		jobs_w.push_back(stod(dummy));
	}
	input >> dummy; // "];"

	input >> dummy;	// "f=["
	for (size_t i = 0; i < n; ++i) {
		input >> dummy;
		jobs_f.push_back(stoi(dummy));
	}
	input >> dummy; // "];"

	input >> dummy;	// "s=["
	for (size_t i = 0; i < n; ++i) {
		input >> dummy;
		jobs_s.push_back(stoi(dummy));
	}
	input >> dummy; // "];"

	input >> dummy;	// "rts=[{
	routes = vector<vector<int> >(F);
	for (size_t i = 0; i < F; ++i) {
		routes[i] = vector<int>();
		do {
			input >> dummy;
			int value;
			try {
				value = stoi(dummy);
				routes[i].push_back(value);
			}
			catch (invalid_argument const& exc) {
				break;
			}
		} while (dummy != "}{" && dummy != "}];");
	}

	for (size_t f = 0; f < F; ++f) {
		products.push_back(Product((f+1), routes[f]));
	}


	pTimes = vector < vector<double> >(F);
	for (size_t i = 0; i < F; ++i) {
		pTimes[i] = vector<double>(routes[i].size());
	}

	input >> dummy; // "p=[["
	for (size_t p = 0; p < F; ++p) {
		for (size_t s = 0; s < routes[p].size(); ++s) {
			input >> dummy;
			pTimes[p][s] = stod(dummy);
		}
		input >> dummy; // "]["/"]];"
	}

	tc = vector<vector<vector<double> > >(F);
	for (size_t p = 0; p < F; ++p) {
		tc[p] = vector<vector<double> >(routes[p].size());
		for (size_t s = 0; s < routes[p].size(); ++s) {
			tc[p][s] = vector<double>(routes[p].size());
		}
	}

	input >> dummy;	// "tc=[[[..."
	for (size_t p = 0; p < F; ++p) {
		// 1st dim: Product
		for (size_t s1 = 0; s1 < routes[p].size(); ++s1) {
			// 2nd dim: Stage
			for (size_t s2 = 0; s2 < routes[p].size(); ++s2) {
				// 3rd dim: Stage
				input >> dummy;
				tc[p][s1][s2] = stod(dummy);
			}
			input >> dummy;		// "]["
		}
	}

	for (size_t f = 0; f < F; ++f) {
		products[f].setProcessingTimes(pTimes[f]);
		for (size_t s1 = 0; s1 < tc[f].size(); ++s1) {
			for (size_t s2 = 0; s2 < tc[f][s1].size(); ++s2) {
				products[f].addTcMax(s1, s2, tc[f][s1][s2]);
			}
		}
	}

	// instantiate entities Jobs + Ops including time constraints
	for (size_t j = 0; j < n; ++j) {
		unscheduledJobs.push_back(unique_ptr<Job>(new Job((j + 1), jobs_s[j], &(products)[jobs_f[j] - 1], jobs_r[j], jobs_d[j], jobs_w[j])));
		for (size_t o = 0; o < products[jobs_f[j]-1].size(); ++o) {
			auto newOp = make_unique<Operation>(unscheduledJobs[j].get(), static_cast<int>(o + 1));
			unscheduledJobs[j]->addOp(move(newOp));
		}

		for (size_t o = 0; o < unscheduledJobs[j]->size(); ++o) {
			if (o > 0) {
				(*unscheduledJobs[j])[o].setPred(&(*unscheduledJobs[j])[o - 1]);
			}
			if (o < (*unscheduledJobs[j]).size() - 1) {
				(*unscheduledJobs[j])[o].setSucc(&(*unscheduledJobs[j])[o + 1]);
			}
		}
	}
	
	TCB::logger.Log(Info, "Problem initialized from " + filename + ".");
	input.close();
}

void Problem::_setG() {
	// TODO implement based on problem properties
	G = 999999;
}

pair<int, int> Problem::_tokenizeTupel(string tupel) {
	pair<int, int> ret = pair<int, int>();
	size_t first = 1;					// "<"
	size_t last = tupel.length() - 3;	// ">:"
	size_t comma = tupel.find(",");		// ","
	ret.first = stoi(tupel.substr(first, comma - first));
	ret.second = stoi(tupel.substr(comma + 1, last));
	return ret;
}

unique_ptr<Schedule> Problem::getSchedule() {
	unique_ptr<Schedule> newSchedule = make_unique<Schedule>();
	// setup machine environment
	for (size_t wc = 0; wc < stgs; ++wc) {
		unique_ptr<Workcenter> newWorkcenter = make_unique<Workcenter>(wc + 1, newSchedule.get());
		for (size_t m = 0; m < m_o[wc]; ++m) {
			unique_ptr<Machine> newMachine = make_unique<Machine>(m+1, m_B[wc], newWorkcenter.get());
			newWorkcenter->addMachine(move(newMachine));
		}
		newSchedule->addWorkcenter(move(newWorkcenter));
	}

	// add jobs
	for (size_t j = 0; j < n; ++j) {
		newSchedule->addJob(unscheduledJobs[j]->clone());
	}

	newSchedule->setProblemRef(this);
	return newSchedule;
}

