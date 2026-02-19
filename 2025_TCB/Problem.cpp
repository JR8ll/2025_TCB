#include "Problem.h"

#include<algorithm>
#include<fstream>
#include<iomanip>
#include<iostream>
#include<numeric>
#include<random>
#include<sstream>
#include <Windows.h>

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
Problem::Problem(ProbParams& params, bool bInteger) {
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
			double pTemp = pDist(TCB::rng);
			if (bInteger) {
				pTemp = floor(pTemp);
			}
			pTimes[i][o] = pTemp;
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
				if (bInteger) {
					tcLength = floor(tcLength);
				}
				tc[i][low][high] = tcLength;
			}
		}
		break;
	case 2:		// from any lower stage there is at most one constraint connecting it to a single higher stage
		for (int i = 0; i < F; ++i) {
			set<pair<int, int> > constraints = set<pair<int, int> >();
			set<int> firstStages = set<int>();
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
				while (low == high || firstStages.find(low) != firstStages.end()) {
					high = stepDist(TCB::rng);
					low = stepDist(TCB::rng);
					if (high < low) {
						int temp = high;
						high = low;
						low = temp;
					}
				}
				constraints.insert(make_pair(low, high));
				firstStages.insert(low);

				// define tc length 
				double p = 0;
				for (int o = low; o < high; ++o) {
					p += pTimes[i][o];
				}
				double tcLength = p * params.tcFlowFactor;
				if (bInteger) {
					tcLength = floor(tcLength);
				}
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
			double rTemp = rDist(TCB::rng);
			if (bInteger) {
				rTemp = floor(rTemp);
			}
			jobs_r[j] = rTemp;
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
		if (bInteger) {
			myD = floor(myD);
		}
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

Job& Problem::operator[](size_t idx) {
	return *unscheduledJobs[idx];
}
Job& Problem::operator[](size_t idx) const {
	return *unscheduledJobs[idx];
}

int Problem::size() const {
	return unscheduledJobs.size();
}

string Problem::getFilename() { return filename; }

int Problem::getN() const { return n; }
int Problem::getStgs() const { return stgs; }
int Problem::getF() const { return F; }

double Problem::getG() const {
	return G;
}

double Problem::getUpperBoundMSP() const {
	double msp = 0.0;
	double maxR = 0.0;
	for (size_t i = 0; i < jobs_r.size(); ++i) {
		for (size_t o = 0; o < stgs; ++o) {
			msp += (*unscheduledJobs[i])[o].getP();
		}
		if (jobs_r[i] > maxR) {
			maxR = jobs_r[i];
		}
	}
	msp += maxR;
	return msp;
}

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
		for (size_t o = 0; o < products[jobs_f[j] - 1].size(); ++o) {
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

	_setG();
	
	TCB::logger.Log(Info, "Problem initialized from " + filename + ".");
	input.close();
}

void Problem::saveToDat(string filename, Schedule* sched, ProbParams* params) {
	bool success = CreateDirectory(L".\\exp", NULL);		// exp directory usually contains problem instances and parameters
	
	string subfolder = filename;
	size_t extensionStart = subfolder.find_last_of('.');
	if (extensionStart != string::npos) {
		subfolder.erase(extensionStart);
	}
	

	ofstream output;
	output.open(".\\exp\\" + filename);
	output << "//seed= " << seed << endl;
	if (sched != nullptr) {
		output << "n= " << n - sched->getNumberOfScheduledJobs() << " ;" << endl;
	} else {
		output << "n= " << n << " ;" << endl;
	}
	output << "m= " << stgs << " ;" << endl;
	output << "F= " << F << " ;" << endl;
	output << endl;
	output << "omega= " << omega << " ;" << endl;
	output << endl;
	output << "Stages_1={ ";
	for (int i = 0; i < stages_1.size(); ++i) {
		output << stages_1[i] << " ";
	}
	output << "};" << endl;
	output << "Stages_b={ ";
	for (int i = 0; i < stages_b.size(); ++i) {
		output << stages_b[i] << " ";
	}
	output << "};" << endl;
	output << endl;
	output << "m_o=[ ";
	for (int i = 0; i < stgs; ++i) {
		output << m_o[i] << " ";
	}
	output << "];" << endl;

	// MACHINE RELEASE
	output << "rm=#[ ";
	if (sched != nullptr) {
		for (size_t wc = 0; wc < sched->size(); ++wc) {
			for (size_t m = 0; m < (*sched)[wc].size(); ++m) {
				output << "<" << (wc + 1) << "," << (m + 1) << ">: " << (*sched)[wc][m].getMSP();
				if (wc < (stgs - 1) || m < (m_o[wc] - 1)) {
					output << ", ";
				}
				else {
					output << " ";
				}
			}
		}
	}
	else {
		for (int i = 0; i < stgs; ++i) {
			for (int j = 0; j < m_o[i]; ++j) {
				output << "<" << (i + 1) << "," << (j + 1) << ">: " << rm[i][j];
				if (i < (stgs - 1) || j < (m_o[i] - 1)) {
					output << ", ";
				}
				else {
					output << " ";
				}
			}
		}
	}	
	output << "]#;" << endl;
	output << endl;
	output << endl;

	if (sched != nullptr) {
		// B_io is vector<vector<set>>> with [F][stgs] empty sets
		// B_iob is vector<vector<vector<int>>> with [F][stgs][n] (zero) int values
		// S_iob is vector<vector<vector<int>>> with [F][stgs][n] (zero) double values
		for (size_t f = 0; f < F; ++f) {
			for (size_t wc = 0; wc < sched->size(); ++wc) {
				int batchIdx = 1;
				for (size_t m = 0; m < (*sched)[wc].size(); ++m) {
					for (size_t b = 0; b < (*sched)[wc][m].size(); ++b) {
						// B_io
						if ((*sched)[wc][m][b].getF() == (f + 1)) {
							B_io[f][wc].insert(batchIdx);
							B_iob[f][wc][batchIdx - 1] = (*sched)[wc][m][b].getCap() - (*sched)[wc][m][b].getAvailableCap();
							S_iob[f][wc][batchIdx - 1] = (*sched)[wc][m][b].getStart();
						}
						++batchIdx;
					}
				}
			}
		}
	}

	// BATCH INDICES
	output << "B_io=[[{ ";
	for (int i = 0; i < B_io.size(); ++i) {
		for (int o = 0; o < B_io[i].size(); ++o) {
			for (auto it = B_io[i][o].begin(); it != B_io[i][o].end(); ++it) {
				output << *it << " ";
			}
			if (o < B_io[i].size() - 1) {
				output << "}{ ";
			}
		}
		if (i < B_io.size() - 1) {
			output << "}][{ ";
		}
	}	
	output << "}]]; " << endl;

	// BATCH OCCUPIED CAPACITIES
	output << "B_iob=[[[ ";
	for (int i = 0; i < B_iob.size(); ++i) {
		for (int o = 0; o < B_iob[i].size(); ++o) {
			for (auto it = B_iob[i][o].begin(); it != B_iob[i][o].end(); ++it) {
				output << *it << " ";
			}
			if (o < B_iob[i].size() - 1) {
				output << "][ ";
			}
		}
		if (i < B_iob.size() - 1) {
			output << "]][[ ";
		}
	}
	output << "]]]; " << endl;

	// BATCH START TIMES
	output << "S_iob=[[[ ";
	for (int i = 0; i < S_iob.size(); ++i) {
		for (int o = 0; o < S_iob[i].size(); ++o) {
			for (auto it = S_iob[i][o].begin(); it != S_iob[i][o].end(); ++it) {
				output << *it << " ";
			}
			if (o < S_iob[i].size() - 1) {
				output << "][ ";
			}
		}
		if (i < S_iob.size() - 1) {
			output << "]][[ ";
		}
	}
	output << "]]]; " << endl;

	output << "B=[ ";
	for (int i = 0; i < stgs; ++i) {
		if (find(stages_b.begin(), stages_b.end(), (i + 1)) != stages_b.end()) {
			output << m_B[i] << " ";
		}
	}
	output << "];" << endl;
	output << endl;
	output << endl;

	// TODO if sched != nullptr, only consider unscheduled jobs
	set<int> scheduledIndices = set<int>();
	if (sched != nullptr) {
		for (size_t j = 0; j < sched->getNumberOfScheduledJobs(); ++j) {
			scheduledIndices.insert(sched->getScheduledJob(j)->getId() - 1);
		}
	}

	output << "d=[ ";
	for (int i = 0; i < n; ++i) {
		if (scheduledIndices.find(i) != scheduledIndices.end()) {
			// this job is already scheduled
		}
		else {
			output << jobs_d[i] << " ";
		}
	}
	output << "];" << endl;

	output << "r=[ ";
	for (int i = 0; i < n; ++i) {
		if (scheduledIndices.find(i) != scheduledIndices.end()) {
			// this job is already scheduled
		}
		else {
			output << jobs_r[i] << " ";
		}
	}
	output << "];" << endl;

	output << "w=[ ";
	for (int i = 0; i < n; ++i) {
		if (scheduledIndices.find(i) != scheduledIndices.end()) {
			// this job is already scheduled
		}
		else {
			output << jobs_w[i] << " ";
		}
	}
	output << "];" << endl;

	output << "f=[ ";
	for (int i = 0; i < n; ++i) {
		if (scheduledIndices.find(i) != scheduledIndices.end()) {
			// this job is already scheduled
		}
		else {
			output << jobs_f[i] << " ";
		}
	}
	output << "];" << endl;

	output << "s=[ ";
	for (int i = 0; i < n; ++i) {
		if (scheduledIndices.find(i) != scheduledIndices.end()) {
			// this job is already scheduled
		}
		else {
			output << jobs_s[i] << " ";
		}
	}
	output << "];" << endl;
	output << endl;


	output << "rts=[{ ";
	for (int i = 0; i < F; ++i) {
		for (int step = 0; step < routes[i].size(); ++step) {
			output << routes[i][step] << " ";
			if (step == (routes[i].size() - 1) && i < (F - 1)) {
				output << "}{ ";
			}
		}
	}
	output << "}];" << endl;

	output << "p=[[ ";
	for (int i = 0; i < F; ++i) {
		for (int j = 0; j < stgs; ++j) {
			output << pTimes[i][j] << " ";
			if (i < (F - 1) && j == (stgs - 1)) {
				output << "][ ";
			}
		}
	}
	output << "]];" << endl;


	output << endl;
	output << endl;




	output << "tc=[[[ ";
	for (int p = 0; p < F; ++p) {
		for (int s = 0; s < stgs; ++s) {
			for (int t = 0; t < stgs; ++t) {
				output << tc[p][s][t] << " ";
				if (s < (stgs - 1) && t == (stgs - 1)) {
					output << "][ ";
				}
				if (p < (F - 1) && s == (stgs - 1) && t == (stgs - 1)) {
					output << "]][[ ";
				}
			}
		}
	}
	output << "]]];" << endl;
	output.close();
	if (params != nullptr) {
		createAutoSchedModelFiles("exp", subfolder, *params);
	}
}

//#include <tchar.h> // Include for _T macro

void Problem::createAutoSchedModelFiles(string topfolder, string subfolder, ProbParams& params) const {
	int lotSize = 25;
	double pTimeSpread = 0.2;	// processing time is uniformly distributed between p - (p * pTimeSpread) and p + (p * pTimeSpread)
	int mttf = 10080;			// mean time to failure (downcal.txt)
	int mttr = 120;				// mean time to repair (downcal.txt)
	double leadTimeFF = 1.3;	// raw processing time * leadTimeFF = lead time (due date is set := release + lead time in simulation)

	string fullpath = topfolder + "\\" + subfolder;
    wstring wideFilename = wstring(fullpath.begin(), fullpath.end());
    CreateDirectory(wideFilename.c_str(), NULL);

	string toolTxt = topfolder + "\\" + subfolder + "\\tool.txt";
	ofstream toolStream(toolTxt);
	toolStream << "STNFAM\t" << "STN\t" << "RULE\t" << "FWLRANK\t" << "WAKERESRANK\t" << "BATCHCRITF\t" << "BATCHPER\t" << "LTIME\t" << "LTUNITS\t" << "ULTIME\t" << "ULTUNITS\t" << "STNCAP\t" << "STNQTY\t" << "STNGRP\t" << "STNFAMSTEP_ACTLIST\t" << "STNFAMLOC\t" << "PRERULERWL\t" << "SETUPGRP" << endl;
	for (size_t o = 1; o <= stgs; ++o) {
		toolStream << "STAGE_" << to_string(o) << "\t";											//STNFAM
		toolStream << "STN_" << to_string(o) << "_" << to_string(1) << "\t";					//STN
		toolStream << "rule_FIRST" << "\t";														//RULE				TODO: custom rule 'rule_tcb'
		toolStream << "rank_FIFO" << "\t";														//FWLRANK			TODO: custom rank 'rank_tcb'
		toolStream << "\t";																		//WAKERESRANK
		toolStream << "crit_samepart" << "\t";													//BATCHCRITF		TODO: custom crit 'crit_tcb'
		toolStream << "lot" << "\t";															//BATCHPER
		toolStream << "1" << "\t";																//LTIME
		toolStream << "min" << "\t";															//LTNUNITS
		toolStream << "1" << "\t";																//ULTIME
		toolStream << "min" << "\t";															//ULTUNITS
		toolStream << to_string(m_B[o-1]) << "\t";												//STNCAP
		toolStream << to_string(m_o[o - 1]) << "\t";											//STNQTY
		toolStream << "GRP_" << to_string(o) << "\t";											//STNGRP
		toolStream << "Custom_actlist_ASISemiOpersDuringSetupAndAdditionalLoadUnload" << "\t";	//STNFAMSTEP_ATCLIST
		toolStream << "Fab" << "\t";															//STNFAMLOC
		toolStream << "no" << "\t";																//PRERULERWL
		toolStream << "\t";																		//SETUPGRP
		toolStream << endl;
	}
	toolStream.close();

	string partTxt = topfolder + "\\" + subfolder + "\\part.txt";
	ofstream partStream(partTxt);
	partStream << "PARTGRP\t" << "PARTFAM\t" << "PART\t" << "ROUTEFILE\t" << "ROUTE\t" << "DFLTLD\t" << "DFLTLDU" << endl;
	for (size_t i = 1; i <= F; ++i) {
		partStream << "Saleable\t";
		partStream << "product_" << to_string(i) << "\t";
		partStream << "part_" << to_string(i) << "\t";
		partStream << "route_" << to_string(i) << ".txt" << "\t";
		partStream << "r_" << to_string(i) << "\t";
		int leadTime = 0;
		uniform_real_distribution<double> ddFFDist(params.dueDateFF.first, params.dueDateFF.second);
		double tempP = 0;
		for (size_t o = 0; o < stgs; ++o) {
			tempP += pTimes[i - 1][o];
		}
		double dueDateFF = ddFFDist(TCB::rng);
		leadTime = tempP * dueDateFF;
		partStream << to_string(leadTime) << "\t";
		partStream << "min";
		partStream << endl;
	}
	partStream.close();

	string orderTxt = topfolder + "\\" + subfolder + "\\order.txt";
	ofstream orderStream(orderTxt);
	orderStream << "LOT\t" << "PART\t" << "PRIOR\t" << "PIECES\t" << "START\t" << "RDIST\t" << "REPEAT\t" << "RUNITS\t" << "RPT#\t" << "LOTSPERRPT\t" << "ORDER\t" << "HOTLOT\t" << endl;
	for (size_t i = 1; i <= F; ++i) {
		orderStream << "Lot_" << to_string(i) << "\t";
		orderStream << "part_" << to_string(i) << "\t";
		orderStream << "10" << "\t";
		orderStream << to_string(lotSize) << "\t";
		orderStream << "01/01/2018 00:00:00" << "\t";		// MM / DD / YYYY HH : MM:SS
		orderStream << "exponential" << "\t";					// TODO: define interarrival distribution
		orderStream << "50" << "\t";						// TODO: define interarrival mean
		orderStream << "min" << "\t";
		orderStream << "200000" << "\t";
		orderStream << "1" << "\t";
		//orderStream << "01/15/2018 10:05:00" << "\t";		// MM/DD/YYYY HH:MM:SS		due dates are defined be default lead time set in part.txt
		orderStream << "O_Lot_" << to_string(i) << "\t";
		orderStream << "no" << "\t";
		orderStream << endl;
	}
	orderStream.close();

	string downcalTxt = topfolder + "\\" + subfolder + "\\downcal.txt";
	ofstream downcalStream(downcalTxt);
	downcalStream << "DOWNCALNAME\t" << "DOWNCALTYPE\t" << "MTTFDIST\t" << "MTTF\t" << "MTTFUNITS\t" << "MTTRDIST\t" << "MTTR\t" << "MTTRUNITS\t" << endl;
	for (size_t i = 1; i <= stgs; ++i) {
		downcalStream << "BREAK_" << to_string(i) << "\t";
		downcalStream << "mttf_by_cal" << "\t";
		downcalStream << "exponential" << "\t";
		downcalStream << to_string(mttf) << "\t";
		downcalStream << "min" << "\t";
		downcalStream << "exponential" << "\t";
		downcalStream << to_string(mttr) << "\t";
		downcalStream << "min" << "\t";
		downcalStream << endl;
	}
	downcalStream.close();

	string attachTxt = topfolder + "\\" + subfolder + "\\attach.txt";
	ofstream attachStream(attachTxt);
	attachStream << "CALNAME\t" << "CALTYPE\t" << "RESTYPE\t" << "RESNAME\t" << "FOADIST\t" << "FOA\t" << "FOAUNITS\t" << endl;
	for (size_t i = 1; i <= stgs; ++i) {
		attachStream << "BREAK_" << to_string(i) << "\t";
		attachStream << "down" << "\t";
		attachStream << "stngrp" << "\t";
		attachStream << "GRP_" << to_string(i) << "\t"; 
		attachStream << "exponential" << "\t";
		attachStream << to_string(mttf) << "\t";
		attachStream << "min" << "\t";
		attachStream << endl;
	}
	attachStream.close();

	for (size_t i = 1; i <= F; ++i) {
		string routeTxt = topfolder + "\\" + subfolder + "\\route_" + to_string(i) + ".txt";
		ofstream routeStream(routeTxt);
		routeStream << "ROUTE\t" << "STEP\t" << "DESC\t" << "STNFAM\t" << "PDIST\t" << "PTIME\t" << "PTIME2\t" << "PTUNITS\t" << "PTPER\t" << "BATCHMN\t" << "BATCHMX\t" << "SETUP\t" << "WHEN\t" << "STIME\t" << "STUNITS\t" << "SVESTN\t" << "FORSTEP\t" << "BatchInterval\t" << "BatchIntUnits\t" << "PartInterval\t" << "PartIntUnits\t" << "RWKSTEP\t" << "REWORK\t" << "RWKTYPE\t" << "StepPercent\t" << "STEP_CQT\t" << "CQT\t" << endl;
		for (size_t o = 1; o <= stgs; ++o) {
			routeStream << "r_" << to_string(i) << "\t";								//ROUTE	
			routeStream << to_string(o) << "\t";										// STEP	
			routeStream << "stage_" << to_string(o) << "\t";							// DESC	
			routeStream << "STAGE_" << to_string(o) << "\t";							// STNFAM	
			routeStream << "uniform" << "\t";											// PDIST	
			routeStream << to_string(pTimes[i-1][o-1]) << "\t";							// PTIME		
			routeStream << to_string(pTimes[i - 1][o - 1] * pTimeSpread) << "\t";		// PTIME2		
			routeStream << "min" << "\t";												// PTUNITS	
			bool batchingStage = m_B[o - 1] > 1;
			if (batchingStage) {
				routeStream << "per_batch" << "\t";										// PTPER
				routeStream << to_string(1) << "\t";									// BATCHMN	
				routeStream << to_string(m_B[o-1]) << "\t";								// BATCHMX	
			}
			else {
				routeStream << "per_lot" << "\t";										// PTPER
				routeStream << "\t";													// BATCHMN	
				routeStream << "\t";													// BATCHMX	
			}	
			routeStream << "\t";	// SETUP	
			routeStream << "\t";	// WHEN	
			routeStream << "\t";	// STIME	
			routeStream << "\t";	// STUNITS	
			routeStream << "\t";	// SVESTN	
			routeStream << "\t";	// FORSTEP	
			routeStream << "\t";	// BatchInterval	
			routeStream << "\t";	// BatchIntUnits	
			routeStream << "\t";	// PartInterval	
			routeStream << "\t";	// PartIntUnits	
			routeStream << "\t";	// RWKSTEP	
			routeStream << "\t";	// REWORK	
			routeStream << "\t";	// RWKTYPE	
			routeStream << "\t";	// StepPercent	

			//std::vector<std::vector<std::vector<double>>> tc;	// time constraints [product][stage1][stage2]
			bool tcExists = false;
			size_t nTc = tc[i - 1][o - 1].size();
			vector<size_t> seq(nTc);
			iota(seq.begin(), seq.end(), 0);
			shuffle(seq.begin(), seq.end(), TCB::rng);
			for (size_t step = 0; step < nTc; ++step) {
				double tcValue = tc[i - 1][o - 1][seq[step]];
				if (tcValue != 999999) {
					routeStream << to_string(seq[step] + 1) << "\t";	// STEP_CQT	
					routeStream << to_string(tcValue) << "\t";	// CQT
					tcExists = true;
					break;
				}
			}
			if (!tcExists) {
				routeStream << "\t";	// STEP_CQT	
				routeStream << "\t";	// CQT
			}
			
			routeStream << endl;
		}
		routeStream.close();
	}
}

double Problem::getAvgWorkloadPerMachine(int stageIdx) {
	if (stageIdx < 0 || stageIdx >= stgs) throw out_of_range("Problem::getAvgWorkloadPerMachine(...) out of range.");
	double totalProcessingTime = 0.0;
	for (size_t j = 0; j < jobs_f.size(); ++j) {
		totalProcessingTime += pTimes[jobs_f[j] - 1][stageIdx] * (double) jobs_s[j] / (double) m_B[stageIdx];	// CONSIDERING JOB AND BATCH SIZE
	}
	return totalProcessingTime / (double)m_o[stageIdx];
}

void Problem::configureBottleneck(pair<int, int> bottleneckStageIndices, double bottleneckCriticality, ProbParams& params, bool integerValues) {
	if (bottleneckStageIndices.first < 0 || bottleneckStageIndices.first >= stgs || bottleneckStageIndices.second < 0 || bottleneckStageIndices.second >= stgs) throw out_of_range("Problem::configureBottleneck(...) out of range.");
	uniform_int_distribution<int> bnStageDist = uniform_int_distribution<int>(bottleneckStageIndices.first, bottleneckStageIndices.second);
	int bottleneckStageIdx = bnStageDist(TCB::rng);

	// get basic stage (= stage with lowest workload)
	int basicStageIdx = 0;
	int maxDistance = 0;
	for (int i = 0; i < stgs; ++i) {	// arg max(o in 1...stgs) |o - o_bottleneck|
		int tempDistance = abs(bottleneckStageIdx - i);
		if (tempDistance > maxDistance) {
			maxDistance = tempDistance;
			basicStageIdx = i;
		}
	}

	double baseWorkload = getAvgWorkloadPerMachine(basicStageIdx);

	// adapt workload of all stages (but the basic stage)
	for (int i = 0; i < stgs; ++i) {
		if (i != basicStageIdx) {
			double myWorkload = getAvgWorkloadPerMachine(i);
			double actualWorkloadRatio = myWorkload / baseWorkload;
			double targetWorkloadRatio = ((maxDistance - abs(bottleneckStageIdx - i))  * bottleneckCriticality) + 1;
			for (int f = 0; f < F; ++f) {
				pTimes[f][i] *= (targetWorkloadRatio / actualWorkloadRatio);
				if (integerValues) {
					pTimes[f][i] = round(pTimes[f][i]);
				}
			}
		}
	}

	// adapt maximal time lags to adapted processing times
	for (size_t i = 0; i < tc.size(); ++i) {
		for (size_t lo = 0; lo < tc[i].size(); ++lo) {
			for (size_t hi = lo + 1; hi < tc[i][lo].size(); ++hi) {
				if (tc[i][lo][hi] != 999999) {
				double pRaw = 0.0;
				for (size_t k = lo; k < hi; ++k) {
					pRaw += pTimes[i][k];
				}
					tc[i][lo][hi] = pRaw * params.tcFlowFactor;
				}
			}
		}
	}
	if (!assertFeasibility()) throw ExcSched("INFEASIBLE MAXIMAL TIME LAGS");

	// adapt due dates to adapted processing times
	uniform_real_distribution<double> ddFFDist(params.dueDateFF.first, params.dueDateFF.second);
	for (size_t j = 0; j < n; ++j) {
		int nSteps = routes[jobs_f[j] - 1].size();
		double myD = jobs_r[j];
		double tempP = 0;
		for (size_t o = 0; o < nSteps; ++o) {
			tempP += pTimes[jobs_f[j] - 1][o];
		}
		double dueDateFF = ddFFDist(TCB::rng);
		myD += dueDateFF * tempP;				// Klemmt & Mönch: r_j + 2 x raw_processing_time
		if (integerValues) {
			myD = floor(myD);
		}
		jobs_d[j] = myD;
		unscheduledJobs[j]->setD(myD);
	}
	

	// DEBUGGING: output average workload of all stages
	//for (int i = 0; i < stgs; ++i) {
	//	cout << "Stage " << (i + 1) << ": Workload/machine = " << fixed << setprecision(1) << getAvgWorkloadPerMachine(i);
	//	if (i == basicStageIdx) cout << " (basic stage)";
	//	if (i == bottleneckStageIdx) cout << " (bottleneck stage)";
	//	cout << endl;
	//}
	/*int debugger = 666;*/
}

void Problem::_setG() {
	double bigInteger = 0;
	for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
		for (size_t o = 0; o < (*unscheduledJobs[j]).size(); ++o) {
			bigInteger += (*unscheduledJobs[j])[o].getP();
		}
	}

	G = bigInteger;
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

bool Problem::assertFeasibility() {
	// check 1 - maximal time lags > raw processing times
	for (size_t i = 0; i < tc.size(); ++i) {
		for (size_t lo = 0; lo < tc[i].size(); ++lo) {
			for (size_t hi = lo + 1; hi < tc[i][lo].size(); ++hi) {
				double pRaw = 0.0;
				for (size_t k = lo; k < hi; ++k) {
					pRaw += pTimes[i][k];
				}
				if (pRaw > tc[i][lo][hi]) {
					cout << "ERROR: infeasible problem instance, tc[" << i << "][" << lo << "][" << hi << "] is smaller than the raw procssing time" << endl;
					return false;
				}
			}
		}
	}

	return true;
}

void Problem::genInstancesTCB25_Feb25_exact() {
	ProbParams params;
	params.omega = 9;
	params.F = 2;
	params.stgs = 5;
	params.n = 10;
	params.m_oIntervals = make_pair(2, 3);
	params.m_BIntervals = make_pair(1, 3);
	//params.m_BValues = vector<int>({ 3, 3, 3, 3, 3 });
	params.pInterval = make_pair(10, 25);
	params.tcScenario = 1;
	params.tcFlowFactor = 1.5;
	params.rInterval = make_pair(0, 0.75);
	params.sInterval = make_pair(1, 1);	// uniform job sizes
	params.wInterval = make_pair(1.0, 3.0);
	params.dueDateFF = make_pair(1.0, 1.3);

	int nmax_tc = 0;	// maximum possible number of timeconstraints = sum(i in 0..stgs) i
	for (int i = 1; i < params.stgs; ++i) {
		nmax_tc += i;
	}
	int nmin_tc = params.stgs - 1; // minimum number of timeconstraints = number of stages (-1)
	params.nTcInterval = make_pair(nmin_tc, nmax_tc);

	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}

	int nInstances = 10;
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params);
		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "ProbI_DS3TC_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "_exact.dat";
		prob.saveToDat(fileName);
	}

	params.F = 3;
	params.n = 9;
	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params);
		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "ProbI_DS3TC_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "_exact.dat";
		prob.saveToDat(fileName);
	}
}

void Problem::genInstancesTCB25_Mar25_discr() {
	ProbParams params;
	params.omega = 9;
	params.F = 5;
	params.stgs = 5;
	params.n = 75;
	params.m_oIntervals = make_pair(2, 4);
	params.m_BIntervals = make_pair(4, 4);
	params.m_BValues = vector<int>({ 4, 1, 1, 1, 1 });
	params.pInterval = make_pair(10, 25);
	params.tcScenario = 1;
	params.tcFlowFactor = 1.5;
	params.rInterval = make_pair(0, 0.75);
	params.sInterval = make_pair(1, 1);	// uniform job sizes
	params.wInterval = make_pair(1.0, 3.0);
	params.dueDateFF = make_pair(1.0, 1.3);
	params.pReadyAtZero = 0.25;
	int nmax_tc = 0;	// maximum possible number of timeconstraints = sum(i in 0..stgs) i
	for (int i = 1; i < params.stgs; ++i) {
		nmax_tc += i;
	}
	int nmin_tc = params.stgs - 1; // minimum number of timeconstraints = number of stages (-1)
	params.nTcInterval = make_pair(nmin_tc, nmax_tc);
	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}
	int nInstances = 10;
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);	// true for discrete time values
		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "ProbI_MISTA_DISCR_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + ".dat";
		prob.saveToDat(fileName);
	}
}

void Problem::genInstancesEURO25() {
	ProbParams params;
	params.omega = 9;
	params.F = 5;
	params.stgs = 10;
	params.n = 50;
	params.m_oIntervals = make_pair(2, 4);
	params.m_BIntervals = make_pair(1, 3);
	//params.m_BValues = vector<int>({ 3, 3, 3, 3, 3 });
	params.pInterval = make_pair(10, 25);
	params.tcScenario = 1;
	params.tcFlowFactor = 1.5;
	params.rInterval = make_pair(0, 0.75);
	params.sInterval = make_pair(1, 1);	// uniform job sizes
	params.wInterval = make_pair(1.0, 3.0);
	params.dueDateFF = make_pair(1.0, 1.3);
	params.pReadyAtZero = 0.25;

	int nmax_tc = 0;	// maximum possible number of timeconstraints = sum(i in 0..stgs) i
	for (int i = 1; i < params.stgs; ++i) {
		nmax_tc += i;
	}
	int nmin_tc = params.stgs - 1; // minimum number of timeconstraints = number of stages (-1)
	params.nTcInterval = make_pair(nmin_tc, nmax_tc);

	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}

	int nInstances = 10;
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params);
		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "ProbI_EURO_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + ".dat";
		prob.saveToDat(fileName, nullptr, &params);
	}

	// 2nd factor combination

	params.F = 10;
	params.n = 100;
	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params);
		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "ProbI_EURO_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + ".dat";
		prob.saveToDat(fileName, nullptr, &params);
	}

	// 3rd factor combination

	params.F = 5;
	params.n = 50;
	params.tcFlowFactor = 3.0;
	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params);
		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "ProbI_EURO_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + ".dat";
		prob.saveToDat(fileName, nullptr, &params);
	}


	// 4th factor combination

	params.F = 10;
	params.n = 100;
	params.tcFlowFactor = 3.0;
	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params);
		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "ProbI_EURO_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + ".dat";
		prob.saveToDat(fileName, nullptr, &params);
	}
}

void Problem::genInstancesEURO25_integer() {
	ProbParams params;
	params.omega = 9;
	params.F = 5;
	params.stgs = 10;
	params.n = 50;
	params.m_oIntervals = make_pair(2, 4);
	params.m_BIntervals = make_pair(1, 3);
	//params.m_BValues = vector<int>({ 3, 3, 3, 3, 3 });
	params.pInterval = make_pair(10, 25);
	params.tcScenario = 1;
	params.tcFlowFactor = 1.5;
	params.rInterval = make_pair(0, 0.75);
	params.sInterval = make_pair(1, 1);	// uniform job sizes
	params.wInterval = make_pair(1.0, 3.0);
	params.dueDateFF = make_pair(1.0, 1.3);
	params.pReadyAtZero = 0.25;

	int nmax_tc = 0;	// maximum possible number of timeconstraints = sum(i in 0..stgs) i
	for (int i = 1; i < params.stgs; ++i) {
		nmax_tc += i;
	}
	int nmin_tc = params.stgs - 1; // minimum number of timeconstraints = number of stages (-1)
	params.nTcInterval = make_pair(nmin_tc, nmax_tc);

	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}

	int nInstances = 10;
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "ProbI_EURO_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + ".dat";
		prob.saveToDat(fileName, nullptr, &params);
	}

	// 2nd factor combination

	params.F = 10;
	params.n = 100;
	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "ProbI_EURO_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + ".dat";
		prob.saveToDat(fileName, nullptr, &params);
	}

	// 3rd factor combination

	params.F = 5;
	params.n = 50;
	params.tcFlowFactor = 3.0;
	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "ProbI_EURO_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + ".dat";
		prob.saveToDat(fileName, nullptr, &params);
	}


	// 4th factor combination

	params.F = 10;
	params.n = 100;
	params.tcFlowFactor = 3.0;
	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "ProbI_EURO_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + ".dat";
		prob.saveToDat(fileName, nullptr, &params);
	}
}

void Problem::genInstancesEURO25_exact(){
	ProbParams params;
	params.omega = 9;
	params.F = 2;
	params.stgs = 4;
	params.n = 10;
	params.m_oIntervals = make_pair(1, 2);
	//params.m_BIntervals = make_pair(1, 3);
	//params.m_BValues = vector<int>({ 3, 3, 3, 3, 3 });
	params.pInterval = make_pair(10, 25);
	params.tcScenario = 1;
	params.tcFlowFactor = 1.5;
	params.rInterval = make_pair(0, 0.75);
	params.sInterval = make_pair(1, 1);	// uniform job sizes
	params.wInterval = make_pair(1.0, 3.0);
	params.dueDateFF = make_pair(1.0, 1.3);
	params.pReadyAtZero = 0.25;

	int nmax_tc = 0;	// maximum possible number of timeconstraints = sum(i in 0..stgs) i
	for (int i = 1; i < params.stgs; ++i) {
		nmax_tc += i;
	}
	int nmin_tc = params.stgs - 1; // minimum number of timeconstraints = number of stages (-1)
	params.nTcInterval = make_pair(nmin_tc, nmax_tc);

	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}

	// just one randomly chosen batching stage
	uniform_int_distribution<int> batchingStageDist(0, params.stgs - 1);
	int nInstances = 30;
	for (int i = 0; i < nInstances; ++i) {
		int batchingStageIdx = batchingStageDist(TCB::rng);
		params.m_BValues = vector<int>();
		for (int stg = 0; stg < params.stgs; ++stg) {
			if (stg == batchingStageIdx) {
				params.m_BValues.push_back(3);
			} else {
				params.m_BValues.push_back(1);
			}
		}

		Problem prob = Problem(params);
		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "ProbI_EURO_exact_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + ".dat";
		prob.saveToDat(fileName);
	}
}

void Problem::genInstancesTCB25_Jun25_exactMILPvsCP() {
	ProbParams params;
	params.omega = 9;
	params.F = 2;
	params.stgs = 3;
	params.n = 8;
	params.m_oIntervals = make_pair(2, 3);
	params.m_BIntervals = make_pair(1, 3);
	//params.m_BValues = vector<int>({ 3, 3, 3, 3, 3 });
	params.pInterval = make_pair(10, 25);
	params.tcScenario = 1;
	params.tcFlowFactor = 1.5;
	params.rInterval = make_pair(0, 0.75);
	params.sInterval = make_pair(1, 1);	// uniform job sizes
	params.wInterval = make_pair(1.0, 3.0);
	params.dueDateFF = make_pair(1.0, 1.3);

	int nmax_tc = 0;	// maximum possible number of timeconstraints = sum(i in 0..stgs) i
	for (int i = 1; i < params.stgs; ++i) {
		nmax_tc += i;
	}
	int nmin_tc = params.stgs - 1; // minimum number of timeconstraints = number of stages (-1)
	params.nTcInterval = make_pair(nmin_tc, nmax_tc);

	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}

	int nInstances = 10;
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "ProbI_TCB_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "_exact.dat";
		prob.saveToDat(fileName);
	}

	params.F = 3;
	params.n = 9;
	params.stgs = 4;
	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "ProbI_TCB_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "_exact.dat";
		prob.saveToDat(fileName);
	}

	params.F = 2;
	params.n = 10;
	params.stgs = 5;
	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "ProbI_TCB_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "_exact.dat";
		prob.saveToDat(fileName);
	}
}

void Problem::genInstancesTCB26_Testing() {
	ProbParams params;
	params.omega = 9;
	params.F = 5;
	params.stgs = 7;
	params.n = 100;
	params.m_oIntervals = make_pair(4, 10);
	params.m_BIntervals = make_pair(2, 6);
	//params.m_BValues = vector<int>({ 3, 3, 3, 3, 3 });
	params.pInterval = make_pair(10, 25);
	params.tcScenario = 2;
	params.tcFlowFactor = 1.5;
	params.rInterval = make_pair(0, 0.75);
	params.sInterval = make_pair(1, 1);	// uniform job sizes
	params.wInterval = make_pair(1.0, 3.0);
	params.dueDateFF = make_pair(1.0, 1.3);

	
	int nmin_tc = params.stgs / 3.0; 
	int nmax_tc = params.stgs - 1; // minimum number of timeconstraints = number of stages (-1)
	params.nTcInterval = make_pair(nmin_tc, nmax_tc);

	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}

	int nInstances = 10;
	double bottleneckCriticality = 0.125;

	// n = 100, 7 stages, tcFF = 1.5, bottlenek 1st stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(0, 0), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckFirst.dat";
		prob.saveToDat(fileName);
	}

	// n = 100, 7 stages, tcFF = 1.5, bottleneck middle stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(params.stgs-1, params.stgs - 1), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckLast.dat";
		prob.saveToDat(fileName);
	}

	// n = 100, 9 stages, tcFF = 1.5, bottleneck last stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair((params.stgs - 1) / 2, (params.stgs - 1) / 2), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckMdl.dat";
		prob.saveToDat(fileName);
	}


	params.tcFlowFactor = 3.0;
	// n = 100, 7 stages, tcFF = 3.0, bottlenek 1st stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(0, 0), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckFirst.dat";
		prob.saveToDat(fileName);
	}

	//  n = 100, 7 stages, tcFF = 3.0, bottleneck middle stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(params.stgs - 1, params.stgs - 1), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckLast.dat";
		prob.saveToDat(fileName);
	}

	//  n = 100, 7 stages, tcFF = 3.0, bottleneck last stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair((params.stgs - 1) / 2, (params.stgs - 1) / 2), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckMdl.dat";
		prob.saveToDat(fileName);
	}

	// 15 stages
	params.stgs = 15;
	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}
	params.tcFlowFactor = 1.5;
	// n = 100, 15 stages, tcFF = 1.5, bottlenek 1st stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(0, 0), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckFirst.dat";
		prob.saveToDat(fileName);
	}

	// n = 100, 15 stages, tcFF = 1.5, bottleneck middle stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(params.stgs - 1, params.stgs - 1), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckLast.dat";
		prob.saveToDat(fileName);
	}

	// n = 100, 15 stages, tcFF = 1.5, bottleneck last stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair((params.stgs - 1) / 2, (params.stgs - 1) / 2), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckMdl.dat";
		prob.saveToDat(fileName);
	}

	
	params.tcFlowFactor = 3.0;
	// n = 100, 15 stages, tcFF = 3.0, bottlenek 1st stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(0, 0), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckFirst.dat";
		prob.saveToDat(fileName);
	}

	// n = 100, 15 stages, tcFF = 3.0, bottleneck middle stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(params.stgs - 1, params.stgs - 1), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckLast.dat";
		prob.saveToDat(fileName);
	}

	// n = 100, 15 stages, tcFF = 3.0, bottleneck last stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair((params.stgs - 1) / 2, (params.stgs - 1) / 2), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckMdl.dat";
		prob.saveToDat(fileName);
	}

	params.n = 200;
	params.tcFlowFactor = 1.5;
	params.stgs = 7;
	// n = 200, 7 stages, tcFF = 1.5, bottlenek 1st stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(0, 0), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckFirst.dat";
		prob.saveToDat(fileName);
	}

	// n = 200, 7 stages, tcFF = 1.5, bottleneck middle stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(params.stgs - 1, params.stgs - 1), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckLast.dat";
		prob.saveToDat(fileName);
	}

	// n = 200, 7 stages, tcFF = 1.5, bottleneck last stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair((params.stgs - 1) / 2, (params.stgs - 1) / 2), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckMdl.dat";
		prob.saveToDat(fileName);
	}


	params.tcFlowFactor = 3.0;
	// n = 200, 7 stages, tcFF = 3.0, bottlenek 1st stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(0, 0), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckFirst.dat";
		prob.saveToDat(fileName);
	}

	// n = 200, 7 stages, tcFF = 3.0, bottleneck middle stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(params.stgs - 1, params.stgs - 1), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckLast.dat";
		prob.saveToDat(fileName);
	}

	// n = 200, 7 stages, tcFF = 3.0, bottleneck last stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair((params.stgs - 1) / 2, (params.stgs - 1) / 2), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckMdl.dat";
		prob.saveToDat(fileName);
	}

	// 15 stages
	params.stgs = 15;
	params.tcFlowFactor = 1.5;
	// n = 200, 15 stages, tcFF = 1.5, bottlenek 1st stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(0, 0), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckFirst.dat";
		prob.saveToDat(fileName);
	}

	// n = 200, 15 stages, tcFF = 1.5, bottleneck middle stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(params.stgs - 1, params.stgs - 1), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckLast.dat";
		prob.saveToDat(fileName);
	}

	// n = 200, 15 stages, tcFF = 1.5, bottleneck last stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair((params.stgs - 1) / 2, (params.stgs - 1) / 2), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckMdl.dat";
		prob.saveToDat(fileName);
	}


	params.tcFlowFactor = 3.0;
	// n = 200, 15 stages, tcFF = 3.0, bottlenek 1st stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(0, 0), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckFirst.dat";
		prob.saveToDat(fileName);
	}

	// n = 200, 15 stages, tcFF = 3.0, bottleneck middle stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(params.stgs - 1, params.stgs - 1), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckLast.dat";
		prob.saveToDat(fileName);
	}

	// n = 200, 15 stages, tcFF = 3.0, bottleneck last stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair((params.stgs - 1) / 2, (params.stgs - 1) / 2), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckMdl.dat";
		prob.saveToDat(fileName);
	}
}

void Problem::genInstancesTCB26small_Testing() {
	ProbParams params;
	params.omega = 9;
	params.F = 1;
	params.stgs = 5;
	params.n = 10;
	params.m_oIntervals = make_pair(1, 2);
	params.m_BIntervals = make_pair(1, 3);
	//params.m_BValues = vector<int>({ 3, 3, 3, 3, 3 });
	params.pInterval = make_pair(10, 25);
	params.tcScenario = 2;
	params.tcFlowFactor = 1.5;
	params.rInterval = make_pair(0, 0.75);
	params.sInterval = make_pair(1, 1);	// uniform job sizes
	params.wInterval = make_pair(1.0, 3.0);
	params.dueDateFF = make_pair(1.0, 1.3);


	int nmin_tc = params.stgs / 3.0;
	int nmax_tc = params.stgs - 1; // minimum number of timeconstraints = number of stages (-1)
	params.nTcInterval = make_pair(nmin_tc, nmax_tc);

	params.routes = vector<vector<int> >(params.F);
	for (int i = 0; i < params.F; ++i) {
		params.routes[i] = vector<int>(params.stgs);
		for (int o = 0; o < params.stgs; ++o) {
			params.routes[i][o] = o + 1;	// flow-shop
		}
	}

	int nInstances = 10;
	double bottleneckCriticality = 0.125;

	// n = 100, 7 stages, tcFF = 1.5, bottlenek 1st stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(0, 0), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckFirst.dat";
		prob.saveToDat(fileName);
	}

	// n = 100, 7 stages, tcFF = 1.5, bottleneck middle stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair(params.stgs - 1, params.stgs - 1), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckLast.dat";
		prob.saveToDat(fileName);
	}

	// n = 100, 9 stages, tcFF = 1.5, bottleneck last stage
	for (int i = 0; i < nInstances; ++i) {
		Problem prob = Problem(params, true);
		prob.configureBottleneck(make_pair((params.stgs - 1) / 2, (params.stgs - 1) / 2), bottleneckCriticality, params, true);

		stringstream tcFFstream;
		tcFFstream << fixed << setprecision(2) << params.tcFlowFactor;
		string fileName = "Inst_TCB2026_F" + to_string(params.F) + "m" + to_string(params.stgs) + "n" + to_string(params.n)
			+ "tcSc" + to_string(params.tcScenario) + "tcFF" + tcFFstream.str() + "_" + to_string(i + 1) + "btlNckMdl.dat";
		prob.saveToDat(fileName);
	}
}




			

