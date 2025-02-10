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

string Problem::getFilename() { return filename; }

int Problem::getN() const { return n; }
int Problem::getStgs() const { return stgs; }
int Problem::getF() const { return F; }

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

pair<int, int> Problem::_tokenizeTupel(string tupel) {
	pair<int, int> ret = pair<int, int>();
	size_t first = 1;					// "<"
	size_t last = tupel.length() - 3;	// ">:"
	size_t comma = tupel.find(",");		// ","
	ret.first = stoi(tupel.substr(first, comma - first));
	ret.second = stoi(tupel.substr(comma + 1, last));
	return ret;
}

unique_ptr<Schedule> Problem::getSchedule() const {
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

	return newSchedule;
}

