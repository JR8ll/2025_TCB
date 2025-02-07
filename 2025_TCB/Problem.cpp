#include "Problem.h"

#include<fstream>

#include "Functions.h"

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
	for (int i = 0; i < stgs; ++i) {
		input >> dummy;
		m_o.push_back(stoi(dummy));
	}
	input >> dummy;	// "];"

	m_B = vector<int>(stages_bSize);

	rm = vector<vector<double> >(stgs);
	for (int i = 0; i < stgs; ++i) {
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
	for (int i = 0; i < F; ++i) {
		B_io[i] = vector<set<int> >(stgs);
		for (int j = 0; j < stgs; ++j) {
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
	for (int i = 0; i < F; ++i) {
		B_iob[i] = vector<vector<int> >(stgs);
		for (int j = 0; j < stgs; ++j) {
			B_iob[i][j] = vector<int>(n);
			for (int k = 0; k < n; ++k) {
				input >> dummy;
				B_iob[i][j][k] = stoi(dummy);
			}
			input >> dummy;		// "][" / "]][[" / "]]];
		}
	}

	input >> dummy; // "S_iob=["
	S_iob = vector<vector<vector<double> > >(F);
	for (int i = 0; i < F; ++i) {
		S_iob[i] = vector<vector<double> >(stgs);
		for (int j = 0; j < stgs; ++j) {
			S_iob[i][j] = vector<double>(n);
			for (int k = 0; k < n; ++k) {
				input >> dummy;
				S_iob[i][j][k] = stod(dummy);
			}
			input >> dummy;		// "][" / "]][[" / "]]];
		}
	}

	input >> dummy;	// "B=["
	for (int i = 0; i < stages_b.size(); ++i) {
		input >> dummy;
		m_B[i] = stoi(dummy);
	}
	input >> dummy; // "];"

	input >> dummy;	// "d=["
	for (int i = 0; i < n; ++i) {
		input >> dummy;
		jobs_d.push_back(stod(dummy));
	}
	input >> dummy; // "];"

	input >> dummy;	// "r=["
	for (int i = 0; i < n; ++i) {
		input >> dummy;
		jobs_r.push_back(stod(dummy));
	}
	input >> dummy; // "];"

	input >> dummy;	// "w=["
	for (int i = 0; i < n; ++i) {
		input >> dummy;
		jobs_w.push_back(stod(dummy));
	}
	input >> dummy; // "];"

	input >> dummy;	// "f=["
	for (int i = 0; i < n; ++i) {
		input >> dummy;
		jobs_f.push_back(stoi(dummy));
	}
	input >> dummy; // "];"

	input >> dummy;	// "s=["
	for (int i = 0; i < n; ++i) {
		input >> dummy;
		jobs_s.push_back(stoi(dummy));
	}
	input >> dummy; // "];"

	input >> dummy;	// "rts=[{
	routes = vector<vector<int> >(F);
	for (int i = 0; i < F; ++i) {
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

	pTimes = vector < vector<double> >(F);
	for (int i = 0; i < F; ++i) {
		pTimes[i] = vector<double>(routes[i].size());
	}

	input >> dummy; // "p=[["
	for (int p = 0; p < F; ++p) {
		for (int s = 0; s < routes[p].size(); ++s) {
			input >> dummy;
			pTimes[p][s] = stod(dummy);
		}
		input >> dummy; // "]["/"]];"
	}

	tc = vector<vector<vector<double> > >(F);
	for (int p = 0; p < F; ++p) {
		tc[p] = vector<vector<double> >(routes[p].size());
		for (int s = 0; s < routes[p].size(); ++s) {
			tc[p][s] = vector<double>(routes[p].size());
		}
	}

	input >> dummy;	// "tc=[[[..."
	for (int p = 0; p < F; ++p) {
		// 1st dim: Product
		for (int s1 = 0; s1 < routes[p].size(); ++s1) {
			// 2nd dim: Stage
			for (int s2 = 0; s2 < routes[p].size(); ++s2) {
				// 3rd dim: Stage
				input >> dummy;
				tc[p][s1][s2] = stod(dummy);
			}
			input >> dummy;		// "]["
		}
	}

	// instantiate entities Jobs + Ops including time constraints
	jobs.resize(n);
	for (int j = 0; j < n; ++j) {
		jobs[j].setD(jobs_d[j]);
		jobs[j].setR(jobs_r[j]);
		jobs[j].setW(jobs_w[j]);
		jobs[j].setF(jobs_f[j]);
		jobs[j].setS(jobs_s[j]);
		int nSteps = routes[jobs_f[j] - 1].size();
		for (int o = 0; o < nSteps; ++o) {
			shared_ptr<Operation> newOp = make_shared<Operation>(Operation(jobs[j].getId(), (o + 1), routes[jobs_f[j] - 1][o], jobs_f[j], pTimes[jobs_f[j] - 1][o], jobs[j]));
			jobs[j].addOp(newOp);
		}
		for (int o = 0; o < nSteps - 1; ++o) {
			jobs[j].at(o)->setSucc(*jobs[j].at(o + 1));
		}

		// add time constraints
		for (int o = 0; o < nSteps; ++o) {
			pOp opPtr = jobs[j].at(o);
			int pIdx = jobs_f[j] - 1;
			for (int s = 0; s < tc[pIdx].size(); ++s) {
				double time = tc[pIdx][s][o];
				if (time != 999999) {
					opPtr->addTcMax(s, time);
				}

				double timeFwd = tc[pIdx][o][s];
				if (timeFwd != 999999) {
					opPtr->addTcMaxFwd(s, timeFwd);
				}
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

