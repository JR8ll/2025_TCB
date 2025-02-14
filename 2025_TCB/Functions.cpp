#include <algorithm>
#include <ctime>
#include <iomanip>
#include <string>
#include <Windows.h>
#include <fstream>
#include <sstream>

#include "Functions.h"
#include "Solver_GA.h"
#include "Solver_MILP.h"
#include "Problem.h"

using namespace std;

void processCmd(int argc, char* argv[], int& iSolver, int& iTilimSeconds, bool& bConsole, GA_params& gaParams, DECOMPMILP_params& decompParams) {
	// argv[1] filename of problem instance to be solved
	if (argc > 1) {
		TCB::prob->loadFromDat(argv[1]);
	}
	else {
		TCB::prob->loadFromDat("debug_data.dat");
	}
	// argv[2] seed for pseudo random-number generator
	if (argc > 2) {
		TCB::seed = atoi(argv[2]);
		TCB::rng = mt19937(TCB::seed);
	}
	
	// argv[3] int describing the solving method to be used
	if (argc > 3) {
		iSolver = atoi(argv[3]);
	}
	
	// argv[4] time limit in seconds
	if (argc > 4) {
		iTilimSeconds = atoi(argv[4]);
	}
	
	// argv[5] console output on(=1)/off(=0)
	if (argc > 5) {
		if (atoi(argv[5]) == 1) bConsole = true;
	}
	
	// argv[6] filename of ga parameters
	if (argc > 6) {
		if (iSolver == ALG_BRKGALISTSCH || iSolver == ALG_BRKGALS2MILP) {
			try {
				loadGaParams(gaParams, argv[6]);
			}
			catch (...) {
				string warning = string("GA parameters file not found: ") + string(argv[1]) + string(" (default parameters initialized.)");
				gaParams = Solver_GA::getDefaultParams();
				TCB::logger.Log(Warning, warning);
			}
		}
	} else {
		string warning = string("GA parameters filename not provided (default parameters initialized.)");
		gaParams = Solver_GA::getDefaultParams();
		TCB::logger.Log(Warning, warning);
	}

	if (argc > 7) {
		if (iSolver == ALG_ITERATEDMILP || iSolver == ALG_BRKGALS2MILP) {
			try {
				loadDecompParams(decompParams, argv[7]);
			}
			catch (...) {
				string warning = string("DecompMILP parameters file not found: ") + string(argv[1]) + string(" (default parameters initialized.)");
				decompParams = Solver_MILP::getDefaultParams();
				TCB::logger.Log(Warning, warning);
			}
		}
	} else {
		string warning = string("DECOMPMILP parameters filename not provided (default parameters initialized.)");
		decompParams = Solver_MILP::getDefaultParams();
		TCB::logger.Log(Warning, warning);
	}
}
void writeSolutions(Schedule* solution, string solverName, string objectiveName, int prescribedTime, int usedTime, GA_params* gaParams, DECOMPMILP_params* decompParams) {
	bool success = CreateDirectory(L".\\results", NULL);
	string fullPath = ".\\results\\TCB_results.csv";
	ifstream checkFile(fullPath);
	bool fileExists = checkFile.good();
	checkFile.close();

	ofstream file(fullPath, ios::app);
	if (file.is_open()) {
		if (!fileExists) {
			// headings
			file << "Problem\t" << "Solver\t" << "Seed\t" << "Objective\t" << "ObjectiveValue\t" << "TimeLimit\t" << "TimeUsed\t" << "GA_params\t" << "MILPCP_params\t" << "CreatedOn" << endl;
		}
		file << solution->getProblem()->getFilename() << "\t" << solverName << "\t" << to_string(TCB::seed) << "\t" << objectiveName << "\t";
		file << to_string(solution->getTWT()) << "\t" << to_string(prescribedTime) << "\t" << to_string(usedTime) << "\t";

		ostringstream gaParamsString;
		ostringstream decompParamsString;
		if (gaParams != nullptr) {
			gaParamsString << gaParams->nPop << "|" << gaParams->pElt << "|" << gaParams->pRpM << "|" << gaParams->rhoe  << "|" << gaParams->K << "|" << gaParams->maxThreads;
		} else {
			gaParamsString << "n/a";
		}
		if (decompParams != nullptr) {
			decompParamsString << decompParams->nDash << "|" << decompParams->cplexTilim << "|" << decompParams->method << "|" << decompParams->initPrioRule << "|";
			decompParamsString << decompParams->kappaLow << "|" << decompParams->kappaHigh << "|" << decompParams->kappaStep;
		} else {
			decompParamsString << "n/a";
		}
		file << gaParamsString.str() << "\t" << decompParamsString.str() << "\t";
		time_t now = time(nullptr);
		tm* localTime = localtime(&now);
		ostringstream dateString;
		dateString << put_time(localTime, "%Y-%m-%d %H:%M:%S");
		file << dateString.str() << endl;
		file.close();
	} else {
		TCB::logger.Log(Error, "writeSolutions() results-file could not be opened");
	}
}

void sortJobsByD(vector<pJob>& unscheduledJobs) {
	sort(unscheduledJobs.begin(), unscheduledJobs.end(), compJobsByD);
}
void sortJobsByR(vector<pJob>& unscheduledJobs) {
	sort(unscheduledJobs.begin(), unscheduledJobs.end(), compJobsByR);
}
void sortJobsByGATC(vector<pJob>& unscheduledJobs, double t, double kappa) {
	double avgP = getAvgP(unscheduledJobs);
	sort(unscheduledJobs.begin(), unscheduledJobs.end(), CompJobsByGATC(avgP, t, kappa));
}

void sortJobsByRK(vector<pJob>& unscheduledJobs, const vector<double>& chr) {
	if (unscheduledJobs.size() != chr.size()) throw ExcSched("sortJobsByRK error different sizes of jobs and chromosome");

	vector<pJob> sortedJobs = vector<pJob>(chr.size());
	vector<pair<double, size_t>> ranking(chr.size());
	for (size_t j = 0; j < chr.size(); ++j) {
		ranking[j] = make_pair(chr[j], j);
	}
	sort(ranking.begin(), ranking.end());
	size_t pos = 0;
	for (auto rank = ranking.begin(); rank != ranking.end(); ++rank) {
		sortedJobs[pos] = move(unscheduledJobs[rank->second]);	// JENS assure that unscheduled size remains unchanged
		++pos;
	}

	for (size_t j = 0; j < chr.size(); ++j) {
		unscheduledJobs[j] = move(sortedJobs[j]);
	}
}

bool compJobsByD(const unique_ptr<Job>& a, const unique_ptr<Job>& b) {
	if (a->getD() == b->getD()) {
		return a->getId() < b->getId();
	}
	return a->getD() < b->getD();
}
bool compJobsByR(const unique_ptr<Job>& a, const unique_ptr<Job>& b) {
	if (a->getR() == b->getR()) {
		return a->getId() < b->getId();
	}
	return a->getR() < b->getR();
}

void shiftJobFromVecToVec(vector<pJob>& source, vector<pJob>& target, size_t sourceIdx) {
	if (sourceIdx >= source.size()) throw out_of_range("shiftJobFromVecToVec() out of range");
	target.push_back(move(source[sourceIdx]));
	source.erase(source.begin() + sourceIdx);
}

double getAvgP(const vector<pJob>& unscheduledJobs) {
	double totalP = 0.0;
	for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
		totalP += unscheduledJobs[j]->getTotalP();
	}
	return totalP / (double)unscheduledJobs.size();	// TODO check if to be divided by nJobs or nOps
}

void loadGaParams(GA_params& gaParams, string filename) {
	ifstream input(filename);
	string line;
	if (!input) throw ExcSched("loadGaParams file not found");
	while (getline(input, line)) {
		istringstream iss(line);
		string key;
		getline(iss, key, ':');
		if (key == "nPop") iss >> gaParams.nPop;
		else if (key == "pElt") iss >> gaParams.pElt;
		else if (key == "pRpM") iss >> gaParams.pRpM;
		else if (key == "rhoe") iss >> gaParams.rhoe;
		else if (key == "K") iss >> gaParams.K;
		else if (key == "maxThreads") iss >> gaParams.maxThreads;
	}
	input.close();
}
void loadDecompParams(DECOMPMILP_params& decompParams, string filename) {
	ifstream input(filename);
	string line;
	if (!input) throw ExcSched("loadDecompParams file not found");

	double low = -1;
	double high = -1;
	double step = -1;

	while (getline(input, line)) {
		istringstream iss(line);
		string key;
		getline(iss, key, ':');
		if (key == "nDash") iss >> decompParams.nDash;
		else if (key == "cplexTilim") iss >> decompParams.cplexTilim;
		else if (key == "method") iss >> decompParams.method;
		else if (key == "initPrioRule") iss >> decompParams.initPrioRule;
		else if (key == "kappasLow") iss >> decompParams.kappaLow;
		else if (key == "kappasHigh") iss >> decompParams.kappaHigh;
		else if (key == "kappasStep") iss >> decompParams.kappaStep;
	}
	input.close();
}

double getObjectiveTWT(const Schedule* sched) {
	return sched->getTWT();
}

vector<double> getDoubleGrid(double low, double high, double step) {
	vector<double> grid = vector<double>();
	while (low <= high + TCB::precision) {
		grid.push_back(low);
		low = low + step;
	}
	return grid;
}

