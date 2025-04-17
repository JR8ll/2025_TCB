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

void processCmd(int argc, char* argv[], int& iSolver, int& iTilimSeconds, bool& bConsole, Sched_params& schedParams, GA_params& gaParams, DECOMPMILP_params& decompParams) {
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

	// argv[6] filename of scheduling parameters
	if (argc > 6) {
		try {
			loadSchedParams(schedParams, argv[6]);
		}
		catch (...) {
			string warning = string("Scheduling parameters file not found: ") + string(argv[6]) + string(" (default parameters initialized.)");
			schedParams = getDefaultParams();
			TCB::logger.Log(Warning, warning);
		}
	}
	else {
		string warning = string("Scheduling parameters filename not provided (default parameters initialized.)");
		schedParams = getDefaultParams();
		TCB::logger.Log(Warning, warning);
	}
	
	// argv[7] filename of ga parameters
	if (argc > 7) {
		if (iSolver == ALG_BRKGALISTSCH || iSolver == ALG_BRKGALS2MILP) {
			try {
				loadGaParams(gaParams, argv[7]);
			}
			catch (...) {
				string warning = string("GA parameters file not found: ") + string(argv[7]) + string(" (default parameters initialized.)");
				gaParams = Solver_GA::getDefaultParams();
				TCB::logger.Log(Warning, warning);
			}
		}
	} else {
		string warning = string("GA parameters filename not provided (default parameters initialized.)");
		gaParams = Solver_GA::getDefaultParams();
		TCB::logger.Log(Warning, warning);
	}

	// argv[8] filename of decomp MILP parameters
	if (argc > 8) {
		if (iSolver == ALG_ITERATEDMILP || iSolver == ALG_BRKGALS2MILP || iSolver == ALG_ITMILPLSHIFT) {
			try {
				loadDecompParams(decompParams, argv[8]);
			}
			catch (...) {
				string warning = string("DecompMILP parameters file not found: ") + string(argv[8]) + string(" (default parameters initialized.)");
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
void writeSolutions(Schedule* solution, int solverType, string solverName, string objectiveName, int prescribedTime, int usedTime, Sched_params* schedParams, GA_params* gaParams, DECOMPMILP_params* decompParams) {
	bool success = CreateDirectory(L".\\results", NULL);
	string fullPath = ".\\results\\TCB_results.csv";
	ifstream checkFile(fullPath);
	bool fileExists = checkFile.good();
	checkFile.close();

	ofstream file(fullPath, ios::app);
	if (file.is_open()) {
		if (!fileExists) {
			// headings
			file << "Problem\t" << "Solver\t" << "Seed\t" << "Objective\t" << "ObjectiveValue\t" << "TimeLimit\t" << "TimeUsed\t" << "SchedParams\t" << "GA_params\t" << "MILPCP_params\t" << "miscReporting\t" << "CreatedOn" << endl;
		}
		
		double objectiveValue = 999999;	// invalid
		if (solution->isValid()) {
			objectiveValue = solution->getTWT();
		}

		file << solution->getProblem()->getFilename() << "\t" << solverName << "\t" << to_string(TCB::seed) << "\t" << objectiveName << "\t";
		file << to_string(objectiveValue) << "\t" << to_string(prescribedTime) << "\t" << to_string(usedTime) << "\t";

		ostringstream schedParamsString;
		ostringstream gaParamsString;
		ostringstream decompParamsString;
		if (schedParams != nullptr) {
			schedParamsString << "kap{" << schedParams->kappaLow << "/" << schedParams->kappaHigh << "/" << schedParams->kappaStep << "}pW{" << schedParams->pWaitLow << "/" << schedParams->pWaitHigh << "/" << schedParams->pWaitStep << "}";
		} else {
			schedParamsString << "n/a";
		}
		if (gaParams != nullptr) {
			gaParamsString << gaParams->nPop << "|" << gaParams->pElt << "|" << gaParams->pRpM << "|" << gaParams->rhoe  << "|" << gaParams->K << "|" << gaParams->maxThreads;
		} else {
			gaParamsString << "n/a";
		}
		if (decompParams != nullptr) {
			decompParamsString << decompParams->nDash << "|" << decompParams->cplexTilim << "|" << decompParams->method << "|" << decompParams->initPrioRule;
		} else {
			decompParamsString << "n/a";
		}
		file << schedParamsString.str() << "\t" << gaParamsString.str() << "\t" << decompParamsString.str() << "\t";
		time_t now = time(nullptr);
		tm* localTime = localtime(&now);
		ostringstream dateString;
		dateString << put_time(localTime, "%Y-%m-%d %H:%M:%S");

		// misc reporting
		if (solverType == ALG_BRKGALISTSCH || solverType == ALG_BRKGALS2MILP) {
			file << "nGen=" << gaParams->iterations << "\t";
		}
		else if (solverType == ALG_LISTSCHEDATC || solverType == ALG_ITMILPLSHIFT) {
			file << "leftShImpr=" << schedParams->leftShiftImprovement << "\t";
		} else {
			file << "n/a\t";
		}

		file << dateString.str() << endl;
		file.close();
	} else {
		TCB::logger.Log(Error, "writeSolutions() results-file could not be opened");
	}
}

void sortJobsByC(std::vector<pJob>& jobs) {
	sort(jobs.begin(), jobs.end(), compJobsByC);
}
void sortJobsByStart(std::vector<pJob>& jobs) {
	sort(jobs.begin(), jobs.end(), compJobsByStart);
}
void sortJobsByWaitingTimeDecr(std::vector<pJob>& jobs) {
	sort(jobs.begin(), jobs.end(), compJobsByWaitingTimeDecr);
}
void sortJobsByD(vector<pJob>& jobs) {
	sort(jobs.begin(), jobs.end(), compJobsByD);
}
void sortJobsByR(vector<pJob>& jobs) {
	sort(jobs.begin(), jobs.end(), compJobsByR);
}
void sortJobsByGATC(vector<pJob>& jobs, double t, double kappa) {
	double avgP = getAvgP(jobs);
	sort(jobs.begin(), jobs.end(), CompJobsByGATC(avgP, t, kappa));
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

void sortJobsDebugging(std::vector<pJob>& jobs) {
	sort(jobs.begin(), jobs.end(), compJobsDebugging);
}

bool compJobsByC(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b) {
	double cA = a->getC();
	double cB = b->getC();
	if (cA == cB) {
		return a->getId() < b->getId();
	}
	return cA < cB;
}
bool compJobsByStart(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b) {
	double startA = a->getStart();
	double startB = b->getStart();
	if (startA == startB) {
		return a->getId() < b->getId();
	}
	return startA < startB;
}
bool compJobsByWaitingTimeDecr(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b) {
	double waitA = a->getWait();
	double waitB = b->getWait();
	if (waitA == waitB) {
		return a->getId() < b->getId();
	}
	return waitA > waitB;
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

bool compJobsDebugging(const std::unique_ptr<Job>& a, const std::unique_ptr<Job>& b) {
	switch (a->getId()) {
	case 1:
		return b->getId() != 8 && b->getId() != 6 && b->getId() != 5 && b->getId() != 9 && b->getId() != 2 && b->getId() != 7 && b->getId() != 3;
	case 2:
		return b->getId() != 8 && b->getId() != 6 && b->getId() != 5 && b->getId() != 9;
	case 3:
		return b->getId() != 8 && b->getId() != 6 && b->getId() != 5 && b->getId() != 9 && b->getId() != 2 && b->getId() != 7;
	case 4:
		return false;
	case 5:
		return b->getId() != 8 && b->getId() != 6;
	case 6:
		return true;
	case 7:
		return b->getId() != 8 && b->getId() != 6 && b->getId() != 5 && b->getId() != 9 && b->getId() != 2;
	case 8:
		return b->getId() != 6;
	case 9:
		return b->getId() != 8 && b->getId() != 6 && b->getId() != 5;
	}
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

void loadSchedParams(Sched_params& schedParams, std::string filename) {
	ifstream input(filename);
	string line;
	if (!input) throw ExcSched("loadSchedParams file not found");
	while (getline(input, line)) {
		istringstream iss(line);
		string key;
		getline(iss, key, ':');
		if (key == "pWaitLow") iss >> schedParams.pWaitLow;
		else if (key == "pWaitHigh") iss >> schedParams.pWaitHigh;
		else if (key == "pWaitStep") iss >> schedParams.pWaitStep;
		else if (key == "kappaLow") iss >> schedParams.kappaLow;
		else if (key == "kappaHigh") iss >> schedParams.kappaHigh;
		else if (key == "kappaStep") iss >> schedParams.kappaStep;
	}
	input.close();
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
	}
	input.close();
}

Sched_params getDefaultParams() {
	Sched_params defaultParams = Sched_params();
	defaultParams.pWaitLow = 0.0;
	defaultParams.pWaitHigh = 0.0;
	defaultParams.pWaitStep = 0.0;

	defaultParams.kappaLow = 0.1;
	defaultParams.kappaHigh = 5.0;
	defaultParams.kappaStep = 0.1;

	return Sched_params();
}

double getObjectiveTWT(const Schedule* sched) {
	return sched->getTWT();
}

vector<double> getDoubleGrid(double low, double high, double step) {
	vector<double> grid = vector<double>();
	if (low >= high) {
		return { low };
	}
	while (low <= high + TCB::precision) {
		grid.push_back(low);
		low = low + step;
	}
	return grid;
}

string extractFileName(const string& fullPath) {
	size_t lastSlash = fullPath.find_last_of("/\\");
	if (lastSlash != std::string::npos) {
		return fullPath.substr(lastSlash + 1);
	}
	return fullPath; // Falls kein Slash gefunden wurde, geben wir den gesamten String zurück
}
void replaceWindowsSpecialCharsWithUnderscore(string& input) {
	replace(input.begin(), input.end(), '*', '_');
	replace(input.begin(), input.end(), '"', '_');
	replace(input.begin(), input.end(), '/', '_');
	replace(input.begin(), input.end(), '\\', '_');
	replace(input.begin(), input.end(), '<', '_');
	replace(input.begin(), input.end(), '>', '_');
	replace(input.begin(), input.end(), ':', '_');
	replace(input.begin(), input.end(), '|', '_');
	replace(input.begin(), input.end(), '?', '_');
}

