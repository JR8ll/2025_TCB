
#include <algorithm>
#include <thread>
#include <iostream>
#include <numeric>
#include <random>
#include "Solver_ILS.h"
#include "Schedule.h"


using namespace std;

Solver_ILS::Solver_ILS() : schedParams(nullptr), params(nullptr) {}
Solver_ILS::Solver_ILS(Sched_params& schedParams, ILS_params& params) {
    this->schedParams = &schedParams;
    this->params = &params;
}

double Solver_ILS::solveILS(Schedule& sched, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, double pWait) {
    auto start = chrono::high_resolution_clock::now();
    chrono::seconds usedTime;
    chrono::time_point<chrono::high_resolution_clock> stop;

    double bestTWT = DBL_MAX;
    unique_ptr<Schedule> bestSched = sched.clone();
    
    params->multiStartIterations = 0;
    
    do {// MULTISTART-LOOP
        // INITIALIZE
        unique_ptr<Schedule> tempSched = sched.clone();
        (tempSched.get()->*init)(rule, *schedParams);

            params->ilsIterations.push_back(0);
        do {// ILS LOOP
            // LOCAL SEARCH
            
           // DEBUGGING
           //cout << "ILS iteration " << params->multiStartIterations + 1 << "." << params->ilsIterations[params->multiStartIterations] + 1 << " TWT: " << bestTWT << endl;
           /*TCB::logger.Log(Info, to_string(params->ilsIterations[params->multiStartIterations] + 1));
           if (params->multiStartIterations + 1 == 1 && params->ilsIterations[params->multiStartIterations] + 1 == 28) {
                tempSched->saveJsonFactory("debugging");
                int debugger = 666;
            }*/
            
            // [JR-2026-Jan-12] wrapped local search in do-while-loop
            bool bLeftShiftApplied = false;
            bool bJobSwapApplied = false;
            bool bBatchConsolidationApplied = false;
            do {  
                bLeftShiftApplied = tempSched->localSearchJobLeftShifting(&sortJobsRandomly, params->applyBestFit);
                bJobSwapApplied = tempSched->localSearchJobSwapping(&sortJobsRandomly, params->applyBestFit);
                bBatchConsolidationApplied = tempSched->localSearchBatchConsolidation(params->applyBestFit);
            } while (bLeftShiftApplied || bJobSwapApplied || bBatchConsolidationApplied);




            double tempTWT = tempSched->getTWT();
            if (tempTWT < bestTWT) {
                bestTWT = tempTWT;
                bestSched = tempSched->clone();
                stop = chrono::high_resolution_clock::now();
                usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
                params->bestAfterSeconds = usedTime.count();
            }
            // PERTURBATION
            uniform_real_distribution<> perturbDistrib(0, 1);
            for (size_t i = 0; i < params->nPerturbationSteps; ++i) {
                double perturbChoice = perturbDistrib(TCB::rng);
                if (perturbChoice < 0.5) {
                    //TCB::logger.Log(Info, "perturbRandomJobSwap");
                    tempSched->perturbRandomJobSwap();
                }
                else {
                    //TCB::logger.Log(Info, "perturbRandomJobRightShifting");
                    tempSched->perturbRandomJobRightShifting();  
                }
            }

            ++params->ilsIterations[params->multiStartIterations];
            stop = chrono::high_resolution_clock::now();
            usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
        } while (usedTime.count() < ((double)iTilimSeconds / (double)params->nStarts) * (params->multiStartIterations + 1));   // MULTISTART
        stop = chrono::high_resolution_clock::now();
        usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
        ++params->multiStartIterations;
    } while (usedTime.count() < iTilimSeconds);

    sched._reconstruct(bestSched.get());
    return bestTWT;
}

double Solver_ILS::solveILSparallelized(Schedule& sched, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, double pWait) {
    unsigned int nCores = thread::hardware_concurrency();

    vector<DWORD> pCores = GetPCoreIndices();
    //cout << "Using " << pCores.size() << " cores..." << endl;

    auto start = chrono::high_resolution_clock::now();
    chrono::seconds usedTime;
    chrono::time_point<chrono::high_resolution_clock> stop;


    vector<ILS_Thread> ILS_threads(nCores);
    vector<thread> threads;
    for (unsigned int core = 0; core < pCores.size(); ++core) {    // [JR-2026-Jan-22] replaced nCores with pCores.size()
        DWORD coreIndex = pCores[core];
        ILS_threads[core].bestSched = sched.clone();
        threads.emplace_back(workerILS, coreIndex, move(ILS_threads[core].bestSched), schedParams, params, init, rule, iTilimSeconds, start, &ILS_threads[core], pWait);
    }

    for (auto& t : threads) t.join();

    int globalBestIdx = 0;
    double globalBestTWT = ILS_threads[0].bestTWT;
    for (int i = 0; i < pCores.size(); ++i) {  // [JR-2026-Jan-22] replaced nCores with pCores.size()
        if (ILS_threads[i].bestTWT < globalBestTWT) {
            globalBestTWT = ILS_threads[i].bestTWT;
            globalBestIdx = i;
            params->bestAfterSeconds = ILS_threads[i].bestAfterSeconds;
            params->ilsIterations = ILS_threads[i].ilsIterations;
            params->multiStartIterations = ILS_threads[i].multiStartIterations;
        }
    }

    sched._reconstruct(ILS_threads[globalBestIdx].bestSched.get());
    return globalBestTWT;
}

double Solver_ILS::solveILSonJobSequence(Schedule& sched, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, double pWait) {
    cout << "Solver_ILS::solveILSonJobSequence(...) not yet implemented." << endl;
    
    return -1.0;
}

void Solver_ILS::workerILS(DWORD coreIndex, unique_ptr<Schedule>& sched, Sched_params* schedParams, ILS_params* ilsParams, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, chrono::time_point<chrono::high_resolution_clock> start, ILS_Thread* localBest, double pWait) {
    DWORD_PTR mask = 1ULL << coreIndex;
    SetThreadAffinityMask(GetCurrentThread(), mask);
    
    double bestTWT = DBL_MAX;
    unique_ptr<Schedule> bestSched = sched->clone();
    chrono::seconds usedTime;
    chrono::time_point<chrono::high_resolution_clock> stop;

    localBest->multiStartIterations = 0;

    do {// MULTISTART-LOOP
        // INITIALIZE
        unique_ptr<Schedule> tempSched = sched->clone();
        (tempSched.get()->*init)(rule, *schedParams);

        // DEBUGGING
        //TCB::logger.Log(Info, "Thread started.");

        localBest->ilsIterations.push_back(0);
        do {// ILS LOOP
            // LOCAL SEARCH
            bool bLeftShiftApplied = false;
            bool bJobSwapApplied = false;
            bool bBatchConsolidationApplied = false;
            do {
                bLeftShiftApplied = tempSched->localSearchJobLeftShifting(&sortJobsRandomly, ilsParams->applyBestFit);
                bJobSwapApplied = tempSched->localSearchJobSwapping(&sortJobsRandomly, ilsParams->applyBestFit);
                bBatchConsolidationApplied = tempSched->localSearchBatchConsolidation(ilsParams->applyBestFit);
            } while (bLeftShiftApplied || bJobSwapApplied || bBatchConsolidationApplied);

            double tempTWT = tempSched->getTWT();
            if (tempTWT < localBest->bestTWT) {
                localBest->bestTWT = tempTWT;
                localBest->bestSched = tempSched->clone();
                stop = chrono::high_resolution_clock::now();
                usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
                localBest->bestAfterSeconds = usedTime.count();
            }
            // PERTURBATION
            uniform_real_distribution<> perturbDistrib(0, 1);
            for (size_t i = 0; i < ilsParams->nPerturbationSteps; ++i) {
                double perturbChoice = perturbDistrib(TCB::rng);
                if (perturbChoice < 0.5) {
                    tempSched->perturbRandomJobSwap();
                }
                else {
                    tempSched->perturbRandomJobRightShifting();
                }
            }

            ++localBest->ilsIterations[localBest->multiStartIterations];
            stop = chrono::high_resolution_clock::now();
            usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
        } while (usedTime.count() < ((double)iTilimSeconds / (double)ilsParams->nStarts) * (ilsParams->multiStartIterations + 1));   // MULTISTART

        //TCB::logger.Log(Info, "Thread ended.");

        stop = chrono::high_resolution_clock::now();
        usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
        ++localBest->multiStartIterations;
    } while (usedTime.count() < iTilimSeconds);
}

ILS_params Solver_ILS::getDefaultParams() {
    ILS_params ilsParams = ILS_params();
    ilsParams.nStarts = 1;
    ilsParams.nPerturbationSteps = 5;
    ilsParams.applyBestFit = true;
    ilsParams.randomizedLocalSearchSequence = false;
    ilsParams.multiStartIterations = 0;
    ilsParams.ilsIterations = vector<size_t>();
    return ilsParams;
}

std::vector<DWORD> Solver_ILS::GetPCoreIndices() {
    cout << "Optimized for 12th Gen Intel(R) Core(TM) i7-12700: using 8 performance cores..." << endl;
    return { 0, 1, 2, 3, 4, 5, 6, 7 };
    
    SYSTEM_INFO si;
    GetSystemInfo(&si);

    // Typisch: P-Cores sind die ersten logischen Kerne
    DWORD numPCores = si.dwNumberOfProcessors / 2;
    std::vector<DWORD> pcores;

    for (DWORD i = 0; i < numPCores; i++) {
        pcores.push_back(i);
    }
    return pcores;
}


Solver_Sequence_ILS::Solver_Sequence_ILS(Schedule* schedule, Sched_params* schedParameters, ILS_params* ilsParameters) : decoder(schedule, schedParameters, nullptr), chr(schedule->getN()) {
    schedParams = schedParameters;
    params = ilsParameters;
}
Solver_Sequence_ILS::~Solver_Sequence_ILS(){}

double Solver_Sequence_ILS::solveILS(Schedule& sched, int iTilimSeconds) {
    auto start = chrono::high_resolution_clock::now();
    chrono::seconds usedTime;
    chrono::time_point<chrono::high_resolution_clock> stop;

    double bestTWT = DBL_MAX;

    params->multiStartIterations = 0;

    do {    // MULTISTART-LOOP
            // INITIALIZE
            // DEFAULT: random initializer (TODO: other initializers to be defined)
        initRandomPermutation();

        params->ilsIterations.push_back(0);
        do {// ILS LOOP
            // LOCAL SEARCH

            cout << "ILS iteration " << params->multiStartIterations + 1 << "." << params->ilsIterations[params->multiStartIterations] + 1 << "currTWT: " << currentTWT << ", bestTWT" <<  bestTWT << endl;

            bool bJobInsertApplied = false;
            bool bJobSwapApplied = false;
            do {
                bJobInsertApplied = localSearchInsertJob(params->applyBestFit);
                bJobSwapApplied = localSearchSwapJob(params->applyBestFit);
            } while (bJobInsertApplied || bJobSwapApplied);

            if (currentTWT < bestTWT) {
                bestTWT = currentTWT;
                bestChr = chr;
                stop = chrono::high_resolution_clock::now();
                usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
                params->bestAfterSeconds = usedTime.count();
            }
            
            // PERTURBATION
            uniform_real_distribution<> perturbDistrib(0, 1);
            for (size_t i = 0; i < params->nPerturbationSteps; ++i) {
                double perturbChoice = perturbDistrib(TCB::rng);
                if (perturbChoice < 1.0 / 3.0) {
                    //TCB::logger.Log(Info, "perturbRandomJobSwap");
                    perturbJobInsert();
                }
                else if(perturbChoice < (1.0 / 3.0) * 2.0) {
                    //TCB::logger.Log(Info, "perturbRandomJobRightShifting");
                    perturbJobSwap();
                }
                else {
                    perturbRandomKey();
                }
            }

            ++params->ilsIterations[params->multiStartIterations];
            stop = chrono::high_resolution_clock::now();
            usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
        } while (usedTime.count() < ((double)iTilimSeconds / (double)params->nStarts) * (params->multiStartIterations + 1));   // MULTISTART
        stop = chrono::high_resolution_clock::now();
        usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
        ++params->multiStartIterations;
    } while (usedTime.count() < iTilimSeconds);

    decoder.formSchedule(bestChr);
    return bestTWT;
}

void Solver_Sequence_ILS::insertJob(vector<double>& perm, size_t jobIdx, size_t posIdx) {
    if (jobIdx < 0 || jobIdx >= perm.size() || posIdx < 0 || posIdx >= perm.size()) throw out_of_range("Solver_Sequence_ILS::insertJob(...) out of range");
    if (jobIdx < posIdx) {
        rotate(perm.rbegin() + (perm.size() - 1 - posIdx), perm.rbegin() + (perm.size() - 1 - jobIdx), perm.rend());
    }
    if (posIdx < jobIdx) {
        rotate(perm.begin() + posIdx, perm.begin() + jobIdx, perm.begin() + jobIdx + 1);
    }
}

void Solver_Sequence_ILS::swapJob(vector<double>& perm, size_t firstIdx, size_t secondIdx){
    if (firstIdx < 0 || firstIdx >= perm.size() || secondIdx < 0 || secondIdx >= perm.size()) throw out_of_range("Solver_Sequence_ILS::swapJob(...) out of range");
    if (firstIdx != secondIdx) {
        swap(perm[firstIdx], perm[secondIdx]);
    }
}

double Solver_Sequence_ILS::evaluateJobInsert(size_t jobIdx, size_t posIdx) {
    vector<double> tempChr = chr;
    insertJob(tempChr, jobIdx, posIdx);
    return currentTWT - decoder.decode(tempChr);
}

double Solver_Sequence_ILS::evaluateJobSwap(size_t firstIdx, size_t secondIdx) {
    vector<double> tempChr = chr;
    swapJob(tempChr, firstIdx, secondIdx);
    return currentTWT - decoder.decode(tempChr);
}

bool Solver_Sequence_ILS::localSearchInsertJob(bool bestFit) {
    // DEFAULT: random order (TODO: other sorting orders to be defined)
    vector<size_t> indices(chr.size());
    iota(indices.begin(), indices.end(), 0);
    shuffle(indices.begin(), indices.end(), TCB::rng);

    double bestImprovement = 0;
    size_t bestJobIdx = 0;
    size_t bestPosIdx = 0;
    bool bImproved = true;
    bool bEventuallyImproved = false;
    while (bImproved) {
        bImproved = false;

        for (size_t i = 0; i < indices.size(); ++i) {
            for (size_t j = 0; j < indices.size(); ++j) {
                if (i != j) {
                    double tempImprovement = evaluateJobInsert(indices[i], indices[j]);
                    if (tempImprovement > bestImprovement) {
                        bestImprovement = tempImprovement;
                        bestJobIdx = indices[i];
                        bestPosIdx = indices[j];
                        if (!bestFit && bestImprovement > 0) {     // FIRST FIT
                            insertJob(chr, bestJobIdx, bestPosIdx);
                            currentTWT = currentTWT - bestImprovement;  // decoding unnecessary 
                            bImproved = true;
                            bEventuallyImproved = true;
                            bestImprovement = 0;
                            break;
                        }
                    
                    }
                }
                if (bImproved) {
                    break;  // FIRST FIT: continue in while loop
                }
            }
            if (bImproved) {
                break;  // FIRST FIT: continue in while loop
            }
        }

        if (bestFit && bestImprovement > 0) {       // BEST FIT
            insertJob(chr, bestJobIdx, bestPosIdx);
            currentTWT = currentTWT - bestImprovement;
            bImproved = true;
            bEventuallyImproved = true;
            bestImprovement = 0;
        }
    }

    return bEventuallyImproved;
}
bool Solver_Sequence_ILS::localSearchSwapJob(bool bestFit) {
    // DEFAULT: random order (TODO: other sorting orders to be defined)
    vector<size_t> indices(chr.size());
    iota(indices.begin(), indices.end(), 0);
    shuffle(indices.begin(), indices.end(), TCB::rng);

    double bestImprovement = 0;
    size_t bestFirstIdx = 0;
    size_t bestSecondIdx = 0;
    bool bImproved = true;
    bool bEventuallyImproved = false;
    while (bImproved) {
        bImproved = false;
        for (size_t i = 0; i < indices.size(); ++i) {
            for (size_t j = 0; j < indices.size(); ++j) {
                if (i != j) {
                    double tempImprovement = evaluateJobSwap(indices[i], indices[j]);
                    if (tempImprovement > bestImprovement) {
                        bestImprovement = tempImprovement;
                        bestFirstIdx = indices[i];
                        bestSecondIdx = indices[j];
                        if (!bestFit && bestImprovement > 0) {     // FIRST FIT
                            insertJob(chr, bestFirstIdx, bestSecondIdx);
                            currentTWT = currentTWT - bestImprovement;  // decoding unnecessary 
                            bImproved = true;
                            bEventuallyImproved = true;
                            bestImprovement = 0;
                            break;
                        }

                    }
                }
                if (bImproved) {
                    break;  // FIRST FIT: continue in while loop
                }
            }
            if (bImproved) {
                break;  // FIRST FIT: continue in while loop
            }
        }

        if (bestFit && bestImprovement > 0) {       // BEST FIT
            insertJob(chr, bestFirstIdx, bestSecondIdx);
            currentTWT = currentTWT - bestImprovement;
            bestImprovement = 0;
            bEventuallyImproved = true;
            bImproved = true;
        }
    }

    return bEventuallyImproved;
}
void Solver_Sequence_ILS::perturbJobInsert(){
    uniform_int_distribution<> idxDistrib(0, chr.size() - 1);
    size_t i = idxDistrib(TCB::rng);
    size_t j = i;
    while (j == i) {
        j = idxDistrib(TCB::rng);
    }
    insertJob(chr, i, j);
    currentTWT = decoder.decode(chr);
}
void Solver_Sequence_ILS::perturbJobSwap() {
    uniform_int_distribution<> idxDistrib(0, chr.size() - 1);
    size_t i = idxDistrib(TCB::rng);
    size_t j = i;
    while (j == i) {
        j = idxDistrib(TCB::rng);
    }
    swapJob(chr, i, j);
    currentTWT = decoder.decode(chr);
}
void Solver_Sequence_ILS::perturbRandomKey() {
    uniform_real_distribution<> rkDistrib(0.0, 1.0);
    uniform_int_distribution<> idxDistrib(0, chr.size() - 1);
    size_t idx = idxDistrib(TCB::rng);
    double newRkValue = rkDistrib(TCB::rng);
    chr[idx] = newRkValue;
    currentTWT = decoder.decode(chr);
}
void Solver_Sequence_ILS::initRandomPermutation(){
    uniform_real_distribution<> rkDistrib(0.0, 1.0);
    generate(chr.begin(), chr.end(), [&] { return rkDistrib(TCB::rng); });
    currentTWT = decoder.decode(chr);
}



