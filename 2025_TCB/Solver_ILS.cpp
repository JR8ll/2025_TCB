
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
        threads.emplace_back(workILS, coreIndex, move(ILS_threads[core].bestSched), schedParams, params, init, rule, iTilimSeconds, start, &ILS_threads[core], pWait);
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

double Solver_ILS::solveILSparallelized(Schedule& sched, initializerRK<pJob> init, const std::vector<double>& randomKeys, int iTilimSeconds, double pWait) {
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
        threads.emplace_back(workILSrk, coreIndex, move(ILS_threads[core].bestSched), schedParams, params, init, randomKeys, iTilimSeconds, start, &ILS_threads[core], pWait);
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

void Solver_ILS::workILS(DWORD coreIndex, unique_ptr<Schedule>& sched, Sched_params* schedParams, ILS_params* ilsParams, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, chrono::time_point<chrono::high_resolution_clock> start, ILS_Thread* localBest, double pWait) {
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

void Solver_ILS::workILSrk(DWORD coreIndex, unique_ptr<Schedule>& sched, Sched_params* schedParams, ILS_params* ilsParams, initializerRK<pJob> init, const vector<double>& randomKeys, int iTilimSeconds, std::chrono::time_point<std::chrono::high_resolution_clock> start, ILS_Thread* localBest, double pWait) {
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
        (tempSched.get()->*init)(sortJobsByRK, randomKeys, *schedParams);

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

vector<DWORD> Solver_ILS::GetPCoreIndices() {
    cout << "Optimized for 12th Gen Intel(R) Core(TM) i7-12700: using 8 performance cores..." << endl;
    return { 0, 1, 2, 3, 4, 5, 6, 7 };
    
    //SYSTEM_INFO si;
    //GetSystemInfo(&si);

    //// Typisch: P-Cores sind die ersten logischen Kerne
    //DWORD numPCores = si.dwNumberOfProcessors / 2;
    //vector<DWORD> pcores;

    //for (DWORD i = 0; i < numPCores; i++) {
    //    pcores.push_back(i);
    //}
    //return pcores;
}



// +++++++ SEQUENCE BASED ILS +++++++++++++
Solver_Sequence_ILS::Solver_Sequence_ILS(Schedule* schedule, Sched_params* schedParameters, ILS_params* ilsParameters) : bestChr(schedule->getN()), currentChr(schedule->getN()) {
    masterSched = schedule;
    schedParams = schedParameters;
    params = ilsParameters;
}
Solver_Sequence_ILS::~Solver_Sequence_ILS(){}

double Solver_Sequence_ILS::solveILSseq(Schedule& sched, int iTilimSeconds) {
    auto start = chrono::high_resolution_clock::now();
    chrono::seconds usedTime;
    chrono::time_point<chrono::high_resolution_clock> stop;

    double bestTWT = DBL_MAX;
    vector<double> bestChr = vector<double>(sched.getN());

    params->multiStartIterations = 0;

    do {    // MULTISTART-LOOP
            // INITIALIZE
            // DEFAULT: random initializer (TODO: other initializers to be defined)
        vector<double> tempChr = bestChr;
        initRandomPermutation(tempChr);

        params->ilsIterations.push_back(0);
        do {// ILS LOOP
            // LOCAL SEARCH
            bool bJobInsertApplied = false;
            bool bJobSwapApplied = false;
           
            do { 
                bJobInsertApplied = localSearchInsertJob(tempChr, params->applyBestFit);
                bJobSwapApplied = localSearchSwapJob(tempChr, params->applyBestFit);
            } while (bJobInsertApplied || bJobSwapApplied);

            double tempTWT = decodeAndGetTWT(tempChr);
            if (tempTWT < bestTWT) {
                bestTWT = tempTWT;
                bestChr = tempChr;
                stop = chrono::high_resolution_clock::now();
                usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
                params->bestAfterSeconds = usedTime.count();
            }

            // DEBÚGGING
            cout << "ILS iteration " << params->multiStartIterations + 1 << "." << params->ilsIterations[params->multiStartIterations] + 1 << "tempTWT: " << tempTWT << ", bestTWT" << bestTWT << endl;
            
            // PERTURBATION
        
            uniform_real_distribution<> perturbDistrib(0, 1);
            for (size_t i = 0; i < params->nPerturbationSteps; ++i) {
                double perturbChoice = perturbDistrib(TCB::rng);
                if (perturbChoice < 1.0 / 3.0) {
                    //TCB::logger.Log(Info, "perturbRandomJobSwap");
                    perturbJobInsert(tempChr);
                }
                else if(perturbChoice < (1.0 / 3.0) * 2.0) {
                    //TCB::logger.Log(Info, "perturbRandomJobRightShifting");
                    perturbJobSwap(tempChr);
                }
                else {
                    perturbRandomKey(tempChr);
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

    formMasterSchedule(bestChr);
    return bestTWT;
}

double Solver_Sequence_ILS::solveILSseqParallelized(Schedule& sched, int iTilimSeconds) {
    unsigned int nCores = thread::hardware_concurrency();

    vector<DWORD> pCores = GetPCoreIndices();
    ////cout << "Using " << pCores.size() << " cores..." << endl;

    auto start = chrono::high_resolution_clock::now();
    chrono::seconds usedTime;
    chrono::time_point<chrono::high_resolution_clock> stop;

    vector<ILSseq_Thread> ILSseq_threads(nCores);
    vector<thread> threads;
    for (unsigned int core = 0; core < pCores.size(); ++core) {   
        DWORD coreIndex = pCores[core];
        ILSseq_threads[core].bestChr = vector<double>(sched.getN());
        threads.emplace_back(workILSseq, coreIndex, sched.clone(), schedParams, params, iTilimSeconds, start, &ILSseq_threads[core]);
    }

    for (auto& t : threads) t.join();

    int globalBestIdx = 0;
    double globalBestTWT = ILSseq_threads[0].bestTWT;
    for (int i = 0; i < pCores.size(); ++i) { 
        
        // DEBUGGING
        //cout << "Thread " << i << ": TWT= " << ILSseq_threads[i].bestTWT << ", ";
        //coutDblVec(ILSseq_threads[i].bestChr);

        if (ILSseq_threads[i].bestTWT < globalBestTWT) {
            globalBestTWT = ILSseq_threads[i].bestTWT;
            globalBestIdx = i;
            params->bestAfterSeconds = ILSseq_threads[i].bestAfterSeconds;
            params->ilsIterations = ILSseq_threads[i].ilsIterations;
            params->multiStartIterations = ILSseq_threads[i].multiStartIterations;
        }
    }

    formMasterSchedule(ILSseq_threads[globalBestIdx].bestChr);
    return globalBestTWT;
}

void Solver_Sequence_ILS::insertJob(size_t jobIdx, size_t posIdx) {
    if (jobIdx < 0 || jobIdx >= currentChr.size() || posIdx < 0 || posIdx >= currentChr.size()) throw out_of_range("Solver_Sequence_ILS::insertJob(...) out of range");
    if (jobIdx < posIdx) {
        rotate(currentChr.rbegin() + (currentChr.size() - 1 - posIdx), currentChr.rbegin() + (currentChr.size() - 1 - jobIdx), currentChr.rend());
    }
    if (posIdx < jobIdx) {
        rotate(currentChr.begin() + posIdx, currentChr.begin() + jobIdx, currentChr.begin() + jobIdx + 1);
    }
}
void Solver_Sequence_ILS::swapJob(size_t firstIdx, size_t secondIdx) {
    if (firstIdx < 0 || firstIdx >= currentChr.size() || secondIdx < 0 || secondIdx >= currentChr.size()) throw out_of_range("Solver_Sequence_ILS::swapJob(...) out of range");
    if (firstIdx != secondIdx) {
        swap(currentChr[firstIdx], currentChr[secondIdx]);
    }
}
double Solver_Sequence_ILS::evaluateJobInsert(size_t jobIdx, size_t posIdx) {
    vector<double> tempChr = currentChr; // bestChr;
    insertJob(tempChr, jobIdx, posIdx);
    return currentTWT - decodeAndGetTWT(tempChr);
}
double Solver_Sequence_ILS::evaluateJobSwap(size_t firstIdx, size_t secondIdx) {
    vector<double> tempChr = currentChr; // bestChr;
    swapJob(tempChr, firstIdx, secondIdx);
    return currentTWT - decodeAndGetTWT(tempChr);
}
bool Solver_Sequence_ILS::localSearchInsertJob(bool bestFit) {
    // DEFAULT: random order (TODO: other sorting orders to be defined)
    vector<size_t> indices(currentChr.size());
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
                    double tempImprovement = evaluateJobInsert(currentChr, indices[i], indices[j]);
                    if (tempImprovement > bestImprovement) {
                        bestImprovement = tempImprovement;
                        bestJobIdx = indices[i];
                        bestPosIdx = indices[j];
                        if (!bestFit && bestImprovement > 0) {     // FIRST FIT
                            insertJob(currentChr, bestJobIdx, bestPosIdx);
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
            insertJob(currentChr, bestJobIdx, bestPosIdx);
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
    vector<size_t> indices(currentChr.size());
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
                    double tempImprovement = evaluateJobSwap(currentChr, indices[i], indices[j]);
                    if (tempImprovement > bestImprovement) {
                        bestImprovement = tempImprovement;
                        bestFirstIdx = indices[i];
                        bestSecondIdx = indices[j];
                        if (!bestFit && bestImprovement > 0) {     // FIRST FIT
                            swapJob(currentChr, bestFirstIdx, bestSecondIdx);
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
            swapJob(currentChr, bestFirstIdx, bestSecondIdx);
            currentTWT = currentTWT - bestImprovement;
            bestImprovement = 0;
            bEventuallyImproved = true;
            bImproved = true;
        }
    }

    return bEventuallyImproved;
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
double Solver_Sequence_ILS::evaluateJobInsert(vector<double>& chromosome, size_t jobIdx, size_t posIdx) {
    vector<double> tempChr = chromosome;
    insertJob(tempChr, jobIdx, posIdx);
    return currentTWT - decodeAndGetTWT(tempChr);
}

double Solver_Sequence_ILS::evaluateJobSwap(vector<double>& chromosome, size_t firstIdx, size_t secondIdx) {
    vector<double> tempChr = chromosome;
    swapJob(tempChr, firstIdx, secondIdx);
    return currentTWT - decodeAndGetTWT(tempChr);
}

bool Solver_Sequence_ILS::localSearchInsertJob(vector<double>& chromosome, bool bestFit) {
    // DEFAULT: random order (TODO: other sorting orders to be defined)
    vector<size_t> indices(chromosome.size());
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
                    double tempImprovement = evaluateJobInsert(chromosome, indices[i], indices[j]);
                    if (tempImprovement > bestImprovement) {
                        bestImprovement = tempImprovement;
                        bestJobIdx = indices[i];
                        bestPosIdx = indices[j];
                        if (!bestFit && bestImprovement > 0) {     // FIRST FIT
                            insertJob(chromosome, bestJobIdx, bestPosIdx);
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
            insertJob(chromosome, bestJobIdx, bestPosIdx);
            currentTWT = currentTWT - bestImprovement;
            bImproved = true;
            bEventuallyImproved = true;
            bestImprovement = 0;
        }
    }

    return bEventuallyImproved;
}
bool Solver_Sequence_ILS::localSearchSwapJob(vector<double>& chromosome, bool bestFit) {
    // DEFAULT: random order (TODO: other sorting orders to be defined)
    vector<size_t> indices(chromosome.size());
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
                    double tempImprovement = evaluateJobSwap(chromosome, indices[i], indices[j]);
                    if (tempImprovement > bestImprovement) {
                        bestImprovement = tempImprovement;
                        bestFirstIdx = indices[i];
                        bestSecondIdx = indices[j];
                        if (!bestFit && bestImprovement > 0) {     // FIRST FIT
                            swapJob(chromosome, bestFirstIdx, bestSecondIdx);
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
            swapJob(chromosome, bestFirstIdx, bestSecondIdx);
            currentTWT = currentTWT - bestImprovement;
            bestImprovement = 0;
            bEventuallyImproved = true;
            bImproved = true;
        }
    }

    return bEventuallyImproved;
}

std::vector<double> Solver_Sequence_ILS::getBestChr() {
    return bestChr;
}

void Solver_Sequence_ILS::perturbJobInsert() {
    uniform_int_distribution<> idxDistrib(0, currentChr.size() - 1);
    size_t i = idxDistrib(TCB::rng);
    size_t j = i;
    while (j == i) {
        j = idxDistrib(TCB::rng);
    }
    insertJob(currentChr, i, j);
    currentTWT = decodeAndGetTWT(currentChr);
}
void Solver_Sequence_ILS::perturbJobSwap() {
    uniform_int_distribution<> idxDistrib(0, currentChr.size() - 1);
    size_t i = idxDistrib(TCB::rng);
    size_t j = i;
    while (j == i) {
        j = idxDistrib(TCB::rng);
    }
    swapJob(currentChr, i, j);
    currentTWT = decodeAndGetTWT(currentChr);
}
void Solver_Sequence_ILS::perturbRandomKey() {
    uniform_real_distribution<> rkDistrib(0.0, 1.0);
    uniform_int_distribution<> idxDistrib(0, currentChr.size() - 1);
    size_t idx = idxDistrib(TCB::rng);
    double newRkValue = rkDistrib(TCB::rng);
    currentChr[idx] = newRkValue;
    currentTWT = decodeAndGetTWT(currentChr);
}

void Solver_Sequence_ILS::perturbJobInsert(vector<double>& chromosome){
    uniform_int_distribution<> idxDistrib(0, chromosome.size() - 1);
    size_t i = idxDistrib(TCB::rng);
    size_t j = i;
    while (j == i) {
        j = idxDistrib(TCB::rng);
    }
    insertJob(chromosome, i, j);
    currentTWT = decodeAndGetTWT(chromosome);
}
void Solver_Sequence_ILS::perturbJobSwap(vector<double>& chromosome) {
    uniform_int_distribution<> idxDistrib(0, chromosome.size() - 1);
    size_t i = idxDistrib(TCB::rng);
    size_t j = i;
    while (j == i) {
        j = idxDistrib(TCB::rng);
    }
    swapJob(chromosome, i, j);
    currentTWT = decodeAndGetTWT(chromosome);
}
void Solver_Sequence_ILS::perturbRandomKey(vector<double>& chromosome) {
    uniform_real_distribution<> rkDistrib(0.0, 1.0);
    uniform_int_distribution<> idxDistrib(0, chromosome.size() - 1);
    size_t idx = idxDistrib(TCB::rng);
    double newRkValue = rkDistrib(TCB::rng);
    chromosome[idx] = newRkValue;
    currentTWT = decodeAndGetTWT(chromosome);
}
void Solver_Sequence_ILS::decode(Schedule* sched, std::vector<double>& chromosome) {
    sched->lSchedJobsWithRandomKeySorting(sortJobsByRK, chromosome);
}
void Solver_Sequence_ILS::formMasterSchedule(std::vector<double>& chromosome) {
    decode(masterSched, chromosome);
}
double Solver_Sequence_ILS::decodeAndGetTWT(std::vector<double>& chr) {
    unique_ptr<Schedule> mySched = masterSched->clone();
    mySched->lSchedJobsWithRandomKeySorting(sortJobsByRK, chr);
    return mySched->getTWT();
}
double Solver_Sequence_ILS::staticDecodeAndGetTWT(Schedule* sched, vector<double>& chr) {
    unique_ptr<Schedule> mySched = sched->clone();
    mySched->lSchedJobsWithRandomKeySorting(sortJobsByRK, chr);
    return mySched->getTWT();
}
void Solver_Sequence_ILS::initRandomPermutation(){
    uniform_real_distribution<> rkDistrib(0.0, 1.0);
    generate(currentChr.begin(), currentChr.end(), [&] { return rkDistrib(TCB::rng); });
    currentTWT = decodeAndGetTWT(currentChr);
}

void Solver_Sequence_ILS::initRandomPermutation(vector<double>& chromosome) {
    uniform_real_distribution<> rkDistrib(0.0, 1.0);
    generate(chromosome.begin(), chromosome.end(), [&] { return rkDistrib(TCB::rng); });
}

void Solver_Sequence_ILS::workILSseq(DWORD coreIndex, unique_ptr<Schedule>& sched, Sched_params* schedParams, ILS_params* ilsParams, int iTilimSeconds, chrono::time_point<chrono::high_resolution_clock> start, ILSseq_Thread* localBest) {
    DWORD_PTR mask = 1ULL << coreIndex;
    SetThreadAffinityMask(GetCurrentThread(), mask);

    double bestTWT = DBL_MAX;
    //vector<double> bestChr = vector<double>(sched->getN());

    chrono::seconds usedTime;
    chrono::time_point<chrono::high_resolution_clock> stop;

    localBest->multiStartIterations = 0;

    Solver_Sequence_ILS worker = Solver_Sequence_ILS(sched.get(), schedParams, ilsParams);

    do {// MULTISTART-LOOP
        // INITIALIZE
        worker.initRandomPermutation();
        vector<double> tempChr = worker.currentChr;

        localBest->ilsIterations.push_back(0);
        do {// ILS LOOP
            // LOCAL SEARCH
            bool bJobInsertApplied = false;
            bool bJobSwapApplied = false;
            do {
                bJobInsertApplied = worker.localSearchInsertJob(ilsParams->applyBestFit);
                bJobSwapApplied = worker.localSearchSwapJob(ilsParams->applyBestFit);
            } while (bJobInsertApplied || bJobSwapApplied);

            double tempTWT = worker.currentTWT;
            if (tempTWT < localBest->bestTWT) {
                localBest->bestTWT = tempTWT;
                localBest->bestChr = worker.currentChr;
                stop = chrono::high_resolution_clock::now();
                usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
                localBest->bestAfterSeconds = usedTime.count();
            }
            // PERTURBATION
            uniform_real_distribution<> perturbDistrib(0, 1);
            for (size_t i = 0; i < ilsParams->nPerturbationSteps; ++i) {
                double perturbChoice = perturbDistrib(TCB::rng);
                if (perturbChoice < 1.0 / 3.0) {
                    //TCB::logger.Log(Info, "perturbRandomJobSwap");
                    worker.perturbJobInsert();
                }
                else if (perturbChoice < (1.0 / 3.0) * 2.0) {
                    //TCB::logger.Log(Info, "perturbRandomJobRightShifting");
                    worker.perturbJobSwap();
                }
                else {
                    worker.perturbRandomKey();
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

Solver_Hybrid_ILS::Solver_Hybrid_ILS(Schedule* schedule, Sched_params* schedParameters, ILS_params* ilsParameters) : phase1(schedule, schedParameters, ilsParameters), phase2(*schedParameters, *ilsParameters) {
    // phase1 = Solver_Sequence_ILS(schedule, schedParameters, ilsParameters);
    // phase2 = Solver_ILS(*schedParameters, *ilsParameters);
    schedParams = schedParameters;
    ilsParams = ilsParameters; 
}

Solver_Hybrid_ILS::~Solver_Hybrid_ILS() {}

double Solver_Hybrid_ILS::solveILShybrid(Schedule& sched, int iTilimTotal) {
    cout << "Solver_Hybrid_ILS::solveILShybrid(...) not yet implemented." << endl;
    double bestTWT = DBL_MAX;
    int iTilimPhase1 = (int)(ilsParams->firstPhaseTimeLimitAllocation * (double)iTilimTotal);
    int iTilimPhase2 = iTilimTotal - iTilimPhase1;

    // PHASE 1 - ILS on sequence/permutation of jobs
    bestTWT = phase1.solveILSseqParallelized(sched, iTilimPhase1);
    //TODO report twt after phase1

    vector<double> chromosome = phase1.getBestChr();    // TODO return different permutations from each core
    sched.reset();

    // PHASE 2 - ILS on schedule (initial schedule from phase 1)
    initializerRK<pJob> init2 = &Schedule::lSchedJobsWithRandomKeySorting;
    bestTWT = phase2.solveILSparallelized(sched, init2, chromosome, iTilimPhase2);
    return bestTWT;
}
