#include <iostream>
#include <thread>

#include "Solver_ILS.h"
#include "Schedule.h"

using namespace std;

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
           //TCB::logger.Log(Info, to_string(params->ilsIterations[params->multiStartIterations] + 1));
           /*if (params->multiStartIterations + 1 == 1 && params->ilsIterations[params->multiStartIterations] + 1 == 5) {
                tempSched->saveJsonFactory("DEBUGGING");
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
                    tempSched->perturbRandomJobSwap();
                }
                else {
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
    auto start = chrono::high_resolution_clock::now();
    chrono::seconds usedTime;
    chrono::time_point<chrono::high_resolution_clock> stop;


    vector<ILS_Thread> ILS_threads(nCores);
    vector<thread> threads;
    for (unsigned int core = 0; core < nCores; ++core) {
        ILS_threads[core].bestSched = sched.clone();
        threads.emplace_back(workerILS, move(ILS_threads[core].bestSched), schedParams, params, init, rule, iTilimSeconds, start, &ILS_threads[core], pWait);
    }

    for (auto& t : threads) t.join();

    int globalBestIdx = 0;
    double globalBestTWT = ILS_threads[0].bestTWT;
    for (int i = 0; i < nCores; ++i) {
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

void Solver_ILS::workerILS(unique_ptr<Schedule>& sched, Sched_params* schedParams, ILS_params* ilsParams, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, chrono::time_point<chrono::high_resolution_clock> start, ILS_Thread* localBest, double pWait) {
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
        /* cout << "Thread id " << this_thread::get_id() << endl;
        cout << *tempSched;*/

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
