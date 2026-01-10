#include <chrono>
#include <iostream>

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
           cout << "ILS iteration " << params->multiStartIterations + 1 << "." << params->ilsIterations[params->multiStartIterations] + 1 << " TWT: " << bestTWT << endl;
           TCB::logger.Log(Info, to_string(params->ilsIterations[params->multiStartIterations] + 1));
           if (params->multiStartIterations + 1 == 1 && params->ilsIterations[params->multiStartIterations] + 1 == 20) {
               //tempSched->saveJsonFactory("DEBUGGING");
               int debugger = 666;
                
            }
           if (!tempSched->isValid()) {
               int debugger = 666;
           }
            
            //cout << "ILS LS JOB SHIFT" << endl;
            tempSched->localSearchJobLeftShifting(&sortJobsRandomly, params->applyBestFit); 
            if (!tempSched->isValid()) {
                int debugger = 666;
            }
            //cout << "ILS LS JOB SWAP" << endl;
            tempSched->localSearchJobSwapping(&sortJobsRandomly, params->applyBestFit);
            //cout << "ILS LS BATCH CON" << endl;
            tempSched->localSearchBatchConsolidation(params->applyBestFit);

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
                    //cout << "ILS Perturbation JOB SWAP" << endl;
                    tempSched->perturbRandomJobSwap();
                }
                else {
                    //cout << "ILS Perturbation JOB SHIFT" << endl;
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

ILS_params Solver_ILS::getDefaultParams() {
    ILS_params ilsParams = ILS_params();
    ilsParams.nStarts = 1;
    ilsParams.nPerturbationSteps = 5;
    ilsParams.applyBestFit = true;
    ilsParams.randomizedLocalSearchSequence = false;
    ilsParams.multiStartIterations = 0;
    ilsParams.ilsIterations = vector<int>();
    return ilsParams;
}
