#include <chrono>

#include "Solver_ILS.h"
#include "Schedule.h"

using namespace std;

Solver_ILS::Solver_ILS(Sched_params& schedParams, ILS_params& params) {
    this->schedParams = &schedParams;
    this->params = &params;
}

double Solver_ILS::solveILS(Schedule& sched, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, double pWait) {
    // TODO ILS PARAMETERS (BESTFIT/FIRSTFIT, JOBSORTINGORDER, etc.)
    auto start = chrono::high_resolution_clock::now();
    chrono::seconds usedTime;
    chrono::time_point<chrono::high_resolution_clock> stop;

    double bestTWT = DBL_MAX;
    unique_ptr<Schedule> bestSched = sched.clone();
    
    // TODO: ILS
    int iterationCounter = 0;
    do {
        // INITIALIZE
        unique_ptr<Schedule> tempSched = sched.clone();
        (tempSched.get()->*init)(rule, *schedParams);
        do {
            // LOCAL SEARCH
            tempSched->saveJsonFactory("BEFORE_LOCAL_SEARCH");   // DEBUGGING
            tempSched->localSearchJobLeftShifting(&sortJobsRandomly, params->applyBestFit);      
            tempSched->saveJsonFactory("AFTER_JOBLEFTSHIFTING");   // DEBUGGING
            tempSched->localSearchJobSwapping(&sortJobsRandomly, params->applyBestFit);
            tempSched->saveJsonFactory("AFTER_JOBSWAPPING");   // DEBUGGING
            tempSched->localSearchBatchConsolidation(params->applyBestFit);
            tempSched->saveJsonFactory("AFTER_BATCHCONSOLIDATION");   // DEBUGGING

            double tempTWT = tempSched->getTWT();
            if (tempTWT < bestTWT) {
                bestTWT = tempTWT;
                bestSched = tempSched->clone();
            }
            // TODO: PERTURBATION
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

            ++iterationCounter;
            stop = chrono::high_resolution_clock::now();
            usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
        } while (usedTime.count() < (double)iTilimSeconds / (double)params->nStarts);   // MULTISTART
        stop = chrono::high_resolution_clock::now();
        usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
    } while (usedTime.count() < iTilimSeconds); 


    return bestTWT;
}

ILS_params Solver_ILS::getDefaultParams() {
    ILS_params ilsParams = ILS_params();
    ilsParams.nStarts = 1;
    ilsParams.nPerturbationSteps = 5;
    ilsParams.applyBestFit = true;
    ilsParams.randomizedLocalSearchSequence = false;
    return ilsParams;
}
