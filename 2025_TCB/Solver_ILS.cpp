#include <chrono>

#include "Solver_ILS.h"
#include "Schedule.h"

using namespace std;

Solver_ILS::Solver_ILS(Sched_params& params) : schedParams(&params) {}

double Solver_ILS::solveILS(Schedule& sched, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, int iStarts, double pWait) {
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
            // TODO: LOCAL SEARCH
            tempSched->saveJsonFactory("BEFORE_LOCAL_SEARCH");   // DEBUGGING
            tempSched->localSearchJobLeftShifting(&sortJobsRandomly, true);      
            tempSched->saveJsonFactory("AFTER_JOBLEFTSHIFTING");   // DEBUGGING
            tempSched->localSearchJobSwapping(&sortJobsRandomly, true);
            tempSched->saveJsonFactory("AFTER_JOBSWAPPING");   // DEBUGGING
            tempSched->localSearchBatchConsolidation(true);
            tempSched->saveJsonFactory("AFTER_BATCHCONSOLIDATION");   // DEBUGGING

            double tempTWT = tempSched->getTWT();
            if (tempTWT < bestTWT) {
                bestTWT = tempTWT;
                bestSched = tempSched->clone();
            }
            // TODO: PERTURBATION


            ++iterationCounter;
            stop = chrono::high_resolution_clock::now();
            usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
        } while (usedTime.count() < (double)iTilimSeconds / (double)iStarts);   // MULTISTART
        stop = chrono::high_resolution_clock::now();
        usedTime = chrono::duration_cast<chrono::seconds>(stop - start);
    } while (usedTime.count() < iTilimSeconds); 


    return bestTWT;
}
