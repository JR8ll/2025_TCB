#include "GaDecoder.h"
#include "Functions.h"
#include "Schedule.h"
#include "Solver_Ga.h"

using namespace std;

GaDecoderJobListSched::GaDecoderJobListSched(Schedule* schedule, Sched_params* schedParameters, GA_params* gaParameters) : masterSched(schedule), schedParams(schedParameters), gaParams(gaParameters) {}
GaDecoderJobListSched::~GaDecoderJobListSched() {}
double GaDecoderJobListSched::decode(const vector<double>& chr) const {
	unique_ptr<Schedule> mySched = masterSched->clone();
	mySched->lSchedJobsWithRandomKeySorting(sortJobsByRK, chr, *schedParams);
	return mySched->getTWT();
}
void GaDecoderJobListSched::applyPersistentLocalSearch(unique_ptr<Schedule> sched, vector<double>& chr) {
    bool bJobSwapApplied = false;
    pair<size_t, size_t> swap = make_pair(0, 0);
    do {
        bJobSwapApplied = sched->localSearchJobSwapping(&sortJobsRandomly, gaParams->applyLocalSearchBestFit);
    } while (bJobSwapApplied);

}
double GaDecoderJobListSched::applyNonPersitentLocalSearch(unique_ptr<Schedule> sched) const {
    bool bLeftShiftApplied = false;
    bool bJobSwapApplied = false;
    bool bBatchConsolidationApplied = false;
    do {
        bLeftShiftApplied = sched->localSearchJobLeftShifting(&sortJobsRandomly, gaParams->applyLocalSearchBestFit);
        bJobSwapApplied = sched->localSearchJobSwapping(&sortJobsRandomly, gaParams->applyLocalSearchBestFit);
        bBatchConsolidationApplied = sched->localSearchBatchConsolidation(gaParams->applyLocalSearchBestFit);
    } while (bLeftShiftApplied || bJobSwapApplied || bBatchConsolidationApplied);
	
    return sched->getTWT();
}
void GaDecoderJobListSched::formSchedule(const std::vector<double>& chr)
{
	masterSched->lSchedJobsWithRandomKeySorting(sortJobsByRK, chr);
}
