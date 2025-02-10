#include "GaDecoder.h"
#include "Functions.h"
#include "Schedule.h"

using namespace std;


GaDecoderJobLs::GaDecoderJobLs(Schedule* schedule) : masterSched(schedule) {}
GaDecoderJobLs::~GaDecoderJobLs() {}
double GaDecoderJobLs::decode(const vector<double>& chr) const {
	unique_ptr<Schedule> mySched = masterSched->clone();
	mySched->lSchedJobsWithRandomKeySorting(sortJobsByRK, chr);
	return mySched->getTWT();
}
void GaDecoderJobLs::formSchedule(const std::vector<double>& chr)
{
	masterSched->lSchedJobsWithRandomKeySorting(sortJobsByRK, chr);
}
