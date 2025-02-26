#include "Solver_ILS.h"
#include "Schedule.h"

using namespace std;

Solver_ILS::Solver_ILS(Sched_params& params) : schedParams(&params) {}

double Solver_ILS::solveILS(Schedule& sched, initializer<pJob> init, prioRule<pJob> rule, int iTilimSeconds, double pWait) {
    // initialize
    unique_ptr<Schedule> tempSched = sched.clone();
    (tempSched.get()->*init)(rule, *schedParams);

    // TODO: ILS


    return 0.0;
}
