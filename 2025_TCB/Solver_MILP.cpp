#include<ilcplex/ilocplex.h>
#include<ilconcert/ilotupleset.h>
#include<map>
#include<set>

#include "Solver_MILP.h"
#include "Schedule.h"

ILOSTLBEGIN

Solver_MILP::Solver_MILP() {}
Solver_MILP::~Solver_MILP() {}

double Solver_MILP::solveDecompositionMILP(Schedule* schedule, int nDash) {
	//IloEnv env;

	//// PARAMETERS
	//IloNum omega = 9;	// TODO set omega

	//IloInt G = 99999;	// TODO compute
	//IloInt m = schedule->size();

	//IloIntTupleSet stageMachine(env, 2);
	//for (int o = 0; o < o; ++o) {
	//	for (int i = 0; i < (*schedule)[o].size(); ++i) {
	//		stageMachine.add(IloIntArray(env, 2, o, i));
	//	}
	//}
	//map<pair<int, int>, int> mapStgMac = map<pair<int, int>, int>();	// for access to tuple "stageMachine"
	//int index = 0;
	//for (IloIntTupleSetIterator it(env, stageMachine); it.ok(); ++it) {
	//	mapStgMac.insert(make_pair(make_pair((*it)[0], (*it)[1]), index));
	//	index++;
	//}

	//IloIntArray B(env);	// capacity[stage]
	//for (int o = 0; o < m; ++o) {
	//	int tempCap = schedule->getCapAtStageIdx(o);	// assumption: parallel identical machines
	//	if (tempCap > 1) {
	//		B.add(tempCap);
	//	}
	//}

	//IloIntArray3 B_iob(env);	// number of jobs assigned to preformed Batch[Products][Stages][Batches]
	//IloNumArray3 S_iob(env);	// start times of preformed Batch[Products][Stages][Batches]	
	//vector<vector<set<int> > > B_io_preformed = vector<vector<set<int> > >(problem.getF());


	// TODO
	return 666;
}




