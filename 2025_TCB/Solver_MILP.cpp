#include<ilcplex/ilocplex.h>
#include<ilconcert/ilotupleset.h>
#include<algorithm>
#include<map>
#include<numeric>
#include<set>

#include "Solver_MILP.h"
#include "Functions.h"
#include "Problem.h"
#include "Schedule.h"

ILOSTLBEGIN

using IloBoolVarArray2 = IloArray<IloBoolVarArray>;
using IloBoolVarArray3 = IloArray<IloBoolVarArray2>;
using IloNumVarArray2 = IloArray<IloNumVarArray>;
using IloBoolArray2 = IloArray<IloBoolArray>;

Solver_MILP::Solver_MILP(Sched_params& sched_params, DECOMPMILP_params& decompParams, double omega) : schedParams(&sched_params), params(&decompParams), omega(omega) {}
Solver_MILP::~Solver_MILP() {}

double Solver_MILP::solveJobBasedMILP(Schedule* schedule, int nDash, int cplexTilim) {
	if (schedule->getProblem() == nullptr) throw ExcSched("solveDecompositionMILP missing problem reference");

	vector<pJob> consideredJobs = vector<pJob>();
	for (size_t j = 0; j < nDash; ++j) {
		try {
			consideredJobs.push_back(move(schedule->get_pJob(0)));
			cout << consideredJobs[0]->getId();	// DEBUGGING
		} catch (const out_of_range& ex) {
			break;	// all jobs considered
		}
	}

	IloEnv env;

	// PARAMETERS
	IloNum omega = this->omega;	// TODO set omega

	IloNum G = schedule->getProblem()->getG();	
	IloInt m = schedule->size();

	int nProducts = schedule->getProblem()->getF();
	int nStages = schedule->size();
	int nJobs = consideredJobs.size();
	int nBatches = nJobs + schedule->getNumberOfScheduledJobs();	// new and (max. number of) already formed batches

	IloIntTupleSet stageMachine(env, 2);
	for (int o = 0; o < nStages; ++o) {
		for (int l = 0; l < (*schedule)[o].size(); ++l) {
			stageMachine.add(IloIntArray(env, 2, o, l));
		}
	}
	map<pair<int, int>, int> mapStgMac = map<pair<int, int>, int>();	// for access to tuple "stageMachine"
	int index = 0;
	for (IloIntTupleSetIterator it(env, stageMachine); it.ok(); ++it) {
		mapStgMac.insert(make_pair(make_pair((*it)[0], (*it)[1]), index));
		index++;
	}
	const map<pair<int, int>, int> constMapStgMac = mapStgMac;

	IloIntArray B(env);	// capacity[stage]
	for (int o = 0; o < m; ++o) {
		int tempCap = schedule->getCapAtStageIdx(o);	// assumption: parallel identical machines
		if (tempCap > 1) {
			B.add(tempCap);
		}
	}

	IloIntArray3 B_iob(env);	// number of jobs assigned to preformed Batch[Products][Stages][Batches]
	IloNumArray3 S_iob(env);	// start times of preformed Batch[Products][Stages][Batches]	
	vector<vector<set<int>>> B_io_preformed = vector<vector<set<int>>>(nProducts);
	for (int i = 0; i < nProducts; ++i) {
		B_io_preformed[i] = vector<set<int>>(nStages);
		B_iob.add(IloIntArray2(env));
		S_iob.add(IloNumArray2(env));
		for (int o = 0; o < nStages; ++o) {
			B_io_preformed[i][o] = set<int>();
			B_iob[i].add(IloIntArray(env, nBatches));
			S_iob[i].add(IloNumArray(env, nBatches));
		}
	}

	vector<map<pair<int, int>, pair<int, int>>> batchIdxMac = vector<map<pair<int, int>, pair<int, int>>>(nStages);
	for (int o = 0; o < nStages; ++o) {
		batchIdxMac[o] = map<pair<int, int>, pair<int, int> >();	// first pair: model indices, second pair: schedule indices
		int batchIdx = 0;
		for (int l = 0; l < (*schedule)[o].size(); ++l) {
			for (int b = 0; b < (*schedule)[o][l].size(); ++b) {
				int productIdx = (*schedule)[o][l][b].getF() - 1;
				if (productIdx >= 0) {
					B_io_preformed[productIdx][o].insert(batchIdx);
					batchIdxMac[o].insert(make_pair(make_pair(productIdx, batchIdx), make_pair(l, b)));

					int debugB_iob = (*schedule)[o][l][b].size();
					double debugS_iob = (*schedule)[o][l][b].getStart();

					B_iob[productIdx][o][batchIdx] = (*schedule)[o][l][b].size();
					S_iob[productIdx][o][batchIdx] = (*schedule)[o][l][b].getStart();
					batchIdx++;
				}
			}
		}
	}
	
	// ready times of machines
	IloNumArray2 rm(env);	// rm[Stages][Machines] machine ready times
	for (int o = 0; o < nStages; ++o) {
		int m = (*schedule)[o].size();
		rm.add(IloNumArray(env, m));
		for (int i = 0; i < m; ++i) {
			rm[o][i] = (*schedule)[o][i].getMSP();
		}
	}

	IloNumArray d(env);		// d[Jobs] job due dates
	IloNumArray r(env);		// r[Jobs] job ready times
	IloNumArray w(env);		// w[Jobs] job weights
	IloIntArray f(env);		// f[Jobs] jobs´ products or families

	for (int j = 0; j < nJobs; ++j) {
		d.add(consideredJobs[j]->getD());
		r.add(consideredJobs[j]->getR());
		w.add(consideredJobs[j]->getW());
		f.add(consideredJobs[j]->getF());
	}

	IloNumArray2 p(env);	// p[Products][Stages] processing times
	for (int i = 0; i < nProducts; ++i) {
		p.add(IloNumArray(env, nStages));
		for (int o = 0; o < nStages; ++o) {
			p[i][o] = schedule->getProblem()->getProduct(i)->getP(o); 
		}
	}

	IloNumArray3 tc(env);	// tc[Products][Stages][Stages] time constraints
	for (int i = 0; i < nProducts; ++i) {
		tc.add(IloNumArray2(env));
		for (int o1 = 0; o1 < nStages; ++o1) {
			tc[i].add(IloNumArray(env, nStages));
			const vector<pair<int, double>>& tcMaxFwd = schedule->getProblem()->getProduct(i)->getTcMaxFwd(o1);
			for (size_t o2 = 0; o2 < nStages; ++o2) {
				double tcDuration = 999999;	// TODO define INF
				for (auto it = tcMaxFwd.begin(); it != tcMaxFwd.end(); ++it) {
					if (it->first == o2) {
						tcDuration = it->second;
					}
				}
				tc[i][o1][o2] = (IloNum)tcDuration;
			}
		}
	}

	IloModel mod(env);

	// DECISION VARIABLES
	IloBoolVarArray3 x(env);	// x[Jobs][Batches][<StageMachine>]
	for (int j = 0; j < nJobs; ++j) {
		x.add(IloBoolVarArray2(env));
		for (int b = 0; b < nBatches; ++b) {
			x[j].add(IloBoolVarArray(env, stageMachine.getCardinality()));
			for (int k = 0; k < stageMachine.getCardinality(); ++k) {
				x[j][b][k] = IloBoolVar(env);
				string varName = "x_" + to_string(j) + "," + to_string(b) + "," + to_string(k);
				x[j][b][k].setName(varName.c_str());
			}
		}
	}

	IloBoolVarArray3 y(env);	// y[Products][Batches][<StageMachine>]
	for (int i = 0; i < nProducts; ++i) {
		y.add(IloBoolVarArray2(env));
		for (int b = 0; b < nBatches; ++b) {
			y[i].add(IloBoolVarArray(env, stageMachine.getCardinality()));
			for (int k = 0; k < stageMachine.getCardinality(); ++k) {
				y[i][b][k] = IloBoolVar(env);
				string varName = "y_" + to_string(i) + "," + to_string(b) + "," + to_string(k);
				y[i][b][k].setName(varName.c_str());
			}
		}
	}

	IloBoolVarArray3 z(env);	// z[Jobs][Stages][Batches]	
	for (int j = 0; j < nJobs; ++j) {
		z.add(IloBoolVarArray2(env));
		for (int o = 0; o < nStages; ++o) {
			z[j].add(IloBoolVarArray(env, nBatches));		
			for (int b = 0; b < nBatches; ++b) {
				z[j][o][b] = IloBoolVar(env);
				string varName = "z_" + to_string(j) + "," + to_string(o) + "," + to_string(b);
				z[j][o][b].setName(varName.c_str());
			}
		}
	}

	IloBoolVarArray2 u(env);	// u[Jobs][<StageMachine>]
	for (int j = 0; j < nJobs; ++j) {
		u.add(IloBoolVarArray(env, stageMachine.getCardinality()));
		for (int k = 0; k < stageMachine.getCardinality(); ++k) {
			u[j][k] = IloBoolVar(env);
			string varName = "u_" + to_string(j) + "," + to_string(k);
			u[j][k].setName(varName.c_str());
		}
	}

	IloBoolVarArray3 e(env);	// e[Jobs][Jobs][Stages]
	for (int j1 = 0; j1 < nJobs; ++j1) {
		e.add(IloBoolVarArray2(env));
		for (int j2 = 0; j2 < nJobs; ++j2) {
			e[j1].add(IloBoolVarArray(env, nStages));
			for (int o = 0; o < nStages; ++o) {
				e[j1][j2][o] = IloBoolVar(env);
				string varName = "e_" + to_string(j1) + "," + to_string(j2) + "," + to_string(o);
				e[j1][j2][o].setName(varName.c_str());
			}
		}
	}

	IloNumVarArray2 s(env);		// s[Jobs][Stages]
	for (int j = 0; j < nJobs; ++j) {
		s.add(IloNumVarArray(env, nStages));
		for (int o = 0; o < nStages; ++o) {
			s[j][o] = IloNumVar(env);
			string varName = "s_" + to_string(j) + "," + to_string(o);
			s[j][o].setName(varName.c_str());
		}
	}

	IloNumVarArray2 S(env);		// S[Batches][<StageMachine>]
	for (int b = 0; b < nBatches; ++b) {
		S.add(IloNumVarArray(env, stageMachine.getCardinality()));
		for (int k = 0; k < stageMachine.getCardinality(); ++k) {
			S[b][k] = IloNumVar(env);
			string varName = "S_" + to_string(b) + "," + to_string(k);
			S[b][k].setName(varName.c_str());
		}
	}

	IloNumVarArray T(env, nJobs);		// T[Jobs]
	for (int j = 0; j < nJobs; ++j) {
		T[j] = IloNumVar(env);
		string varName = "T_" + to_string(j);
		T[j].setName(varName.c_str());
	}

	// OBJECTIVE FUNCTION
	IloExpr objTwt(env);
	IloExpr objJobStart(env);
	for (int j = 0; j < nJobs; ++j) {
		objTwt += w[j] * T[j];
		for (int o = 0; o < nStages; ++o) {
			objJobStart += s[j][o];
		}
	}
	IloExpr obj(env);
	obj = omega * objTwt + (1.0 - omega) * objJobStart;		// [JR-2026-Feb-23] added to 2nd summand: "(1.0 - omega) * "

	mod.add(IloMinimize(env, obj));

	// CONSTRAINTS

	// R03: each job is assigned to a fresh or to one of the preformed batches			// forall..
	for (int j = 0; j < nJobs; ++j) {													// j in Jobs
		for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {				// o in Stages_b
			int oBatch = schedule->getBatchingStages()[o] - 1;
			for (int i = 0; i < nProducts; ++i) {
				if (f[j] == (i + 1)) {													// i in Products : i == f[j]
					IloExpr cR03_SumX(env);

					for (int b = 0; b < nJobs; ++b) {									// sum(b in Batches, l in 1..m[o])
						for (int l = 0; l < (*schedule)[oBatch].size(); ++l) {			
							cR03_SumX += x[j][b][constMapStgMac.at({oBatch, l})];
						}
					}
					IloExpr cR03_SumZ(env);
					std::set<int>& B_io = B_io_preformed[i][oBatch];
					for (auto it = B_io.begin(); it != B_io.end(); ++it) {				// sum(b in B_io[i][o])
						cR03_SumZ += z[j][oBatch][*it];
					}
					IloExpr cR03(env);
					cR03 = cR03_SumX + cR03_SumZ; // == 1;
					string conLabel = "R03 (j" + to_string(j) + "o" + to_string(oBatch) + "i" + to_string(i) + ")";
					mod.add(cR03 == 1).setName(conLabel.c_str());
				}
			}
		}
	}

	// R04: capacity																	// forall..
	for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {					// o in Stages_b
		int oBatch = schedule->getBatchingStages()[o] - 1;
		for (int b = 0; b < nBatches; ++b) {											// b in Batches
			for (int l = 0; l < (*schedule)[oBatch].size(); ++l) {						// l in 1..m[o]
				IloExpr cR04_SumX(env);
				for (int j = 0; j < nJobs; ++j) {										// sum(j in Jobs)
					cR04_SumX += x[j][b][constMapStgMac.at({oBatch, l})];
				}
				IloExpr cR04(env);
				cR04 = cR04_SumX; // <= B[o];
				string conLabel = "R04 (o" + to_string(oBatch) + "b" + to_string(b) + "l" + to_string(l) + ")";
				mod.add(cR04 <= B[o]).setName(conLabel.c_str());
			}
		}
	}

	// R05: capacity of preformed batches												// forall..
	for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {					// o in Stages_b
		int oBatch = schedule->getBatchingStages()[o] - 1;

		for (auto it = batchIdxMac[oBatch].begin(); it != batchIdxMac[oBatch].end(); ++it) {
			IloExpr cR05_SumZ(env);
			for (int j = 0; j < nJobs; ++j) {											// sum(j in Jobs)

				cR05_SumZ += z[j][oBatch][it->first.second];
			}
			IloExpr cR05(env);
			cR05 = cR05_SumZ + B_iob[it->first.first][oBatch][it->first.second];

			string conLabel = "R05 (o" + to_string(oBatch) + "b" + to_string(it->first.second) + "i" + to_string(it->first.first) + ")";
			mod.add(cR05 <= B[o]).setName(conLabel.c_str());
		}
	}

	// R06: unique family of a batch													// forall..
	for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {					// o in Stages_b
		int oBatch = schedule->getBatchingStages()[o] - 1;
		for (int b = 0; b < nBatches; ++b) {												// b in Batches
			for (int l = 0; l < (*schedule)[oBatch].size(); ++l) {						// l in 1..m[o]
				IloExpr cR06_SumY(env);
				for (int i = 0; i < nProducts; ++i) {									// sum(i in Products)
					cR06_SumY += y[i][b][constMapStgMac.at({oBatch, l})];
				}
				IloExpr cR06(env);
				cR06 = cR06_SumY; // == 1;
				string conLabel = "R06 (o" + to_string(oBatch) + "b" + to_string(b) + "l" + to_string(l) + ")";
				mod.add(cR06 == 1).setName(conLabel.c_str());
			}
		}
	}

	// R07: matching families of jobs and batches										// forall..
	for (int j = 0; j < nJobs; ++j) {													// j in Jobs
		for (int b = 0; b < nBatches; ++b) {											// b in Batches
			for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {			// o in Stages_b
				int oBatch = schedule->getBatchingStages()[o] - 1;
				for (int l = 0; l < (*schedule)[oBatch].size(); ++l) {					// l in 1..m_o[o]
					for (int i = 0; i < nProducts; ++i) {
						if (f[j] == (i + 1)) {											// i in Products : i == f[j]
							IloExpr cR07(env);
							cR07 = x[j][b][constMapStgMac.at({oBatch, l})]; // <= y[i][b][oBatch, l];
							string conLabel = "R07 (j" + to_string(j) + "b" + to_string(b) + "o" + to_string(oBatch) + "l" + to_string(l) + "i" + to_string(i) + ")";
							mod.add(cR07 <= y[i][b][constMapStgMac.at({oBatch, l})]).setName(conLabel.c_str());

						}
					}
				}
			}
		}
	}

	// R08: each operation is assigned to one machine										// forall..
	for (int j = 0; j < nJobs; ++j) {														// j in Jobs
		for (int o = 0; o < schedule->getDiscreteStages().size(); ++o) {					// o in Stages_1
			int oSingle = schedule->getDiscreteStages()[o] - 1;
			IloExpr cR08_SumU(env);
			for (int l = 0; l < (*schedule)[oSingle].size(); ++l) {							// sum(l in 1..m_o[o])
				cR08_SumU += u[j][constMapStgMac.at({oSingle, l})];
			}
			IloExpr cR08(env);
			cR08 = cR08_SumU; // == 1;
			string conLabel = "R08 (j" + to_string(j) + "o" + to_string(oSingle) + ")";
			mod.add(cR08 == 1).setName(conLabel.c_str());
		}
	}

	// R09: release dates of jobs									// forall..
	for (int j = 0; j < nJobs; ++j) {								// j in Jobs
		IloExpr cR09(env);
		cR09 = s[j][0]; // >= r[j];
		string conLabel = "R09 (j" + to_string(j) + ")";
		mod.add(cR09 >= r[j]).setName(conLabel.c_str());
	}

	// R10: release dates of machines													// forall..
	for (int j = 0; j < nJobs; ++j) {													// j in Jobs
		for (int o = 0; o < schedule->getDiscreteStages().size(); ++o) {				// o in Stages_1
			int oSingle = schedule->getDiscreteStages()[o] - 1;
			for (int l = 0; l < (*schedule)[oSingle].size(); ++l) {					// l in 1..m_o[o]
				IloExpr cR10(env);
				cR10 = s[j][oSingle] + (rm[oSingle][l] * (1 - u[j][constMapStgMac.at({oSingle, l})])); // >= rm[oSingle][l];
				string conLabel = "R10 (j" + to_string(j) + "o" + to_string(oSingle) + "l" + to_string(l) + ")";
				mod.add(cR10 >= rm[oSingle][l]).setName(conLabel.c_str());
			}
		}
	}

	// R11: release dates of batching machines											// forall..
	for (int j = 0; j < nJobs; ++j) {													// j in Jobs
		for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {				// o in Stages_b
			int oBatch = schedule->getBatchingStages()[o] - 1;
			for (int b = 0; b < nBatches; ++b) {											// b in Batches
				for (int l = 0; l < (*schedule)[oBatch].size(); ++l) {					// l in 1..m_o[o]
					IloExpr cR11(env);
					cR11 = S[b][constMapStgMac.at({oBatch, l})]; // >= rm[oBatch][l];
					string conLabel = "R11 (j" + to_string(j) + "o" + to_string(oBatch) + "b" + to_string(b) + "l" + to_string(l) + ")";
					mod.add(cR11 >= rm[oBatch][l]).setName(conLabel.c_str());
				}
			}
		}
	}

	// R12: processing times of batches													// forall..
	for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {					// o in Stages_b
		int oBatch = schedule->getBatchingStages()[o] - 1;
		for (int b = 0; b < (nBatches - 1); ++b) {										// b in 1..(n-1)
			for (int l = 0; l < (*schedule)[oBatch].size(); ++l) {						// l in 1..m_o[o]
				IloExpr cR12_SumY(env);
				for (int i = 0; i < nProducts; ++i) {									// sum(i in Products)
					cR12_SumY += p[i][oBatch] * y[i][b][constMapStgMac.at({oBatch, l})];
				}
				IloExpr cR12(env);
				cR12 = S[b][constMapStgMac.at({oBatch, l})] + cR12_SumY; // <= S[b + 1][oBatch, l];
				string conLabel = "R12 (o" + to_string(oBatch) + "b" + to_string(b) + "l" + to_string(l) + ")";
				mod.add(cR12 <= S[b + 1][constMapStgMac.at({oBatch, l})]).setName(conLabel.c_str());
			}
		}
	}

	// R13: precedence of jobs (single processing)										// forall..
	for (int j = 0; j < nJobs; ++j) {
		for (int k = 0; k < nJobs; ++k) {
			if (j != k) {																// j, k : j != k
				for (int i = 0; i < nProducts; ++i) {
					if (f[j] == (i + 1)) {												// i in Products : i == f[j]
						for (int o = 0; o < schedule->getDiscreteStages().size(); ++o) {
							int oSingle = schedule->getDiscreteStages()[o] - 1;
							IloExpr cR13(env);
							cR13 = s[j][oSingle] + p[i][oSingle] + (G * (e[j][k][oSingle] - 1)); // <= s[k][oSingle];
							string conLabel = "R13 (j" + to_string(j) + "k" + to_string(k) + "i" + to_string(i) + ")";
							mod.add(cR13 <= s[k][oSingle]).setName(conLabel.c_str());
						}
					}
				}
			}
		}
	}

	// R14: precedence of operations of the same job (across stages)					// forall..
	for (int j = 0; j < nJobs; ++j) {													// j in Jobs
		for (int i = 0; i < nProducts; ++i) {
			if (f[j] == (i + 1)) {														// i in Products : i == f[j]
				for (int o = 0; o < (nStages - 1); ++o) {								// o in 1..(m-1)
					IloExpr cR14(env);
					cR14 = s[j][o] + p[i][o]; 
					string conLabel = "R14 (j" + to_string(j) + "i" + to_string(i) + "o" + to_string(o) + ")";
					mod.add(cR14 <= s[j][o + 1]).setName(conLabel.c_str());
				}
			}
		}
	}

	// R15: precedence of different jobs on the same machine and stage					// forall..
	for (int j = 0; j < nJobs; ++j) {
		for (int k = 0; k < nJobs; ++k) {
			if (j != k) {																// j, k in Jobs : j != k
				for (int o = 0; o < schedule->getDiscreteStages().size(); ++o) {		// o in Stages_1
					int oSingle = schedule->getDiscreteStages()[o] - 1;
					for (int l = 0; l < (*schedule)[oSingle].size(); ++l) {				// l in 1..m_o[o]
						IloExpr cR15(env);
						cR15 = u[j][constMapStgMac.at({oSingle, l})] + u[k][constMapStgMac.at({oSingle, l})] - e[j][k][oSingle] - e[k][j][oSingle]; // <= 1;
						string conLabel = "R15 (j" + to_string(j) + "k" + to_string(k) + "o" + to_string(oSingle) + "l" + to_string(l) + ")";
						mod.add(cR15 <= 1).setName(conLabel.c_str());
					}
				}
			}
		}
	}

	// R16: only one of two jobs can precede the other								// forall..
	for (int j = 0; j < nJobs; ++j) {
		for (int k = 0; k < nJobs; ++k) {											// j, k in Jobs
			for (int o = 0; o < schedule->getDiscreteStages().size(); ++o) {		// o in Stages_1
				int oSingle = schedule->getDiscreteStages()[o] - 1;
				for (int l = 0; l < (*schedule)[oSingle].size(); ++l) {				// l in 1..m_o[o]
					IloExpr cR16(env);
					cR16 = e[j][k][oSingle] + e[k][j][oSingle]; 
					string conLabel = "R16 (j" + to_string(j) + "k" + to_string(k) + "o" + to_string(oSingle) + "l" + to_string(l) + ")";
					mod.add(cR16 <= 1).setName(conLabel.c_str());
				}
			}
		}
	}

	// R17: matching start times of jobs and batches								// forall..
	for (int j = 0; j < nJobs; ++j) {												// j in Jobs
		for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {						// o in Stages_b
			int oBatch = schedule->getBatchingStages()[o] - 1;
			for (int b = 0; b < nBatches; ++b) {									// b in Batches
				for (int l = 0; l < (*schedule)[oBatch].size(); ++l) {					// l in 1..m_o[o]
					IloExpr cR17(env);
					cR17 = s[j][oBatch] + (G * (1 - x[j][b][constMapStgMac.at({oBatch, l})])); // >= S[b][oBatch, l];
					string conLabel = "R17 (j" + to_string(j) + "o" + to_string(oBatch) + "b" + to_string(b) + "l" + to_string(l) + ")";
					mod.add(cR17 >= S[b][constMapStgMac.at({oBatch, l})]).setName(conLabel.c_str());
				}
			}
		}
	}

	// R18: matching start times of jobs and batches									// forall..
	for (int j = 0; j < nJobs; ++j) {												// j in Jobs
		for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {						// o in Stages_b
			int oBatch = schedule->getBatchingStages()[o] - 1;
			for (int b = 0; b < nBatches; ++b) {									// b in Batches
				for (int l = 0; l < (*schedule)[oBatch].size(); ++l) {					// l in 1..m_o[o]
					IloExpr cR18(env);
					cR18 = s[j][oBatch] + (G * (x[j][b][constMapStgMac.at({oBatch, l})] - 1)); // <= S[b][oBatch, l];
					string conLabel = "R18 (j" + to_string(j) + "o" + to_string(oBatch) + "b" + to_string(b) + "l" + to_string(l) + ")";
					mod.add(cR18 <= S[b][constMapStgMac.at({oBatch, l})]).setName(conLabel.c_str());
				}
			}
		}
	}

	// R19: matching start times of jobs and batches for preformed batches				// forall..
	for (int j = 0; j < nJobs; ++j) {												// j in Jobs
		for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {						// o in Stages_b
			int oBatch = schedule->getBatchingStages()[o] - 1;
			for (int i = 0; i < nProducts; ++i) {
				if (f[j] == (i + 1)) {													// i in Products : f[j] == i
					std::set<int>& B_io = B_io_preformed[i][oBatch];
					for (auto it = B_io.begin(); it != B_io.end(); ++it) {				// b in B_io[i][o]
						IloExpr cR19(env);
						cR19 = s[j][oBatch] + G * (1 - z[j][oBatch][*it]); // >= S_iob[i][oBatch][*it];
						string conLabel = "R19 (j" + to_string(j) + "o" + to_string(oBatch) + "i" + to_string(i) + "b" + to_string(*it - 1) + ")";
						mod.add(S_iob[i][oBatch][*it] <= cR19).setName(conLabel.c_str());
					}
				}
			}
		}
	}

	// R20: matching start times of jobs and batches for preformed batches				// forall..
	for (int j = 0; j < nJobs; ++j) {												// j in Jobs
		for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {						// o in Stages_b
			int oBatch = schedule->getBatchingStages()[o] - 1;
			for (int i = 0; i < nProducts; ++i) {
				if (f[j] == (i + 1)) {													// i in Products : f[j] == i
					std::set<int>& B_io = B_io_preformed[i][oBatch];
					for (auto it = B_io.begin(); it != B_io.end(); ++it) {				// b in B_io[i][o]
						IloExpr cR20(env);
						cR20 = s[j][oBatch] + G * (z[j][oBatch][*it] - 1);
						string conLabel = "R20 (j" + to_string(j) + "o" + to_string(oBatch) + "i" + to_string(i) + "b" + to_string(*it - 1) + ")";
						mod.add(cR20 <= S_iob[i][oBatch][*it]).setName(conLabel.c_str());
					}
				}
			}
		}
	}

	// R21: time constraints														// forall..
	for (int j = 0; j < nJobs; ++j) {												// j in Jobs
		for (int i = 0; i < nProducts; ++i) {
			if (f[j] == (i + 1)) {													// i in Products, f[j] == i
				for (int o = 0; o < nStages; ++o) {
					for (int p = 0; p < nStages; ++p) {								// o, p in Stages
						IloExpr cR21(env);
						cR21 = s[j][p] - s[j][o];
						string conLabel = "R21 (j" + to_string(j) + "i" + to_string(i) + "o" + to_string(o) + "p" + to_string(p) + ")";
						mod.add(cR21 <= tc[i][o][p]).setName(conLabel.c_str());
					}
				}
			}
		}
	}

	// R22: tardiness																	// forall..
	for (int j = 0; j < nJobs; ++j) {												// j in Jobs
		for (int i = 0; i < nProducts; ++i) {
			if (f[j] == (i + 1)) {														// i in Products : f[j] == i
				IloExpr cR22(env);
				cR22 = s[j][m - 1] + p[i][m - 1] - d[j];
				string conLabel = "R22 (j" + to_string(j) + "i" + to_string(i) + ")";
				mod.add(cR22 <= T[j]).setName(conLabel.c_str());
			}
		}
	}

	// set timelimit and solve
	IloCplex cplex(mod);
	cplex.setOut(env.getNullStream());				// GO-LIVE: uncomment
	//cplex.exportModel("milp.lp");					// GO-LIVE: comment
	cplex.setParam(IloCplex::EpGap, params->presetMipGap);	// [JR-2026-Jan-23]
	cplex.setParam(IloCplex::TiLim, cplexTilim);	

	try {
		if (!cplex.solve()) {
			std::stringstream errorMsg;
			IloAlgorithm::Status status = cplex.getStatus();
			errorMsg << "CPLEX failed to solve the MILP. Status: " << cplex.getStatus() << ". ";
			if (status == IloAlgorithm::Infeasible) {
				errorMsg << "The model is infeasible. ";
				// TODO Use conflict refiner to identify conflicting constraints
			}
			else if (status == IloAlgorithm::Unbounded) {
				errorMsg << "The model is unbounded. ";
			}
			else if (status == IloAlgorithm::InfeasibleOrUnbounded) {
				errorMsg << "The model is infeasible or unbounded. ";
			}
			else if (status == IloAlgorithm::Error) {
				errorMsg << "An error occurred during optimization. Error code: " << cplex.getCplexStatus() << ". ";
			}
			TCB::logger.Log(Error, errorMsg.str());
			cout << errorMsg.str() << endl;
		}
	}
	catch (IloException& e) {
		TCB::logger.Log(Error, e.getMessage());
		env.end();
	}

	params->relMipGap = cplex.getMIPRelativeGap();	// [JR-2026-Jan-21] output rel MIP gap of last iteration (only really interesting for exact solutions with nDash == n

	cout << ".";
	// construct schedule from solution
	for (int stg = 0; stg < nStages; ++stg) {
		if (schedule->getCapAtStageIdx(stg) <= 1) {
			// no batching at this stages
			for (int j = 0; j < nJobs; ++j) {
				pBat newSingleBatch = make_unique<Batch>(schedule->getCapAtStageIdx(stg));	// assumption: parallel identical machines
				newSingleBatch->addOp(consideredJobs[j]->getOpPtr(stg));
				int macIdx;
				for (int mac = 0; mac < (*schedule)[stg].size(); ++mac) {
					IloNum val = cplex.getValue(u[j][constMapStgMac.at({stg, mac})]);		// u[j][stgMac]	: machine assignment
					if (val >= 1 - TCB::precision) {
						macIdx = mac;
					}
				}
				IloNum start = cplex.getValue(s[j][stg]);
				if (!(*schedule)[stg][macIdx].addBatch(move(newSingleBatch), start)) { 
					TCB::logger.Log(Error, "Problem in SolverMILP::solve(): infeasible op assignment " + to_string((j + 1)) + " at stage " + to_string((stg + 1)) + ": " + to_string(start));
				}
			}
		}
		else {
			// add to preformed batches (dvar z[j][o][b])
			for (int j = 0; j < nJobs; ++j) {
				for (auto it = batchIdxMac[stg].begin(); it != batchIdxMac[stg].end(); ++it) { //vector[stg], map< pair<product, batchIdx>, macIdx>
					try {
						IloNum zVal = cplex.getValue(z[j][stg][it->first.second]);
						if (zVal >= 1 - TCB::precision) {
							if (!(*schedule)[stg][it->second.first][it->second.second].addOp(&(*consideredJobs[j])[stg])) {
								stringstream errorMsg;
								errorMsg << "Problem in SolverMILP::solve(): infeasible op assignment " << consideredJobs[j]->getId() << " at stage " << (stg + 1);
								TCB::logger.Log(Error, errorMsg.str());
							}
						}
					}
					catch (const IloCplex::Exception exc) {
						std::stringstream errorMsg;
						errorMsg << "CPLEX exception caught: " << exc.getMessage();
						TCB::logger.Log(Error, errorMsg.str());
					}
				}
			}

			// form fresh batches (dvar x[j][b][o][l])
			//int batchingStageIndex = 0;
			//const vector<int>& stagesB = schedule->getBatchingStages();
			//for (int i = 0; i < stagesB.size(); ++i) {
			//	if (stagesB[i] == (stg + 1)) {
			//		batchingStageIndex = i;
			//		break;
			//	}
			//}
			int capacity = schedule->getCapAtStageIdx(stg);
			vector<vector<pBat>> newBatches = vector<vector<pBat>>((*schedule)[stg].size());
			for (size_t l = 0; l < ((*schedule)[stg].size()); ++l) {
				for (size_t b = 0; b < nBatches; ++b) {
					newBatches[l].push_back(make_unique<Batch>(capacity));
				}
			}
				
			/*for (int o = 0; o < (*schedule)[stg].size(); ++o) {
				newBatches[o].resize(nBatches);
			}*/
			/*for (int mac = 0; mac < (*schedule)[stg].size(); ++mac) {
				for (int b = 0; b < nBatches; ++b) {
					newBatches[mac][b] = make_unique<Batch>(Batch(capacity));
				}
			}*/

			for (int b = 0; b < nJobs; ++b) {
				for (int j = 0; j < nJobs; ++j) {
					for (int mac = 0; mac < (*schedule)[stg].size(); ++mac) {
						IloNum xJBOL = cplex.getValue(x[j][b][constMapStgMac.at({stg, mac})]);
						if (xJBOL >= 1 - TCB::precision) {
							if (!newBatches[mac][b]->addOp(consideredJobs[j]->getOpPtr(stg))) {
								throw ExcSched("Could not add op to batch (invalid solution?), SolverMILP::initModel()");
							}
						}
					}
				}
			}

			// assign batches to machines
			for (int b = 0; b < nJobs; ++b) {
				for (int mac = 0; mac < (*schedule)[stg].size(); ++mac) {
					for (int j = 0; j < nJobs; ++j) {
						IloNum xJBOL = cplex.getValue(x[j][b][constMapStgMac.at({stg, mac})]);
						if (xJBOL >= 1 - TCB::precision) {
							IloNum sBOL = cplex.getValue(S[b][mapStgMac[{stg, mac}]]);
							IloNum sJO = cplex.getValue(s[j][stg]);
							try {
								if (newBatches[mac][b] != nullptr) {
									if (!(*schedule)[stg][mac].addBatch(move(newBatches[mac][b]), (double)sBOL)) {
										throw ExcSched("");
									}	// sBOL and sJO should be equal but the tolerance is too high (big G too large?)
								}	 
							}
							catch (ExcSched& e) {
								TCB::logger.Log(Error, e.getMessage());
								break; // batch must only be assigned once
							}
						}
					}
				}
			}
		}
	}

	while (consideredJobs.size() > 0) {
		schedule->markAsScheduled(move(consideredJobs[0]));
		consideredJobs.erase(consideredJobs.begin());
	}

	//double bestObjectiveValue = cplex.getBestObjValue();
	env.end();
	return schedule->getTWT();	// bestObjectiveValue; // integrated objective function
}
double Solver_MILP::solveDecompJobBasedMILP(Schedule* schedule, int nDash, int cplexTilim) {
	double twt = FLT_MAX;
	int iteration = 1;	// DEBUGGING
	while (schedule->getN() > 0) {
		// SOLVE with nDash Jobs
		twt = solveJobBasedMILP(schedule, nDash, cplexTilim);
		//schedule->saveJsonFactory("MILP_it" + to_string(iteration));	// DEBUGGING
		iteration++;	// DEBUGGING
	}
	return twt;
}
double Solver_MILP::solveDecompJobBasedMILP(Schedule* schedule, int nDash, int cplexTilim, prioRuleKappa<pJob> rule, double kappa)
{	
	double t = schedule->getMinMSP(0);
	double twt = numeric_limits<double>::max();
	
	while (schedule->getN() > 0) {
		// DYNAMIC SORTING
		if (rule != nullptr) {
			schedule->sortUnscheduled(rule, kappa);
		}

		// SOLVE
		twt = solveJobBasedMILP(schedule, nDash, cplexTilim);

		// UPDATE PARAMETERS
		t = schedule->getMinMSP(0);
	}
	return twt;
}
double Solver_MILP::solveDecompJobBasedDynamicSortingMILP(Schedule* schedule, int nDash, int cplexTilim, prioRule<pJob> initRule, prioRuleKappa<pJob> dynRule, double kappa) {
	// OBTAIN INITIAL WAITING TIMES FROM A SCHEDULE WITH LIST SCHEDULING AND SORTING BY initRule (DEFAULT: EDD)
	unique_ptr<Schedule> initSched = schedule->clone();
	initSched->lSchedJobsWithSorting(initRule);
	initSched->updateWaitingTimes();
	schedule->mimicWaitingTimes(initSched.get());

	while (schedule->getN() > 0) {
		schedule->sortUnscheduled(dynRule, kappa);
		solveJobBasedMILP(schedule, nDash, cplexTilim);

		// DYNAMIC SORTING AND UPDATE WAITING TIMES
		unique_ptr<Schedule> tempSched = schedule->clone();
		tempSched->lSchedJobsWithSorting(dynRule, kappa);
		schedule->mimicWaitingTimes(tempSched.get());
	}

	return schedule->getTWT();
}
double Solver_MILP::solveDecompJobBasedDynamicSortingGridMILP(Schedule* schedule, int nDash, int cplexTilim, prioRule<pJob> initRule, prioRuleKappa<pJob> dynRule) {
	vector<double> kappas = getDoubleGrid(schedParams->kappaLow, schedParams->kappaHigh, schedParams->kappaStep);
	// OBTAIN INITIAL WAITING TIMES FROM A SCHEDULE WITH LIST SCHEDULING AND SORTING BY initRule (DEFAULT: EDD)
	unique_ptr<Schedule> initSched = schedule->clone();
	initSched->lSchedJobsWithSorting(initRule);
	initSched->updateWaitingTimes();

	// DEBUGGING
	//initSched->saveJsonFactory("INITIAL_SCHEDULE_EDD");
	//schedule->sortUnscheduled(sortJobsDebugging);	// DEBUGGING 23.06.2026


	schedule->mimicWaitingTimes(initSched.get());

	while (schedule->getN() > 0) {
		// DYNAMIC SORTING AND UPDATE WAITING TIMES 
		unique_ptr<Schedule> tempSched = schedule->clone();
		double bestKappa = tempSched->lSchedJobsWithSorting(dynRule, kappas);	// GRID SEARCH
		tempSched->updateWaitingTimes();
		schedule->mimicWaitingTimes(tempSched.get());	
		schedule->sortUnscheduled(dynRule, bestKappa);

		// ASSIGN THE NEXT NDASH JOBS

		solveJobBasedMILP(schedule, nDash, cplexTilim);

		// DEBUGGING
		//schedule->saveJsonFactory("PARTIAL_SCHEDULE_AFTER_ITERATION");
		//int breakpoint = 666;

	}

	
	return schedule->getTWT();
}

double Solver_MILP::solveGapAwareMILP(Schedule* schedule, int nDash, int cplexTilim)
{
	//cout << "DEVELOPMENT IN PROGRESS: CONSIDER [JR-2026-Okt-26] TAG" << endl;

	if (schedule->getProblem() == nullptr) throw ExcSched("solveDecompositionMILP missing problem reference");

	vector<pJob> consideredJobs = vector<pJob>();
	for (size_t j = 0; j < nDash; ++j) {
		try {
			consideredJobs.push_back(move(schedule->get_pJob(0)));
			//cout << consideredJobs[0]->getId();	// DEBUGGING
		}
		catch (const out_of_range& ex) {
			break;	// all jobs considered
		}
	}


	

	IloEnv env;

	// PARAMETERS
	IloNum omega = this->omega;	// TODO set omega

	IloNum G = schedule->getProblem()->getG();
	IloInt m = schedule->size();

	int nProducts = schedule->getProblem()->getF();
	int nStages = schedule->size();
	int nJobs = consideredJobs.size();
	int nBatches = (*schedule).getNumberOfScheduledJobs() + consideredJobs.size();	// new and (max. number of) already formed batches
	
	
	
	IloIntTupleSet stageMachine(env, 2);
	for (int o = 0; o < nStages; ++o) {
		for (int l = 0; l < (*schedule)[o].size(); ++l) {
			stageMachine.add(IloIntArray(env, 2, o, l));
		}
	}
	map<pair<int, int>, int> mapStgMac = map<pair<int, int>, int>();	// for access to tuple "stageMachine"
	int index = 0;
	for (IloIntTupleSetIterator it(env, stageMachine); it.ok(); ++it) {
		mapStgMac.insert(make_pair(make_pair((*it)[0], (*it)[1]), index));
		index++;
	}
	const map<pair<int, int>, int> constMapStgMac = mapStgMac;

	
	IloIntArray3 B(env);		// capacity[stageMachine][batches]
	// the following parameters consider batches formed in previous iterations
	IloNumArray3 batStart(env);	// earliest start of a batch
	IloNumArray3 batEnd(env);	// latest start of a batch
	IloNumArray3 batDur(env);	// processing time of a batch
	IloNumArray3 batFam(env);	// family (product) of a batch


	// CALCULATE UPPER BOUND ON BATCH COMPLETION
	int rMax = 0;
	int pMax = 0;
	int rawP = 0;
	for (int j = 0; j < consideredJobs.size(); j++) {
		if (consideredJobs[j]->getR() > rMax) {
			rMax = consideredJobs[j]->getR();
		}


		for (int o = 0; o < nStages; ++o) {
			if (consideredJobs[j]->getP(o) > pMax) {
				pMax = consideredJobs[j]->getP(o);
			}
			rawP += consideredJobs[j]->getP(o);
		}
	}
	int maxMSP = 0;
	for (int l = 0; l < (*schedule)[nStages - 1].size(); ++l) {
		if ((*schedule)[nStages - 1][l].getMSP() > maxMSP) {
			maxMSP = (*schedule)[nStages - 1][l].getMSP();
		}
	}
	int latestBatchStart = max(rMax, maxMSP) + rawP + pMax;	 

	// variables for schedule construction
	vector<vector<vector<Batch*>>> existingBatches = vector<vector<vector<Batch*>>>(nStages);

	for (int o = 0; o < m; ++o) {
		B.add(IloIntArray2(env));
		batStart.add(IloNumArray2(env));
		batEnd.add(IloNumArray2(env));
		batDur.add(IloNumArray2(env));
		batFam.add(IloNumArray2(env));
		existingBatches[o] = vector < vector<Batch*>>((*schedule)[o].size());

		for (int m = 0; m < (*schedule)[o].size(); ++m) {
			B[o].add(IloIntArray(env));
			batStart[o].add(IloNumArray(env));
			batEnd[o].add(IloNumArray(env));
			batDur[o].add(IloNumArray(env));
			batFam[o].add(IloNumArray(env));

			int nAlreadyScheduledBatches = 0;
			Machine* mac = &(*schedule)[o][m];
			// previously scheduled batches
			for (int b = 0; b < mac->size(); b++) {
				Batch* bat = &(*mac)[b];
				B[o][m].add(IloInt(bat->getAvailableCap()));
				batStart[o][m].add(IloNum(bat->getStart()));
				batEnd[o][m].add(IloNum(bat->getC()));
				batDur[o][m].add(IloNum(bat->getP()));
				batFam[o][m].add(IloNum(bat->getF()));

				existingBatches[o][m].push_back(bat);

				++nAlreadyScheduledBatches;
			}
			
			// additional batches 
			for (int b = nAlreadyScheduledBatches; b < nBatches; ++b) {
				B[o][m].add(IloInt(schedule->getCapAtStageIdx(o)));
				batStart[o][m].add(IloInt(0));
				batEnd[o][m].add(IloInt(latestBatchStart));
				batDur[o][m].add(IloInt(0));
				batFam[o][m].add(IloInt(0));
			}
		}
	}

	IloNumArray d(env);		// d[Jobs] job due dates
	IloNumArray r(env);		// r[Jobs] job ready times
	IloNumArray w(env);		// w[Jobs] job weights
	IloIntArray f(env);		// f[Jobs] jobs´ products or families

	for (int j = 0; j < nJobs; ++j) {
		d.add(consideredJobs[j]->getD());
		r.add(consideredJobs[j]->getR());
		w.add(consideredJobs[j]->getW());
		f.add(consideredJobs[j]->getF());
	}

	IloNumArray2 p(env);	// p[Products][Stages] processing times
	for (int i = 0; i < nProducts; ++i) {
		p.add(IloNumArray(env, nStages));
		for (int o = 0; o < nStages; ++o) {
			p[i][o] = schedule->getProblem()->getProduct(i)->getP(o);
		}
	}

	IloNumArray3 tc(env);	// tc[Products][Stages][Stages] time constraints
	for (int i = 0; i < nProducts; ++i) {
		tc.add(IloNumArray2(env));
		for (int o1 = 0; o1 < nStages; ++o1) {
			tc[i].add(IloNumArray(env, nStages));
			const vector<pair<int, double>>& tcMaxFwd = schedule->getProblem()->getProduct(i)->getTcMaxFwd(o1);
			for (size_t o2 = 0; o2 < nStages; ++o2) {
				double tcDuration = 999999;	// TODO define INF
				for (auto it = tcMaxFwd.begin(); it != tcMaxFwd.end(); ++it) {
					if (it->first == o2) {
						tcDuration = it->second;
					}
				}
				tc[i][o1][o2] = (IloNum)tcDuration;
			}
		}
	}
	IloModel mod(env);

	// DECISION VARIABLES
	IloBoolVarArray3 x(env);	// x[<StageMachine>][Jobs][Batches]		: JOB 2 BATCH ASSIGNMENT
	for (int k = 0; k < stageMachine.getCardinality(); ++k) {
		x.add(IloBoolVarArray2(env));
		for (int j = 0; j < nJobs; ++j) {
			x[k].add(IloBoolVarArray(env, nBatches));
			for (int b = 0; b < nBatches; ++b) {
				x[k][j][b] = IloBoolVar(env);
				string varName = "x_" + to_string(k) + "," + to_string(j) + "," + to_string(b);
				x[k][j][b].setName(varName.c_str());
			}
		}
	}

	IloBoolVarArray3 y(env);	// y[<StageMachine>][Batches][Batches]	: PRECEDENCE OF BATCHES
	for (int k = 0; k < stageMachine.getCardinality(); ++k) {	
		y.add(IloBoolVarArray2(env));
		for (int b = 0; b < nBatches; ++b) {
			y[k].add(IloBoolVarArray(env, nBatches));
			for (int c = 0; c < nBatches; ++c) {	
				y[k][b][c] = IloBoolVar(env);
				string varName = "y_" + to_string(k) + "," + to_string(b) + "," + to_string(c);
				y[k][b][c].setName(varName.c_str());
			}
		}
	}

	IloBoolVarArray3 z(env);	// z[<StageMachine>][Batches][Products]	: BATCH FAMILY (PRODUCT)
	for (int k = 0; k < stageMachine.getCardinality(); ++k) {
		z.add(IloBoolVarArray2(env));
		for (int b = 0; b < nBatches; ++b) {
			z[k].add(IloBoolVarArray(env, nProducts));
			for (int i = 0; i < nProducts; ++i) {
				z[k][b][i] = IloBoolVar(env);
				string varName = "z_" + to_string(k) + "," + to_string(b) + "," + to_string(i);
				z[k][b][i].setName(varName.c_str());
			}
		}
	}

	IloNumVarArray2 S(env);		// S[<StageMachine>][Batches]			: BATCH START
	for (int k = 0; k < stageMachine.getCardinality(); ++k) {
		S.add(IloNumVarArray(env, nBatches));
		for (int b = 0; b < nBatches; ++b) {
			S[k][b] = IloNumVar(env);
			string varName = "S_" + to_string(k) + "," + to_string(b);
			S[k][b].setName(varName.c_str());
		}
	}

	IloNumVarArray2 C(env);		// C[Jobs][Stages]			: OPERATION COMPLETION
	for (int j = 0; j < nJobs; ++j) {
		C.add(IloNumVarArray(env, nStages));
		for (int o = 0; o < nStages; ++o) {
			C[j][o] = IloNumVar(env);
			string varName = "C_" + to_string(j) + "," + to_string(o);
			C[j][o].setName(varName.c_str());
		}
	}

	IloNumVarArray T(env, nJobs);		// T[Jobs]						: JOB TARDINESS
	for (int j = 0; j < nJobs; ++j) {
		T[j] = IloNumVar(env);
		string varName = "T_" + to_string(j);
		T[j].setName(varName.c_str());
	}

	// OBJECTIVE FUNCTION
	IloExpr objTwt(env);
	IloExpr objOpComp(env);
	for (int j = 0; j < nJobs; ++j) {
		objTwt += w[j] * T[j];
		for (int o = 0; o < nStages; ++o) {
			objOpComp += C[j][o];
		}
	}
	IloExpr obj(env);
	obj = omega * objTwt + (1.0 - omega) * objOpComp;		

	mod.add(IloMinimize(env, obj));

// CONSTRAINTS

//R02: Every operation is assigned once
//R02: forall(j in jobs, o in stages) sum(b in batches, l in 1..m_o[o]) X[<o, l>][j][b] == 1;
	for (int j = 0; j < nJobs; ++j) {									// forall(j in jobs, 
		for (int o = 0; o < nStages; ++o) {								//   o in stages)
			IloExpr cR02(env);											//   sum(
			for (int b = 0; b < nBatches; ++b) {						//     b in batches,
				for (int l = 0; l < (*schedule)[o].size(); ++l) {		//	   l in machines)
					cR02 += x[constMapStgMac.at({ o, l })][j][b];
				}
			}
			string conLabel = "R02 (j" + to_string(j) + "o" + to_string(o) + ")";
			mod.add(cR02 == 1).setName(conLabel.c_str());
		}
	}

//R03: Batch capacity
//R03: forall(sm in stageMachine, b in batches) sum(j in jobs) X[sm][j][b] <= B[sm][b];
	for (int o = 0; o < nStages; ++o) {									// forall(
		for (int l = 0; l < (*schedule)[o].size(); ++l) {				//   sm in stagemachine,
			for (int b = 0; b < nBatches; ++b) {						//	 b in batches)
				IloExpr cR03(env);
				for (int j = 0; j < nJobs; ++j) {						//   sum(j in jobs)
					cR03 += x[constMapStgMac.at({ o, l })][j][b];
				}
				string conLabel = "R03 (o" + to_string(o) + "l" + to_string(l) + "b" + to_string(b) + ")";
				mod.add(cR03 <= B[o][l][b]);							// B[sm][b]
			}
		}
	}

//R04: Compatibility: batch formation
//R04: forall(j in jobs, i in families, sm in stageMachine, b in batches) e[j][i] * X[sm][j][b] <= Z[sm][b][i];
	for (int j = 0; j < nJobs; ++j) {												// forall(j in jobs,
		for (int i = 0; i < nProducts; ++i) {										//   i in products,
			for (int o = 0; o < nStages; ++o) {										//   
				for (int l = 0; l < (*schedule)[o].size(); ++l) {					//   sm in stageMachine
					for (int b = 0; b < nBatches; ++b) {							//   b in batches
						int e = i == f[j]-1;
						IloExpr cR04(env);
						cR04 = e * x[constMapStgMac.at({ o, l })][j][b];	
						mod.add(cR04 <= z[constMapStgMac.at({ o, l })][b][i]);	// e[j][i] * x <= z
					}
				}
			}
		}
	}

//R05: Compatibility: batch family
	for (int o = 0; o < nStages; ++o) {
		for (int l = 0; l < (*schedule)[o].size(); ++l) {
			for (int b = 0; b < nBatches; ++b) {
				IloExpr cR05(env);
				for (int i = 0; i < nProducts; ++i) {
					cR05 += z[constMapStgMac.at({ o, l })][b][i];
				}
				mod.add(cR05 <= 1);
			}
		}
	}

//R06: job release dates
	for (int j = 0; j < nJobs; ++j) {
		IloExpr cR06(env);
		int productIdx = f[j] - 1;
		cR06 = C[j][0];
		mod.add(cR06 >= r[j] + p[productIdx][0]);
	}
 
//R07: (fixed) starting times of batches (from partial schedule)	
//R07: forall(sm in stageMachine, b in batches) gamma[sm][b] <= S[sm][b];
	for (int o = 0; o < nStages; ++o) {
		for (int l = 0; l < (*schedule)[o].size(); ++l) {
			for (int b = 0; b < nBatches; ++b) {
				IloExpr cR07(env);
				cR07 = S[constMapStgMac.at({ o, l })][b];
				mod.add(cR07 >= batStart[o][l][b]);
			}
		}
	}

//R08: (fixed) starting times of batches (from partial schedule)	
//R08: forall(sm in stageMachine, b in batches) S[sm][b] <= delta[sm][b];
	for (int o = 0; o < nStages; ++o) {
		for (int l = 0; l < (*schedule)[o].size(); ++l) {
			for (int b = 0; b < nBatches; ++b) {
				IloExpr cR08(env);
				cR08 = S[constMapStgMac.at({ o, l })][b];
				mod.add(cR08 <= batEnd[o][l][b]);
			}
		}
	}

//R09: Precedence of one job´s operations
//R09: forall(j in jobs, o in 1..O - 1) C[j][o] + p[f[j]][o] <= C[j][o + 1];	
	for (int j = 0; j < nJobs; ++j) {
		for (int o = 0; o < nStages - 1; ++o) {
			IloExpr cR09(env);
			int productIdx = f[j] - 1;
			cR09 = C[j][o + 1] - C[j][o];
			mod.add(cR09 >= p[productIdx][o + 1]);
		}
	}

//R10: time constraints
//R10: forall(j in jobs, s in stages, t in stages : s < t) (C[j][t] - p[f[j]][t]) >= (C[j][s] - p[f[j]][s]) - tc[f[j]][s][t];
	for (int j = 0; j < nJobs; ++j) {
		for (int o1 = 0; o1 < nStages; ++o1) {
			for (int o2 = o1 + 1; o2 < nStages; ++o2) {
				IloExpr cR10(env);
				int productIdx = f[j] - 1;
				cR10 = (C[j][o2] - p[productIdx][o2]) - (C[j][o1] - p[productIdx][o1]);
				mod.add(cR10 <= tc[productIdx][o1][o2]);
			}
		}
	}


//R11: Link job completion and batch start
//R11: forall(j in jobs, o in stages, b in batches, l in 1..m_o[o]) S[<o, l>][b] + p[f[j]][o] <= C[j][o] + M * (1 - X[<o, l>][j][b]);
	for (int j = 0; j < nJobs; ++j) {
		for (int o = 0; o < nStages; ++o) {
			for (int b = 0; b < nBatches; ++b) {
				for (int l = 0; l < (*schedule)[o].size(); ++l) {
					IloExpr cR11(env);
					int productIdx = f[j] - 1;
					cR11 = C[j][o] - S[constMapStgMac.at({ o, l })][b] + G * (1 - x[constMapStgMac.at({ o, l })][j][b]);
					mod.add(cR11 >= p[productIdx][o]);
				}
			}
		}
	}


//R12: Link job completion and batch start
//R12: forall(j in jobs, o in stages, b in batches, l in 1..m_o[o]) C[j][o] <= S[<o, l>][b] + p[f[j]][o] + M * (1 - X[<o, l>][j][b]);
	for (int j = 0; j < nJobs; ++j) {
		for (int o = 0; o < nStages; ++o) {
			for (int b = 0; b < nBatches; ++b) {
				for (int l = 0; l < (*schedule)[o].size(); ++l) {
					IloExpr cR12(env);
					int productIdx = f[j] - 1;
					cR12 = C[j][o] - S[constMapStgMac.at({ o, l })][b] - G * (1 - x[constMapStgMac.at({ o, l })][j][b]);
					mod.add(cR12 <= p[productIdx][o]);
				}
			}
		}
	}

//R13: Batch processing does not overlap at a machine
//R13: forall(j in jobs, o in stages, l in 1..m_o[o], b in batches, c in batches : b != c) S[<o, l>][b] + p[f[j]][o] * X[<o, l>][j][b] <= S[<o, l>][c] + M * (1 - Y[<o, l>][b][c]);	// TODO
	for (int j = 0; j < nJobs; ++j) {
		for (int o = 0; o < nStages; ++o) {
			for (int l = 0; l < (*schedule)[o].size(); ++l) {
				for (int b1 = 0; b1 < nBatches; ++b1) {
					for (int b2 = 0; b2 < nBatches; ++b2) {
						if (b1 != b2) {
							IloExpr cR13(env);
							cR13 = S[constMapStgMac.at({ o, l })][b1] + p[f[j]-1][o] * x[constMapStgMac.at({ o, l })][j][b1];
							mod.add(cR13 <= S[constMapStgMac.at({ o, l })][b2] + G * (1 - y[constMapStgMac.at({ o, l })][b1][b2]));
						}
					}
				}
			}
		}
	}

//R14: Batch processing does not overlap at a machine
//R14: forall(j in jobs, o in stages, l in 1..m_o[o], b in batches, c in batches : b != c) S[<o, l>][b] + p[f[j]][o] * X[<o, l>][j][b] <= S[<o, l>][c] + M * (1 - Y[<o, l>][b][c]);	// TODO
	for (int j = 0; j < nJobs; ++j) {
		for (int o = 0; o < nStages; ++o) {
			for (int l = 0; l < (*schedule)[o].size(); ++l) {
				for (int b1 = 0; b1 < nBatches; ++b1) {
					for (int b2 = 0; b2 < nBatches; ++b2) {
						if (b1 != b2) {
							IloExpr cR14(env);
							cR14 = S[constMapStgMac.at({ o, l })][b2] + x[constMapStgMac.at({ o, l })][j][b2] * p[f[j]-1][o] - (G * y[constMapStgMac.at({ o, l })][b1][b2]);
							mod.add(cR14 <= S[constMapStgMac.at({ o, l })][b1]);
						}
					}
				}
			}
		}
	}

//R15: No overlaps with previously scheduled batches
//R15: forall(o in stages, l in machines, b in batches, c in batches : b!=c) Y[<o,l>][b][c] * gamma[b] <= S[<o,l>][c];
	for (int o = 0; o < nStages; ++o) {
		for (int l = 0; l < (*schedule)[o].size(); ++l) {
			for (int b1 = 0; b1 < nBatches; ++b1) {
				for (int b2 = 0; b2 < nBatches; ++b2) {
					if (b1 != b2) {
						IloExpr cR15(env);
						cR15 = y[constMapStgMac.at({ o, l })][b1][b2] * (batStart[o][l][b1] + batDur[o][l][b1]);
						mod.add(cR15 <= S[constMapStgMac.at({ o, l })][b2]);
					}
				}
			}
		}
	}

//R16: No overlaps with previously scheduled batches
//R16: forall(o in stages, l in machines, b in batches, c in batches : b!=c) S[<o,l>][b] <= (1- Y[<o,l>][b][c]) * M + delta[<o,l>][c];
	for (int o = 0; o < nStages; ++o) {
		for (int l = 0; l < (*schedule)[o].size(); ++l) {
			for (int b1 = 0; b1 < nBatches; ++b1) {
				for (int b2 = 0; b2 < nBatches; ++b2) {
					if (b1 != b2) {
						for (int j = 0; j < nJobs; ++j) {
							int productIdx = f[j] - 1;
							IloExpr cR16(env);
							cR16 = (1 - y[constMapStgMac.at({ o, l })][b1][b2]) * G + batEnd[o][l][b2] - batDur[o][l][b2];
							mod.add(cR16 >= S[constMapStgMac.at({ o, l })][b1] + p[productIdx][o] * x[constMapStgMac.at({ o, l })][j][b1]);
						}
					}
				}
			}
		}
	}

//R17: Link job completion with previously scheduled batches
//R17: forall(o in stages, l in machines, b in batches, j in jobs) X[sm][j][b] * gamma[sm][b] + p[f[j]][o] <= C[j][o];	
	for (int o = 0; o < nStages; ++o) {
		for (int l = 0; l < (*schedule)[o].size(); ++l) {
			for (int b = 0; b < nBatches; ++b) {
				for (int j = 0; j < nJobs; ++j) {
					IloExpr cR17(env);
					cR17 = x[constMapStgMac.at({ o, l })][j][b] * batStart[o][l][b] + p[f[j] - 1][o];
					mod.add(cR17 <= C[j][o]);
				}
			}
		}
	}

//R18: Link job completion with previously scheduled batches
//R18: forall(o in stages, l in machines, b in batches, j in jobs) X[sm][j][b] * gamma[sm][b] + p[f[j]][o] <= C[j][o];	
	for (int o = 0; o < nStages; ++o) {
		for (int l = 0; l < (*schedule)[o].size(); ++l) {
			for (int b = 0; b < nBatches; ++b) {
				for (int j = 0; j < nJobs; ++j) {
					IloExpr cR18(env);
					cR18 = batEnd[o][l][b] + G * (1 - x[constMapStgMac.at({ o, l })][j][b]);
					mod.add(C[j][o] <= cR18);
				}
			}
		}
	}

//R19: Compatibility: batch compatiblity of previously scheduled batches
//R19: forall(sm in stageMachine, b in batches, j in jobs) X[sm][j][b] <= e[j][f[j]] + (1 - eB[sm][b][f[j]]);
	for (int j = 0; j < nJobs; ++j) {												// forall(j in jobs,
		for (int i = 0; i < nProducts; ++i) {										//   i in products,
			for (int o = 0; o < nStages; ++o) {										//   
				for (int l = 0; l < (*schedule)[o].size(); ++l) {					//   sm in stageMachine
					for (int b = 0; b < nBatches; ++b) {							//   b in batches
						int e = i == f[j] - 1;
						int eB = i == batFam[o][l][b] - 1;
						IloExpr cR19(env);
						cR19 = x[constMapStgMac.at({ o, l })][j][b];
						mod.add(cR19 <= e + (1 - eB));
					}
				}
			}
		}
	}

//R20: Tardiness
//R20: forall(j in jobs) C[j][O] - d[j] <= T[j];
	for (int j = 0; j < nJobs; ++j) {
		IloExpr cR20(env);
		cR20 = C[j][nStages - 1] - T[j];
		mod.add(cR20 <= d[j]);
	}


// SET TIMELIMIT AND SOLVE

	IloCplex cplex(mod);
	cplex.setOut(env.getNullStream());									// GO-LIVE: uncomment
	//cplex.exportModel("gapAware_milp.lp");							// GO-LIVE: comment
	cplex.setParam(IloCplex::EpGap, params->presetMipGap);	
	cplex.setParam(IloCplex::TiLim, cplexTilim);

	try {
		if (!cplex.solve()) {
			std::stringstream errorMsg;
			IloAlgorithm::Status status = cplex.getStatus();
			errorMsg << "CPLEX failed to solve the MILP. Status: " << cplex.getStatus() << ". ";
			if (status == IloAlgorithm::Infeasible) {
				errorMsg << "The model is infeasible. ";
				// TODO Use conflict refiner to identify conflicting constraints
			}
			else if (status == IloAlgorithm::Unbounded) {
				errorMsg << "The model is unbounded. ";
			}
			else if (status == IloAlgorithm::InfeasibleOrUnbounded) {
				errorMsg << "The model is infeasible or unbounded. ";
			}
			else if (status == IloAlgorithm::Error) {
				errorMsg << "An error occurred during optimization. Error code: " << cplex.getCplexStatus() << ". ";
			}
			TCB::logger.Log(Error, errorMsg.str());
			cout << errorMsg.str() << endl;
		}
	}
	catch (IloException& e) {
		TCB::logger.Log(Error, e.getMessage());
		env.end();
	}

	params->relMipGap = cplex.getMIPRelativeGap();	// output rel MIP gap of last iteration (only really interesting for exact solutions with nDash == n

	cout << ".";


// DEBUGGING
/*for (int j = 0; j < nJobs; ++j) {
	for (int stg = 0; stg < nStages; ++stg) {
		IloNum cJ = cplex.getValue(C[j][stg]);
		cout << "J" << to_string(consideredJobs[j]->getId()) << "." << to_string(stg + 1) << ": " << to_string(cJ) << ", ";
	}
	cout << endl;
}

cout << "BatStart:" << endl;
for (int o = 0; o < nStages; ++o) {
	cout << "Stage " << to_string(o) << ": " << endl;
	for (int l = 0; l < (*schedule)[o].size(); ++l) {
		cout << "M" << to_string(l + 1) << ": ";
		for (int b = 0; b < nBatches; ++b) {
			cout << "B" << to_string(b) << ": " << to_string(batStart[o][l][b]) << ", ";
		}
		cout << endl;
	}
	cout << endl;
}

cout << "BatDur:" << endl;
for (int o = 0; o < nStages; ++o) {
	cout << "Stage " << to_string(o) << ": " << endl;
	for (int l = 0; l < (*schedule)[o].size(); ++l) {
		cout << "M" << to_string(l + 1) << ": ";
		for (int b = 0; b < nBatches; ++b) {
			cout << "B" << to_string(b) << ": " << to_string(batDur[o][l][b]) << ", ";
		}
		cout << endl;
	}
	cout << endl;
}

cout << "BatEnd:" << endl;
for (int o = 0; o < nStages; ++o) {
	cout << "Stage " << to_string(o) << ": " << endl;
	for (int l = 0; l < (*schedule)[o].size(); ++l) {
		cout << "M" << to_string(l + 1) << ": ";
		for (int b = 0; b < nBatches; ++b) {
			cout << "B" << to_string(b) << ": " << to_string(batEnd[o][l][b]) << ", ";
		}
		cout << endl;
	}
	cout << endl;
}*/


// CONSTRUCT SCHEDULE FROM SOLUTION	
	for (int stg = 0; stg < nStages; ++stg) {
		for (int mac = 0; mac < (*schedule)[stg].size(); ++mac) {
			// derive sequence from y-variables
			vector<int> batchSequence = vector<int>(nBatches);
			vector<int> precedence = vector<int>(nBatches, 0);
			iota(batchSequence.begin(), batchSequence.end(), 0);

			for (int b1 = 0; b1 < nBatches; ++b1) {
				for (int b2 = 0; b2 < nBatches; ++b2) {
					if (b1 != b2) {
						IloBool yB1B2 = cplex.getValue(y[constMapStgMac.at({ stg, mac })][b1][b2]);
						if (yB1B2 >= 1 - TCB::precision) {
							// b1 is processed before b2
							++precedence[b2];
						}
					}
				}
			}
			sort(batchSequence.begin(), batchSequence.end(), [&](int a, int b) {
				return precedence[a] < precedence[b];
				});

			
			for (int bat = 0; bat < nBatches; ++bat) {
				int b = batchSequence[bat];

				// IDENTITFY OPERATIONS TO BE ASSIGNED TO THIS BATCH
				vector<int> opsToAssign = vector<int>();
				for (int j = 0; j < nJobs; ++j) {
					IloBool xJB = cplex.getValue(x[constMapStgMac.at({ stg, mac })][j][b]);
					if (xJB >= 1 - TCB::precision) {
						opsToAssign.push_back(j);
					}
				}
				// GET BATCH START TIME
				IloNum sB = cplex.getValue(S[constMapStgMac.at({ stg, mac })][b]);

				if (!opsToAssign.empty()) {
					if (B[stg][mac][b] < (*schedule)[stg].getCap()) {
						// FIND PREVIOUSLY SCHEDULED BATCH
						for (int op = 0; op < opsToAssign.size(); ++op) {
							if (!existingBatches[stg][mac][b]->addOp(consideredJobs[opsToAssign[op]]->getOpPtr(stg))) {
								throw ExcSched("Could not add op to batch (invalid solution?), SolverMILP::initModel()");
							}
						}

						//int start = batStart[stg][mac][b];
						//for (int schedB = 0; schedB < (*schedule)[mac].size(); ++schedB) {
						//	if ((*schedule)[stg][mac].getBatch(schedB)->getStart() == start) {
						//		// ASSIGN OPERATIONS
						//		for (int op = 0; op < opsToAssign.size(); ++op) {
						//			if (!(*schedule)[stg][mac].getBatch(schedB)->addOp(consideredJobs[opsToAssign[op]]->getOpPtr(stg))) {
						//				throw ExcSched("Could not add op to batch (invalid solution?), SolverMILP::initModel()");
						//			}
						//		}
						//	}
						//}	
					}
					else {
						// CREATE NEW BATCH
						pBat newBatch = make_unique<Batch>((*schedule)[stg].size());
						// ASSIGN OPERATIONS
						for (int op = 0; op < opsToAssign.size(); ++op) {
							if (!newBatch->addOp(consideredJobs[opsToAssign[op]]->getOpPtr(stg))) {
								throw ExcSched("Could not add op to batch (invalid solution?), SolverMILP::initModel()");
							}
						}
						// ASSIGN BATCH
						if (!(*schedule)[stg][mac].addBatch(move(newBatch), (double) sB)) {
							throw ExcSched("BATCH ASSIGNMENT FAILED FROM MILP SOLUTION");
						}
					}
				}
			}				
		}
	}

	//schedule->saveJsonFactory("TEST_NEWMILP");	// DEBUGGING

	while (consideredJobs.size() > 0) {
		schedule->markAsScheduled(move(consideredJobs[0]));
		consideredJobs.erase(consideredJobs.begin());
	}

	env.end();
	return schedule->getTWT();	// bestObjectiveValue; // integrated objective function

}
double Solver_MILP::solveDecompGapAwareMILP(Schedule* schedule, int nDash, int cplexTilim) {
	double twt = FLT_MAX;
	int iteration = 1;	// DEBUGGING
	while (schedule->getN() > 0) {
		// SOLVE with nDash Jobs
		twt = solveGapAwareMILP(schedule, nDash, cplexTilim);
		iteration++;	// DEBUGGING
	}
	return twt;
}
double Solver_MILP::solveDecompGapAwareDynamicSortingGridMILP(Schedule* schedule, int nDash, int cplexTilim, prioRule<pJob> initRule, prioRuleKappa<pJob> dynRule) {
	vector<double> kappas = getDoubleGrid(schedParams->kappaLow, schedParams->kappaHigh, schedParams->kappaStep);
	// OBTAIN INITIAL WAITING TIMES FROM A SCHEDULE WITH LIST SCHEDULING AND SORTING BY initRule (DEFAULT: EDD)
	unique_ptr<Schedule> initSched = schedule->clone();
	initSched->lSchedJobsWithSorting(initRule);
	initSched->updateWaitingTimes();

	// DEBUGGING
	//initSched->saveJsonFactory("INITIAL_SCHEDULE_EDD");
	//schedule->sortUnscheduled(sortJobsDebugging);	// DEBUGGING 23.06.2026


	schedule->mimicWaitingTimes(initSched.get());

	while (schedule->getN() > 0) {
		// DYNAMIC SORTING AND UPDATE WAITING TIMES 
		
		// DEBUGGING
		schedule->saveJsonFactory("DEBUGGING");

		unique_ptr<Schedule> tempSched = schedule->clone();
		double bestKappa = tempSched->lSchedJobsWithSorting(dynRule, kappas);	// GRID SEARCH
		tempSched->updateWaitingTimes();
		schedule->mimicWaitingTimes(tempSched.get());
		schedule->sortUnscheduled(dynRule, bestKappa);

		// ASSIGN THE NEXT NDASH JOBS

		solveGapAwareMILP(schedule, nDash, cplexTilim);

		// DEBUGGING
		//schedule->saveJsonFactory("PARTIAL_SCHEDULE_AFTER_ITERATION");
		//int breakpoint = 666;

	}


	return schedule->getTWT();
}

DECOMPMILP_params Solver_MILP::getDefaultParams() {
	DECOMPMILP_params decompParams = DECOMPMILP_params();
	decompParams.nDash = 10;
	decompParams.cplexTilim = 60;
	decompParams.method = DECOMP_SOLVER_MILP;
	decompParams.initPrioRule = PRIORULE_EDD;
	decompParams.presetMipGap = 0.0001;
	return decompParams;
}





