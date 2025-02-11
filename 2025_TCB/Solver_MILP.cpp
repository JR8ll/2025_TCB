#include<ilcplex/ilocplex.h>
#include<ilconcert/ilotupleset.h>
#include<map>
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

Solver_MILP::Solver_MILP() {}
Solver_MILP::~Solver_MILP() {}

double Solver_MILP::solveDecompositionMILP(Schedule* schedule, int nDash, int cplexTilim) {
	if (schedule->getProblem() == nullptr) throw ExcSched("solveDecompositionMILP missing problem reference");

	IloEnv env;

	// PARAMETERS
	IloNum omega = 9;	// TODO set omega

	IloInt G = 99999;	// TODO compute
	IloInt m = schedule->size();

	IloIntTupleSet stageMachine(env, 2);
	for (int o = 0; o < o; ++o) {
		for (int i = 0; i < (*schedule)[o].size(); ++i) {
			stageMachine.add(IloIntArray(env, 2, o, i));
		}
	}
	map<pair<int, int>, int> mapStgMac = map<pair<int, int>, int>();	// for access to tuple "stageMachine"
	int index = 0;
	for (IloIntTupleSetIterator it(env, stageMachine); it.ok(); ++it) {
		mapStgMac.insert(make_pair(make_pair((*it)[0], (*it)[1]), index));
		index++;
	}

	IloIntArray B(env);	// capacity[stage]
	for (int o = 0; o < m; ++o) {
		int tempCap = schedule->getCapAtStageIdx(o);	// assumption: parallel identical machines
		if (tempCap > 1) {
			B.add(tempCap);
		}
	}

	int nProducts = schedule->getProblem()->getF();
	int nStages = schedule->size();
	int nJobs = schedule->getN();

	IloIntArray3 B_iob(env);	// number of jobs assigned to preformed Batch[Products][Stages][Batches]
	IloNumArray3 S_iob(env);	// start times of preformed Batch[Products][Stages][Batches]	
	vector<vector<set<int> > > B_io_preformed = vector<vector<set<int> > >(nProducts);

	for (int i = 0; i < nProducts; ++i) {
		B_io_preformed[i] = vector<set<int> >(nStages);
		B_iob.add(IloIntArray2(env));
		S_iob.add(IloNumArray2(env));
		for (int o = 0; o < nStages; ++o) {
			B_io_preformed[i][o] = set<int>();
			B_iob[i].add(IloIntArray(env, nJobs));
			S_iob[i].add(IloNumArray(env, nJobs));
		}
	}
	vector<map<pair<int, int>, pair<int, int> > > batchIdxMac = vector<map<pair<int, int>, pair<int, int> > >(nStages);
	for (int o = 0; o < nStages; ++o) {
		batchIdxMac[o] = map<pair<int, int>, pair<int, int> >();
		int batchIdx = 0;
		for (int l = 0; l < (*schedule)[o].size(); ++l) {
			for (int b = 0; b < (*schedule)[o][l].size(); ++b) {
				int productIdx = (*schedule)[o][l][b].getF() - 1;
				if (productIdx >= 0) {
					B_io_preformed[productIdx][o].insert(batchIdx);
					batchIdxMac[o].insert(make_pair(make_pair(productIdx, batchIdx), make_pair(l, b)));
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
		d.add(schedule->getJob(j).getD());
		r.add(schedule->getJob(j).getR());
		w.add(schedule->getJob(j).getW());
		f.add(schedule->getJob(j).getF());
	}

	IloNumArray2 p(env);	// p[Products][Stages] processing times
	for (int f = 0; f < nProducts; ++f) {
		p.add(IloNumArray(env, nStages));
		for (int o = 0; o < nStages; ++o) {
			p[f][o] = schedule->getProblem()->getProduct(f)->getP(o); 
		}
	}

	IloNumArray3 tc(env);	// tc[Products][Stages][Stages] time constraints
	for (int f = 0; f < nProducts; ++f) {
		tc.add(IloNumArray2(env));
		for (int o1 = 0; o1 < nStages; ++o1) {
			tc[f].add(IloNumArray(env, nStages));
			const vector<pair<int, double>>& tcMaxFwd = schedule->getProblem()->getProduct(f)->getTcMaxFwd(o1);
			for (size_t o2 = 0; o2 < nStages; ++o2) {
				double tcDuration = 999999;	// TODO define INF
				for (auto it = tcMaxFwd.begin(); it != tcMaxFwd.end(); ++it) {
					if (it->first == o2) {
						tcDuration = it->second;
					}
				}
				tc[f][o1][o2] = (IloNum)tcDuration;
			}
		}
	}

	IloModel mod(env);

	// DECISION VARIABLES
	IloBoolVarArray3 x(env);	// x[Jobs][Batches][<StageMachine>]
	for (int j = 0; j < nJobs; ++j) {
		x.add(IloBoolVarArray2(env));
		for (int b = 0; b < nJobs; ++b) {
			x[j].add(IloBoolVarArray(env, stageMachine.getCardinality()));
			for (int k = 0; k < stageMachine.getCardinality(); ++k) {
				x[j][b][k] = IloBoolVar(env);
				string varName = "x_" + to_string(j) + "," + to_string(b) + "," + to_string(k);
				x[j][b][k].setName(varName.c_str());
			}
		}
	}

	IloBoolVarArray3 y(env);	// y[Products][Batches][<StageMachine>]
	for (int f = 0; f < nProducts; ++f) {
		y.add(IloBoolVarArray2(env));
		for (int b = 0; b < nJobs; ++b) {
			y[f].add(IloBoolVarArray(env, stageMachine.getCardinality()));
			for (int k = 0; k < stageMachine.getCardinality(); ++k) {
				y[f][b][k] = IloBoolVar(env);
				string varName = "y_" + to_string(f) + "," + to_string(b) + "," + to_string(k);
				y[f][b][k].setName(varName.c_str());
			}
		}
	}

	IloBoolVarArray3 z(env);	// z[Jobs][Stages][Batches]	
	for (int j = 0; j < nJobs; ++j) {
		z.add(IloBoolVarArray2(env));
		for (int o = 0; o < nStages; ++o) {
			z[j].add(IloBoolVarArray(env, nJobs));		
			for (int b = 0; b < nJobs; ++b) {
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
	for (int b = 0; b < nJobs; ++b) {
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
	obj = omega * objTwt + objJobStart;

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
							cR03_SumX += x[j][b][mapStgMac[{oBatch, l}]];
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
		for (int b = 0; b < nJobs; ++b) {												// b in Batches
			for (int l = 0; l < (*schedule)[oBatch].size(); ++l) {						// l in 1..m[o]
				IloExpr cR04_SumX(env);
				for (int j = 0; j < nJobs; ++j) {										// sum(j in Jobs)
					cR04_SumX += x[j][b][mapStgMac[{oBatch, l}]];
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
			cR05 = cR05_SumZ + B_iob[it->first.first][o][it->first.second];
			string conLabel = "R05 (o" + to_string(oBatch) + "b" + to_string(it->first.second) + "i" + to_string(it->first.first) + ")";
			mod.add(cR05 <= B[o]).setName(conLabel.c_str());
		}
	}

	// R06: unique family of a batch													// forall..
	for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {					// o in Stages_b
		int oBatch = schedule->getBatchingStages()[o] - 1;
		for (int b = 0; b < nJobs; ++b) {												// b in Batches
			for (int l = 0; l < (*schedule)[oBatch].size(); ++l) {						// l in 1..m[o]
				IloExpr cR06_SumY(env);
				for (int i = 0; i < nProducts; ++i) {									// sum(i in Products)
					cR06_SumY += y[i][b][mapStgMac[{oBatch, l}]];
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
		for (int b = 0; b < nJobs; ++b) {												// b in Batches
			for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {			// o in Stages_b
				int oBatch = schedule->getBatchingStages()[o] - 1;
				for (int l = 0; l < (*schedule)[oBatch].size(); ++l) {					// l in 1..m_o[o]
					for (int i = 0; i < nProducts; ++i) {
						if (f[j] == (i + 1)) {											// i in Products : i == f[j]
							IloExpr cR07(env);
							cR07 = x[j][b][mapStgMac[{oBatch, l}]]; // <= y[i][b][oBatch, l];
							string conLabel = "R07 (j" + to_string(j) + "b" + to_string(b) + "o" + to_string(oBatch) + "l" + to_string(l) + "i" + to_string(i) + ")";
							mod.add(cR07 <= y[i][b][mapStgMac[{oBatch, l}]]).setName(conLabel.c_str());

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
				cR08_SumU += u[j][mapStgMac[{oSingle, l}]];
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
		for (int o = 0; o < schedule->getDiscreteStages().size(); ++o) {						// o in Stages_1
			int oSingle = schedule->getDiscreteStages()[o] - 1;
			for (int l = 0; l < (*schedule)[oSingle].size(); ++l) {					// l in 1..m_o[o]
				IloExpr cR10(env);
				cR10 = s[j][oSingle] + (rm[oSingle][l] * (1 - u[j][mapStgMac[{oSingle, l}]])); // >= rm[oSingle][l];
				string conLabel = "R10 (j" + to_string(j) + "o" + to_string(oSingle) + "l" + to_string(l) + ")";
				mod.add(cR10 >= rm[oSingle][l]).setName(conLabel.c_str());
			}
		}
	}

	// R11: release dates of batching machines											// forall..
	for (int j = 0; j < nJobs; ++j) {													// j in Jobs
		for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {				// o in Stages_b
			int oBatch = schedule->getDiscreteStages()[o] - 1;
			for (int b = 0; b < nJobs; ++b) {											// b in Batches
				for (int l = 0; l < (*schedule)[oBatch].size(); ++l) {					// l in 1..m_o[o]
					IloExpr cR11(env);
					cR11 = S[b][mapStgMac[{oBatch, l}]]; // >= rm[oBatch][l];
					string conLabel = "R11 (j" + to_string(j) + "o" + to_string(oBatch) + "b" + to_string(b) + "l" + to_string(l) + ")";
					mod.add(cR11 >= rm[oBatch][l]).setName(conLabel.c_str());
				}
			}
		}
	}

	// R12: processing times of batches													// forall..
	for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {					// o in Stages_b
		int oBatch = schedule->getBatchingStages()[o] - 1;
		for (int b = 0; b < (nJobs - 1); ++b) {											// b in 1..(n-1)
			for (int l = 0; l < (*schedule)[oBatch].size(); ++l) {						// l in 1..m_o[o]
				IloExpr cR12_SumY(env);
				for (int i = 0; i < nProducts; ++i) {									// sum(i in Products)
					cR12_SumY += p[i][oBatch] * y[i][b][mapStgMac[{oBatch, l}]];
				}
				IloExpr cR12(env);
				cR12 = S[b][mapStgMac[{oBatch, l}]] + cR12_SumY; // <= S[b + 1][oBatch, l];
				string conLabel = "R12 (o" + to_string(oBatch) + "b" + to_string(b) + "l" + to_string(l) + ")";
				mod.add(cR12 <= S[b + 1][mapStgMac[{oBatch, l}]]).setName(conLabel.c_str());
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
					cR14 = s[j][o] + p[i][o]; // -s[j][o + 1]; // <= p[i][o];
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
						cR15 = u[j][mapStgMac[{oSingle, l}]] + u[k][mapStgMac[{oSingle, l}]] - e[j][k][oSingle] - e[k][j][oSingle]; // <= 1;
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

	// R17: matching start times of jobs and batches									// forall..
	for (int j = 0; j < nJobs; ++j) {												// j in Jobs
		for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {						// o in Stages_b
			int oBatch = schedule->getBatchingStages()[o] - 1;
			for (int b = 0; b < nJobs; ++b) {									// b in Batches
				for (int l = 0; l < (*schedule)[oBatch].size(); ++l) {					// l in 1..m_o[o]
					IloExpr cR17(env);
					cR17 = s[j][oBatch] + (G * (1 - x[j][b][mapStgMac[{oBatch, l}]])); // >= S[b][oBatch, l];
					string conLabel = "R17 (j" + to_string(j) + "o" + to_string(oBatch) + "b" + to_string(b) + "l" + to_string(l) + ")";
					mod.add(cR17 >= S[b][mapStgMac[{oBatch, l}]]).setName(conLabel.c_str());
				}
			}
		}
	}

	// R18: matching start times of jobs and batches									// forall..
	for (int j = 0; j < nJobs; ++j) {												// j in Jobs
		for (int o = 0; o < schedule->getBatchingStages().size(); ++o) {						// o in Stages_b
			int oBatch = schedule->getBatchingStages()[o] - 1;
			for (int b = 0; b < nJobs; ++b) {									// b in Batches
				for (int l = 0; l < (*schedule)[oBatch].size(); ++l) {					// l in 1..m_o[o]
					IloExpr cR18(env);
					cR18 = s[j][oBatch] + (G * (x[j][b][mapStgMac[{oBatch, l}]] - 1)); // <= S[b][oBatch, l];
					string conLabel = "R18 (j" + to_string(j) + "o" + to_string(oBatch) + "b" + to_string(b) + "l" + to_string(l) + ")";
					mod.add(cR18 <= S[b][mapStgMac[{oBatch, l}]]).setName(conLabel.c_str());
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
	cplex.setParam(IloCplex::TiLim, cplexTilim);	// TODO extract timelimit parameter

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

	// construct schedule from solution
	for (int stg = 0; stg < nStages; ++stg) {
		if (schedule->getCapAtStageIdx(stg) <= 1) {
			// no batching at this stages
			for (int j = 0; j < nJobs; ++j) {
				pBat newSingleBatch = make_unique<Batch>(schedule->getCapAtStageIdx(stg));	// assumption: parallel identical machines
				newSingleBatch->addOp(&schedule->getJob(j)[stg]);
				int macIdx;
				for (int mac = 0; mac < (*schedule)[stg].size(); ++mac) {
					IloNum val = cplex.getValue(u[j][mapStgMac[{stg, mac}]]);		// u[j][stgMac]	: machine assignment
					if (val >= 1 - TCB::precision) {
						macIdx = mac;
					}
				}
				IloNum start = cplex.getValue(s[j][stg]);
				if (!(*schedule)[stg][macIdx].addBatch(move(newSingleBatch), start)) { 
					std::cout << "Problem in SolverMILP::solve(): infeasible op assignment " << (j + 1) << " at stage " << (stg + 1) << ": " << start << endl;
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
							if (!(*schedule)[stg][it->second.first][it->second.second].addOp(&schedule->getJob(j)[stg])) {
								stringstream errorMsg;
								errorMsg << "Problem in SolverMILP::solve(): infeasible op assignment " << schedule->getJob(j).getId() << " at stage " << (stg + 1);
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
			int batchingStageIndex = 0;
			const vector<int>& stagesB = schedule->getBatchingStages();
			for (int i = 0; i < stagesB.size(); ++i) {
				if (stagesB[i] == (stg + 1)) {
					batchingStageIndex = i;
					break;
				}
			}
			int capacity = schedule->getCapAtStageIdx(batchingStageIndex);
			vector<vector<pBat>> newBatches = vector<vector<pBat>>((*schedule)[stg].size());
			for (int o = 0; o < (*schedule)[stg].size(); ++o) {
				newBatches[o].resize(nJobs);
			}
			for (int mac = 0; mac < (*schedule)[stg].size(); ++mac) {
				for (int b = 0; b < nJobs; ++b) {
					newBatches[mac][b] = make_unique<Batch>(Batch(capacity));
				}
			}

			for (int b = 0; b < nJobs; ++b) {
				for (int j = 0; j < nJobs; ++j) {
					for (int mac = 0; mac < (*schedule)[stg].size(); ++mac) {
						IloNum xJBOL = cplex.getValue(x[j][b][mapStgMac[{stg, mac}]]);
						if (xJBOL >= 1 - TCB::precision) {
							if (!newBatches[mac][b]->addOp(&schedule->getJob(j)[stg])) {
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
						IloNum xJBOL = cplex.getValue(x[j][b][mapStgMac[{stg, mac}]]);
						if (xJBOL >= 1 - TCB::precision) {
							IloNum sBOL = cplex.getValue(S[b][mapStgMac[{stg, mac}]]);
							try {
								(*schedule)[stg][mac].addBatch(move(newBatches[mac][b]), sBOL); //  preformed.at(stg).at(mac).addBatch(newBatches[mac][b], sBOL);
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
	double bestObjectiveValue = cplex.getBestObjValue();
	env.end();
	return bestObjectiveValue;
}




