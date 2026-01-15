#include <algorithm>
#include <iostream>
#include <queue>
#include <map>
#include <numeric>
#include <Windows.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "Schedule.h"
#include "Functions.h"
#include "Problem.h"
#include "Solver_MILP.h"	// DECOMPMILP_params

namespace pt = boost::property_tree;
using namespace std;

Schedule::Schedule() {
	workcenters = vector<pWc>();
	unscheduledJobs = vector<pJob>();
	scheduledJobs = vector<pJob>();
	problem = nullptr;
}

ostream& operator<<(ostream& os, const Schedule& sched) {
	os << "---------SCHEDULE ----------" << endl;
	for (size_t wc = 0; wc < sched.size(); ++wc) {
		os << sched[wc] << endl;
	}
	return os;
}

Workcenter& Schedule::operator[](size_t idx) { return *workcenters[idx]; }
Workcenter& Schedule::operator[](size_t idx) const { return *workcenters[idx]; }

Job& Schedule::getJob(size_t idx) { return *unscheduledJobs[idx]; };
Job& Schedule::getJob(size_t idx) const { return *unscheduledJobs[idx]; }

pJob Schedule::get_pJob(size_t idx) {
	if (idx >= unscheduledJobs.size()) throw out_of_range("Schedule::get_pJob() out of range");
	pJob returnJob = move(unscheduledJobs[idx]);
	unscheduledJobs.erase(unscheduledJobs.begin() + idx);
	return returnJob;
}

std::unique_ptr<Schedule> Schedule::clone() const {
	auto newSchedule = make_unique<Schedule>();
	newSchedule->setProblemRef(problem);
	
	for (const auto& wc : workcenters) {
		newSchedule->addWorkcenter(wc->clone(newSchedule.get()));
	}
	for (const auto& job : unscheduledJobs) {
		newSchedule->addJob(move(job->clone()));
	}
	for (const auto& job : scheduledJobs) {
		newSchedule->markAsScheduled(move(job->clone()));
	}
	newSchedule->_reconstruct(this);
	return newSchedule;
}

void Schedule::_reconstruct(const Schedule* orig) {
	for (size_t wc = 0; wc < (*orig).size(); ++wc) {
		for (size_t m = 0; m < (*orig)[wc].size(); ++m) {
			for (size_t b = 0; b < (*orig)[wc][m].size(); ++b) {
				(*this)[wc][m].addBatch(move((*orig)[wc][m][b].clone()), (*orig)[wc][m][b].getStart());
				for (size_t j = 0; j < (*orig)[wc][m][b].size(); ++j) {
					Operation* remoteOp = &(*orig)[wc][m][b][j];
					Operation* localOp = findInScheduledJobs(remoteOp);
					if (localOp == nullptr) {
						localOp = findInUnscheduledJobs(remoteOp);
					}
					if (localOp == nullptr) {
						TCB::logger.Log(Error, "Exception thrown in Schedule::_reconstruct(...)");
						throw ExcSched("Schedule::_reconstruct() operation not found");
					}
					(*this)[wc][m][b].addOp(localOp);
				}
			}
		}
	}
	for (size_t wc = 0; wc < size(); ++wc) {
		for (size_t m = 0; m < workcenters[wc]->size(); ++m) {
			for (size_t b = 0; b < (*workcenters[wc])[m].size(); ++b) {
				(*this)[wc][m][b].setStart((*this)[wc][m][b].getStart());
			}
		}
	}
}

int Schedule::size() const { return workcenters.size();  }
size_t Schedule::getN() const { return unscheduledJobs.size(); }

bool Schedule::contains(Operation* op) const {
	for (size_t wc = 0; wc < size(); ++wc) {
		for (size_t m = 0; m < (*workcenters[wc]).size(); ++m) {
			for (size_t b = 0; b < (*workcenters[wc])[m].size(); ++b) {
				for (size_t o = 0; o < (*workcenters[wc])[m][b].size(); ++o) {
					if (op->getId() == (*workcenters[wc])[m][b][o].getId() && op->getStg() == (*workcenters[wc])[m][b][o].getStg()) {
						return true;
					}
				}
			}
		}
	}
	return false;
}

int Schedule::getCapAtStageIdx(size_t stgIdx) const {
	if (stgIdx >= size()) throw out_of_range("Schedule::getCapAtStageIdx() out of range");
	return (*workcenters[stgIdx])[0].getCap();
}

const vector<int> Schedule::getBatchingStages() const {
	vector<int> batchingStages = vector<int>();
	for (size_t o = 0; o < size(); ++o) {
		if ((*workcenters[o])[0].getCap() > 1) {	// assumption: parallel identical machines
			batchingStages.push_back(o + 1);
		}
	}
	return batchingStages;
}
const std::vector<int> Schedule::getDiscreteStages() const {
	vector<int> discreteStages = vector<int>();
	for (size_t o = 0; o < size(); ++o) {
		if ((*workcenters[o])[0].getCap() <= 1) {	// assumption: parallel identical machines
			discreteStages.push_back(o + 1);
		}
	}
	return discreteStages;
}

const std::vector<pWc>& Schedule::getWorkcenters() const {
	return workcenters;
}
void Schedule::addWorkcenter(pWc wc) {
	workcenters.push_back(move(wc));
}
void Schedule::addJob(pJob job) {
	unscheduledJobs.push_back(move(job));
}

void Schedule::schedOp(Operation* op, double pWait) {
	int wcIdx = op->getWorkcenterId() - 1;
 	workcenters[wcIdx]->schedOp(op, pWait);
}
void Schedule::schedOp(Operation* op, double pWait, double inflation, bool batchinStageInflationOnly, bool opsWithoutTcInflationOnly) {
	int wcIdx = op->getWorkcenterId() - 1;
	workcenters[wcIdx]->schedOp(op, pWait, inflation, batchinStageInflationOnly, opsWithoutTcInflationOnly);
}
void Schedule::schedOpDelayed(Operation* op, double startingAt) {
	int wcIdx = op->getWorkcenterId() - 1;
	workcenters[wcIdx]->schedOpDelayed(op, startingAt);
}

Problem* Schedule::getProblem() const {
	return problem;
}

void Schedule::setProblemRef(Problem* prob) {
	problem = prob;
}

void Schedule::reset() {
	for (size_t wc = 0; wc < size(); ++wc) {
		for (size_t m = 0; m < workcenters[wc]->size(); ++m) {
			(*workcenters[wc])[m].removeAllBatches();
		}
	}
	while (!scheduledJobs.empty()) {
		shiftJobFromVecToVec(scheduledJobs, unscheduledJobs, 0);
	}

	for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
		for (size_t o = 0; o < (*unscheduledJobs[j]).size(); ++o) {
			(*unscheduledJobs[j])[o].setWait(0);
		}
	}
}
void Schedule::clearJobs() {
	for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
		for (size_t o = 0; o < (*unscheduledJobs[j]).size(); ++o) {	
			(*unscheduledJobs[j])[o].setPred(nullptr);
			(*unscheduledJobs[j])[o].setSucc(nullptr);
		}
	}

	for (size_t j = 0; j < scheduledJobs.size(); ++j) {
		for (size_t o = 0; o < (*scheduledJobs[j]).size(); ++o) {
			(*scheduledJobs[j])[o].setPred(nullptr);
			(*scheduledJobs[j])[o].setSucc(nullptr);
		}
	}
	unscheduledJobs.clear();
	scheduledJobs.clear();
}

void Schedule::sortUnscheduled(prioRule<pJob> rule) {
	rule(unscheduledJobs);
}

void Schedule::sortUnscheduled(prioRuleKappa<pJob> rule, double kappa) {
	double t = getMinMSP(0);
	rule(unscheduledJobs, t, kappa);
}

void Schedule::sortUnscheduled(prioRuleKeySet<pJob> rule, std::vector<double>& chr) {
	rule(unscheduledJobs, chr);
}

void Schedule::sortScheduled(prioRule<pJob> rule) {
	updateWaitingTimes();
	rule(scheduledJobs);
}

void Schedule::updateWaitingTimes() {
	for (size_t wc = 0; wc < size(); ++wc) {
		workcenters[wc]->updateWaitingTimes();
	}
}

void Schedule::mimicWaitingTimes(const Schedule* wtSched) {
	map<pair<int, int>, double> mapWait = map<pair<int, int>, double>();	// <job id, stage id>, waiting time 
	for (size_t wc = 0; wc < wtSched->size(); ++wc) {
		for (size_t m = 0; m < (*wtSched)[wc].size(); ++m) {
			for (size_t b = 0; b < (*wtSched)[wc][m].size(); ++b) {
				for (size_t op = 0; op < (*wtSched)[wc][m][b].size(); ++op) {
					mapWait.insert(make_pair(make_pair((*wtSched)[wc][m][b][op].getId(), (*wtSched)[wc][m][b][op].getStg()), (*wtSched)[wc][m][b][op].getWait()));
				}
			}
		}
	}

	for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
		for (size_t o = 0; o < unscheduledJobs[j]->size(); ++o) {
			(*unscheduledJobs[j])[o].setWait(mapWait[make_pair((*unscheduledJobs[j])[o].getId(), (*unscheduledJobs[j])[o].getStg())]);
		}
	}

	for (size_t j = 0; j < scheduledJobs.size(); ++j) {
		for (size_t o = 0; o < scheduledJobs[j]->size(); ++o) {
			(*scheduledJobs[j])[o].setWait(mapWait[make_pair((*scheduledJobs[j])[o].getId(), (*scheduledJobs[j])[o].getStg())]);
		}
	}
}

void Schedule::markAsScheduled(size_t jobIdx) {
	if (jobIdx >= unscheduledJobs.size()) throw out_of_range("Schedule::markAsScheduled() out of range");
	shiftJobFromVecToVec(unscheduledJobs, scheduledJobs, jobIdx);
}
void Schedule::markAsScheduled(pJob scheduledJob) {
	scheduledJobs.push_back(move(scheduledJob));
}

int Schedule::getNumberOfScheduledJobs() const {
	return scheduledJobs.size();
}

const Job* Schedule::getScheduledJob(size_t idx) const {
	if (idx >= scheduledJobs.size()) throw out_of_range("Schedule::getScheduledJob() out of range");
	return scheduledJobs[idx].get();
}

Operation* Schedule::findInScheduledJobs(Operation* remoteOp) const {
	for (size_t j = 0; j < scheduledJobs.size(); ++j) {
		if (scheduledJobs[j]->getId() == remoteOp->getId()) {
			for (size_t o = 0; o < scheduledJobs[j]->size(); ++o) {
				if ((*scheduledJobs[j])[o].getStg() == remoteOp->getStg()) {
					return &(*scheduledJobs[j])[o];	// local op
				}
			}
		}
	}
	return nullptr;
}
Operation* Schedule::findInUnscheduledJobs(Operation* remoteOp) const {
	for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
		if (unscheduledJobs[j]->getId() == remoteOp->getId()) {
			for (size_t o = 0; o < unscheduledJobs[j]->size(); ++o) {
				if ((*unscheduledJobs[j])[o].getStg() == remoteOp->getStg()) {
					return &(*unscheduledJobs[j])[o];	// local op
				}
			}
		}
	}
	return nullptr;
}

void Schedule::lSchedFirstJob(double pWait) {
	for (size_t op = 0; op < (*unscheduledJobs.begin())->size(); ++op) {
		schedOp(&(**unscheduledJobs.begin())[op], pWait);
	}
	shiftJobFromVecToVec(unscheduledJobs, scheduledJobs, 0);
}

void Schedule::lSchedFirstJobInflated(double pWait, double inflation, bool batchinStageInflationOnly, bool opsWithoutTcInflationOnly) {
	for (size_t op = 0; op < (*unscheduledJobs.begin())->size(); ++op) {
		schedOp(&(**unscheduledJobs.begin())[op], pWait, inflation, batchinStageInflationOnly, opsWithoutTcInflationOnly);
	}
	shiftJobFromVecToVec(unscheduledJobs, scheduledJobs, 0);
}
void Schedule::lSchedJobsStageWise(double pWait) {
	for (size_t i = 0; i < size(); ++i) {	
		
		for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
			Operation* op = &(*unscheduledJobs[j])[i];

			int maxLookAhead = workcenters[i]->getCap() - op->getS();	// looking ahead, ops may be batched if their common start time thus is earlier than their latest start if not batched
			int lookAheadCapRqrmt = op->getS();
			
			if (!op->isScheduled()) {
				vector<Operation*> lookingAhead = vector<Operation*>();
				double earliest = op->getEarliestStart();
				double earliestC = earliest + op->getP();
				double latest = min(earliest + op->getP(), op->getLatestStartConsideringBwdTc());
				double commonStart = earliest;

				if(maxLookAhead > 0) {
					// this operation should not wait if a full batch is ready to start at earliestC
					bool bFullBatchWaiting = false;
					/*int remainingCap = workcenters[i]->getCap();
					for (size_t k = j + 1; k < unscheduledJobs.size(); ++k) {
						Operation* lookAheadOp = &(*unscheduledJobs[k])[i];
						if (lookAheadOp->getF() == op->getF()) {
							if (!lookAheadOp->isScheduled()) {
								double tempEarliest = lookAheadOp->getEarliestStart();
								if (tempEarliest > earliest && tempEarliest <= earliestC) {
									remainingCap -= lookAheadOp->getS();
								}
							}
						}
						if (remainingCap <= 0) {
							bFullBatchWaiting = true;
							break;
						}
					}*/

					// look for operations worthy to wait for
					if (!bFullBatchWaiting) {
						for (size_t k = j + 1; k < unscheduledJobs.size(); ++k) {
							Operation* lookAheadOp = &(*unscheduledJobs[k])[i];
							if (lookAheadOp->getF() == op->getF()) {
								if (!lookAheadOp->isScheduled()) {
									double tempEarliest = lookAheadOp->getEarliestStart();
									if (tempEarliest <= latest) {
										if (workcenters[i]->getCap() >= lookAheadCapRqrmt + lookAheadOp->getS()) {
											// CANDIDATE FOUND
											lookAheadCapRqrmt += lookAheadOp->getS();
											// TODO: TRY MORE SOPHISTICATED DECISIONS (e.g. based on weights or anticipated wT)
											commonStart = max(earliest, tempEarliest);
											lookingAhead.push_back(lookAheadOp);
										}
									}
								}
							}
							if (lookingAhead.size() >= maxLookAhead || lookAheadCapRqrmt >= workcenters[i]->getCap()) break;
						}
					}
				}	
				schedOpDelayed(op, commonStart);
				for (size_t i = 0; i < lookingAhead.size(); ++i) {
					schedOpDelayed(lookingAhead[i], commonStart);
				}
			}
		}
	}
	while (!unscheduledJobs.empty()) {
		shiftJobFromVecToVec(unscheduledJobs, scheduledJobs, 0);
	}
	/*if (!this->isValid()) {
		TCB::logger.Log(Error, "Exception thrown in Schedule::lSchedJobsStageWise(...)");
		throw ExcSched("ERROR: invalid schedule after Schedule::lSchedStages.");
	}*/
}
void Schedule::lSchedJobsStageWiseWithSorting(prioRule<pJob> rule, double pWait) {
	rule(unscheduledJobs);
	lSchedJobsStageWise(pWait);
}
void Schedule::lSchedJobsStageWiseBackward(double pWait) {
	// temporarily change all jobs´ release time to an upper bound value
	/*double tZero = problem->getUpperBoundMSP();
	for (size_t i = 0; i < unscheduledJobs.size(); ++i) {
		unscheduledJobs[i]->setR(unscheduledJobs[i]->getR() + tZero);
	}*/
	
	// +++ A) last stage: forward order at last stage, as early as possible +++
	size_t last = size() - 1;
	for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
		Operation* op = &(*unscheduledJobs[j])[last];
		int maxLookAhead = workcenters[last]->getCap() - op->getS();	// looking ahead, ops may be batched if their common start time thus is earlier than their latest start if not batched
		int lookAheadCapRqrmt = op->getS();
		if (!op->isScheduled()) {
			vector<Operation*> lookingAhead = vector<Operation*>();
			double earliest = op->getEarliestStartForBackwardScheduling(); // op->getEarliestStart();
			double earliestC = earliest + op->getP();
			double latest = earliest + op->getP(); // , op->getLatestStartConsideringFwdTc());
			if (op->getSucc() != nullptr) {
				if (op->getSucc()->isScheduled()) {
					latest = op->getSucc()->getStart() - op->getP();
				}
			}
			double commonStart = earliest;
			if (maxLookAhead > 0) {
				// this operation should not wait if a full batch is ready to start at earliestC
				bool bFullBatchWaiting = false;
				/*int remainingCap = workcenters[i]->getCap();
				for (size_t k = j + 1; k < unscheduledJobs.size(); ++k) {
					Operation* lookAheadOp = &(*unscheduledJobs[k])[i];
					if (lookAheadOp->getF() == op->getF()) {
						if (!lookAheadOp->isScheduled()) {
							double tempEarliest = lookAheadOp->getEarliestStart();
							if (tempEarliest > earliest && tempEarliest <= earliestC) {
								remainingCap -= lookAheadOp->getS();
							}
						}
					}
					if (remainingCap <= 0) {
						bFullBatchWaiting = true;
						break;
					}
				}*/

				// look for operations worthy to wait for
				if (!bFullBatchWaiting) {
					for (size_t k = j + 1; k < unscheduledJobs.size(); ++k) {
						Operation* lookAheadOp = &(*unscheduledJobs[k])[last];
						if (lookAheadOp->getF() == op->getF()) {
							if (!lookAheadOp->isScheduled()) {
								double tempEarliest = lookAheadOp->getEarliestStart();
								if (tempEarliest <= latest) {
									if (workcenters[last]->getCap() >= lookAheadCapRqrmt + lookAheadOp->getS()) {
										// CANDIDATE FOUND
										lookAheadCapRqrmt += lookAheadOp->getS();
										// TODO: TRY MORE SOPHISTICATED DECISIONS (e.g. based on weights or anticipated wT)
										commonStart = max(earliest, tempEarliest);
										lookingAhead.push_back(lookAheadOp);
									}
								}
							}
						}
						if (lookingAhead.size() >= maxLookAhead || lookAheadCapRqrmt >= workcenters[last]->getCap()) break;
					}
				}
			}
			schedOpDelayed(op, commonStart);
			for (size_t i = 0; i < lookingAhead.size(); ++i) {
				schedOpDelayed(lookingAhead[i], commonStart);
			}
		}
	}

	// reset all jobs´ release time to their original values
/*	for (size_t i = 0; i < unscheduledJobs.size(); ++i) {
		unscheduledJobs[i]->setR(unscheduledJobs[i]->getR() - tZero);
	}*/	

	// +++ B) from 2nd-to-last to 1st stage: backward order, as late as possible +++
	for (int i = size() - 2; i >= 0; --i) {
		for (int j = unscheduledJobs.size() - 1; j >= 0; --j) {
			Operation* op = &(*unscheduledJobs[j])[i];
			int maxLookAhead = workcenters[i]->getCap() - op->getS();	// looking ahead, ops may be batched if their common start time thus is earlier than their latest start if not batched
			int lookAheadCapRqrmt = op->getS();
			bool bSuccessorMovable = false;
			if (!op->isScheduled()) {

				if (op->getId() == 5 && op->getStg() == 2) {
					int stop = 666;
				}

				vector<Operation*> lookingAhead = vector<Operation*>();
				double earliest = op->getEarliestStartForBackwardScheduling(); // op->getEarliestStart();
				double latest = op->getSucc()->getStart() - op->getP();
				double commonStart = latest;

				double latestAtWc = workcenters[i]->findLatestAvailableTimeSlotBefore(latest, op->getP());
				if (latestAtWc < latest) {
					double necessaryLeftShift = latest - latestAtWc;
					double maximalLeftShift = latest - earliest;
					// TODO Try to postpone successors without turmoil
					vector<pair<double, double>> possibleLeftShift = getLeftShiftOptions(op->getSucc());
					for (size_t p = 0; p < possibleLeftShift.size(); ++p) {
						if ((possibleLeftShift[p].first >= necessaryLeftShift && possibleLeftShift[p].second <= necessaryLeftShift) 
							|| (possibleLeftShift[p].first == possibleLeftShift[p].second && possibleLeftShift[p].first >= necessaryLeftShift)) {
							if (possibleLeftShift[p].second <= maximalLeftShift) {
								// left shift successor
								Operation* leftShiftedOp = op->getSucc();
								size_t leftShiftOpIdx = leftShiftedOp->getIdxInBatch();
								Batch* leftShiftBatch = leftShiftedOp->getBatch();
								size_t leftShiftBatchIdx = leftShiftBatch->getIdx();
								Machine* leftShiftMac = leftShiftBatch->getMachine();
								leftShiftBatch->removeOp(leftShiftOpIdx);
								if (leftShiftBatch->isEmpty()) {
									leftShiftMac->removeBatch(leftShiftBatchIdx);
								}
								schedOpDelayed(leftShiftedOp, leftShiftBatch->getStart() - necessaryLeftShift);

								// update common start
								commonStart = latestAtWc;
								break;
							}
						}
					}
				}

				if (commonStart < earliest) {
					int stop = 666;
				}

				if (maxLookAhead > 0) {
					// this operation should not wait if a full batch is ready to start at earliestC
					bool bFullBatchWaiting = false;
					/*int remainingCap = workcenters[i]->getCap();
					for (int k = j - 1; j >= 0; --k) {
						Operation* lookAheadOp = &(*unscheduledJobs[k])[i];
						if (lookAheadOp->getF() == op->getF()) {
							if (!lookAheadOp->isScheduled()) {
								double tempEarliest = lookAheadOp->getEarliestStart();
								if (tempEarliest > earliest && tempEarliest <= earliestC) {
									remainingCap -= lookAheadOp->getS();
								}
							}
						}
						if (remainingCap <= 0) {
							bFullBatchWaiting = true;
							break;
						}
					}*/
					if (!bFullBatchWaiting) {
						for (int k = j - 1; k >= 0; --k) {
							Operation* lookAheadOp = &(*unscheduledJobs[k])[i];
							if (lookAheadOp->getF() == op->getF()) {
								if (!lookAheadOp->isScheduled()) {
									double tempLatest = lookAheadOp->getSucc()->getStart() - lookAheadOp->getP();

									// if the lookAhead operation was not considered would that lead to a tc violation?
									//bool bConsiderTemp = false;
									//double tempStartIfNotConsidered = commonStart - op->getP();	// MISTAKE: this is on valid for single machines, TODO find more general test
									//const vector<pair<int, double>> tempFwdTc = lookAheadOp->getTcMaxFwd();
									//for (size_t tc = 0; tc < tempFwdTc.size(); ++tc) {
									//	if (tempFwdTc[tc].second < 999999 && tempStartIfNotConsidered + tempFwdTc[tc].second < (*lookAheadOp->getJob())[tempFwdTc[tc].first].getStart()) {
									//		bConsiderTemp = true;
									//	}
									//}
									bool bConsiderTemp = true;
									
									if (bConsiderTemp && tempLatest >= earliest) {
										if (workcenters[i]->getCap() >= lookAheadCapRqrmt + lookAheadOp->getS()) {
											// CANDIDATE FOUND
											lookAheadCapRqrmt += lookAheadOp->getS();
											// TODO: TRY MORE SOPHISTICATED DECISIONS (e.g. based on weights or anticipated wT)
											commonStart = min(commonStart, tempLatest);
											lookingAhead.push_back(lookAheadOp);
										}
									}
								}
							}
							else {
								// stop looking for batch candidates outside of immediate successors in the sequence
								break;
							}
							if (lookingAhead.size() >= maxLookAhead || lookAheadCapRqrmt >= workcenters[last]->getCap()) break;
						}
					}
				}
				schedOpDelayed(op, commonStart);
				for (size_t i = 0; i < lookingAhead.size(); ++i) {
					schedOpDelayed(lookingAhead[i], commonStart);
				}
			}
		}
	}

	/*if (!this->isValid()) {
		TCB::logger.Log(Error, "Exception thrown in Schedule::lSchedStages(...)");
		throw ExcSched("ERROR: invalid schedule after Schedule::lSchedStages.");
	}*/


	// probably left shifting local search in the end  
}
void Schedule::lSchedJobsStageWiseBackwardWithSorting(prioRule<pJob> rule, double pWait) {
	rule(unscheduledJobs);
	lSchedJobsStageWiseBackward(pWait);	
}
void Schedule::lSchedJobs(double pWait) {
	while(!unscheduledJobs.empty()) {
		lSchedFirstJob(pWait);
	}
}
void Schedule::lSchedJobs(vector<double> pWaitVec) {
	double bestTWT = DBL_MAX; // numeric_limits<double>::max();
	double bestWait = pWaitVec[0];
	for (size_t w = 0; w < pWaitVec.size(); ++w) {
		unique_ptr<Schedule> copySched = this->clone();
		copySched->lSchedJobs(pWaitVec[w]);
		double tempTWT = copySched->getTWT();
		if (tempTWT < bestTWT) {
			bestTWT = tempTWT;
			bestWait = pWaitVec[w];
		}
	}
	//TCB::logger.Log(Info, "List Scheduling with best pWait " + to_string(bestWait));
	lSchedJobs(bestWait);
}
void Schedule::lSchedJobsInflated(double pWait, double inflation, bool batchinStageInflationOnly, bool opsWithoutTcInflationOnly) {
	while (!unscheduledJobs.empty()) {
		lSchedFirstJobInflated(pWait, inflation, batchinStageInflationOnly, opsWithoutTcInflationOnly);
	}
	leftShiftBatches();	// deflate
}

void Schedule::lSchedJobsWithSorting(prioRule<pJob> rule, double pWait) {
	rule(unscheduledJobs);
	lSchedJobs(pWait);
}
void Schedule::lSchedJobsWithSorting(prioRule<pJob> rule, Sched_params& sched_params) {
	vector<double> pWaitVec = getDoubleGrid(sched_params.pWaitLow, sched_params.pWaitHigh, sched_params.pWaitStep);
	if (pWaitVec.size() == 1) {
		lSchedJobsWithSorting(rule, pWaitVec[0]);
	} else {
		double bestTWT = DBL_MAX; // numeric_limits<double>::max();
		double bestWait = pWaitVec[0];
		for (size_t w = 0; w < pWaitVec.size(); ++w) {
			unique_ptr<Schedule> copySched = this->clone();
			copySched->lSchedJobsWithSorting(rule, pWaitVec[w]);
			double tempTWT = copySched->getTWT();
			if (tempTWT < bestTWT) {
				bestTWT = tempTWT;
				bestWait = pWaitVec[w];
			}
		}
		lSchedJobsWithSorting(rule, bestWait);
	}
}

void Schedule::lSchedJobsWithSorting(prioRuleKappa<pJob> rule, double kappa, double pWait) {
	double t = 0.0;	// dynamic computation of priority index (increase t)
	while (!unscheduledJobs.empty()) {
		t = getMinMSP(0);
		rule(unscheduledJobs, t, kappa);
		lSchedFirstJob(pWait);
	}
}

double Schedule::lSchedJobsWithSorting(prioRuleKappa<pJob> rule, const vector<double>& kappaGrid, double pWait, objectiveFunction objectiveFunction) {
	if (kappaGrid.size() == 1) {
		lSchedJobsWithSorting(rule, kappaGrid[0], pWait);
		return kappaGrid[0];
	}
	
	double bestObjectiveValue = DBL_MAX; //  numeric_limits<double>::max();
	double bestKappa = 0.0;
	for (size_t kappa = 0; kappa < kappaGrid.size(); ++kappa) {
		lSchedJobsWithSorting(rule, kappaGrid[kappa], pWait);
		double tempObjectiveValue = objectiveFunction(this);
		if (tempObjectiveValue < bestObjectiveValue) {
			bestObjectiveValue = tempObjectiveValue;
			bestKappa = kappaGrid[kappa];
		}
		reset();
	}
	lSchedJobsWithSorting(rule, bestKappa, pWait);
	//TCB::logger.Log(Info, "Found a schedule with best kappa value = " + to_string(bestKappa));
	return bestKappa;
}

double Schedule::lSchedJobsWithSorting(prioRuleKappa<pJob> rule, Sched_params& sched_params, objectiveFunction objectiveFunction) {
	vector<double> pWaitVec = getDoubleGrid(sched_params.pWaitLow, sched_params.pWaitHigh, sched_params.pWaitStep);
	vector<double> kappas = getDoubleGrid(sched_params.kappaLow, sched_params.kappaHigh, sched_params.kappaStep);

	double bestTWT = DBL_MAX; //  numeric_limits<double>::max();
	double bestWait = pWaitVec[0];
	for (size_t w = 0; w < pWaitVec.size(); ++w) {
		unique_ptr<Schedule> copySched = this->clone();
		copySched->lSchedJobsWithSorting(rule, kappas, pWaitVec[w]);
		double tempTWT = copySched->getTWT();
		if (tempTWT < bestTWT) {
			bestTWT = tempTWT;
			bestWait = pWaitVec[w];
		}
	}
	return lSchedJobsWithSorting(rule, kappas, bestWait);
}

void Schedule::lSchedJobsWithRandomKeySorting(prioRuleKeySet<pJob> rule, const std::vector<double>& keys, double pWait) {
	rule(unscheduledJobs, keys);
	lSchedJobs(pWait);
}
void Schedule::lSchedJobsWithRandomKeySorting(prioRuleKeySet<pJob> rule, const std::vector<double>& keys, Sched_params& sched_params) {
	vector<double> pWaitVec = getDoubleGrid(sched_params.pWaitLow, sched_params.pWaitHigh, sched_params.pWaitStep);
	double bestTWT = DBL_MAX; // numeric_limits<double>::max();
	double bestWait = pWaitVec[0];
	for (size_t w = 0; w < pWaitVec.size(); ++w) {
		unique_ptr<Schedule> copySched = this->clone();
		copySched->lSchedJobsWithRandomKeySorting(rule, keys, pWaitVec[w]);
		double tempTWT = copySched->getTWT();
		if (tempTWT < bestTWT) {
			bestTWT = tempTWT;
			bestWait = pWaitVec[w];
		}
	}
	lSchedJobsWithRandomKeySorting(rule, keys, bestWait);
}

void Schedule::lSchedGifflerThompson(prioRule<pJob> rule, double pWait) {
	// 1) consider all "next" operations of jobs to be processed
	vector<queue<Operation*>> unscheduledOps = vector<queue<Operation*>>(unscheduledJobs.size());
	for (size_t j = 0; j < unscheduledJobs.size(); ++j) {
		for (size_t o = 0; o < (*unscheduledJobs[j]).size(); ++o) {
			unscheduledOps[j].push(&(*unscheduledJobs[j])[o]);
		}
	}

	// 2) get earliest completion time at any suitable workcenter => op*: operation with min C, wc*: corresponding workcenter 
	bool allOpsScheduled = false;
	double earliestC = DBL_MAX; // numeric_limits<double>::max();
	while (!allOpsScheduled) {
		size_t nextOp = 0;
		allOpsScheduled = true;
		for (size_t j = 0; j < unscheduledOps.size(); ++j) {
			if (!unscheduledOps[j].empty()) {
				allOpsScheduled = false;

				size_t bestMacIdx = 0;
				size_t bestBatIdx = 0;
				double bestStart = 0.0;
				//double tempC getEarliestC(unscheduledOps[j].front());

			}
		}
	} 
	
	// 3) further consider all operations which can be started before the completion of op* at wc* (Overlapping set)

	// Original Giffler & Thompson algorithm: branch (construct all possible schedules with any of the ops of 3) to be processed next and loop back to 1))
	// List scheduling adaptation: assign op from set 3) by priority index, loop back to 1)


}

void Schedule::leftShiftBatches() {
	bool bImproved = true;
	while (bImproved) {
		bImproved = false;
		for (size_t wc = 0; wc < size(); ++wc) {
			for (size_t i = 0; i < workcenters[wc]->size(); ++i) {
				for (size_t b = 0; b < (*workcenters[wc])[i].size(); ++b) {
					if (workcenters[wc]->leftShift(i, b)) {
						bImproved = true;
					}
				}
			}
		}
	}
}

void Schedule::localSearchBatchLeftShifting() {
	cout << "Schedule::localSearchBatchLeftShifting() not yet implemented. " << endl;

}

double Schedule::localSearchEvaluateBatchConsolidation(size_t idxWc, size_t tgtMac, size_t tgtBatch, size_t srcMac, size_t srcBatch, size_t& opIdx) {
	double change = 0.0;
	if (idxWc < 0 || idxWc >= size()) throw out_of_range("ERROR: Schedule::localSearchEvaluateBatchConsolidation(...) idxWc out of bounds!");
	if (tgtMac < 0 || tgtMac >= workcenters[idxWc]->size()) throw out_of_range("ERROR: Schedule::localSearchEvaluateBatchConsolidation(...) idx1stMac out of bounds!");
	if (srcMac < 0 || srcMac >= workcenters[idxWc]->size()) throw out_of_range("ERROR: Schedule::localSearchEvaluateBatchConsolidation(...) idx2ndMac out of bounds!");
	if (tgtBatch < 0 || tgtBatch >= (*workcenters[idxWc])[tgtMac].size()) throw out_of_range("ERROR: Schedule::localSearchEvaluateBatchConsolidation(...) idx1stBatch out of bounds!");
	if (srcBatch < 0 || srcBatch >= (*workcenters[idxWc])[srcMac].size()) throw out_of_range("ERROR: Schedule::localSearchEvaluateBatchConsolidation(...) idx2ndBatch out of bounds!");

	Workcenter* wc = &(*workcenters[idxWc]);
	Machine* mac1 = &(*wc)[tgtMac];
	Machine* mac2 = &(*wc)[srcMac];
	Batch* bat1 = &(*mac1)[tgtBatch];
	Batch* bat2 = &(*mac2)[srcBatch];

	// [JR-2025-Dec-16] Assessment of batch consolidation´s implications is too cumbersome, instead execute change to a copy of the schedule
	//if (bat1->getStart() >= bat2->getStart() || bat1->getF() != bat2->getF()) return 0.0; // symmetry braking and disregarding of incompatible batches
	//for (size_t j = 0; j < bat2->size(); ++j) {
	//	Operation* op = &(*bat2)[j];
	//	double opR = op->getAvailability();
	//	// this move is only relevant if op can be left shifted (< current start) and it cannot simply be inserted into previous batch (otherwise an insertion/shift move would do the trick)
	//	if (opR < bat2->getStart() && opR > bat1->getStart()) {
	//		
	//		if (op->getS() <= bat1->getAvailableCap()) {
	//			double leftShiftAtThisStage = bat2->getStart() - opR;
	//			vector<pair<Operation*, double>> rightShifts = vector <pair<Operation*, double>>();		
	//		}
	//		// WORK IN PROGRESS
	//	}
	//}

	double myTWT = getTWT();
	double bestTWT = myTWT;

	for (size_t movingOpIdx = 0; movingOpIdx < bat2->size(); ++movingOpIdx) {
		Operation* op = &(*bat2)[movingOpIdx];
		double opR = op->getAvailability();
		if (opR < bat2->getStart() && opR > bat1->getStart()) {	// this move is only relevant if op can be left shifted (< current start) and it cannot simply be inserted into previous batch (otherwise an insertion/shift move would do the trick)
			if (op->getS() <= bat1->getAvailableCap()) {
				unique_ptr<Schedule> copySched = clone();
				if (copySched->localSearchConsolidateBatch(idxWc, tgtMac, tgtBatch, srcMac, srcBatch, movingOpIdx)) {
					double tempTWT = copySched->getTWT();
					if (tempTWT < bestTWT) {
						bestTWT = tempTWT;
						opIdx = movingOpIdx;
					}
				}	
			}
		}
	}

	return myTWT - bestTWT;
}

bool Schedule::localSearchConsolidateBatch(size_t wcIdx, size_t tgtMacIdx, size_t tgtBatchIdx, size_t srcMacIdx, size_t srcBatchIdx, size_t opIdx) {
	Workcenter* wc = workcenters[wcIdx].get();
	Machine* srcMac = &(*wc)[srcMacIdx];
	Machine* tgtMac = &(*wc)[tgtMacIdx];
	Batch* srcBat = &(*srcMac)[srcBatchIdx];
	Batch* tgtBat = &(*tgtMac)[tgtBatchIdx];
	Operation* op = &(*srcBat)[opIdx];
	Job* job = op->getJob();

	// right shift target batch and transfer operation to target batch
	double rightShiftTgt = op->getAvailability() - tgtBat->getStart();
	double currentStart = (*tgtMac)[tgtBatchIdx].getStart();
	double originalCompletion = (*tgtMac)[tgtBatchIdx].getC();
	(*tgtMac)[tgtBatchIdx].setStart(currentStart + rightShiftTgt, false);
	srcBat->removeOp(opIdx);
	if (!tgtBat->addOp(op)) {
		return false;
	}

	// there are only right shifts at this stage if the source batch carrys only one operation!
	if(srcBat->size() > 0) {
		// right shifting of target Batch and successive batches at target machine
		size_t nNumberOfBatchesAtTarget = tgtMac->size();
		for (size_t b = tgtBatchIdx + 1; b < nNumberOfBatchesAtTarget; ++b) {
			double rightShiftOp = (*tgtMac)[b - 1].getC() - (*tgtMac)[b].getStart();
			if (rightShiftOp <= 0) break;																// [JR-2026-Jan-12] changed rightShiftTgt to rightShiftOp
			double currentStart = (*tgtMac)[b].getStart();
			(*tgtMac)[b].setStart(currentStart + rightShiftOp, false);	// no validity check here		// [JR-2026-Jan-12] changed rightShiftTgt to rightShiftOp
		}

		// right shifting at successive stages
		for (size_t o = (wcIdx + 1); o < size(); ++o) {
			Workcenter* currentWc = workcenters[o].get();
			for (size_t m = 0; m < currentWc->size(); ++m) {
				Machine* currentMac = &(*currentWc)[m];
				for (size_t b = 0; b < currentMac->size(); ++b) {
					Batch* currentBat = &(*currentMac)[b];
					if (currentBat->getStart() < originalCompletion) {
						// if a batch was started before the left shifted batch, it cannot have been affected
						continue;
					}
					// Batches from here on may have been affected, successive batches may still be affected even if 
					
					double rightShiftOp = currentBat->getR() - currentBat->getStart();	// at the same machine
					// TODO (MAYBE?): instead of stupid right shifting, parallel machines should be considered

					if (b > 0) {
						rightShiftOp = max(rightShiftOp, (*currentMac)[b - 1].getC() - currentBat->getStart());
					}
					if (rightShiftOp > 0) {
						double currentBatStart = currentBat->getStart();
						currentBat->setStart(currentBatStart + rightShiftOp, false);
					}
				}
			}
		}
	} else {
		// delete srcbatch after job insertion into target
		size_t srcBatIdx = srcBat->getIdx();
		srcMac->removeBatch(srcBatIdx);
		// TODO this may have enabled additional left shifts
	}

	// left shifting of the successors of the left shifted operation
	for (size_t s = 0; s < size() - wcIdx; ++s) {
		Operation* succOp = op->getSucc(s);
		if (succOp != nullptr) {
			double maxLeftShift = max(0, succOp->getStart() - succOp->getAvailability());
			if (maxLeftShift > 0) {
				Batch* succBatch = succOp->getBatch();
				double originalSuccStart = succBatch->getStart();
				Machine* succMac = succBatch->getMachine();
				size_t tempBatIdx = succBatch->getIdx();
				bool bOnlyOp = succBatch->size() <= 1;
				bool bDiscreteMove;	// if true, succOp is inserted into existing batch
				vector<pair<double, double>> possibleLeftShifts = getLeftShiftOptions(succOp);
				for (size_t opt = 0; opt < possibleLeftShifts.size(); ++opt) {
					if ((possibleLeftShifts[opt].first == possibleLeftShifts[opt].second && possibleLeftShifts[opt].first <= maxLeftShift + TCB::precision)
						||(possibleLeftShifts[opt].first >= possibleLeftShifts[opt].second && possibleLeftShifts[opt].second <= maxLeftShift + TCB::precision) ) {
						bDiscreteMove = executeLeftShiftOption(succOp, possibleLeftShifts[opt]);
						break;
					}
				}

				if ((bDiscreteMove && bOnlyOp) || !bDiscreteMove) {
					// in this case left shifting succeeding batches at this machine may have become feasible
					for (size_t b = 0; b < succMac->size(); ++b) {
						// consider all batches starting later than original succBatch start
						double tempStart = (*succMac)[b].getStart();
						if (tempStart < originalSuccStart) {
							continue;
						}
						double tempR = (*succMac)[b].getR();
						double tempLeftShift = tempStart - tempR;
						if (b > 0) {
							tempLeftShift = min(tempLeftShift, tempStart - (*succMac)[b - 1].getC());
						}
						if (tempLeftShift > 0) {
							(*succMac)[b].setStart(tempStart - tempLeftShift, false);
						}
					}
				}		
			}	
		}
	}

	// Eventually, right shifts may have caused tc violations => further right shifts may be necessary
	for (int xWc = size() - 1; xWc >= 0; --xWc) {
		Workcenter* xWorkcenter = workcenters[xWc].get();
		for (size_t xM = 0; xM < xWorkcenter->size(); ++xM) {
			Machine* xMachine = &(*xWorkcenter)[xM];
			for (int xB = xMachine->size() - 1; xB >= 0; --xB) {
				Batch* xBatch = &(*xMachine)[xB];
				for (size_t xOp = 0; xOp < xBatch->size(); ++xOp) {
					xWorkcenter->ensureValidityFixedBatchFormation(&(*xBatch)[xOp]);
				}
			}
		}
	}

	return isValid();
}

void Schedule::localSearchOpLeftShifting(prioRule<pJob> rule, double pWait) {
	bool bImproved = true;
	while (bImproved) {
		bImproved = false;
		rule(scheduledJobs);
		for (size_t j = 0; j < scheduledJobs.size(); ++j) {
			for (size_t o = 0; o < (*scheduledJobs[j]).size(); ++o) {
				size_t wcIdx = (*scheduledJobs[j])[o].getWorkcenterId() - 1;
				size_t mIdx = 0;
				size_t bIdx = 0;
				size_t jIdx = 0;
				if(workcenters[wcIdx]->locateOp(&(*scheduledJobs[j])[o], mIdx, bIdx, jIdx)) {
					if (workcenters[wcIdx]->leftShift(mIdx, bIdx, jIdx, pWait)) {
						bImproved = true;
					}
				}
			}
			/*for (int o = (*scheduledJobs[j]).size() - 1; o >= 0; --o) {
				size_t wcIdx = (*scheduledJobs[j])[o].getWorkcenterId() - 1;
				size_t mIdx = 0;
				size_t bIdx = 0;
				size_t jIdx = 0;
				if (workcenters[wcIdx]->locateOp(&(*scheduledJobs[j])[o], mIdx, bIdx, jIdx)) {
					if (workcenters[wcIdx]->leftShift(mIdx, bIdx, jIdx, pWait)) {
						bImproved = true;
					}
				}
			}*/
		}
	}
}
void Schedule::perturbRandomJobSwap() {
	uniform_int_distribution<> distrib(0, scheduledJobs.size()-1);
	int max = 1000;
	int i = 0;
	int j = distrib(TCB::rng);
	int k = distrib(TCB::rng);
	while (i < max && (j == k || scheduledJobs[j]->getF() != scheduledJobs[k]->getF() || scheduledJobs[j]->getStart() < scheduledJobs[k]->getR() || scheduledJobs[k]->getStart() < scheduledJobs[j]->getR())) {
		j = distrib(TCB::rng);
		k = distrib(TCB::rng);
		++i;
	}
	locSearchSwapJobs(j, k);
}
void Schedule::perturbRandomBatchRightShifting(){
	// 1) randomly choose batch at first stage (TODO: maybe also consider laters stages)
	int wcIdx = 0;	
	uniform_int_distribution<> mDistrib(0, workcenters[0]->size() - 1);
	int mIdx = mDistrib(TCB::rng);
	Machine* machine = &(*workcenters[0])[mIdx];
	uniform_int_distribution<> bDistrib(0, machine->size() - 1);
	int bIdx = bDistrib(TCB::rng);

	// 2) define duration of right-shift: equivalent to a randomly chosen product´s processing time
	uniform_int_distribution<> pDistrib(0, problem->getF() - 1);
	int productIdx = pDistrib(TCB::rng);
	double rightShiftAmount = problem->getProduct(productIdx)->getP(wcIdx);

	// 3) Actually right-shift
	double tempRightShift = rightShiftAmount;
	for (int b = bIdx; b < machine->size(); ++b) {
		Batch* batch = &(*machine)[b];
		batch->setStart(batch->getStart() + tempRightShift, false);
		if (b < machine->size() - 1) {
			tempRightShift = min(tempRightShift, max(batch->getC() - (*machine)[b+1].getStart(), 0));
		}
	}


	// 4) ensure validity for all shifted operations
	for (int b = machine->size() - 1; b >= bIdx; --b) {
		Batch* batch = &(*machine)[b];
		for (size_t op = 0; op < batch->size(); ++op) {
			Operation* operation = &(*batch)[op];
			workcenters[0]->ensureValidityFixedBatchFormation(operation);
		}
	}

	//if (!isValid()) {
	//	TCB::logger.Log(Error, "Exception thrown in Schedule::perturbRandomRightShifting()");
	//	throw ExcSched("ERROR in Schedule::perturbRandomRightShifting()");
	//}

}
void Schedule::perturbRandomJobRightShifting() {
	// 1) randomly choose batch and job at first stage
	int wcIdx = 0;
	uniform_int_distribution<> mDistrib(0, workcenters[0]->size() - 1);
	int mIdx = mDistrib(TCB::rng);
	Machine* machine = &(*workcenters[0])[mIdx];	// machines may be empty!
	while (machine->size() == 0) {
		mIdx = mDistrib(TCB::rng);
		machine = &(*workcenters[0])[mIdx];
	}
	uniform_int_distribution<> bDistrib(0, machine->size() - 1);
	int bIdx = bDistrib(TCB::rng);
	Batch* batch = &(*machine)[bIdx];
	uniform_int_distribution<> opDistrib(0, batch->size() - 1);
	int opIdx = opDistrib(TCB::rng);
	Operation* movingOp = &(*batch)[opIdx];

	uniform_int_distribution<> productDistrib(0, problem->getF()-1);
	int productIdx = productDistrib(TCB::rng);
	double minDelay = problem->getProduct(productIdx)->getP(0);	// right-shift by the amount of a randomly chosen product´s processing time (to facilitate an operation of such product to fill in the gap during next local search)

	double newStart = movingOp->getStart() + minDelay;
	batch->removeOp(opIdx);
	if (batch->isEmpty()) {
		machine->removeBatch(bIdx);
	}
	schedOpDelayed(movingOp, newStart);
}
bool Schedule::localSearchJobSwapping(prioRule<pJob> rule, bool bestFit) {
	bool bImprovedOnce = false;
	int debug = 0;
	bool bImproved = true;
	int nJobs = scheduledJobs.size();
	while (bImproved) {

		debug++;

		bImproved = false;
		rule(scheduledJobs);
		bool swapFeasible = true;
		double bestImprovement = 0;
		size_t best1 = 0;
		size_t best2 = 0;
		for (size_t j = 0; j < nJobs; ++j) {
			for (int k = nJobs - 1; k >= 0; --k) {
				if (j != k) {
					double tempImprovement = locSearchEvaluateJobSwap(j, k, swapFeasible);
					if (swapFeasible && tempImprovement > bestImprovement) {
						if (!bestFit) { // FIRST FIT
							locSearchSwapJobs(j, k);
							bImprovedOnce = true;
							bImproved = true;
							break;
						} else {			// BEST FIT
							bestImprovement = tempImprovement;
							best1 = j;
							best2 = k;
						}
					}
				}
			}
			if (bImproved) {
				break; // continue with while loop
			}
		}
		if (bestFit && bestImprovement > 0) {
			locSearchSwapJobs(best1, best2);
			bImprovedOnce = true;
			bImproved = true;
		}
	}
	return bImprovedOnce;
}
bool Schedule::localSearchJobLeftShifting(prioRule<pJob> rule, bool bestFit) {
	int maxConsecutiveInnerStageImprovement = scheduledJobs.size();	// to avoid infinite loop (there may be infinite shifts improving completion of inner operations only (w/o effect on TWT)
	int nConsecutiveInnerStageImprovement = 0;
	bool bImproved = true;
	bool bImprovedOnce = false;
	int nJobs = scheduledJobs.size();
	
	while (bImproved) {
		bImproved = false;
		rule(scheduledJobs);
		pair<double, double> bestImprovement = make_pair(0, 0);
		double bestTwtImprovement = 0.0;
		size_t best = 0;
		vector<vector<pair<double, double>>> bestPossibleLeftShifts = vector<vector<pair<double, double>>>(size());
		for (size_t j = 0; j < nJobs; ++j) {
			vector<vector<pair<double, double>>> tempPossibleLeftShifts = vector<vector<pair<double, double>>>(size());
			pair<double, double> tempImprovement = locSearchEvaluateJobLeftShift(j, tempPossibleLeftShifts);
			double tempTwtImprovement = tempImprovement.first * scheduledJobs[j]->getW();
			if (tempTwtImprovement > bestTwtImprovement || (tempTwtImprovement == bestTwtImprovement && tempImprovement.second > bestImprovement.second)) {
				if (tempTwtImprovement <= 0) {
					++nConsecutiveInnerStageImprovement;
				} else {
					nConsecutiveInnerStageImprovement = 0;
				}
				if (!bestFit) {
					locSearchLeftShiftJob(j, tempPossibleLeftShifts);
					if (nConsecutiveInnerStageImprovement < maxConsecutiveInnerStageImprovement) {
						bImprovedOnce = true;
						bImproved = true;
					}
					break;
				}
				else {
					bestImprovement = tempImprovement;
					bestTwtImprovement = tempTwtImprovement;
					bestPossibleLeftShifts = tempPossibleLeftShifts;
					best = j;
				}
			}
		}
		if (bestFit && (bestImprovement.first > 0 || bestImprovement.second > 0)) {
			locSearchLeftShiftJob(best, bestPossibleLeftShifts);
			if (nConsecutiveInnerStageImprovement < maxConsecutiveInnerStageImprovement) {
				bImprovedOnce = true;
				bImproved = true;
			}
		}
	}
	return bImprovedOnce;
}
bool Schedule::localSearchBatchConsolidation(bool bestFit) {
	bool bImproving = true;
	bool bImprovedOnce = false;
	while (bImproving) {
		bImproving = false;
		double bestImprovement = 0.0;
		size_t bestWc = 0;
		size_t bestM = 0;
		size_t bestB = 0;
		size_t bestOp = 0;

		for (size_t o = 0; o < size(); ++o) {
			Workcenter* workcenter = &(*workcenters[o]);
			for (size_t m = 0; m < workcenter->size(); ++m) {
				Machine* machine = &(*workcenter)[m];

				for (int b = 0; (b + 1) < machine->size(); ++b) {		// beware of trap: ->size() returns size_t, so ->size() - 1 would become a big value for ->size() == 0 !!!
						Batch* batch = &(*machine)[b];
						Batch* succB = &(*machine)[b + 1];	

						if (batch->getF() == succB->getF()) {
							for (size_t op = 0; op < succB->size(); ++op) {
								size_t opIdx = op;	// by value			
								double tempImprovement = localSearchEvaluateBatchConsolidation(o, m, b, m, b+1, opIdx);	// opIdx by reference
								
								if (tempImprovement > bestImprovement) {
									if (!bestFit) {
										// FIRST FIT => execute if improvement > 0
										localSearchConsolidateBatch(o, m, b, m, b+1, opIdx);
										bImprovedOnce = true;
										bImproving = true;
										break;
									} else {
										// BEST FIT (GREEDY) => keep looking for best move
										bestWc = o;
										bestM = m;
										bestB = b;
										bestOp = opIdx;
										bestImprovement = tempImprovement;
									}
								}
							}
						}	
						
						if (bImproving) break;
					}
					if (bImproving) break;				
			}
			if (bImproving) break;
		}

		if (bestFit && bestImprovement > 0) {
			// BEST FIT => execute best move if improvement > 0
			localSearchConsolidateBatch(bestWc, bestM, bestB, bestM, bestB+1, bestOp);		

			bImprovedOnce = true;
			bImproving = true;
		}
	}
	return bImprovedOnce;
}

double Schedule::locSearchEvaluateJobSwap(size_t idxFirst, size_t idxSecond, bool& feasible) {
	Job* job1 = scheduledJobs[idxFirst].get();
	Job* job2 = scheduledJobs[idxSecond].get();

	if (job1->getF() != job2->getF()) {
		feasible = false;
		return -1.0;
	}

	if (job1->getStart() < job2->getR()
		|| job2->getStart() < job1->getR()) {
		feasible = false;
		return -1.0;
	}

	feasible = true;
	
	double grossBetter = (max(job1->getC() - job1->getD(), 0) * job1->getW())
		+ (max(job2->getC() - job2->getD(), 0) * job2->getW());

	double grossWorse = (max(job1->getC() - job2->getD(), 0) * job2->getC())
		+ (max(job2->getC() - job1->getD(), 0) * job1->getW());

	
	return grossBetter - grossWorse;
}

pair<double, double> Schedule::locSearchEvaluateJobLeftShift(size_t idxFirst, vector<vector<pair<double, double>>>& possibleLeftShift) {
	Job* job = scheduledJobs[idxFirst].get();
	possibleLeftShift = vector<vector<pair<double, double>>>(size());		// [stage][shift option] true => continuous value (not depending on batching)

	// 1) at each stage check the potential for left shifting disregarding maximal time lags (time constraints) and previous stages
	double earliest = job->getR();
	for (size_t o = 0; o < size(); ++o) {
		possibleLeftShift[o] = this->getLeftShiftOptions(job->getOpPtr(o));
	}

	

	// 2) remove mutually exclusive options 	
	vector<double> currentLeeway = getLeeway(job);	// as of now, how much time between end and start of succeeding operations? => if job[i] is left shifted by t, then job[i-1] must be left shifted by at least max((t - leeway[i]), 0)
	vector<vector<pair<size_t, double>>> currentTcSlack = getTcSlack(job);
	bool possibleOverlaps = false;
	bool possibleTcViolations = false;
	do {
		possibleOverlaps = constrainLeftShiftOptionsFromOverlaps(possibleLeftShift, currentLeeway);
		if (possibleLeftShift[0].empty()) {
			return make_pair(-1, 0);
		}
		possibleTcViolations = constrainLeftShiftOptionsFromTimeConstraints(possibleLeftShift, currentTcSlack);
		if (possibleLeftShift[0].empty()) {
			return make_pair(-1, 0);
		}
	} while (possibleOverlaps || possibleTcViolations);


	// 4) evaluate best options
	pair<double, double> evaluation = make_pair(possibleLeftShift[possibleLeftShift.size() - 1][0].first, 0);
	for (size_t i = 0; i < (possibleLeftShift.size() - 1); ++i) {
		evaluation.second += possibleLeftShift[i][0].first;
	}

	return evaluation; // first = actual left shift of job (=last operation), second = secondary left shifts of intermediate operations
}
pair<double, double> Schedule::locSearchEvaluateJobRightShift(size_t idxJob, size_t idxStg, double time, std::vector<std::vector<std::pair<double, double>>>& possibleRightShifts)
{
	cout << "Schedule::locSearchEvaluateJobRightShift(...) not yet implemented." << endl;
	Job* job = scheduledJobs[idxJob].get();
	possibleRightShifts = vector<vector<pair<double, double>>>(size());
	// 1) at each stage check the possible right shifting disregarding maximal time lags (time constraints) and overlaps
	for (size_t o = 0; o < job->size(); ++o) {
		double delay = 0;
		if (o == idxStg) {
			delay = time;
		}
		possibleRightShifts[o] = getRightShiftOptions(&(*job)[o], delay);
	}

	// 2) remove mutually exclusive options
	vector<double> currentLeeway = getLeeway(job);
	vector<vector<pair<size_t, double>>> currentTcSlack = getTcSlack(job);
	bool possibleOverlaps = false;
	bool possibleTcViolations = false;

	// TODO eine Schleife statt zweier

	do {
		possibleOverlaps = constrainRightShiftOptionsFromOverlaps(possibleRightShifts, currentLeeway);
		//possibleTcViolations = constrainRightShiftOptionsFromTimeConstraints(possibleRightShifts, currentTcSlack);
	} while (possibleOverlaps);
		
	do {
		possibleTcViolations = constrainRightShiftOptionsFromTimeConstraints(possibleRightShifts, currentTcSlack);
	} while(possibleTcViolations);

	// 3) evaluate best option (== smallest right shift)
	
	return std::pair<double, double>();
}


pair<double, double> Schedule::locSearchEvaluateBatchLeftShift(Batch* batch, double time, vector<pair<double, double>>& possibleLeftShifts) {
	cout << "Schedule::locSearchEvaluateBatchLeftShift not yet implemented." << endl;
	pair<double, double> evaluation = make_pair(0, 0);
	possibleLeftShifts = vector<pair<double, double>>(size());
	
	// 1) check the potential for left shifting disregarding maximal time lags (time constraints) and previous stages
	possibleLeftShifts = getLeftShiftOptions(batch);
	
	// TODO


	return evaluation;
}
pair<double, double> Schedule::locSearchEvaluateBatchRightShift(Batch* batch, double time, bool& possible) {
	cout << "Schedule::locSearchEvaluateBatchRightShift not yet implemented." << endl;
	pair<double, double> evaluation = make_pair(0, 0);
	for (size_t i = 0; i < batch->size(); ++i) {
		Operation* op = &(*batch)[i];
		
	
	}
	
	return evaluation;
}

double Schedule::locSearchEvaluateOpConsolidation(size_t idxJob, size_t idxStg, bool& feasible) {
	cout << "Schedule::locSearchEvaluateOpConsolidation() not yet implemented." << endl;
	double evaluation = 0.0;
	Job* job = scheduledJobs[idxJob].get();
	Operation* myOp = &(*job)[idxStg];
	Workcenter* myWc = workcenters[idxStg].get();
	if (myWc->getCap() > 2) {
		Batch* myBatch = myOp->getBatch();
		Machine* myMachine = myBatch->getMachine();
		size_t myBatchIdx = myBatch->getIdx();
		int requiredCap = myBatch->getCap() - myBatch->getAvailableCap();

		// consider consolidation of ops at the same machine
		if (myBatchIdx < myMachine->size() - 1) {
			Batch* succBatch = &(*myMachine)[myBatchIdx + 1];
			if (succBatch->getF() == myOp->getF() && succBatch->getAvailableCap() >= requiredCap) {
				bool bMovePossible = false;
				double succBatchR = succBatch->getR();
				double succBatchW = succBatch->getW();
				double delayOfMyOp = max(0, succBatchR - myOp->getStart());
				if (delayOfMyOp > 0) {
					// CONSOLIDATION POSSIBLE (LEAVE CASE delay <= 0 TO LOCAL SEARCH OPERATION LEFT SHIFT)
					double leftShift = succBatch->getStart() - succBatchR;
					pair<double, double> addedTime = locSearchEvaluateBatchRightShift(myBatch, delayOfMyOp, bMovePossible);
					if (bMovePossible) {
						// TODO WORK IN PROGRESS
						//pair<double, double> reducedTime = locSearchEvaluateBatchLeftShift(succBatch, leftShift, bMovePossible);
					}
	
				}
			}
		}

		// TODO consolidate ops at different machines
	}


	// new plan 15.12.2025

	


	return evaluation;
}
bool Schedule::locSearchSwapJobs(size_t idxFirst, size_t idxSecond) {
	Job* job1 = scheduledJobs[idxFirst].get();
	Job* job2 = scheduledJobs[idxSecond].get();
	size_t mIdx1 = 0;
	size_t bIdx1 = 0;
	size_t jIdx1 = 0;
	size_t mIdx2 = 0;
	size_t bIdx2 = 0;
	size_t jIdx2 = 0;
	for (size_t o = 0; o < workcenters.size(); ++o) {
		Workcenter* wc = workcenters[o].get();
		Operation* op1 = job1->getOpPtr(o);
		Operation* op2 = job2->getOpPtr(o);
		if (!wc->locateOp(op1, mIdx1, bIdx1, jIdx1) || !wc->locateOp(op2, mIdx2, bIdx2, jIdx2)) {
			return false;
		}

		wc->swapOps(mIdx1, bIdx1, jIdx1, mIdx2, bIdx2, jIdx2);
	}
}
bool Schedule::locSearchLeftShiftJob(size_t jobIdx, vector<vector<pair<double, double>>>& options) {
	Job* job = scheduledJobs[jobIdx].get();
	for (size_t o = 0; o < job->size(); ++o) {
		executeLeftShiftOption(jobIdx, o, options[o][0]);	// first option must be largest shift
	}
	return true;
}

vector<pair<double, double>> Schedule::getLeftShiftOptions(Operation* op) {
	vector<pair<double, double>> options = vector<pair<double, double>>();
	options.push_back(make_pair(0.0, 0.0));
	
	double earliest = op->getRconsideringRawP();

	int stageIdx = op->getStg() - 1;
	int nMachines = this->getWorkcenters()[stageIdx]->size();
	for (size_t i = 0; i < nMachines; ++i) {
		int nBatches = (*this->getWorkcenters()[stageIdx])[i].size();
			// +++ options in existing batches (DISCRETE => option.first == option.second) +++
			for (size_t j = 0; j < nBatches; ++j) {
			double currentBatchStart = (*this->getWorkcenters()[stageIdx])[i][j].getStart();

			if (currentBatchStart < earliest - TCB::precision) {
				continue;	// op cannot go into this batch, continue with next batches
			}

			if (currentBatchStart >= op->getStart()) {
				break;		// no favorable batches to insert the batch at this machine, continue with next machine
			}

			if((*this->getWorkcenters()[stageIdx])[i][j].getF() == op->getF() && (*this->getWorkcenters()[stageIdx])[i][j].getAvailableCap() >= op->getS()) {
				double tempOption = op->getStart() - (*this->getWorkcenters()[stageIdx])[i][j].getStart();
				options.push_back(make_pair(tempOption, tempOption));
			}
		}

			bool bOnlyOpInBatch = op->getBatch()->size() == 1;

		// +++ free time slots (CONTINUOUS => option.first <= option.second) +++
		if (nBatches == 0) {
			double tempOptionFrom = op->getStart() - earliest;
			double tempOptionTill = 0.0;
			options.push_back(make_pair(tempOptionFrom, tempOptionTill));
		}
		
		
		for (size_t j = 0; j < nBatches; ++j) {
			double currentBatchStart = (*this->getWorkcenters()[stageIdx])[i][j].getStart();
			if (j == 0) {
				// before 1st batch
				if (earliest + op->getP() <= currentBatchStart || (bOnlyOpInBatch && (*this->getWorkcenters()[stageIdx])[i][j].contains(op))) {
					double tempOptionFrom = op->getStart() - earliest;
					double tempOptionTill = op->getStart() - (currentBatchStart - op->getP());
					if ((bOnlyOpInBatch && (*this->getWorkcenters()[stageIdx])[i][j].contains(op))) {
						tempOptionTill = op->getStart() - currentBatchStart;
					}
					options.push_back(make_pair(tempOptionFrom, tempOptionTill));
				}
			}
			else {
				// inbetween batches
				double gapFrom = max((*this->getWorkcenters()[stageIdx])[i][j - 1].getC(), earliest);
				if (gapFrom + op->getP() <= currentBatchStart || (bOnlyOpInBatch && (*this->getWorkcenters()[stageIdx])[i][j].contains(op))) {
					double tempOptionFrom = op->getStart() - gapFrom;
					double tempOptionTill = op->getStart() - (currentBatchStart - op->getP());
					if ((bOnlyOpInBatch && (*this->getWorkcenters()[stageIdx])[i][j].contains(op))) {
						tempOptionTill = op->getStart() - currentBatchStart;
					}
					if (tempOptionFrom > 0.0) {
						options.push_back(make_pair(tempOptionFrom, tempOptionTill));
					}
				}
			}
		}
		// after last batch
		if (nBatches > 0) {
			if ((*this->getWorkcenters()[stageIdx])[i][nBatches - 1].getC() < op->getStart()) {
				double gapFrom = max(earliest, (*this->getWorkcenters()[stageIdx])[i][nBatches - 1].getC());
				double tempOptionFrom = op->getStart() - gapFrom;
				double tempOptionTill = 0.0;
				options.push_back(make_pair(tempOptionFrom, tempOptionTill));
			}
		}
	}

	// NON-INCREASING ORDER
	sort(options.begin(), options.end(), [&](pair<double, double> a, pair<double, double> b) {
		return a.first > b.first;
		});
	
	return options;
}
vector<pair<double, double>> Schedule::getLeftShiftOptions(Batch* batch) {
	vector<pair<double, double>> options = vector<pair<double, double>>();
	options.push_back(make_pair(0.0, 0.0));
	if (batch->size() == 0) {
		return options;
	}

	double earliest = batch->getRconsideringRawP();
	size_t myBatIdx = batch->getIdx();
	size_t myMacIdx = batch->getMachine()->getIdx();
	size_t myStgIdx = (*batch)[0].getStg() - 1;
	int nMachines = this->getWorkcenters()[myStgIdx]->size();
	for (size_t i = 0; i < nMachines; ++i) {
		int nBatches = (*workcenters[myStgIdx])[i].size();
		if (nBatches == 0) {
			double tempOptionFrom = batch->getStart() - earliest;
			double tempOptionTill = 0.0;
			options.push_back(make_pair(tempOptionFrom, tempOptionTill));
		}

		for (size_t j = 0; j < nBatches; ++j) {
			if (i != myMacIdx || j != myBatIdx) {
				double currentBatchStart = (*this->getWorkcenters()[myStgIdx])[i][j].getStart();
				if (j == 0) {
					// before 1st batch
					if (earliest + batch->getP() <= currentBatchStart) {
						double tempOptionFrom = batch->getStart() - earliest;
						double tempOptionTill = batch->getStart() - (currentBatchStart - batch->getP());
						options.push_back(make_pair(tempOptionFrom, tempOptionTill));
					}
				}
				else {
					// inbetween batches
					double gapFrom = max((*this->getWorkcenters()[myStgIdx])[i][j - 1].getC(), earliest);
					if (gapFrom + batch->getP() <= currentBatchStart) {
						double tempOptionFrom = batch->getStart() - gapFrom;
						double tempOptionTill = batch->getStart() - (currentBatchStart - batch->getP());
						if (tempOptionFrom > 0.0) {
							options.push_back(make_pair(tempOptionFrom, tempOptionTill));
						}
					}
				}
			}
		}

		// after last batch
		if (nBatches > 0) {
			if ((*this->getWorkcenters()[myStgIdx])[i][nBatches - 1].getC() < batch->getStart()) {
				double gapFrom = max(earliest, (*this->getWorkcenters()[myStgIdx])[i][nBatches - 1].getC());
				double tempOptionFrom = batch->getStart() - gapFrom;
				double tempOptionTill = 0.0;
				options.push_back(make_pair(tempOptionFrom, tempOptionTill));
			}
		}
	}
	return options;
}
vector<pair<double, double>> Schedule::getRightShiftOptions(Operation* op, double minDelay) {
	double lastStartAddendum = 10;
	vector<pair<double, double>> options = vector<pair<double, double>>();
	if (minDelay == 0) {
		options.push_back(make_pair(0.0, 0.0));
	}
	double newStart = op->getStart() + minDelay;
	int stageIdx = op->getStg() - 1;
	int nMachines = this->getWorkcenters()[stageIdx]->size();
	for (size_t i = 0; i < nMachines; ++i) {
		int nBatches = (*this->getWorkcenters()[stageIdx])[i].size();
		// +++ options in existing batches (DISCRETE => option.first == option.second) +++
		for (size_t j = 0; j < nBatches; ++j) {
			double currentBatchStart = (*this->getWorkcenters()[stageIdx])[i][j].getStart();

			if (currentBatchStart < newStart) {
				continue;	
			}

			if ((*this->getWorkcenters()[stageIdx])[i][j].getF() == op->getF() && (*this->getWorkcenters()[stageIdx])[i][j].getAvailableCap() >= op->getS()) {
				double tempOption = (*this->getWorkcenters()[stageIdx])[i][j].getStart() - op->getStart();
				options.push_back(make_pair(tempOption, tempOption));
			}
		}

		bool bOnlyOpInBatch = op->getBatch()->size() == 1;

		// +++ free time slots (CONTINUOUS => option.first <= option.second) +++
		if (nBatches == 0) {
			double tempOptionFrom = minDelay;
			double tempOptionTill = minDelay + (lastStartAddendum * op->getP());	// to infinity and beyond!
			options.push_back(make_pair(tempOptionFrom, tempOptionTill));
		}

		for (size_t j = 0; j < nBatches; ++j) {
			double currentBatchStart = (*this->getWorkcenters()[stageIdx])[i][j].getStart();
			if (j == 0) {
				// before 1st batch
				if (newStart + op->getP() <= currentBatchStart) {
					double tempOptionFrom = minDelay;
					double tempOptionTill = (currentBatchStart - op->getP()) - op->getStart();
					options.push_back(make_pair(tempOptionFrom, tempOptionTill));
				}
			} else {
				// inbetween batches
				double gapFrom = max((*this->getWorkcenters()[stageIdx])[i][j].getC(), newStart);
				if (gapFrom + op->getP() <= currentBatchStart) {
					double tempOptionFrom = gapFrom - op->getStart();
					double tempOptionTill = (currentBatchStart - op->getP()) - op->getStart();
					options.push_back(make_pair(tempOptionFrom, tempOptionTill));
				}
 			}
		}
		
		// after last batch
		if (nBatches > 0) {
			double tempOptionFrom = max((*this->getWorkcenters()[stageIdx])[i][nBatches-1].getC(), newStart) - op->getStart();
			double tempOptionTill = tempOptionFrom + (lastStartAddendum * op->getP());	// to infinity and beyond!
			options.push_back(make_pair(tempOptionFrom, tempOptionTill));
		}
	}
	
	// NON-DECREASING ORDER
	sort(options.begin(), options.end(), [&](pair<double, double> a, pair<double, double> b) {
		return a.first < b.first;
		});
	return options;
}
vector<double> Schedule::getLeeway(Job* job) {
	vector<double> leeway = vector<double>(job->size());
	leeway[0] = (*job)[0].getStart() - job->getR();
	for (size_t i = 1; i < job->size(); ++i) {
		leeway[i] = (*job)[i].getStart() - (*job)[i - 1].getC();
	}
	return leeway;
}
vector<double> Schedule::getRightSideLeeway(Job* job) {
	vector<double> rsLeeway = vector<double>(job->size());
	for (size_t i = 0; i < job->size() - 1; ++i) {
		rsLeeway[i] = (*job)[i + 1].getStart() - (*job)[i].getC();
	}
	rsLeeway[job->size() - 1] = 100 * (*job)[job->size() - 1].getP(); // infinity!
	return rsLeeway;
}
vector<vector<pair<size_t, double>>> Schedule::getTcSlack(Job* job) {
	vector<vector<pair<size_t, double>>> tcSlack = vector<vector<pair<size_t, double>>>(job->size());
	for (size_t i = 0; i < job->size(); ++i) {
		vector<pair<int, double>> tcFwd = job->getTcMaxFwd(i);
		tcSlack[i] = vector<pair<size_t, double>>();
		for (size_t j = 0; j < tcFwd.size(); ++j) {
			if (tcFwd[j].first > i && tcFwd[j].second < 999999) {
				double tempTcSlack = ((*job)[i].getStart() + tcFwd[j].second) - (*job)[tcFwd[j].first].getStart();
				tcSlack[i].push_back(make_pair(tcFwd[j].first, tempTcSlack));
			}
		}
	}
	return tcSlack;
}
bool Schedule::constrainLeftShiftOptionsFromOverlaps(std::vector<std::vector<std::pair<double, double>>>& options, std::vector<double>& leeway) {
	bool changeApplied = false;
	for (int o = options.size() - 1; o >= 1; --o) {
		for (int opt = options[o].size() - 1; opt >= 0; --opt) {
			double optFrom = options[o][opt].first;
			double optTill = options[o][opt].second;
						
			// 1) delete discrete options if infeasible
			if (optFrom <= (optTill + TCB::precision) && optFrom >= (optTill - TCB::precision)) {
				bool feasible = false;
				for (size_t prevOpt = 0; prevOpt < options[o - 1].size(); ++prevOpt) {
					if (optFrom <= options[o - 1][prevOpt].first + leeway[o]) {
						feasible = true;
						break;
					}
				}
				if (!feasible) {
					changeApplied = true;
					options[o].erase(options[o].begin() + opt);
					if (options[o].empty()) {	// if not even the zero shift option is available for an operation, then no shifts are feasible
						for (size_t i = 0; i < options.size(); ++i) {
							options[i].clear();
						}
						return true;
					}
				}
			} else if (optFrom > optTill - TCB::precision) {
				// 2) reduce continuous options if necessary
				double maxLeftShift = 0.0;
				for (size_t prevOpt = 0; prevOpt < options[o - 1].size(); ++prevOpt) {
					maxLeftShift = max(maxLeftShift, options[o - 1][prevOpt].first + leeway[o]);
				}
				
				if (maxLeftShift < options[o][opt].first) {
					changeApplied = true;
					if (maxLeftShift < options[o][opt].second) {
						options[o].erase(options[o].begin() + opt);
						if (options[o].empty()) {	// if not even the zero shift option is available for an operation, then no shifts are feasible
							for (size_t i = 0; i < options.size(); ++i) {
								options[i].clear();
							}
							return true;
						}
					} else {
						options[o][opt].first = maxLeftShift;
					}
				}
			}	
		}
	}
	return changeApplied;
}
bool Schedule::constrainLeftShiftOptionsFromTimeConstraints(std::vector<std::vector<std::pair<double, double>>>& options, vector<vector<pair<size_t, double>>>& tcSlack) {
	bool changeApplied = false;
	for (int o = 0; o < options.size() - 1; ++o) {
		for (int opt = 0; opt < options[o].size(); ++opt) {
			double optFrom = options[o][opt].first;
			double optTill = options[o][opt].second;

			// 1) reduce discrete options if infeasible
			if (optFrom <= optTill + TCB::precision && optFrom >= optTill - TCB::precision) {
				bool allFeasible = true;
				bool tempFeasible = false;
				if (tcSlack[o].size() < 1) {
					tempFeasible = true;
				}
				for (size_t tc = 0; tc < tcSlack[o].size(); ++tc) {
					tempFeasible = false;
					size_t succO = tcSlack[o][tc].first;	// successive stage connected by a time constraint
					for (size_t succOpt = 0; succOpt < options[succO].size(); ++succOpt) {
						if (options[o][opt].first <= options[succO][succOpt].first + tcSlack[o][tc].second + TCB::precision) {
							tempFeasible = true;
							break;
						}
					}
					if (!tempFeasible) {
						allFeasible = false;
					} 
				}
				
				if (!allFeasible) {
					changeApplied = true;
					options[o].erase(options[o].begin() + opt);
					if (options[o].empty()) {	// if not even the zero shift option is available for an operation, then no shifts are feasible
						for (size_t i = 0; i < options.size(); ++i) {
							options[i].clear();
						}
						return true;
					}
				}
			}

			// 2) reduce continuous options if necessary
			if (optFrom > optTill) {	
				double minMaxLeftShift = DBL_MAX;	// the minimum of the maximal feasible left shifts of all stages	[JR-2026-Jan-10]
				for (size_t tc = 0; tc < tcSlack[o].size(); ++tc) {
					double maxLeftShift = 0;			// the maximal feasible left shift bound by one successing stage [JR-2026-Jan-10] moved this down into the for(tc<tcSlack[o])-loop
					size_t succO = tcSlack[o][tc].first;
					for (size_t succOpt = 0; succOpt < options[succO].size(); ++succOpt) {
						maxLeftShift = max(maxLeftShift, options[succO][succOpt].first + tcSlack[o][tc].second);	// the max possible option remains
					}

					if (maxLeftShift < minMaxLeftShift) {				// [JR-2026-Jan-10] consider the minimum of the maximal feasible shifts
						minMaxLeftShift = maxLeftShift;
					}
				}
				if (minMaxLeftShift < options[o][opt].first) {			// [JR-2026-Jan-10] moved this block out of the for(tc < tcSlack[o])-block
					changeApplied = true;
					if (minMaxLeftShift < options[o][opt].second) {		// [JR-2026-Jan-10] changed maxLeftShift to minMaxLeftShift
						options[o].erase(options[o].begin() + opt);
						if (options[o].empty()) {	// if not even the zero shift option is available for an operation, then no shifts are feasible
							for (size_t i = 0; i < options.size(); ++i) {
								options[i].clear();
							}
							return true;
						}
					} else {
						options[o][opt].first = minMaxLeftShift;		// [JR-2026-Jan-10] changed maxLeftShift to minMaxLeftShift
					}
				}
			}

		}
	}
	return changeApplied;
}

bool Schedule::constrainRightShiftOptionsFromOverlaps(std::vector<std::vector<std::pair<double, double>>>& options, std::vector<double>& leeway) {
	bool changeApplied = false;
	for (int o = options.size() - 1; o >= 1; --o) {
		for (int opt = options[o].size() - 1; opt >= 0; --opt) {
			double optFrom = options[o][opt].first;
			double optTill = options[o][opt].second;

			// 1) delete discrete options if infeasible
			if (optFrom == optTill) {
				bool feasible = false;
				for (size_t prevOpt = 0; prevOpt < options[o - 1].size(); ++prevOpt) {
					if (optFrom >= options[o - 1][prevOpt].first + leeway[o]) {
						feasible = true;
						break;
					}
				}
				if (!feasible) {
					changeApplied = true;
					options[o].erase(options[o].begin() + opt);
				}
			}

			// 2) reduce continuous options if necessary
			if (optFrom < optTill) {
				double minRightShift = DBL_MAX;
				for (size_t prevOpt = 0; prevOpt < options[o - 1].size(); ++prevOpt) {
					minRightShift = min(minRightShift, options[o - 1][prevOpt].first + leeway[o]);
				}

				if (minRightShift > options[o][opt].first) {
					changeApplied = true;
					if (minRightShift > options[o][opt].second) {
						options[o].erase(options[o].begin() + opt);
					}
					else {
						options[o][opt].first = minRightShift;
					}
				}
			}
		}
	}
	return changeApplied;
}
bool Schedule::constrainRightShiftOptionsFromTimeConstraints(std::vector<std::vector<std::pair<double, double>>>& options, std::vector<std::vector<std::pair<size_t, double>>>& tcSlack) {
	bool changeApplied = false;
	for (int o = 0; o < options.size() - 1; ++o) {
		for (int opt = 0; opt < options[o].size(); ++opt) {
			double optFrom = options[o][opt].first;
			double optTill = options[o][opt].second;

			// 1) reduce discrete options if infeasible
			if (optFrom == optTill) {
				bool allFeasible = true;
				bool tempFeasible = false;
				if (tcSlack[o].size() < 1) {
					tempFeasible = true;
				}
				for (size_t tc = 0; tc < tcSlack[o].size(); ++tc) {
					tempFeasible = false;
					size_t succO = tcSlack[o][tc].first;	// successive stage connected by a time constraint
					for (size_t succOpt = 0; succOpt < options[succO].size(); ++succOpt) {
						if (options[o][opt].first >= options[succO][succOpt].first - tcSlack[o][tc].second) {
							tempFeasible = true;
							break;
						}
					}
					if (!tempFeasible) {
						allFeasible = false;
					}
				}

				if (!allFeasible) {
					changeApplied = true;
					options[o].erase(options[o].begin() + opt);
				}
			}

			// 2) reduce continuous options if necessary
			if (optFrom < optTill) {
				double minRightShift = DBL_MAX;
				for (size_t tc = 0; tc < tcSlack[o].size(); ++tc) {
					size_t succO = tcSlack[o][tc].first;
					for (size_t succOpt = 0; succOpt < options[succO].size(); ++succOpt) {
						minRightShift = min(minRightShift, options[succO][succOpt].first - tcSlack[o][tc].second);
					}
					if (minRightShift > options[o][opt].first) {
						changeApplied = true;
						if (minRightShift > options[o][opt].second) {
							options[o].erase(options[o].begin() + opt);
						}
						else {
							options[o][opt].first = minRightShift;
						}
					}
				}
			}

		}
	}
	return changeApplied;
}



void Schedule::executeLeftShiftOption(size_t jobIdx, size_t stgIdx, std::pair<double, double>& option) {
	if (option.first > 0) {
		Operation* operation = &(*scheduledJobs[jobIdx])[stgIdx];
		bool bIntoExistingBatch = option.first == option.second;
		double newStart = operation->getStart() - option.first;

		if (!workcenters[stgIdx]->moveOpDisregardingTc(operation, newStart, bIntoExistingBatch)) {
			TCB::logger.Log(Error, "Exception thrown in Schedule::executeLeftShiftOption(...)");
			throw(ExcSched("ERROR: Schedule::executeLeftShiftOption(...) invalid."));
		}
	}
}
bool Schedule::executeLeftShiftOption(Operation* operation, std::pair<double, double>& option) {
	bool bIntoExistingBatch = true;
	if (option.first > 0) {
		bIntoExistingBatch = option.first == option.second;
		double newStart = operation->getStart() - option.first;
		if (!workcenters[operation->getStg()-1]->moveOpDisregardingTc(operation, newStart, bIntoExistingBatch)) {
			TCB::logger.Log(Error, "Exception thrown in Schedule::executeLeftShiftOption(...)");
			throw(ExcSched("ERROR: Schedule::executeLeftShiftOption(...) invalid."));
		}
	}
	return bIntoExistingBatch;
}

bool Schedule::isValid() const {
	// ALL OPERATIONS OF ALL JOBS ARE ASSIGNED
	for (size_t j = 0; j < problem->getN(); ++j) {
		for (size_t o = 0; o < (*problem)[j].size(); ++o) {
			if (!contains(&(*problem)[j][o])) {
				//TCB::logger.Log(Error, "missing operation " + to_string((*problem)[j][o].getId()) + "." + to_string((*problem)[j][o].getStg()));
				return false;
			}
		}
	}

	// NO OVERLAPPING PROCESSING ON ANY MACHINE
	for (size_t wc = 0; wc < size(); ++wc) {
		for (size_t m = 0; m < (*workcenters[wc]).size(); ++m) {
			if ((*workcenters[wc])[m].hasOverlaps()) {
				//TCB::logger.Log(Error, "overlapping processing of batches at machine " + to_string(workcenters[wc]->getId()) + "." + to_string((*workcenters[wc])[m].getId()));
				return false;
			}
		}
	}

	// NO OPERATION IS STARTED BEFORE ITS ROUTE PREDECESSOR IS COMPLETED + TIME CONSTRAINTS ARE MET
	for (size_t wc = 0; wc < size(); ++wc) {
		for (size_t m = 0; m < (*workcenters[wc]).size(); ++m) {
			for (size_t b = 0; b < (*workcenters[wc])[m].size(); ++b) {
				for (size_t o = 0; o < (*workcenters[wc])[m][b].size(); ++o) {
					if (!(*workcenters[wc])[m][b][o].checkProcessingOrder()) {
						//TCB::logger.Log(Error, "Processing order violated");
						return false;
					}
					if (!(*workcenters[wc])[m][b][o].checkTimeConstraints()) {
						//TCB::logger.Log(Error, "Time constraint violated");
						return false;
					}
				}
			}
		}
	}
	return true;
}

double Schedule::getTWT() const {
	double twt = 0;
	for (size_t wc = 0; wc < size(); ++wc) {
		twt += workcenters[wc]->getTWT();
	}
	return twt;
}

double Schedule::getMinMSP(size_t stgIdx) const {
	if (stgIdx >= size()) throw out_of_range("Schedule::getMSP() out of range");
	return workcenters[stgIdx]->getMinMSP();
}

void Schedule::saveJson(std::string solver) {
	pt::ptree treeFile;

	treeFile.put("Problem", TCB::prob->getFilename());
	treeFile.put("Solver", solver);
	stringstream ssTWT;
	ssTWT << setprecision(3) << fixed << getTWT();
	treeFile.put("TWT", ssTWT.str());

	pt::ptree treeSchedule;
	for (size_t wc = 0; wc < workcenters.size(); ++wc) {
		Workcenter* WC = &(*workcenters[wc]);
		pt::ptree treeWorkcenter;
		treeWorkcenter.put("stage", WC->getId());
		for (size_t m = 0; m < WC->size(); ++m) {
			Machine* MAC = &(*WC)[m];
			pt::ptree treeMachine;
			treeMachine.put("id", MAC->getId());
			treeMachine.put("capacity", MAC->getCap());
			for (size_t b = 0; b < MAC->size(); ++b) {
				Batch* BAT = &(*MAC)[b];
				pt::ptree treeBatch;
				stringstream ssStart;
				stringstream ssCompletion;
				ssStart << setprecision(3) << fixed << BAT->getStart();
				ssCompletion << setprecision(3) << fixed << BAT->getC();
				treeBatch.put("start", ssStart.str());
				treeBatch.put("completion", ssCompletion.str());
				for (size_t j = 0; j < BAT->size(); ++j) {
					Operation* OP = &(*BAT)[j];
					pt::ptree treeOp;
					treeOp.put("id", OP->getId());
					treeOp.put("stage", OP->getStg());
					treeBatch.add_child("operations.Operation", treeOp);
				}
				treeMachine.add_child("batches.Batch", treeBatch);
			}
			treeWorkcenter.add_child("machines.Machine", treeMachine);
		}
		treeSchedule.add_child("workcenters.Workcenter", treeWorkcenter);
	}
	treeFile.add_child("Schedule", treeSchedule);

	string probFileName = extractFileName(TCB::prob->getFilename());
	string probFileNameWithoutExtension = probFileName.substr(0, probFileName.find(".dat"));
	string solverName = solver;
	replaceWindowsSpecialCharsWithUnderscore(solverName);

	string sSeed = to_string(TCB::seed);
	string schedFileName = probFileNameWithoutExtension + "_" + solverName + sSeed + ".json";

	bool success = CreateDirectory(L".\\results", NULL);
	string pathAndFilename = string(".\\results\\").append(schedFileName);
	pt::write_json(pathAndFilename, treeFile);
}
void Schedule::saveJsonFactory(std::string solver) const {
	pt::ptree treeFile;

	treeFile.put("Problem", TCB::prob->getFilename());
	treeFile.put("Solver", solver);
	stringstream ssTWT;
	ssTWT << setprecision(3) << fixed << getTWT();
	treeFile.put("TWT", ssTWT.str());

	pt::ptree treeFactory;
	pt::ptree arrayWorkareas;
	pt::ptree treeWorkarea;
	treeWorkarea.put("name", "Flow Shop");
	pt::ptree arrayWorkcenters;

	for (size_t wc = 0; wc < workcenters.size(); ++wc) {
		Workcenter* WC = &(*workcenters[wc]);
		pt::ptree treeWorkcenter;
		pt::ptree arrayResources;
		treeWorkcenter.put("name", "Stage " + to_string(WC->getId()));

		for (size_t m = 0; m < WC->size(); ++m) {
			Machine* MAC = &(*WC)[m];
			pt::ptree treeMachine;
			pt::ptree arrayLoad;
			treeMachine.put("name", "M" + to_string(MAC->getId()));
			
			for (size_t b = 0; b < MAC->size(); ++b) {
				Batch* BAT = &(*MAC)[b];
				pt::ptree treeBatch;
				pt::ptree arrayContent;
				stringstream ssStart;
				stringstream ssCompletion;
				stringstream ssP;
				ssStart << setprecision(3) << fixed << BAT->getStart();
				ssCompletion << setprecision(3) << fixed << BAT->getC();
				ssP << setprecision(3) << fixed << BAT->getP();
				treeBatch.put("id", "B" + to_string((b + 1)));
				treeBatch.put("type", "Batch");
				treeBatch.put("start", ssStart.str());
				treeBatch.put("C", ssCompletion.str());
				treeBatch.put("f", BAT->getF());
				treeBatch.put("p", ssP.str());
				treeBatch.put("capacity", BAT->getCap());
				for (size_t j = 0; j < BAT->size(); ++j) {
					Operation* OP = &(*BAT)[j];
					stringstream ssD;
					stringstream ssR;
					stringstream ssW;	
					ssD << setprecision(3) << fixed << OP->getD();
					ssR << setprecision(3) << fixed << OP->getR();
					ssW << setprecision(3) << fixed << OP->getW();
					pt::ptree treeOp;
					treeOp.put("id", to_string(OP->getId()) + "." + to_string(OP->getStg()));
					treeOp.put("type", "Operation");
					treeOp.put("d", ssD.str());
					treeOp.put("r", ssR.str());
					treeOp.put("p", to_string(OP->getP()));
					treeOp.put("f", to_string(OP->getF()));
					treeOp.put("s", to_string(OP->getS()));
					treeOp.put("w", ssW.str());
					arrayContent.push_back(make_pair("", treeOp));
				}
				treeBatch.add_child("content", arrayContent);
				arrayLoad.push_back(make_pair("", treeBatch));
			}
			treeMachine.add_child("load", arrayLoad);
			arrayResources.push_back(std::make_pair("", treeMachine));
		}
		treeWorkcenter.add_child("resources", arrayResources);
		arrayWorkcenters.push_back(make_pair("", treeWorkcenter));
		
	}
	treeWorkarea.add_child("workcenters", arrayWorkcenters);
	arrayWorkareas.push_back(make_pair("", treeWorkarea));
	treeFactory.add_child("workareas", arrayWorkareas);
	treeFile.add_child("Factory", treeFactory);

	string probFileName = extractFileName(TCB::prob->getFilename());
	string probFileNameWithoutExtension = probFileName.substr(0, probFileName.find(".dat"));
	string solverName = solver;
	replaceWindowsSpecialCharsWithUnderscore(solverName);

	string sSeed = to_string(TCB::seed);
	string schedFileName = probFileNameWithoutExtension + "_" + solverName + sSeed + ".json";

	bool success = CreateDirectory(L".\\results", NULL);
	string pathAndFilename = string(".\\results\\").append(schedFileName);
	pt::write_json(pathAndFilename, treeFile);
}

void Schedule::debugSetR(size_t scheduledJobIdx, double newR) {
	cout << "+++ WARNING: INFEASIBLE OPERATION EXECUTED IN COURSE OF DEBUGGING +++" << endl;
	scheduledJobs[scheduledJobIdx]->setR(newR);
}
void Schedule::debugAddMachine(size_t stgIdx) {
	unique_ptr<Machine> newMac = make_unique<Machine>((*workcenters[stgIdx])[workcenters[stgIdx]->size() - 1].getId() + 1, workcenters[stgIdx]->getCap(), &(*workcenters[stgIdx]));
	workcenters[stgIdx]->addMachine(move(newMac));
}

void Schedule::debugAllBatchesNotEmptyAndWithMachineReference() {
	for (size_t wcIdx = 0; wcIdx < size(); ++wcIdx) {
		Workcenter* wc = workcenters[wcIdx].get();
		for (size_t mIdx = 0; mIdx < wc->size(); ++mIdx) {
			Machine* mac = &(*wc)[mIdx];
			for (size_t bIdx = 0; bIdx < mac->size(); ++bIdx) {
				if ((*mac)[bIdx].size() <= 0) {
					cout << "Batch at stage " << (wcIdx + 1) << ", M" << (mIdx + 1) << " is empty!" << endl;
					throw ExcSched("EMPTY BATCH");
				}

				if ((*mac)[bIdx].getMachine() == nullptr) {					
					cout << "Stage " << (wcIdx + 1) << ", M" << (mIdx + 1) << ", Batch " << (bIdx + 1) << ": no machine reference!" << endl;
					throw ExcSched("MISSING MACHINE REFERENCE");
				}
			}
		}
	}
}

