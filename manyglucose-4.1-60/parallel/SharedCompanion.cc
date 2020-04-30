/***************************************************************************************[SharedCompanion.cc]
 Glucose -- Copyright (c) 2009-2014, Gilles Audemard, Laurent Simon
                                CRIL - Univ. Artois, France
                                LRI  - Univ. Paris Sud, France (2009-2013)
                                Labri - Univ. Bordeaux, France

 Syrup (Glucose Parallel) -- Copyright (c) 2013-2014, Gilles Audemard, Laurent Simon
                                CRIL - Univ. Artois, France
                                Labri - Univ. Bordeaux, France

Glucose sources are based on MiniSat (see below MiniSat copyrights). Permissions and copyrights of
Glucose (sources until 2013, Glucose 3.0, single core) are exactly the same as Minisat on which it
is based on. (see below).

Glucose-Syrup sources are based on another copyright. Permissions and copyrights for the parallel
version of Glucose-Syrup (the "Software") are granted, free of charge, to deal with the Software
without restriction, including the rights to use, copy, modify, merge, publish, distribute,
sublicence, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

- The above and below copyrights notices and this permission notice shall be included in all
copies or substantial portions of the Software;
- The parallel version of Glucose (all files modified since Glucose 3.0 releases, 2013) cannot
be used in any competitive event (sat competitions/evaluations) without the express permission of
the authors (Gilles Audemard / Laurent Simon). This is also the case for any competitive event
using Glucose Parallel as an embedded SAT engine (single core or not).


--------------- Original Minisat Copyrights

Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include "../core/Solver.h"
#include "../parallel/ParallelSolver.h"
#include "../core/SolverTypes.h"
//#include "../parallel/ClausesBuffer.h"
#include "../parallel/SharedCompanion.h"
#include "../parallel/PeriodClausesBuffer.h"

using namespace Glucose;

SharedCompanion::SharedCompanion(int _nbThreads) :
    clausesMgr(_nbThreads),
    nbThreads(_nbThreads),
    numLiveThreads(0),          // added by nabesima
    bjobFinished(false),
    jobFinishedBy(NULL),
    panicMode(false), // The bug in the SAT2014 competition :)
    jobStatus(l_Undef),
    jobFinishedPeriod(0),    // added by gotou
    random_seed(9164825) {

    pthread_mutex_init(&mutexSharedClauseCompanion,NULL); // This is the shared companion lock
    //pthread_mutex_init(&mutexSharedUnitCompanion,NULL); // This is the shared companion lock    modified by nabesima
    //pthread_mutex_init(&mutexSharedCompanion,NULL); // This is the shared companion lock
    pthread_mutex_init(&mutexJobFinished,NULL); // This is the shared companion lock
    pthread_barrier_init(&barrierMargin0AfterImport, NULL, nbThreads); // added by gotou
    if (_nbThreads> 0)  {
        setNbThreads(_nbThreads);
//        fprintf(stdout,"c Shared companion initialized: handling of clauses of %d threads.\nc %d ints for the sharing clause buffer (not expandable) .\n", _nbThreads, clausesBuffer.maxSize());
    }
}

// added by nabesima for incremental SAT solving
void SharedCompanion::init() {
    bjobFinished      = false;
    jobFinishedBy     = NULL;
    panicMode         = false;
    jobStatus         = l_Undef;
    jobFinishedPeriod = 0;
}

void SharedCompanion::setNbThreads(int _nbThreads) {
    nbThreads = _nbThreads;
    clausesMgr.setNumThreads(nbThreads);
    clausesBuffer.setNbThreads(nbThreads);
    pthread_barrier_init(&barrierMargin0AfterImport, NULL, nbThreads); // added by gotou
    // added by nabesima
    last_propagations.growTo(nbThreads);
    pthread_rwlock_init(&last_propagations_rwlock, NULL);
}

void SharedCompanion::printStats() {
}

// No multithread safe
bool SharedCompanion::addSolver(ParallelSolver* s) {
    watchedSolvers.push(s);
    pthread_mutex_t* mu = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(mu,NULL);
    assert(s->thn == watchedSolvers.size()-1); // all solvers must have been registered in the good order
    nextUnit.push(0);

    return true;
}
void SharedCompanion::newVar(bool sign) {
    isUnary .push(l_Undef);
}

void SharedCompanion::addLearnt(ParallelSolver *s, Lit unary) {
    // modified by nabesima
    PrdClauses& prdClauses = clausesMgr.get(s->thn).last();
    assert(prdClauses.period() == s->periods);
    prdClauses.pushClause(unary);
}

Lit SharedCompanion::getUnary(ParallelSolver *s) {
    int sn = s->thn;
    Lit ret = lit_Undef;

    // modified by nabesima
    //pthread_mutex_lock(&mutexSharedUnitCompanion);
    pthread_mutex_lock(&mutexSharedClauseCompanion);
    if (nextUnit[sn] < unitLit.size())
        ret = unitLit[nextUnit[sn]++];
    //pthread_mutex_unlock(&mutexSharedUnitCompanion);
    pthread_mutex_unlock(&mutexSharedClauseCompanion);
    return ret;
}

// Specialized functions for this companion
// must be multithread safe
// Add a clause to the threads-wide clause database (all clauses, through)
bool SharedCompanion::addLearnt(ParallelSolver *s, Clause & c) {
    // modified by nabesima
    PrdClauses& prdClauses = clausesMgr.get(s->thn).last();
    assert(prdClauses.period() == s->periods);
    prdClauses.pushClause(c);
    return true;

//    int sn = s->thn; // thread number of the solver
//    bool ret = false;
//    assert(watchedSolvers.size()>sn);
//
//    pthread_mutex_lock(&mutexSharedClauseCompanion);
//    ret = clausesBuffer.pushClause(sn, c);
//    pthread_mutex_unlock(&mutexSharedClauseCompanion);
//    return ret;
}


bool SharedCompanion::getNewClause(ParallelSolver *s, int & threadOrigin, vec<Lit>& newclause) { // gets a new interesting clause for solver s
    int sn = s->thn;

    // First, let's get the clauses on the big blackboard
    pthread_mutex_lock(&mutexSharedClauseCompanion);
    // modified by gotou
//    bool b = clausesBuffer.getClause(sn, threadOrigin, newclause);
    bool b = clausesBuffer.getClause(sn, s->periods, s->margin, threadOrigin, newclause);
    pthread_mutex_unlock(&mutexSharedClauseCompanion);

    return b;
}

bool SharedCompanion::jobFinished(int64_t periods) {
    bool ret = false;
    pthread_mutex_lock(&mutexJobFinished);
    // modified by nabesima
    //ret = bjobFinished;
    ret = bjobFinished && (jobStatus == l_Undef || jobFinishedPeriod + margin < periods);
    pthread_mutex_unlock(&mutexJobFinished);
    return ret;
}

// modified by nabesima
//bool SharedCompanion::IFinished(ParallelSolver *s) {
bool SharedCompanion::IFinished(ParallelSolver *s, lbool status) {
    bool found_better_one = false;

    pthread_mutex_lock(&mutexJobFinished);

    //     status != l_Undef && jobStatus != l_Undef -> status == jobStatus
    assert(status == l_Undef || jobStatus == l_Undef || status == jobStatus);

    // modified by nabesima
//    if (!bjobFinished) {
    if (!bjobFinished
            || (status != l_Undef && s->periods <  jobFinishedPeriod)
            || (status != l_Undef && s->periods == jobFinishedPeriod  &&  s->thn < jobFinishedBy->thn)
       ) {
        found_better_one  = true;
        bjobFinished      = true;
        jobFinishedBy     = s;
        jobFinishedPeriod = s->periods;         // added by gotou
        jobStatus         = status;             // added by nabesima
    }

    pthread_mutex_unlock(&mutexJobFinished);
    return found_better_one;
}


// added by gotou
bool SharedCompanion::moveToNextPeriod(int threadId, int margin) {
    bool b;
    pthread_mutex_lock(&mutexSharedClauseCompanion);
    b = clausesBuffer.nextPeriodProcess(threadId, margin);
    pthread_mutex_unlock(&mutexSharedClauseCompanion);
    return b;
}
