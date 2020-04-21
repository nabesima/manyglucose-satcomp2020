/***********************************************************************************[LocalVals.h]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson
Copyright (c) 2017-2018, Hidetomo Nabeshima

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

#include "PrdClauses.h"
#include "ParallelSolver.h"

#include <assert.h>

using namespace Glucose;

PrdClauses::PrdClauses(int thread_id, int period) :
        thn(thread_id),
        prd(period),
        num_clauses(0),
        num_exported_threads(0),
        completed(false)
{
    pthread_mutex_init(&lock_num_exported_threads, NULL);
    pthread_mutex_init(&lock_completed, NULL);
    pthread_cond_init(&is_completed, NULL);
}

bool PrdClauses::pushClause(const Lit lit) {
    assert(!completed);
    lits.push(1);             // size
    lits.push(toInt(lit));
    lits.push(ParallelSolver::calculateHash(lit)); // added by kanbara.
    num_clauses++;
    return true;
}

bool PrdClauses::pushClause(const Clause& c) {
    assert(!completed);
    lits.push(c.size());      // size
    for (int i=0; i < c.size(); i++)
        lits.push(toInt(c[i]));
    lits.push(c.getHash()); // added by kanbara
    num_clauses++;
    return true;
}

// When the period of the thread is finished (it means that the addition of clauses is completed),
// then this method is called by the thread. This method notifies waiting threads to be completed.
void PrdClauses::completeAddition() {
    assert(completed == false);
    pthread_mutex_lock(&lock_completed);
    completed = true;
    // Send signals to threads which are waiting to be completed.
    pthread_cond_broadcast(&is_completed);
    pthread_mutex_unlock(&lock_completed);
}

// Wait the addition of clauses to be completed.
void PrdClauses::waitAdditionCompleted(void) {
    pthread_mutex_lock(&lock_completed);
    if (!completed)
        pthread_cond_wait(&is_completed, &lock_completed);
    pthread_mutex_unlock(&lock_completed);
}

// When exporting to the specified thread is finished, then this method is called.
void PrdClauses::completeExportation(int thread_id) {
    pthread_mutex_lock(&lock_num_exported_threads);
    num_exported_threads++;
    pthread_mutex_unlock(&lock_num_exported_threads);
}

int PrdClauses::getNumExportedThreads(void) {
    pthread_mutex_lock(&lock_num_exported_threads);
    int num = num_exported_threads;
    pthread_mutex_unlock(&lock_num_exported_threads);
    return num;
}

