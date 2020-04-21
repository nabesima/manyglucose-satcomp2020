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

#include "PrdClausesQueue.h"

#include <assert.h>

using namespace Glucose;

PrdClausesQueue::PrdClausesQueue(int thread_id, int nb_threads) :
        thn(thread_id),
        num_threads(nb_threads),
        next_period(nb_threads)
{
    pthread_rwlock_init(&rwlock, NULL);

    // Add an empty set of clauses to which clauses acquired at period 0 are stored.
    queue.insert(new PrdClauses(thn, 0));
}

PrdClausesQueue::~PrdClausesQueue()
{
    while (queue.size() > 0) {
        delete queue.peek();
        queue.pop();
    }
}

// When the current period of the thread is finished, then this method is called by the thread.
// This method notifies waiting threads to be completed.
void PrdClausesQueue::completeAddtion()
{
    assert(queue.size() > 0);
    PrdClauses& last = *queue[queue.size() - 1];

    pthread_rwlock_wrlock(&rwlock);

    // Add an empty set of clauses to which clauses acquired at the next period are stored.
    queue.insert(new PrdClauses(thn, last.period() + 1));

    // DEBUG
    //printf("T%d@p%" PRIu64 ": %dclauses, %dlits\n", thn, last.period(), last.numClauses(), last.size());

    // Complete and notify it to all waiting threads
    last.completeAddition();

    // Remove a set of clauses that were sent to other threads
    while (queue.size() > 1) {
        PrdClauses* head = queue.peek();
        if (head->getNumExportedThreads() != num_threads - 1)
            break;
        queue.pop();
        delete head;
    }

    pthread_rwlock_unlock(&rwlock);
}

// When exporting to the specified thread is finished, then this method is called.
void PrdClausesQueue::completeExportation(int thn, PrdClauses& prdClauses)
{
    assert(next_period[thn] == prdClauses.period());
    next_period[thn]++;
    prdClauses.completeExportation(thn);
}


// Get a set of clauses which are generated at the specified period.
PrdClauses* PrdClausesQueue::get(int thread, int period) {
    int p = next_period[thread];
    if (period < p) return NULL;    // 'period' is already exported to 'thn'

    pthread_rwlock_rdlock(&rwlock);
    assert(queue.size() > 0);

    //printf("queue[0]->period() = %d, period = %d\n", queue[0]->period(), period);
    assert(queue[0]->period() <= p);

    int index = p - queue[0]->period();
    assert(0 <= index);
    assert(index < queue.size());
    PrdClauses *prdClauses = queue[index];
    //printf("> queue[%d]->period() = %d, period = %d\n", index, queue[index]->period(), period);
    assert(prdClauses->period() == p);
    pthread_rwlock_unlock(&rwlock);

    return prdClauses;
}





