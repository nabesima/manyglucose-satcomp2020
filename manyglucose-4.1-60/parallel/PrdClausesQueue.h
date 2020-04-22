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

#ifndef Glucose_PrdClausesQueue_h
#define Glucose_PrdClausesQueue_h

#include "../mtl/Queue.h"
#include "../parallel/PrdClauses.h"

namespace Glucose {

//=================================================================================================
// A set of clauses acquired at a certain thread
class PrdClausesQueue {
private:
    int thn;                       // thread number
    int num_threads;               // the number of threads
    vec<int> next_period;          // the next period for exporting to the specified thread
    Queue<PrdClauses *> queue;     // a list of sets of clauses.
    pthread_rwlock_t    rwlock;    // read/write lock of this object.

public:
    PrdClausesQueue(int thread_id, int nb_threads);
    ~PrdClausesQueue();

    // When the current period of the thread is finished, then this method is called by the thread.
    // This method notifies waiting threads to be completed.
    void completeAddtion(void);

    // When exporting to the specified thread is finished, then this method is called.
    void completeExportation(int thread_id, PrdClauses& prdClauses);

    // Get a set of clauses which are generated at the specified period.
    PrdClauses* get(int thread, int period);

    // Return the last set of clauses
    PrdClauses& last() { assert(queue.size() > 0); return *queue[queue.size() - 1]; }
};
//=================================================================================================

}

#endif
