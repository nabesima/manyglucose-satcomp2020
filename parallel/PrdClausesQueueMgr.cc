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

#include "PrdClausesQueueMgr.h"

#include <assert.h>

using namespace Glucose;

PrdClausesQueueMgr::PrdClausesQueueMgr(int num_threads) {
    setNumThreads(num_threads);
}

void PrdClausesQueueMgr::setNumThreads(int num_threads) {
    assert(queues.size() == 0);
    for (int i=0; i < num_threads; i++)
        queues.push(new PrdClausesQueue(i, num_threads));
}

PrdClausesQueueMgr::~PrdClausesQueueMgr() {
    for (int i=0; i < queues.size(); i++)
        delete queues[i];
    queues.clear();
}

PrdClausesQueue& PrdClausesQueueMgr::get(int thread_id) const {
    assert(0 <= thread_id);
    //printf("thread_id = %d, queues.size = %d\n", thread_id, queues.size());
    assert(thread_id < queues.size());
    return *queues[thread_id];
}
