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

#ifndef Glucose_PrdClauses_h
#define Glucose_PrdClauses_h

#include "../mtl/Vec.h"
#include "../core/SolverTypes.h"

namespace Glucose {

//=================================================================================================
// A set of clauses acquired with a certain period of a thread
class PrdClauses {
private:
    int thn;                     // thread number
    int64_t prd;                 // period number
    int num_clauses;             // the number of clauses
    vec<uint32_t> lits;          // concatenated literals

    int num_exported_threads;    // the number of threads to which these clauses are exported.
    pthread_mutex_t lock_num_exported_threads;    // mutex on the variable "num_exported_threads"

    bool completed;                     // whether the addition of clauses from the thread is finished
    pthread_mutex_t lock_completed;     // mutex on the variable "completed"
    pthread_cond_t  is_completed;       // condition variable that says that this set of clauses is completed

public:
    PrdClauses(int thread_id, int period);

    bool pushClause(const Lit lit);
    bool pushClause(const Clause& c);

    // When the period of the thread is finished (it means that the addition of clauses is completed),
    // then this method is called by the thread. This method notifies waiting threads to be completed.
    void completeAddition();

    // Wait the addition of clauses to be completed.
    void waitAdditionCompleted(void);

    // Methods for exportation
    int size(void) const { return lits.size(); };
    uint32_t operator [] (int index) const { return lits[index]; }

    // When exporting to the specified thread is finished, then this method is called.
    void completeExportation(int thread_id);
    int getNumExportedThreads(void);


    // Misc
    int64_t  period(void)          const { return prd; }
    int      numClauses(void)      const { return num_clauses; }
};

//=================================================================================================

}

#endif
