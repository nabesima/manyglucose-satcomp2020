/**********************************************************************************[ClausesBuffer.cc]
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

/* ClausesBuffer
 *
 * This class is responsible for exchanging clauses between threads.
 * It is based on a fixed-length FIFO array of literals.
 * If the FIFO is full, then old clauses are removed (even if it was not yet sent to all threads)
 *
 * a clause " l1 l2 l3" is pushed in the FIFO with the following 6 unsigned integers
 * 3 nseen origin l1 l2 l3
 * + 3 is the size of the pushed clause
 * + nseen is the number of thread which imported this clause (initialized with nthreads-1)
 *       (when set to 0, the clause is removed from the fifo)
 * + origin is the thread id of the thread which added this clause to the fifo
 * + l1 l2 l3 are the literals of the clause
 *
 * **********************************************************************************************
 * **CAREFUL** This class is not thread-safe. In glucose-syrup, the SharedCompanion is
 * responsible for ensuring atomicity of main functions
 * **********************************************************************************************
 *
 * */

#include "../parallel/PeriodClausesBuffer.h"

//=================================================================================================

using namespace Glucose;


PeriodClausesBuffer::PeriodClausesBuffer(int _nbThreads) : nbThreads(_nbThreads) {
    lastOfThread.growTo(_nbThreads);
    for (int i = 0; i < nbThreads; i++) lastOfThread[i] = {i?0:1, 0};
    elems.growTo(nbThreads);
    elems_lim.growTo(nbThreads);
    for (int i = 0; i < nbThreads; i++) elems_lim[i].insert( {0, 0} );
}

PeriodClausesBuffer::PeriodClausesBuffer() : nbThreads(0) {}

void PeriodClausesBuffer::setNbThreads(int _nbThreads) {
    nbThreads = _nbThreads;
    lastOfThread.growTo(_nbThreads);
    for (int i = 0; i < nbThreads; i++) lastOfThread[i] = {i?0:1, 0};
    elems.growTo(nbThreads);
    elems_lim.growTo(nbThreads);
    for (int i = 0; i < nbThreads; i++) elems_lim[i].insert( {0, 0} );
}


bool PeriodClausesBuffer::pushClause(int threadId, Clause & c) {
    assert( 0 <= threadId && threadId < nbThreads );

    // 追加
    // 0：節サイズ
    // 1: lit0
    // 2: lit2
    // ...
    elems[threadId].insert( c.size() );
    for (int i = 0; i < c.size(); ++i)
        elems[threadId].insert( toInt(c[i]) );

    return true;
}

bool PeriodClausesBuffer::pushClause(int threadId, Lit lit) {
    assert( 0 <= threadId && threadId < nbThreads );

    // 追加
    elems[threadId].insert( 1 );
    elems[threadId].insert( toInt(lit) );

    return true;
}

// threadId : import先のスレッド番号
// marginedPeriod : マージン分減算済みのimport元ピリオド

bool PeriodClausesBuffer::getClause(int threadId, int period, int margin, int & threadOrigin, vec<Lit> & resultClause) {
    int marginedPeriod = period - margin;

    if ( marginedPeriod < 0 ) return false;
    assert(marginedPeriod >= 0);

    LastRead thislast = lastOfThread[threadId];    // どのスレッドのどの節まで受け取ったかを表す

    // 節がない時の処理
    thislast = ifEmptyClauseCurrentPeriod( threadId, period, margin, thislast);

    if ( isEndOfCurrentThread(thislast) ) {
        lastOfThread[threadId] = {threadId?0:1, 0};
        return false;
    }

    assert( threadId  != thislast.thread );
    assert( nbThreads != thislast.thread );

    threadOrigin = thislast.thread;

    // 節抽出
    int              bias  = getPeriodIndex (thislast.thread, marginedPeriod);
    Queue<uint32_t>& elem  = elems[thislast.thread];
    int              csize = elem[bias + thislast.index];  // ここで添え字が size を超えることがある

    resultClause.clear();
    for (int i = 0; i < csize; ++i) {
        resultClause.push( toLit( elem[bias + thislast.index + i + 1]));
    }

    // 次のインデックスに移動
    lastOfThread[threadId] = nextIndex(threadId, period, margin, thislast);

    return true;
}


bool PeriodClausesBuffer::nextPeriodProcess(int threadId, int margin) {
    ElemsLim last_lim = elems_lim[threadId][ elems_lim[threadId].size()-1 ];
    elems_lim[threadId].insert( ElemsLim{last_lim.period+1, elems[threadId].size()} );

    assert( elems_lim[threadId].size() > 1 );
    if ( last_lim.period - 10*margin == elems_lim[threadId].peek().period ) {
        ElemsLim second_lim = elems_lim[threadId][1];
        int bias = second_lim.index;

        for (int i = 0; i < bias; ++i) elems[threadId].pop();

        elems_lim[threadId].pop();
        for (int i = 0; i < elems_lim[threadId].size(); ++i)  elems_lim[threadId][i].index -= bias;
    }
    return true;
}


PeriodClausesBuffer::LastRead PeriodClausesBuffer::nextIndex (int threadId, int period, int margin, LastRead lastRead) {
    int marginedPeriod = period - margin;
    int bias  = getPeriodIndex(lastRead.thread, marginedPeriod);
    int csize = elems[lastRead.thread][bias + lastRead.index]; // csize
    lastRead.index += csize + 1;

    assert( bias + lastRead.index <= getPeriodIndex(lastRead.thread, marginedPeriod+1));

    // 現スレッドピリオド分の節は終わった
    if ( bias + lastRead.index == getPeriodIndex(lastRead.thread, marginedPeriod+1) ) {
        // 次のスレッドに移動
        lastRead.thread ++;
        // 自分のスレッドだった場合次のスレッドにさらに移動
        if ( lastRead.thread == threadId ) lastRead.thread ++;

        // lastRead.threadが lastRead.thread == nbThreadsとなることがあるが，正常
        // import終了判定はlastRead.thread == nbThreadsで行なう

        // ピリオドの先頭節に移動
        lastRead.index = 0;
    }

    return lastRead;
}

PeriodClausesBuffer::LastRead PeriodClausesBuffer::ifEmptyClauseCurrentPeriod(int threadId, int period, int margin, LastRead lastRead) {
    int marginedPeriod = period - margin;
    while ( lastRead.thread < nbThreads  &&
            getPeriodIndex(lastRead.thread, marginedPeriod) == getPeriodIndex(lastRead.thread, marginedPeriod+1 ) ){
        lastRead.thread ++;
        if (lastRead.thread == threadId) lastRead.thread ++;
    }
    return lastRead;
}


bool PeriodClausesBuffer::isEndOfCurrentThread (LastRead lastRead) {
    return lastRead.thread == nbThreads;
}


int PeriodClausesBuffer::getPeriodIndex (int threadId, int period) {
    assert( elems_lim[threadId].size() > 0 );

    for (int i = 0; i < elems_lim[threadId].size(); ++i )
        if (elems_lim[threadId][i].period == period)
            return elems_lim[threadId][i].index;

    return elems[threadId].size();    // この添え字の場所にはデータは存在しない
}

void PeriodClausesBuffer::debug_print() {
//    printf("  elems[i][j] \n");
//    for (int i = 0; i < nbThreads; ++i) {
//        for (int j = 0; j < elems[i].size(); ++j) {
//            printf("%6u", elems[i][j]);
//        }
//        printf("\n");
//    }
    for (int i = 0; i < nbThreads; ++i) {
      printf("elems[%d].size() = %5d\n", i, elems[i].size());
    }
    printf("\n");

    for (int i = 0; i < nbThreads; ++i) {
      printf("elems_lim[%d]: ", i);
        for (int j = 0; j < elems_lim[i].size(); ++j) {
            printf("{%3d, %5d},", elems_lim[i][j].period, elems_lim[i][j].index);
        }                                                                                                                                                                                              printf("\n");
    }
    for (int i = 0; i < nbThreads; ++i) {
      printf("lastOfThread[%d] = {%d, %5d}\n", i, lastOfThread[i].thread, lastOfThread[i].index);
    }
    printf("\n\n\n");
}

//=================================================================================================

