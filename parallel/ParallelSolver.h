/**************************************************************************************[ParallelSolver.h]
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

#ifndef PARALLELSOLVER_H
#define	PARALLELSOLVER_H

#include "pthread_barrier.h"    // added by nabesima for mac os x
#include "../core/SolverTypes.h"
#include "../core/Solver.h"
#include "../simp/SimpSolver.h"
#include "../parallel/SharedCompanion.h"
#include "../mtl/IntSet.h"  // added by kanbara

namespace Glucose {

    enum ParallelStats{
        nbexported=coreStatsSize,
        nbimported,
        nbexportedunit,
        nbimportedunit,
        nbimportedInPurgatory,
        nbImportedGoodClauses
    } ;

#define parallelStatsSize (coreStatsSize + 6)

// added by nabesima
#define DBG_EXCHANGE        1
#define DBG_SPEED           2
#define DBG_SEARCH          4
#define DBG_PERIOD          8
#define DBG_PERIOD_STATS    16

//=================================================================================================
    //class MultiSolvers;
    //class SolverCompanion;
    //   class MultiSolvers;

    class ParallelSolver : public SimpSolver {
        friend class MultiSolvers;
        friend class SolverCompanion;
        friend class SharedCompanion;
//    friend class ReasoningCompanion;
//    friend class SolverConfiguration;

    protected :
        // Multithread :
        int		thn; // internal thread number
        //MultiSolvers* belongsto; // Not working (due to incomplete types)
        SharedCompanion *sharedcomp;

        // modified by nabesima
        //bool coreFUIP; // true if one core is specialized for branching on all FUIP
        //bool ImTheSolverFUIP;
        pthread_mutex_t *pmfinished; // mutex on which main process may wait for... As soon as one process finishes it release the mutex
        pthread_cond_t *pcfinished; // condition variable that says that a thread as finished

        // added by nabesima for incremental SAT solving
        vec<Lit> ps_assumptions;    // 'ps' means parallel solver
        bool     ps_do_simp;
        bool     ps_turn_off_simp;

    public:
        // Constructor/Destructor:
        //
        ParallelSolver(int threadId);
        ParallelSolver(const ParallelSolver &s);
        ~ParallelSolver();

        /**
         * Clone function
         */
        virtual Clone* clone() const {
            return  new ParallelSolver(*this);
        }

        int  threadNumber  ()      const;
        void setThreadNumber(int i);
        void reportProgress();
        void reportProgressArrayImports(vec<unsigned int> &totalColumns);
        virtual void reduceDB();
        virtual lbool         solve_                   (bool do_simp = true, bool turn_off_simp = false);

        // added by nabesima
        void setAssumptions     (vec<Lit>& assumps)  { assumps.copyTo(ps_assumptions); }
        void setDoSimp          (bool do_simp)       { ps_do_simp = do_simp; }
        void setTurnOffSimp     (bool turn_off_simp) { ps_turn_off_simp = turn_off_simp; }
        vec<Lit>& getAssumptions()                   { return ps_assumptions; }
        bool getDoSimp          () const             { return ps_do_simp; }
        bool getTurnOffSimp     () const             { return ps_turn_off_simp; }

        // added by nabesima
        void incNumLiveThreads() { sharedcomp->incNumLiveThreads(); }
        void decNumLiveThreads() { sharedcomp->decNumLiveThreads(); }
        bool hasNoLiveThreads()  { return sharedcomp->hasNoLiveThreads(); }

        vec<Lit>    importedClause; // Temporary clause used to copy each imported clause
        unsigned int    goodlimitlbd; // LBD score of the "good" clauses, locally
        int    goodlimitsize;
        bool purgatory; // mode of operation
        bool shareAfterProbation; // Share any none glue clause only after probation (seen 2 times in conflict analysis)
        bool plingeling; // plingeling strategy for sharing clauses (experimental)
        int nbTimesSeenBeforeExport;
        // Stats front end
//    uint64_t   getNbExported() { return nbexported;}
        //   uint64_t   getNbImported() { return nbimported;}
        //   uint64_t   getNbExportedUnit() {return nbexportedunit;}

        uint32_t  firstSharing, limitSharingByGoodLBD, limitSharingByFixedLimitLBD, limitSharingByFixedLimitSize;
        //uint32_t  probationByFollowingRoads, probationByFriend;    // modified by nabesima
        //uint32_t  survivorLayers; // Number of layers for a common clause to survive  // modified by nabesima
        bool dontExportDirectReusedClauses ; // When true, directly reused clauses are not exported
        uint64_t nbNotExportedBecauseDirectlyReused;

        vec<uint32_t> goodImportsFromThreads; // Stats of good importations from other threads

        // added by nabesima for margin strategy
        int64_t   margin;
        int64_t   periods;
        int       prd_type;
        uint64_t  base_conflicts;
        uint64_t  next_conflicts;
        uint64_t  base_lit_scans;
        uint64_t  num_lit_scans;
        uint64_t  next_lit_scans;
        void incLitScans(int d);

        // added by nabesima for block-based period
        double      prd_unit_time;
        double      prd_next_time;
        struct PrdCoeff {
            CoreStats index;
            double    value;
            PrdCoeff() : index(static_cast<CoreStats>(0)), value(0) {}
            PrdCoeff(CoreStats i, double v) : index(i), value(v) {}
        };
        vec<PrdCoeff>         prd_coeff;
        int                   prd_stats_qsize;
        vec<bqueue<uint64_t>> prd_stats_qs;

        // added by nabesima for forced applying imported clauses
        int64_t last_applied_periods;
        bool   fapp_imp_clauses;
        double fapp_prd_rate;
        int    fapp_cla_lb;

        // added by nabesima for exporting to outside of the system
        vec<Lit> unariesExpOutside;

        // added by nabesima for debug
        FILE   *dbg_log;
        int     dbg_type;
        double  dbg_start_time, dbg_next_time;
        double  last_period_time;

        // added by kanbara
        vec<double> learnts_used_rate, import_used_rate, overlap_rate;
        vec<Lit> sort_tmp;
        IntSet table;
        bool div_strategy;
        int  div_periods;
        double used_rate_weight, overlap_rate_threshold, delta_rnd_freq, min_rnd_freq, max_rnd_freq;

        // added by nabesima
        virtual Lit  pickBranchLit();
        virtual bool addClause_(vec<Lit>& ps);

        virtual void parallelExportClauseDuringConflictAnalysis(Clause &c,CRef confl);
        virtual bool parallelImportClauses(); // true if the empty clause was received
        virtual void parallelImportUnaryClauses();
        virtual void parallelExportUnaryClause(Lit p);
        virtual void parallelExportClauseDuringSearch(Clause &c);
        //virtual bool parallelJobIsFinished();  // modified by nabesima
        virtual bool panicModeIsEnabled();

        bool shareClause(Clause & c); // true if the clause was succesfully sent


        // added by gotou
        Queue<uint32_t> imported_clauses_before_applying;
        virtual void moveToNextPeriod();
        virtual bool applyImportedClauses();
        // added by nabesima
        virtual bool shouldApplyImportedClauses();
        virtual void completeCurrPeriod();
        virtual bool shouldFinish();
                void initPeriodCoeffs();
        virtual bool updatePeriod();
        virtual void searchLoopPrehook();
        virtual void analyzePreHook(CRef confl);
        virtual void analyzePostHook(vec<Lit>& learnt, unsigned int lbd);
        virtual void preprocessNewPeriodHook();
        virtual void debugHook(const char* msg);



        // add by kanbara
        virtual bool calculateClausesOverlapRate();
        virtual bool calculateClausesUsedRate();
        virtual void adjustStrategy(bool plan);
        virtual void printConfData();

        uint32_t calculateHash(Clause &c);
        static uint32_t calculateHash(Lit l);
    };


    inline int      ParallelSolver::threadNumber  ()      const   {return thn;}
    inline void     ParallelSolver::setThreadNumber (int i)       {thn = i;}
}
#endif	/* PARALLELSOLVER_H */
