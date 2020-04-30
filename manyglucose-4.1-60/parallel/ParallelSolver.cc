/***************************************************************************************[ParallelSolver.cc]
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

#include "../parallel/ParallelSolver.h"
#include "../mtl/Sort.h"
#include "../utils/System.h"
#include "../mtl/IntSet.h"

using namespace Glucose;

//=====================================================================
// == Options

const char* _cunstable = "CORE/PARALLEL -- UNSTABLE FEATURES";
const char* _parallel = "PARALLEL";

extern BoolOption opt_dontExportDirectReusedClauses; // (_cunstable, "reusedClauses",    "Don't export directly reused clauses", false);
extern BoolOption opt_plingeling; // (_cunstable, "plingeling",    "plingeling strategy for sharing clauses (exploratory feature)", false);

// added by nabesima
const char* _det = "DETERMINISTIC";
extern IntOption    opt_margin;
extern IntOption    opt_prd_type;
extern Int64Option  opt_prd_confs;
extern Int64Option  opt_prd_lit_scans;
extern DoubleOption opt_prd_unit_time;
extern IntOption    opt_prd_stats_qsize;
extern BoolOption   opt_fapp_imp_cla;
extern DoubleOption opt_fapp_prd_rate;
extern IntOption    opt_fapp_cla_lb;

// added by kanbara
extern BoolOption   opt_div_strategy;
extern IntOption    opt_div_periods;
extern DoubleOption opt_used_rate_wt;
extern DoubleOption opt_overlap_rate_th;
extern DoubleOption opt_delta_rnd_freq;
extern DoubleOption opt_max_rnd_freq;
extern DoubleOption opt_min_rnd_freq;

//=====================================================================

//=====================================================================


ParallelSolver::ParallelSolver(int threadId) :
    SimpSolver()
    , thn(threadId) // The thread number of this solver
    , sharedcomp(NULL)        // added by nabesima
    , pmfinished(NULL)        // added by nabesima
    , pcfinished(NULL)        // added by nabesima
    , ps_do_simp(true)        // added by nabesima
    , ps_turn_off_simp(false) // added by nabesima
    , goodlimitlbd(7)
    , goodlimitsize(25)
    , purgatory(true)
    , shareAfterProbation(!opt_plingeling) // only share clauses after probation
    , plingeling(opt_plingeling)
    , nbTimesSeenBeforeExport(2)
    , firstSharing(5000) // Strong limit : do not share anything (except unary clauses) before this number of conflicts
    , limitSharingByGoodLBD(true) // Moving limit of what a good LBD is (median value of last learnt clauses set)
    , limitSharingByFixedLimitLBD(0) // No fixed bound (like 8 in plingeling)
    , limitSharingByFixedLimitSize(0) // No fixed boud (like 40 in plingeling)
    , dontExportDirectReusedClauses(opt_dontExportDirectReusedClauses)
    , nbNotExportedBecauseDirectlyReused(0)
    // added by nabesima
    , margin(opt_margin)
    , periods(0)
    , prd_type(opt_prd_type)
    , base_conflicts(opt_prd_confs)
    , next_conflicts(UINT64_MAX)
    , base_lit_scans(opt_prd_lit_scans)
    , num_lit_scans(0)
    , next_lit_scans(UINT64_MAX)
    , prd_unit_time(opt_prd_unit_time)
    , prd_next_time(opt_prd_unit_time)
    , prd_stats_qsize(opt_prd_stats_qsize)
    , last_applied_periods(0)
    , fapp_imp_clauses(opt_fapp_imp_cla)
    , fapp_prd_rate(opt_fapp_prd_rate)
    , fapp_cla_lb(opt_fapp_cla_lb)
    , dbg_log(NULL)
    , dbg_type(0)
    , dbg_start_time(0.0)
    , dbg_next_time(0.0)
    , last_period_time(0.0)

    // added by kanbara
    , div_strategy(opt_div_strategy)
    , div_periods(opt_div_periods)
    , used_rate_weight(opt_used_rate_wt)
    , overlap_rate_threshold(opt_overlap_rate_th)
    , delta_rnd_freq(opt_delta_rnd_freq)
    , min_rnd_freq(opt_min_rnd_freq)
    , max_rnd_freq(opt_max_rnd_freq)
{
    useUnaryWatched = true; // We want to use promoted clauses here !
    stats.growTo(parallelStatsSize,0);

    imported_clauses_before_applying.clear();   // added by gotou
    initPeriodCoeffs();                         // added by nabesima
}

ParallelSolver::~ParallelSolver() {
    printf("c Solver of thread %d ended.\n", thn);
    fflush(stdout);
    // added by nabesima
    if (dbg_log) { fflush(dbg_log); fclose(dbg_log); }
}

ParallelSolver::ParallelSolver(const ParallelSolver &s) :
    SimpSolver(s)
    , thn(-1)           // added by nabesima
    , sharedcomp(s.sharedcomp)
    , pmfinished(NULL)        // added by nabesima
    , pcfinished(NULL)        // added by nabesima
    , ps_do_simp(true)        // added by nabesima
    , ps_turn_off_simp(false) // added by nabesima
    , goodlimitlbd(s.goodlimitlbd)
    , goodlimitsize(s.goodlimitsize)
    , purgatory(s.purgatory)
    , shareAfterProbation(s.shareAfterProbation) // only share clauses after probation
    , plingeling(s.plingeling)
    , nbTimesSeenBeforeExport(2)
    , firstSharing(s.firstSharing) // Strong limit : do not share anything (except unary clauses) before this number of conflicts
    , limitSharingByGoodLBD(s.limitSharingByGoodLBD) // Moving limit of what a good LBD is (median value of last learnt clauses set)
    , limitSharingByFixedLimitLBD(s.limitSharingByFixedLimitLBD) // No fixed bound (like 8 in plingeling)
    , limitSharingByFixedLimitSize(s.limitSharingByFixedLimitSize) // No fixed boud (like 40 in plingeling)
    , dontExportDirectReusedClauses(s.dontExportDirectReusedClauses)
    , nbNotExportedBecauseDirectlyReused(s.nbNotExportedBecauseDirectlyReused)
    // added by nabesima
    , margin(s.margin)
    , periods(s.periods)
    , prd_type(s.prd_type)
    , base_conflicts(s.base_conflicts)
    , next_conflicts(s.next_conflicts)
    , base_lit_scans(s.base_lit_scans)
    , num_lit_scans(0)
    , next_lit_scans(s.next_lit_scans)
    , prd_unit_time(s.prd_unit_time)
    , prd_next_time(s.prd_unit_time)
    , prd_stats_qsize(s.prd_stats_qsize)
    , last_applied_periods(s.last_applied_periods)
    , fapp_imp_clauses(s.fapp_imp_clauses)
    , fapp_prd_rate(s.fapp_prd_rate)
    , fapp_cla_lb(s.fapp_cla_lb)
    , dbg_log(NULL)
    , dbg_type(0)
    , dbg_start_time(0.0)
    , dbg_next_time(0.0)
    , last_period_time(0.0)

    // added by kanbara
    , div_strategy(opt_div_strategy)
    , div_periods(opt_div_periods)
    , used_rate_weight(opt_used_rate_wt)
    , overlap_rate_threshold(opt_overlap_rate_th)
    , delta_rnd_freq(opt_delta_rnd_freq)
    , min_rnd_freq(opt_min_rnd_freq)
    , max_rnd_freq(opt_max_rnd_freq)
{
    s.goodImportsFromThreads.memCopyTo(goodImportsFromThreads);
    useUnaryWatched = s.useUnaryWatched;
    s.stats.copyTo(stats);
    s.elimclauses.copyTo(elimclauses); // This should be done more efficiently some day

    imported_clauses_before_applying.clear();   // added by gotou

    // added by nabesima
    s.ps_assumptions.copyTo(ps_assumptions);
    s.prd_coeff.copyTo(prd_coeff);
    prd_stats_qs.growTo(s.prd_stats_qs.size());
    for (int i=0; i < prd_stats_qs.size(); i++)
        prd_stats_qs[i].initSize(prd_stats_qsize);
    num_lit_scans = s.num_lit_scans;
}


// modified by nabesima
// Strategy to reduce unary watches list
//struct reduceDB_oneWatched_lt {
//    ClauseAllocator& ca;
//
//    reduceDB_oneWatched_lt(ClauseAllocator& ca_) : ca(ca_) {
//    }
//
//    bool operator()(CRef x, CRef y) {
//
//        // Main criteria... Like in MiniSat we keep all binary clauses
//        if (ca[x].size() > 2 && ca[y].size() == 2) return 1;
//
//        if (ca[y].size() > 2 && ca[x].size() == 2) return 0;
//        if (ca[x].size() == 2 && ca[y].size() == 2) return 0;
//
//        // Second one  based on literal block distance
//        if (ca[x].size() > ca[y].size()) return 1;
//        if (ca[x].size() < ca[y].size()) return 0;
//
//        if (ca[x].lbd() > ca[y].lbd()) return 1;
//        if (ca[x].lbd() < ca[y].lbd()) return 0;
//
//        // Finally we can use old activity or size, we choose the last one
//        return ca[x].activity() < ca[y].activity();
//        //return x->size() < y->size();
//
//        //return ca[x].size() > 2 && (ca[y].size() == 2 || ca[x].activity() < ca[y].activity()); }
//    }
//};
struct reduceDB_oneWatched_lt {
    ClauseAllocator& ca;
    int& num_comps;

    reduceDB_oneWatched_lt(ClauseAllocator& ca_, int& num_comps_) : ca(ca_), num_comps(num_comps_) {
    }

    bool operator()(CRef x, CRef y) {
        num_comps++;

        // Main criteria... Like in MiniSat we keep all binary clauses
        if (ca[x].size() > 2 && ca[y].size() == 2) return 1;

        if (ca[y].size() > 2 && ca[x].size() == 2) return 0;
        if (ca[x].size() == 2 && ca[y].size() == 2) return 0;

        // Second one  based on literal block distance
        if (ca[x].size() > ca[y].size()) return 1;
        if (ca[x].size() < ca[y].size()) return 0;

        if (ca[x].lbd() > ca[y].lbd()) return 1;
        if (ca[x].lbd() < ca[y].lbd()) return 0;

        // Finally we can use old activity or size, we choose the last one
        return ca[x].activity() < ca[y].activity();
        //return x->size() < y->size();

        //return ca[x].size() > 2 && (ca[y].size() == 2 || ca[x].activity() < ca[y].activity()); }
    }
};

// @overide
void ParallelSolver::reduceDB() {
    // added by nabesima
    if (dbg_log && (dbg_type & DBG_SEARCH)) {
        fprintf(dbg_log, "P%lu:C%lu:L%lu REDUCEDB START\n", periods, conflicts, num_lit_scans);
    }
    stats[RedLearntCompCount]++;
    MTIME(seqchrono.start(RedLearntComp));

    int i, j;
    stats[nbReduceDB]++;

    int limit;

    int red_learnt_comp_inside = 0;    // added by nabesima
    if (chanseokStrategy)
        sort(learnts, reduceDBAct_lt(ca, red_learnt_comp_inside));
    else
        sort(learnts, reduceDB_lt(ca, red_learnt_comp_inside));

    // added by nabesima
    MTIME(seqchrono.stop(RedLearntComp));
    stats[sumReductionLitScans] += learnts.size() * log2(learnts.size());
    stats[RedLearntCompInside ] += red_learnt_comp_inside;

    if (!chanseokStrategy && !panicModeIsEnabled()) {
        // We have a lot of "good" clauses, it is difficult to compare them. Keep more !
        if (ca[learnts[learnts.size() / RATIOREMOVECLAUSES]].lbd() <= 3) nbclausesbeforereduce += specialIncReduceDB;
        // Useless :-)
        if (ca[learnts.last()].lbd() <= 5) nbclausesbeforereduce += specialIncReduceDB;
    }
    // Don't delete binary or locked clauses. From the rest, delete clauses from the first half
    // Keep clauses which seem to be usefull (their lbd was reduce during this sequence)

    if (!panicModeIsEnabled()) {
        limit = learnts.size() / 2;
    } else {
        limit = panicModeLastRemoved;
    }
    panicModeLastRemoved = 0;

    // added by nabesima
    stats[RedLearntLoopCount]++;
    MTIME(seqchrono.start(RedLearntLoop));

    uint64_t sumsize = 0;
    for (i = j = 0; i < learnts.size(); i++) {

        Clause& c = ca[learnts[i]];

        if (i == learnts.size() / 2)
            goodlimitlbd = c.lbd();
        sumsize += c.size();
        if (c.lbd() > 2 && c.size() > 2 && c.canBeDel() && !locked(c) && (i < limit)) {
            removeClause(learnts[i]);
            stats[nbRemovedClauses]++;
            panicModeLastRemoved++;
        } else {
            if (!c.canBeDel()) limit++; //we keep c, so we can delete an other clause
            c.setCanBeDel(true); // At the next step, c can be delete
            learnts[j++] = learnts[i];
        }
    }

    // added by nabesima
    MTIME(seqchrono.stop(RedLearntLoop));
    stats[sumReductionLitScans] += learnts.size();
    stats[RedLearntLoopInside ] += learnts.size();

    learnts.shrink(i - j);

    if (learnts.size() > 0)
        goodlimitsize = 1 + (double) sumsize / (double) learnts.size();

    // Special treatment for imported clauses
    if (!panicModeIsEnabled())
        limit = unaryWatchedClauses.size() - (learnts.size() * (chanseokStrategy?4:2));
    else
        limit = panicModeLastRemovedShared;
    panicModeLastRemovedShared = 0;
    if ((unaryWatchedClauses.size() > 100) && (limit > 0)) {

        // added by nabesima
        stats[RedUWLearntCompCount]++;
        MTIME(seqchrono.start(RedUWLearntComp));
        int red_uw_learnt_comp_inside = 0;

        sort(unaryWatchedClauses, reduceDB_oneWatched_lt(ca, red_uw_learnt_comp_inside));

        // added by nabesima
        stats[RedUWLearntCompInside] += red_uw_learnt_comp_inside;
        stats[RedUWLearntLoopCount]++;
        MTIME(seqchrono.toggle(RedUWLearntComp, RedUWLearntLoop));

        for (i = j = 0; i < unaryWatchedClauses.size(); i++) {
            Clause& c = ca[unaryWatchedClauses[i]];
            if (c.lbd() > 2 && c.size() > 2 && c.canBeDel() && !locked(c) && (i < limit)) {
                removeClause(unaryWatchedClauses[i], c.getOneWatched()); // remove from the purgatory (or not)
                stats[nbRemovedUnaryWatchedClauses]++;
                panicModeLastRemovedShared++;
            } else {
                if (!c.canBeDel()) limit++; //we keep c, so we can delete an other clause
                c.setCanBeDel(true); // At the next step, c can be delete
                unaryWatchedClauses[j++] = unaryWatchedClauses[i];
            }
        }
        // added by nabesima
        MTIME(seqchrono.stop(RedUWLearntLoop));
        stats[RedUWLearntLoopInside] += learnts.size();

        unaryWatchedClauses.shrink(i - j);
    }

    checkGarbage();

    if (dbg_log && (dbg_type & DBG_SEARCH)) {
        fprintf(dbg_log, "P%lu:C%lu:L%lu REDUCEDB END\n", periods, conflicts, num_lit_scans);
    }
}

bool ParallelSolver::addClause_(vec<Lit>& ps)
{
#ifndef NDEBUG
    // DEBUG
    for (int i = 0; i < ps.size(); i++)
        if (isEliminated(var(ps[i]))) {
            printf("TH%d: ", thn); printLits(ps); printf("\n");
            printf("TH%d: Var %d is already eliminated!\n", thn, var(ps[i]));
            break;
        }
#endif
    return SimpSolver::addClause_(ps);
}

// These Two functions are useless here !!
void ParallelSolver::reportProgress() {
    // modified by nabesima
    //printf("c | %2d | %6d | %10d | %10d | %8d | %8d | %8d | %8d | %8d | %6.3f |\n",(int)thn,(int)starts,(int)decisions,(int)conflicts,(int)stats[originalClausesSeen],(int)learnts.size(),(int)stats[nbexported],(int)stats[nbimported],(int)stats[nbPromoted],progressEstimate()*100);
    printf("c | %2d | %6d | %10d | %10d | %8d | %8d | %8d | %8d | %8d | %8d | %3d | %6.3f |\n",(int)thn,(int)starts,(int)decisions,(int)conflicts,(int)stats[originalClausesSeen],(int)learnts.size(),(int)stats[nbexported],(int)stats[nbimported],(int)stats[nbPromoted],(int)periods,(int)margin,progressEstimate()*100);
    //printf("c thread=%d confl=%lld starts=%llu reduceDB=%llu learnts=%d broadcast=%llu  blockedReuse=%lld imported=%llu promoted=%llu limitlbd=%llu limitsize=%llu\n", thn, conflicts, starts, nbReduceDB, learnts.size(), nbexported, nbNotExportedBecauseDirectlyReused, nbimported, nbPromoted, goodlimitlbd, goodlimitsize);
}

void ParallelSolver::reportProgressArrayImports(vec<unsigned int> &totalColumns) {
    return ; // TODO : does not currently work
    unsigned int totalImports = 0;
    printf("c %3d | ", thn);
    for (int i = 0; i <  sharedcomp->nbThreads; i++) {
        totalImports += goodImportsFromThreads[i];
        totalColumns[i] += goodImportsFromThreads[i];
        printf(" %8d", goodImportsFromThreads[i]);
    }
    printf(" | %8d\n", totalImports);

}

/*_________________________________________________________________________________________________
  |
  |  panicModeIsEnabled : ()   ->  [bool]
  |
  |  Description:
  |  is panic mode (save memory) is enabled ?
  |________________________________________________________________________________________________@*/

bool ParallelSolver::panicModeIsEnabled() {
    return sharedcomp->panicMode;
}

/*_________________________________________________________________________________________________
  |
  |  parallelImportUnaryClauses : ()   ->  [void]
  |
  |  Description:
  |  import all unary clauses from other cores
  |________________________________________________________________________________________________@*/

void ParallelSolver::parallelImportUnaryClauses() {
    Lit l;
    while ((l = sharedcomp->getUnary(this)) != lit_Undef) {
        if (value(var(l)) == l_Undef) {
            uncheckedEnqueue(l);
            stats[nbimportedunit]++;
        }
    }
}

/*_________________________________________________________________________________________________
  |
  |  parallelImportClauses : ()   ->  [bool]
  |
  |  Description:
  |  import all clauses from other cores
  |  Output : if there is a final conflict
  |________________________________________________________________________________________________@*/

bool ParallelSolver::parallelImportClauses() {

    if (periods < margin) return false;

    stats[ImpThreadLoopCount]++;
    MTIME(seqchrono.start(ImpThreadLoop));
    for (int i=1; i < sharedcomp->nbThreads; i++) {        // i=0 denote the current thread

        int target = (thn + i) % sharedcomp->nbThreads;    // target thread number from which clauses are imported
        assert(target != thn);

        PrdClausesQueue& queue = sharedcomp->getPrdClausesQueue(target);

        if (dbg_log && (dbg_type & DBG_EXCHANGE))
            fprintf(dbg_log, "S%" PRIu64 ",P%" PRIu64 ",C%" PRIu64 ",IMP from thread %d up to period %" PRIu64 "\n", num_lit_scans, periods, conflicts, target, periods - margin);

        stats[ImpQueueLoopCount]++;
        MTIME(seqchrono.start(ImpQueueLoop));
        PrdClauses* p = NULL;
        while ((p = queue.get(thn, periods - margin)) != NULL) {

            PrdClauses& prdClauses = *p;

            if (dbg_log && (dbg_type & DBG_EXCHANGE))
                fprintf(dbg_log, "S%" PRIu64 ",P%" PRIu64 ",C%" PRIu64 ",waitAdditionCompleted thread %d, period %" PRIu64 "\n", num_lit_scans, periods, conflicts, target, prdClauses.period());

            MTIME(seqchrono.start(WaitingTime));    // to avoid adding wait time to ImpQueueLoopTime
            parchrono.start(WaitingTime);

            prdClauses.waitAdditionCompleted();

            parchrono.stop(WaitingTime);
            MTIME(seqchrono.stop(WaitingTime));

            if (dbg_log && (dbg_type & DBG_EXCHANGE))
                fprintf(dbg_log, "S%" PRIu64 ",P%" PRIu64 ",C%" PRIu64 ",IMP from thread %d, period %" PRIu64 "\n", num_lit_scans, periods, conflicts, target, prdClauses.period());

            stats[ImpClauseLoopCount]++;
            MTIME(seqchrono.start(ImpClauseLoop));
            vec<Lit> lits;
            int      idx = 0;
            while (idx < prdClauses.size()) {

                int size = prdClauses[idx++];    // clause size

                stats[ImpClauseLitLoopCount]++;
                MTIME(seqchrono.start(ImpClauseLitLoop));

                // 0 : clause size
                // 1 : imported from thread number
                // 2 : 0th literal of clause
                // 3 : 1st literal of clause
                // ...
                // size + 1 : last literal of clause
                // size + 2 : hash of clause
                imported_clauses_before_applying.insert(size);
                imported_clauses_before_applying.insert(target);
                for (int j=0; j < size; j++)
                    imported_clauses_before_applying.insert(prdClauses[idx++]);
                imported_clauses_before_applying.insert(prdClauses[idx++]);    // added by kanbara for storing hash code

                MTIME(seqchrono.stop(ImpClauseLitLoop));
                num_lit_scans += size;
                stats[sumSharingLitScans    ] += size;
                stats[ImpClauseLitLoopInside] += size;

                // added by nabesima
//                if (dbg_log && (dbg_type & DBG_EXCHANGE)) {
//                    fprintf(dbg_log, "S%" PRIu64 ",P%" PRIu64 ",C%" PRIu64 ",IMP: ", stats[sumLitScans], periods, conflicts);
//                    for (int j=0; j < size; j++) {
//                        Lit p = toLit(prdClauses[idx - size + j]);
//                        fprintf(dbg_log, "%s%d ", sign(p) ? "-" : "", var(p)+1);
//                    }
//                    fprintf(dbg_log, "\n");
//                }
            }
            MTIME(seqchrono.stop(ImpClauseLoop));

            queue.completeExportation(thn, prdClauses);
        }
        MTIME(seqchrono.stop(ImpQueueLoop));
    }
    MTIME(seqchrono.stop(ImpThreadLoop));

    return false;
}


/*_________________________________________________________________________________________________
  |
  |  parallelExportUnaryClause : (Lit p)   ->  [void]
  |
  |  Description:
  |  export unary clauses to other cores
  |________________________________________________________________________________________________@*/

void ParallelSolver::parallelExportUnaryClause(Lit p) {
    // added by nabesima
    stats[ExpUnaryClauseCount]++;
    MTIME(seqchrono.start(ExpUnaryClause));

    // Multithread
    sharedcomp->addLearnt(this, p); // TODO: there can be a contradiction here (two theads proving a and -a)

    // added by nabesima
    MTIME(seqchrono.stop(ExpUnaryClause));

    // added by nabesima
    unariesExpOutside.push(p);

    stats[nbexportedunit]++;

    // added by nabesima
    if (dbg_log && (dbg_type & DBG_EXCHANGE)) {
        fprintf(dbg_log, "UNRY,S%" PRIu64 ",P%" PRIu64 ",C%" PRIu64 ",EXP: ", num_lit_scans, periods, conflicts);
        fprintf(dbg_log, "%s%d ", sign(p) ? "-" : "", var(p)+1);
        fprintf(dbg_log, "\n");
    }
}


/*_________________________________________________________________________________________________
  |
  |  parallelExportClauseDuringSearch : (Clause &c)   ->  [void]
  |
  |  Description:
  |  Verify if a new learnt clause is useful for export
  |  @see search
  |
  |________________________________________________________________________________________________@*/

void ParallelSolver::parallelExportClauseDuringSearch(Clause &c) {
    //
    // Multithread
    // Now I'm sharing the clause if seen in at least two conflicts analysis shareClause(ca[cr]);
    if ((plingeling && !shareAfterProbation && c.lbd() < 8 && c.size() < 40) ||
        (c.lbd() <= 2)) { // For this class of clauses, I'm sharing them asap (they are Glue CLauses, no probation for them)
        if (dbg_log && (dbg_type & DBG_EXCHANGE)) fprintf(dbg_log, "SRCH,");
        stats[ExpClauseDuringSearchCount]++;
        MTIME(seqchrono.start(ExpClauseDuringSearch));
        shareClause(c);
        MTIME(seqchrono.stop(ExpClauseDuringSearch));
        stats[ExpClauseDuringSearchInside] += c.size();
        c.setExported(2);
    }

}

/*_________________________________________________________________________________________________
  |
  |  parallelExportClauseDuringConflictAnalysis : (Clause &c,CRef confl)   ->  [void]
  |
  |  Description:
  |    Verify if the clause using during conflict analysis is good for export
  |    @see : analyze
  |  Output:
  |________________________________________________________________________________________________@*/


void ParallelSolver::parallelExportClauseDuringConflictAnalysis(Clause &c,CRef confl) {
//    if (dbg_log && (dbg_type & DBG_EXCHANGE)) {
//        fprintf(dbg_log, "CCHK,");
//        fprintf(dbg_log, "S%" PRIu64 ",P%" PRIu64 ",C%" PRIu64 ",EXP: ", stats[sumLitScans], periods, conflicts);
//        for (int i=0; i < c.size(); i++)
//            fprintf(dbg_log, "%s%d ", sign(c[i]) ? "-" : "", var(c[i])+1);
//        fprintf(dbg_log, ", wasimp=%d, exported=%d, lbd=%d/%d, size=%d/%d, 1st=%d, 2nd=%d, 3rd=%d\n",
//                c.wasImported(), c.getExported(), c.lbd(), goodlimitlbd, c.size(), goodlimitsize,
//                (shareAfterProbation && c.getExported() != (unsigned int)nbTimesSeenBeforeExport && conflicts > firstSharing),
//                (!c.wasImported() && c.getExported() + 1 == (unsigned int)nbTimesSeenBeforeExport),
//                (c.lbd() == 2 || (c.size() < goodlimitsize && c.lbd() <= goodlimitlbd)));
//    }
    if (dontExportDirectReusedClauses && (confl == lastLearntClause) && (c.getExported() < (unsigned int)nbTimesSeenBeforeExport)) { // Experimental stuff  // modified by gotou (add cast)
        c.setExported(nbTimesSeenBeforeExport);
        nbNotExportedBecauseDirectlyReused++;
    } else if (shareAfterProbation && c.getExported() != (unsigned int)nbTimesSeenBeforeExport && conflicts > firstSharing) {                               // modified by gotou (add cast)
        c.setExported(c.getExported() + 1);
        if (!c.wasImported() && c.getExported() == (unsigned int)nbTimesSeenBeforeExport) { // It's a new interesting clause:                               // modified by gotou (add cast)
            if (c.lbd() == 2 || (c.size() < goodlimitsize && c.lbd() <= goodlimitlbd)) {
                //if (dbg_log && (dbg_type & DBG_EXCHANGE)) fprintf(dbg_log, "CANL,");
                stats[ExpClauseDuringAnalysisCount]++;
                MTIME(seqchrono.start(ExpClauseDuringAnalysis));
                shareClause(c);
                MTIME(seqchrono.stop(ExpClauseDuringAnalysis));
                stats[ExpClauseDuringAnalysisInside] += c.size();
            }
        }
    }

}

/*_________________________________________________________________________________________________
  |
  |  shareClause : (Clause &c)   ->  [bool]
  |
  |  Description:
  |  share a clause to other cores
  | @see : analyze
  |  Output: true if the clause is indeed sent
  |________________________________________________________________________________________________@*/

bool ParallelSolver::shareClause(Clause & c) {

    // added by kanbara
    if(c.getHash() == 0)
      c.setHash(calculateHash(c));

    bool sent = sharedcomp->addLearnt(this, c);
    if (sent) {
        stats[nbexported]++;
        // added by nabesima
        num_lit_scans += c.size();
        stats[sumSharingLitScans] += c.size();
        // added by nabesima
//        if (dbg_log && (dbg_type & DBG_EXCHANGE)) {
//            fprintf(dbg_log, "S%" PRIu64 ",P%" PRIu64 ",C%" PRIu64 ",EXP: ", stats[sumLitScans], periods, conflicts);
//            for (int i=0; i < c.size(); i++)
//                fprintf(dbg_log, "%s%d ", sign(c[i]) ? "-" : "", var(c[i])+1);
//            fprintf(dbg_log, "\n");
//        }
        if (dbg_log && (dbg_type & DBG_SEARCH)) {
            fprintf(dbg_log, "P%lu:C%lu:L%lu EXP ", periods, conflicts, num_lit_scans);
            fprintClause(dbg_log, c);
            fprintf(dbg_log, " LBD=%d\n", c.lbd());
        }
    }
    return sent;
}

/*_________________________________________________________________________________________________
  |
  |  parallelJobIsFinished : ()   ->  [bool]
  |
  |  Description:
  |  Is a core already finish the search
  |
  |________________________________________________________________________________________________@*/

// modified by nabesima
//bool ParallelSolver::parallelJobIsFinished() {
//    // modified by gotou
////    // Parallel: another job has finished let's quit
////    return (sharedcomp->jobFinished());
//    if (sharedcomp->jobFinished() && ( margin == 0 || periods > sharedcomp->jobFinishedPeriod))
//        return true;
//    return false;
//}

// @overide
lbool ParallelSolver::solve_(bool do_simp, bool turn_off_simp) {
    vec<Var> extra_frozen;
    lbool    result = l_True;
    do_simp &= use_simplification;
    if (do_simp){
        // Assumptions must be temporarily frozen to run variable elimination:
        for (int i = 0; i < assumptions.size(); i++){
            Var v = var(assumptions[i]);

            // If an assumption has been eliminated, remember it.
            assert(!isEliminated(v));

            if (!frozen[v]){
                // Freeze and store.
                setFrozen(v, true);
                extra_frozen.push(v);
            }
        }
        result = lbool(eliminate(turn_off_simp));
    }

    model.clear();
    conflict.clear();
    if (!ok) return l_False;

    solves++;

    lbool status = l_Undef;

    // added by nabesima
    next_conflicts = conflicts + base_conflicts;
    next_lit_scans = num_lit_scans + base_lit_scans;
    if (prd_type == prd_type_blocks)
        for (int i=0; i < prd_coeff.size(); i++)
            prd_stats_qs[i].push(stats[prd_coeff[i].index]);
    seqchrono.clearAll();            // to ignore time in preprocessing
    parchrono.start(RunningTime);
    stats[numFreeVars] = nFreeVars();
    stats[numClauses ] = nClauses();
    stats[numLiterals] = stats[clauses_literals];

    // added by nabesima
    if (dbg_log && (dbg_type & DBG_SPEED)) {
        dbg_start_time = realTime();
        dbg_next_time  = dbg_start_time + 0.5;
    }

    // Search:
    int curr_restarts = 0;
    // modified by gotou
//    while (status == l_Undef && !sharedcomp->jobFinished()) {
    //while (status == l_Undef && !parallelJobIsFinished()) {
    // modified by nabesima
    while (status == l_Undef && !shouldFinish()) {
        // added by nabesima
        stats[SearchLoopCount]++;
        MTIME(seqchrono.start(SearchLoop));

        status = search(luby_restart?luby(restart_inc, curr_restarts)*luby_restart_factor:0);  // the parameter is useless in glucose, kept to allow modifications

        // added by nabesima
        MTIME(seqchrono.stop(SearchLoop));

        if(!withinBudget())
            break;

        curr_restarts++;
    }

    // DEBUG
    printf("Thread %d is exited. withinBudget() = %d, shouldFinish() = %d\n", thn, withinBudget(), shouldFinish());
    fflush(stdout);

    completeCurrPeriod();    // added by nabesima

    if (verbosity >= 1)
        printf("c =========================================================================================================\n");

/*
  if (do_simp)
  // Unfreeze the assumptions that were frozen:
  for (int i = 0; i < extra_frozen.size(); i++)
  setFrozen(extra_frozen[i], false);
*/


    // modified by nabesima
//    bool firstToFinish = false;
//    if (status != l_Undef)
//        firstToFinish = sharedcomp->IFinished(this);
//    if (firstToFinish) {
//        printf("c Thread %d is 100%% pure glucose! First thread to finish! (%s answer).\n", threadNumber(), status == l_True ? "SAT" : status == l_False ? "UNSAT" : "UNKOWN");
//        sharedcomp->jobStatus = status;
//    }
    bool foundBetterAns = sharedcomp->IFinished(this, status);
    //printf("c Thread %d (period %" PRIu64 ") finish! (%s answer%s)\n", threadNumber(), periods, status == l_True ? "SAT" : status == l_False ? "UNSAT" : "UNKOWN", foundBetterAns ? ", better one" : "");

    if (foundBetterAns && status == l_True) {
        extendModel();

        // Extend & copy model:
        model.growTo(nVars());
        for (int i = 0; i < nVars(); i++) model[i] = value(i);
    } else if (status == l_False && conflict.size() == 0)
        ok = false;

    pthread_cond_signal(pcfinished);

    parchrono.stop(RunningTime);    // added by nabesima

    //cancelUntil(0);

    return status;
}

// added by nabesima
void ParallelSolver::incLitScans(int d) {
    num_lit_scans += d;
}

// added by nabesima for debug
Lit ParallelSolver::pickBranchLit() {
    Lit p = Solver::pickBranchLit();
    if (dbg_log && (dbg_type & DBG_SEARCH)) {
        fprintf(dbg_log, "P%lu:C%lu:L%lu DEC ", periods, conflicts, num_lit_scans);
        fprintLit(dbg_log, p);
        fprintf(dbg_log, "\n");
    }
    return p;
}


// added by gotou
/*_________________________________________________________________________________________________
  |
  |  nextPeriod : ()   ->  [void]
  |
  |  Description:
  |  Move to the next period.
  |
  |________________________________________________________________________________________________@*/

void ParallelSolver::moveToNextPeriod() {
//    nextperiodprocessing_time.start();
//    if ( margin == 0 )
//        pthread_barrier_wait(&(sharedcomp->barrierMargin0AfterImport));
//
//    sharedcomp->moveToNextPeriod(thn, margin);
//    nextperiodprocessing_time.stop();

    periods ++;

    if (dbg_log && (dbg_type & DBG_SEARCH)) {
        fprintf(dbg_log, "P%lu:C%lu: NEW_PERIOD ", periods, conflicts);
        fprintf(dbg_log, " litScans=%lu, nextLitScans=%lu\n", num_lit_scans, next_lit_scans);
    }
    if (dbg_log && (dbg_type & DBG_PERIOD)) {
        double curr =
                parchrono.getTime(RunningTime) +
                parchrono.getTime(ExchangingTime) +
                parchrono.getTime(PeriodUpdateTime);

        if (periods == 1)
            last_period_time = curr;

        fprintf(dbg_log, "Thn,%d,Prd,%lu,Time,%f\n",
                thn,
                periods,
                curr - last_period_time);
        last_period_time = curr;
    }
    if (dbg_log && (dbg_type & DBG_PERIOD_STATS) && periods % 10 == 0) {
        double waiting_time       = parchrono.getTime(WaitingTime);
        parchrono.start(WaitingTime);
        double running_time       = parchrono.getTime(RunningTime);
        double exchanging_time    = parchrono.getTime(ExchangingTime);
        double period_update_time = parchrono.getTime(PeriodUpdateTime);
        double no_waiting_time    = running_time + exchanging_time + period_update_time;

        if (periods == 10) {
            fprintf(dbg_log,
                    "Thn,Prd,NoWaitingTime,RunningTime,ExchangingTime,PeriodUpdateTime,WaitingTime,"
                    "SearchLoopCount,SearchLoopInside,SearchConflictCount,SearchDecisionCount,"
                    "PropCleanWatchLoopCount,PropCleanWatchLoopInside,PropQueueLoopCount,"
                    "PropBinWatchLoopCount,PropBinWatchLoopInside,PropBinWatchLoopInsideConf,PropBinWatchLoopInsideEnq,PropBinWatchLoopInsideNext,"
                    "PropWatchLoopCount,PropWatchLoopInside,PropWatchLoopInsideNext,PropWatchLoopInsideBlocker,PropClauseLoopCount,PropClauseLoopInside,PropClauseLoopBreakCount,"
                    "UnaryPropWatchLoopCount,UnaryPropWatchLoopInside,UnaryPropWatchLoopInsideBlocker,UnaryPropClauseLoopCount,UnaryPropClauseLoopInside,UnaryPropClauseLoopBreakCount,UnaryPropConfClauseLoopCount,UnaryPropConfClauseLoopInside,"
                    "AnalyzeLoopCount,AnalyzeClauseLoopCount,AnalyzeClauseLoopInside,AnalyzeOutMinLoopCount,AnalyzeOutMinLoopInside,AnalyzeOutLBDLoopCount,AnalyzeOutLBDLoopInside,AnalyzeLastDecLoopCount,AnalyzeLastDecLoopInside,"
                    "LitRedLoopCount,LitRedClauseLoopCount,LitRedClauseLoopInside,MinBinResCount,"
                    "RedLearntCompCount,RedLearntCompInside,RedLearntLoopCount,RedLearntLoopInside,RedUWLearntCompCount,RedUWLearntCompInside,RedUWLearntLoopCount,RedUWLearntLoopInside,"
                    "SimpRmSatClausesLoopCount,SimpClauseLoopCount,SimpClauseLoopInside,SimpVarHeapLoopCount,SimpVarHeapLoopInside,"
                    "ExpUnaryClauseCount,ExpClauseDuringSearchCount,ExpClauseDuringSearchInside,ExpClauseDuringAnalysisCount,ExpClauseDuringAnalysisInside,"
                    "ImpThreadLoopCount,ImpQueueLoopCount,ImpClauseLoopCount,ImpClauseLitLoopCount,ImpClauseLitLoopInside,"
                    "ApplyImportedClausesLoopCount,ApplyImportedClausesLoopInside,"
                    "RelocCleanWatchLoopCount,RelocCleanWatchLoopInside,RelocLoopCount,RelocLoopInside,DetachLoopCount,DetachLoopInside\n");
        }

        fprintf(dbg_log, "%d,%lu,%f,%f,%f,%f,%f",
                thn,
                periods,
                no_waiting_time, running_time, exchanging_time, period_update_time, waiting_time);
        for (int i=SearchLoopCount; i <= DetachLoopInside; i++)
            fprintf(dbg_log, ",%" PRIu64, stats[i]);
        fprintf(dbg_log, "\n");

        parchrono.stop(WaitingTime);
    }

//    waitingNextPeriod = false;
}

bool ParallelSolver::applyImportedClauses() {
    int queue_index = 0;


    if (decisionLevel() > 0) {
        cancelUntil(0);
        stats[numForcedImports]++;
    }

    assert(decisionLevel() == 0);

    stats[ApplyImportedClausesLoopCount]++;
    MTIME(seqchrono.start(ApplyImportedClausesLoop));

    while ( queue_index < imported_clauses_before_applying.size() ) {
        stats[ApplyImportedClausesLoopInside]++;

        int size = imported_clauses_before_applying[queue_index];
        int importedFromThread = imported_clauses_before_applying[queue_index+1];
        importedClause.clear();
        for ( int i = 0; i < size; ++i )
            importedClause.push( toLit(imported_clauses_before_applying[queue_index + 2 + i]) );

        if (importedClause.size() == 0) {
            MTIME(seqchrono.stop(ApplyImportedClausesLoop));
            return true;
        }

        num_lit_scans += size;    // added by nabesima
        stats[sumSharingLitScans] += size;

        // added by gotou
        if (importedClause.size() == 1) {
            if (value(var(importedClause[0])) == l_Undef) {
                uncheckedEnqueue(importedClause[0]);
                stats[nbimportedunit]++;
                if (dbg_log && (dbg_type & DBG_SEARCH)) {
                    fprintf(dbg_log, "P%lu:C%lu: IMP ", periods, conflicts);
                    fprintLits(dbg_log, importedClause);
                    fprintf(dbg_log, "\n");
            }
            }
            // modified by kanbara
            // queue_index += size + 2;
            queue_index += size + 3;
            continue;
        }

        //printf("Thread %d imports clause from thread %d\n", threadNumber(), importedFromThread);
        CRef cr = ca.alloc(importedClause, true, true);
        ca[cr].setLBD(importedClause.size());
        if (plingeling) // 0 means a broadcasted clause (good clause), 1 means a survivor clause, broadcasted
            ca[cr].setExported(2); // A broadcasted clause (or a survivor clause) do not share it anymore
        else {
            ca[cr].setExported(1); // next time we see it in analyze, we share it (follow route / broadcast depending on the global strategy, part of an ongoing experimental stuff: a clause in one Watched will be set to exported 2 when promotted.
        }

        // added by kanbara
        ca[cr].setHash(imported_clauses_before_applying[queue_index + 2 + size]);

        ca[cr].setImportedFrom(importedFromThread);
        if(useUnaryWatched)
            unaryWatchedClauses.push(cr);
        else
            learnts.push(cr);

        if (plingeling || ca[cr].size() <= 2) {//|| importedRoute == 0) { // importedRoute == 0 means a glue clause in another thread (or any very good clause)
            ca[cr].setOneWatched(false); // Warning: those clauses will never be promoted by a conflict clause (or rarely: they are propagated!)
            attachClause(cr);
            stats[nbImportedGoodClauses]++;
        } else {
            if(useUnaryWatched) {
                attachClausePurgatory(cr); //
                ca[cr].setOneWatched(true);
            } else {
                attachClause(cr);
                ca[cr].setOneWatched(false);
            }
            stats[nbimportedInPurgatory]++;
        }
        assert(ca[cr].learnt());
        stats[nbimported]++;

        // added by nabesima
        if (dbg_log && (dbg_type & DBG_SEARCH)) {
            fprintf(dbg_log, "P%lu:C%lu: IMP ", periods, conflicts);
            fprintClause(dbg_log, cr);
            fprintf(dbg_log, "\n");
        }

        // modified by kanbara
        // queue_index += size + 2;
        queue_index += size + 3;
    }
    MTIME(seqchrono.stop(ApplyImportedClausesLoop));

    imported_clauses_before_applying.clear();

    last_applied_periods = periods;    // added by nabesima

    return false;
}

// added by nabesima
bool ParallelSolver::shouldApplyImportedClauses() {
    if (decisionLevel() == 0) return true;
    if (fapp_imp_clauses
            && periods - last_applied_periods >= margin * fapp_prd_rate
            && imported_clauses_before_applying.size() / sharedcomp->nbThreads >= fapp_cla_lb)
        return true;
    return false;
}
void ParallelSolver::completeCurrPeriod() {
    sharedcomp->completeAdditionToCurrPeriod(thn, propagations);
}
bool ParallelSolver::shouldFinish() {
    return sharedcomp->jobFinished(periods);
}
void ParallelSolver::initPeriodCoeffs() {
    // Coefficients are tuned by 2016-2017 instances (4.1-41b)
    // - based on properties whose max running time ratio >= 10%
    prd_coeff.push(PrdCoeff(SearchDecisionCount, 9.710020606567300e-07));
    prd_coeff.push(PrdCoeff(PropBinWatchLoopInsideNext, 2.214945109571396e-08));
    prd_coeff.push(PrdCoeff(PropWatchLoopCount, 2.857831824512144e-07));
    prd_coeff.push(PrdCoeff(PropWatchLoopInsideNext, 6.090353982242565e-09));
    prd_coeff.push(PrdCoeff(PropClauseLoopCount, 2.803365253553293e-07));
    prd_coeff.push(PrdCoeff(PropClauseLoopInside, 2.396763104145166e-08));
    prd_coeff.push(PrdCoeff(UnaryPropWatchLoopCount, 6.568775555448704e-07));
    prd_coeff.push(PrdCoeff(UnaryPropWatchLoopInside, 4.092636578956352e-09));
    prd_coeff.push(PrdCoeff(UnaryPropWatchLoopInsideBlocker, 3.162052910627721e-07));
    prd_coeff.push(PrdCoeff(UnaryPropClauseLoopCount, 2.191431783220884e-19));
    prd_coeff.push(PrdCoeff(UnaryPropClauseLoopInside, 4.845600922167139e-10));
    prd_coeff.push(PrdCoeff(UnaryPropClauseLoopBreakCount, 2.100872711380898e-08));
    prd_coeff.push(PrdCoeff(UnaryPropConfClauseLoopCount, 8.894903073539904e-04));
    prd_coeff.push(PrdCoeff(AnalyzeClauseLoopCount, 2.358954361882080e-07));
    prd_coeff.push(PrdCoeff(LitRedLoopCount, 5.479807925553783e-07));
    prd_coeff.push(PrdCoeff(LitRedClauseLoopCount, 2.027371838581897e-07));
    prd_coeff.push(PrdCoeff(LitRedClauseLoopInside, 1.680037403557589e-08));

    // Coefficients are tuned by 2016-2017 instances (4.1-47)
    // - based on 1% of every period
//    prd_coeff.push(PrdCoeff(SearchDecisionCount, 1.3951181828357395e-06));
//    prd_coeff.push(PrdCoeff(PropBinWatchLoopCount, 8.135371664136694e-07));
//    prd_coeff.push(PrdCoeff(PropBinWatchLoopInsideNext, 1.8730892554223123e-08));
//    prd_coeff.push(PrdCoeff(PropWatchLoopCount, 7.577918838423548e-08));
//    prd_coeff.push(PrdCoeff(PropWatchLoopInside, 9.076478229881346e-10));
//    prd_coeff.push(PrdCoeff(PropWatchLoopInsideNext, 1.925407221305776e-08));
//    prd_coeff.push(PrdCoeff(PropWatchLoopInsideBlocker, 6.085915445895534e-08));
//    prd_coeff.push(PrdCoeff(PropClauseLoopCount, 2.411319168032802e-07));
//    prd_coeff.push(PrdCoeff(PropClauseLoopInside, 2.3072959807639104e-08));
//    prd_coeff.push(PrdCoeff(UnaryPropWatchLoopCount, 2.020543576885977e-08));
//    prd_coeff.push(PrdCoeff(UnaryPropWatchLoopInside, 7.700404788345862e-09));
//    prd_coeff.push(PrdCoeff(UnaryPropWatchLoopInsideBlocker, 3.458215990519878e-07));
//    prd_coeff.push(PrdCoeff(UnaryPropClauseLoopCount, 4.035284562962225e-15));
//    prd_coeff.push(PrdCoeff(UnaryPropClauseLoopInside, 2.289169681349136e-08));
//    prd_coeff.push(PrdCoeff(UnaryPropConfClauseLoopInside, 6.744112813012007e-06));
//    prd_coeff.push(PrdCoeff(AnalyzeClauseLoopCount, 1.1474611778823935e-07));
//    prd_coeff.push(PrdCoeff(AnalyzeLastDecLoopInside, 1.064574011312786e-06));
//    prd_coeff.push(PrdCoeff(LitRedLoopCount, 4.897558305943611e-07));
//    prd_coeff.push(PrdCoeff(LitRedClauseLoopCount, 1.6882172120692342e-07));
//    prd_coeff.push(PrdCoeff(LitRedClauseLoopInside, 1.578702091294888e-08));
//    prd_coeff.push(PrdCoeff(RedLearntCompInside, 4.721851579353707e-08));
//    prd_coeff.push(PrdCoeff(RedUWLearntCompInside, 7.76551129373079e-08));
//    prd_coeff.push(PrdCoeff(SimpClauseLoopCount, 2.3925112776446383e-08));
//    prd_coeff.push(PrdCoeff(SimpClauseLoopInside, 2.2522623491131393e-08));
//    prd_coeff.push(PrdCoeff(SimpVarHeapLoopInside, 4.84554697126839e-08));
//    prd_coeff.push(PrdCoeff(ApplyImportedClausesLoopInside, 1.0577335087191132e-06));
//    prd_coeff.push(PrdCoeff(RelocCleanWatchLoopInside, 1.1245077697681242e-07));
//    prd_coeff.push(PrdCoeff(RelocLoopInside, 1.0675157677349713e-07));

    // Coefficients are tuned by 2017 instances (4.1-48)
    // - based on 50% of every period
//    prd_coeff.push(PrdCoeff(SearchDecisionCount, 2.552682353603104e-06));
//    prd_coeff.push(PrdCoeff(PropBinWatchLoopCount, 6.051850371801265e-07));
//    prd_coeff.push(PrdCoeff(PropBinWatchLoopInsideNext, 5.537462949944591e-08));
//    prd_coeff.push(PrdCoeff(PropWatchLoopCount, 6.614764754310713e-08));
//    prd_coeff.push(PrdCoeff(PropWatchLoopInsideNext, 2.7799315964616305e-08));
//    prd_coeff.push(PrdCoeff(PropClauseLoopCount, 2.9093481210001306e-07));
//    prd_coeff.push(PrdCoeff(PropClauseLoopInside, 1.9262365703339007e-08));
//    prd_coeff.push(PrdCoeff(UnaryPropWatchLoopCount, 1.8462028206993465e-07));
//    prd_coeff.push(PrdCoeff(UnaryPropWatchLoopInsideBlocker, 4.155204583397585e-07));
//    prd_coeff.push(PrdCoeff(UnaryPropClauseLoopCount, 4.4827333452441817e-14));
//    prd_coeff.push(PrdCoeff(UnaryPropClauseLoopInside, 5.027674181252701e-09));
//    prd_coeff.push(PrdCoeff(UnaryPropConfClauseLoopInside, 1.2975277926677068e-05));
//    prd_coeff.push(PrdCoeff(AnalyzeClauseLoopCount, 2.293670149371837e-07));
//    prd_coeff.push(PrdCoeff(LitRedLoopCount, 4.501774495608066e-07));
//    prd_coeff.push(PrdCoeff(LitRedClauseLoopCount, 3.564442409960364e-08));
//    prd_coeff.push(PrdCoeff(LitRedClauseLoopInside, 9.536087793514327e-08));
//    prd_coeff.push(PrdCoeff(RedLearntCompInside, 2.8853476712447138e-08));
//    prd_coeff.push(PrdCoeff(RedUWLearntCompInside, 8.451138311827189e-08));
//    prd_coeff.push(PrdCoeff(SimpClauseLoopInside, 1.8789005628199493e-08));
//    prd_coeff.push(PrdCoeff(SimpVarHeapLoopInside, 1.5251795705963834e-07));
//    prd_coeff.push(PrdCoeff(ApplyImportedClausesLoopInside, 1.114227037065348e-06));
//    prd_coeff.push(PrdCoeff(RelocCleanWatchLoopInside, 1.7337777018241134e-08));
//    prd_coeff.push(PrdCoeff(RelocLoopInside, 1.5356689906349226e-07));

    // Coefficients are tuned by 2017 instances (4.1-49)
    // - based on 10% of every period
    // - exclude periods >= 1.2 sec
    // - reduced the number of periods < 1.0 sec
//    prd_coeff.push(PrdCoeff(SearchDecisionCount, 2.5065649563822687e-06));
//    prd_coeff.push(PrdCoeff(PropBinWatchLoopCount, 7.579655540470209e-07));
//    prd_coeff.push(PrdCoeff(PropBinWatchLoopInsideNext, 4.867457986582924e-08));
//    prd_coeff.push(PrdCoeff(PropWatchLoopCount, 6.653002197758043e-08));
//    prd_coeff.push(PrdCoeff(PropWatchLoopInsideNext, 2.579797646948263e-08));
//    prd_coeff.push(PrdCoeff(PropClauseLoopCount, 3.114517639097563e-07));
//    prd_coeff.push(PrdCoeff(PropClauseLoopInside, 1.784031002755822e-08));
//    prd_coeff.push(PrdCoeff(UnaryPropWatchLoopCount, 2.7997389385243658e-08));
//    prd_coeff.push(PrdCoeff(UnaryPropWatchLoopInsideBlocker, 4.02563429888713e-07));
//    prd_coeff.push(PrdCoeff(UnaryPropClauseLoopCount, 9.946711408311292e-15));
//    prd_coeff.push(PrdCoeff(UnaryPropClauseLoopInside, 8.964304505475538e-09));
//    prd_coeff.push(PrdCoeff(UnaryPropClauseLoopBreakCount, 7.042356149070837e-10));
//    prd_coeff.push(PrdCoeff(UnaryPropConfClauseLoopInside, 1.1971367324300188e-05));
//    prd_coeff.push(PrdCoeff(AnalyzeClauseLoopCount, 2.3112003917131086e-07));
//    prd_coeff.push(PrdCoeff(LitRedLoopCount, 4.411682824259703e-07));
//    prd_coeff.push(PrdCoeff(LitRedClauseLoopCount, 3.505878936092569e-08));
//    prd_coeff.push(PrdCoeff(LitRedClauseLoopInside, 9.811612678394326e-08));
//    prd_coeff.push(PrdCoeff(RedLearntCompInside, 5.6015654830724756e-08));
//    prd_coeff.push(PrdCoeff(RedLearntLoopInside, 2.4355448648453826e-07));
//    prd_coeff.push(PrdCoeff(RedUWLearntCompInside, 7.386281179842156e-08));
//    prd_coeff.push(PrdCoeff(RedUWLearntLoopInside, 3.4451678768411686e-07));
//    prd_coeff.push(PrdCoeff(SimpClauseLoopInside, 1.87926165539707e-08));
//    prd_coeff.push(PrdCoeff(SimpVarHeapLoopInside, 1.3808291844754817e-07));
//    prd_coeff.push(PrdCoeff(ApplyImportedClausesLoopInside, 1.1322609119189405e-06));
//    prd_coeff.push(PrdCoeff(RelocLoopInside, 1.7757497947515505e-07));

    // Coefficients are tuned by 2016 & 2017 instances (4.1-51)
    // - based on the running time of 10 consecutive periods.
//    prd_coeff.push(PrdCoeff(SearchDecisionCount, 9.983206498672585e-07));
//    prd_coeff.push(PrdCoeff(PropBinWatchLoopCount, 6.395135209401361e-07));
//    prd_coeff.push(PrdCoeff(PropWatchLoopCount, 1.0076737317601754e-07));
//    prd_coeff.push(PrdCoeff(PropClauseLoopCount, 2.2958445326502038e-07));
//    prd_coeff.push(PrdCoeff(UnaryPropWatchLoopCount, 1.6239649746506312e-07));
//    prd_coeff.push(PrdCoeff(UnaryPropWatchLoopInsideBlocker, 2.988726316765735e-07));
//    prd_coeff.push(PrdCoeff(UnaryPropClauseLoopCount, 2.2044918939206413e-17));
//    prd_coeff.push(PrdCoeff(UnaryPropClauseLoopBreakCount, 4.39009952599934e-08));
//    prd_coeff.push(PrdCoeff(UnaryPropConfClauseLoopCount, 0.00035535900717009224));
//    prd_coeff.push(PrdCoeff(AnalyzeClauseLoopCount, 9.284141308407916e-08));
//    prd_coeff.push(PrdCoeff(AnalyzeLastDecLoopInside, 1.3610150466130722e-06));
//    prd_coeff.push(PrdCoeff(LitRedLoopCount, 1.6274077531544594e-07));
//    prd_coeff.push(PrdCoeff(LitRedClauseLoopCount, 1.7121744652931586e-07));
//    prd_coeff.push(PrdCoeff(LitRedClauseLoopInside, 7.260101342576074e-09));
//    prd_coeff.push(PrdCoeff(RedUWLearntCompInside, 7.908140232368078e-08));
//    prd_coeff.push(PrdCoeff(SimpClauseLoopInside, 1.7241883463473122e-09));
//    prd_coeff.push(PrdCoeff(SimpVarHeapLoopInside, 2.510111704306389e-07));
//    prd_coeff.push(PrdCoeff(ImpQueueLoopCount, 0.0008366756600355268));
//    prd_coeff.push(PrdCoeff(ImpClauseLoopCount, 2.0259350544252956e-07));
//    prd_coeff.push(PrdCoeff(ApplyImportedClausesLoopInside, 9.536739315067232e-07));
//    prd_coeff.push(PrdCoeff(RelocCleanWatchLoopInside, 1.6377078442017202e-07));
//    prd_coeff.push(PrdCoeff(RelocLoopInside, 9.083620106739547e-08));


    // Ensure the size of period statistics array
    prd_stats_qs.growTo(prd_coeff.size());
    for (int i=0; i < prd_stats_qs.size(); i++)
        prd_stats_qs[i].initSize(prd_stats_qsize);
}
bool ParallelSolver::updatePeriod() {
    if (prd_type == prd_type_conf) {
        if (conflicts >= next_conflicts) {
            next_conflicts = conflicts + base_conflicts;
            return true;
        }
        return false;
    }
    else if (prd_type == prd_type_lit_scan) {
        if (num_lit_scans >= next_lit_scans) {
            stats[sumLitScanFractinos] += num_lit_scans - next_lit_scans;
            next_lit_scans = num_lit_scans + base_lit_scans;
            return true;
        }
        return false;
    }
    else if (prd_type == prd_type_blocks) {
        double time = 0.0;
        for (int i=0; i < prd_coeff.size(); i++)
            time += (stats[prd_coeff[i].index] - prd_stats_qs[i].peek()) * prd_coeff[i].value;
        if (time >= prd_next_time) {
            if (prd_stats_qs[0].size() < prd_stats_qsize)
                prd_next_time += prd_unit_time;
            for (int i=0; i < prd_coeff.size(); i++)
                prd_stats_qs[i].push(stats[prd_coeff[i].index]);
            return true;
        }

        return false;
    }
    assert(false);
    return false;
}

void ParallelSolver::searchLoopPrehook() {
    if (dbg_log && (dbg_type & DBG_SPEED)) {
        double curr = realTime() - parchrono.getTime(WaitingTime);
        if (curr > dbg_next_time) {
            fprintf(dbg_log, "thread:%d\t"
                             "time:%f\t"
                             "period:%" PRIu64 "\t"
                             "propagations:%" PRIu64 "\t"
                             "conflicts:%" PRIu64 "\t"
                             "litfullscans:%" PRIu64 "\t"
                             "len:%f\t"
                             "\n",
                             thn,
                             curr - dbg_start_time,
                             periods,
                             propagations,
                             conflicts,
                             num_lit_scans,
                             (double)stats[learnts_literals] / (learnts.size() + permanentLearnts.size() + unaryWatchedClauses.size())
            );
            dbg_next_time = curr + 0.5;
        }
    }
}
void ParallelSolver::analyzePreHook(CRef confl) {
  if (dbg_log && (dbg_type & DBG_SEARCH)) {
    fprintf(dbg_log, "P%lu:C%lu:L%lu CONFLICT ", periods, conflicts, num_lit_scans);
    fprintClause(dbg_log, ca[confl]);
    fprintf(dbg_log, "\n");
  }
}
void ParallelSolver::analyzePostHook(vec<Lit>& learnt, unsigned int lbd) {
  if (dbg_log && (dbg_type & DBG_SEARCH)) {
    fprintf(dbg_log, "P%lu:C%lu:L%lu LEARNT ", periods, conflicts, num_lit_scans);
    fprintLits(dbg_log, learnt);
    fprintf(dbg_log, " LBD=%d\n", lbd);
  }
}
void ParallelSolver::preprocessNewPeriodHook() {
    // added by kanbara
    if (div_strategy && unaryWatchedClauses.size() > 0 && periods % div_periods == 0) {
        usedrate_time.start();
        bool calcOverlap = calculateClausesUsedRate();
        usedrate_time.stop();

        if (calcOverlap) {
            overlaprate_time.start();
            bool plan = calculateClausesOverlapRate();
            overlaprate_time.stop();

            adjustStrategy(plan);
        }

        rnd_freq.push(random_var_freq);
    }
}

void ParallelSolver::debugHook(const char* msg) {
    if (dbg_log && (dbg_type & DBG_SEARCH))
        fprintf(dbg_log, "P%lu:C%lu:L%lu %s\n", periods, conflicts, num_lit_scans, msg);
}


// added by kanbara
// return true : check overlap rate
bool ParallelSolver::calculateClausesUsedRate() {

  int numUsedLearnts = 0;
  for(int i = 0; i < learnts.size(); i++) {
    Clause &c = ca[learnts[i]];
    if(c.isUsed()) {
      numUsedLearnts++;
      c.setUsed(false);
    }
  }

  double learntsUsedRate = (double)numUsedLearnts / (double)learnts.size();
  learnts_used_rate.push(learntsUsedRate);

  int numUsedImports= 0;
  for(int i = 0; i < unaryWatchedClauses.size(); i++) {
    Clause &c = ca[unaryWatchedClauses[i]];
    if(c.isUsed()) {
        numUsedImports++;
      c.setUsed(false);
    }
  }
  double importsUsedRate = (double)numUsedImports / (double)unaryWatchedClauses.size();
  import_used_rate.push(importsUsedRate);

  if (used_rate_weight == 0) {
    return true;
  }

  return (learntsUsedRate * used_rate_weight >= importsUsedRate);
}

// added by kanbara
// return true  : give diversity
// return false : give intensity
bool ParallelSolver::calculateClausesOverlapRate() {

    table.clear();

    for (int i = 0; i < learnts.size(); i++) {
        Clause &c = ca[learnts[i]];
        // if hash is 0, hash has not calculated yet.
        if (c.getHash() == 0)
            c.setHash(calculateHash(c));

        table.add(c.getHash());
    }

    int numOfMatchClauses = 0;
    for (int i = 0; i < unaryWatchedClauses.size(); i++) {
        Clause &c = ca[unaryWatchedClauses[i]];
        if (c.getHash() == 0)
            c.setHash(calculateHash(c));

        if (table.has(c.getHash()))
            numOfMatchClauses++;
    }

    double rate = (double) numOfMatchClauses / (double) unaryWatchedClauses.size();
    overlap_rate.push(rate);
    return (rate >= overlap_rate_threshold);
}

// added by kanbara
uint32_t ParallelSolver::calculateHash(Clause &c) {

  sort_tmp.clear();
  for(int i = 0; i < c.size(); i++)
    sort_tmp.push(c[i]);

  sort(sort_tmp);

  uint32_t hash = (uint32_t)sort_tmp[0].x;
  for(int i = 1; i < sort_tmp.size(); i++)
    hash = hash * 31 + sort_tmp[i].x;

  return hash;
}
uint32_t ParallelSolver::calculateHash(Lit l) {
  return l.x + 31;
}

// added by kanbara
// plan == true  : give diversity
// plan == false : give intensity
void ParallelSolver::adjustStrategy(bool plan) {

  if(plan)
    random_var_freq + delta_rnd_freq < max_rnd_freq ? random_var_freq += delta_rnd_freq : random_var_freq = max_rnd_freq;
  else
    random_var_freq - delta_rnd_freq > min_rnd_freq ? random_var_freq -= delta_rnd_freq : random_var_freq = min_rnd_freq;
}

void ParallelSolver::printConfData() {
  printf("Thread:%d Period:%" PRIu64 " Conflict:%" PRIu64 " random_seed: %f ",
     thn,
     periods,
     conflicts,
     random_seed
     );
}
