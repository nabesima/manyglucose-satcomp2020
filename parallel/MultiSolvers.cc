/***************************************************************************************[MultiSolvers.cc]
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

#include <pthread.h>
#include "../parallel/MultiSolvers.h"
#include "../mtl/Sort.h"
#include "../utils/System.h"
#include "../simp/SimpSolver.h"
#include <errno.h>
#include <string.h>
#include <string>
#include "../parallel/SolverConfiguration.h"
#include <map>      // added by nabesima
#include <set>      // added by nabesima

using namespace Glucose;

extern const char *_parallel;
extern const char *_cunstable;
// Options at the parallel solver level
static IntOption opt_nbsolversmultithreads(_parallel, "nthreads", "Number of core threads for syrup (0 for automatic)", 0);
static IntOption opt_maxnbsolvers(_parallel, "maxnbthreads", "Maximum number of core threads to ask for (when nbthreads=0)", 4);
// modified by nabesima
//static IntOption opt_maxmemory(_parallel, "maxmemory", "Maximum memory to use (in Mb, 0 for no software limit)", 20000);
//static IntOption opt_statsInterval(_parallel, "statsinterval", "Seconds (real time) between two stats reports", 5);
static IntOption opt_maxmemory(_parallel, "maxmemory", "Maximum memory to use (in Mb, 0 for no software limit)", 0);
static IntOption opt_statsInterval(_parallel, "statsinterval", "Seconds (real time) between two stats reports", 30);

// Configuration of solver threads. added by nabesima
IntOption opt_ps_conf(_parallel, "ps-conf", "Parallel solvers configuration (0=same, 1=syrup, 2=syrup + diff-seed)", 2);

//
// Shared with ClausesBuffer.cc
BoolOption opt_whenFullRemoveOlder(_parallel, "removeolder", "When the FIFO for exchanging clauses between threads is full, remove older clauses", false);
IntOption opt_fifoSizeByCore(_parallel, "fifosize", "Size of the FIFO structure for exchanging clauses between threads, by threads", 100000);
//
// Shared options with Solver.cc
BoolOption opt_dontExportDirectReusedClauses(_cunstable, "reusedClauses", "Don't export directly reused clauses", false);
BoolOption opt_plingeling(_cunstable, "plingeling", "plingeling strategy for sharing clauses (exploratory feature)", false);

// added by nabesima
BoolOption   opt_fapp_imp_cla   (_cunstable, "fapp-imp-cla",    "Forced apply import clauses", false);
DoubleOption opt_fapp_prd_rate  (_cunstable, "fapp-prd-rate",   "Rate of margin after which imported clauses are applied forcedly", 5, DoubleRange(0.0, true, HUGE_VAL, false));
IntOption    opt_fapp_cla_lb    (_cunstable, "fapp-cla-lb",     "Lower-bound of unapplied imported clauses for forced application", 1000000, IntRange(0,INT32_MAX));
extern const char *_det;
IntOption    opt_margin         (_det, "margin",          "Delayed period of learnt clause exchange", 20, IntRange(0, INT32_MAX));
IntOption    opt_prd_type       (_det, "prd-type",        "the type of period (0=conflicts, 1=lit-scans, 2=blocks)", 2, IntRange(0, 2));
Int64Option  opt_prd_confs      (_det, "prd-confs",       "the number of conflicts period", 100, Int64Range(0, INT64_MAX));
Int64Option  opt_prd_lit_scans  (_det, "prd-lit-scans",   "the number of scanned literals per period", 2000000, Int64Range(0, INT64_MAX));
DoubleOption opt_prd_unit_time  (_det, "prd-unit-time",   "Unit time of blocked-base period", 0.5, DoubleRange(0.0, false, HUGE_VAL, false));
IntOption    opt_prd_stats_qsize(_det, "prd-stats-qsize", "Size of bounded queue of stats to compute block-based period", 10000, IntRange(1, INT32_MAX));
StringOption opt_dbg_log_dir    (_det, "dbg-log-dir",     "debug log directory name");
IntOption    opt_dbg_log_type   (_det, "dbg-type",        "debug log types (1=exhange, 2=speed, 4=search, 8=period)", 0);

// added by kanbara
BoolOption   opt_div_strategy   (_cunstable, "div-strategy",    "Use diversity control strategy", false);
IntOption    opt_div_periods    (_cunstable, "div-periods",     "Base periods of diversity control strategy", 100, IntRange(1, INT32_MAX));
DoubleOption opt_used_rate_wt   (_cunstable, "used-rate-wt",    "Weight of the rate of used imported clauses", 0.0, DoubleRange(0.0, true, 1.0, true));
DoubleOption opt_overlap_rate_th(_cunstable, "overlap-rate-th", "Threshold of the rate of overlapped imported clauses", 0.0, DoubleRange(0.0, true, 1.0, true));
DoubleOption opt_delta_rnd_freq (_cunstable, "delta-rnd-freq",  "Increments or decrements of random variable selection probability", 0.001, DoubleRange(0.0, true, 1.0, true));
DoubleOption opt_min_rnd_freq   (_cunstable, "min-rnd-freq",    "Minimum value of random variable selection probability", 0.00, DoubleRange(0.0, true, 1.0, true));
DoubleOption opt_max_rnd_freq   (_cunstable, "max-rnd-freq",    "Maximum value of random variable selection probability", 0.03, DoubleRange(0.0, true, 1.0, true));

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


// modified by nabesima
//static inline double cpuTime(void) {
//    struct rusage ru;
//    getrusage(RUSAGE_SELF, &ru);
//    return (double) ru.ru_utime.tv_sec + (double) ru.ru_utime.tv_usec / 1000000;
//}


void MultiSolvers::informEnd(lbool res) {
    result = res;
    pthread_cond_broadcast(&cfinished);
}


MultiSolvers::MultiSolvers(ParallelSolver *s) :
    use_simplification(true), ok(true), /* maxnbthreads(4), */ nbthreads(opt_nbsolversmultithreads), nbsolvers(opt_nbsolversmultithreads), nbcompanions(4), nbcompbysolver(2), 
    allClonesAreBuilt(false), showModel(false), winner(-1), var_decay(1 / 0.95), clause_decay(1 / 0.999), cla_inc(1), var_inc(1), random_var_freq(0.02), restart_first(100),
    restart_inc(1.5), learntsize_factor((double) 1 / (double) 3), learntsize_inc(1.1), expensive_ccmin(true), polarity_mode(polarity_false), maxmemory(opt_maxmemory),
    maxnbsolvers(opt_maxnbsolvers), ps_conf(opt_ps_conf), verb(0), verbEveryConflicts(10000), numvar(0), numclauses(0),
    // added by nabesima
    start_real_time(0.0), input_file_name(NULL), termCallbackState(NULL), termCallback(NULL)
{
    result = l_Undef;
    SharedCompanion *sc = new SharedCompanion();
    this->sharedcomp = sc;

    // Generate only solver 0.
    // It loads the formula
    // All others solvers are clone of this one
    solvers.push(s);
    s->verbosity = 0; // No reportf in solvers... All is done in MultiSolver
    s->setThreadNumber(0);
    //s->belongsto = this;
    s->sharedcomp = sc;
    sc->addSolver(s);
    assert(solvers[0]->threadNumber() == 0);

    // added by nabesima
    sc->margin = s->margin;

    pthread_mutex_init(&m, NULL);  //PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_init(&mfinished, NULL); //PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_init(&cfinished, NULL);

    if(nbsolvers > 0)
        fprintf(stdout, "c %d solvers engines and 1 companion as a blackboard created.\n", nbsolvers);
}


MultiSolvers::MultiSolvers() : MultiSolvers(new ParallelSolver(-1)) {

}


MultiSolvers::~MultiSolvers() { }


/**
 * Generate All solvers
 */

void MultiSolvers::generateAllSolvers() {
    assert(solvers[0] != NULL);
    assert(!allClonesAreBuilt);

    // Reset num of code block executions. added by nabesima
    for (int i=SearchLoopCount; i <= DetachLoopInside; i++)
        solvers[0]->stats[i] = 0;

    for(int i = 1; i < nbsolvers; i++) {
        ParallelSolver *s = (ParallelSolver *) solvers[0]->clone();
        solvers.push(s);
        s->verbosity = 0; // No reportf in solvers... All is done in MultiSolver
        s->setThreadNumber(i);
        s->sharedcomp = this->sharedcomp;
        this->sharedcomp->addSolver(s);
        assert(solvers[i]->threadNumber() == i);

        // added by nabesima for ipasir interface
        s->setTermCallback(termCallbackState, termCallback);

    /* modified by kanbara
       move to SolverConfiguration.cc::configure();
        // added by gotou in order to ensure diversity of solvers
        s->random_seed += i;
    */
    }

    // added by nabesima
    if (opt_dbg_log_dir != NULL && opt_dbg_log_type != 0) {
        std::string dirname(opt_dbg_log_dir);
        for(int i = 0; i < nbsolvers; i++) {
            std::string instance = std::string(input_file_name);
            std::string basename = instance.substr(instance.find_last_of('/') + 1);
            std::string dbg_file = dirname + "/" + basename + ".thn" + std::to_string(i) + ".csv";
            FILE *fp = fopen(dbg_file.c_str(), "wb");
            solvers[i]->dbg_log = fp;
            solvers[i]->dbg_type = opt_dbg_log_type;
        }
    }

    adjustParameters();

    allClonesAreBuilt = true;
}


/**
 * Choose solver for threads i (if no given in command line see above)
 */


ParallelSolver *MultiSolvers::retrieveSolver(int i) {
    return new ParallelSolver(i);
}


Var MultiSolvers::newVar(bool sign, bool dvar) {
    assert(solvers[0] != NULL);
    //int v;    // modified by nabesima
    sharedcomp->newVar(sign);
    if(!allClonesAreBuilt) { // At the beginning we want to generate only solvers 0
        // modified by nabesima
        //v = solvers[0]->newVar(sign, dvar);
        //assert(numvar == v + 1); // Just a useless check
#ifndef NDEBUG
        int v = solvers[0]->newVar(sign, dvar);
        assert(numvar == v); // Just a useless check
#else
        solvers[0]->newVar(sign, dvar);
#endif
    } else {
        for(int i = 0; i < nbsolvers; i++) {
            //v = solvers[i]->newVar(sign, dvar);    // modified by nabesima
            solvers[i]->newVar(sign, dvar);
        }
    }
    return numvar++;
}

// modified by nabesima
//bool MultiSolvers::addClause_(vec<Lit> &ps) {
//    assert(solvers[0] != NULL); // There is at least one solver.
//    // Check if clause is satisfied and remove false/duplicate literals:
//    if(!okay()) return false;
//
//    sort(ps);
//    Lit p;
//    int i, j;
//    for(i = j = 0, p = lit_Undef; i < ps.size(); i++)
//        if(solvers[0]->value(ps[i]) == l_True || ps[i] == ~p)
//            return true;
//        else if(solvers[0]->value(ps[i]) != l_False && ps[i] != p)
//            ps[j++] = p = ps[i];
//    ps.shrink(i - j);
//
//
//    if(ps.size() == 0) {
//        return ok = false;
//    }
//    else if(ps.size() == 1) {
//        assert(solvers[0]->value(ps[0]) == l_Undef); // TODO : Passes values to all threads
//        solvers[0]->uncheckedEnqueue(ps[0]);
//        if(!allClonesAreBuilt) {
//            return ok = ((solvers[0]->propagate()) == CRef_Undef); // checks only main solver here for propagation constradiction
//        }
//
//        // Here, all clones are built.
//        // Gives the unit clause to everybody
//        //for(int i = 0; i < nbsolvers; i++)
//        for(int i = 1; i < nbsolvers; i++)    // modified by nabesima
//            solvers[i]->uncheckedEnqueue(ps[0]);
//        return ok = ((solvers[0]->propagate()) == CRef_Undef); // checks only main solver here for propagation constradiction
//    } else {
//        //		printf("Adding clause %0xd for solver %d.\n",(void*)c, thn);
//        // At the beginning only solver 0 load the formula
//        solvers[0]->addClause(ps);
//
//        if(!allClonesAreBuilt) {
//            numclauses++;
//            return true;
//        }
//        // Clones are built, need to pass the clause to all the threads
//        for(int i = 1; i < nbsolvers; i++) {
//            solvers[i]->addClause(ps);
//        }
//        numclauses++;
//    }
//    return true;
//}
bool MultiSolvers::addClause_(vec<Lit> &ps) {

    if(!okay()) return false;

    if (allClonesAreBuilt) {
        bool ret = true;
        for (int i=0; i < nbsolvers; i++)
            ret &= solvers[i]->addClause(ps);
        return ret;
    }
    else {
        assert(getPrimarySolver() != NULL); // There is at least one solver.
        return getPrimarySolver()->addClause(ps);
    }
}

bool MultiSolvers::simplify() {
    assert(solvers[0] != NULL); // There is at least one solver.

    if(!okay()) return false;
    return ok = solvers[0]->simplify();
}


bool MultiSolvers::eliminate() {
    // TODO allow variable elimination when all threads are built!
    assert(!allClonesAreBuilt);

    SimpSolver *s = (SimpSolver *)getPrimarySolver();
    s->use_simplification = use_simplification;
    if(!use_simplification) return true;

    return s->eliminate(true);
}

// TODO: Use a template here
void *localLaunch(void *arg) {
    ParallelSolver *s = (ParallelSolver *) arg;

    // added by nabesima 
    s->incNumLiveThreads();

    // modified by nabesima
    //(void) s->solve();
    (void)s->solveLimited(s->getAssumptions(), s->getDoSimp(), s->getTurnOffSimp());

    // added by nabesima 
    s->decNumLiveThreads();

    pthread_exit(NULL);
}


// modified by nabesima
//#define MAXIMUM_SLEEP_DURATION 5


void MultiSolvers::printStats() {
    static int nbprinted = 1;
    double cpu_time = cpuTime();
    double real_time = realTime();  // added by nabesima
    printf("c\n");

    // modified by nabesima
//    printf("c |-------------------------------------------------------------------------------------------------------|\n");
//    printf("c | id | starts | decisions  |  confls    |  Init T  |  learnts | exported | imported | promoted |    %%   | \n");
//    printf("c |-------------------------------------------------------------------------------------------------------|\n");
    printf("c |------------------------------------------------------------------------------------------------------------------------|\n");
    printf("c | id | starts | decisions  |  confls    |  Init T  |  learnts | exported | imported | promoted |  period  | mrg |    %%   | \n");
    printf("c |------------------------------------------------------------------------------------------------------------------------|\n");

    //printf("%.0fs | ",cpu_time);
    for(int i = 0; i < solvers.size(); i++) {
        solvers[i]->reportProgress();
        //printf(" %2d: %12ld confl. |", i,  (long int) solvers[i]->conflicts);
    }
    long long int totalconf = 0;
    long long int totalprop = 0;
    for(int i = 0; i < solvers.size(); i++) {
        totalconf += (long int) solvers[i]->conflicts;
        totalprop += solvers[i]->propagations;
    }
    printf("c \n");

    printf("c synthesis %11lld conflicts %11lld propagations %8.0f conflicts/sec %8.0f propagations/sec\n",
           totalconf, totalprop, (double) totalconf / cpu_time, (double) totalprop / cpu_time);

    // modified by nabesima
    printf("c real: %g s, cpu: %g s, memory: %.2f Mb\n", real_time - start_real_time, cpu_time, memUsed());
    fflush(stdout);

    nbprinted++;
}


// Still a ugly function... To be rewritten with some statistics class some day
void MultiSolvers::printFinalStats() {
    sharedcomp->printStats();
    printf("c\nc\n");
    printf("c\n");
    printf("c |---------------------------------------- FINAL STATS --------------------------------------------------|\n");
    printf("c\n");

    printf("c |---------------|-----------------");
    for(int i = 0; i < solvers.size(); i++)
        printf("|------------");
    printf("|\n");

    printf("c | Threads       |      Total      ");
    for(int i = 0; i < solvers.size(); i++) {
        printf("| %10d ", i);
    }
    printf("|\n");

    printf("c |---------------|-----------------");
    for(int i = 0; i < solvers.size(); i++)
        printf("|------------");
    printf("|\n");


//--
    printf("c | Conflicts     ");
    long long int totalconf = 0;
    for(int i = 0; i < solvers.size(); i++)
        totalconf += solvers[i]->conflicts;
    printf("| %15lld ", totalconf);

    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->conflicts);
    printf("|\n");

    //--
    printf("c | Decisions     ");
    long long int totaldecs = 0;
    for(int i = 0; i < solvers.size(); i++)
        totaldecs += solvers[i]->decisions;
    printf("| %15lld ", totaldecs);

    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->decisions);
    printf("|\n");

    //--
    printf("c | Propagations  ");
    long long int totalprops = 0;
    for(int i = 0; i < solvers.size(); i++)
        totalprops += solvers[i]->propagations;
    printf("| %15lld ", totalprops);

    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->propagations);
    printf("|\n");

    // added by nabesima
    //--
    printf("c | ReduceDBs     ");
    long long int totalreducedbs = 0;
    for(int i = 0; i < solvers.size(); i++)
        totalreducedbs += solvers[i]->stats[nbReduceDB];
    printf("| %15lld ", totalreducedbs);

    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[nbReduceDB]);
    printf("|\n");

    // added by nabesima
    //--
    printf("c | SimpDBs       ");
    long long int totalsimpdbs = 0;
    for(int i = 0; i < solvers.size(); i++)
        totalsimpdbs += solvers[i]->stats[nbSimplifyDB];
    printf("| %15lld ", totalsimpdbs);

    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[nbSimplifyDB]);
    printf("|\n");

    // added by nabesima
    //--
    printf("c | Restarts      ");
    long long int totalrestarts = 0;
    for(int i = 0; i < solvers.size(); i++)
        totalrestarts += solvers[i]->starts;
    printf("| %15lld ", totalrestarts);

    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->starts);
    printf("|\n");

    printf("c | RestartBlocks ");
    long long int totalrestartblocks = 0;
    for(int i = 0; i < solvers.size(); i++)
        totalrestartblocks += solvers[i]->stats[nbstopsrestarts];
    printf("| %15lld ", totalrestartblocks);

    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[nbstopsrestarts]);
    printf("|\n");

    // added by nabesima
    //--
    printf("c | ForcedImports ");
    long long int totalforcedimports = 0;
    for(int i = 0; i < solvers.size(); i++)
        totalforcedimports += solvers[i]->stats[numForcedImports];
    printf("| %15lld ", totalforcedimports);

    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[numForcedImports]);
    printf("|\n");

    // added by gotou
    //--
    printf("c | LitFullScans  ");
    long long int totalprop_lit_full_scans = 0;
    for(int i = 0; i < solvers.size(); i++)
        totalprop_lit_full_scans += solvers[i]->num_lit_scans;
    printf("| %15lld ", totalprop_lit_full_scans);

    for(int i = 0; i < solvers.size(); i++)
        printf("|%11" PRIu64" ", solvers[i]->num_lit_scans);
    printf("|\n");

    //--
    printf("c | LitScanFracts ");
    long long int totalprop_lit_scan_fracts = 0;
    for(int i = 0; i < solvers.size(); i++)
        totalprop_lit_scan_fracts += solvers[i]->stats[sumLitScanFractinos];
    printf("| %15lld ", totalprop_lit_scan_fracts);

    for(int i = 0; i < solvers.size(); i++)
        printf("|%11" PRIu64" ", solvers[i]->stats[sumLitScanFractinos]);
    printf("|\n");

    //--
    printf("c | Periods       ");
    long long int totalperiods = 0;
    for(int i = 0; i < solvers.size(); i++)
        totalperiods += solvers[i]->periods;
    printf("| %15lld ", totalperiods);

    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->periods);
    printf("|\n");


    //--
    printf("c | Margin        ");
    long long int totalmargin = 0;
    for(int i = 0; i < solvers.size(); i++)
        totalmargin += solvers[i]->margin;
    printf("| %15lld ", totalmargin);

    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->margin);
    printf("|\n");


    printf("c | Avg_Trail     ");
    printf("|                 ");
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->conflicts==0 ? 0 : solvers[i]->stats[sumTrail] / solvers[i]->conflicts);
    printf("|\n");

    //--
    printf("c | Avg_DL        ");
    printf("|                 ");
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->conflicts==0 ? 0 : solvers[i]->stats[sumDecisionLevels] / solvers[i]->conflicts);
    printf("|\n");

    //--
    printf("c | Avg_Res       ");
    printf("|                 ");
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->conflicts==0 ? 0 : solvers[i]->stats[sumRes] / solvers[i]->conflicts);
    printf("|\n");

    //--
    printf("c | Avg_Res_Seen  ");
    printf("|                 ");
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->conflicts==0 ? 0 : solvers[i]->stats[sumResSeen] / solvers[i]->conflicts);
    printf("|\n");

    //--

    printf("c |---------------|-----------------");
    for(int i = 0; i < solvers.size(); i++)
        printf("|------------");
    printf("|\n");

    printf("c | Exported      ");
    uint64_t exported = 0;
    for(int i = 0; i < solvers.size(); i++)
        exported += solvers[i]->stats[nbexported];
    printf("| %15" PRIu64" ", exported);

    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[nbexported]);
    printf("|\n");
//--
    printf("c | Imported      ");
    uint64_t imported = 0;
    for(int i = 0; i < solvers.size(); i++)
        imported += solvers[i]->stats[nbimported];
    printf("| %15" PRIu64" ", imported);
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[nbimported]);
    printf("|\n");
//--

    printf("c | Good          ");
    uint64_t importedGood = 0;
    for(int i = 0; i < solvers.size(); i++)
        importedGood += solvers[i]->stats[nbImportedGoodClauses];
    printf("| %15" PRIu64" ", importedGood);
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[nbImportedGoodClauses]);
    printf("|\n");
//--

    printf("c | Purge         ");
    uint64_t importedPurg = 0;
    for(int i = 0; i < solvers.size(); i++)
        importedPurg += solvers[i]->stats[nbimportedInPurgatory];
    printf("| %15" PRIu64" ", importedPurg);
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[nbimportedInPurgatory]);
    printf("|\n");
//--

    printf("c | Promoted      ");
    uint64_t promoted = 0;
    for(int i = 0; i < solvers.size(); i++)
        promoted += solvers[i]->stats[nbPromoted];
    printf("| %15" PRIu64" ", promoted);
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[nbPromoted]);
    printf("|\n");
//--

    printf("c | Remove_Imp    ");
    uint64_t removedimported = 0;
    for(int i = 0; i < solvers.size(); i++)
        removedimported += solvers[i]->stats[nbRemovedUnaryWatchedClauses];
    printf("| %15" PRIu64" ", removedimported);
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[nbRemovedUnaryWatchedClauses]);
    printf("|\n");
//--

    printf("c | Blocked_Reuse ");
    uint64_t blockedreused = 0;
    for(int i = 0; i < solvers.size(); i++)
        blockedreused += solvers[i]->nbNotExportedBecauseDirectlyReused;
    printf("| %15" PRIu64" ", blockedreused);
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->nbNotExportedBecauseDirectlyReused);
    printf("|\n");
//--
    printf("c |---------------|-----------------");
    for(int i = 0; i < solvers.size(); i++)
        printf("|------------");
    printf("|\n");

    printf("c | Unaries       ");
    printf("|                 ");
    for(int i = 0; i < solvers.size(); i++) {
        printf("| %10" PRIu64" ", solvers[i]->stats[nbUn]);
    }
    printf("|\n");
//--

    printf("c | Binaries      ");
    printf("|                 ");
    for(int i = 0; i < solvers.size(); i++) {
        printf("| %10" PRIu64" ", solvers[i]->stats[nbBin]);
    }
    printf("|\n");
//--


    printf("c | Glues         ");
    printf("|                 ");
    for(int i = 0; i < solvers.size(); i++) {
        printf("| %10" PRIu64" ", solvers[i]->stats[nbDL2]);
    }
    printf("|\n");
//--

    printf("c |---------------|-----------------");
    for(int i = 0; i < solvers.size(); i++)
        printf("|------------");
    printf("|\n");

    printf("c | Orig_Seen     ");
    uint64_t origseen = 0;

    for(int i = 0; i < solvers.size(); i++) {
        origseen += solvers[i]->stats[originalClausesSeen];
    }
    printf("| %13" PRIu64" %% ", origseen * 100 / nClauses() / solvers.size());

    for(int i = 0; i < solvers.size(); i++) {
        printf("| %10" PRIu64" ", solvers[i]->stats[originalClausesSeen]);
    }

    printf("|\n");


    int winner = -1;
    for(int i = 0; i < solvers.size(); i++) {
        if(sharedcomp->winner() == solvers[i])
            winner = i;
    }

//--
    if(winner != -1) {
        printf("c | Diff Orig seen");
        printf("|                 ");

        for(int i = 0; i < solvers.size(); i++) {
            if(i == winner) {
                printf("|      X     ");
                continue;
            }
            if(solvers[i]->stats[originalClausesSeen] > solvers[winner]->stats[originalClausesSeen])
                printf("| %10" PRIu64" ", solvers[i]->stats[originalClausesSeen] - solvers[winner]->stats[originalClausesSeen]);
            else
                printf("| -%9" PRIu64" ", solvers[winner]->stats[originalClausesSeen] - solvers[i]->stats[originalClausesSeen]);

        }

        printf("|\n");
    }


//--

    if(winner != -1) {
        int sum = 0;
        printf("c | Hamming       ");
        for(int i = 0; i < solvers.size(); i++) {
            if(i == winner)
                continue;
            int nb = 0;
            for(int j = 0; j < nVars(); j++) {
                if(solvers[i]->valuePhase(j) != solvers[winner]->valuePhase(j)) nb++;
            }
            sum += nb;

        }
        sum = sum / (solvers.size() > 1 ? solvers.size() - 1 : 1);

        printf("| %13d %% ", sum * 100 / nVars());

        for(int i = 0; i < solvers.size(); i++) {
            if(i == winner) {
                printf("|      X     ");
                continue;
            }
            int nb = 0;
            for(int j = 0; j < nVars(); j++) {
                if(solvers[i]->valuePhase(j) != solvers[winner]->valuePhase(j)) nb++;
            }
            printf("| %10d ", nb);
            sum += nb;

        }
        printf("|\n");
    }

    printf("c |---------------|-----------------");
    for(int i = 0; i < solvers.size(); i++)
        printf("|------------");
    printf("|\n");

    printThreadParameters();

    // added by gotou
    printf("c\nc\n");
    printf("c [Basic stats]\n");

    // added by nabesima
    printf("c Threads : %d\n", solvers.size());

    if ( winner != -1 ) printf("c Winner : %d\n", winner);

    for(int i = 0; i < solvers.size(); i++)
        printf("c FreeVars_%d : %" PRIu64"\n", i, solvers[i]->stats[numFreeVars]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c Clauses_%d : %" PRIu64"\n", i, solvers[i]->stats[numClauses]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c Literals_%d : %" PRIu64"\n", i, solvers[i]->stats[numLiterals]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c Unaries_%d : %" PRIu64"\n", i, solvers[i]->stats[nbUn]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c Binaries_%d : %" PRIu64"\n", i, solvers[i]->stats[nbBin]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c Glues_%d : %" PRIu64"\n", i, solvers[i]->stats[nbDL2]);

    for(int i = 0; i < solvers.size(); i++)
        printf("c Conflicts_%d : %" PRIu64"\n", i, solvers[i]->conflicts);
    printf("c Conflicts_total : %lld\n", totalconf);
    for(int i = 0; i < solvers.size(); i++)
        printf("c Decisions_%d : %" PRIu64"\n", i, solvers[i]->decisions);
    printf("c Decisions_total : %lld\n", totaldecs);
    for(int i = 0; i < solvers.size(); i++)
        printf("c Propagations_%d : %" PRIu64"\n", i, solvers[i]->propagations);
    printf("c Propagations_total : %lld\n", totalprops);
    for(int i = 0; i < solvers.size(); i++)
        printf("c ReduceDBs_%d : %" PRIu64"\n", i, solvers[i]->stats[nbReduceDB]);
    printf("c ReduceDBs_total : %lld\n", totalreducedbs);
    for(int i = 0; i < solvers.size(); i++)
        printf("c SimpDBs_%d : %" PRIu64"\n", i, solvers[i]->stats[nbSimplifyDB]);
    printf("c SimpDBs_total : %lld\n", totalsimpdbs);
    for(int i = 0; i < solvers.size(); i++)
        printf("c Restarts_%d : %" PRIu64"\n", i, solvers[i]->starts);
    printf("c Restarts_total : %lld\n", totalrestarts);
    for(int i = 0; i < solvers.size(); i++)
        printf("c RestartBlockings_%d : %" PRIu64"\n", i, solvers[i]->stats[nbstopsrestarts]);

    // added by nabesima
    for(int i = 0; i < solvers.size(); i++)
        printf("c ForcedImports_%d : %" PRIu64"\n", i, solvers[i]->stats[numForcedImports]);
    printf("c ForcedImports_total : %lld\n", totalforcedimports);

    for(int i = 0; i < solvers.size(); i++)
        printf("c PropsLitFullScans_%d : %" PRIu64"\n", i, solvers[i]->num_lit_scans);
    printf("c PropsLitFullScans_total : %lld\n", totalprop_lit_full_scans);

    for(int i = 0; i < solvers.size(); i++)
        printf("c LitScanFractions_%d : %" PRIu64"\n", i, solvers[i]->stats[sumLitScanFractinos]);
    printf("c LitScanFractions_total : %lld\n", totalprop_lit_scan_fracts);

    for(int i = 0; i < solvers.size(); i++)
        printf("c PropagationLitScans_%d : %" PRIu64"\n", i, solvers[i]->stats[sumPropagationLitScans]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c AnalysisLitScans_%d : %" PRIu64"\n", i, solvers[i]->stats[sumAnalysisLitScans]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c ReductionLitScans_%d : %" PRIu64"\n", i, solvers[i]->stats[sumReductionLitScans]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c SimplificationLitScans_%d : %" PRIu64"\n", i, solvers[i]->stats[sumSimplificationLitScans]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c SharingLitScans_%d : %" PRIu64"\n", i, solvers[i]->stats[sumSharingLitScans]);

    for(int i = 0; i < solvers.size(); i++)
        printf("c Periods_%d : %" PRIu64"\n", i, solvers[i]->periods);
    printf("c Periods_total : %lld\n", totalperiods);

    for(int i = 0; i < solvers.size(); i++)
        printf("c Margin_%d : %" PRIu64"\n", i, solvers[i]->margin);
    printf("c Margin_total : %lld\n", totalmargin);

    for(int i = 0; i < solvers.size(); i++)
        printf("c Exported_%d : %" PRIu64"\n", i, solvers[i]->stats[nbexported]);
    printf("c Exported_total : %" PRIu64"\n", exported);

    for(int i = 0; i < solvers.size(); i++)
        printf("c Imported_%d : %" PRIu64"\n", i, solvers[i]->stats[nbimported]);
    printf("c Imported_total : %" PRIu64"\n", imported);

    for(int i = 0; i < solvers.size(); i++)
        printf("c Good_%d : %" PRIu64"\n", i, solvers[i]->stats[nbImportedGoodClauses]);
    printf("c Good_total : %" PRIu64"\n", importedGood);

    for(int i = 0; i < solvers.size(); i++)
        printf("c Purge_%d : %" PRIu64"\n", i, solvers[i]->stats[nbimportedInPurgatory]);
    printf("c Purge_total : %" PRIu64"\n", importedPurg);

    for(int i = 0; i < solvers.size(); i++)
        printf("c Promoted_%d : %" PRIu64"\n", i, solvers[i]->stats[nbPromoted]);
    printf("c Promoted_total : %" PRIu64"\n", promoted);

    for(int i = 0; i < solvers.size(); i++)
        printf("c Remove_Imp_%d : %" PRIu64"\n", i, solvers[i]->stats[nbRemovedUnaryWatchedClauses]);
    printf("c Remove_Imp_total : %" PRIu64"\n", removedimported);

    for(int i = 0; i < solvers.size(); i++)
        printf("c Blocked_Reuse_%d : %" PRIu64"\n", i, solvers[i]->nbNotExportedBecauseDirectlyReused);
    printf("c Blocked_Reuse_total : %" PRIu64"\n", blockedreused);

    // added by nabesima
    printf("c\n");
    printf("c [Executed number of each blocks]\n");

    // [search]
    for(int i = 0; i < solvers.size(); i++)
        printf("c SearchLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[SearchLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c SearchLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[SearchLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c SearchConflictCount_%d : %" PRIu64"\n", i, solvers[i]->stats[SearchConflictCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c SearchDecisionCount_%d : %" PRIu64"\n", i, solvers[i]->stats[SearchDecisionCount]);
    // [propagation]
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropCleanWatchLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[PropCleanWatchLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropCleanWatchLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[PropCleanWatchLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropQueueLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[PropQueueLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropBinWatchLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[PropBinWatchLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropBinWatchLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[PropBinWatchLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropBinWatchLoopInsideConf_%d : %" PRIu64"\n", i, solvers[i]->stats[PropBinWatchLoopInsideConf]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropBinWatchLoopInsideEnq_%d : %" PRIu64"\n", i, solvers[i]->stats[PropBinWatchLoopInsideEnq]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropBinWatchLoopInsideNext_%d : %" PRIu64"\n", i, solvers[i]->stats[PropBinWatchLoopInsideNext]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropWatchLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[PropWatchLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropWatchLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[PropWatchLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropWatchLoopInsideNext_%d : %" PRIu64"\n", i, solvers[i]->stats[PropWatchLoopInsideNext]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropWatchLoopInsideBlocker_%d : %" PRIu64"\n", i, solvers[i]->stats[PropWatchLoopInsideBlocker]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropClauseLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[PropClauseLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropClauseLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[PropClauseLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropClauseLoopBreakCount_%d : %" PRIu64"\n", i, solvers[i]->stats[PropClauseLoopBreakCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c UnaryPropWatchLoopCount_%d : %" PRIu64 "\n", i, solvers[i]->stats[UnaryPropWatchLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c UnaryPropWatchLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[UnaryPropWatchLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c UnaryPropWatchLoopInsideBlocker_%d : %" PRIu64"\n", i, solvers[i]->stats[UnaryPropWatchLoopInsideBlocker]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c UnaryPropClauseLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[UnaryPropClauseLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c UnaryPropClauseLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[UnaryPropClauseLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c UnaryPropClauseLoopBreakCount_%d : %" PRIu64"\n", i, solvers[i]->stats[UnaryPropClauseLoopBreakCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c UnaryPropConfClauseLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[UnaryPropConfClauseLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c UnaryPropConfClauseLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[UnaryPropConfClauseLoopInside]);
    // [analysis]
    for(int i = 0; i < solvers.size(); i++)
        printf("c AnalyzeLoopCount_%d : %" PRIu64 "\n", i, solvers[i]->stats[AnalyzeLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c AnalyzeClauseLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[AnalyzeClauseLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c AnalyzeClauseLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[AnalyzeClauseLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c AnalyzeOutMinLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[AnalyzeOutMinLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c AnalyzeOutMinLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[AnalyzeOutMinLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c AnalyzeOutLBDLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[AnalyzeOutLBDLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c AnalyzeOutLBDLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[AnalyzeOutLBDLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c AnalyzeLastDecLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[AnalyzeLastDecLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c AnalyzeLastDecLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[AnalyzeLastDecLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c LitRedLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[LitRedLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c LitRedClauseLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[LitRedClauseLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c LitRedClauseLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[LitRedClauseLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c MinBinResCount_%d : %" PRIu64"\n", i, solvers[i]->stats[MinBinResCount]);
    // [reduction]
    for(int i = 0; i < solvers.size(); i++)
        printf("c RedLearntCompCount_%d : %" PRIu64"\n", i, solvers[i]->stats[RedLearntCompCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c RedLearntCompInside_%d : %" PRIu64"\n", i, solvers[i]->stats[RedLearntCompInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c RedLearntLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[RedLearntLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c RedLearntLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[RedLearntLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c RedUWLearntCompCount_%d : %" PRIu64"\n", i, solvers[i]->stats[RedUWLearntCompCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c RedUWLearntCompInside_%d : %" PRIu64"\n", i, solvers[i]->stats[RedUWLearntCompInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c RedUWLearntLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[RedUWLearntLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c RedUWLearntLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[RedUWLearntLoopInside]);
    // [simplification]
    for(int i = 0; i < solvers.size(); i++)
        printf("c SimpRmSatClausesLoopCount_%d : %" PRIu64 "\n", i, solvers[i]->stats[SimpRmSatClausesLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c SimpClauseLoopCount_%d : %" PRIu64 "\n", i, solvers[i]->stats[SimpClauseLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c SimpClauseLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[SimpClauseLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c SimpVarHeapLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[SimpVarHeapLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c SimpVarHeapLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[SimpVarHeapLoopInside]);
    // [sharing]
    for(int i = 0; i < solvers.size(); i++)
        printf("c ExpUnaryClauseCount_%d : %" PRIu64"\n", i, solvers[i]->stats[ExpUnaryClauseCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c ExpClauseDuringSearchCount_%d : %" PRIu64"\n", i, solvers[i]->stats[ExpClauseDuringSearchCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c ExpClauseDuringSearchInside_%d : %" PRIu64"\n", i, solvers[i]->stats[ExpClauseDuringSearchInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c ExpClauseDuringAnalysisCount_%d : %" PRIu64"\n", i, solvers[i]->stats[ExpClauseDuringAnalysisCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c ExpClauseDuringAnalysisInside_%d : %" PRIu64"\n", i, solvers[i]->stats[ExpClauseDuringAnalysisInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c ImpThreadLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[ImpThreadLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c ImpQueueLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[ImpQueueLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c ImpClauseLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[ImpClauseLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c ImpClauseLitLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[ImpClauseLitLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c ImpClauseLitLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[ImpClauseLitLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c ApplyImportedClausesLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[ApplyImportedClausesLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c ApplyImportedClausesLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[ApplyImportedClausesLoopInside]);
    // [reloc]
    for(int i = 0; i < solvers.size(); i++)
        printf("c RelocCleanWatchLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[RelocCleanWatchLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c RelocCleanWatchLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[RelocCleanWatchLoopInside]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c RelocLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[RelocLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c RelocLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[RelocLoopInside]);
    // [detach]
    for(int i = 0; i < solvers.size(); i++)
        printf("c DetachLoopCount_%d : %" PRIu64"\n", i, solvers[i]->stats[DetachLoopCount]);
    for(int i = 0; i < solvers.size(); i++)
        printf("c DetachLoopInside_%d : %" PRIu64"\n", i, solvers[i]->stats[DetachLoopInside]);

    printf("c\n");
#ifdef MEASURE_TIME
    printf("c [Running time for each blocks]\n");

    // [search]
    for(int i = 0; i < solvers.size(); i++)
        printf("c SearchLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(SearchLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c SearchConflictTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(SearchConflict));
    for(int i = 0; i < solvers.size(); i++)
        printf("c SearchDecisionTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(SearchDecision));
    // [propagation]
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropCleanWatchLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(PropCleanWatchLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropQueueLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(PropQueueLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropBinWatchLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(PropBinWatchLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropWatchLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(PropWatchLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropClauseLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(PropClauseLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c PropClauseLoopBreakTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(PropClauseLoopBreak));
    for(int i = 0; i < solvers.size(); i++)
        printf("c UnaryPropWatchLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(UnaryPropWatchLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c UnaryPropClauseLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(UnaryPropClauseLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c UnaryPropClauseLoopBreakTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(UnaryPropClauseLoopBreak));
    for(int i = 0; i < solvers.size(); i++)
        printf("c UnaryPropConfClauseLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(UnaryPropConfClauseLoop));
    // [analysis]
    for(int i = 0; i < solvers.size(); i++)
        printf("c AnalyzeLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(AnalyzeLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c AnalyzeClauseLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(AnalyzeClauseLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c AnalyzeOutMinLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(AnalyzeOutMinLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c AnalyzeOutLBDLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(AnalyzeOutLBDLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c AnalyzeLastDecLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(AnalyzeLastDecLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c LitRedLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(LitRedLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c LitRedClauseLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(LitRedClauseLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c MinBinResTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(MinBinRes));
    // [reduction]
    for(int i = 0; i < solvers.size(); i++)
        printf("c RedLearntCompTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(RedLearntComp));
    for(int i = 0; i < solvers.size(); i++)
        printf("c RedLearntLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(RedLearntLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c RedUWLearntCompTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(RedUWLearntComp));
    for(int i = 0; i < solvers.size(); i++)
        printf("c RedUWLearntLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(RedUWLearntLoop));
    // [simplification]
    for(int i = 0; i < solvers.size(); i++)
        printf("c SimpClauseLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(SimpClauseLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c SimpRmSatClausesLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(SimpRmSatClausesLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c SimpVarHeapLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(SimpVarHeapLoop));
    // [sharing]
    for(int i = 0; i < solvers.size(); i++)
        printf("c ExpUnaryClauseTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(ExpUnaryClause));
    for(int i = 0; i < solvers.size(); i++)
        printf("c ExpClauseDuringSearchTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(ExpClauseDuringSearch));
    for(int i = 0; i < solvers.size(); i++)
        printf("c ExpClauseDuringAnalysisTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(ExpClauseDuringAnalysis));
    for(int i = 0; i < solvers.size(); i++)
        printf("c ImpThreadLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(ImpThreadLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c ImpQueueLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(ImpQueueLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c ImpClauseLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(ImpClauseLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c ImpClauseLitLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(ImpClauseLitLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c ApplyImportedClausesLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(ApplyImportedClausesLoop));
    // [reloc]
    for(int i = 0; i < solvers.size(); i++)
        printf("c RelocCleanWatchLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(RelocCleanWatchLoop));
    for(int i = 0; i < solvers.size(); i++)
        printf("c RelocLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(RelocLoop));
    // [detach]
    for(int i = 0; i < solvers.size(); i++)
        printf("c DetachLoopTime_%d : %f\n", i, solvers[i]->seqchrono.getTime(DetachLoop));
    // ALL
    for(int i = 0; i < solvers.size(); i++) {
        double total = 0.0;
        for (ProcName j=SearchLoop; j < NumSeqProcTypes; j = static_cast<ProcName>(j + 1))
            total += solvers[i]->seqchrono.getTime(j);
        printf("c CumBlockRunningTime_%d : %f\n", i, total);
    }
#endif

    // added by kanbara
    if (solvers[0]->div_strategy) {
        for(int i = 0; i < solvers.size(); i++)
            printf("c UsedRateTime_%d : %f\n", i, solvers[i]->usedrate_time.getSumTime());
        for(int i = 0; i < solvers.size(); i++)
            printf("c OverlapRateTime_%d : %f\n", i, solvers[i]->overlaprate_time.getSumTime());

        for(int i = 0; i < solvers.size(); i++) {
          for(int j = 0; j < solvers[i]->learnts_used_rate.size(); j++) {
              printf("c LearntsUsedRate_%d_%d : %f\n", i, j, solvers[i]->learnts_used_rate[j]);
          }
        }

        for(int i = 0; i < solvers.size(); i++) {
          for(int j = 0; j < solvers[i]->import_used_rate.size(); j++) {
              printf("c ImportUsedRate_%d_%d : %f\n", i, j, solvers[i]->import_used_rate[j]);
          }
        }

        for(int i = 0; i < solvers.size(); i++) {
          for(int j = 0; j < solvers[i]->overlap_rate.size(); j++) {
              printf("c OverlapRate_%d_%d : %f\n", i, j, solvers[i]->overlap_rate[j]);
          }
        }

        for(int i = 0; i < solvers.size(); i++) {
          for(int j = 0; j < solvers[i]->rnd_freq.size(); j++) {
              printf("c random_var_freq_%d_%d : %f\n", i, j, solvers[i]->rnd_freq[j]);
          }
        }
    }

    printf("c\n");
    if (solvers[0]->prd_type == prd_type_blocks)
        printPrdCoeffStats();

    printf("c\n");
    printf("c [Time stats for parallel proccessing]\n");

    for(int i = 0; i < solvers.size(); i++)
           printf("c RunningTime_%d : %f\n", i, solvers[i]->parchrono.getTime(RunningTime));
    double runnings = 0;
    for(int i = 0; i < solvers.size(); i++)
        runnings += solvers[i]->parchrono.getTime(RunningTime);
    printf("c RunningTime_total : %f\n", runnings);

    for(int i = 0; i < solvers.size(); i++)
        printf("c WaitingTime_%d : %f\n", i, solvers[i]->parchrono.getTime(WaitingTime));
    double waitings = 0;
    for(int i = 0; i < solvers.size(); i++)
        waitings += solvers[i]->parchrono.getTime(WaitingTime);
    printf("c WaitingTime_total : %f\n", waitings);

    for(int i = 0; i < solvers.size(); i++)
        printf("c ExchangingTime_%d : %f\n", i, solvers[i]->parchrono.getTime(ExchangingTime));
    double exchangings = 0;
    for(int i = 0; i < solvers.size(); i++)
        exchangings += solvers[i]->parchrono.getTime(ExchangingTime);
    printf("c ExchangingTime_total : %f\n", exchangings);

    for(int i = 0; i < solvers.size(); i++)
        printf("c PeriodProcTime_%d : %f\n", i, solvers[i]->parchrono.getTime(PeriodUpdateTime));
    double nextperiods = 0;
    for(int i = 0; i < solvers.size(); i++)
        nextperiods += solvers[i]->parchrono.getTime(PeriodUpdateTime);
    printf("c PeriodProcTime_total : %f\n", nextperiods);

    vec<double> sum_non_waiting_time(solvers.size());
    vec<double> sum_solve_time(solvers.size());
    for(int i = 0; i < solvers.size(); i++) {
        sum_non_waiting_time[i] =
                  solvers[i]->parchrono.getTime(RunningTime)
                + solvers[i]->parchrono.getTime(ExchangingTime)
                + solvers[i]->parchrono.getTime(PeriodUpdateTime);
        sum_solve_time[i] =
                  sum_non_waiting_time[i]
                + solvers[i]->parchrono.getTime(WaitingTime);
    }

    for(int i = 0; i < solvers.size(); i++)
        printf("c NonWaitingTime_%d : %f\n", i, sum_non_waiting_time[i]);
    double tot_non_waiting = 0;
    for(int i = 0; i < solvers.size(); i++)
        tot_non_waiting += sum_non_waiting_time[i];
    printf("c NonWaitingTime_total : %f\n", tot_non_waiting);
    for(int i = 0; i < solvers.size(); i++)
        printf("c SolvingTime_%d : %f\n", i, sum_solve_time[i]);
    double tot_solvings = 0;
    for(int i = 0; i < solvers.size(); i++)
        tot_solvings += sum_solve_time[i];
    printf("c SolvingTime_total : %f\n", tot_solvings);
}

void MultiSolvers::printThreadParameters() {
    printf("c\n");
    printf("c |------------------------------------------ THREAD PARAMETERS ------------------------------------------|\n");
    printf("c\n");
    printf("c |-----------------------");
    for (int i = 0; i < solvers.size(); i++)
        printf("|------------");
    printf("|\n");
    printf("c | Threads               ");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10d ", i);
    printf("|\n");
    printf("c |-----------------------");
    for (int i = 0; i < solvers.size(); i++)
        printf("|------------");
    printf("|\n");
    printf("c | randomizeFirstDescent ");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10s ", solvers[i]->randomizeFirstDescent ? "True" : "False");
    printf("|\n");
    printf("c | adaptStrategeis       ");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10s ", solvers[i]->adaptStrategies ? "True" : "False");
    printf("|\n");
    printf("c | forceUnsatOnNewDescent");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10s ", solvers[i]->forceUnsatOnNewDescent ? "True" : "False");
    printf("|\n");
    printf("c | Initial random_seed   ");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10.0f ", solvers[i]->random_seed);
    printf("|\n");
    if (solvers.size() > 8) {
        printf("c | goodlimitlbd          ");
        for (int i = 0; i < solvers.size(); i++)
            printf("| %10d ", solvers[i]->goodlimitlbd);
        printf("|\n");
        printf("c | goodlimitsize         ");
        for (int i = 0; i < solvers.size(); i++)
            printf("| %10d ", solvers[i]->goodlimitsize);
        printf("|\n");
    }

    printf("c | chanseokStrategy      ");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10s ", solvers[i]->chanseokStrategy ? "True" : "False");
    printf("|\n");
    printf("c | coLBDBound            ");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10d ", solvers[i]->coLBDBound);
    printf("|\n");
    printf("c | glureduce             ");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10s ", solvers[i]->glureduce ? "True" : "False");
    printf("|\n");
    printf("c | firstReduceDB         ");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10d ", solvers[i]->firstReduceDB);
    printf("|\n");
    printf("c | incReduceDB           ");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10d ", solvers[i]->incReduceDB);
    printf("|\n");
    printf("c | luby_restart          ");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10s ", solvers[i]->luby_restart ? "True" : "False");
    printf("|\n");
    printf("c | luby_restart_factor   ");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10d ", solvers[i]->luby_restart_factor);
    printf("|\n");
    printf("c | var_decay             ");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10.4f ", solvers[i]->var_decay);
    printf("|\n");
    printf("c | max_var_decay         ");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10.4f ", solvers[i]->max_var_decay);
    printf("|\n");
    printf("c | randomize_on_restarts ");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10s ", solvers[i]->randomize_on_restarts ? "True" : "False");
    printf("|\n");
    printf("c | remove_satisfied      ");
    for (int i = 0; i < solvers.size(); i++)
        printf("| %10s ", solvers[i]->remove_satisfied ? "True" : "False");
    printf("|\n");


    printf("c |-----------------------");
    for (int i = 0; i < solvers.size(); i++)
        printf("|------------");
    printf("|\n");
}

void MultiSolvers::printPrdCoeffStats() {

    std::map<int, std::string> stat2str;
    // search
    stat2str[SearchLoopCount] = "SearchLoopCount";
    stat2str[SearchLoopInside] = "SearchLoopInside";
    stat2str[SearchConflictCount] = "SearchConflictCount";
    stat2str[SearchDecisionCount] = "SearchDecisionCount";
    // propagation
    stat2str[PropCleanWatchLoopCount] = "PropCleanWatchLoopCount";
    stat2str[PropCleanWatchLoopInside] = "PropCleanWatchLoopInside";
    stat2str[PropQueueLoopCount] = "PropQueueLoopCount";
    stat2str[PropBinWatchLoopCount] = "PropBinWatchLoopCount";
    stat2str[PropBinWatchLoopInside] = "PropBinWatchLoopInside";
    stat2str[PropBinWatchLoopInsideConf] = "PropBinWatchLoopInsideConf";
    stat2str[PropBinWatchLoopInsideEnq] = "PropBinWatchLoopInsideEnq";
    stat2str[PropBinWatchLoopInsideNext] = "PropBinWatchLoopInsideNext";
    stat2str[PropWatchLoopCount] = "PropWatchLoopCount";
    stat2str[PropWatchLoopInside] = "PropWatchLoopInside";
    stat2str[PropWatchLoopInsideNext] = "PropWatchLoopInsideNext";
    stat2str[PropWatchLoopInsideBlocker] = "PropWatchLoopInsideBlocker";
    stat2str[PropClauseLoopCount] = "PropClauseLoopCount";
    stat2str[PropClauseLoopInside] = "PropClauseLoopInside";
    stat2str[PropClauseLoopBreakCount] = "PropClauseLoopBreakCount";
    stat2str[UnaryPropWatchLoopCount] = "UnaryPropWatchLoopCount";
    stat2str[UnaryPropWatchLoopInside] = "UnaryPropWatchLoopInside";
    stat2str[UnaryPropWatchLoopInsideBlocker] = "UnaryPropWatchLoopInsideBlocker";
    stat2str[UnaryPropClauseLoopCount] = "UnaryPropClauseLoopCount";
    stat2str[UnaryPropClauseLoopInside] = "UnaryPropClauseLoopInside";
    stat2str[UnaryPropClauseLoopBreakCount] = "UnaryPropClauseLoopBreakCount";
    stat2str[UnaryPropConfClauseLoopCount] = "UnaryPropConfClauseLoopCount";
    stat2str[UnaryPropConfClauseLoopInside] = "UnaryPropConfClauseLoopInside";
    // conflict analysis
    stat2str[AnalyzeLoopCount] = "AnalyzeLoopCount";
    stat2str[AnalyzeClauseLoopCount] = "AnalyzeClauseLoopCount";
    stat2str[AnalyzeClauseLoopInside] = "AnalyzeClauseLoopInside";
    stat2str[AnalyzeOutMinLoopCount] = "AnalyzeOutMinLoopCount";
    stat2str[AnalyzeOutMinLoopInside] = "AnalyzeOutMinLoopInside";
    stat2str[AnalyzeOutLBDLoopCount] = "AnalyzeOutLBDLoopCount";
    stat2str[AnalyzeOutLBDLoopInside] = "AnalyzeOutLBDLoopInside";
    stat2str[AnalyzeLastDecLoopCount] = "AnalyzeLastDecLoopCount";
    stat2str[AnalyzeLastDecLoopInside] = "AnalyzeLastDecLoopInside";
    stat2str[LitRedLoopCount] = "LitRedLoopCount";
    stat2str[LitRedClauseLoopCount] = "LitRedClauseLoopCount";
    stat2str[LitRedClauseLoopInside] = "LitRedClauseLoopInside";
    stat2str[MinBinResCount] = "MinBinResCount";
    // clause reduction
    stat2str[RedLearntCompCount] = "RedLearntCompCount";
    stat2str[RedLearntCompInside] = "RedLearntCompInside";
    stat2str[RedLearntLoopCount] = "RedLearntLoopCount";
    stat2str[RedLearntLoopInside] = "RedLearntLoopInside";
    stat2str[RedUWLearntCompCount] = "RedUWLearntCompCount";
    stat2str[RedUWLearntCompInside] = "RedUWLearntCompInside";
    stat2str[RedUWLearntLoopCount] = "RedUWLearntLoopCount";
    stat2str[RedUWLearntLoopInside] = "RedUWLearntLoopInside";
    // simplification
    stat2str[SimpRmSatClausesLoopCount] = "SimpRmSatClausesLoopCount";
    stat2str[SimpClauseLoopCount] = "SimpClauseLoopCount";
    stat2str[SimpClauseLoopInside] = "SimpClauseLoopInside";
    stat2str[SimpVarHeapLoopCount] = "SimpVarHeapLoopCount";
    stat2str[SimpVarHeapLoopInside] = "SimpVarHeapLoopInside";
    // sharing
    stat2str[ExpUnaryClauseCount] = "ExpUnaryClauseCount";
    stat2str[ExpClauseDuringSearchCount] = "ExpClauseDuringSearchCount";
    stat2str[ExpClauseDuringSearchInside] = "ExpClauseDuringSearchInside";
    stat2str[ExpClauseDuringAnalysisCount] = "ExpClauseDuringAnalysisCount";
    stat2str[ExpClauseDuringAnalysisInside] = "ExpClauseDuringAnalysisInside";
    stat2str[ImpThreadLoopCount] = "ImpThreadLoopCount";
    stat2str[ImpQueueLoopCount] = "ImpQueueLoopCount";
    stat2str[ImpClauseLoopCount] = "ImpClauseLoopCount";
    stat2str[ImpClauseLitLoopCount] = "ImpClauseLitLoopCount";
    stat2str[ImpClauseLitLoopInside] = "ImpClauseLitLoopInside";
    stat2str[ApplyImportedClausesLoopCount] = "ApplyImportedClausesLoopCount";
    stat2str[ApplyImportedClausesLoopInside] = "ApplyImportedClausesLoopInside";
    // reloc
    stat2str[RelocCleanWatchLoopCount] = "RelocCleanWatchLoopCount";
    stat2str[RelocCleanWatchLoopInside] = "RelocCleanWatchLoopInside";
    stat2str[RelocLoopCount] = "RelocLoopCount";
    stat2str[RelocLoopInside] = "RelocLoopInside";
    // detach
    stat2str[DetachLoopCount] = "DetachLoopCount";
    stat2str[DetachLoopInside] = "DetachLoopInside";

    printf("c [Estimated ratio of running time]\n");
    for(int i=0; i < solvers.size(); i++) {
        ParallelSolver &s = *solvers[i];

        double total = 0.0;
        for (int j=0; j < s.prd_coeff.size(); j++)
            total += s.stats[s.prd_coeff[j].index] * s.prd_coeff[j].value;

        for (int j=0; j < s.prd_coeff.size(); j++) {
            double time = s.stats[s.prd_coeff[j].index] * s.prd_coeff[j].value;
            double ratio = time / total;
            printf("c %sRatio_%d : %f\n", stat2str[s.prd_coeff[j].index].c_str(), s.thn, ratio);
        }
    }
}

// Well, all those parameteres are just naive guesses... No experimental evidences for this.
void MultiSolvers::adjustParameters() {
    SolverConfiguration::configure(this, nbsolvers);
}


void MultiSolvers::adjustNumberOfCores() {
    // added by nabesima
    if (maxmemory == 0)
        maxmemory = memSize();

    double mem = memUsed();
    if(nbthreads == 0) { // Automatic configuration
        if(verb >= 1)
            printf("c |  Automatic Adjustement of the number of solvers. MaxMemory=%5d, MaxCores=%3d.                       |\n", maxmemory, maxnbsolvers);
        unsigned int tmpnbsolvers = maxmemory * 4 / 10 / mem;
        if(tmpnbsolvers > maxnbsolvers) tmpnbsolvers = maxnbsolvers;
        if(tmpnbsolvers < 1) tmpnbsolvers = 1;
        if(verb >= 1)
            printf("c |  One Solver is taking %9.2fMb... Let's take %3d solvers for this run (max 40%% of the maxmemory).  |\n", mem, tmpnbsolvers);
        nbsolvers = tmpnbsolvers;
        nbthreads = nbsolvers;
    }
    // added by nabesima
    int pos_threads = std::max((int)(maxmemory * 4 / 10 / mem), 1);
    if (nbthreads > pos_threads) {
        nbthreads = nbsolvers = pos_threads;
        if (verb >= 1)
            printf("c |  One Solver is taking %9.2fMb... Let's take %3d solvers for this run (max 40%% of the maxmemory).  |\n", mem, pos_threads);
    }

    assert(nbthreads == nbsolvers);
}


lbool MultiSolvers::solve_(bool do_simp, bool turn_off_simp) {
    pthread_attr_t thAttr;
    int i;

    if (!allClonesAreBuilt) {    // added by nabesima
        adjustNumberOfCores();
        sharedcomp->setNbThreads(nbsolvers);
        if(verb >= 1) {
            printf("c |  Generating clones                                                                                    |\n");
            fflush(stdout);    // added by nabesima
        }
        generateAllSolvers();
        if(verb >= 1) {
            printf("c |  all clones generated. Memory = %8.2fMb.                                                           |\n", memUsed());
            printf("c =========================================================================================================\n");

            printThreadParameters(); // added by kanbara
            fflush(stdout);    // added by nabesima
        }
    }
    else { // added by nabesima
        if(verb >= 1) {
            printf("c =========================================================================================================\n");
            fflush(stdout);
        }
    }

    model.clear();
    conflict.clear();

    /* Initialize and set thread detached attribute */
    pthread_attr_init(&thAttr);
    pthread_attr_setdetachstate(&thAttr, PTHREAD_CREATE_JOINABLE);

    sharedcomp->init();    // added by nabesima for incremental SAT solving
    threads.clear();       // added by nabesima for incremental SAT solving

//    // DEBUG
//    for(i = 0; i < nbsolvers; i++)
//        printf("STARTS thread %d period %lu, prdClausesQueue last period %lu\n", i, solvers[i]->periods, sharedcomp->getPrdClausesQueue(i).last().period());

    // added by nabesima
    for(i = 0; i < nbsolvers; i++)
        solvers[i]->periods = sharedcomp->getPrdClausesQueue(i).last().period();

    // modified by nabesima to avoid a deadlock when incremental SAT solving
    // // Launching all solvers
    // for(i = 0; i < nbsolvers; i++) {
    //     pthread_t *pt = (pthread_t *) malloc(sizeof(pthread_t));
    //     threads.push(pt);
    //     solvers[i]->pmfinished = &mfinished;
    //     solvers[i]->pcfinished = &cfinished;
    //     solvers[i]->setAssumptions(assumptions);
    //     solvers[i]->setDoSimp(do_simp);
    //     solvers[i]->setTurnOffSimp(turn_off_simp);
    //     pthread_create(threads[i], &thAttr, &localLaunch, (void *) solvers[i]);
    // }

    bool launched = false;
    bool done = false;
    //bool adjustedlimitonce = false;    // modified by nabesima for determinisity

    (void) pthread_mutex_lock(&m);
    while(!done) {

        // added by nabesima
        if (!launched) {
            // Launching all solvers
            for(i = 0; i < nbsolvers; i++) {
                pthread_t *pt = (pthread_t *) malloc(sizeof(pthread_t));
                threads.push(pt);
                solvers[i]->pmfinished = &mfinished;
                solvers[i]->pcfinished = &cfinished;
                solvers[i]->setAssumptions(assumptions);
                solvers[i]->setDoSimp(do_simp);
                solvers[i]->setTurnOffSimp(turn_off_simp);
                pthread_create(threads[i], &thAttr, &localLaunch, (void *) solvers[i]);
            }
            launched = true;
        }
        // added by nabesima
        if (sharedcomp->hasNoLiveThreads())
            break;

        struct timespec timeout;
        time(&timeout.tv_sec);
        //timeout.tv_sec += MAXIMUM_SLEEP_DURATION;  // modified by nabesima
        timeout.tv_sec += opt_statsInterval;
        timeout.tv_nsec = 0;
        if(pthread_cond_timedwait(&cfinished, &mfinished, &timeout) != ETIMEDOUT)
            done = true;
        //else
        else if (verb >= 1)    // modified by nabesima
            printStats();

#ifdef NONDETERMINISTIC    // added by nabesima
        float mem = memUsed();
        if((maxmemory > 0) && (mem > maxmemory) && !sharedcomp->panicMode)
            printf("c ** reduceDB switching to Panic Mode due to memory limitations !\n"), sharedcomp->panicMode = true;

        if(!done && !adjustedlimitonce) {
            uint64_t sumconf = 0;
            uint64_t sumimported = 0;
            for(int i = 0; i < nbsolvers; i++) {
                sumconf += solvers[i]->conflicts;
                sumimported += solvers[i]->stats[nbimported];
            }
            if(sumconf > 10000000 && sumimported > 4 * sumconf) { // too many many imported clauses (after a while)
                for(int i = 0; i < nbsolvers; i++) { // we have like 32 threads, so we need to export just very good clauses
                    solvers[i]->goodlimitlbd -= 2;
                    solvers[i]->goodlimitsize -= 4;
                }
                adjustedlimitonce = true;
                printf("c adjusting (once) the limits to send fewer clauses.\n");
            }
        }
#endif
    }

    (void) pthread_mutex_unlock(&m);

    for(i = 0; i < nbsolvers; i++) { // Wait for all threads to finish
        pthread_join(*threads[i], NULL);
        free(threads[i]);        // added by nabesima
    }

    assert(sharedcomp != NULL);
    result = sharedcomp->jobStatus;
    if(result == l_True) {
        sharedcomp->jobFinishedBy->extendModel();
        int n = sharedcomp->jobFinishedBy->nVars();
        model.growTo(n);
        for(int i = 0; i < n; i++) {
            model[i] = sharedcomp->jobFinishedBy->model[i];
            assert(model[i] != l_Undef);
        }
    }
    // added by nabesima
    else if (result == l_False) {
        sharedcomp->jobFinishedBy->conflict.copyTo(conflict);
    }

    // added by nabesima
    for(i = 0; i < nbsolvers; i++)
        solvers[i]->cancelUntil(0);

//    // DEBUG
//    for(i = 0; i < nbsolvers; i++)
//        printf("FINISHED thread %d period %lu, prdClausesQueue last period %lu\n", i, solvers[i]->periods, sharedcomp->getPrdClausesQueue(i).last().period());


    return result;
    /*
      for(int i=0;i<NBTHREADS;i++)
      pthread_join(*threads[i],&status);
    */

}

// added by gotou
void MultiSolvers::processInterruptedTermination() {
    for (int i = 0; i < solvers.size(); ++i) {
        solvers[i]->parchrono.stopAll();
        // added by nabesima
        if (solvers[i]->dbg_log) { fflush(solvers[i]->dbg_log); fclose(solvers[i]->dbg_log); }
    }
}

// added by nabesima (for naps interface)
void MultiSolvers::setFrozen(Var v, bool b) {
    if (allClonesAreBuilt) {
        for(int i = 0; i < nbsolvers; i++) {
            //v = solvers[i]->newVar(sign, dvar);    // modified by nabesima
            solvers[i]->setFrozen(v, b);
        }
    }
    else
        getPrimarySolver()->setFrozen(v, b);
}
bool MultiSolvers::isEliminated(Var v) const {
    if (allClonesAreBuilt) {
        bool ret = false;
        for(int i = 0; i < nbsolvers; i++)
            ret |= solvers[i]->isEliminated(v);    // a variable is eliminated if at least one thread eliminates it.
        return ret;
    }
    else {
        assert(getPrimarySolver() != NULL); // There is at least one solver.
        return getPrimarySolver()->isEliminated(v);
    }
}
void MultiSolvers::toDimacs(const char* file, int tid) {
    if (allClonesAreBuilt) {
        assert(tid < nbsolvers);
        assert(solvers[tid] != NULL);
        solvers[tid]->toDimacs(file);
    }
    else {
        assert(getPrimarySolver() != NULL); // There is at least one solver.
        getPrimarySolver()->toDimacs(file);
    }
}
Lit MultiSolvers::getTrail(int idx) const {
    assert(getPrimarySolver() != NULL); // There is at least one solver.
    return getPrimarySolver()->trail[idx];
}
int MultiSolvers::getTrailSize() const {
    assert(getPrimarySolver() != NULL); // There is at least one solver.
    return getPrimarySolver()->trail.size();
}
uint64_t MultiSolvers::getNumConflicts() const {
    if (allClonesAreBuilt) {
        uint64_t num = 0;
        for(int i = 0; i < nbsolvers; i++)
            num += solvers[i]->conflicts;
        return num;
    }
    else {
        assert(getPrimarySolver() != NULL); // There is at least one solver.
        return getPrimarySolver()->conflicts;
    }
}
uint64_t MultiSolvers::getNumDecisions() const {
    if (allClonesAreBuilt) {
        uint64_t num = 0;
        for(int i = 0; i < nbsolvers; i++)
            num += solvers[i]->decisions;
        return num;
    }
    else {
        assert(getPrimarySolver() != NULL); // There is at least one solver.
        return getPrimarySolver()->decisions;
    }
}
uint64_t MultiSolvers::getNumPropagations() const {
    if (allClonesAreBuilt) {
        uint64_t num = 0;
        for(int i = 0; i < nbsolvers; i++)
            num += solvers[i]->propagations;
        return num;
    }
    else {
        assert(getPrimarySolver() != NULL); // There is at least one solver.
        return getPrimarySolver()->propagations;
    }
}
uint64_t MultiSolvers::getNumRestarts() const {
    if (allClonesAreBuilt) {
        uint64_t num = 0;
        for(int i = 0; i < nbsolvers; i++)
            num += solvers[i]->starts;
        return num;
    }
    else {
        assert(getPrimarySolver() != NULL); // There is at least one solver.
        return getPrimarySolver()->starts;
    }
}
void MultiSolvers::getUnaryAndBinLearntClauses(std::vector<std::vector<int> > &out) {
    int num = allClonesAreBuilt ? nbsolvers : 1;
    
    std::set<std::vector<int> > set;

    // Collects unary learnt clauses
    for (int i=0; i < num; i++) {
        vec<Lit> &unaries = solvers[i]->unariesExpOutside;
        for (int j=0; j < unaries.size(); j++) {
            std::vector<int> unary;
            Lit l = unaries[j];
            int n = sign(l) ? -(var(l) + 1) : +(var(l) + 1);
            unary.push_back(n);
            set.insert(unary);
        }
        unaries.clear();
    }

    // Collects binary learnt clauses
    for (int i=0; i < num; i++) {
        vec<CRef> &ls = solvers[i]->learnts;
        for (int j=0; j < ls.size(); j++) {
            Clause &c = solvers[i]->ca[ls[j]];
            if (c.size() != 2) continue;
            if (c.getExpOutside()) continue;
            
            std::vector<int> bin;
            for (int k=0; k < c.size(); k++) {
                Lit l = c[k];
                int n = sign(l) ? -(var(l) + 1) : +(var(l) + 1);
                bin.push_back(n);
            }
            if (bin[0] > bin[1]) std::swap(bin[0], bin[1]);

            set.insert(bin);
            c.setExpOutside(true);
        }
    }

    // Copy elements in set to out
    out.clear();
    std::copy(set.begin(), set.end(), std::back_inserter(out));
}