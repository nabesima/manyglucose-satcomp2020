//#define __STDC_LIMIT_MACROS
//#define __STDC_FORMAT_MACROS
#include "parallel/MultiSolvers.h"
//#include "utils/System.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>

using namespace std;
using namespace Glucose;

extern "C" {
    static const char * sig = "manyglucose-4.1.27";
};

class IPAsirManyGlucose: public MultiSolvers {
    vec<Lit> assumptions, clause;
    int szfmap;
    unsigned char * fmap;
    bool nomodel;
    unsigned long long calls;
    double total_time;
    void reset() {
        if (fmap)
            delete[] fmap, fmap = 0, szfmap = 0;
    }
    Lit import(int lit) {
        while (abs(lit) > nVars())
            (void) newVar();
        return mkLit(Var(abs(lit) - 1), (lit > 0));
    }
    void ana() {
        fmap = new unsigned char[szfmap = nVars()];
        memset(fmap, 0, szfmap);
        for (int i = 0; i < conflict.size(); i++) {
            int tmp = var(conflict[i]);
            assert(0 <= tmp && tmp < szfmap);
            fmap[tmp] = 1;
        }
    }
    double ps(double s, double t) {
        return t ? s / t : 0;
    }
public:
    IPAsirManyGlucose() :
            szfmap(0), fmap(0), nomodel(false), calls(0), total_time(0.0) {
        // MiniSAT by default produces non standard conforming messages.
        // So either we have to set this to '0' or patch the sources.
        setVerbosity(1);
    }
    ~IPAsirManyGlucose() {
        reset();
    }
    void printStatistics() {
        printStats();
    }
    void add(int lit) {
        reset();
        nomodel = true;
        if (lit)
            clause.push(import(lit));
        else
            addClause(clause), clause.clear();
    }
    void assume(int lit) {
        reset();
        nomodel = true;
        assumptions.push(import(lit));
    }
    int solve() {
        double start = realTime();
        setvbuf(stdout, (char *) NULL, _IONBF, 0);
        calls++;
        reset();
        lbool res = MultiSolvers::solve(assumptions, false);
        assumptions.clear();
        nomodel = (res != l_True);
        printStatistics();
        total_time += realTime() - start;
        return (res == l_Undef) ? 0 : (res == l_True ? 10 : 20);
    }
    int val(int lit) {
        if (nomodel)
            return 0;
        lbool res = getModelValue(import(lit));
        return (res == l_True) ? lit : -lit;
    }
    int failed(int lit) {
        if (!fmap)
            ana();
        int tmp = var(import(lit));
        assert(0 <= tmp && tmp < nVars());
        return fmap[tmp] != 0;
    }
    double get_time() {
        return total_time;
    }
    void stats() {
       printStats();
    }
};

extern "C" {
#include "ipasir.h"
static IPAsirManyGlucose * import(void * s) {
    return (IPAsirManyGlucose*) s;
}
const char * ipasir_signature() {
    return sig;
}
void * ipasir_init() {
    return new IPAsirManyGlucose();
}
void ipasir_release(void * s) {
    import(s)->stats();
    delete import(s);
}
int ipasir_solve(void * s) {
    return import(s)->solve();
}
void ipasir_add(void * s, int l) {
    import(s)->add(l);
}
void ipasir_assume(void * s, int l) {
    import(s)->assume(l);
}
int ipasir_val(void * s, int l) {
    return import(s)->val(l);
}
int ipasir_failed(void * s, int l) {
    return import(s)->failed(l);
}
void ipasir_set_terminate(void * s, void * state,
        int (*callback)(void * state)) {
    import(s)->setTermCallback(state, callback);
}
}
