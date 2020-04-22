#include "../api/isatlib.h"
#include "../api/ipasirmanyglucose.cc"

using namespace std;
using namespace Glucose;

class IPAsirManyGlucoseExtended : public IPAsirManyGlucose {

  typedef IPAsirManyGlucose base;

public:
  IPAsirManyGlucoseExtended () {
  }

  void new_var(int howmany){
     int lit = howmany + nVars();
     while (abs(lit) > nVars()) (void) newVar();
  }

  void add_clauses(int * clauses, int length){
     for (int i=0; i<length; i++){
        //cout << clauses[i] << " ";
        //if (!clauses[i]) cout << endl;
        base::add(clauses[i]);
     }
  }

  void add_assumptions(int * literals, int length){
     for (int i=0; i<length; i++)
        base::assume(literals[i]);
  }

  int * model(){
     int vars = nVars();
     int * ansModel = (int *)malloc(sizeof(int)*vars);
     for(int i=0; i<vars; i++)
        ansModel[i] = val(i+1);
     return ansModel;
  }

  int * get_failed(int * literals, int length){
     int * failed = (int *)malloc(sizeof(int)*length);
     for(int i=0; i<length; i++)
        failed[i] = base::failed(literals[i]);
     return failed;
  }

  void add_blocking_clause(int * literals, int length){
     // all variable
     if (literals == NULL){
        int * ansModel = model();
        for(int i=0; i<nVars(); i++)
          base::add(-ansModel[i]);
        base::add(0);
     }
     else {
        for(int i=0; i<length; i++)
          base::add(-literals[i]);
        base::add(0);
     }
  }

  void release(int * p){
     free(p);
  }

  void print_clauses(){
  }

  void set_frozen(int literal){
     base::setFrozen(literal, true);
  }

  int is_eliminated(int literal){
      return base::isEliminated(literal) ? 1 : 0;
  }

};

static IPAsirManyGlucoseExtended * importex (void * s) { return (IPAsirManyGlucoseExtended*) s; }
void isat_set_verbosity(void * s, int v){ importex(s) -> setVerbosity(v); }
int isat_vars(void * s){ return importex(s) -> nVars(); }
int isat_clauses(void * s){ return importex(s) -> nClauses(); }
void isat_new_var(void * s, int howmany){importex(s) -> new_var(howmany);}
void isat_add_clauses(void * s, int * clauses, int length){ importex(s) -> add_clauses(clauses,length);}
void isat_add_assumptions(void * s, int * literals, int length){ importex(s) -> add_assumptions(literals, length);}
int * isat_model(void * s){ return importex(s) -> model();}
int * isat_failed(void * s, int * literals, int length){ return importex(s) -> get_failed(literals,length);}
void isat_add_blocking_clause(void * s, int * literals, int length){ importex(s) -> add_blocking_clause(literals,length);}
void isat_clean(void * s, int * p){ importex(s) -> release(p);}
double isat_get_time(void * s){ return importex(s) -> get_time();}
// int isat_add_PB(void * s, int * weights, int * literals, int length, int k, bool isgeq){ return importex(s) -> add_PB_constraint(weights, literals, length, k, isgeq);}
// int isat_add_PB_with_assumption(void * s, int * weights, int * literals, int length, int k, bool isgeq, int assumption){ return importex(s) -> add_PB_constraint_with_assumption(weights, literals, length, k, isgeq, assumption);}
void isat_print_clauses(void * s){ importex(s) -> print_clauses();}
// int * isat_MUS(void * s);
void isat_freeze(void * s, int literal){ importex(s) -> set_frozen(literal);}
int isat_is_eliminated(void * s, int literal){ return importex(s) -> is_eliminated(literal);}
// void ipasir_delete_learnts(void * s);
// void ipasir_set_polarity(void * s, int * polarity);
// int * ipasir_polarity(void * s);
