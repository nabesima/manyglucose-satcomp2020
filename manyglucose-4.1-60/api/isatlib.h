
#ifndef isatlib_h_INCLUDED
#define isatlib_h_INCLUDED

#include "core/Solver.h"
//#include "simp/SimpSolver.h"

#ifdef __cplusplus
extern "C" {
#endif
  
  void isat_set_verbosity(void * s, int v);
  int isat_vars(void * s);
  int isat_clauses(void * s);
  void isat_new_var(void * s, int howmany);
  void isat_add_clauses(void * s, int * clauses, int length);
  void isat_add_assumptions(void * s, int * literals, int length);
  int * isat_model(void * s);
  int * isat_failed(void * s, int * literals, int length);
  void isat_add_blocking_clause(void * s, int * literals, int length);
  void isat_clean(void * s, int * p);
  double isat_get_time(void * s);
  int isat_add_PB(void * s, int * weights, int * literals, int length, int k, bool isgeq);
  int isat_add_PB_with_assumption(void * s, int * weights, int * literals, int length, int k, bool isgeq, int assumption);
  void isat_print_clauses(void * s);
  // int * isat_MUS(void * s);
  void isat_freeze(void * s, int literal);
  int isat_is_eliminated(void * s, int literal);
  // void ipasir_delete_learnts(void * s);
  // void ipasir_set_polarity(void * s, int * polarity);
  // int * ipasir_polarity(void * s);
  
#ifdef __cplusplus
}
#endif
  
#endif
