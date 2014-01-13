#ifndef M2_PROBLEMS_HEADER
#define M2_PROBLEMS_HEADER

#include "m2.h"

void mwn_initial_data(m2vol *V);
void mwn_boundary_conditions(m2sim *m2);
void mwn_analysis(m2sim *m2);
void mwn_initialize_problem(m2sim *m2);


#endif /* M2_PROBLEMS_HEADER */
