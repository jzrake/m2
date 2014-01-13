#include <math.h>
#include "m2.h"
#include "hydro.h"

void m2_self_test()
{
  m2aux A;
  m2prim P;
  double U[5], UfromP[5];
  double B[4];

  U[DDD] = 1.0;
  U[TAU] = 8.0;
  U[S11] = 0.1;
  U[S22] = 0.2;
  U[S33] = 0.3;

  B[1] = 1.0;
  B[2] = 2.0;
  B[3] = 3.0;

  srmhd_from_conserved(NULL, U, B, NULL, 1.0, &A, &P);
  srmhd_from_primitive(NULL, &P, NULL, NULL, 1.0, UfromP, NULL);

  printf("state:\n------\n"); m2_print_state(&P, &A, U);
  printf("U out:\n------\n"); m2_print_state(NULL, NULL, UfromP);

  ASSERTEQF(U[TAU], UfromP[TAU]);
  ASSERTEQF(U[DDD], UfromP[DDD]);
  ASSERTEQF(U[S11], UfromP[S11]);
  ASSERTEQF(U[S22], UfromP[S22]);
  ASSERTEQF(U[S33], UfromP[S33]);
}
