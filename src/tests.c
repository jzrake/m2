#include <math.h>
#include "m2.h"
#include "hydro.h"

void test_srmhd_c2p()
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

void test_srmhd_waves()
{
  m2prim P;
  m2aux A;
  double evals1[8];
  double evals2[8];
  double n[4] = { 0.0, 1.0, 0.0, 0.0 };

  P.v1 = 0.1;
  P.v2 = 0.2;
  P.v3 = 0.3;
  P.B1 = 0.0;
  P.B2 = 0.0;
  P.B3 = 0.0;
  P.d = 1.5;
  P.p = 4.3;

  srmhd_from_primitive(NULL, &P, NULL, NULL, 1.0, NULL, &A);
  srmhd_eigenvalues(&A, n, evals1);
  srhyd_eigenvalues(&A, n, evals2);

  ASSERTEQF(evals1[0], evals2[0]);
  ASSERTEQF(evals1[7], evals2[7]);


  P.v1 = 0.00001;
  P.v2 = 0.00002;
  P.v3 = 0.00003;
  P.B1 = 0.0001;
  P.B2 = 0.0001;
  P.B3 = 0.0;
  P.d = 1.0;
  P.p = 0.0001;

  srmhd_from_primitive(NULL, &P, NULL, NULL, 1.0, NULL, &A);
  srmhd_eigenvalues(&A, n, evals1);
  nrmhd_from_primitive(NULL, &P, NULL, NULL, 1.0, NULL, &A);
  nrmhd_eigenvalues(&A, n, evals2);

  printf("SRMHD:\n------\n");
  printf("slow magnetosonic: %+8.6e %+8.6e\n", evals1[2], evals1[5]);
  printf("Alfven waves     : %+8.6e %+8.6e\n", evals1[1], evals1[6]);
  printf("fast magnetosonic: %+8.6e %+8.6e\n", evals1[0], evals1[7]);

  printf("NRMHD:\n------\n");
  printf("slow magnetosonic: %+8.6e %+8.6e\n", evals2[2], evals2[5]);
  printf("Alfven waves     : %+8.6e %+8.6e\n", evals2[1], evals2[6]);
  printf("fast magnetosonic: %+8.6e %+8.6e\n", evals2[0], evals2[7]);
}

void m2_self_test()
{
  test_srmhd_c2p();
  test_srmhd_waves();
}
