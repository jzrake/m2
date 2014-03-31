#include <math.h>
#include "m2.h"
#include "hydro.h"


void test_nrmhd_c2p()
{
  printf("\n\n********* %s *********\n", __FUNCTION__);
  m2aux A;
  m2prim P;
  double U[5], UfromP[5];
  double B[4];

  A.m2 = NULL;

  U[DDD] = 1.0;
  U[TAU] = 8.0;
  U[S11] = 0.1;
  U[S22] = 0.2;
  U[S33] = 0.3;

  B[1] = 1.0;
  B[2] = 2.0;
  B[3] = 3.0;

  nrmhd_from_conserved(NULL, U, B, NULL, 1.0, &A, &P);
  nrmhd_from_primitive(NULL, &P, NULL, NULL, 1.0, UfromP, NULL);

  printf("\nstate:\n------\n"); m2_print_state(&P, &A, U);
  printf("\nU out:\n------\n"); m2_print_state(NULL, NULL, UfromP);

  ASSERTEQF(U[TAU], UfromP[TAU]);
  ASSERTEQF(U[DDD], UfromP[DDD]);
  ASSERTEQF(U[S11], UfromP[S11]);
  ASSERTEQF(U[S22], UfromP[S22]);
  ASSERTEQF(U[S33], UfromP[S33]);
}

void test_srmhd_c2p()
{
  printf("\n\n********* %s *********\n", __FUNCTION__);
  m2aux A;
  m2prim P;
  double U[5], UfromP[5];
  double B[4];

  A.m2 = NULL;

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

  printf("\nstate:\n------\n"); m2_print_state(&P, &A, U);
  printf("\nU out:\n------\n"); m2_print_state(NULL, NULL, UfromP);

  ASSERTEQF(U[TAU], UfromP[TAU]);
  ASSERTEQF(U[DDD], UfromP[DDD]);
  ASSERTEQF(U[S11], UfromP[S11]);
  ASSERTEQF(U[S22], UfromP[S22]);
  ASSERTEQF(U[S33], UfromP[S33]);
}


void test_nonrelativistic_limit()
{
  printf("\n\n********* %s *********\n", __FUNCTION__);
  m2sim m2;
  m2prim P;
  m2aux A;
  double flux_nr[8];
  double flux_sr[8];
  double n[4] = { 0.0, 1.0, 0.0, 0.0 };
  double fac = 1e-3;

  A.m2 = &m2;
  P.v1 = 0.1 * fac;
  P.v2 = 0.2 * fac;
  P.v3 = 0.3 * fac;
  P.B1 = 0.1 * sqrt(fac);
  P.B2 = 0.2 * sqrt(fac);
  P.B3 = 0.3 * sqrt(fac);
  P.p = 4.3 * fac;
  P.d = 1.5;

  m2.magnetized = 1;
  m2.relativistic = 0;
  m2sim_from_primitive(&m2, &P, NULL, NULL, 1.0, NULL, &A);
  m2aux_fluxes(&A, n, flux_nr);
  //m2_print_state(NULL, &A, NULL);

  m2.magnetized = 1;
  m2.relativistic = 1;
  m2sim_from_primitive(&m2, &P, NULL, NULL, 1.0, NULL, &A);
  m2aux_fluxes(&A, n, flux_sr);
  //m2_print_state(NULL, &A, NULL);

  printf("nrhyd flux: [%8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e]\n",
	 flux_nr[0], flux_nr[1], flux_nr[2], flux_nr[3],
	 flux_nr[4], flux_nr[5], flux_nr[6], flux_nr[7]);
  printf("srhyd flux: [%8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e]\n",
	 flux_sr[0], flux_sr[1], flux_sr[2], flux_sr[3],
	 flux_sr[4], flux_sr[5], flux_sr[6], flux_sr[7]);
}


void test_srmhd_waves()
{
  printf("\n\n********* %s *********\n", __FUNCTION__);
  m2prim P;
  m2aux A;
  double evals1[8];
  double evals2[8];
  double n[4] = { 0.0, 1.0, 0.0, 0.0 };

  A.m2 = NULL;
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


  P.v1 = 0.00010;
  P.v2 = 0.00000;
  P.v3 = 0.00000;
  P.B1 = 0.00000;
  P.B2 = 0.00000;
  P.B3 = 0.00001;
  P.d  = 1.00000;
  P.p  = 0.00100;

  srmhd_from_primitive(NULL, &P, NULL, NULL, 1.0, NULL, &A);
  srmhd_eigenvalues(&A, n, evals1);
  nrmhd_from_primitive(NULL, &P, NULL, NULL, 1.0, NULL, &A);
  nrmhd_eigenvalues(&A, n, evals2);

  printf("\nSRMHD (with non-relativistic conditions):\n------\n");
  printf("slow magnetosonic: %+8.6e %+8.6e\n", evals1[2], evals1[5]);
  printf("Alfven waves     : %+8.6e %+8.6e\n", evals1[1], evals1[6]);
  printf("fast magnetosonic: %+8.6e %+8.6e\n", evals1[0], evals1[7]);

  printf("\nNRMHD:\n------\n");
  printf("slow magnetosonic: %+8.6e %+8.6e\n", evals2[2], evals2[5]);
  printf("Alfven waves     : %+8.6e %+8.6e\n", evals2[1], evals2[6]);
  printf("fast magnetosonic: %+8.6e %+8.6e\n", evals2[0], evals2[7]);


  P.v1 = 0.0;
  P.v2 = 0.0;
  P.v3 = 0.0;
  P.B1 = 1.0;
  P.B2 = 0.0;
  P.B3 = 0.0;
  P.d = 1.0;
  P.p = 1.0;
  srmhd_from_primitive(NULL, &P, NULL, NULL, 1.0, NULL, &A);
  srmhd_eigenvalues(&A, n, evals1);
  printf("\nSRMHD (longitudinal field):\n------\n");
  printf("slow magnetosonic: %+8.6e %+8.6e\n", evals1[2], evals1[5]);
  printf("Alfven waves     : %+8.6e %+8.6e\n", evals1[1], evals1[6]);
  printf("fast magnetosonic: %+8.6e %+8.6e\n", evals1[0], evals1[7]);


  P.v1 = 0.0;
  P.v2 = 0.0;
  P.v3 = 0.0;
  P.B1 = 0.0;
  P.B2 = 1.0;
  P.B3 = 0.0;
  P.d = 1.0;
  P.p = 1.0;
  srmhd_from_primitive(NULL, &P, NULL, NULL, 1.0, NULL, &A);
  srmhd_eigenvalues(&A, n, evals1);
  printf("\nSRMHD (transverse field):\n------\n");
  printf("slow magnetosonic: %+8.6e %+8.6e\n", evals1[2], evals1[5]);
  printf("Alfven waves     : %+8.6e %+8.6e\n", evals1[1], evals1[6]);
  printf("fast magnetosonic: %+8.6e %+8.6e\n", evals1[0], evals1[7]);
}

void test_quartic()
{
  printf("\n\n********* %s *********\n", __FUNCTION__);
  double roots[4];
  double A1[5] = {+1, -1, +1, -2, +1};
  double A2[5] = {+11166.1589491647,
		  -49480.3771267108,
		  +82221.6612223221,
		  -60722.5167663058,
		  +16816.5669665401};
  double A3[5] = {+11782.9065783065,
		  -48988.7329357975,
		  +76377.5655986056,
		  -52923.3676374477,
		  +13751.6553530872};
  double A4[5] = {+00.0335111570,
		  -01.1655211505,
		  +12.7687728282,
		  -56.1746104789,
		  +85.6898329342};
  double A5[5] = {+0751.8716976283,
		  -3465.6791436275,
		  +5989.3459074759,
		  -4599.4258609677,
		  +1324.2662465491};
  int i, n, nr;
  double y, *A, *As[5] = {A1, A2, A3, A4, A5};

  for (i=0; i<5; ++i) {
    A = As[i];
    nr = m2_solve_quartic_equation(A[4], A[3], A[2], A[1], A[0], roots);
    printf("\nthere are %d roots:\n", nr);
    for (n=0; n<4; ++n) {
      y = (A[4]*pow(roots[n], 4) +
	   A[3]*pow(roots[n], 3) +
	   A[2]*pow(roots[n], 2) +
	   A[1]*pow(roots[n], 1) +
	   A[0]*pow(roots[n], 0));
      printf("r%d = %f, y(r%d) = %8.6e\n", n, roots[n], n, y);
    }
  }
}

void m2_self_test()
{
  test_nonrelativistic_limit();
  test_nrmhd_c2p();
  test_srmhd_c2p();
  test_srmhd_waves();
  test_quartic();
}
