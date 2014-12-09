#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "jsw_rand.h"


typedef complex double Complex;



typedef struct fourier_mode
{
  double k[4];
  Complex A[4];
} fourier_mode;



double m2_force_free_vector_potential(double x[4], double n[4], int model, int k2)
{
#define RAND jsw_random_double(&rand, -1, 1)
  int m,i,j,k;
  jsw_rand_t rand;
  jsw_seed(&rand, model);


  int k2_sphere = k2;//26;
  int k_cube = floor(sqrt(k2_sphere)) + 1;
  fourier_mode *modes = NULL;
  int num_modes = 0;
  Complex A[4] = {0, 0, 0, 0};
  double amp;


  for (i=-k_cube; i<=k_cube; ++i) {
    for (j=-k_cube; j<=k_cube; ++j) {
      for (k=-k_cube; k<=k_cube; ++k) {
  	fourier_mode M;

  	if (i*i + j*j + k*k != k2_sphere) {
  	  continue;
  	}
  	else {
  	  M.k[0] = 0.0;
  	  M.k[1] = i;
  	  M.k[2] = j;
  	  M.k[3] = k;
  	  M.A[0] = 0.0;
  	  M.A[1] = RAND + RAND*I;
  	  M.A[2] = RAND + RAND*I;
  	  M.A[3] = RAND + RAND*I;
  	  amp = sqrt(M.A[1]*conj(M.A[1]) + M.A[2]*conj(M.A[2]) + M.A[3]*conj(M.A[3]));
  	  M.A[1] /= amp;
  	  M.A[2] /= amp;
  	  M.A[3] /= amp;
  	  num_modes += 1;
  	  modes = (fourier_mode *) realloc(modes, num_modes * sizeof(fourier_mode));
  	  modes[num_modes-1] = M;

  	  /* printf("k[%d] = [%d %d %d] is on shell\n", num_modes, i, j, k); */
  	}
      }
    }
  }


  for (m=0; m<num_modes; ++m) {
    fourier_mode M = modes[m];
    double a = sqrt(M.k[1]*M.k[1] + M.k[2]*M.k[2] + M.k[3]*M.k[3]);
    Complex K[4] = {0, I*M.k[1], I*M.k[2], I*M.k[3]};
    Complex Ikx  = (K[1]*x[1] + K[2]*x[2] + K[3]*x[3]) * M_PI;
    Complex P[4] = {0, M.A[1], M.A[2], M.A[3]}; /* a times psi */

    Complex T[4] = {0, /* T = K cross (a psi) */
		    K[2]*P[3] - K[3]*P[2],
		    K[3]*P[1] - K[1]*P[3],
		    K[1]*P[2] - K[2]*P[1]};

    Complex S[4] = {0, /* S = K cross T / alpha */
		    (K[2]*T[3] - K[3]*T[2])/a,
		    (K[3]*T[1] - K[1]*T[3])/a,
		    (K[1]*T[2] - K[2]*T[1])/a};

    A[1] += (S[1] + T[1])/a * cexp(Ikx);
    A[2] += (S[2] + T[2])/a * cexp(Ikx);
    A[3] += (S[3] + T[3])/a * cexp(Ikx);
  }


  free(modes);
  return creal(A[1]*n[1] + A[2]*n[2] + A[3]*n[3]);

#undef RAND
}


double m2_force_free_magnetic_field(double x[4], double n[4], int model, int k2)
{
  /* correct up to factor of alpha */
  return m2_force_free_vector_potential(x, n, model, k2);
}


double m2_magnetic_rope_vector_potential(double x[4], double n[4])
{
  double m[4] = {0, 1/sqrt(3), 1/sqrt(3), 1/sqrt(3)};
  //double m[4] = {0, 0, 0, 1};
  double Ap = 1.0; /* poloidal A -> toroidal B */
  double At = 0.0; /* toroidal A -> poloidal B */
  double L = 0.05;
  double f = sqrt(x[1]*x[1] + x[2]*x[2]); /* cylindrical radius */
  double R = sqrt(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
  //double D = 1 / (1 + pow(R/L, 3));
  double D = pow(R/L,2) * exp(-pow(R/L, 4));
  double mdotrhat = (m[1]*x[1] + m[2]*x[2] + m[3]*x[3]) / R;
  double phihat[4] = {0,
		      phihat[1] = x[2]/R*m[3] - x[3]/R*m[2],
		      phihat[2] = x[3]/R*m[1] - x[1]/R*m[3],
		      phihat[3] = x[1]/R*m[2] - x[2]/R*m[1]};
  double A[4] = {0,
		 Ap * (3*x[1]/R*mdotrhat - m[1]) * D,
		 Ap * (3*x[2]/R*mdotrhat - m[2]) * D,
		 Ap * (3*x[3]/R*mdotrhat - m[3]) * D};

  A[1] += At * phihat[1] * (f/L) * D;
  A[2] += At * phihat[2] * (f/L) * D;
  A[3] += At * phihat[3] * (f/L) * D;

  return A[1]*n[1] + A[2]*n[2] + A[3]*n[3];
}
