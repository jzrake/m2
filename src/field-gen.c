#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "jsw_rand.h"
#include "spheromak.h"


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


double m2_spheromak(double x[4], double n[4], int N, int L, int M)
{
  double R = 1e-16 + sqrt(x[1]*x[1] + x[2]*x[2]);
  double r = 1e-16 + sqrt(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
  double t = acos(x[3] / r);
  double p = atan2(x[2], x[1]);
  double rhat[4] = {0, x[1]/r, x[2]/r, x[3]/r};
  double that[4] = {0, x[3]/r * x[1]/R, x[3]/r * x[2]/R, -R/r};
  double phat[4] = {0,-x[2]/R, x[1]/R, 0};

  double r0 = SphericalBesselJZero(N, L);
  double Br, Bt, Bp, B[4];

  if (r < r0) {
    SpheromakInternal(N, L, M, r, t, p, &Br, &Bt, &Bp);
  }
  else {
    SpheromakExternal(N, L, M, r, t, p, &Br, &Bt, &Bp);
  }

  B[1] = Br * rhat[1] + Bt * that[1] + Bp * phat[1];
  B[2] = Br * rhat[2] + Bt * that[2] + Bp * phat[2];
  B[3] = Br * rhat[3] + Bt * that[3] + Bp * phat[3];

  return B[1]*n[1] + B[2]*n[2] + B[3]*n[3];
}



double m2_tubomak_vector_potential(double x[4], double n[4])
{
#define J CylindricalBesselJ
  double a = 3.83170597021; /* first zero of A-phi will be at r=1 */
  double r = sqrt(x[1]*x[1] + x[2]*x[2]);
  if (fabs(r) < 1e-12) r = 1e-12;
  double Af =   a*J(1, a*r);
  double Az = ((a*J(1, a*r))/r + (pow(a,2)*(J(0, a*r) - J(2, a*r)))/2.)/a;
  double Ax = Af * (-x[2]/r);
  double Ay = Af * (+x[1]/r);
  if (r < 1) {
    return Ax*n[1] + Ay*n[2] + Az*n[3];
  }
  else {
    r = 1;
    double Az = ((a*J(1, a*r))/r + (pow(a,2)*(J(0, a*r) - J(2, a*r)))/2.)/a;
    return Az * n[3];
  }
#undef J
}

