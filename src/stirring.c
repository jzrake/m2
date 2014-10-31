#include <stdlib.h>
#include <math.h>
#include "m2.h"


#define DOT(a,b) (a[1]*b[1] + a[2]*b[2] + a[3]*b[3])



double m2stirring_get_field(m2stirring *S, double x[4], double n[4])
{
  int m;
  double Atot = 0.0;

  for (m=0; m<S->num_waves; ++m) {
    double *k = &S->wavenumbers[4*m];
    double *A = &S->fourrieramp[4*m];
    double d = DOT(k, x) + A[0];
    Atot += DOT(A, n) * sin(d);
  }

  return Atot;
}


void m2stirring_next_realization(m2stirring *S)
{
  int m;
  double kfreq = 1.0 / S->wavelength;

  S->wavenumbers = realloc(S->wavenumbers, S->num_waves*4*sizeof(double));
  S->fourrieramp = realloc(S->fourrieramp, S->num_waves*4*sizeof(double));

  for (m=0; m<S->num_waves; ++m) {

    double *k = &S->wavenumbers[4*m];
    double *A = &S->fourrieramp[4*m];

    k[1] = jsw_random_double(&S->random, -1, 1);
    k[2] = jsw_random_double(&S->random, -1, 1);
    k[3] = jsw_random_double(&S->random, -1, 1);
    k[0] = sqrt(DOT(k, k));

    k[1] = (int) (kfreq * k[1] / k[0]) * 2 * M_PI;
    k[2] = (int) (kfreq * k[2] / k[0]) * 2 * M_PI;
    k[3] = (int) (kfreq * k[3] / k[0]) * 2 * M_PI;
    k[0] = sqrt(DOT(k, k));

    A[1] = jsw_random_double(&S->random, -1, 1);
    A[2] = jsw_random_double(&S->random, -1, 1);
    A[3] = jsw_random_double(&S->random, -1, 1);
    A[0] = sqrt(DOT(A, A));

    double Adotkhat = DOT(A, k) / k[0];

    A[1] -= Adotkhat * A[1] / A[0];
    A[2] -= Adotkhat * A[2] / A[0];
    A[3] -= Adotkhat * A[3] / A[0];
    A[0] = S->amplitude * S->wavelength / sqrt(S->num_waves * DOT(A, A));
    /* B0 := S->amplitude ~ sqrt(N) * k0 * A[0] */

    A[1] *= A[0];
    A[2] *= A[0];
    A[3] *= A[0];
    A[0] = jsw_random_double(&S->random, -M_PI, M_PI);
    /* zero-component now holds the phase */
  }
}


void m2stirring_init(m2stirring *S)
{
  S->wavenumbers = NULL;
  S->fourrieramp = NULL;
  S->num_waves = 0;
  S->amplitude = 0.0;
  S->wavelength = 1.0;
  S->stirring_type = M2_STIRRING_NONE;
}


void m2stirring_del(m2stirring *S)
{
  free(S->wavenumbers);
  free(S->fourrieramp);
}
