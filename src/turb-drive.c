#include <math.h>
#include "jsw_rand.h"


#define NWAVES 10
#define DOT(a,b) (a[1]*b[1] + a[2]*b[2] + a[3]*b[3])

#define random_real jsw_random_double

static double VectorPotentialAmplitudes[NWAVES][4];
static double VectorPotentialWavenumbers[NWAVES][4];


void init_amplitudes()
{
  int m;
  double kfreq = 16;
  for (m=0; m<NWAVES; ++m) {

    double *A = VectorPotentialAmplitudes[m];
    double *k = VectorPotentialWavenumbers[m];

    k[1] = random_real(-1, 1);
    k[2] = random_real(-1, 1);
    k[3] = random_real(-1, 1);
    k[0] = sqrt(DOT(k, k));

    k[1] = (int) (kfreq * k[1] / k[0]) * 2 * M_PI;
    k[2] = (int) (kfreq * k[2] / k[0]) * 2 * M_PI;
    k[3] = (int) (kfreq * k[3] / k[0]) * 2 * M_PI;
    k[0] = sqrt(DOT(k, k));

    A[1] = random_real(-1, 1);
    A[2] = random_real(-1, 1);
    A[3] = random_real(-1, 1);
    A[0] = sqrt(DOT(A, A));

    double Adotkhat = DOT(A, k) / k[0];

    A[1] -= Adotkhat * A[1] / A[0];
    A[2] -= Adotkhat * A[2] / A[0];
    A[3] -= Adotkhat * A[3] / A[0];
    A[0] = sqrt(DOT(A, A));

    A[1] /= A[0];
    A[2] /= A[0];
    A[3] /= A[0];
    A[0] = random_real(-M_PI, M_PI); /* zero-component now holds the phase */
  }
}


double m2_stochastic_vector_potential(double x[4], double n[4])
{
  static int num_calls = 0;

  if (num_calls == 0) {
    jsw_seed(1);
    init_amplitudes();
  }

  num_calls +=1 ;

  int m;
  double Atot = 0.0;

  for (m=0; m<NWAVES; ++m) {
    double *A = VectorPotentialAmplitudes[m];
    double *k = VectorPotentialWavenumbers[m];
    double d = DOT(k, x) + A[0];
    Atot += 0.1 * DOT(A, n) * sin(d);
  }

  return Atot;
}
