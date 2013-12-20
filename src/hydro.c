#include <math.h>
#include "hydro.h"


#define SYSTEM_CHOICE(m2,E)						\
  do {									\
    if (0) { }								\
    else if (m2->physics == (M2_RELATIVISTIC | M2_MAGNETIZED)) {	\
      MSG(FATAL, "srmhd not yet implemented");				\
      return -1;							\
    }									\
    else if (m2->physics == (M2_RELATIVISTIC | M2_UNMAGNETIZED)) {	\
      return srhyd_##E;							\
    }									\
    else if (m2->physics == (M2_NONRELATIVISTIC | M2_MAGNETIZED)) {	\
      return nrmhd_##E;							\
    }									\
    else if (m2->physics == (M2_NONRELATIVISTIC | M2_UNMAGNETIZED)) {	\
      return nrhyd_##E;							\
    }									\
    else {								\
      MSG(FATAL, "unknown physics");					\
      return -1;							\
    }									\
  } while (0)								\


int m2sim_from_primitive(m2sim *m2, m2prim *P, double *B, double *X, double dV,
			 double *U, m2aux *aux)
{
  SYSTEM_CHOICE(m2, from_primitive(m2, P, B, X, dV, U, aux));
}
int m2sim_from_conserved(m2sim *m2, double *U, double *B, double *X, double dV,
			 m2aux *aux, m2prim *P)
{
  SYSTEM_CHOICE(m2, from_conserved(m2, U, B, X, dV, aux, P));
}
int m2sim_from_auxiliary(m2sim *m2, m2aux *aux, double *X, double dV,
			 m2prim *P, double *U)
{
  SYSTEM_CHOICE(m2, from_auxiliary(m2, aux, X, dV, P, U));
}
int m2aux_eigenvalues(m2aux *aux, double n[4], double *evals)
{
  SYSTEM_CHOICE(aux->m2, eigenvalues(aux, n, evals));
}
int m2aux_add_geometrical_source_terms(m2aux *aux, double x0[4], double x1[4],
				       double *U)
{
  if (aux->m2->geometry == M2_CARTESIAN) {
    return 0;
  }
  else if (aux->m2->geometry == M2_SPHERICAL) {
    if (aux->m2->physics & M2_MAGNETIZED) {
      MSG(FATAL, "MHD source terms not yet added");
    }
    const double p  = aux->gas_pressure;
    const double dt = x1[0] - x0[0];
    const double r0 = x0[1];
    const double r1 = x1[1];
    const double t0 = x0[2];
    const double t1 = x1[2];
    const double f0 = x0[3];
    const double f1 = x1[3];
    const double u0 = aux->velocity_four_vector[0];
    const double v2 = aux->velocity_four_vector[2] / u0;
    const double v3 = aux->velocity_four_vector[3] / u0;
    const double S1 = aux->momentum_density[1];
    const double S2 = aux->momentum_density[2];
    const double S3 = aux->momentum_density[3];
    U[S11] += dt * 0.5 * (-1.0) *
      (r1 - r0)*(r1 + r0)*(2*p + S2*v2 + S3*v3)*(f1 - f0)*(cos(t1) - cos(t0));
    U[S22] += dt * 0.5 *
      (r1*r1 - r0*r0)*(f1 - f0)*((0 + S1*v2)*(cos(t1) - cos(t0)) +
				 (p + S3*v3)*(sin(t1) - sin(t0)));
    U[S33] += dt * 0.5 * v3 *
      (r1*r1 - r0*r0)*(f1 - f0)*(S1*(cos(t1) - cos(t0)) -
				 S2*(sin(t1) - sin(t0)));
    return 0;
  }
  else {
    MSG(FATAL, "invalid geometry");
    return 0;
  }
}
int m2aux_fluxes(m2aux *aux, double n[4], double *F)
{
  const double SR = aux->m2->physics & M2_RELATIVISTIC ? 1.0 : 0.0;
  const double u0 = aux->velocity_four_vector[0];
  const double u1 = aux->velocity_four_vector[1];
  const double u2 = aux->velocity_four_vector[2];
  const double u3 = aux->velocity_four_vector[3];
  const double b0 = aux->magnetic_four_vector[0];
  const double b1 = aux->magnetic_four_vector[1];
  const double b2 = aux->magnetic_four_vector[2];
  const double b3 = aux->magnetic_four_vector[3];
  const double S0 = aux->momentum_density[0]; /* T^{0,0} = tau + D */
  const double S1 = aux->momentum_density[1];
  const double S2 = aux->momentum_density[2];
  const double S3 = aux->momentum_density[3];
  const double d0 = aux->comoving_mass_density;
  const double pg = aux->gas_pressure;
  const double pb = aux->magnetic_pressure;
  const double D0 = d0 * u0;
  const double n1 = n[1];
  const double n2 = n[2];
  const double n3 = n[3];
  const double T0 = S0 - D0; /* tau */
  const double v1 = u1 / u0;
  const double v2 = u2 / u0;
  const double v3 = u3 / u0;
  const double B1 = b1 * u0 - b0 * u1 * SR;
  const double B2 = b2 * u0 - b0 * u2 * SR;
  const double B3 = b3 * u0 - b0 * u3 * SR;
  const double Bn = B1*n1 + B2*n2 + B3*n3;
  const double vn = v1*n1 + v2*n2 + v3*n3;
  F[DDD] = D0 * vn;
  F[TAU] = T0 * vn + (pg + pb) * vn - b0 * Bn / u0;
  F[S11] = S1 * vn + (pg + pb) * n1 - b1 * Bn / u0;
  F[S22] = S2 * vn + (pg + pb) * n2 - b2 * Bn / u0;
  F[S33] = S3 * vn + (pg + pb) * n3 - b3 * Bn / u0;
  F[B11] = B1 * vn - Bn * v1;
  F[B22] = B2 * vn - Bn * v2;
  F[B33] = B3 * vn - Bn * v3;
  return 0;
}
double m2aux_maximum_wavespeed(m2aux *aux)
{
  double evals[8];
  double n[4] = { 0.0,
		  aux->velocity_four_vector[1],
		  aux->velocity_four_vector[2],
		  aux->velocity_four_vector[3] };
  double N = sqrt(n[1]*n[1] + n[2]*n[2] + n[3]*n[3]);
  if (fabs(N) < 1e-10) { /* if no velocity, then just use sound speed */
    n[1] = 0.0;
    n[2] = 0.0;
    n[3] = 0.0;
  }
  else {
    n[1] /= N;
    n[2] /= N;
    n[3] /= N;
  }
  m2aux_eigenvalues(aux, n, evals);
  return fabs(evals[0]) > fabs(evals[7]) ? fabs(evals[0]) : fabs(evals[7]);
}
