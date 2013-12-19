#include <math.h>
#include "hydro.h"


#define gamma_law_index (4./3.)


int nrmhd_from_primitive(m2sim *m2, m2prim *P, double *B, double *X, double dV,
			 double *U, m2aux *aux)
{
  double v1 = P->v1;
  double v2 = P->v2;
  double v3 = P->v3;
  double d = P->d;
  double p = P->p;
  double u = p / (gamma_law_index - 1.0);
  double Eb = 0.5 * (P->B1*P->B1 + P->B2*P->B2 + P->B3*P->B3);

  if (aux) {
    aux->velocity_four_vector[0] = 1.0;
    aux->velocity_four_vector[1] = v1;
    aux->velocity_four_vector[2] = v2;
    aux->velocity_four_vector[3] = v3;
    aux->magnetic_four_vector[0] = 0.0;
    aux->magnetic_four_vector[1] = P->B1;
    aux->magnetic_four_vector[2] = P->B2;
    aux->magnetic_four_vector[3] = P->B3;
    aux->momentum_density[0] = d * 0.5 * (v1*v1 + v2*v2 + v3*v3) + u + Eb + d;
    aux->momentum_density[1] = d * v1;
    aux->momentum_density[2] = d * v2;
    aux->momentum_density[3] = d * v3;
    aux->comoving_mass_density = d;
    aux->gas_pressure = p;
    aux->magnetic_pressure = Eb;
    aux->m2 = m2;
  }
  if (U) {
    U[DDD] = dV * (d);
    U[TAU] = dV * (d * 0.5 * (v1*v1 + v2*v2 + v3*v3) + u + Eb);
    U[S11] = dV * (d * v1);
    U[S22] = dV * (d * v2);
    U[S33] = dV * (d * v3);
  }

  return 0;
}


int nrmhd_from_conserved(m2sim *m2, double *U, double *B, double *X, double dV,
			 m2aux *aux, m2prim *P)
{
  double T0 = U[TAU] / dV; /* total energy density */
  double S1 = U[S11] / dV;
  double S2 = U[S22] / dV;
  double S3 = U[S33] / dV;
  double D0 = U[DDD] / dV;

  double Ek = 0.5 * (S1*S1 + S2*S2 + S3*S3) / D0;
  double Eb = 0.5 * (B[1]*B[1] + B[2]*B[2] + B[3]*B[3]);
  double pg = (T0 - Ek - Eb) * (gamma_law_index - 1.0);
  double v1 = S1 / D0;
  double v2 = S2 / D0;
  double v3 = S3 / D0;

  if (pg < 0.0) {
    MSGF(FATAL, "got negative p: %f %f %f", T0, Ek, Eb);
  }

  if (aux) {
    aux->velocity_four_vector[0] = 1.0;
    aux->velocity_four_vector[1] = v1;
    aux->velocity_four_vector[2] = v2;
    aux->velocity_four_vector[3] = v3;
    aux->magnetic_four_vector[0] = 0.0;
    aux->magnetic_four_vector[1] = B[1];
    aux->magnetic_four_vector[2] = B[2];
    aux->magnetic_four_vector[3] = B[3];
    aux->momentum_density[0] = T0 + D0; /* odd for Newtonian, see fluxes */
    aux->momentum_density[1] = S1;
    aux->momentum_density[2] = S2;
    aux->momentum_density[3] = S3;
    aux->comoving_mass_density = D0;
    aux->gas_pressure = pg;
    aux->magnetic_pressure = Eb;
    aux->m2 = m2;
  }

  if (P) {
    P->v1 = v1;
    P->v2 = v2;
    P->v3 = v3;
    P->d = D0;
    P->p = pg;
  }

  return 0;
}


int nrmhd_from_auxiliary(m2sim *m2, m2aux *aux, double *X, double dV,
			 m2prim *P, double *U)
{
  return 0;
}


int nrmhd_eigenvalues(m2aux *aux, double n[4], double *evals)
{
  const double d = aux->comoving_mass_density;
  const double p = aux->gas_pressure;
  const double u0 = aux->velocity_four_vector[0];
  const double v1 = aux->velocity_four_vector[1] / u0;
  const double v2 = aux->velocity_four_vector[2] / u0;
  const double v3 = aux->velocity_four_vector[3] / u0;
  const double vn = v1*n[1] + v2*n[2] + v3*n[3];
  const double cs = sqrt(gamma_law_index * p / d);
  evals[0] = vn - cs;
  evals[1] = vn;
  evals[2] = vn;
  evals[3] = vn;
  evals[4] = vn + cs;
  return 0;
}
