#include <math.h>
#include "hydro.h"


#define gamma_law_index (4./3.)


int srhyd_from_primitive(m2sim *m2, m2prim *P, double *B, double *X, double dV,
			 double *U, m2aux *aux)
{
  double u0 = 1.0 / sqrt(1.0 - (P->v1*P->v1 + P->v2*P->v2 + P->v3*P->v3));
  double u1 = u0 * P->v1;
  double u2 = u0 * P->v2;
  double u3 = u0 * P->v3;
  double d = P->d;
  double p = P->p;
  double u = p / (gamma_law_index - 1.0);
  double H = d + u + p;

  if (aux) {
    aux->velocity_four_vector[0] = u0;
    aux->velocity_four_vector[1] = u1;
    aux->velocity_four_vector[2] = u2;
    aux->velocity_four_vector[3] = u3;
    aux->magnetic_four_vector[0] = 0.0;
    aux->magnetic_four_vector[1] = 0.0;
    aux->magnetic_four_vector[2] = 0.0;
    aux->magnetic_four_vector[3] = 0.0;
    aux->momentum_density[0] = H * u0 * u0;
    aux->momentum_density[1] = H * u0 * u1;
    aux->momentum_density[2] = H * u0 * u2;
    aux->momentum_density[3] = H * u0 * u3;
    aux->comoving_mass_density = P->d;
    aux->gas_pressure = P->p;
    aux->magnetic_pressure = 0.0;
    aux->m2 = m2;
  }
  if (U) {
    U[DDD] = dV * (d * u0);
    U[TAU] = dV * (H * u0 * u0 - p - d * u0);
    U[S11] = dV * (H * u0 * u1);
    U[S22] = dV * (H * u0 * u2);
    U[S33] = dV * (H * u0 * u3);
  }

  return 0;
}


int srhyd_from_conserved(m2sim *m2, double *U, double *B, double *X, double dV,
			 m2aux *aux, m2prim *P)
/*
 * This function transforms the conserved variables into a full state
 * vector. `U` has dimensions of mass, not mass density.
 *
 * The equations are based directly on Marti & Muller section 9.2. The gas
 * pressure is used as the independent variable, from which the rest of
 * primitive variables may be derived. That means that a guess value for the
 * pressure must be obtained directly from the conserved variables. This
 * function uses the exact relations
 *
 *            S^2 / D^2 = h^2 * (W^2 - 1)
 *                    p = D * (h*W - 1) - tau
 *
 * which are easy to verify. When h->1, they can be solved for the pressure and
 * evaluated using only the conserved variables.
 *
 *                    p = D * (sqrt(1 + S^2/D^2) - 1) - tau
 */
{
  static const double ERROR_TOLR = 1e-8;
  static const int NEWTON_MAX_ITER = 50;

  const double D    =  U[DDD] / dV;
  const double Tau  =  U[TAU] / dV;
  const double S1   =  U[S11] / dV;
  const double S2   =  U[S22] / dV;
  const double S3   =  U[S33] / dV;
  const double SS   =  S1*S1 + S2*S2 + S3*S3;

  /* This expression for p becomes exact when h->1 */
  double p = D * (sqrt(1.0 + SS/(D*D)) - 1.0) - Tau;
  double f = 1.0, g = 1.0;
  int n_iter = 0;

  while (fabs(f) > ERROR_TOLR) {

    const double vv   =  SS / pow(Tau + D + p, 2);
    const double W2   =  1.0 / (1.0 - vv);
    const double W    =  sqrt(W2);
    const double Rho  =  D / W;
    const double e    =  (Tau + D*(1.0 - W) + p*(1.0 - W2)) / (D*W);
    const double h    =  1.0 + e + p/Rho;
    const double cs2  =  gamma_law_index * p / (Rho * h);

    /*
     * The Newton-Rapheson step is contained here:
     * 
     * f := function whose zero is being sought
     * g := derivative of f
     * 
     * The expression for g is Marti & Muller Eqn 9.2.78 and is approximate, but
     * it goes asymptotically to the exact derivative as f->0. The iteration
     * step is p_{n+1} = p_{n} - f/g.
     */

    f = Rho * e * (gamma_law_index - 1.0) - p;
    g = vv*cs2 - 1.0; /* This is Marti & Muller Eqn 9.2.78 */
    p -= f/g;
    
    if (n_iter++ == NEWTON_MAX_ITER) {
      MSGF(FATAL,
	   "root-finder failed to locate a root: "
	   "U = %+8.6e %+8.6e %+8.6e %+8.6e %+8.6e",
	   U[DDD], U[TAU], U[S11], U[S22], U[S33]);
    }
  }

  if (p < 0.0) {
    MSG(FATAL, "got a negative pressure");
  }
  if (D < 0.0) {
    MSG(FATAL, "got a negative density");
  }

  const double vv  =  SS / pow(Tau + D + p, 2);
  const double W2  =  1.0 / (1.0 - vv);
  const double u0  =  sqrt(W2);
  const double u1  =  u0 * S1 / (Tau + D + p);
  const double u2  =  u0 * S2 / (Tau + D + p);
  const double u3  =  u0 * S3 / (Tau + D + p);
  const double d   =  D / u0;

  if (aux) {
    aux->velocity_four_vector[0] = u0;
    aux->velocity_four_vector[1] = u1;
    aux->velocity_four_vector[2] = u2;
    aux->velocity_four_vector[3] = u3;
    aux->magnetic_four_vector[0] = 0.0;
    aux->magnetic_four_vector[1] = 0.0;
    aux->magnetic_four_vector[2] = 0.0;
    aux->magnetic_four_vector[3] = 0.0;
    aux->momentum_density[0] = Tau + D;
    aux->momentum_density[1] = S1;
    aux->momentum_density[2] = S2;
    aux->momentum_density[3] = S3;
    aux->comoving_mass_density = d;
    aux->gas_pressure = p;
    aux->magnetic_pressure = 0.0;
    aux->m2 = m2;
  }

  if (P) {
    P->v1 = u1 / u0;
    P->v2 = u2 / u0;
    P->v3 = u3 / u0;
    P->B1 = 0.0;
    P->B2 = 0.0;
    P->B3 = 0.0;
    P->d = d;
    P->p = p;
  }

  return 0;
}


int srhyd_from_auxiliary(m2sim *m2, m2aux *aux, double *X, double dV,
			 m2prim *P, double *U)
{
  return 0;
}


int srhyd_eigenvalues(m2aux *aux, double n[4], double *evals)
{
  const double d = aux->comoving_mass_density;
  const double p = aux->gas_pressure;
  const double u = aux->gas_pressure / (gamma_law_index - 1.0);
  const double u0 = aux->velocity_four_vector[0];
  const double v1 = aux->velocity_four_vector[1] / u0;
  const double v2 = aux->velocity_four_vector[2] / u0;
  const double v3 = aux->velocity_four_vector[3] / u0;
  const double vv = v1*v1 + v2*v2 + v3*v3;
  const double vn = v1*n[1] + v2*n[2] + v3*n[3];
  const double vn2 = vn*vn;
  const double cs2 = gamma_law_index * p / (d + u + p);
  evals[0] = (vn*(1-cs2) - sqrt(cs2*(1-vv)*(1-vv*cs2-vn2*(1-cs2))))/(1-vv*cs2);
  evals[1] = vn;
  evals[2] = vn;
  evals[3] = vn;
  evals[4] = vn;
  evals[5] = vn;
  evals[6] = vn;
  evals[7] = (vn*(1-cs2) + sqrt(cs2*(1-vv)*(1-vv*cs2-vn2*(1-cs2))))/(1-vv*cs2);
  return 0;
}

double srhyd_measure(m2aux *aux, int flag)
{
  return 0.0;
}
