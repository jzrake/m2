#include <math.h>
#include "hydro.h"
#include "srmhd-c2p.h"

#define gamma_law_index (m2 ? m2->gamma_law_index : 5./3.)


int srmhd_from_primitive(m2sim *m2, m2prim *P, double *B, double *X, double dV,
			 double *U, m2aux *aux)
{
  double B1 = P->B1;
  double B2 = P->B2;
  double B3 = P->B3;
  double v1 = P->v1;
  double v2 = P->v2;
  double v3 = P->v3;
  double u0 = 1.0 / sqrt(1.0 - (v1*v1 + v2*v2 + v3*v3));
  double u1 = u0 * v1;
  double u2 = u0 * v2;
  double u3 = u0 * v3;
  double b0 = B1*u1 + B2*u2 + B3*u3;
  double b1 = (P->B1 + b0 * u1) / u0;
  double b2 = (P->B2 + b0 * u2) / u0;
  double b3 = (P->B3 + b0 * u3) / u0;
  double bb = b1*b1 + b2*b2 + b3*b3 - b0*b0;
  double dg = P->d;
  double pg = P->p;
  double ug = pg / (gamma_law_index - 1.0);
  double pb = 0.5 * bb;
  double ub = 0.5 * bb;
  double H0 = dg + (ug + ub) + (pg + pb);

  if (aux) {
    aux->velocity_four_vector[0] = u0;
    aux->velocity_four_vector[1] = u1;
    aux->velocity_four_vector[2] = u2;
    aux->velocity_four_vector[3] = u3;
    aux->magnetic_four_vector[0] = b0;
    aux->magnetic_four_vector[1] = b1;
    aux->magnetic_four_vector[2] = b2;
    aux->magnetic_four_vector[3] = b3;
    aux->momentum_density[0] = H0 * u0 * u0 - b0 * b0;
    aux->momentum_density[1] = H0 * u0 * u1 - b0 * b1;
    aux->momentum_density[2] = H0 * u0 * u2 - b0 * b2;
    aux->momentum_density[3] = H0 * u0 * u3 - b0 * b3;
    aux->comoving_mass_density = dg;
    aux->gas_pressure = pg;
    aux->magnetic_pressure = pb;
    aux->m2 = m2;
  }
  if (U) {
    U[DDD] = dV * (dg * u0);
    U[TAU] = dV * (H0 * u0 * u0 - b0 * b0 - (pg + pb) - dg * u0);
    U[S11] = dV * (H0 * u0 * u1 - b0 * b1);
    U[S22] = dV * (H0 * u0 * u2 - b0 * b2);
    U[S33] = dV * (H0 * u0 * u3 - b0 * b3);
  }

  return 0;
}


int srmhd_from_conserved(m2sim *m2, double *U, double *B, double *X, double dV,
			 m2aux *aux, m2prim *P)
{
  int return_code = 0;
  double Uin[8];
  double Pin[8];
  double B1 = B[1];
  double B2 = B[2];
  double B3 = B[3];
  int error = 1;
  srmhd_c2p *c2p = srmhd_c2p_new();

  Uin[0] = U[DDD] / dV;
  Uin[1] = U[TAU] / dV;
  Uin[2] = U[S11] / dV;
  Uin[3] = U[S22] / dV;
  Uin[4] = U[S33] / dV;
  Uin[5] = B1;
  Uin[6] = B2;
  Uin[7] = B3;

  /* if user set pressure_floor negative then it's disabled */
  srmhd_c2p_set_pressure_floor(c2p, m2->pressure_floor * U[DDD]);
  srmhd_c2p_set_gamma(c2p, gamma_law_index);
  srmhd_c2p_new_state(c2p, Uin);
  srmhd_c2p_estimate_from_cons(c2p);

  if (error) {
    error = srmhd_c2p_solve_noble1dw(c2p, Pin);
  }
  if (error) {
    error = srmhd_c2p_solve_anton2dzw(c2p, Pin);
    if (error == SRMHD_C2P_SUCCESS) {
      return_code |= 1<<M2_BIT_NEARLY_FAILED_CONSERVED_INVERSION;
      /*MSG(WARNING, "got a success with anton2dzw after noble1dw failed");*/
    }
  }
  if (error) {
    /*MSGF(WARNING, "%s", srmhd_c2p_get_error(c2p, error));*/
    /*m2_print_state(NULL, NULL, U);*/
    srmhd_c2p_del(c2p);
    return 1<<M2_BIT_FAILED_CONSERVED_INVERSION;
  }

  double v1 = Pin[2];
  double v2 = Pin[3];
  double v3 = Pin[4];
  double u0 = 1.0 / sqrt(1.0 - (v1*v1 + v2*v2 + v3*v3));
  double u1 = u0 * v1;
  double u2 = u0 * v2;
  double u3 = u0 * v3;
  double b0 = B1*u1 + B2*u2 + B3*u3;
  double b1 = (B1 + b0 * u1) / u0;
  double b2 = (B2 + b0 * u2) / u0;
  double b3 = (B3 + b0 * u3) / u0;
  double bb = -b0*b0 + b1*b1 + b2*b2 + b3*b3;
  double dg = Pin[0];
  double pg = Pin[1];
  double ug = pg / (gamma_law_index - 1.0);
  double pb = 0.5 * bb;
  double ub = 0.5 * bb;
  double H0 = dg + (ug + ub) + (pg + pb);

  /* in case a pressure floor was used, the input conserved state is made
     consistent with it */
  if (srmhd_c2p_put_pressure_floor(c2p)) {
    U[DDD] = dV * (dg * u0);
    U[TAU] = dV * (H0 * u0 * u0 - b0 * b0 - (pg + pb) - dg * u0);
    U[S11] = dV * (H0 * u0 * u1 - b0 * b1);
    U[S22] = dV * (H0 * u0 * u2 - b0 * b2);
    U[S33] = dV * (H0 * u0 * u3 - b0 * b3);
  }

  if (aux) {
    aux->velocity_four_vector[0] = u0;
    aux->velocity_four_vector[1] = u1;
    aux->velocity_four_vector[2] = u2;
    aux->velocity_four_vector[3] = u3;
    aux->magnetic_four_vector[0] = b0;
    aux->magnetic_four_vector[1] = b1;
    aux->magnetic_four_vector[2] = b2;
    aux->magnetic_four_vector[3] = b3;
    aux->momentum_density[0] = H0 * u0 * u0 - b0 * b0;
    aux->momentum_density[1] = H0 * u0 * u1 - b0 * b1;
    aux->momentum_density[2] = H0 * u0 * u2 - b0 * b2;
    aux->momentum_density[3] = H0 * u0 * u3 - b0 * b3;
    aux->comoving_mass_density = dg;
    aux->gas_pressure = pg;
    aux->magnetic_pressure = pb;
    aux->m2 = m2;
  }
  if (P) {
    P->v1 = v1;
    P->v2 = v2;
    P->v3 = v3;
    P->B1 = B1;
    P->B2 = B2;
    P->B3 = B3;
    P->d = dg;
    P->p = pg;
  }

  if (srmhd_c2p_put_pressure_floor(c2p)) {
    return_code |= 1<<M2_BIT_NEGATIVE_PRESSURE;
  }

  srmhd_c2p_del(c2p);
  return return_code;
}


int srmhd_from_auxiliary(m2sim *m2, m2aux *aux, double *X, double dV,
			 m2prim *P, double *U)
{
  return 0;
}


int srmhd_eigenvalues(m2aux *aux, double n[4], double *evals)
{
  m2sim *m2 = aux->m2;
  int quartic_solver = m2 ? m2->quartic_solver : M2_QUARTIC_FULL_COMPLEX;
  if (quartic_solver == M2_QUARTIC_NONE) {
    evals[0] = -1.0;
    evals[1] =  0.0;
    evals[2] =  0.0;
    evals[3] =  0.0;
    evals[4] =  0.0;
    evals[5] =  0.0;
    evals[6] =  0.0;
    evals[7] = +1.0;
    return 0;
  }
  const double n1 = n[1];
  const double n2 = n[2];
  const double n3 = n[3];
  const double u0 = aux->velocity_four_vector[0];
  const double v1 = aux->velocity_four_vector[1] / u0;
  const double v2 = aux->velocity_four_vector[2] / u0;
  const double v3 = aux->velocity_four_vector[3] / u0;
  const double b0 = aux->magnetic_four_vector[0];
  const double b1 = aux->magnetic_four_vector[1];
  const double b2 = aux->magnetic_four_vector[2];
  const double b3 = aux->magnetic_four_vector[3];
  const double bb = b1*b1 + b2*b2 + b3*b3 - b0*b0;
  const double bn = b1*n1 + b2*n2 + b3*n3;
  const double vn = v1*n1 + v2*n2 + v3*n3;
  const double dg = aux->comoving_mass_density;
  const double pg = aux->gas_pressure;
  const double ug = pg / (gamma_law_index - 1.0);
  const double Hg = dg + ug + pg; /* gas enthalpy density */
  const double cs2 = gamma_law_index * pg / Hg; /* sound speed squared */
  const double C = Hg + bb; /* constant in Alfven wave expression */

  const double W2 = u0*u0;
  const double W4 = W2*W2;
  const double V2 = vn*vn; /* Vx := vn^x */
  const double V3 = vn*V2;
  const double V4 = vn*V3;
  const double K  =    Hg * (1.0/cs2-1)  * W4;
  const double L  =  -(Hg +   bb/cs2)    * W2;
  const double A4 =    K    - L          -   b0*b0;
  const double A3 = -4*K*vn + L*vn*2     + 2*b0*bn;
  const double A2 =  6*K*V2 + L*(1.0-V2) +   b0*b0 - bn*bn;
  const double A1 = -4*K*V3 - L*vn*2     - 2*b0*bn;
  const double A0 =    K*V4 + L*V2       +   bn*bn;

  double roots[4];
  int nr = 0;

  switch (quartic_solver) {
  case M2_QUARTIC_ALGORITHMIC:
    nr = m2_solve_quartic_equation1(A4, A3, A2, A1, A0, roots);
    break;
  case M2_QUARTIC_FULL_COMPLEX:
    nr = m2_solve_quartic_equation2(A4, A3, A2, A1, A0, roots);
    break;
  }

  if (nr == 4) {
    evals[0] = roots[0];
    evals[1] = (bn - sqrt(C) * vn * u0) / (b0 - sqrt(C) * u0);
    evals[2] = roots[1];
    evals[3] = vn;
    evals[4] = vn;
    evals[5] = roots[2];
    evals[6] = (bn + sqrt(C) * vn * u0) / (b0 + sqrt(C) * u0);
    evals[7] = roots[3];
    return 0;
  }
  else {
    /*m2_print_state(NULL, aux, NULL);*/
    /*MSGF(INFO, "A = [%+12.10f %+12.10f %+12.10f %+12.10f %+12.10f]",*/
        /*A0, A1, A2, A3, A4);*/
    /*MSGF(INFO, "x = [%+12.10f %+12.10f %+12.10f %+12.10f]",*/
        /*roots[0], roots[1], roots[2], roots[3]);*/
    /*MSGF(INFO, "magnetosonic polynomial N4=0 has %d real roots", nr);*/
    return 1;
  }
}


double srmhd_measure(m2aux *aux, int flag)
{
  m2sim *m2 = aux->m2;
  double u0 = aux->velocity_four_vector[0];
  double v1 = aux->velocity_four_vector[1] / u0;
  double v2 = aux->velocity_four_vector[2] / u0;
  double v3 = aux->velocity_four_vector[3] / u0;
  double dg = aux->comoving_mass_density;
  double pg = aux->gas_pressure;
  double pb = aux->magnetic_pressure;
  double vv = v1*v1 + v2*v2 + v3*v3;
  double v0 = sqrt(vv);
  double ug = pg / (gamma_law_index - 1.0);
  double Hg = dg + ug + pg; /* gas enthalpy density */
  double eg = Hg*u0*u0 - pg - u0*dg; /* gas energy density */
  double ek = dg*u0*(u0 - 1.0);
  double cs = sqrt(gamma_law_index * pg / Hg);

  switch (flag) {
  case M2_SIGMA: return pb / eg;
  case M2_SOUND_SPEED: return cs;
  case M2_MACH_NUMBER: return (v0/sqrt(1.0 - v0*v0)) / (cs/sqrt(1.0 - cs*cs));
  case M2_OBSERVER_MASS_DENSITY: return dg * u0;
  case M2_TOTAL_ENERGY_DENSITY: return aux->momentum_density[0] - (pg + pb) - dg * u0;
  case M2_INTERNAL_ENERGY_DENSITY: return ug;
  case M2_KINETIC_ENERGY_DENSITY: return ek;
  case M2_MAGNETIC_ENERGY_DENSITY: return pb;
  default:
    return 0.0;
  }
}
