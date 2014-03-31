#include <math.h>
#include "hydro.h"


#define SYSTEM_CHOICE(m2,E)						\
  do {									\
    if (0) { }								\
    else if (m2->relativistic == 1 && m2->magnetized == 1) {		\
      return srmhd_##E;							\
    }									\
    else if (m2->relativistic == 1 && m2->magnetized == 0) {		\
      return srhyd_##E;							\
    }									\
    else if (m2->relativistic == 0 && m2->magnetized == 1) {		\
      return nrmhd_##E;							\
    }									\
    else if (m2->relativistic == 0 && m2->magnetized == 0) {		\
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
double m2aux_measure(m2aux *aux, int flag)
{
  SYSTEM_CHOICE(aux->m2, measure(aux, flag));
}
int m2vol_from_primitive(m2vol *V)
{
  m2sim *m2 = V->m2;
  m2prim *P = &V->prim;
  m2aux  *A = &V->aux;
  double *B = &V->prim.B1 - 1;
  double *U = V->consA;
  double dV = V->volume;
  SYSTEM_CHOICE(m2, from_primitive(m2, P, B, NULL, dV, U, A));
}
int m2aux_add_geometrical_source_terms(m2aux *aux, double x0[4], double x1[4],
				       double *U)
{
  if (aux->m2->geometry == M2_CARTESIAN) {
    return 0;
  }
  else if (aux->m2->geometry == M2_CYLINDRICAL) {
    const double dt = x1[0] - x0[0];
    const double r0 = x0[1];
    const double r1 = x1[1];
    const double f0 = x0[2];
    const double f1 = x1[2];
    const double z0 = x0[3];
    const double z1 = x1[3];
    const double S1 = aux->momentum_density[1];
    const double S2 = aux->momentum_density[2];
    const double u0 = aux->velocity_four_vector[0];
    const double u1 = aux->velocity_four_vector[1];
    const double u2 = aux->velocity_four_vector[2];
    const double b0 = aux->magnetic_four_vector[0];
    const double b1 = aux->magnetic_four_vector[1];
    const double b2 = aux->magnetic_four_vector[2];
    const double v2 = aux->velocity_four_vector[2] / u0;
    const double B1 = b1 * u0 - b0 * u1;
    const double B2 = b2 * u0 - b0 * u2;
    const double p  = aux->gas_pressure + aux->magnetic_pressure;
    U[S11] += dt * ((-f0 + f1)*(-r0 + r1)*
		    (-(b2*B2) + u0*(p + S2*v2))*(-z0 + z1))/u0;
    U[S22] += dt * ((f0 - f1)*(r0 - r1)*(-(B1*b2) + S1*u0*v2)*(z0 - z1))/u0;
    U[S33] += dt * 0.0;
    return 0;
  }
  else if (aux->m2->geometry == M2_SPHERICAL) {
    const double dt = x1[0] - x0[0];
    const double r0 = x0[1];
    const double r1 = x1[1];
    const double t0 = x0[2];
    const double t1 = x1[2];
    const double f0 = x0[3];
    const double f1 = x1[3];
    const double S1 = aux->momentum_density[1];
    const double S2 = aux->momentum_density[2];
    const double S3 = aux->momentum_density[3];
    const double u0 = aux->velocity_four_vector[0];
    const double u1 = aux->velocity_four_vector[1];
    const double u2 = aux->velocity_four_vector[2];
    const double u3 = aux->velocity_four_vector[3];
    const double b0 = aux->magnetic_four_vector[0];
    const double b1 = aux->magnetic_four_vector[1];
    const double b2 = aux->magnetic_four_vector[2];
    const double b3 = aux->magnetic_four_vector[3];
    const double v2 = aux->velocity_four_vector[2] / u0;
    const double v3 = aux->velocity_four_vector[3] / u0;
    const double B1 = b1 * u0 - b0 * u1;
    const double B2 = b2 * u0 - b0 * u2;
    const double B3 = b3 * u0 - b0 * u3;
    const double p  = aux->gas_pressure + aux->magnetic_pressure;
    /* ------------------------------------ */
    /* Mathematica made the following mess: */
    /* ------------------------------------ */
    U[S11] += dt * ((B2*b2 + B3*b3 - (2*p + S2*v2 + S3*v3)*u0)*(f0 - f1)*
		    (cos(t0) - cos(t1))*(-(r0*r0) + r1*r1))/(2.*u0);
    U[S22] += dt * ((r0 - r1)*(r0 + r1)*(f0 - f1)*
		    ((B1*b2 - S1*v2*u0)*cos(t0) +
		     (-(B1*b2) + S1*v2*u0)*cos(t1) + (B3*b3 - (p + S3*v3)*u0)*
		     (sin(t0) - sin(t1))))/(2.*u0);
    U[S33] += dt * ((f0 - f1)*((-(B1*b3) + S1*v3*u0)*cos(t0) +
			       (B1*b3 - S1*v3*u0)*cos(t1) +
			       (B2*b3 - S2*v3*u0)*(sin(t0) - sin(t1)))*
		    (-(r0*r0) + r1*r1))/(2.*u0);
    return 0;
  }
  else if (aux->m2->geometry == M2_PARABOLIC) {
    const double dt = x1[0] - x0[0];
    const double s0 = x0[1];
    const double s1 = x1[1];
    const double t0 = x0[2];
    const double t1 = x1[2];
    const double f0 = x0[3];
    const double f1 = x1[3];
    const double S1 = aux->momentum_density[1];
    const double S2 = aux->momentum_density[2];
    const double S3 = aux->momentum_density[3];
    const double u0 = aux->velocity_four_vector[0];
    const double u1 = aux->velocity_four_vector[1];
    const double u2 = aux->velocity_four_vector[2];
    const double u3 = aux->velocity_four_vector[3];
    const double b0 = aux->magnetic_four_vector[0];
    const double b1 = aux->magnetic_four_vector[1];
    const double b2 = aux->magnetic_four_vector[2];
    const double b3 = aux->magnetic_four_vector[3];
    const double v1 = aux->velocity_four_vector[1] / u0;
    const double v2 = aux->velocity_four_vector[2] / u0;
    const double v3 = aux->velocity_four_vector[3] / u0;
    const double B1 = b1 * u0 - b0 * u1;
    const double B2 = b2 * u0 - b0 * u2;
    const double B3 = b3 * u0 - b0 * u3;
    const double p  = aux->gas_pressure + aux->magnetic_pressure;
    U[S11] += dt * 
(((-f0 + f1)*(sqrt(s0*s0 + t0*t0)*(3*t0*(b1*B2 - S2*u0*v1)*(s0*s0) + 
s0*(-3*b2*B2 - 5*b3*B3 + u0*(8*p + 3*S2*v2 + 5*S3*v3))*(t0*t0) - 
2*(3*b2*B2 + b3*B3 - u0*(4*p + 3*S2*v2 + S3*v3))*(s0*s0*s0) + 
6*(b1*B2 - S2*u0*v1)*(t0*t0*t0)) + sqrt(s1*s1 + 
t0*t0)*(3*t0*(-(b1*B2) + S2*u0*v1)*(s1*s1) + s1*(3*b2*B2 + 5*b3*B3 - 
u0*(8*p + 3*S2*v2 + 5*S3*v3))*(t0*t0) + 2*(3*b2*B2 + b3*B3 - u0*(4*p 
+ 3*S2*v2 + S3*v3))*(s1*s1*s1) + 6*(-(b1*B2) + S2*u0*v1)*(t0*t0*t0)) 
+ sqrt(s1*s1 + t1*t1)*(3*t1*(b1*B2 - S2*u0*v1)*(s1*s1) + s1*(-3*b2*B2 
- 5*b3*B3 + u0*(8*p + 3*S2*v2 + 5*S3*v3))*(t1*t1) - 2*(3*b2*B2 + 
b3*B3 - u0*(4*p + 3*S2*v2 + S3*v3))*(s1*s1*s1) + 6*(b1*B2 - 
S2*u0*v1)*(t1*t1*t1)) + sqrt(s0*s0 + t1*t1)*(3*t1*(-(b1*B2) + 
S2*u0*v1)*(s0*s0) + s0*(3*b2*B2 + 5*b3*B3 - u0*(8*p + 3*S2*v2 + 
5*S3*v3))*(t1*t1) + 2*(3*b2*B2 + b3*B3 - u0*(4*p + 3*S2*v2 + 
S3*v3))*(s0*s0*s0) + 6*(-(b1*B2) + S2*u0*v1)*(t1*t1*t1)) + 
3*(-(b1*B2) + S2*u0*v1)*log(t0 + sqrt(s0*s0 + t0*t0))*(s0*s0*s0*s0) + 
3*(b1*B2 - S2*u0*v1)*log(t1 + sqrt(s0*s0 + t1*t1))*(s0*s0*s0*s0) + 
3*(b1*B2 - S2*u0*v1)*log(t0 + sqrt(s1*s1 + t0*t0))*(s1*s1*s1*s1) + 
3*(-(b1*B2) + S2*u0*v1)*log(t1 + sqrt(s1*s1 + t1*t1))*(s1*s1*s1*s1) + 
3*(b2*B2 - b3*B3 - S2*u0*v2 + S3*u0*v3)*log(s0 + sqrt(s0*s0 + 
t0*t0))*(t0*t0*t0*t0) + 3*(-(b2*B2) + b3*B3 + S2*u0*v2 - 
S3*u0*v3)*log(s1 + sqrt(s1*s1 + t0*t0))*(t0*t0*t0*t0) + 3*(-(b2*B2) + 
b3*B3 + S2*u0*v2 - S3*u0*v3)*log(s0 + sqrt(s0*s0 + 
t1*t1))*(t1*t1*t1*t1) + 3*(b2*B2 - b3*B3 - S2*u0*v2 + 
S3*u0*v3)*log(s1 + sqrt(s1*s1 + t1*t1))*(t1*t1*t1*t1)))/(24.*u0));

    U[S22] += dt * 
(((f0 - f1)*(5*b3*B3*t0*(s0*s0)*sqrt(s0*s0 + t0*t0) - 
8*p*t0*u0*(s0*s0)*sqrt(s0*s0 + t0*t0) - 
5*S3*t0*u0*v3*(s0*s0)*sqrt(s0*s0 + t0*t0) - 
5*b3*B3*t0*(s1*s1)*sqrt(s1*s1 + t0*t0) + 8*p*t0*u0*(s1*s1)*sqrt(s1*s1 
+ t0*t0) + 5*S3*t0*u0*v3*(s1*s1)*sqrt(s1*s1 + t0*t0) - 
5*b3*B3*t1*(s0*s0)*sqrt(s0*s0 + t1*t1) + 8*p*t1*u0*(s0*s0)*sqrt(s0*s0 
+ t1*t1) + 5*S3*t1*u0*v3*(s0*s0)*sqrt(s0*s0 + t1*t1) + 
5*b3*B3*t1*(s1*s1)*sqrt(s1*s1 + t1*t1) - 8*p*t1*u0*(s1*s1)*sqrt(s1*s1 
+ t1*t1) - 5*S3*t1*u0*v3*(s1*s1)*sqrt(s1*s1 + t1*t1) + 
2*b3*B3*sqrt(s0*s0 + t0*t0)*(t0*t0*t0) - 8*p*u0*sqrt(s0*s0 + 
t0*t0)*(t0*t0*t0) - 2*S3*u0*v3*sqrt(s0*s0 + t0*t0)*(t0*t0*t0) - 
2*b3*B3*sqrt(s1*s1 + t0*t0)*(t0*t0*t0) + 8*p*u0*sqrt(s1*s1 + 
t0*t0)*(t0*t0*t0) + 2*S3*u0*v3*sqrt(s1*s1 + t0*t0)*(t0*t0*t0) - 
2*b3*B3*sqrt(s0*s0 + t1*t1)*(t1*t1*t1) + 8*p*u0*sqrt(s0*s0 + 
t1*t1)*(t1*t1*t1) + 2*S3*u0*v3*sqrt(s0*s0 + t1*t1)*(t1*t1*t1) + 
2*b3*B3*sqrt(s1*s1 + t1*t1)*(t1*t1*t1) - 8*p*u0*sqrt(s1*s1 + 
t1*t1)*(t1*t1*t1) - 2*S3*u0*v3*sqrt(s1*s1 + t1*t1)*(t1*t1*t1) - 
3*B1*(b2*(s0*(t0*t0*sqrt(s0*s0 + t0*t0) - t1*t1*sqrt(s0*s0 + t1*t1)) 
+ s1*(-(t0*t0*sqrt(s1*s1 + t0*t0)) + sqrt(s1*s1 + t1*t1)*(2*(s1*s1) + 
t1*t1)) + 2*(sqrt(s0*s0 + t0*t0) - sqrt(s0*s0 + t1*t1))*(s0*s0*s0) - 
2*sqrt(s1*s1 + t0*t0)*(s1*s1*s1)) + b1*(t0*(s1*s1)*sqrt(s1*s1 + 
t0*t0) - t1*sqrt(s1*s1 + t1*t1)*(s1*s1 + 2*(t1*t1)) + 
s0*s0*(-(t0*sqrt(s0*s0 + t0*t0)) + t1*sqrt(s0*s0 + t1*t1)) + 
2*(-sqrt(s0*s0 + t0*t0) + sqrt(s1*s1 + t0*t0))*(t0*t0*t0) + 
2*sqrt(s0*s0 + t1*t1)*(t1*t1*t1))) + 
3*S1*u0*(v2*(s0*(t0*t0*sqrt(s0*s0 + t0*t0) - t1*t1*sqrt(s0*s0 + 
t1*t1)) + s1*(-(t0*t0*sqrt(s1*s1 + t0*t0)) + sqrt(s1*s1 + 
t1*t1)*(2*(s1*s1) + t1*t1)) + 2*(sqrt(s0*s0 + t0*t0) - sqrt(s0*s0 + 
t1*t1))*(s0*s0*s0) - 2*sqrt(s1*s1 + t0*t0)*(s1*s1*s1)) + 
v1*(t0*(s1*s1)*sqrt(s1*s1 + t0*t0) - t1*sqrt(s1*s1 + t1*t1)*(s1*s1 + 
2*(t1*t1)) + s0*s0*(-(t0*sqrt(s0*s0 + t0*t0)) + t1*sqrt(s0*s0 + 
t1*t1)) + 2*(-sqrt(s0*s0 + t0*t0) + sqrt(s1*s1 + t0*t0))*(t0*t0*t0) + 
2*sqrt(s0*s0 + t1*t1)*(t1*t1*t1))) + 3*((-(b1*B1) + b3*B3 + S1*u0*v1 
- S3*u0*v3)*log(t0 + sqrt(s0*s0 + t0*t0))*(s0*s0*s0*s0) + 
b1*B1*log(t1 + sqrt(s0*s0 + t1*t1))*(s0*s0*s0*s0) - b3*B3*log(t1 + 
sqrt(s0*s0 + t1*t1))*(s0*s0*s0*s0) - S1*u0*v1*log(t1 + sqrt(s0*s0 + 
t1*t1))*(s0*s0*s0*s0) + S3*u0*v3*log(t1 + sqrt(s0*s0 + 
t1*t1))*(s0*s0*s0*s0) + b1*B1*log(t0 + sqrt(s1*s1 + 
t0*t0))*(s1*s1*s1*s1) - b3*B3*log(t0 + sqrt(s1*s1 + 
t0*t0))*(s1*s1*s1*s1) - S1*u0*v1*log(t0 + sqrt(s1*s1 + 
t0*t0))*(s1*s1*s1*s1) + S3*u0*v3*log(t0 + sqrt(s1*s1 + 
t0*t0))*(s1*s1*s1*s1) + (-(b1*B1) + b3*B3 + S1*u0*v1 - 
S3*u0*v3)*log(t1 + sqrt(s1*s1 + t1*t1))*(s1*s1*s1*s1) + (B1*b2 - 
S1*u0*v2)*log(s0 + sqrt(s0*s0 + t0*t0))*(t0*t0*t0*t0) - B1*b2*log(s1 
+ sqrt(s1*s1 + t0*t0))*(t0*t0*t0*t0) + S1*u0*v2*log(s1 + sqrt(s1*s1 + 
t0*t0))*(t0*t0*t0*t0) - B1*b2*log(s0 + sqrt(s0*s0 + 
t1*t1))*(t1*t1*t1*t1) + S1*u0*v2*log(s0 + sqrt(s0*s0 + 
t1*t1))*(t1*t1*t1*t1) + B1*b2*log(s1 + sqrt(s1*s1 + 
t1*t1))*(t1*t1*t1*t1) - S1*u0*v2*log(s1 + sqrt(s1*s1 + 
t1*t1))*(t1*t1*t1*t1))))/(24.*u0));

      U[S33] += dt * 
(-((f0 - f1)*(S1*u0*v3*(s0*(-5*(t0*t0)*sqrt(s0*s0 + t0*t0) + 
5*(t1*t1)*sqrt(s0*s0 + t1*t1)) + s1*(5*(t0*t0)*sqrt(s1*s1 + t0*t0) - 
5*(t1*t1)*sqrt(s1*s1 + t1*t1) + 2*(s1*s1)*(sqrt(s1*s1 + t0*t0) - 
sqrt(s1*s1 + t1*t1))) + 2*(-sqrt(s0*s0 + t0*t0) + sqrt(s0*s0 + 
t1*t1))*(s0*s0*s0)) + B1*b3*(5*s0*(t0*t0*sqrt(s0*s0 + t0*t0) - 
t1*t1*sqrt(s0*s0 + t1*t1)) + s1*(-5*(t0*t0)*sqrt(s1*s1 + t0*t0) + 
sqrt(s1*s1 + t1*t1)*(2*(s1*s1) + 5*(t1*t1))) + 2*(sqrt(s0*s0 + t0*t0) 
- sqrt(s0*s0 + t1*t1))*(s0*s0*s0) - 2*sqrt(s1*s1 + t0*t0)*(s1*s1*s1)) 
+ (B2*b3 - S2*u0*v3)*(-5*t0*(s1*s1)*sqrt(s1*s1 + t0*t0) + 
t1*sqrt(s1*s1 + t1*t1)*(5*(s1*s1) + 2*(t1*t1)) + 
5*(s0*s0)*(t0*sqrt(s0*s0 + t0*t0) - t1*sqrt(s0*s0 + t1*t1)) + 
2*(sqrt(s0*s0 + t0*t0) - sqrt(s1*s1 + t0*t0))*(t0*t0*t0) - 
2*sqrt(s0*s0 + t1*t1)*(t1*t1*t1)) + 3*((B2*b3 - S2*u0*v3)*log(t0 + 
sqrt(s0*s0 + t0*t0))*(s0*s0*s0*s0) - B2*b3*log(t1 + sqrt(s0*s0 + 
t1*t1))*(s0*s0*s0*s0) + S2*u0*v3*log(t1 + sqrt(s0*s0 + 
t1*t1))*(s0*s0*s0*s0) - B2*b3*log(t0 + sqrt(s1*s1 + 
t0*t0))*(s1*s1*s1*s1) + S2*u0*v3*log(t0 + sqrt(s1*s1 + 
t0*t0))*(s1*s1*s1*s1) + (B2*b3 - S2*u0*v3)*log(t1 + sqrt(s1*s1 + 
t1*t1))*(s1*s1*s1*s1) + (B1*b3 - S1*u0*v3)*log(s0 + sqrt(s0*s0 + 
t0*t0))*(t0*t0*t0*t0) - B1*b3*log(s1 + sqrt(s1*s1 + 
t0*t0))*(t0*t0*t0*t0) + S1*u0*v3*log(s1 + sqrt(s1*s1 + 
t0*t0))*(t0*t0*t0*t0) - B1*b3*log(s0 + sqrt(s0*s0 + 
t1*t1))*(t1*t1*t1*t1) + S1*u0*v3*log(s0 + sqrt(s0*s0 + 
t1*t1))*(t1*t1*t1*t1) + B1*b3*log(s1 + sqrt(s1*s1 + 
t1*t1))*(t1*t1*t1*t1) - S1*u0*v3*log(s1 + sqrt(s1*s1 + 
t1*t1))*(t1*t1*t1*t1))))/(24.*u0));

    return 0;
  }
  else {
    MSG(FATAL, "invalid geometry");
    return 0;
  }
}
int m2aux_fluxes(m2aux *aux, double n[4], double *F)
{
  const double SR = aux->m2->relativistic != 0;
  const double u0 = aux->velocity_four_vector[0];
  const double u1 = aux->velocity_four_vector[1];
  const double u2 = aux->velocity_four_vector[2];
  const double u3 = aux->velocity_four_vector[3];
  const double b0 = aux->magnetic_four_vector[0];
  const double b1 = aux->magnetic_four_vector[1];
  const double b2 = aux->magnetic_four_vector[2];
  const double b3 = aux->magnetic_four_vector[3];
  const double S0 = aux->momentum_density[0]; /* T^{0,0} = tau + D + pg + pb */
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
  const double T0 = S0 - (D0 + pg + pb); /* tau */
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
double m2aux_maximum_wavespeed(m2aux *aux, int *err)
{
  double evals[8];
  double n[4] = { 0.0,
		  aux->velocity_four_vector[1],
		  aux->velocity_four_vector[2],
		  aux->velocity_four_vector[3] };
  double N = sqrt(n[1]*n[1] + n[2]*n[2] + n[3]*n[3]);
  double a;
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
  a = fabs(evals[0]) > fabs(evals[7]) ? fabs(evals[0]) : fabs(evals[7]);
  if (a != a) {
    *err = 1;
    return 0.0;
  }
  *err = 0;
  return a;
}
