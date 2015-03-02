#include <string.h>
#include "m2.h"
#include "hydro.h"
#include "riemann.h"

/* #define RIEMANN_DEBUG */

static void _hllc_nrhyd(m2riemann_problem *R);
static void _hllc_nrmhd(m2riemann_problem *R);
static void _hlld_srmhd(m2riemann_problem *R);

void riemann_solver_hlle(m2riemann_problem *R)
{
  int q;
  double *F_hll = R->F_hll;
  double *U_hll = R->U_hll;
  double Sl = R->Sl;
  double Sr = R->Sr;
  double *Ul = R->Ul;
  double *Ur = R->Ur;
  double *Fl = R->Fl;
  double *Fr = R->Fr;

  for (q=0; q<8; ++q) {
    U_hll[q] = (Sr*Ur[q] - Sl*Ul[q] +       (Fl[q] - Fr[q])) / (Sr - Sl);
    F_hll[q] = (Sr*Fl[q] - Sl*Fr[q] + Sr*Sl*(Ur[q] - Ul[q])) / (Sr - Sl);
  }

  memcpy(R->F_hat, F_hll, 8 * sizeof(double));
}

void riemann_solver_hllc(m2riemann_problem *R)
{
  if (R->m2->relativistic) {
    MSG(FATAL, "relativistic HLLC Riemann solver not implemented");
  }
  else if (R->m2->magnetized) {
    _hllc_nrmhd(R);
  }
  else {
    _hllc_nrhyd(R);
  }
}

void riemann_solver_hlld(m2riemann_problem *R)
{
  if (R->m2->relativistic && R->m2->magnetized) {
    _hlld_srmhd(R);
  }
  else {
    MSG(FATAL, "HLLD Riemann solver only implemented for srmhd");
  }
}






void _hllc_nrhyd(m2riemann_problem *R)
/*------------------------------------------------------------------------------
 *
 * PURPOSE: Obtain the intercell flux between two constant states using a 3-wave
 *          approximation.
 *
 * REFERENCES: Toro (1999) Riemann Solvers and Numerical Methods for Fluid
 *             Dynamics
 *
 *------------------------------------------------------------------------------
 */
{
  int q;
  double *F_hll = R->F_hll;
  double *U_hll = R->U_hll;
  double Sl = R->Sl;
  double Sr = R->Sr;
  double *Ul = R->Ul;
  double *Ur = R->Ur;
  double *Fl = R->Fl;
  double *Fr = R->Fr;
  double v1L = R->Pl.v1;
  double v2L = R->Pl.v2;
  double v3L = R->Pl.v3;
  double v1R = R->Pr.v1;
  double v2R = R->Pr.v2;
  double v3R = R->Pr.v3;
  double pgL = R->Pl.p;
  double pgR = R->Pr.p;
  double dgL = R->Pl.d;
  double dgR = R->Pr.d;
  double *n = R->nhat;
  double qL = v1L*n[1] + v2L*n[2] + v3L*n[3];
  double qR = v1R*n[1] + v2R*n[2] + v3R*n[3];
  double Sm = (dgR*qR*(Sr - qR) -
	       dgL*qL*(Sl - qL) + (pgL - pgR)) / (dgR*(Sr - qR) - dgL*(Sl - qL));
  double dgl = dgL * (Sl - qL) / (Sl - Sm);
  double dgr = dgR * (Sr - qR) / (Sr - Sm);
  double vparL[4] = {0.0,
		     v1L - qL * n[1],
		     v2L - qL * n[2],
		     v3L - qL * n[3]};
  double vparR[4] = {0.0,
		     v1R - qR * n[1],
		     v2R - qR * n[2],
		     v3R - qR * n[3]};

  R->Uls[TAU] = dgl * (Ul[TAU]/dgL + (Sm - qL)*(Sm + pgL/(dgL*(Sl - qL))));
  R->Uls[S11] = dgl * (Sm * n[1] + vparL[1]);
  R->Uls[S22] = dgl * (Sm * n[2] + vparL[2]);
  R->Uls[S33] = dgl * (Sm * n[3] + vparL[3]);
  R->Uls[DDD] = dgl;

  R->Urs[TAU] = dgr * (Ur[TAU]/dgR + (Sm - qR)*(Sm + pgR/(dgR*(Sr - qR))));
  R->Urs[S11] = dgr * (Sm * n[1] + vparR[1]);
  R->Urs[S22] = dgr * (Sm * n[2] + vparR[2]);
  R->Urs[S33] = dgr * (Sm * n[3] + vparR[3]);
  R->Urs[DDD] = dgr;

  double *U = R->U_hat;
  double *F = R->F_hat;
  double *Uls = R->Uls;
  double *Urs = R->Urs;
  double s = 0.0; /* x/t */

  for (q=0; q<8; ++q) {
    U_hll[q] = (Sr*Ur[q] - Sl*Ul[q] +       (Fl[q] - Fr[q])) / (Sr - Sl);
    F_hll[q] = (Sr*Fl[q] - Sl*Fr[q] + Sr*Sl*(Ur[q] - Ul[q])) / (Sr - Sl);
  }

  if      (        s<=Sl) for (q=0; q<8; ++q) U[q] = Ul [q];
  else if (Sl<s && s<=Sm) for (q=0; q<8; ++q) U[q] = Uls[q];
  else if (Sm<s && s<=Sr) for (q=0; q<8; ++q) U[q] = Urs[q];
  else if (Sr<s         ) for (q=0; q<8; ++q) U[q] = Ur [q];

  if      (        s<=Sl) for (q=0; q<8; ++q) F[q] = Fl[q];
  else if (Sl<s && s<=Sm) for (q=0; q<8; ++q) F[q] = Fl[q] + Sl*(Uls[q]-Ul[q]);
  else if (Sm<s && s<=Sr) for (q=0; q<8; ++q) F[q] = Fr[q] + Sr*(Urs[q]-Ur[q]);
  else if (Sr<s         ) for (q=0; q<8; ++q) F[q] = Fr[q];
}

void _hllc_nrmhd(m2riemann_problem *R)
/*------------------------------------------------------------------------------
 *
 * PURPOSE: Obtain the intercell flux between two constant states using a 3-wave
 *          approximation.
 *
 * REFERENCES: S. Li - Journal of Computational Physics 203 (2005) 344-257
 *
 *------------------------------------------------------------------------------
 */
{
  int q;
  double *F_hll = R->F_hll;
  double *U_hll = R->U_hll;
  double Sl = R->Sl;
  double Sr = R->Sr;
  double *Ul = R->Ul;
  double *Ur = R->Ur;
  double *Fl = R->Fl;
  double *Fr = R->Fr;
  double *U = R->U_hat;
  double *F = R->F_hat;
  double *n = R->nhat;
  double *Uls = R->Uls;
  double *Urs = R->Urs;
  double v1L = R->Pl.v1;
  double v2L = R->Pl.v2;
  double v3L = R->Pl.v3;
  double v1R = R->Pr.v1;
  double v2R = R->Pr.v2;
  double v3R = R->Pr.v3;
  double pL = R->Al.gas_pressure + R->Al.magnetic_pressure;
  double pR = R->Ar.gas_pressure + R->Ar.magnetic_pressure;
  double dL = R->Pl.d;
  double dR = R->Pr.d;
  double qL = v1L*n[1] + v2L*n[2] + v3L*n[3];
  double qR = v1R*n[1] + v2R*n[2] + v3R*n[3];
  double Sm = (dR*qR*(Sr - qR) -
	       dL*qL*(Sl - qL) + (pL - pR)) / (dR*(Sr - qR) - dL*(Sl - qL));
  double dl = dL * (Sl - qL) / (Sl - Sm);
  double dr = dR * (Sr - qR) / (Sr - Sm);
  double nt[4] = { 1.0 - n[0], 1.0 - n[1], 1.0 - n[2], 1.0 - n[3] };
  double B1L = R->Ul[B11];
  double B2L = R->Ul[B22];
  double B3L = R->Ul[B33];
  double B1R = R->Ur[B11];
  double B2R = R->Ur[B22];
  double B3R = R->Ur[B33];
  double BnL = B1L*n[1] + B2L*n[2] + B3L*n[3];
  double BnR = B1R*n[1] + B2R*n[2] + B3R*n[3];

  for (q=0; q<8; ++q) {
    U_hll[q] = (Sr*Ur[q] - Sl*Ul[q] +       (Fl[q] - Fr[q])) / (Sr - Sl);
    F_hll[q] = (Sr*Fl[q] - Sl*Fr[q] + Sr*Sl*(Ur[q] - Ul[q])) / (Sr - Sl);
  }
  nrmhd_from_conserved(R->m2, U_hll, &U_hll[B11], NULL, 1.0, &R->As, NULL);

  double u1L = R->Al.velocity_four_vector[1];
  double u2L = R->Al.velocity_four_vector[2];
  double u3L = R->Al.velocity_four_vector[3];
  double u1R = R->Ar.velocity_four_vector[1];
  double u2R = R->Ar.velocity_four_vector[2];
  double u3R = R->Ar.velocity_four_vector[3];
  double u1s = R->As.velocity_four_vector[1];
  double u2s = R->As.velocity_four_vector[2];
  double u3s = R->As.velocity_four_vector[3];
  double B1s = U_hll[B11];
  double B2s = U_hll[B22];
  double B3s = U_hll[B33];
  double Bns = B1s*n[1] + B2s*n[2] + B3s*n[3];
  double BuL = B1L*u1L + B2L*u2L + B3L*u3L;
  double BuR = B1R*u1R + B2R*u2R + B3R*u3R;
  double Bus = B1s*u1s + B2s*u2s + B3s*u3s;
  double ps = dL * (Sl - qL)*(Sm - qL) + pL; /* same for right/left */
  double fL = (ps*Sm - pL*qL) / (Sl - qL);
  double fR = (ps*Sm - pR*qR) / (Sr - qR);
  double s = 0.0; /* x/t */

  Uls[TAU] = ((fL        +Ul[TAU]*nt[0])*(Sl-qL)-(Bns*Bus-BnL*BuL))/(Sl-Sm);
  Uls[S11] = ((dL*Sm*n[1]+Ul[S11]*nt[1])*(Sl-qL)-(Bns*B1s-BnL*B1L))/(Sl-Sm);
  Uls[S22] = ((dL*Sm*n[2]+Ul[S22]*nt[2])*(Sl-qL)-(Bns*B2s-BnL*B2L))/(Sl-Sm);
  Uls[S33] = ((dL*Sm*n[3]+Ul[S33]*nt[3])*(Sl-qL)-(Bns*B3s-BnL*B3L))/(Sl-Sm);
  Uls[DDD] = dl;
  Uls[B11] = B1s;
  Uls[B22] = B2s;
  Uls[B33] = B3s;

  Urs[TAU] = ((fR        +Ur[TAU]*nt[0])*(Sr-qR)-(Bns*Bus-BnR*BuR))/(Sr-Sm);
  Urs[S11] = ((dR*Sm*n[1]+Ur[S11]*nt[1])*(Sr-qR)-(Bns*B1s-BnR*B1R))/(Sr-Sm);
  Urs[S22] = ((dR*Sm*n[2]+Ur[S22]*nt[2])*(Sr-qR)-(Bns*B2s-BnR*B2R))/(Sr-Sm);
  Urs[S33] = ((dR*Sm*n[3]+Ur[S33]*nt[3])*(Sr-qR)-(Bns*B3s-BnR*B3R))/(Sr-Sm);
  Urs[DDD] = dr;
  Urs[B11] = B1s;
  Urs[B22] = B2s;
  Urs[B33] = B3s;

#ifdef RIEMANN_DEBUG
  for (q=0; q<8; ++q) {
    double lhs = (Sm-Sl)/(Sr-Sl)*Uls[q] - (Sm-Sr)/(Sr-Sl)*Uls[q];
    double rhs = (Sr*Ur[q] - Sl*Ul[q] - (R->Fr[q] - R->Fl[q])) / (Sr - Sl);
    printf("%f == %f\n", lhs, rhs);
  }
  double psL = dL * (Sl - qL)*(Sm - qL) + pL;
  double psR = dR * (Sr - qR)*(Sm - qR) + pR;
  printf("psL == psR ? %f == %f\n", psL, psR);
#endif

  if      (        s<=Sl) for (q=0; q<8; ++q) U[q] = Ul [q];
  else if (Sl<s && s<=Sm) for (q=0; q<8; ++q) U[q] = Uls[q];
  else if (Sm<s && s<=Sr) for (q=0; q<8; ++q) U[q] = Urs[q];
  else if (Sr<s         ) for (q=0; q<8; ++q) U[q] = Ur [q];

  if      (        s<=Sl) for (q=0; q<8; ++q) F[q] = Fl[q];
  else if (Sl<s && s<=Sm) for (q=0; q<8; ++q) F[q] = Fl[q] + Sl*(Uls[q]-Ul[q]);
  else if (Sm<s && s<=Sr) for (q=0; q<8; ++q) F[q] = Fr[q] + Sr*(Urs[q]-Ur[q]);
  else if (Sr<s         ) for (q=0; q<8; ++q) F[q] = Fr[q];
}

void _hlld_srmhd(m2riemann_problem *R)
{

}

