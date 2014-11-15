#ifndef M2_RIEMANN_HEADER
#define M2_RIEMANN_HEADER
#include "m2.h"

typedef struct
{
  m2sim *m2;
  double nhat[4];
  double Fl[8];
  double Fr[8];
  double Ul[8];
  double Ur[8];
  double lamL[8];
  double lamR[8];
  m2prim Pl;
  m2prim Pr;
  m2aux Al;
  m2aux Ar;
  m2aux As;
  double Uls[8]; /* U left star */
  double Urs[8]; /* U right star */
  double F_hll[8];
  double U_hll[8];
  double F_hat[8];
  double U_hat[8];
  double Sl; /* left-most wavespeed */
  double Sr; /* right-most wavespeed */
  double Sc; /* contact wavespeed */
} m2riemann_problem;

void riemann_solver_hlle(m2riemann_problem *R);
void riemann_solver_hllc(m2riemann_problem *R);

#endif /* M2_RIEMANN_HEADER */
