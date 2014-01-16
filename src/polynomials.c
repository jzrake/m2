/*------------------------------------------------------------------------------
 * FILE: polynomials.c
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP (adapted from C++ code by Weiqun Zhang)
 *
 * DESCRIPTION: Implementation of exact solution to the quartic polynomial
 *
 * REFERENCES:
 *------------------------------------------------------------------------------
 */

#include <math.h>
#include <stdlib.h>
#include "polynomials.h"

static int compare(const void *x, const void *y)
{
  double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return +1;
  return 0;
}

int solve_quartic_equation(double d4, double d3,
			   double d2, double d1, double d0,
			   double roots[4])
{
  int nr, nr12, nr34;
  double r1, r2, r3, r4;
  double a3 = d3/d4;
  double a2 = d2/d4;
  double a1 = d1/d4;
  double a0 = d0/d4;
  double au2 = -a2;
  double au1 = (a1*a3 - 4.0*a0) ;
  double au0 = 4.0*a0*a2 - a1*a1 - a0*a3*a3;
  double x1, x3, roots3[3];
  double u1;

  nr = solve_cubic_equation(1.0, au2, au1, au0, roots3);
  x1 = roots3[0];
  x3 = roots3[2];

  if (nr == 1) {
    u1 = x1;
  }
  else {
    u1 = (x1 > x3) ? x1 : x3;
  }

  double R2 = 0.25*a3*a3 + u1 - a2;
  double R = (R2>0.0) ? sqrt(R2) : 0.0;
  double D2, E2;

  if (R != 0.0) {
    double foo1 = 0.75*a3*a3 - R2 - 2.0*a2;
    double foo2 = 0.25*(4.0*a3*a2 - 8.0*a1 - a3*a3*a3) / R;
    D2 = foo1 + foo2;
    E2 = foo1 - foo2;
  }
  else {
    double foo1 = 0.75*a3*a3 - 2.0*a2;
    double foo2 = 2.0 * sqrt(u1*u1 - 4.0*a0);
    D2 = foo1 + foo2;
    E2 = foo1 - foo2;
  }

  if (D2 >= 0.0) {
    double D = sqrt(D2);
    r1 = -0.25*a3 + 0.5*R - 0.5*D;
    r2 = -0.25*a3 + 0.5*R + 0.5*D;
    nr12 = 2;
  }
  else {
    r1 = r2 = -0.25*a3 + 0.5*R;
    nr12 = 0;
  }

  if (E2 >= 0.0) {
    double E = sqrt(E2);
    r3 = -0.25*a3 - 0.5*R - 0.5*E;
    r4 = -0.25*a3 - 0.5*R + 0.5*E;
    nr34 = 2;
  }
  else {
    r3 = r4 = -0.25*a3 - 0.5*R;
    nr34 = 0;
  }

  int i = 0;
  if (nr12 != 0) {
    roots[i++] = r1;
    roots[i++] = r2;
  }
  if (nr34 != 0) {
    roots[i++] = r3;
    roots[i++] = r4;
  }

  qsort(roots, nr12 + nr34, sizeof(double), compare);

  return nr12 + nr34;
}


int solve_cubic_equation(double c3, double c2, double c1, double c0,
			 double roots[3])
{
  double a2 = c2/c3;
  double a1 = c1/c3;
  double a0 = c0/c3;
  double q = a1/3.0 - a2*a2/9.0;
  double r = (a1*a2 - 3.0*a0)/6.0 - a2*a2*a2 / 27.0;
  double delta = q*q*q + r*r;
  double *x1 = &roots[0];
  double *x2 = &roots[1];
  double *x3 = &roots[2];
  int nr;

  if (delta>0.0) {
    double s1 = r + sqrt(delta);
    s1 = (s1>=0.0) ? pow(s1,1./3.) : -pow(-s1,1./3.);

    double s2 = r - sqrt(delta);
    s2 = (s2>=0.0) ? pow(s2,1./3.) : -pow(-s2,1./3.);

    *x1 = (s1+s2) - a2/3.0;
    *x2 = *x3 = -0.5 * (s1+s2) - a2/3.0;
    nr = 1;
  }
  else if (delta < 0.0) {
    double theta = acos(r/sqrt(-q*q*q)) / 3.0;
    double costh = cos(theta);
    double sinth = sin(theta);
    double sq = sqrt(-q);

    *x1 = 2.0*sq*costh - a2/3.0;
    *x2 = -sq*costh - a2/3.0 - sqrt(3.) * sq * sinth;
    *x3 = -sq*costh - a2/3.0 + sqrt(3.) * sq * sinth;

    nr = 3;
  }
  else {
    double s = (r>=0.0) ? pow(r,1./3.) : -pow(-r,1./3.);
    *x1 = 2.0*s - a2/3.0;
    *x2 = *x3 = -s - a2/3.0;

    nr = 3;
  }
  return nr;
}
