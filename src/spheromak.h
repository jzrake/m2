#ifndef SPHEROMAK_H
#define SPHEROMAK_H

#include <complex.h>
typedef complex double Complex;


double LegendreP(int n, int m, double x); /* uses spherical normalization */
double SphericalBesselJZero(int N, int L);
double SphericalBesselJ(int L, double r);
double CylindricalBesselJ(int N, double X);
Complex SphericalHarmonicY(int L, int M, double t, double p);


void SpheromakExternal(int N, int L, int M, double r, double t, double p, double *Br, double *Bt, double *Bp);
void SpheromakInternal(int N, int L, int M, double r, double t, double p, double *Br, double *Bt, double *Bp);


#endif /* SPHEROMAK_H */
