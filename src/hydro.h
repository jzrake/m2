#ifndef M2_HYDRO_HEADER
#define M2_HYDRO_HEADER

#include "m2.h"


/* -----------------------------------------------------------------------------
 * NRHYD functions
 * ---------------------------------------------------------------------------*/
int nrhyd_from_primitive(m2sim *m2, m2prim *P, double *B, double *X, double dV,
			 double *U, m2aux *aux);
int nrhyd_from_conserved(m2sim *m2, double *U, double *B, double *X, double dV,
			 m2aux *aux, m2prim *P);
int nrhyd_from_auxiliary(m2sim *m2, m2aux *aux, double *X, double dV,
			 m2prim *P, double *U);
int nrhyd_eigenvalues(m2aux *aux, double n[4], double *evals);
double nrhyd_measure(m2aux *aux, int flag);


/* -----------------------------------------------------------------------------
 * SRHYD functions
 * ---------------------------------------------------------------------------*/
int srhyd_from_primitive(m2sim *m2, m2prim *P, double *B, double *X, double dV,
			 double *U, m2aux *aux);
int srhyd_from_conserved(m2sim *m2, double *U, double *B, double *X, double dV,
			 m2aux *aux, m2prim *P);
int srhyd_from_auxiliary(m2sim *m2, m2aux *aux, double *X, double dV,
			 m2prim *P, double *U);
int srhyd_eigenvalues(m2aux *aux, double n[4], double *evals);
double srhyd_measure(m2aux *aux, int flag);


/* -----------------------------------------------------------------------------
 * NRMHD functions
 * ---------------------------------------------------------------------------*/
int nrmhd_from_primitive(m2sim *m2, m2prim *P, double *B, double *X, double dV,
			 double *U, m2aux *aux);
int nrmhd_from_conserved(m2sim *m2, double *U, double *B, double *X, double dV,
			 m2aux *aux, m2prim *P);
int nrmhd_from_auxiliary(m2sim *m2, m2aux *aux, double *X, double dV,
			 m2prim *P, double *U);
int nrmhd_eigenvalues(m2aux *aux, double n[4], double *evals);
double nrmhd_measure(m2aux *aux, int flag);


#endif // M2_HYDRO_HEADER
