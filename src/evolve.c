#include <string.h>
#include "m2.h"


static void riemann_hll(m2vol *VL, m2vol *VR, double n[4], double *F)
{
  int q;
  double UL[5];
  double UR[5];
  double FL[5];
  double FR[5];
  double lamL[5];
  double lamR[5];
  m2aux AL;
  m2aux AR;

  m2sim_from_primitive(VL->m2, &VL->prim, NULL, NULL, 1.0, UL, &AL);
  m2sim_from_primitive(VR->m2, &VR->prim, NULL, NULL, 1.0, UR, &AR);
  m2aux_eigenvalues(&AL, n, lamL);
  m2aux_eigenvalues(&AR, n, lamR);
  m2aux_fluxes(&AL, n, FL);
  m2aux_fluxes(&AR, n, FR);

  double am = M2_MIN3(lamL[0], lamR[0], 0.0);
  double ap = M2_MAX3(lamL[4], lamR[4], 0.0);

  for (q=0; q<5; ++q) {
    F[q] = (ap*FL[q] - am*FR[q] + ap*am*(UR[q] - UL[q])) / (ap - am);
  }
}


void m2sim_calculate_flux(m2sim *m2)
{
  int i, j, k, q;
  int *L = m2->local_grid_size;
  double n[4] = { 0, 0, 0, 0 };
  m2vol *VL, *VR;
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {
	if (i != L[1] - 1) {
	  VL = m2->volumes + M2_IND(i+0,j,k);
	  VR = m2->volumes + M2_IND(i+1,j,k);
	  n[1] = 1.0;
	  riemann_hll(VL, VR, n, VL->flux1);
	  n[1] = 0.0;
	}
	else {
	  for (q=0; q<5; ++q) VL->flux1[q] = 0.0;	  
	}
	if (j != L[2] - 1) {
	  VL = m2->volumes + M2_IND(i,j+0,k);
	  VR = m2->volumes + M2_IND(i,j+1,k);
	  n[2] = 1.0;
	  riemann_hll(VL, VR, n, VL->flux2);
	  n[2] = 0.0;
	}
	else {
	  for (q=0; q<5; ++q) VL->flux2[q] = 0.0;	  
	}
	if (k != L[3] - 1) {
	  VL = m2->volumes + M2_IND(i,j,k+0);
	  VR = m2->volumes + M2_IND(i,j,k+1);
	  n[3] = 1.0;
	  riemann_hll(VL, VR, n, VL->flux3);
	  n[3] = 0.0;
	}
	else {
	  for (q=0; q<5; ++q) VL->flux3[q] = 0.0;	  
	}
      }
    }
  }
}


void m2sim_exchange_flux(m2sim *m2, double dt)
{
  int i, j, k, q;
  int *L = m2->local_grid_size;
  m2vol *VL, *VR;
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {
	if (i != L[1] - 1) {
	  VL = m2->volumes + M2_IND(i+0,j,k);
	  VR = m2->volumes + M2_IND(i+1,j,k);
	  for (q=0; q<5; ++q) {
	    VL->consA[q] -= dt * VL->area1 * VL->flux1[q];
	    VR->consA[q] += dt * VL->area1 * VL->flux1[q];
	  }
	}
	if (j != L[2] - 1) {
	  VL = m2->volumes + M2_IND(i,j+0,k);
	  VR = m2->volumes + M2_IND(i,j+1,k);
	  for (q=0; q<5; ++q) {
	    VL->consA[q] -= dt * VL->area2 * VL->flux2[q];
	    VR->consA[q] += dt * VL->area2 * VL->flux2[q];
	  }
	}
	if (k != L[3] - 1) {
	  VL = m2->volumes + M2_IND(i,j,k+0);
	  VR = m2->volumes + M2_IND(i,j,k+1);
	  for (q=0; q<5; ++q) {
	    VL->consA[q] -= dt * VL->area3 * VL->flux3[q];
	    VR->consA[q] += dt * VL->area3 * VL->flux3[q];
	  }
	}
      }
    }
  }
}


void m2sim_cache_conserved(m2sim *m2)
{
  int n;
  int *L = m2->local_grid_size;
  m2vol *V;
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    memcpy(V->consB, V->consA, 5 * sizeof(double));
  }
}


int m2sim_from_conserved_all(m2sim *m2)
{
  int n;
  int *L = m2->local_grid_size;
  m2vol *V;
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    m2sim_from_conserved(m2, V->consA, NULL, NULL, V->volume,
			 NULL, &V->prim);
  }
  return 0;
}


double m2sim_minimum_courant_time(m2sim *m2)
{
  int n;
  int *L = m2->local_grid_size;
  double dt, mindt=-1.0;
  m2vol *V;
  m2aux aux;
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    m2sim_from_primitive(m2, &V->prim, NULL, NULL, V->volume, NULL, &aux);
    dt = m2vol_minimum_dimension(V) / m2aux_maximum_wavespeed(&aux);
    if (dt < mindt || mindt < 0.0) {
      mindt = dt;
    }
  }
  return mindt;
}


void m2sim_enforce_boundary_condition(m2sim *m2)
{
  int *L = m2->local_grid_size;
  m2->volumes[     0] = m2->volumes[     1];
  m2->volumes[L[1]-2] = m2->volumes[L[1]-1];
}
