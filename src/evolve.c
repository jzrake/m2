#include <string.h>
#include "m2.h"


static void riemann_hll(m2vol *VL, m2vol *VR, double n[4], double *F)
{
  int q;
  double UL[8];
  double UR[8];
  double FL[8];
  double FR[8];
  double lamL[8];
  double lamR[8];
  m2aux AL;
  m2aux AR;

  m2sim_from_primitive(VL->m2, &VL->prim, NULL, NULL, 1.0, UL, &AL);
  m2sim_from_primitive(VR->m2, &VR->prim, NULL, NULL, 1.0, UR, &AR);
  m2aux_eigenvalues(&AL, n, lamL);
  m2aux_eigenvalues(&AR, n, lamR);
  m2aux_fluxes(&AL, n, FL);
  m2aux_fluxes(&AR, n, FR);

  UL[B11] = VL->prim.B1;
  UL[B22] = VL->prim.B2;
  UL[B33] = VL->prim.B3;
  UR[B11] = VR->prim.B1;
  UR[B22] = VR->prim.B2;
  UR[B33] = VR->prim.B3;

  double am = M2_MIN3(lamL[0], lamR[0], 0.0);
  double ap = M2_MAX3(lamL[7], lamR[7], 0.0);

  if (ap != ap || am != am) {
    MSG(FATAL, "got NAN wave speeds");
  }

  for (q=0; q<8; ++q) {
    F[q] = (ap*FL[q] - am*FR[q] + ap*am*(UR[q] - UL[q])) / (ap - am);
  }
}


void m2sim_calculate_flux(m2sim *m2)
{
  int i, j, k, q;
  int *L = m2->local_grid_size;
  double n1[4] = { 0, 1, 0, 0 };
  double n2[4] = { 0, 0, 1, 0 };
  double n3[4] = { 0, 0, 0, 1 };
  m2vol *V0, *V1, *V2, *V3;
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {

	V0 = M2_VOL(i+0, j+0, k+0);
	V1 = M2_VOL(i+1, j+0, k+0);
	V2 = M2_VOL(i+0, j+1, k+0);
	V3 = M2_VOL(i+0, j+0, k+1);

	if (V1) riemann_hll(V0, V1, n1, V0->flux1);
	else for (q=0; q<8; ++q) V0->flux1[q] = 0.0;

	if (V2) riemann_hll(V0, V2, n2, V0->flux2);
	else for (q=0; q<8; ++q) V0->flux2[q] = 0.0;

	if (V3) riemann_hll(V0, V3, n3, V0->flux3);
	else for (q=0; q<8; ++q) V0->flux3[q] = 0.0;
      }
    }
  }
}


void m2sim_calculate_emf(m2sim *m2)
{
  int i, j, k;
  int *G = m2->domain_resolution;
  int *L = m2->local_grid_size;
  m2vol *V0, *V1, *V2;
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {

	/* ---------------------- */
	/* EMF in the x-direction */
	/* ---------------------- */
	V0 = M2_VOL(i, j+0, k+0);
	V1 = M2_VOL(i, j+1, k+0); V1 = V1 ? V1 : V0;
	V2 = M2_VOL(i, j+0, k+1); V2 = V2 ? V2 : V0;
	if (G[2] == 1 && G[3] == 1) V0->emf1 = 0.0;
	else if (G[2] == 1) V0->emf1 = 0.50 * (+V0->flux3[B22]*V0->line1
					       +V1->flux3[B22]*V0->line1);
	else if (G[3] == 1) V0->emf1 = 0.50 * (-V0->flux2[B33]*V0->line1
					       -V2->flux2[B33]*V0->line1);
	else V0->emf1 = 0.25 * (-V0->flux2[B33]*V0->line1
				-V2->flux2[B33]*V0->line1
				+V0->flux3[B22]*V0->line1
				+V1->flux3[B22]*V0->line1);


	/* ---------------------- */
	/* EMF in the y-direction */
	/* ---------------------- */
	V0 = M2_VOL(i+0, j, k+0);
	V1 = M2_VOL(i+0, j, k+1); V1 = V1 ? V1 : V0;
	V2 = M2_VOL(i+1, j, k+0); V2 = V2 ? V2 : V0;
	if (G[3] == 1 && G[1] == 1) V0->emf2 = 0.0;
	else if (G[3] == 1) V0->emf2 = 0.50 * (+V0->flux1[B33]*V0->line2
					       +V1->flux1[B33]*V0->line2);
	else if (G[1] == 1) V0->emf2 = 0.50 * (-V0->flux3[B11]*V0->line2
					       -V2->flux3[B11]*V0->line2);
	else V0->emf2 = 0.25 * (-V0->flux3[B11]*V0->line2
				-V2->flux3[B11]*V0->line2
				+V0->flux1[B33]*V0->line2
				+V1->flux1[B33]*V0->line2);


	/* ---------------------- */
	/* EMF in the z-direction */
	/* ---------------------- */
	V0 = M2_VOL(i+0, j+0, k);
	V1 = M2_VOL(i+1, j+0, k); V1 = V1 ? V1 : V0;
	V2 = M2_VOL(i+0, j+1, k); V2 = V2 ? V2 : V0;
	if (G[1] == 1 && G[2] == 1) V0->emf3 = 0.0;
	else if (G[1] == 1) V0->emf3 = 0.50 * (+V0->flux2[B11]*V0->line3
					       +V1->flux2[B11]*V0->line3);
	else if (G[2] == 1) V0->emf3 = 0.50 * (-V0->flux1[B22]*V0->line3
					       -V2->flux1[B22]*V0->line3);
	else V0->emf3 = 0.25 * (-V0->flux1[B22]*V0->line3
				-V2->flux1[B22]*V0->line3
				+V0->flux2[B11]*V0->line3
				+V1->flux2[B11]*V0->line3);

	/* if (V0->global_index[2] == 0) { */
	/*   printf("on the North pole: theta-EMF=%f phi-EMF=%f\n", */
	/* 	 V0->emf2, V0->emf3); */
	/* } */
      }
    }
  }
}


void m2sim_exchange_flux(m2sim *m2, double dt)
{
  int i, j, k, q;
  int *G = m2->domain_resolution;
  int *L = m2->local_grid_size;
  m2vol *V0, *V1, *V2, *V3;

  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {

	if (m2->physics & M2_MAGNETIZED) {

	  /* --------------- */
	  /* x-directed face */
	  /* --------------- */
	  V0 = M2_VOL(i, j+0, k+0);
	  V1 = M2_VOL(i, j-1, k+0); V1 = V1 ? V1 : V0;
	  V2 = M2_VOL(i, j+0, k-1); V2 = V2 ? V2 : V0;
	  if (G[2] == 1 && G[3] == 1) { /* symmetric in y and z */

	  }
	  else if (G[2] == 1) { /* symmetric in y */
	    V0->Bflux1A += dt * V0->emf2;
	    V0->Bflux1A -= dt * V2->emf2;
	  }
	  else if (G[3] == 1) { /* symmetric in z */
	    V0->Bflux1A += dt * V1->emf3;
	    V0->Bflux1A -= dt * V0->emf3;
	  }
	  else {
	    V0->Bflux1A += dt * V0->emf2;
	    V0->Bflux1A -= dt * V2->emf2;
	    V0->Bflux1A += dt * V1->emf3;
	    V0->Bflux1A -= dt * V0->emf3;
	  }

	  /* --------------- */
	  /* y-directed face */
	  /* --------------- */
	  V0 = M2_VOL(i+0, j, k+0);
	  V1 = M2_VOL(i+0, j, k-1); V1 = V1 ? V1 : V0;
	  V2 = M2_VOL(i-1, j, k+0); V2 = V2 ? V2 : V0;
	  if (G[3] == 1 && G[1] == 1) { /* symmetric in z and x */

	  }
	  else if (G[3] == 1) { /* symmetric in z */
	    V0->Bflux2A += dt * V0->emf3;
	    V0->Bflux2A -= dt * V2->emf3;
	  }
	  else if (G[1] == 1) { /* symmetric in x */
	    V0->Bflux2A += dt * V1->emf1;
	    V0->Bflux2A -= dt * V0->emf1;
	  }
	  else {
	    V0->Bflux2A += dt * V0->emf3;
	    V0->Bflux2A -= dt * V2->emf3;
	    V0->Bflux2A += dt * V1->emf1;
	    V0->Bflux2A -= dt * V0->emf1;
	  }

	  /* --------------- */
	  /* z-directed face */
	  /* --------------- */
	  V0 = M2_VOL(i+0, j+0, k);
	  V1 = M2_VOL(i-1, j+0, k); V1 = V1 ? V1 : V0;
	  V2 = M2_VOL(i+0, j-1, k); V2 = V2 ? V2 : V0;
	  if (G[1] == 1 && G[2] == 1) { /* symmetric in x and y */

	  }
	  else if (G[1] == 1) { /* symmetric in x */
	    V0->Bflux3A += dt * V0->emf1;
	    V0->Bflux3A -= dt * V2->emf1;
	  }
	  else if (G[2] == 1) { /* symmetric in y */
	    V0->Bflux3A += dt * V1->emf2;
	    V0->Bflux3A -= dt * V0->emf2;
	  }
	  else {
	    V0->Bflux3A += dt * V0->emf1;
	    V0->Bflux3A -= dt * V2->emf1;
	    V0->Bflux3A += dt * V1->emf2;
	    V0->Bflux3A -= dt * V0->emf2;
	  }
	}



	V0 = M2_VOL(i+0, j+0, k+0);
	V1 = M2_VOL(i+1, j+0, k+0);
	V2 = M2_VOL(i+0, j+1, k+0);
	V3 = M2_VOL(i+0, j+0, k+1);

	for (q=0; q<5; ++q) {
	  /* x-face */
	  if (V0) V0->consA[q] -= dt * V0->area1 * V0->flux1[q];
	  if (V1) V1->consA[q] += dt * V0->area1 * V0->flux1[q];

	  /* y-face */
	  if (V0) V0->consA[q] -= dt * V0->area2 * V0->flux2[q];
	  if (V2) V2->consA[q] += dt * V0->area2 * V0->flux2[q];

	  /* z-face */
	  if (V0) V0->consA[q] -= dt * V0->area3 * V0->flux3[q];
	  if (V3) V3->consA[q] += dt * V0->area3 * V0->flux3[q];
	}
      }
    }
  }
}


void m2sim_add_source_terms(m2sim *m2, double dt)
{
  int n;
  int *L = m2->local_grid_size;
  int *G = m2->domain_resolution;
  int *I;
  m2vol *V;
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    V->x0[0] = 0.0;
    V->x1[0] = dt;
    I = V->global_index;
    if (I[1] < 0 ||
	I[2] < 0 ||
	I[3] < 0 ||
	I[1] >= G[1] ||
	I[2] >= G[2] ||
	I[3] >= G[3]) {
    }
    else {
      m2aux_add_geometrical_source_terms(&V->aux, V->x0, V->x1, V->consA);
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
    V->Bflux1B = V->Bflux1A;
    V->Bflux2B = V->Bflux2A;
    V->Bflux3B = V->Bflux3A;
  }
}


void m2sim_average_runge_kutta(m2sim *m2, double b)
{
  int n, q;
  int *L = m2->local_grid_size;
  m2vol *V;
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    for (q=0; q<5; ++q) {
      V->consA[q] = b * V->consA[q] + (1.0 - b) * V->consB[q];
    }
    V->Bflux1A = b * V->Bflux1A + (1.0 - b) * V->Bflux1B;
    V->Bflux2A = b * V->Bflux2A + (1.0 - b) * V->Bflux2B;
    V->Bflux3A = b * V->Bflux3A + (1.0 - b) * V->Bflux3B;
  }
}


void m2sim_magnetic_flux_to_cell_center(m2sim *m2)
/*
 * Convert the face-centered magnetic flux to a pointwise magnetic field at the
 * cell center. Although the left-most cell along each axis is given a
 * cell-centered field value (using only its +1/2 face), it should never be used
 * since the field value at its (-1/2) face is unknown.
 */
{
  int i, j, k;
  int *L = m2->local_grid_size;
  m2vol *VC, *V1, *V2, *V3;
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {
	VC = m2->volumes + M2_IND(i, j, k);
	V1 = (i == 0) ? VC : m2->volumes + M2_IND(i-1, j, k);
	V2 = (j == 0) ? VC : m2->volumes + M2_IND(i, j-1, k);
	V3 = (k == 0) ? VC : m2->volumes + M2_IND(i, j, k-1);
	VC->prim.B1 = 0.5 * (VC->Bflux1A/VC->area1 + V1->Bflux1A/V1->area1);
	VC->prim.B2 = 0.5 * (VC->Bflux2A/VC->area2 + V2->Bflux2A/V2->area2);
	VC->prim.B3 = 0.5 * (VC->Bflux3A/VC->area3 + V3->Bflux3A/V3->area3);
      }
    }
  }
}


int m2sim_from_conserved_all(m2sim *m2)
{
  int n, error;
  int *L = m2->local_grid_size;
  m2vol *V;
  double B[4];
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    /* assume fields have been averaged from faces to center (prim) */
    B[1] = V->prim.B1;
    B[2] = V->prim.B2;
    B[3] = V->prim.B3;
    error = m2sim_from_conserved(m2, V->consA, B, NULL, V->volume,
				 &V->aux, &V->prim);
    if (error) {
      MSGF(FATAL, "at global index [%d %d %d]",
	   V->global_index[1],
	   V->global_index[2],
	   V->global_index[3]);
    }
  }
  return 0;
}


double m2sim_minimum_courant_time(m2sim *m2)
{
  int n;
  int *L = m2->local_grid_size;
  double dt, mindt = -1.0;
  m2vol *V;
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    dt = m2vol_minimum_dimension(V) / m2aux_maximum_wavespeed(&V->aux);
    if (dt < mindt || mindt < 0.0) {
      mindt = dt;
    }
  }
  return mindt;
}


void initial_data(m2vol *V);
void m2sim_enforce_boundary_condition(m2sim *m2)
{
  int n;
  int *L = m2->local_grid_size;
  int *G = m2->domain_resolution;
  int *I;
  m2vol *V;
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    I = V->global_index;
    if (I[1] < 0 ||
  	I[2] < 0 ||
  	I[3] < 0 ||
  	I[1] >= G[1] ||
  	I[2] >= G[2] ||
  	I[3] >= G[3]) {
      initial_data(V);
      m2sim_from_primitive(m2,
  			   &V->prim, NULL, NULL,
  			   V ->volume,
  			   V ->consA,
  			   &V->aux);
    }
  }
}
