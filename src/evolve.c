#include <string.h>
#include <math.h>
#if _OPENMP
#include <omp.h>
#endif
#include "m2.h"

#define ENABLE_PLM 1

static void riemann_hll(m2vol *VL, m2vol *VR, int axis, double *F)
{
  int q;
  double UL[8];
  double UR[8];
  double FL[8];
  double FR[8];
  double lamL[8];
  double lamR[8];
  double yl[8], yr[8]; /* cell-centered left/right interpolation variables */
  double yL[8], yR[8]; /* face-centered left/right interpolation variables */
  m2prim PL;
  m2prim PR;
  m2aux AL;
  m2aux AR;

  double xL = m2vol_coordinate_centroid(VL, axis);
  double xR = m2vol_coordinate_centroid(VR, axis);
  double n[4] = { 0.0, 0.0, 0.0, 0.0 };

  if (VL->zone_type != M2_ZONE_TYPE_FULL ||
      VR->zone_type != M2_ZONE_TYPE_FULL) {
    return;
  }

  n[axis] = 1.0;
  m2vol_to_interpolated(VL, yl, 1);
  m2vol_to_interpolated(VR, yr, 1);

  switch (axis) {
  case 1:
    for (q=0; q<8; ++q) yL[q] = yl[q] + VL->grad1[q] * (VL->x1[1] - xL);
    for (q=0; q<8; ++q) yR[q] = yr[q] + VR->grad1[q] * (VR->x0[1] - xR);
    break;
  case 2:
    for (q=0; q<8; ++q) yL[q] = yl[q] + VL->grad2[q] * (VL->x1[2] - xL);
    for (q=0; q<8; ++q) yR[q] = yr[q] + VR->grad2[q] * (VR->x0[2] - xR);
    break;
  case 3:
    for (q=0; q<8; ++q) yL[q] = yl[q] + VL->grad3[q] * (VL->x1[3] - xL);
    for (q=0; q<8; ++q) yR[q] = yr[q] + VR->grad3[q] * (VR->x0[3] - xR);
    break;
  }

  m2sim_from_interpolated(VL->m2, yL, &PL);
  m2sim_from_interpolated(VR->m2, yR, &PR);

  UL[B11] = PL.B1;
  UL[B22] = PL.B2;
  UL[B33] = PL.B3;
  UR[B11] = PR.B1;
  UR[B22] = PR.B2;
  UR[B33] = PR.B3;

  switch (axis) {
  case 1: UL[B11] = UR[B11] = PL.B1 = PR.B1 = VL->Bflux1A / VL->area1; break;
  case 2: UL[B22] = UR[B22] = PL.B2 = PR.B2 = VL->Bflux2A / VL->area2; break;
  case 3: UL[B33] = UR[B33] = PL.B3 = PR.B3 = VL->Bflux3A / VL->area3; break;
  }

  m2sim_from_primitive(VL->m2, &PL, NULL, NULL, 1.0, UL, &AL);
  m2sim_from_primitive(VR->m2, &PR, NULL, NULL, 1.0, UR, &AR);
  m2aux_eigenvalues(&AL, n, lamL);
  m2aux_eigenvalues(&AR, n, lamR);
  m2aux_fluxes(&AL, n, FL);
  m2aux_fluxes(&AR, n, FR);

  double am = M2_MIN3(lamL[0], lamR[0], 0.0);
  double ap = M2_MAX3(lamL[7], lamR[7], 0.0);

  if (ap != ap || am != am) {
    MSG(FATAL, "got NAN wave speeds");
  }

  for (q=0; q<8; ++q) {
    F[q] = (ap*FL[q] - am*FR[q] + ap*am*(UR[q] - UL[q])) / (ap - am);
  }
}

static double plm_gradient(double *xs, double *fs, double plm)
{
#define minval3(a,b,c) ((a<b) ? ((a<c) ? a : c) : ((b<c) ? b : c))
#define sign(x) ((x > 0.0) - (x < 0.0))

  double xL = xs[0];
  double x0 = xs[1];
  double xR = xs[2];

  double yL = fs[0];
  double y0 = fs[1];
  double yR = fs[2];

  double a = (yL - y0) / (xL - x0) * plm;
  double b = (yL - yR) / (xL - xR);
  double c = (y0 - yR) / (x0 - xR) * plm;

  double sa = sign(a), sb = sign(b), sc = sign(c);
  double fa = fabs(a), fb = fabs(b), fc = fabs(c);
  double g = 0.25 * fabs(sa + sb) * (sa + sc) * minval3(fa, fb, fc);

  if (g != g) {
    MSGF(FATAL, "got NAN gradient: [%f %f %f] [%f %f %f]",
	 xs[0], xs[1], xs[2], fs[0], fs[1], fs[2]);
  }

  return g;
#undef minval3
#undef sign
}


void m2sim_calculate_flux(m2sim *m2)
{
  int i, j, k, q;
  int *L = m2->local_grid_size;
  int *G = m2->domain_resolution;
  m2vol *V0, *V1, *V2, *V3;
#ifdef _OPENMP
#pragma omp parallel for private(V0,V1,V2,V3,j,k,q)
#endif
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {

	V0 = M2_VOL(i+0, j+0, k+0);
	V1 = M2_VOL(i+1, j+0, k+0);
	V2 = M2_VOL(i+0, j+1, k+0);
	V3 = M2_VOL(i+0, j+0, k+1);


	/* -----------------------------------------
	 * axis 1
	 * ----------------------------------------- */
	if (V0->global_index[1] == -1 ||
	    V0->global_index[1] == G[1] - 1) {
	  if (m2->boundary_conditions_flux1) m2->boundary_conditions_flux1(V0);
	  else for (q=0; q<8; ++q) V0->flux1[q] = 0.0;
	}
	else if (V1 == NULL) for (q=0; q<8; ++q) V0->flux1[q] = 0.0;
	else riemann_hll(V0, V1, 1, V0->flux1);


	/* -----------------------------------------
	 * axis 2
	 * ----------------------------------------- */
	if (V0->global_index[2] == -1 ||
	    V0->global_index[2] == G[2] - 1) {
	  if (m2->boundary_conditions_flux2) m2->boundary_conditions_flux2(V0);
	  else for (q=0; q<8; ++q) V0->flux2[q] = 0.0;
	}
	else if (V2 == NULL) for (q=0; q<8; ++q) V0->flux2[q] = 0.0;
	else riemann_hll(V0, V2, 2, V0->flux2);


	/* -----------------------------------------
	 * axis 3
	 * ----------------------------------------- */
	if (V0->global_index[3] == -1 ||
	    V0->global_index[3] == G[3] - 1) {
	  if (m2->boundary_conditions_flux3) m2->boundary_conditions_flux3(V0);
	  else for (q=0; q<8; ++q) V0->flux3[q] = 0.0;
	}
	else if (V3 == NULL) for (q=0; q<8; ++q) V0->flux3[q] = 0.0;
	else riemann_hll(V0, V3, 3, V0->flux3);
      }
    }
  }
}


void m2sim_calculate_emf(m2sim *m2)
{
  int i, j, k;
  int *G = m2->domain_resolution;
  int *L = m2->local_grid_size;
  int *I;
  m2vol *V0, *V1, *V2;

#define EMF(a1,a2,a3,B1,B2,B3,f1,f2,f3,line,emf)			\
  do {									\
    if (G[a2] == 1 && G[a3] == 1) {					\
      V0->emf = 0.0;							\
    }									\
    else if (G[a2] == 1) {						\
      if (I[a3] == -1 || I[a3] == G[a3] - 1) {				\
	if (m2->boundary_conditions_##emf){m2->boundary_conditions_##emf(V0);} \
	else V0->emf = 0.0;						\
      }									\
      else if (V1 == NULL) V0->emf = 0.0; /* processor boundary */	\
      else {								\
	V0->emf = 0.50 * (+V0->f3[B2]*V0->line				\
			  +V1->f3[B2]*V0->line);			\
      }									\
    }									\
    else if (G[a3] == 1) {						\
      if (I[a2] == -1 || I[a2] == G[a2] - 1) {				\
	if (m2->boundary_conditions_##emf){m2->boundary_conditions_##emf(V0);} \
	else V0->emf = 0.0;						\
      }									\
      else if (V2 == NULL) V0->emf = 0.0; /* processor boundary */	\
      else {								\
	V0->emf = 0.50 * (-V0->f2[B3]*V0->line				\
			  -V2->f2[B3]*V0->line);			\
      }									\
    }									\
    else {								\
      if (I[a2] == -1 || I[a2] == G[a2] - 1 ||				\
	  I[a3] == -1 || I[a3] == G[a3] - 1) {				\
	if (m2->boundary_conditions_##emf){m2->boundary_conditions_##emf(V0);} \
	else V0->emf = 0.0;						\
      }									\
      else if (V1 == NULL || V2 == NULL) { /* processor boundary */	\
	V0->emf = 0.0;							\
      }									\
      else {								\
	V0->emf = 0.25 * (-V0->f2[B3]*V0->line				\
			  -V2->f2[B3]*V0->line				\
			  +V0->f3[B2]*V0->line				\
			  +V1->f3[B2]*V0->line);			\
      }									\
    }									\
  }									\
  while (0)								\

#ifdef _OPENMP
#pragma omp parallel for private(V0,V1,V2,j,k) default(shared)
#endif
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {
	I = M2_VOL(i, j, k)->global_index;

	V0 = M2_VOL(i, j+0, k+0);
	V1 = M2_VOL(i, j+1, k+0);
	V2 = M2_VOL(i, j+0, k+1);
	EMF(1,2,3,B11,B22,B33,flux1,flux2,flux3,line1,emf1);

	V0 = M2_VOL(i+0, j, k+0);
	V1 = M2_VOL(i+0, j, k+1);
	V2 = M2_VOL(i+1, j, k+0);
	EMF(2,3,1,B22,B33,B11,flux2,flux3,flux1,line2,emf2);

	V0 = M2_VOL(i+0, j+0, k);
	V1 = M2_VOL(i+1, j+0, k);
	V2 = M2_VOL(i+0, j+1, k);
	EMF(3,1,2,B33,B11,B22,flux3,flux1,flux2,line3,emf3);
      }
    }
  }
#undef EMF
}


void m2sim_calculate_gradient(m2sim *m2)
{
  int i, j, k, q, n;
  int *L = m2->local_grid_size;
  m2vol *V[3];
  double x[3];
  double ys[24];
#ifdef _OPENMP
#pragma omp parallel for private(V,x,ys,j,k,q,n) default(shared)
#endif
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {


	/* ----------------------- */
	/* gradient in x-direction */
	/* ----------------------- */
	V[0] = M2_VOL(i-1, j, k);
	V[1] = M2_VOL(i+0, j, k);
	V[2] = M2_VOL(i+1, j, k);
	if (V[0] && V[0]->zone_type != M2_ZONE_TYPE_SHELL && 
	    V[2] && V[2]->zone_type != M2_ZONE_TYPE_SHELL && ENABLE_PLM) {
	  for (n=0; n<3; ++n) {
	    x[n] = m2vol_coordinate_centroid(V[n], 1);
	  }
	  m2vol_to_interpolated(V[0], &ys[0], 3);
	  m2vol_to_interpolated(V[1], &ys[1], 3);
	  m2vol_to_interpolated(V[2], &ys[2], 3);
	  for (q=0; q<8; ++q) {
	    V[1]->grad1[q] = plm_gradient(x, &ys[3*q], m2->plm_parameter);
	  }
	}
	else {
	  for (q=0; q<8; ++q) {
	    V[1]->grad1[q] = 0.0;
	  }
	}


	/* ----------------------- */
	/* gradient in y-direction */
	/* ----------------------- */
	V[0] = M2_VOL(i, j-1, k);
	V[1] = M2_VOL(i, j+0, k);
	V[2] = M2_VOL(i, j+1, k);
	if (V[0] && V[0]->zone_type != M2_ZONE_TYPE_SHELL && 
	    V[2] && V[2]->zone_type != M2_ZONE_TYPE_SHELL && ENABLE_PLM) {
	  for (n=0; n<3; ++n) {
	    x[n] = m2vol_coordinate_centroid(V[n], 2);
	  }
	  m2vol_to_interpolated(V[0], &ys[0], 3);
	  m2vol_to_interpolated(V[1], &ys[1], 3);
	  m2vol_to_interpolated(V[2], &ys[2], 3);
	  for (q=0; q<8; ++q) {
	    V[1]->grad2[q] = plm_gradient(x, &ys[3*q], m2->plm_parameter);
	  }
	}
	else {
	  for (q=0; q<8; ++q) {
	    V[1]->grad2[q] = 0.0;
	  }
	}


	/* ----------------------- */
	/* gradient in z-direction */
	/* ----------------------- */
	V[0] = M2_VOL(i, j, k-1);
	V[1] = M2_VOL(i, j, k+0);
	V[2] = M2_VOL(i, j, k+1);
	if (V[0] && V[0]->zone_type != M2_ZONE_TYPE_SHELL && 
	    V[2] && V[2]->zone_type != M2_ZONE_TYPE_SHELL && ENABLE_PLM) {
	  for (n=0; n<3; ++n) {
	    x[n] = m2vol_coordinate_centroid(V[n], 3);
	  }
	  m2vol_to_interpolated(V[0], &ys[0], 3);
	  m2vol_to_interpolated(V[1], &ys[1], 3);
	  m2vol_to_interpolated(V[2], &ys[2], 3);
	  for (q=0; q<8; ++q) {
	    V[1]->grad3[q] = plm_gradient(x, &ys[3*q], m2->plm_parameter);
	  }
	}
	else {
	  for (q=0; q<8; ++q) {
	    V[1]->grad3[q] = 0.0;
	  }
	}
      }
    }
  }

  if (m2->boundary_conditions_gradient) {
    m2->boundary_conditions_gradient(m2);
  }
  else {
    /* guard zones will just have zero gradient */
  }
}


void m2sim_exchange_flux(m2sim *m2, double dt)
{
  int i, j, k, q;
  int *G = m2->domain_resolution;
  int *L = m2->local_grid_size;
  m2vol *V0, *V1, *V2, *V3;
#ifdef _OPENMP
#pragma omp parallel for private(V0,V1,V2,V3,j,k,q) default(shared)
#endif
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {

	if (m2->physics & M2_MAGNETIZED) {
	  /* --------------- */
	  /* x-directed face */
	  /* --------------- */
	  V0 = M2_VOL(i, j+0, k+0);
	  V1 = M2_VOL(i, j-1, k+0);
	  V2 = M2_VOL(i, j+0, k-1);
	  if (G[2] == 1 && G[3] == 1) { /* symmetric in y and z */

	  }
	  else if (G[2] == 1) { /* symmetric in y */
	    if (V0) V0->Bflux1A += dt * V0->emf2;
	    if (V0) V0->Bflux1A -= dt * V2->emf2;
	  }
	  else if (G[3] == 1) { /* symmetric in z */
	    if (V1) V0->Bflux1A += dt * V1->emf3;
	    if (V0) V0->Bflux1A -= dt * V0->emf3;
	  }
	  else {
	    if (V0) V0->Bflux1A += dt * V0->emf2;
	    if (V2) V0->Bflux1A -= dt * V2->emf2;
	    if (V1) V0->Bflux1A += dt * V1->emf3;
	    if (V0) V0->Bflux1A -= dt * V0->emf3;
	  }

	  /* --------------- */
	  /* y-directed face */
	  /* --------------- */
	  V0 = M2_VOL(i+0, j, k+0);
	  V1 = M2_VOL(i+0, j, k-1);
	  V2 = M2_VOL(i-1, j, k+0);
	  if (G[3] == 1 && G[1] == 1) { /* symmetric in z and x */

	  }
	  else if (G[3] == 1) { /* symmetric in z */
	    if (V0) V0->Bflux2A += dt * V0->emf3;
	    if (V2) V0->Bflux2A -= dt * V2->emf3;
	  }
	  else if (G[1] == 1) { /* symmetric in x */
	    if (V1) V0->Bflux2A += dt * V1->emf1;
	    if (V0) V0->Bflux2A -= dt * V0->emf1;
	  }
	  else {
	    if (V0) V0->Bflux2A += dt * V0->emf3;
	    if (V2) V0->Bflux2A -= dt * V2->emf3;
	    if (V1) V0->Bflux2A += dt * V1->emf1;
	    if (V0) V0->Bflux2A -= dt * V0->emf1;
	  }

	  /* --------------- */
	  /* z-directed face */
	  /* --------------- */
	  V0 = M2_VOL(i+0, j+0, k);
	  V1 = M2_VOL(i-1, j+0, k);
	  V2 = M2_VOL(i+0, j-1, k);
	  if (G[1] == 1 && G[2] == 1) { /* symmetric in x and y */

	  }
	  else if (G[1] == 1) { /* symmetric in x */
	    if (V0) V0->Bflux3A += dt * V0->emf1;
	    if (V2) V0->Bflux3A -= dt * V2->emf1;
	  }
	  else if (G[2] == 1) { /* symmetric in y */
	    if (V1) V0->Bflux3A += dt * V1->emf2;
	    if (V0) V0->Bflux3A -= dt * V0->emf2;
	  }
	  else {
	    if (V0) V0->Bflux3A += dt * V0->emf1;
	    if (V2) V0->Bflux3A -= dt * V2->emf1;
	    if (V1) V0->Bflux3A += dt * V1->emf2;
	    if (V0) V0->Bflux3A -= dt * V0->emf2;
	  }
	}



	V0 = M2_VOL(i+0, j+0, k+0);
	V1 = M2_VOL(i-1, j+0, k+0);
	V2 = M2_VOL(i+0, j-1, k+0);
	V3 = M2_VOL(i+0, j+0, k-1);

	for (q=0; q<5; ++q) {
	  V0->consA[q] -= dt * V0->area1 * V0->flux1[q];
	  V0->consA[q] -= dt * V0->area2 * V0->flux2[q];
	  V0->consA[q] -= dt * V0->area3 * V0->flux3[q];
	  if (V1) V0->consA[q] += dt * V1->area1 * V1->flux1[q];
	  if (V2) V0->consA[q] += dt * V2->area2 * V2->flux2[q];
	  if (V3) V0->consA[q] += dt * V3->area3 * V3->flux3[q];
	}
      }
    }
  }
}


void m2sim_add_source_terms(m2sim *m2, double dt)
{
  int n;
  int *L = m2->local_grid_size;
  m2vol *V;
#ifdef _OPENMP
#pragma omp parallel for private(V) default(shared)
#endif
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    if (V->zone_type != M2_ZONE_TYPE_FULL) {
      continue;
    }
    V->x0[0] = 0.0;
    V->x1[0] = dt;
    m2aux_add_geometrical_source_terms(&V->aux, V->x0, V->x1, V->consA);
    if (m2->add_physical_source_terms) {
      m2->add_physical_source_terms(V);
    }
  }
}


void m2sim_cache_conserved(m2sim *m2)
{
  int n;
  int *L = m2->local_grid_size;
  m2vol *V;
#ifdef _OPENMP
#pragma omp parallel for private(V) default(shared)
#endif
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    if (V->zone_type != M2_ZONE_TYPE_FULL) {
      continue;
    }
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
#ifdef _OPENMP
#pragma omp parallel for private(V,q) default(shared)
#endif
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    if (V->zone_type != M2_ZONE_TYPE_FULL) {
      continue;
    }
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
  int *G = m2->domain_resolution;
  int *L = m2->local_grid_size;
  m2vol *VC, *V1, *V2, *V3;
  double A1L, A2L, A3L;
  double A1R, A2R, A3R;
#ifdef _OPENMP
#pragma omp parallel for private(VC,V1,V2,V3,j,k) default(shared)
#endif
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {
	VC = M2_VOL(i, j, k);
	V1 = (G[1] > 1) ? M2_VOL(i-1, j, k) : VC;
	V2 = (G[2] > 1) ? M2_VOL(i, j-1, k) : VC;
	V3 = (G[3] > 1) ? M2_VOL(i, j, k-1) : VC;

	A1R = (VC && VC->area1 > 1e-12) ? VC->area1 : 1.0;
	A2R = (VC && VC->area2 > 1e-12) ? VC->area2 : 1.0;
	A3R = (VC && VC->area3 > 1e-12) ? VC->area3 : 1.0;
	A1L = (V1 && V1->area1 > 1e-12) ? V1->area1 : 1.0;
	A2L = (V2 && V2->area2 > 1e-12) ? V2->area2 : 1.0;
	A3L = (V3 && V3->area3 > 1e-12) ? V3->area3 : 1.0;

	VC->prim.B1 = V1 ? 0.5 * (VC->Bflux1A/A1R + V1->Bflux1A/A1L) : 0.0;
	VC->prim.B2 = V2 ? 0.5 * (VC->Bflux2A/A2R + V2->Bflux2A/A2L) : 0.0;
	VC->prim.B3 = V3 ? 0.5 * (VC->Bflux3A/A3R + V3->Bflux3A/A3L) : 0.0;
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
#ifdef _OPENMP
#pragma omp parallel for private(V,B,error) default(shared)
#endif
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    if (V->zone_type != M2_ZONE_TYPE_FULL) {
      continue;
    }
    /* assume fields have been averaged from faces to center (prim) */
    B[0] = 0.0;
    B[1] = V->prim.B1;
    B[2] = V->prim.B2;
    B[3] = V->prim.B3;
    error = m2sim_from_conserved(m2, V->consA, B, NULL, V->volume,
				 &V->aux, &V->prim);
    if (error) {
      //m2->initial_data(V);
      //m2sim_from_primitive(V->m2, &V->prim, NULL, NULL, V->volume,
      //		   V->consA, &V->aux);
      m2_print_state(&V->prim, &V->aux, V->consA);
      MSGF(FATAL, "at global index [%d %d %d]",
    	   V->global_index[1],
    	   V->global_index[2],
    	   V->global_index[3]);
    }
  }
  return 0;
}


int m2sim_from_primitive_all(m2sim *m2)
{
  int *L = m2->local_grid_size;
  int n;
  m2vol *V;
#ifdef _OPENMP
#pragma omp parallel for private(V) default(shared)
#endif
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    if (V->zone_type != M2_ZONE_TYPE_FULL) {
      continue;
    }
    m2sim_from_primitive(V->m2, &V->prim, NULL, NULL, V->volume,
			 V->consA, &V->aux);
  }
  return 0;
}


double m2sim_minimum_courant_time(m2sim *m2)
{
  int n;
  int *L = m2->local_grid_size;
  double dt, mindt = -1.0;
  m2vol *V;
  double globalmindt = mindt;
#ifdef _OPENMP
#pragma omp parallel private(V,dt) firstprivate(mindt) default(shared)
#endif
  {
#ifdef _OPENMP
#pragma omp for
#endif
    for (n=0; n<L[0]; ++n) {
      V = m2->volumes + n;
      if (V->zone_type != M2_ZONE_TYPE_FULL) {
	continue;
      }
      dt = m2vol_minimum_dimension(V) / m2aux_maximum_wavespeed(&V->aux);
      if (dt < mindt || mindt < 0.0) {
	mindt = dt;
      }
    }
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      if (mindt < globalmindt || globalmindt < 0.0) {
	globalmindt = mindt;
      }
    }
  }
  return globalmindt;
}


void m2sim_enforce_boundary_condition(m2sim *m2)
{
  if (m2->boundary_conditions == NULL) {
    MSG(FATAL, "no boundary conditions supplied");
  }
  else {
    m2->boundary_conditions(m2);
  }
}


void m2sim_runge_kutta_substep(m2sim *m2, double dt, double rkparam)
{
  m2sim_calculate_gradient(m2);
  m2sim_calculate_flux(m2);
  m2sim_calculate_emf(m2);
  m2sim_exchange_flux(m2, dt);
  m2sim_add_source_terms(m2, dt);
  m2sim_average_runge_kutta(m2, rkparam);
  m2sim_enforce_boundary_condition(m2);
  m2sim_magnetic_flux_to_cell_center(m2);
  m2sim_from_conserved_all(m2);
}


void m2sim_drive(m2sim *m2)
{
  double dt;
  int rk_order = m2->rk_order;
  double kzps; /* kilozones per second */
  char checkpoint_name[1024];

#ifdef _OPENMP
  double start_wtime, stop_wtime;
#else
  clock_t start_cycle, stop_cycle;
#endif

  m2sim_run_analysis(m2);
  m2sim_cache_conserved(m2);
  if (m2->status.checkpoint_cadence * m2->status.checkpoint_number <
      m2->status.time_simulation) {
    m2->status.checkpoint_number += 1;
    snprintf(checkpoint_name, sizeof(checkpoint_name), "chkpt.%04d.m2",
	     m2->status.checkpoint_number);
    m2sim_save_checkpoint(m2, checkpoint_name);
  }
  dt = m2->cfl_parameter * m2sim_minimum_courant_time(m2);

#ifdef _OPENMP
  start_wtime = omp_get_wtime();
#else
  start_cycle = clock();
#endif

  switch (rk_order) {
  case 1:
    m2sim_runge_kutta_substep(m2, dt, 1.0);
    break;
  case 2:
    m2sim_runge_kutta_substep(m2, dt, 1.0);
    m2sim_runge_kutta_substep(m2, dt, 0.5);
    break;
  case 3:
    m2sim_runge_kutta_substep(m2, dt, 1.0);
    m2sim_runge_kutta_substep(m2, dt, 1.0/4.0);
    m2sim_runge_kutta_substep(m2, dt, 2.0/3.0);
    break;
  }

#ifdef _OPENMP
  stop_wtime = omp_get_wtime();
  kzps = 1e-3 * m2->local_grid_size[0] / (stop_wtime - start_wtime);
#else
  stop_cycle = clock();
  kzps = 1e-3 * m2->local_grid_size[0] / (stop_cycle - start_cycle) *
    CLOCKS_PER_SEC;
#endif

  /* print status message */
  printf("%08d: t=%5.4f dt=%5.4e kzps=%4.3f\n",
	 m2->status.iteration_number,
	 m2->status.time_simulation,
	 dt,
	 kzps);

  m2->status.iteration_number += 1;
  m2->status.time_simulation += dt;
}
