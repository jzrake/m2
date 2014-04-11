#include <string.h>
#include <math.h>
#if _OPENMP
#include <omp.h>
#endif
#include "m2.h"
#include "riemann.h"

#define ENABLE_PLM 1


static void riemann_solver(m2vol *VL, m2vol *VR, int axis, double *F)
{
  int q;
  m2riemann_problem R = { .m2 = VL->m2 };
  double *UL = R.Ul;
  double *UR = R.Ur;
  double *FL = R.Fl;
  double *FR = R.Fr;
  double *lamL = R.lamL;
  double *lamR = R.lamR;
  double yl[8], yr[8]; /* cell-centered left/right interpolation variables */
  double yL[8], yR[8]; /* face-centered left/right interpolation variables */
  double xL = m2vol_coordinate_centroid(VL, axis);
  double xR = m2vol_coordinate_centroid(VR, axis);
  double *n = R.nhat;
  m2prim *PL = &R.Pl;
  m2prim *PR = &R.Pr;
  m2aux *AL = &R.Al;
  m2aux *AR = &R.Ar;

  n[0] = 0.0;
  n[1] = 0.0;
  n[2] = 0.0;
  n[3] = 0.0;

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

  m2sim_from_interpolated(VL->m2, yL, PL);
  m2sim_from_interpolated(VR->m2, yR, PR);

  UL[B11] = PL->B1;
  UL[B22] = PL->B2;
  UL[B33] = PL->B3;
  UR[B11] = PR->B1;
  UR[B22] = PR->B2;
  UR[B33] = PR->B3;

  switch (axis) {
  case 1: UL[B11] = UR[B11] = PL->B1 = PR->B1 = VL->Bflux1A / VL->area1; break;
  case 2: UL[B22] = UR[B22] = PL->B2 = PR->B2 = VL->Bflux2A / VL->area2; break;
  case 3: UL[B33] = UR[B33] = PL->B3 = PR->B3 = VL->Bflux3A / VL->area3; break;
  }

  m2sim_from_primitive(VL->m2, PL, NULL, NULL, 1.0, UL, AL);
  m2sim_from_primitive(VR->m2, PR, NULL, NULL, 1.0, UR, AR);
  m2aux_eigenvalues(AL, n, lamL);
  m2aux_eigenvalues(AR, n, lamR);
  m2aux_fluxes(AL, n, FL);
  m2aux_fluxes(AR, n, FR);

  R.Sl = M2_MIN3(lamL[0], lamR[0], 0.0);
  R.Sr = M2_MAX3(lamL[7], lamR[7], 0.0);

  switch (VL->m2->riemann_solver) {
  case M2_RIEMANN_HLLE: riemann_solver_hlle(&R); break;
  case M2_RIEMANN_HLLC: riemann_solver_hllc(&R); break;
  default: MSG(FATAL, "invalid riemann solver"); break;
  }
  memcpy(F, R.F_hat, 8 * sizeof(double));
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
  int i, j, k;
  int *L = m2->local_grid_size;
  int *G = m2->domain_resolution;
  m2vol *V0, *V1, *V2, *V3;

  double n1[4] = { 0.0, 1.0, 0.0, 0.0 };
  double n2[4] = { 0.0, 0.0, 1.0, 0.0 };
  double n3[4] = { 0.0, 0.0, 0.0, 1.0 };

#ifdef _OPENMP
#pragma omp parallel for private(V0,V1,V2,V3,j,k)
#endif
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {

	V0 = M2_VOL(i+0, j+0, k+0);
	V1 = M2_VOL(i+1, j+0, k+0);
	V2 = M2_VOL(i+0, j+1, k+0);
	V3 = M2_VOL(i+0, j+0, k+1);

	if (V1 == NULL) m2aux_fluxes(&V0->aux, n1, V0->flux1);
	else riemann_solver(V0, V1, 1, V0->flux1);

	if (V2 == NULL) m2aux_fluxes(&V0->aux, n2, V0->flux2);
	else riemann_solver(V0, V2, 2, V0->flux2);

	if (V3 == NULL) m2aux_fluxes(&V0->aux, n3, V0->flux3);
	else riemann_solver(V0, V3, 3, V0->flux3);
      }
    }
  }

  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {
	V0 = M2_VOL(i, j, k);
	if ((V0->global_index[1] == -1 ||
	     V0->global_index[1] == G[1] - 1) &&
	    (m2->boundary_conditions_flux1)) {
	  m2->boundary_conditions_flux1(V0);
	}
	if ((V0->global_index[2] == -1 ||
	     V0->global_index[2] == G[2] - 1) &&
	    (m2->boundary_conditions_flux2)) {
	  m2->boundary_conditions_flux2(V0);
	}
	if ((V0->global_index[3] == -1 ||
	     V0->global_index[3] == G[3] - 1) &&
	    (m2->boundary_conditions_flux3)) {
	  m2->boundary_conditions_flux3(V0);
	}
      }
    }
  }
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

static void _calculate_emf1(m2sim *m2)
{
  int i;
  int *L = m2->local_grid_size;
  m2vol *V0;
  for (i=0; i<L[1]; ++i) {
    V0 = M2_VOL(i, 0, 0);
    V0->emf1 = 0.0;
    V0->emf2 = +V0->flux1[B33] * V0->line2;
    V0->emf3 = -V0->flux1[B22] * V0->line3;
  }
}

typedef struct
{
  double *F1;
  double *F2;
  double *x0;
  double *x1;
  double v1, v2, B1, B2;
} m2emf_descr;

static double _calculate_emf_single(m2emf_descr emfs[4], int axis)
{
  int a1 =   1 + (axis + 0) % 3;
  int a2 =   1 + (axis + 1) % 3;
  int b1 = B11 + (axis + 0) % 3;
  int b2 = B11 + (axis + 1) % 3;

  m2emf_descr *emf00 = &emfs[0];
  m2emf_descr *emf10 = &emfs[1];
  m2emf_descr *emf01 = &emfs[2];
  m2emf_descr *emf11 = &emfs[3];

  double dx = emf00->x1[a1] - emf00->x0[a1];
  double dy = emf00->x1[a2] - emf00->x0[a2];

  double Eh0 = -emf00->F1[b2]; /* E_{i+1/2,j+0} */
  double Eh1 = -emf01->F1[b2]; /* E_{i+1/2,j+1} */
  double E0h = +emf00->F2[b1]; /* E_{i+0,j+1/2} */
  double E1h = +emf10->F2[b1]; /* E_{i+1,j+1/2} */
  double E00 = -(emf00->v1*emf00->B2 - emf00->v2*emf00->B1); /* E_{i+0,j+0} */
  double E10 = -(emf10->v1*emf10->B2 - emf10->v2*emf10->B1); /* E_{i+1,j+0} */
  double E01 = -(emf01->v1*emf01->B2 - emf01->v2*emf01->B1); /* E_{i+0,j+1} */
  double E11 = -(emf11->v1*emf11->B2 - emf11->v2*emf11->B1); /* E_{i+1,j+1} */

  double dEdx0 = 0.0, dEdx1 = 0.0;
  double dEdy0 = 0.0, dEdy1 = 0.0;

  /* dE/dx @ {i+0,j+1/2} */
  if      (emf00->F2[DDD] == 0.0) dEdx0 = (Eh1 - E01)/dx + (Eh0 - E00)/dx;
  else if (emf00->F2[DDD]  < 0.0) dEdx0 = (Eh1 - E01)/dx * 2;
  else if (emf00->F2[DDD]  > 0.0) dEdx0 =              2 * (Eh0 - E00)/dx;

  /* dE/dx @ {i+1,j+1/2} */
  if      (emf10->F2[DDD] == 0.0) dEdx1 = (E11 - Eh1)/dx + (E10 - Eh0)/dx;
  else if (emf10->F2[DDD]  < 0.0) dEdx1 = (E11 - Eh1)/dx * 2;
  else if (emf10->F2[DDD]  > 0.0) dEdx1 =              2 * (E10 - Eh0)/dx;

  /* dE/dy @ {i+1/2,j+0} */
  if      (emf00->F1[DDD] == 0.0) dEdy0 = (E1h - E10)/dy + (E0h - E00)/dy;
  else if (emf00->F1[DDD]  < 0.0) dEdy0 = (E1h - E10)/dy * 2;
  else if (emf00->F1[DDD]  > 0.0) dEdy0 =              2 * (E0h - E00)/dy;

  /* dE/dy @ {i+1/2,j+1} */
  if      (emf01->F1[DDD] == 0.0) dEdy1 = (E11 - E1h)/dy + (E01 - E0h)/dy;
  else if (emf01->F1[DDD]  < 0.0) dEdy1 = (E11 - E1h)/dy * 2;
  else if (emf01->F1[DDD]  > 0.0) dEdy1 =              2 * (E01 - E0h)/dy;

  return 0.25 * (E0h + 0.5*dx*dEdx0 +
		 E1h - 0.5*dx*dEdx1 +
		 Eh0 + 0.5*dy*dEdy0 +
		 Eh1 - 0.5*dy*dEdy1);
}

static void _calculate_emf2(m2sim *m2)
{
  int i, j, n;
  int *L = m2->local_grid_size;
  int bnd;
  m2vol *vols[4];
  m2emf_descr emfs[4];

  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {

      vols[0] = M2_VOL(i,   j,   0);
      vols[1] = M2_VOL(i+1, j,   0);
      vols[2] = M2_VOL(i,   j+1, 0);
      vols[3] = M2_VOL(i+1, j+1, 0);

      if (vols[1] == NULL ||
	  vols[2] == NULL ||
	  vols[3] == NULL) {
	bnd = 1;
      }
      else {
	for (n=0; n<4; ++n) {
	  emfs[n].F1 = vols[n]->flux1;
	  emfs[n].F2 = vols[n]->flux2;
	  emfs[n].x0 = vols[n]->x0;
	  emfs[n].x1 = vols[n]->x1;
	  emfs[n].v1 = vols[n]->prim.v1;
	  emfs[n].v2 = vols[n]->prim.v2;
	  emfs[n].B1 = vols[n]->prim.B1;
	  emfs[n].B2 = vols[n]->prim.B2;
	}
	bnd = 0;
      }

      vols[0]->emf1 = -vols[0]->flux2[B33] * vols[0]->line1;
      vols[0]->emf2 = +vols[0]->flux1[B33] * vols[0]->line2;
      vols[0]->emf3 = bnd?0.0:_calculate_emf_single(emfs, 3) * vols[0]->line3;
    }
  }
}

static void _calculate_emf3(m2sim *m2)
{
  int i, j, k, n, m;
  int *L = m2->local_grid_size;
  m2vol *vols[4][4];
  m2emf_descr emfs[4][4];

  for (i=0; i<L[1]-1; ++i) {
    for (j=0; j<L[2]-1; ++j) {
      for (k=0; k<L[3]-1; ++k) {

	vols[1][0] = M2_VOL(i, j,   k  );
	vols[1][1] = M2_VOL(i, j+1, k  );
	vols[1][2] = M2_VOL(i, j,   k+1);
	vols[1][3] = M2_VOL(i, j+1, k+1);
	vols[2][0] = M2_VOL(i,   j, k  );
	vols[2][1] = M2_VOL(i,   j, k+1);
	vols[2][2] = M2_VOL(i+1, j, k  );
	vols[2][3] = M2_VOL(i+1, j, k+1);
	vols[3][0] = M2_VOL(i,   j,   k);
	vols[3][1] = M2_VOL(i+1, j,   k);
	vols[3][2] = M2_VOL(i,   j+1, k);
	vols[3][3] = M2_VOL(i+1, j+1, k);

	for (m=1; m<4; ++m) {
	  for (n=0; n<4; ++n) {
	    emfs[m][n].x0 = vols[m][n]->x0;
	    emfs[m][n].x1 = vols[m][n]->x1;
	    switch (m) {
	    case 1:
	      emfs[m][n].F1 = vols[m][n]->flux2;
	      emfs[m][n].F2 = vols[m][n]->flux3;
	      emfs[m][n].v1 = vols[m][n]->prim.v2;
	      emfs[m][n].v2 = vols[m][n]->prim.v3;
	      emfs[m][n].B1 = vols[m][n]->prim.B2;
	      emfs[m][n].B2 = vols[m][n]->prim.B3;
	      break;
	    case 2:
	      emfs[m][n].F1 = vols[m][n]->flux3;
	      emfs[m][n].F2 = vols[m][n]->flux1;
	      emfs[m][n].v1 = vols[m][n]->prim.v3;
	      emfs[m][n].v2 = vols[m][n]->prim.v1;
	      emfs[m][n].B1 = vols[m][n]->prim.B3;
	      emfs[m][n].B2 = vols[m][n]->prim.B1;
	      break;
	    case 3:
	      emfs[m][n].F1 = vols[m][n]->flux1;
	      emfs[m][n].F2 = vols[m][n]->flux2;
	      emfs[m][n].v1 = vols[m][n]->prim.v1;
	      emfs[m][n].v2 = vols[m][n]->prim.v2;
	      emfs[m][n].B1 = vols[m][n]->prim.B1;
	      emfs[m][n].B2 = vols[m][n]->prim.B2;
	      break;
	    }
	  }
	}
	vols[1][0]->emf1 = _calculate_emf_single(emfs[1], 1)*vols[1][0]->line1;
	vols[1][0]->emf2 = _calculate_emf_single(emfs[2], 2)*vols[1][0]->line2;
	vols[1][0]->emf3 = _calculate_emf_single(emfs[3], 3)*vols[1][0]->line3;
      }
    }
  }
}

static void _exchange_flux1(m2sim *m2, double dt)
{
  int i, q;
  int *L = m2->local_grid_size;
  m2vol *V0, *V1;

  for (i=0; i<L[1]; ++i) {
    V0 = M2_VOL(i,   0, 0);
    V1 = M2_VOL(i-1, 0, 0);

    V1 = V1 ? V1 : V0;

    for (q=0; q<5; ++q) {
      V0->consA[q] -= dt * V0->area1 * V0->flux1[q];
      V0->consA[q] += dt * V1->area1 * V1->flux1[q];
    }

    if (m2->magnetized) {
      V0->Bflux2A -= dt * V1->emf3;
      V0->Bflux2A += dt * V0->emf3;
      V0->Bflux3A -= dt * V0->emf2;
      V0->Bflux3A += dt * V1->emf2;
    }
  }
}

static void _exchange_flux2(m2sim *m2, double dt)
{
  int i, j, q;
  int *L = m2->local_grid_size;
  m2vol *V0, *V1, *V2;

  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      V0 = M2_VOL(i,   j, 0);
      V1 = M2_VOL(i-1, j, 0);
      V2 = M2_VOL(i, j-1, 0);

      V1 = V1 ? V1 : V0;
      V2 = V2 ? V2 : V0;

      for (q=0; q<5; ++q) {
	V0->consA[q] -= dt * V0->area1 * V0->flux1[q];
	V0->consA[q] += dt * V1->area1 * V1->flux1[q];
	V0->consA[q] -= dt * V0->area2 * V0->flux2[q];
	V0->consA[q] += dt * V2->area2 * V2->flux2[q];
      }

      if (m2->magnetized) {
	V0->Bflux1A -= dt * V0->emf3;
	V0->Bflux1A += dt * V2->emf3;
	V0->Bflux2A += dt * V0->emf3;
	V0->Bflux2A -= dt * V1->emf3;

	V0->Bflux3A -= dt * V0->emf2;
	V0->Bflux3A += dt * V1->emf2;
	V0->Bflux3A += dt * V0->emf1;
	V0->Bflux3A -= dt * V2->emf1;
      }
    }
  }
}

static void _exchange_flux3(m2sim *m2, double dt)
{
  int i, j, k, q;
  int *L = m2->local_grid_size;
  m2vol *V0, *V1, *V2, *V3;

  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {

	V0 = M2_VOL(i,   j, k  );
	V1 = M2_VOL(i-1, j, k  );
	V2 = M2_VOL(i, j-1, k  );
	V3 = M2_VOL(i, j,   k-1);

	V1 = V1 ? V1 : V0;
	V2 = V2 ? V2 : V0;
	V3 = V3 ? V3 : V0;

	for (q=0; q<5; ++q) {
	  V0->consA[q] -= dt * V0->area1 * V0->flux1[q];
	  V0->consA[q] += dt * V1->area1 * V1->flux1[q];
	  V0->consA[q] -= dt * V0->area2 * V0->flux2[q];
	  V0->consA[q] += dt * V2->area2 * V2->flux2[q];
	  V0->consA[q] -= dt * V0->area3 * V0->flux3[q];
	  V0->consA[q] += dt * V3->area3 * V3->flux3[q];
	}

	if (m2->magnetized) {
	  V0->Bflux1A += dt * V2->emf3;
	  V0->Bflux1A -= dt * V0->emf3;
	  V0->Bflux1A += dt * V0->emf2;
	  V0->Bflux1A -= dt * V3->emf2;

	  V0->Bflux2A += dt * V3->emf1;
	  V0->Bflux2A -= dt * V0->emf1;
	  V0->Bflux2A += dt * V0->emf3;
	  V0->Bflux2A -= dt * V1->emf3;

	  V0->Bflux3A += dt * V1->emf2;
	  V0->Bflux3A -= dt * V0->emf2;
	  V0->Bflux3A += dt * V0->emf1;
	  V0->Bflux3A -= dt * V2->emf1;
	}
      }
    }
  }
}

void m2sim_calculate_emf(m2sim *m2)
{
  if (!m2->magnetized) return;
  int *G = m2->domain_resolution;
  int ndim = 0;
  ndim += G[1] > 1;
  ndim += G[2] > 1;
  ndim += G[3] > 1;
  switch (ndim) {
  case 1: _calculate_emf1(m2); break;
  case 2: _calculate_emf2(m2); break;
  case 3: _calculate_emf3(m2); break;
  }
}

void m2sim_exchange_flux(m2sim *m2, double dt)
{
  int *G = m2->domain_resolution;
  int ndim = 0;
  ndim += G[1] > 1;
  ndim += G[2] > 1;
  ndim += G[3] > 1;
  switch (ndim) {
  case 1: _exchange_flux1(m2, dt); break;
  case 2: _exchange_flux2(m2, dt); break;
  case 3: _exchange_flux3(m2, dt); break;
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

	/* if (V1 && V2 && V3) { */
	/*   double divB = */
	/*     (VC->Bflux1A * VC->area1 - V1->Bflux1A * V1->area1 +  */
	/*      VC->Bflux2A * VC->area2 - V2->Bflux2A * V2->area2 +  */
	/*      VC->Bflux3A * VC->area3 - V3->Bflux3A * V3->area3); */
	/*   if (fabs(divB) > 1e-16) { */
	/*     printf("large divB: %16.12e\n", divB); */
	/*   } */
	/* } */
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
  int n, err;
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
      dt = m2vol_minimum_dimension(V) / m2aux_maximum_wavespeed(&V->aux, &err);
      if (err) {
	MSG(INFO, "got bad wavespeed");
	m2_print_state(NULL, &V->aux, NULL);
	MSGF(FATAL, "at global index [%d %d %d]",
	     V->global_index[1],
	     V->global_index[2],
	     V->global_index[3]);
      }
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
  if (m2->boundary_conditions_cell == NULL) {

  }
  else {
    m2->boundary_conditions_cell(m2);
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
  dt = m2->cfl_parameter * m2sim_minimum_courant_time(m2);
  m2->status.time_step = dt;

  double cad = m2->cadence_checkpoint_tpl;
  double now = m2->status.time_simulation;
  double las = m2->status.time_last_checkpoint_tpl;

  if (cad > 0.0 && now - las > cad) {
    m2->status.checkpoint_number_tpl += 1;
    m2->status.time_last_checkpoint_tpl += cad;
    snprintf(checkpoint_name, sizeof(checkpoint_name), "chkpt.%04d.m2",
	     m2->status.checkpoint_number_tpl);
    m2sim_save_checkpoint_tpl(m2, checkpoint_name);
  }

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

  m2->status.iteration_number += 1;
  m2->status.time_simulation += dt;

  /* print status message */
  snprintf(m2->status.message, 1024,
	   "%08d: t=%5.4f dt=%5.4e kzps=%4.3f",
	   m2->status.iteration_number,
	   m2->status.time_simulation,
	   dt,
	   kzps);
}
