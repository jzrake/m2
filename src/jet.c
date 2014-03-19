#include <math.h>
#include <string.h>
#include "m2.h"

#define AMBIENT_DENSITY  0.010
#define AMBIENT_PRESSURE 0.001
#define NU_PARAMETER 0.75 /* nu parameter */

static void initial_data(m2vol *V)
{
  V->prim.d = AMBIENT_DENSITY;
  V->prim.p = AMBIENT_PRESSURE;
  V->prim.v1 = 0.0;
  V->prim.v2 = 0.0;
  V->prim.v3 = 0.0;
  V->prim.B1 = 0.0;
  V->prim.B2 = 0.0;
  V->prim.B3 = 0.0;
  V->Bflux1A = V->prim.B1 * V->area1;
  V->Bflux2A = V->prim.B2 * V->area2;
  V->Bflux3A = V->prim.B3 * V->area3;
}

static void initial_data_cell(m2sim *m2, double X[4], double N[4], double *prim)
{
  prim[RHO] = AMBIENT_DENSITY;
  prim[PRE] = AMBIENT_PRESSURE;
  prim[V11] = 0.0;
  prim[V22] = 0.0;
  prim[V33] = 0.0;
}

static void initial_data_edge(m2sim *m2, double X[4], double N[4], double *E)
{
  double x[4];
  m2_to_cartesian(X, x, M2_PARABOLIC);
  double n = NU_PARAMETER;
  double R = sqrt(x[1]*x[1] + x[2]*x[2]);
  double r = sqrt(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
  double t = acos(x[3]/r);
  //double f = atan(x[2]/x[1]);
  double P = pow(r, n) * (1.0 - cos(t));
  double Af = 5.0 * P / (R + 0.1); /* A_phi */
  *E = -Af * N[3];
}

static void boundary_conditions(m2sim *m2)
{
  int *L = m2->local_grid_size;
  int *G = m2->domain_resolution;
  int *I;
  int i, j, k;
  m2vol *V0;
#ifdef _OPENMP
#pragma omp parallel for private(V0,j,k) default(shared)
#endif
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {
	V0 = M2_VOL(i, j, k);
	I = V0->global_index;

	if (I[2] == 0) {
	  double X[4], x[4];

	  m2vol_coordinate_centroid_3d(V0, X);
	  m2_to_cartesian(X, x, M2_PARABOLIC);

	  double n = NU_PARAMETER;
	  double r = sqrt(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
	  double cost = x[3]/r;
	  double rdisk = r * pow(1.0 - cost, 1.0/n);
	  double omega = rdisk < 1.0 ? 0.25 : 0.25 * pow(rdisk, -1.5);
	  double vkepl = omega * rdisk;
	  double vr = rdisk < 1.0 ? 0.5 : 0.5 * pow(rdisk, -0.5);

	  V0->prim.v1 = 0.00;
	  V0->prim.v2 = vr;
	  V0->prim.v3 = vkepl;
	  V0->prim.p = 0.01;
	  V0->prim.d = 1.00;
	  m2sim_from_primitive(m2,
			       &V0->prim, NULL, NULL,
			       V0 ->volume,
			       V0 ->consA,
			       &V0->aux);
	}
	else if (I[2] == G[2] - 1) {
	  initial_data(V0);
	  m2sim_from_primitive(m2,
			       &V0->prim, NULL, NULL,
			       V0 ->volume,
			       V0 ->consA,
			       &V0->aux);
	}
	if (I[1] == G[1] - 1) {
	  /* initial_data(V0); */
	  /* m2sim_from_primitive(m2, */
	  /* 		       &V0->prim, NULL, NULL, */
	  /* 		       V0 ->volume, */
	  /* 		       V0 ->consA, */
	  /* 		       &V0->aux); */
	}
      }
    }
  }
}

static void boundary_conditions_emf1(m2vol *V)
{
  V->emf1 = 0.0;
}
static void boundary_conditions_emf2(m2vol *V)
{
  V->emf2 = 0.0;
}
static void boundary_conditions_emf3(m2vol *V)
{
  V->emf3 = 0.0;
}
static void boundary_conditions_flux1(m2vol *V)
{
  double n[4] = {0.0, 1.0, 0.0, 0.0};
  m2vol *V1 = V->global_index[1] == -1 ? m2vol_neighbor(V, 1, 1) : V;
  m2aux_fluxes(&V1->aux, n, V->flux1);
}
static void boundary_conditions_flux2(m2vol *V)
{
  double n[4] = {0.0, 0.0, 1.0, 0.0};
  m2vol *V1 = V->global_index[2] == -1 ? m2vol_neighbor(V, 2, 1) : V;
  m2aux_fluxes(&V1->aux, n, V->flux2);
}

void initialize_problem_jet(m2sim *m2)
{
  m2sim_set_resolution(m2, 128, 128, 1);
  m2sim_set_extent0(m2, 0.01, 4.0, 0.0);
  m2sim_set_extent1(m2, 8.0, 16.0, 2.0*M2_PI);
  m2->number_guard_zones0[1] = 1;
  m2->number_guard_zones0[2] = 1;
  //m2->geometry = M2_SPHERICAL;
  m2->geometry = M2_PARABOLIC;
  //m2->coordinate_scaling1 = M2_LOGARITHMIC;
  //m2->coordinate_scaling2 = M2_LOGARITHMIC;
  m2->physics = M2_RELATIVISTIC | M2_MAGNETIZED;
  m2->rk_order = 2;
  m2->boundary_conditions = boundary_conditions;
  m2->boundary_conditions_flux1 = boundary_conditions_flux1;
  m2->boundary_conditions_flux2 = boundary_conditions_flux2;
  m2->boundary_conditions_emf1 = boundary_conditions_emf1;
  m2->boundary_conditions_emf2 = boundary_conditions_emf2;
  m2->boundary_conditions_emf3 = boundary_conditions_emf3;
  m2->initial_data_cell = initial_data_cell;
  m2->initial_data_edge = initial_data_edge;
  m2->plm_parameter = 1.25;
  m2->cfl_parameter = 0.30;
  m2->simple_eigenvalues = 0;
  m2->interpolation_fields = M2_PRIMITIVE_AND_FOUR_VELOCITY;
  m2->status.checkpoint_cadence = 10;
}
