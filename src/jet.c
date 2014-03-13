#include <math.h>
#include <string.h>
#include "m2.h"

static void initial_data(m2vol *V)
{
  V->prim.d = 0.01;
  V->prim.p = 0.0001;
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
	V0 = M2_VOL(i+0, j, k);
	I = V0->global_index;
	if (I[1] == 0) {
	  /* initial_data(V0); */
	  /* m2sim_from_primitive(m2, */
	  /* 		       &V0->prim, NULL, NULL, */
	  /* 		       V0 ->volume, */
	  /* 		       V0 ->consA, */
	  /* 		       &V0->aux); */
	}
	else if (I[1] == G[1] - 1) {
	  initial_data(V0);
	  m2sim_from_primitive(m2,
			       &V0->prim, NULL, NULL,
			       V0 ->volume,
			       V0 ->consA,
			       &V0->aux);
	}
	if (I[2] == G[2] - 1) {
	  initial_data(V0);
	  m2sim_from_primitive(m2,
			       &V0->prim, NULL, NULL,
			       V0 ->volume,
			       V0 ->consA,
			       &V0->aux);
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
  int q;
  double t = m2vol_coordinate_centroid(V, 2);
  for (q=0; q<8; ++q) V->flux1[q] = 0.0;
  V->flux1[DDD] =  0.1 * (t<0.2);
  V->flux1[S11] = 10.0 * (t<0.2);
  V->flux1[TAU] = 10.0 * (t<0.2);
}
static void boundary_conditions_flux2(m2vol *V)
{
  int q;
  for (q=0; q<8; ++q) V->flux2[q] = 0.0;
}

void initialize_problem_jet(m2sim *m2)
{
  m2sim_set_resolution(m2, 128, 64, 1);
  m2sim_set_extent0(m2, 1.0, 0.0, 0.0);
  m2sim_set_extent1(m2, 10.0, 0.5*M2_PI, 2.0*M2_PI);

  m2->number_guard_zones0[1] = 1;
  m2->number_guard_zones0[2] = 1;

  m2->geometry = M2_SPHERICAL;
  m2->coordinate_scaling1 = M2_LOGARITHMIC;
  m2->physics = M2_RELATIVISTIC | M2_MAGNETIZED;
  m2->rk_order = 2;
  m2->initial_data = initial_data;
  m2->boundary_conditions = boundary_conditions;
  m2->boundary_conditions_flux1 = boundary_conditions_flux1;
  m2->boundary_conditions_flux2 = boundary_conditions_flux2;
  m2->boundary_conditions_emf1 = boundary_conditions_emf1;
  m2->boundary_conditions_emf2 = boundary_conditions_emf2;
  m2->boundary_conditions_emf3 = boundary_conditions_emf3;

  m2->plm_parameter = 1.25;
  m2->cfl_parameter = 0.30;
  m2->simple_eigenvalues = 0;
  m2->interpolation_fields = M2_PRIMITIVE_AND_FOUR_VELOCITY;
}
