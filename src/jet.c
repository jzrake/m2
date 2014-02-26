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
  printf("called BC's on sim\n");
}

static void boundary_conditions_flux1(m2vol *V)
{
  if (V->global_index[1] == -1) {
    printf("called BC's on flux1: %d %d\n",
	   V->global_index[1],
	   V->global_index[2]);
  }
}

static void boundary_conditions_flux2(m2vol *V)
{
  if (V->global_index[1] == -1) {
    printf("called BC's on flux2: %d %d\n",
	   V->global_index[1],
	   V->global_index[2]);
  }
  /* printf("called BC's on flux2: %d %d\n", */
  /* 	 V->global_index[1], */
  /* 	 V->global_index[2]); */
}

void initialize_problem_jet(m2sim *m2)
{
  //m2sim_set_resolution(m2, 64, 512, 1);
  //m2sim_set_extent0(m2, 0.0,  1.0, 0.0);
  //m2sim_set_extent1(m2, 2.0, 16.0, 2.0*M2_PI);
  //m2sim_set_geometry(m2, M2_PARABOLIC);

  m2sim_set_resolution(m2, 256, 128, 1);
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

  m2->plm_parameter = 1.25;
  m2->cfl_parameter = 0.30;
  m2->simple_eigenvalues = 0;
  m2->interpolation_fields = M2_PRIMITIVE_AND_FOUR_VELOCITY;
}
