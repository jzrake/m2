#include <math.h>
#include "m2.h"

static void initial_data(m2vol *V)
{
  double x = m2vol_coordinate_centroid(V, 1);

  if (x < 0.5) {
    V->prim.v1 = 0.0;
    V->prim.v2 = 0.0;
    V->prim.v3 = 0.0;
    V->prim.d  = 1.0;
    V->prim.p  = 1.0;
  }
  else {
    V->prim.v1 = 0.0;
    V->prim.v2 = 0.0;
    V->prim.v3 = 0.0;
    V->prim.d  = 0.1;
    V->prim.p  = 0.125;
  }

  V->Bflux1A = 1.0 * V->area1;
  V->Bflux2A = 1.0 * V->area2;
  V->Bflux3A = 0.0 * V->area3;
}

static void boundary_conditions(m2sim *m2)
{
  int n;
  int *L = m2->local_grid_size;
  int *G = m2->domain_resolution;
  int *I;
  m2vol *V;
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    I = V->global_index;
    if (I[1] < 0 || I[1] >= G[1]) {
      initial_data(V);
      m2sim_from_primitive(m2,
                           &V->prim, NULL, NULL,
                           V ->volume,
                           V ->consA,
                           &V->aux);
    }
  }
}

void initialize_problem_shocktube(m2sim *m2)
{
  m2sim_set_resolution(m2, 256, 1, 1);
  m2sim_set_extent0(m2, 0.0, 0.0, 0.0);
  m2sim_set_extent1(m2, 1.0, 1.0, 1.0);
  m2sim_set_geometry(m2, M2_CARTESIAN);
  m2sim_set_physics(m2, M2_RELATIVISTIC | M2_MAGNETIZED);
  m2sim_set_ct_scheme(m2, M2_CT_FULL3D);
  m2sim_set_analysis(m2, NULL);
  m2sim_set_boundary_conditions(m2, boundary_conditions);
  m2sim_set_initial_data(m2, initial_data);

  m2->number_guard_zones0[1] = 2;
  m2->number_guard_zones1[1] = 2;
  m2->plm_parameter = 2.00;
  m2->cfl_parameter = 0.40;
  m2->simple_eigenvalues = 1;
  m2->interpolation_fields = M2_PRIMITIVE_AND_FOUR_VELOCITY;
}
