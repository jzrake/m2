#include <math.h>
#include <string.h>
#include "m2.h"

static void initial_data(m2vol *V)
{
  double m = 0.0;
  double face1[4] = {0.0,
		     V->global_index[1] + 0.5,
		     V->global_index[2],
		     V->global_index[3]};
  double face2[4] = {0.0,
		     V->global_index[1],
		     V->global_index[2] + 0.5,
		     V->global_index[3]};
  double x1[4]; /* position of r-face */
  double x2[4]; /* position of t-face */
  m2sim_index_to_position(V->m2, face1, x1);
  m2sim_index_to_position(V->m2, face2, x2);
  double R1 = x1[1];
  double T1 = x1[2];
  V->prim.d = 0.1;
  V->prim.p = 0.001;
  V->prim.v1 = 0.0;
  V->prim.v2 = 0.0;
  V->prim.v3 = 0.0;
  V->prim.B1 = m / pow(R1, 2) * ((T1 > M2_PI/2) - (T1 < M2_PI/2));
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
  double r, t, R, f;
  m2vol *V0;
#ifdef _OPENMP
#pragma omp parallel for private(V0,j,k) default(shared)
#endif
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {
	V0 = M2_VOL(i+0, j, k);
	I = V0->global_index;
	if (I[2] == G[2] - 1) {
	  initial_data(V0);
	  m2sim_from_primitive(m2,
			       &V0->prim, NULL, NULL,
			       V0 ->volume,
			       V0 ->consA,
			       &V0->aux);
	}
	else if (I[1] == 0) {
	  r = m2vol_coordinate_centroid(V0, 1);
	  t = m2vol_coordinate_centroid(V0, 2);
	  R = r * sin(t);
	  f = pow(cos(t), 4);
	  initial_data(V0);
	  V0->prim.d = 1.00;
	  V0->prim.p = 0.01;
	  V0->prim.v1 = 0.99 * f;
	  V0->prim.v2 = 0.0;
	  V0->prim.v3 = 0.0 * f * R;
	  V0->prim.B3 = 1.0 * sin(t) * cos(t);
	  V0->Bflux3A = V0->prim.B3 * V0->area3;
	  if (V0->prim.v1*V0->prim.v1 +
	      V0->prim.v2*V0->prim.v2 +
	      V0->prim.v3*V0->prim.v3 > 1.0) {
	    MSG(FATAL, "superluminal");
	  }
	  m2sim_from_primitive(m2,
			       &V0->prim, NULL, NULL,
			       V0 ->volume,
			       V0 ->consA,
			       &V0->aux);
	}
	else if (I[1] == G[1] - 1) {
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

void initialize_problem_jet(m2sim *m2)
{
  m2sim_set_resolution(m2, 512, 128, 1);
  m2sim_set_guard_zones(m2, 0);
  m2sim_set_extent0(m2, 1e0, 0.0, 0.0);
  m2sim_set_extent1(m2, 1e2, M2_PI/2, 2*M2_PI);
  m2sim_set_geometry(m2, M2_SPHERICAL);
  m2sim_set_physics(m2, M2_RELATIVISTIC | M2_MAGNETIZED);
  m2sim_set_ct_scheme(m2, M2_CT_FULL3D);
  m2sim_set_rk_order(m2, 2);
  m2sim_set_boundary_conditions(m2, boundary_conditions);
  m2sim_set_initial_data(m2, initial_data);

  m2->plm_parameter = 1.00;
  m2->cfl_parameter = 0.30;
  m2->simple_eigenvalues = 0;
  m2->interpolation_fields = M2_PRIMITIVE_AND_FOUR_VELOCITY;
  m2->coordinate_scaling1 = M2_LOGARITHMIC;
}
