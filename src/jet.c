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

  double rj = 0.2; /* jet width */
  double Bj = 1.0; /* B-phi in the jet */
  double dj = 1.00; /* density */
  double pj = 0.01; /* gas pressure at the jet edge */
  double gj = 5.00; /* jet Lorentz factor */

#ifdef _OPENMP
#pragma omp parallel for private(V0,j,k) default(shared)
#endif
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {
	V0 = M2_VOL(i+0, j, k);
	I = V0->global_index;
	if (I[1] == G[1] - 1) {
	  initial_data(V0);
	  m2sim_from_primitive(m2,
			       &V0->prim, NULL, NULL,
			       V0 ->volume,
			       V0 ->consA,
			       &V0->aux);
	}
	else if (I[2] == 0) {
	  initial_data(V0);

	  double xp[4], xc[4];
	  m2vol_coordinate_centroid_3d(V0, xp);
	  m2_to_cartesian(xp, xc, M2_PARABOLIC);
	  double rc = sqrt(xc[1]*xc[1] + xc[2]*xc[2]);

	  if (rc < rj) {
	    V0->prim.v2 = sqrt(1.0 - 1.0/(gj*gj));
	    V0->prim.d = dj;
	    V0->prim.p = pj - Bj * Bj * (rc/rj - 1.0);
	    V0->prim.B3 = Bj * sqrt(rc/rj);
	    V0->Bflux3A = V0->prim.B3 * V0->area3;
	  }
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
      }
    }
  }
}

void initialize_problem_jet(m2sim *m2)
{
  m2sim_set_resolution(m2, 64, 512, 1);
  m2sim_set_extent0(m2, 0.0,  1.0, 0.0);
  m2sim_set_extent1(m2, 2.0, 16.0, 2.0*M2_PI);
  m2sim_set_geometry(m2, M2_PARABOLIC);
  m2sim_set_physics(m2, M2_RELATIVISTIC | M2_MAGNETIZED);
  m2sim_set_rk_order(m2, 2);
  m2sim_set_boundary_conditions(m2, boundary_conditions);
  m2sim_set_initial_data(m2, initial_data);

  m2->plm_parameter = 1.25;
  m2->cfl_parameter = 0.30;
  m2->simple_eigenvalues = 0;
  m2->interpolation_fields = M2_PRIMITIVE_AND_FOUR_VELOCITY;
}
