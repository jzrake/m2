#include <math.h>
#include <string.h>
#include "m2.h"

static double Bfield(double r, double t, double f, int a)
{
  double m = 0.1;
  switch (a) {
  case 1: return m / pow(r, 3) * 2 * cos(t);
  case 2: return m / pow(r, 3) * 1 * sin(t);
  case 3: return 0.0;
  }
  return 0.0;
}

static void initial_data(m2vol *V)
{
  double R0 = 0.5*(V->x0[1] + V->x1[1]);
  double T0 = 0.5*(V->x0[2] + V->x1[2]);
  double F0 = 0.5*(V->x0[3] + V->x1[3]);
  double R1 = V->x1[1];
  double T1 = V->x1[2];
  double F1 = V->x1[3];
  V->prim.d = 1.00;
  V->prim.p = 1.00;
  V->prim.v1 = 0.0;
  V->prim.v2 = 0.0;
  V->prim.v3 = 0.0;
  V->prim.B1 = Bfield(R0, T0, F0, 1);
  V->prim.B2 = Bfield(R0, T0, F0, 2);
  V->prim.B3 = Bfield(R0, T0, F0, 3);
  V->Bflux1A = Bfield(R1, T0, F0, 1) * V->area1;
  V->Bflux2A = Bfield(R0, T1, F0, 2) * V->area2;
  V->Bflux3A = Bfield(R0, T0, F1, 3) * V->area3;
}

static void boundary_conditions(m2sim *m2)
{
  int *L = m2->local_grid_size;
  int *G = m2->domain_resolution;
  int *I;
  int i, j, k;
  m2vol *V0, *V1;
#ifdef _OPENMP
#pragma omp parallel for private(V0,V1,j,k) default(shared)
#endif
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {
	V0 = M2_VOL(i+0, j, k);
	I = V0->global_index;
	if (I[1] == 0) {
	  V1 = M2_VOL(i+1, j, k);
	  initial_data(V0);
	  V0->prim.d = V1->prim.d;
	  V0->prim.p = V1->prim.p;
	  V0->prim.v1 *= -1.0;
	  V0->prim.v2 = 0.0;
	  m2sim_from_primitive(m2,
			       &V0->prim, NULL, NULL,
			       V0 ->volume,
			       V0 ->consA,
			       &V0->aux);
	}
	else if (I[1] == G[1] - 1) {
	  V1 = M2_VOL(i-1, j, k);
	  initial_data(V0);
	  V0->prim.d = V1->prim.d;
	  V0->prim.p = V1->prim.p;
	  V0->prim.v1 = V1->prim.v1;
	  V0->prim.v2 = 0.0;
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


static void add_physical_source_terms(m2vol *V)
{
  double r = m2vol_coordinate_centroid(V, 1);
  double dt = V->x1[0] - V->x0[0];
  double dV = V->volume;
  double vr = V->prim.v1;
  double d0 = V->prim.d;
  double Fr = -1.0 / (r*r);
  V->consA[S11] += dt * dV * d0 * Fr;
  V->consA[TAU] += dt * dV * d0 * Fr * vr;
}

void initialize_problem_flux_burial(m2sim *m2)
{
  m2sim_set_resolution(m2, 64*2, 128*2, 1);
  m2sim_set_extent0(m2, 1.0, 0.0, 0.0);
  m2sim_set_extent1(m2, 4.0, M2_PI, 2*M2_PI);
  m2sim_set_geometry(m2, M2_SPHERICAL);
  m2sim_set_rk_order(m2, 2);
  m2sim_set_boundary_conditions(m2, boundary_conditions);
  m2sim_set_initial_data(m2, initial_data);

  m2->relativistic = 0;
  m2->magnetized = 1;
  m2->number_guard_zones0[1] = 1;
  m2->number_guard_zones0[2] = 1;
  m2->add_physical_source_terms = add_physical_source_terms;
  m2->plm_parameter = 1.25;
  m2->cfl_parameter = 0.30;
  m2->simple_eigenvalues = 0;
  m2->interpolation_fields = M2_PRIMITIVE_AND_FOUR_VELOCITY;
  m2->coordinate_scaling1 = M2_LOGARITHMIC;
}
