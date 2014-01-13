#include <math.h>
#include "m2.h"

#define FNAME_VOLUME_INTEGRALS_I "volume-integrals-in.dat"
#define FNAME_VOLUME_INTEGRALS_O "volume-integrals-out.dat"

void mwn_initial_data(m2vol *V)
{
  double R = m2vol_coordinate_centroid(V, 1);
  if (R < 1.0) {
    V->prim.v1 = 0.0;
    V->prim.v2 = 0.0;
    V->prim.v3 = 0.0;
    V->prim.d  = 1.0;
    V->prim.p  = 0.001;

    V->Bflux1A = 0.0 * V->area1;
    V->Bflux2A = 0.0 * V->area2;
    V->Bflux3A = 0.0 * V->area3;
  }
  else {
    V->prim.v1 = 0.0;
    V->prim.v2 = 0.0;
    V->prim.v3 = 0.0;
    V->prim.d  = 1000.0;
    V->prim.p  = 0.001;

    V->Bflux1A = 0.0 * V->area1;
    V->Bflux2A = 0.0 * V->area2;
    V->Bflux3A = 0.0 * V->area3;
  }
}

void mwn_boundary_conditions(m2sim *m2)
{
  int n;
  int *L = m2->local_grid_size;
  int *G = m2->domain_resolution;
  int *I;
  m2vol *V;
  double r, t;
  double a = 1.0;
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    I = V->global_index;

    if (I[1] == 0) {
      r = m2vol_coordinate_centroid(V, 1);
      t = m2vol_coordinate_centroid(V, 2);
      V->prim.v1 = 1.0 * sqrt(a + (1 - a) * sin(t) * sin(t));
      V->prim.v2 = 0.0;
      V->prim.v3 = 0.0;
      V->prim.d = 1.00;
      V->prim.p = 0.001;
      V->Bflux1A = 0.0;
      V->Bflux2A = 0.0;
      V->Bflux3A = sqrt(0.1) * sin(t) * V->area3;
      m2sim_from_primitive(m2,
                           &V->prim, NULL, NULL,
                           V ->volume,
                           V ->consA,
                           &V->aux);
    }
    else if (I[1] == G[1] - 1) {
      mwn_initial_data(V);
      m2sim_from_primitive(m2,
                           &V->prim, NULL, NULL,
                           V ->volume,
                           V ->consA,
                           &V->aux);
    }
    else if (I[2] == 0) {
      if (V->consA[S22] < 0.0) {
	V->consA[S22] *= -1.0;
      }
    }
    else if (I[2] == G[2] - 1) {
      if (V->consA[S22] > 0.0) {
	V->consA[S22] *= -1.0;
      }
    }
  }
}

static int inside_cavity_cut(m2vol *V)
{
  return V->x1[1] < 1.0;
}
static int outside_cavity_cut(m2vol *V)
{
  return V->x1[1] >= 1.0;
}
static void write_log_entry(m2sim *m2, char *fname, int (*cut)(m2vol *V))
{
  FILE *outfile = fopen(fname, "a");
  if (outfile == NULL) {
    MSGF(FATAL, "open file '%s' failed", fname);
  }
  fprintf(outfile, "%+8.6e %+8.6e %+8.6e %+8.6e %+8.6e\n",
  	  m2->status.time_simulation,
  	  m2sim_volume_integral(m2, M2_COMOVING_MASS_DENSITY, cut),
  	  m2sim_volume_integral(m2, M2_KINETIC_ENERGY_DENSITY, cut),
  	  m2sim_volume_integral(m2, M2_INTERNAL_ENERGY_DENSITY, cut),
  	  m2sim_volume_integral(m2, M2_MAGNETIC_PRESSURE, cut));
  fclose(outfile);
}
static void write_equatorial_slice_1d(m2sim *m2, char *fname)
{
  FILE *outfile = fopen(fname, "w");
  m2vol *V;
  int *L = m2->local_grid_size;
  int *G = m2->domain_resolution;
  int i, j = G[2] / 2, k = 0;
  printf("[write equatorial slice 1d: %s]\n", fname);
  fprintf(outfile, "# r d p u1 u2 u3 b1 b2 b3 sigma mach\n");
  for (i=0; i<L[1]; ++i) {
    V = M2_VOL(i, j, k);
    fprintf(outfile,
	    "%+8.6e %+8.6e %+8.6e %+8.6e "
	    "%+8.6e %+8.6e %+8.6e %+8.6e "
	    "%+8.6e %+8.6e %+8.6e\n",
	    m2vol_coordinate_centroid(V, 1),
	    m2aux_get(&V->aux, M2_COMOVING_MASS_DENSITY),
	    m2aux_get(&V->aux, M2_GAS_PRESSURE),
	    m2aux_get(&V->aux, M2_VELOCITY_FOUR_VECTOR1),
	    m2aux_get(&V->aux, M2_VELOCITY_FOUR_VECTOR2),
	    m2aux_get(&V->aux, M2_VELOCITY_FOUR_VECTOR3),
	    m2aux_get(&V->aux, M2_MAGNETIC_FOUR_VECTOR1),
	    m2aux_get(&V->aux, M2_MAGNETIC_FOUR_VECTOR2),
	    m2aux_get(&V->aux, M2_MAGNETIC_FOUR_VECTOR3),
	    m2aux_get(&V->aux, M2_SIGMA),
	    m2aux_get(&V->aux, M2_MACH_NUMBER));
  }
  fclose(outfile);
}
void mwn_analysis(m2sim *m2)
{
  char slice_fname[1024];
  if (m2->status.iteration_number % 100 == 0) {
    write_log_entry(m2, FNAME_VOLUME_INTEGRALS_I, inside_cavity_cut);
    write_log_entry(m2, FNAME_VOLUME_INTEGRALS_O, outside_cavity_cut);
  }
  if (m2->status.iteration_number % 10000 == 0) {
    snprintf(slice_fname, sizeof(slice_fname),
	     "equatorial_slice.%06d.dat", m2->status.iteration_number);
    write_equatorial_slice_1d(m2, slice_fname);
  }
}
void mwn_initialize_problem(m2sim *m2)
{
  m2sim_set_resolution(m2, 96, 64, 1);
  m2sim_set_guard_zones(m2, 0);
  m2sim_set_geometry(m2, M2_SPHERICAL);
  m2sim_set_extent0(m2, 0.1, 0.0  , 0.0    );
  m2sim_set_extent1(m2, 1e1, M2_PI, 2*M2_PI);
  m2sim_set_physics(m2, M2_NONRELATIVISTIC | M2_MAGNETIZED);
  m2sim_set_ct_scheme(m2, M2_CT_OUTOFPAGE3);
  m2sim_set_analysis(m2, mwn_analysis);
  m2sim_set_boundary_conditions(m2, mwn_boundary_conditions);
  m2sim_set_initial_data(m2, mwn_initial_data);

  remove(FNAME_VOLUME_INTEGRALS_I);
  remove(FNAME_VOLUME_INTEGRALS_O);
}
