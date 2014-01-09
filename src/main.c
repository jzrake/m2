#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "m2.h"


void m2sim_write_equatorial_slice_1d(m2sim *m2, char *fname);
void m2sim_write_log_entry(m2sim *m2, char *fname, int (*cut)(m2vol *V));

#define FNAME_VOLUME_INTEGRALS_I "volume-integrals-in.dat"
#define FNAME_VOLUME_INTEGRALS_O "volume-integrals-out.dat"

static int inside_cavity_cut(m2vol *V)
{
  return V->x1[1] < 1.0;
}
static int outside_cavity_cut(m2vol *V)
{
  return V->x1[1] >= 1.0;
}

void initial_data(m2vol *V)
{
  double R = m2vol_coordinate_centroid(V, 1);

  /* double x = m2vol_coordinate_centroid(V, 1); */
  /* double y = m2vol_coordinate_centroid(V, 2); */
  /* double R = sqrt(x*x + y*y); */

  /* if (R < 0.125) { */
  /*   V->prim.v1 = 0.0; */
  /*   V->prim.v2 = 0.0; */
  /*   V->prim.v3 = 0.0; */
  /*   V->prim.d  = 1.0; */
  /*   V->prim.p  = 100.0; */
  /* } */
  /* else { */
  /*   V->prim.v1 = 0.0; */
  /*   V->prim.v2 = 0.0; */
  /*   V->prim.v3 = 0.0; */
  /*   V->prim.d  = 1.0; */
  /*   V->prim.p  = 1.0; */
  /* } */

  if (R < 1.0) {
    V->prim.v1 = 0.0;
    V->prim.v2 = 0.0;
    V->prim.v3 = 0.0;
    V->prim.d  = 1.0;
    V->prim.p  = 0.01;
  }
  else {
    V->prim.v1 = 0.0;
    V->prim.v2 = 0.0;
    V->prim.v3 = 0.0;
    V->prim.d  = 1000.0;
    V->prim.p  = 0.01;
  }

  V->Bflux1A = 0.0 * V->area1;
  V->Bflux2A = 0.0 * V->area2;
  V->Bflux3A = 0.0 * V->area3;
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
  m2sim_from_conserved_all(m2);
  m2sim_magnetic_flux_to_cell_center(m2);
}


void m2sim_drive(m2sim *m2)
{
  double dt;
  int rk_order = 2;

  clock_t start_cycle = 0, stop_cycle = 0;
  double kzps; /* kilozones per second */


  if (m2->status.iteration_number % 100 == 0) {
    m2sim_write_log_entry(m2, FNAME_VOLUME_INTEGRALS_I, inside_cavity_cut);
    m2sim_write_log_entry(m2, FNAME_VOLUME_INTEGRALS_O, outside_cavity_cut);
  }

  dt = 0.5 * m2sim_minimum_courant_time(m2);

  m2sim_cache_conserved(m2);

  start_cycle = clock();
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
  stop_cycle = clock();
  kzps = 1e-3 * m2->local_grid_size[0] /
    (stop_cycle - start_cycle) * CLOCKS_PER_SEC;

  /* print status message */
  printf("%08d: t=%5.4f dt=%5.4e kzps=%4.3f\n",
	 m2->status.iteration_number,
	 m2->status.time_simulation,
	 dt,
	 kzps);

  m2->status.iteration_number += 1;
  m2->status.time_simulation += dt;
}


void m2sim_write_equatorial_slice_1d(m2sim *m2, char *fname)
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


void m2sim_write_log_entry(m2sim *m2, char *fname, int (*cut)(m2vol *V))
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


void m2sim_write_volume_integrals(m2sim *m2, char *fname)
{
  FILE *outfile = fopen(fname, "a");
  m2vol *V;
  int *L = m2->local_grid_size;
  int *G = m2->domain_resolution;
  int i, j = G[2] / 2, k = 0;
  
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


int main(int argc, char **argv)
{
  m2sim *m2 = m2sim_new();


  m2sim_set_resolution(m2, 48, 64, 1);
  m2sim_set_guard_zones(m2, 0);


  if (0) {
    m2sim_set_geometry(m2, M2_CYLINDRICAL);
    m2sim_set_extent0(m2, 0.1, 0.0    , -0.5);
    m2sim_set_extent1(m2, 1.0, 2*M2_PI, +0.5);
    m2sim_set_physics(m2, M2_NONRELATIVISTIC | M2_MAGNETIZED);
  }
  else if (1) {
    m2sim_set_geometry(m2, M2_SPHERICAL);
    m2sim_set_extent0(m2, 0.1, 0    , 0.0    );
    m2sim_set_extent1(m2, 1e1, M2_PI, 2*M2_PI);
    m2sim_set_physics(m2, M2_NONRELATIVISTIC | M2_MAGNETIZED);
  }
  else if (0) {
    m2sim_set_geometry(m2, M2_CARTESIAN);
    m2sim_set_extent0(m2, -0.5, -0.5, 0.0);
    m2sim_set_extent1(m2, +0.5, +0.5, 1.0);
    m2sim_set_physics(m2, M2_NONRELATIVISTIC | M2_MAGNETIZED);
  }
  else if (0) {
    m2sim_set_geometry(m2, M2_CARTESIAN);
    m2sim_set_extent0(m2, 0.0, 0.0, 0.0);
    m2sim_set_extent1(m2, 1.0, 1.0, 1.0);
    m2sim_set_physics(m2, M2_NONRELATIVISTIC | M2_MAGNETIZED);
  }


  m2sim_initialize(m2);
  m2sim_map(m2, initial_data);
  m2sim_magnetic_flux_to_cell_center(m2);
  m2sim_calculate_conserved(m2);


  printf("[m2]: astrophysical MHD code\n");
  printf("[m2]: memory usage %d MB]\n", m2sim_memory_usage(m2));
  m2sim_print(m2);

  remove(FNAME_VOLUME_INTEGRALS_I);
  remove(FNAME_VOLUME_INTEGRALS_O);

  if (1) {
    m2sim_visualize(m2, argc, argv);
  }
  else {
    while (m2->status.time_simulation < 1000.0) {
      m2sim_drive(m2);
    }
  }
  //  m2sim_visualize(m2, argc, argv);
  //  m2sim_write_ascii_2d(m2, "m2.dat");
  m2sim_write_ascii_1d(m2, "m2.dat");
  m2sim_del(m2);


  return 0;
}


