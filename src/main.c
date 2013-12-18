#include <stdio.h>
#include <stdlib.h>
#include "m2.h"


void initial_data(m2vol *V)
{
  double x[4];
  m2vol_coordinate_centroid_3d(V, x);
  if (x[1] < 0.5) {
    V->prim.v1 = 0.0;
    V->prim.v2 = 0.0;
    V->prim.v3 = 0.0;
    V->prim.d = 1.0;
    V->prim.p = 1.0;
  }
  else {
    V->prim.v1 = 0.0;
    V->prim.v2 = 0.0;
    V->prim.v3 = 0.0;
    V->prim.d = 0.125;
    V->prim.p = 0.1;
  }
}


void m2sim_runge_kutta_substep(m2sim *m2, double dt, double rkparam)
{
  m2sim_calculate_flux(m2);
  m2sim_exchange_flux(m2, dt);
  m2sim_add_source_terms(m2, dt);
  m2sim_average_runge_kutta(m2, rkparam);
  m2sim_enforce_boundary_condition(m2);
  m2sim_from_conserved_all(m2);
}


int main()
{
  m2sim *m2 = m2sim_new();


  m2sim_set_resolution(m2, 1024, 1, 1);
  m2sim_set_guard_zones(m2, 1);
  m2sim_set_geometry(m2, M2_CARTESIAN);
  //  m2sim_set_geometry(m2, M2_SPHERICAL);
  m2sim_set_physics(m2, M2_NONRELATIVISTIC | M2_UNMAGNETIZED);
  m2sim_set_extent0(m2, 0.1, 0.5*M2_PI-0.1, -0.1);
  m2sim_set_extent1(m2, 1.0, 0.5*M2_PI+0.1,  0.1);
  m2sim_initialize(m2);
  m2sim_map(m2, initial_data);
  m2sim_calculate_conserved(m2);


  printf("[m2]: astrophysical MHD code\n");
  printf("[m2]: memory usage %d MB]\n", m2sim_memory_usage(m2));

  double time_simulation = 0.0;
  double dt;
  int iteration_number = 0;
  int rk_order = 3;

  clock_t start_cycle = 0, stop_cycle = 0;
  double kzps; /* kilozones per second */


  while (time_simulation < 0.2) {

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
	   iteration_number,
	   time_simulation,
	   dt,
	   kzps);

    iteration_number += 1;
    time_simulation += dt;
  }

  FILE *outfile = fopen("m2.dat", "w");
  int i;
  for (i=0; i<m2->local_grid_size[1]; ++i) {
    fprintf(outfile, "%f %f %f %f\n",
	    m2vol_coordinate_centroid(&m2->volumes[M2_IND(i,0,0)], 1),
	    m2->volumes[M2_IND(i,0,0)].prim.d,
	    m2->volumes[M2_IND(i,0,0)].prim.p,
	    m2->volumes[M2_IND(i,0,0)].prim.v1);
  }

  m2sim_del(m2);


  return 0;
}
