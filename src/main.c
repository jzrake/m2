#include <stdio.h>
#include <stdlib.h>
#include "m2.h"
#include "problems.h"

void m2sim_runge_kutta_substep(m2sim *m2, double dt, double rkparam)
{
  m2sim_calculate_gradient(m2);
  m2sim_calculate_flux(m2);
  m2sim_calculate_emf(m2);
  m2sim_exchange_flux(m2, dt);
  m2sim_add_source_terms(m2, dt);
  m2sim_average_runge_kutta(m2, rkparam);
  m2sim_enforce_boundary_condition(m2);
  m2sim_magnetic_flux_to_cell_center(m2);
  m2sim_from_conserved_all(m2);
}


void m2sim_drive(m2sim *m2)
{
  double dt;
  int rk_order = m2->rk_order;

  clock_t start_cycle = 0, stop_cycle = 0;
  double kzps; /* kilozones per second */

  m2sim_run_analysis(m2);
  m2sim_cache_conserved(m2);

  dt = 0.5 * m2sim_minimum_courant_time(m2);

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


int main(int argc, char **argv)
{
  m2sim *m2 = m2sim_new();


  if (0) {
    m2sim_set_geometry(m2, M2_CYLINDRICAL);
    m2sim_set_extent0(m2, 0.1, 0.0    , -0.5);
    m2sim_set_extent1(m2, 1.0, 2*M2_PI, +0.5);
    m2sim_set_physics(m2, M2_NONRELATIVISTIC | M2_MAGNETIZED);
  }
  else if (1) {
    mwn_initialize_problem(m2);
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
  m2sim_run_initial_data(m2);
  m2sim_magnetic_flux_to_cell_center(m2);
  m2sim_from_primitive_all(m2);

  printf("[m2]: astrophysical MHD code\n");
  printf("[m2]: memory usage %d MB\n", m2sim_memory_usage(m2));
  m2sim_print(m2);

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


