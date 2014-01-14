#include <stdio.h>
#include <stdlib.h>
#include "m2.h"


void initialize_proble_mmwn(m2sim *m2);
void initialize_problem_shocktube(m2sim *m2);


int main(int argc, char **argv)
{
  if (0) {
    m2_self_test();
    return 0;
  }

  m2sim *m2 = m2sim_new();

  //initialize_problem_mwn(m2);
  initialize_problem_shocktube(m2);

  m2sim_initialize(m2);
  m2sim_run_initial_data(m2);
  m2sim_magnetic_flux_to_cell_center(m2);
  m2sim_from_primitive_all(m2);

  printf("[m2]: astrophysical MHD code\n");
  printf("[m2]: memory usage %d MB\n", m2sim_memory_usage(m2));
  m2sim_print(m2);

  if (0) {
    m2sim_visualize(m2, argc, argv);
  }
  else {
    while (m2->status.time_simulation < 0.2) {
      m2sim_drive(m2);
    }
  }
  //  m2sim_visualize(m2, argc, argv);
  //  m2sim_write_ascii_2d(m2, "m2.dat");
  m2sim_write_ascii_1d(m2, "m2.dat");
  m2sim_del(m2);

  return 0;
}
