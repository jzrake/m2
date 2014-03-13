#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "m2.h"

void initialize_problem_shocktube(m2sim *m2);
void initialize_problem_mwn(m2sim *m2);
void initialize_problem_jet(m2sim *m2);
void initialize_problem_flux_burial(m2sim *m2);
void initialize_problem_spheromak(m2sim *m2);

int main(int argc, char **argv)
{
  int opt_vis = 0;
  int opt_self_test = 0;
  int n;
  char *restart_file = NULL;
  for (n=1; n<argc; ++n) {
    if (argv[n][0] == '-' && argv[n][1] == '-') {
      if (strcmp(argv[n], "--self-test") == 0) opt_self_test = 1;
      else if (strcmp(argv[n], "--vis") == 0) opt_vis = 1;
      else if (strncmp(argv[n], "--restart=", 10) == 0) {
	restart_file = &argv[n][10];
      }
      else {
	printf("[m2]: unrecognized option '%s'\n", argv[n]);
	return 0;
      }
    }
  }

  if (opt_self_test) {
    m2_self_test();
    return 0;
  }

  m2sim *m2 = m2sim_new();

  //initialize_problem_spheromak(m2);
  //initialize_problem_flux_burial(m2);
  initialize_problem_jet(m2);

  m2sim_initialize(m2);
  if (restart_file) {
    if (m2sim_load_checkpoint(m2, restart_file)) {
      MSG(FATAL, "checkpoint load failed");
    }
  }
  else {
    m2sim_run_initial_data(m2);
  }

  printf("[m2]: astrophysical MHD code\n");
  printf("[m2]: memory usage %d MB\n", m2sim_memory_usage(m2));
  m2sim_print(m2);

  if (opt_vis) {
    m2sim_visualize(m2, argc, argv);
  }
  else {
    while (m2->status.time_simulation < 500.0) {
      m2sim_drive(m2);
    }
  }
  m2sim_write_ascii_1d(m2, "m2.dat");
  m2sim_del(m2);

  return 0;
}
