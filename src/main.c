#include <stdio.h>
#include <stdlib.h>
#include "m2.h"


void initial_data(m2vol *V)
{
  double x[4];
  double index[4] = { V->global_index[0],
		      V->global_index[1],
		      V->global_index[2],
		      V->global_index[3] };
  m2sim_index_to_position(V->m2, index, x);
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


int main()
{
  m2sim *m2 = m2sim_new();


  m2sim_set_resolution(m2, 400, 1, 1);
  m2sim_set_guard_zones(m2, 2);
  m2sim_set_geometry(m2, M2_CARTESIAN);
  m2sim_set_physics(m2, M2_RELATIVISTIC | M2_UNMAGNETIZED);
  m2sim_set_extent0(m2, 0.0, 0.0, 0.0);
  m2sim_set_extent1(m2, 1.0, 1.0, 1.0);
  m2sim_initialize(m2);
  m2sim_map(m2, initial_data);
  m2sim_calculate_conserved(m2);


  double time_simulation = 0.0;
  double dt;


  while (time_simulation < 0.1) {

    dt = 0.5 * m2sim_minimum_courant_time(m2);

    m2sim_cache_conserved(m2);
    m2sim_calculate_flux(m2);
    m2sim_exchange_flux(m2, dt);
    m2sim_from_conserved_all(m2);
    m2sim_enforce_boundary_condition(m2);
    time_simulation += dt;
  }


  int i;
  for (i=0; i<m2->local_grid_size[1]; ++i) {
    printf("%d %f %f %f\n", i,
  	   m2->volumes[i].prim.d,
  	   m2->volumes[i].prim.p,
  	   m2->volumes[i].prim.v1);
  }

  m2sim_del(m2);


  return 0;
}
