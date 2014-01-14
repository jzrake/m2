#include <string.h>
#include <float.h>
#include "m2.h"


m2sim *m2sim_new()
{
  m2sim *m2 = (m2sim*) malloc(sizeof(m2sim));
  m2->domain_extent_lower[0] = 0.0;
  m2->domain_extent_lower[1] = 0.0;
  m2->domain_extent_lower[2] = 0.0;
  m2->domain_extent_lower[3] = 0.0;
  m2->domain_extent_upper[0] = 1.0;
  m2->domain_extent_upper[1] = 1.0;
  m2->domain_extent_upper[2] = 1.0;
  m2->domain_extent_upper[3] = 1.0;
  m2->domain_resolution[0] = 0;
  m2->domain_resolution[1] = 128;
  m2->domain_resolution[2] = 1;
  m2->domain_resolution[3] = 1;
  m2->local_grid_size[0] = 0;
  m2->local_grid_size[1] = 0;
  m2->local_grid_size[2] = 0;
  m2->local_grid_size[3] = 0;
  m2->number_guard_zones = 2;
  m2->geometry = M2_CARTESIAN;
  m2->physics = M2_NONRELATIVISTIC | M2_UNMAGNETIZED;
  m2->volumes = NULL;
  m2->rk_order = 2;
  m2->ct_scheme = M2_CT_FULL3D;

  /* callback functions */
  m2->analysis = NULL;
  m2->boundary_conditions = NULL;
  m2->boundary_conditions_gradient = NULL;
  m2->initial_data = NULL;

  /* status initializer */
  m2->status.time_simulation = 0.0;
  m2->status.iteration_number = 0;
  return m2;
}
void m2sim_del(m2sim *m2)
{
  free(m2->volumes);
  free(m2);
}
void m2sim_set_resolution(m2sim *m2, int n1, int n2, int n3)
{
  m2->domain_resolution[1] = n1;
  m2->domain_resolution[2] = n2;
  m2->domain_resolution[3] = n3;
}
void m2sim_set_guard_zones(m2sim *m2, int ng)
{
  m2->number_guard_zones = ng;
}
void m2sim_set_extent0(m2sim *m2, double x1, double x2, double x3)
{
  m2->domain_extent_lower[1] = x1;
  m2->domain_extent_lower[2] = x2;
  m2->domain_extent_lower[3] = x3;
}
void m2sim_set_extent1(m2sim *m2, double x1, double x2, double x3)
{
  m2->domain_extent_upper[1] = x1;
  m2->domain_extent_upper[2] = x2;
  m2->domain_extent_upper[3] = x3;
}
void m2sim_set_geometry(m2sim *m2, int geometry)
{
  m2->geometry = geometry;
}
void m2sim_set_physics(m2sim *m2, int modes)
{
  m2->physics = modes;
}
void m2sim_set_ct_scheme(m2sim *m2, int mode)
{
  m2->ct_scheme = mode;
}
void m2sim_set_rk_order(m2sim *m2, int order)
{
  m2->rk_order = order;
}
void m2sim_set_analysis(m2sim *m2, m2sim_operator analysis)
{
  m2->analysis = analysis;
}
void m2sim_set_boundary_conditions(m2sim *m2, m2sim_operator bc)
{
  m2->boundary_conditions = bc;
}
void m2sim_set_boundary_conditions_gradient(m2sim *m2, m2sim_operator bcg)
{
  m2->boundary_conditions_gradient = bcg;
}
void m2sim_set_initial_data(m2sim *m2, m2vol_operator id)
{
  m2->initial_data = id;
}
void m2sim_initialize(m2sim *m2)
{
  int *N = m2->domain_resolution;
  int *L = m2->local_grid_size;
  int ng = m2->number_guard_zones;
  int i, j, k, d;
  double I0[4], I1[4], x0[4], x1[4];
  m2vol *V;
  L[1] = N[1] + 2*ng*(N[1] > 1);
  L[2] = N[2] + 2*ng*(N[2] > 1);
  L[3] = N[3] + 2*ng*(N[3] > 1);
  L[0] = L[1] * L[2] * L[3];
  m2->volumes = (m2vol*) malloc(L[0] * sizeof(m2vol));
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {
	V = m2->volumes + M2_IND(i,j,k);
	V->m2 = m2;
	/* ----------------------- */
	/* cache cell global index */
	/* ----------------------- */
	for (d=1; d<=3; ++d) {
	  V->global_index[0] = M2_IND(i,j,k);
	  V->global_index[1] = i - ng * (N[1] > 1);
	  V->global_index[2] = j - ng * (N[2] > 1);
	  V->global_index[3] = k - ng * (N[3] > 1);
	}
	/* ----------------------- */
	/* cache cell volumes      */
	/* ----------------------- */
	I0[0] = 0.0;
	I0[1] = V->global_index[1] - 0.5;
	I0[2] = V->global_index[2] - 0.5;
	I0[3] = V->global_index[3] - 0.5;
	I1[0] = 0.0;
	I1[1] = V->global_index[1] + 0.5;
	I1[2] = V->global_index[2] + 0.5;
	I1[3] = V->global_index[3] + 0.5;
	m2sim_index_to_position(m2, I0, x0);
	m2sim_index_to_position(m2, I1, x1);
	V->volume = m2_volume_measure(x0, x1, m2->geometry);
	memcpy(V->x0, x0, 4 * sizeof(double));
	memcpy(V->x1, x1, 4 * sizeof(double));
	/* ------------------------ */
	/* cache cell surface areas */
	/* ------------------------ */
	I0[0] = 0.0;
	I0[1] = V->global_index[1] + 0.5;
	I0[2] = V->global_index[2] - 0.5;
	I0[3] = V->global_index[3] - 0.5;
	I1[0] = 0.0;
	I1[1] = V->global_index[1] + 0.5;
	I1[2] = V->global_index[2] + 0.5;
	I1[3] = V->global_index[3] + 0.5;
	m2sim_index_to_position(m2, I0, x0);
	m2sim_index_to_position(m2, I1, x1);
	V->area1 = m2_area_measure(x0, x1, m2->geometry, 1); /* axis 1 */
	I0[0] = 0.0;
	I0[1] = V->global_index[1] - 0.5;
	I0[2] = V->global_index[2] + 0.5;
	I0[3] = V->global_index[3] - 0.5;
	I1[0] = 0.0;
	I1[1] = V->global_index[1] + 0.5;
	I1[2] = V->global_index[2] + 0.5;
	I1[3] = V->global_index[3] + 0.5;
	m2sim_index_to_position(m2, I0, x0);
	m2sim_index_to_position(m2, I1, x1);
	V->area2 = m2_area_measure(x0, x1, m2->geometry, 2); /* axis 2 */
	I0[0] = 0.0;
	I0[1] = V->global_index[1] - 0.5;
	I0[2] = V->global_index[2] - 0.5;
	I0[3] = V->global_index[3] + 0.5;
	I1[0] = 0.0;
	I1[1] = V->global_index[1] + 0.5;
	I1[2] = V->global_index[2] + 0.5;
	I1[3] = V->global_index[3] + 0.5;
	m2sim_index_to_position(m2, I0, x0);
	m2sim_index_to_position(m2, I1, x1);
	V->area3 = m2_area_measure(x0, x1, m2->geometry, 3); /* axis 3 */
	/* ---------------------------- */
	/* cache cell linear dimensions */
	/* ---------------------------- */
	I0[0] = 0.0;
	I0[1] = V->global_index[1] - 0.5;
	I0[2] = V->global_index[2] + 0.5;
	I0[3] = V->global_index[3] + 0.5;
	I1[0] = 0.0;
	I1[1] = V->global_index[1] + 0.5;
	I1[2] = V->global_index[2] + 0.5;
	I1[3] = V->global_index[3] + 0.5;
	m2sim_index_to_position(m2, I0, x0);
	m2sim_index_to_position(m2, I1, x1);
	V->line1 = m2_line_measure(x0, x1, m2->geometry, 1); /* axis 1 */
	I0[0] = 0.0;
	I0[1] = V->global_index[1] + 0.5;
	I0[2] = V->global_index[2] - 0.5;
	I0[3] = V->global_index[3] + 0.5;
	I1[0] = 0.0;
	I1[1] = V->global_index[1] + 0.5;
	I1[2] = V->global_index[2] + 0.5;
	I1[3] = V->global_index[3] + 0.5;
	m2sim_index_to_position(m2, I0, x0);
	m2sim_index_to_position(m2, I1, x1);
	V->line2 = m2_line_measure(x0, x1, m2->geometry, 2); /* axis 2 */
	I0[0] = 0.0;
	I0[1] = V->global_index[1] + 0.5;
	I0[2] = V->global_index[2] + 0.5;
	I0[3] = V->global_index[3] - 0.5;
	I1[0] = 0.0;
	I1[1] = V->global_index[1] + 0.5;
	I1[2] = V->global_index[2] + 0.5;
	I1[3] = V->global_index[3] + 0.5;
	m2sim_index_to_position(m2, I0, x0);
	m2sim_index_to_position(m2, I1, x1);
	V->line3 = m2_line_measure(x0, x1, m2->geometry, 3); /* axis 3 */
      }
    }
  }
}
void m2sim_map(m2sim *m2, m2vol_operator f)
{
  int *L = m2->local_grid_size;
  int i, j, k;
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {
	f(m2->volumes + M2_IND(i,j,k));
      }
    }
  }
}
double m2vol_minimum_dimension(m2vol *V)
{
  double l1 = V->m2->domain_resolution[1] > 1 ? V->line1 : DBL_MAX;
  double l2 = V->m2->domain_resolution[2] > 1 ? V->line2 : DBL_MAX;
  double l3 = V->m2->domain_resolution[3] > 1 ? V->line3 : DBL_MAX;
  return M2_MIN3(l1, l2, l3);
}
double m2vol_coordinate_centroid(m2vol *V, int axis)
{
  return 0.5 * (V->x1[axis] + V->x0[axis]);
}
void m2vol_coordinate_centroid_3d(m2vol *V, double x[4])
{
  x[0] = 0.0;
  x[1] = 0.5 * (V->x1[1] + V->x0[1]);
  x[2] = 0.5 * (V->x1[2] + V->x0[2]);
  x[3] = 0.5 * (V->x1[3] + V->x0[3]);
}
int m2sim_memory_usage(m2sim *m2)
{
  return m2->local_grid_size[0] * sizeof(m2vol) / (1<<20);
}
double m2aux_get(m2aux *aux, int member)
{
  switch (member) {
  case M2_VELOCITY_FOUR_VECTOR0: return aux->velocity_four_vector[0];
  case M2_VELOCITY_FOUR_VECTOR1: return aux->velocity_four_vector[1];
  case M2_VELOCITY_FOUR_VECTOR2: return aux->velocity_four_vector[2];
  case M2_VELOCITY_FOUR_VECTOR3: return aux->velocity_four_vector[3];
  case M2_MAGNETIC_FOUR_VECTOR0: return aux->magnetic_four_vector[0];
  case M2_MAGNETIC_FOUR_VECTOR1: return aux->magnetic_four_vector[1];
  case M2_MAGNETIC_FOUR_VECTOR2: return aux->magnetic_four_vector[2];
  case M2_MAGNETIC_FOUR_VECTOR3: return aux->magnetic_four_vector[3];
  case M2_COMOVING_MASS_DENSITY: return aux->comoving_mass_density;
  case M2_GAS_PRESSURE: return aux->gas_pressure;
  case M2_MAGNETIC_PRESSURE: return aux->magnetic_pressure;
  default: return m2aux_measure(aux, member);
  }
}
double m2sim_volume_integral(m2sim *m2, int member, int (*cut_cb)(m2vol *V))
{
  m2vol *V;
  int *L = m2->local_grid_size;
  int *G = m2->domain_resolution;
  int n;
  double y = 0.0;
  for (n=0; n<L[0]; ++n) {
    V = m2->volumes + n;
    if (V->global_index[1] < 0 || (V->global_index[1] >= G[1] && G[1] > 1) ||
	V->global_index[2] < 0 || (V->global_index[2] >= G[2] && G[2] > 1) ||
	V->global_index[3] < 0 || (V->global_index[3] >= G[3] && G[3] > 1)) {
      continue;
    }
    else if (cut_cb) {
      if (cut_cb(V) == 0) {
	continue;
      }
    }
    y += m2aux_get(&V->aux, member) * V->volume;
  }
  return y;
}
void m2sim_index_global_to_local(m2sim *m2, int global_index[4], int I[4])
{
  int *G = m2->domain_resolution;
  int i = global_index[1];
  int j = global_index[2];
  int k = global_index[3];
  int ng = m2->number_guard_zones;
  I[0] = global_index[0];
  I[1] = i + ng * (G[1] > 1);
  I[2] = j + ng * (G[2] > 1);
  I[3] = k + ng * (G[3] > 1);
}
void m2sim_run_analysis(m2sim *m2)
{
  if (m2->analysis) {
    m2->analysis(m2);
  }
  else {
    /* it's OK, analysis function isn't required */
  }
}
void m2sim_run_initial_data(m2sim *m2)
{
  if (m2->initial_data) {
    m2sim_map(m2, m2->initial_data);
  }
  else {
    MSG(FATAL, "no initial data function supplied");
  }
}
