#include <string.h>
#include <float.h>
#include <math.h>
#include "m2.h"


void m2sim_new(m2sim *m2)
{
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
  m2->local_grid_start[0] = 0;
  m2->local_grid_start[1] = 0;
  m2->local_grid_start[2] = 0;
  m2->local_grid_start[3] = 0;
  m2->number_guard_zones0[0] = 0;
  m2->number_guard_zones0[1] = 0;
  m2->number_guard_zones0[2] = 0;
  m2->number_guard_zones0[3] = 0;
  m2->number_guard_zones1[0] = 0;
  m2->number_guard_zones1[1] = 0;
  m2->number_guard_zones1[2] = 0;
  m2->number_guard_zones1[3] = 0;
  m2->periodic_dimension[0] = 0;
  m2->periodic_dimension[1] = 0;
  m2->periodic_dimension[2] = 0;
  m2->periodic_dimension[3] = 0;

  m2->volumes = NULL;
  m2->user_struct = NULL;
  m2->cart_comm = NULL;
  m2->mpi_types = NULL;

  /* solver and physics configuration */
  m2->coordinate_scaling1 = M2_LINEAR;
  m2->coordinate_scaling2 = M2_LINEAR;
  m2->coordinate_scaling3 = M2_LINEAR;
  m2->geometry = M2_CARTESIAN; 
  m2->relativistic = 0;
  m2->magnetized = 0;
  m2->rk_order = 2;
  m2->gradient_fields = M2_PRIMITIVE;
  m2->gradient_estimation_method = M2_GRADIENT_ESTIMATION_PLM;
  m2->riemann_solver = M2_RIEMANN_HLLE;
  m2->quartic_solver = M2_QUARTIC_FULL_COMPLEX;
  m2->plm_parameter = 1.5;
  m2->cfl_parameter = 0.4;
  m2->gamma_law_index = 4./3;

  /* callback functions */
  m2->analysis = NULL;
  m2->boundary_conditions_cell = NULL;
  m2->boundary_conditions_flux1 = NULL;
  m2->boundary_conditions_flux2 = NULL;
  m2->boundary_conditions_flux3 = NULL;
  m2->boundary_conditions_emf1 = NULL;
  m2->boundary_conditions_emf2 = NULL;
  m2->boundary_conditions_emf3 = NULL;

  m2->boundary_conditions_gradient = NULL;
  m2->add_physical_source_terms = NULL;
  m2->initial_data = NULL;
  m2->initial_data_cell = NULL;
  m2->initial_data_face = NULL;
  m2->initial_data_edge = NULL;

  m2->cadence_checkpoint_hdf5 = 0.0;
  m2->cadence_checkpoint_tpl = 0.0;
  m2->cadence_analysis = 0.0;

  /* status initializer */
  m2->status.iteration_number = 0;
  m2->status.time_simulation = 0.0;
  m2->status.time_last_checkpoint_hdf5 = 0.0;
  m2->status.time_last_checkpoint_tpl = 0.0;
  m2->status.time_last_analysis = 0.0;
  m2->status.checkpoint_number_hdf5 = 0;
  m2->status.checkpoint_number_tpl = 0;
  m2->status.drive_attempt = 0;
  m2->suppress_extrapolation_at_unhealthy_zones = 1;
  m2->pressure_floor = -1.0; /* disabled */

  mesh_new(&m2->mesh);
  jsw_seed(&m2->stirring.random, 12345);
}



void m2sim_del(m2sim *m2)
{
  if (m2->cart_comm) {
    m2sim_delete_mpi_types(m2);
  }
  free(m2->volumes);
  m2->volumes = NULL;
  mesh_deallocate(&m2->mesh);
  m2stirring_del(&m2->stirring);
}



void m2sim_initialize(m2sim *m2)
/*
 * Each axes is either periodic or is bounded by hard walls at both sides. If it
 * is periodic, then the number of guard zones in both directions
 *
 */
{
  int *L = m2->local_grid_size;
  int *S = m2->local_grid_start;
  int n, d;
  double I0[4], I1[4], x0[4], x1[4];

  L[0] = L[1] * L[2] * L[3];

  if (m2->cart_comm) {
    m2sim_create_mpi_types(m2);
  }


  m2->mesh.cells_shape[1] = L[1];
  m2->mesh.cells_shape[2] = L[2];
  m2->mesh.cells_shape[3] = L[3];
  mesh_allocate(&m2->mesh);

  printf("F: %d %d\n", m2->mesh.faces_shape[1][0], m2->mesh.faces_shape[1][1]);
  printf("C: %d %d\n", m2->mesh.cells_shape   [0], m2->mesh.cells_shape   [1]);

  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    int I[4];
    struct mesh_cell *C = m2->mesh.cells + n;
    mesh_linear_to_multi(n, m2->mesh.cells_shape, I);
    C->global_i = I[1] + S[1];
    C->global_j = I[2] + S[2];
    C->global_k = I[3] + S[3];
    I0[1] = C->global_i - 0.5;
    I0[2] = C->global_j - 0.5;
    I0[3] = C->global_k - 0.5;
    I1[1] = C->global_i + 0.5;
    I1[2] = C->global_j + 0.5;
    I1[3] = C->global_k + 0.5;
    m2sim_index_to_position(m2, I0, x0);
    m2sim_index_to_position(m2, I1, x1);
    C->volume = m2_volume_measure(x0, x1, m2->geometry);
    C->x[1] = 0.5 * (x0[1] + x1[1]);
    C->x[2] = 0.5 * (x0[2] + x1[2]);
    C->x[3] = 0.5 * (x0[3] + x1[3]);
  }

  for (d=1; d<=3; ++d) {
    for (n=0; n<m2->mesh.faces_shape[d][0]; ++n) {
      int I[4];
      struct mesh_face *F = m2->mesh.faces[d] + n;
      mesh_linear_to_multi(n, m2->mesh.faces_shape[d], I);
      I0[1] = I[1] + S[1];
      I0[2] = I[2] + S[2];
      I0[3] = I[3] + S[3];
      I1[1] = I[1] + S[1] + 1.0 * (d != 1);
      I1[2] = I[2] + S[2] + 1.0 * (d != 2);
      I1[3] = I[3] + S[3] + 1.0 * (d != 3);
      m2sim_index_to_position(m2, I0, x0);
      m2sim_index_to_position(m2, I1, x1);
      F->area = m2_area_measure(x0, x1, m2->geometry, d);
      F->x[1] = 0.5 * (x0[1] + x1[1]);
      F->x[2] = 0.5 * (x0[2] + x1[2]);
      F->x[3] = 0.5 * (x0[3] + x1[3]);
    }

    for (n=0; n<m2->mesh.edges_shape[d][0]; ++n) {
      int I[4];
      struct mesh_edge *E = m2->mesh.edges[d] + n;
      mesh_linear_to_multi(n, m2->mesh.edges_shape[d], I);
      I0[1] = I[1] + S[1];
      I0[2] = I[2] + S[2];
      I0[3] = I[3] + S[3];
      I1[1] = I[1] + S[1] + 1.0 * (d == 1);
      I1[2] = I[2] + S[2] + 1.0 * (d == 2);
      I1[3] = I[3] + S[3] + 1.0 * (d == 3);
      m2sim_index_to_position(m2, I0, x0);
      m2sim_index_to_position(m2, I1, x1);
      E->length = m2_line_measure(x0, x1, m2->geometry, d);
      E->x[1] = 0.5 * (x0[1] + x1[1]);
      E->x[2] = 0.5 * (x0[2] + x1[2]);
      E->x[3] = 0.5 * (x0[3] + x1[3]);
    }
  }
}



void m2sim_from_interpolated(m2sim *m2, double *y, m2prim *P)
{
  double u0;
  switch (m2->gradient_fields) {
  case M2_PRIMITIVE:
    P->v1 = y[0];
    P->v2 = y[1];
    P->v3 = y[2];
    P->B1 = y[3];
    P->B2 = y[4];
    P->B3 = y[5];
    P->d  = y[6];
    P->p  = y[7];
    break;
  case M2_PRIMITIVE_AND_FOUR_VELOCITY:
    if (m2->relativistic) {
      /* gamma * beta = {y[0], y[1], y[2]} */
      u0 = sqrt(1.0 + y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
    }
    else {
      u0 = 1.0;
    }
    P->v1 = y[0] / u0;
    P->v2 = y[1] / u0;
    P->v3 = y[2] / u0;
    P->B1 = y[3];
    P->B2 = y[4];
    P->B3 = y[5];
    P->d  = y[6];
    P->p  = y[7];
    break;
  default:
    MSG(FATAL, "unknown interpolation fields");
    break;
  }
}



void m2sim_to_interpolated(m2sim *m2, m2prim *P, m2aux *A, double *y, int stride)
{
  switch (m2->gradient_fields) {
  case M2_PRIMITIVE:
    y[0*stride] = P->v1;
    y[1*stride] = P->v2;
    y[2*stride] = P->v3;
    y[3*stride] = P->B1;
    y[4*stride] = P->B2;
    y[5*stride] = P->B3;
    y[6*stride] = P->d;
    y[7*stride] = P->p;
    break;
  case M2_PRIMITIVE_AND_FOUR_VELOCITY:
    /* gamma * beta = {y[0], y[1], y[2]} */
    y[0*stride] = A->velocity_four_vector[1];
    y[1*stride] = A->velocity_four_vector[2];
    y[2*stride] = A->velocity_four_vector[3];
    y[3*stride] = P->B1;
    y[4*stride] = P->B2;
    y[5*stride] = P->B3;
    y[6*stride] = P->d;
    y[7*stride] = P->p;
    break;
  default:
    MSG(FATAL, "unknown interpolation fields");
    break;
  }
}



double m2vol_minimum_dimension(m2vol *V)
{
  double l1 = V->m2->domain_resolution[1] > 1 ? V->line1 : DBL_MAX;
  double l2 = V->m2->domain_resolution[2] > 1 ? V->line2 : DBL_MAX;
  double l3 = V->m2->domain_resolution[3] > 1 ? V->line3 : DBL_MAX;
  if (l1 < 1e-12) l1 = DBL_MAX; /* don't worry about exactly zero cell size */
  if (l2 < 1e-12) l2 = DBL_MAX;
  if (l3 < 1e-12) l3 = DBL_MAX;
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
  const double u0 = aux->velocity_four_vector[0];
  const double u1 = aux->velocity_four_vector[1];
  const double u2 = aux->velocity_four_vector[2];
  const double u3 = aux->velocity_four_vector[3];
  const double b0 = aux->magnetic_four_vector[0];
  const double b1 = aux->magnetic_four_vector[1];
  const double b2 = aux->magnetic_four_vector[2];
  const double b3 = aux->magnetic_four_vector[3];
  switch (member) {
  case M2_MAGNETIC1: return b1 * u0 - b0 * u1;
  case M2_MAGNETIC2: return b2 * u0 - b0 * u2;
  case M2_MAGNETIC3: return b3 * u0 - b0 * u3;
  case M2_VELOCITY1: return u1 / u0;
  case M2_VELOCITY2: return u2 / u0;
  case M2_VELOCITY3: return u3 / u0;
  case M2_VELOCITY_FOUR_VECTOR0: return u0;
  case M2_VELOCITY_FOUR_VECTOR1: return u1;
  case M2_VELOCITY_FOUR_VECTOR2: return u2;
  case M2_VELOCITY_FOUR_VECTOR3: return u3;
  case M2_MAGNETIC_FOUR_VECTOR0: return b0;
  case M2_MAGNETIC_FOUR_VECTOR1: return b1;
  case M2_MAGNETIC_FOUR_VECTOR2: return b2;
  case M2_MAGNETIC_FOUR_VECTOR3: return b3;
  case M2_COMOVING_MASS_DENSITY: return aux->comoving_mass_density;
  case M2_GAS_PRESSURE: return aux->gas_pressure;
  case M2_MAGNETIC_PRESSURE: return aux->magnetic_pressure;
  case M2_PLASMA_BETA: return aux->gas_pressure / aux->magnetic_pressure;
  case M2_ENTROPY: return aux->gas_pressure /
      pow(aux->comoving_mass_density, aux->m2->gamma_law_index);
  default: return m2aux_measure(aux, member);
  }
}



double m2aux_mach(m2aux *aux, double n[4], int mode)
/*
 * Return a signed, directional Mach number for either fast, slow, or Alfven
 * modes. The Mach number is based on the separation and inclination of the left
 * and right going characteristic four-velocities.
 */
{
  double ap = 0.0, am = 0.0;
  double up = 0.0, um = 0.0;
  double evals[8];
  m2aux_eigenvalues(aux, n, evals);
  switch (mode) {
  case M2_MACH_ALFVEN:
    am = evals[1];
    ap = evals[6];
    break;
  case M2_MACH_FAST:
    am = evals[0];
    ap = evals[7];
    break;
  case M2_MACH_SLOW:
    am = evals[2];
    ap = evals[5];
    break;
  }
  if (aux->m2->relativistic) {
    up = ap / sqrt(1.0 - ap*ap);
    um = am / sqrt(1.0 - am*am);
  }
  else {
    up = ap;
    um = am;
  }
  return (up + um) / (up - um);
}



void m2sim_index_global_to_local(m2sim *m2, int global_index[4], int I[4])
{
  I[0] = global_index[0];
  I[1] = global_index[1] - m2->local_grid_start[1];
  I[2] = global_index[2] - m2->local_grid_start[2];
  I[3] = global_index[3] - m2->local_grid_start[3];
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
  int n, d;
  double N[4] = {0,0,0,0}; /* unit vector for faces and edges */
  double prim[5], Bfield, Efield;

  struct mesh_cell *C;
  struct mesh_face *F;
  struct mesh_edge *E;

  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    C = m2->mesh.cells + n;
    m2->initial_data_cell(m2, C->x, NULL, prim);
    C->prim.d  = prim[RHO];
    C->prim.p  = prim[PRE];
    C->prim.v1 = prim[V11];
    C->prim.v2 = prim[V22];
    C->prim.v3 = prim[V33];
    C->zone_type = M2_ZONE_TYPE_FULL; /* TODO: account for guard zones */
  }

  for (d=1; d<=3; ++d) {
    N[d] = 1.0;
    for (n=0; n<m2->mesh.faces_shape[d][0]; ++n) {
      F = m2->mesh.faces[d] + n;
      if (m2->initial_data_face) {
	m2->initial_data_face(m2, F->x, N, &Bfield);
      }
      else {
	Bfield = 0.0;
      }
      F->BfluxA = Bfield * F->area;
    }
    for (n=0; n<m2->mesh.edges_shape[d][0]; ++n) {
      E = m2->mesh.edges[d] + n;
      if (m2->initial_data_edge) {
	m2->initial_data_edge(m2, E->x, N, &Efield);
      }
      else {
	Efield = 0.0;
      }
      E->emf = Efield * E->length;
    }
    N[d] = 0.0;
  }

  /* in case data was set manually on interior zones, synch it first */
  m2sim_synchronize_guard(m2);
  /* append to face-centered B data with curl(A), where A is the emf data on the
     edges */
  m2sim_exchange_flux(m2, 1.0);
  m2sim_magnetic_flux_to_cell_center(m2);
  m2sim_from_primitive_all(m2);
  m2sim_synchronize_guard(m2);
}

m2vol *m2vol_neighbor(m2vol *V, int axis, int dist)
{
  /* TODO: emove this function or change to a cell function */
  m2sim *m2 = V->m2;
  int *L = m2->local_grid_size;
  int I[4];
  m2sim_index_global_to_local(m2, V->global_index, I);
  I[axis] += dist;
  return M2_VOL(I[1], I[2], I[3]);
}




/* -----------------------------------------------------------------------------
* These functions return 1 if the error was fatal, and 0 otherwise. For each of
* the error flags, there is a corresponding response description that includes
* whether that error should be discarded or pushed, or is fatal.
* -----------------------------------------------------------------------------
*/
int m2sim_error_cell(m2sim *m2, struct mesh_cell *C, int error_flag)
{
  return 0;
}
int m2sim_error_face(m2sim *m2, struct mesh_face *F, int error_flag)
{
  return 0;
}
int m2sim_error_edge(m2sim *m2, struct mesh_edge *E, int error_flag)
{
  return 0;
}
int m2sim_error_general(m2sim *m2, const char *message)
{
  return 1;
}
