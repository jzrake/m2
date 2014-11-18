#include <string.h>
#include <math.h>
#include "m2.h"
#if M2_HAVE_MPI
#include <mpi.h>
#endif



static inline void global_sum(m2sim *m2, int *N);
static inline void global_min(m2sim *m2, double *x);
static inline double plm_gradient(double *xs, double *fs, double plm, int *error);
static inline int global_index_cell(m2sim *m2, struct mesh_cell *C, int axis);
/* static inline int global_index_face(m2sim *m2, struct mesh_face *F, int d); */
static inline int global_index_edge(m2sim *m2, struct mesh_edge *E, int d);



int m2sim_calculate_flux(m2sim *m2)
{
  struct mesh_cell *cells[2];
  struct mesh_face *F, *G;
  int n, d;

  /* Solve the riemann problem between all adjacent cells. Faces that have only
     one cell attached are given the flux of the interior cell. */
  for (d=1; d<=3; ++d) {
    for (n=0; n<m2->mesh.faces_shape[d][0]; ++n) {
      m2sim_riemann_solver(m2, &m2->mesh.faces[d][n]);
    }
  }

  /* Here we apply reflecting boundary conditions on the fluxes, if
     necessary. */
  for (d=1; d<=3; ++d) {

    int B1, B2, B3;
    int S1, S2, S3;

    if (m2->domain_resolution[d] == 1) {
      continue;
    }
    else {
      switch (d) {
      case 1: B1 = B11; B2 = B22; B3 = B33; S1 = S11; S2 = S22; S3 = S33; break;
      case 2: B1 = B22; B2 = B33; B3 = B11; S1 = S22; S2 = S33; S3 = S11; break;
      case 3: B1 = B33; B2 = B11; B3 = B22; S1 = S33; S2 = S11; S3 = S22; break;
      }
    }

    for (n=0; n<m2->mesh.faces_shape[d][0]; ++n) {
      G = NULL;
      F = &m2->mesh.faces[d][n];
      mesh_face_cells(&m2->mesh, F, cells);
      if (cells[0] == NULL &&
	  m2->boundary_condition_lower[d] == M2_BOUNDARY_REFLECTING) {
	/* TODO: when shell thing is removed, this should be a zero not a -1 */
	if (global_index_cell(m2, cells[1], d) == -1) {
	  G = F + m2->mesh.faces_stride[d][d];
	}
      }
      if (cells[1] == NULL &&
	  m2->boundary_condition_upper[d] == M2_BOUNDARY_REFLECTING) {
	if (global_index_cell(m2, cells[0], d) == m2->domain_resolution[d]-1) {
	  G = F - m2->mesh.faces_stride[d][d];
	}
      }
      if (G) {
	F->flux[DDD] = -G->flux[DDD];
	F->flux[TAU] = -G->flux[TAU];
	F->flux[S1] = +G->flux[S1];
	F->flux[S2] = -G->flux[S2];
	F->flux[S3] = -G->flux[S3];
	F->flux[B1] = -G->flux[B1];
	F->flux[B2] = +G->flux[B2];
	F->flux[B3] = +G->flux[B3];
      }
    }
  }
  return 0;
}



int m2sim_calculate_emf(m2sim *m2)
{
  if (!m2->magnetized) {
    return 0;
  }
  struct mesh_edge *E;
  struct mesh_face *F[4];
  int n, d, num_faces;
  int B1, B2, B3;
  int a2, a3, g2, g3;
  for (d=1; d<=3; ++d) {
    switch (d) {
    case 1: B1 = B11; B2 = B22; B3 = B33; break;
    case 2: B1 = B22; B2 = B33; B3 = B11; break;
    case 3: B1 = B33; B2 = B11; B3 = B22; break;
    }
    for (n=0; n<m2->mesh.edges_shape[d][0]; ++n) {
      E = &m2->mesh.edges[d][n];
      num_faces = mesh_edge_faces(&m2->mesh, E, F);
      E->emf = 0.0;
      if (F[0]) E->emf -= F[0]->flux[B3];
      if (F[1]) E->emf -= F[1]->flux[B3];
      if (F[2]) E->emf += F[2]->flux[B2];
      if (F[3]) E->emf += F[3]->flux[B2];
      E->emf *= E->length / num_faces;
      /* reflecting boundary condition */
      a2 = (d+1-1)%3+1;
      a3 = (d+2-1)%3+1;
      if ((num_faces != 4) &&
	  (m2->boundary_condition_lower[a2] == M2_BOUNDARY_REFLECTING ||
	   m2->boundary_condition_lower[a3] == M2_BOUNDARY_REFLECTING ||
	   m2->boundary_condition_upper[a2] == M2_BOUNDARY_REFLECTING ||
	   m2->boundary_condition_upper[a3] == M2_BOUNDARY_REFLECTING)) {
	g2 = global_index_edge(m2, E, a2);
	g3 = global_index_edge(m2, E, a3);
	/* If the edge lies in any of the four domain boundary surfaces, and
	   that surface is conducting, then set the emf on that edge to zero. */
	if (g2 == -1 &&
	    m2->boundary_condition_lower[a2] == M2_BOUNDARY_REFLECTING) {
	  E->emf = 0.0;
	}
	if (g3 == -1 &&
	    m2->boundary_condition_lower[a3] == M2_BOUNDARY_REFLECTING) {
	  E->emf = 0.0;
	}
	if (g2 == m2->domain_resolution[a2] &&
	    m2->boundary_condition_upper[a2] == M2_BOUNDARY_REFLECTING) {
	  E->emf = 0.0;
	}
	if (g3 == m2->domain_resolution[a3] &&
	    m2->boundary_condition_upper[a3] == M2_BOUNDARY_REFLECTING) {
	  E->emf = 0.0;
	}
      }
    }
  }
  return 0;
}



int m2sim_exchange_flux(m2sim *m2, double dt)
{
  struct mesh_cell *cells[2];
  struct mesh_face *faces[4], *F;
  struct mesh_edge *E;
  int n, d, q;
  for (d=1; d<=3; ++d) {
    for (n=0; n<m2->mesh.faces_shape[d][0]; ++n) {
      F = &m2->mesh.faces[d][n];
      mesh_face_cells(&m2->mesh, F, cells);
      for (q=0; q<5; ++q) {
        if (cells[0]) cells[0]->consA[q] -= dt * F->area * F->flux[q];
        if (cells[1]) cells[1]->consA[q] += dt * F->area * F->flux[q];
      }
    }
    for (n=0; n<m2->mesh.edges_shape[d][0]; ++n) {
      E = &m2->mesh.edges[d][n];
      mesh_edge_faces(&m2->mesh, E, faces);
      if (faces[0]) faces[0]->BfluxA -= E->emf * dt;
      if (faces[1]) faces[1]->BfluxA += E->emf * dt;
      if (faces[2]) faces[2]->BfluxA += E->emf * dt;
      if (faces[3]) faces[3]->BfluxA -= E->emf * dt;
    }
  }
  return 0;
}



int m2sim_calculate_gradient(m2sim *m2)
{
  struct mesh_cell *C[3];
  int *dims = m2->mesh.cells_shape;
  int n, d, q, i;
  int err = 0;
  int set_grad_to_zero;
  int I[4];
  int m[3];
  double x[3];
  double y[24];
  double plm = m2->plm_parameter;
  double *G;
  for (d=1; d<=3; ++d) {
    if (m2->domain_resolution[d] == 1) {
      continue;
    }
    for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
      mesh_linear_to_multi(n, dims, I);
      I[d] -= 1; m[0] = mesh_multi_to_linear(dims, I); I[d] += 1;
      I[d] += 1; m[2] = mesh_multi_to_linear(dims, I); I[d] -= 1;
      m[1] = n;
      C[1] = m2->mesh.cells + m[1];
      switch (d) {
      case 1: G = C[1]->grad1; break;
      case 2: G = C[1]->grad2; break;
      case 3: G = C[1]->grad3; break;
      }
      if (m2->gradient_estimation_method == M2_GRADIENT_ESTIMATION_PCM) {
	set_grad_to_zero = 1;
      }
      else if (m[0] == -1 || m[2] == -1) { /* bordering the grid edge */
	set_grad_to_zero = 1;
      }
      else if (m2->gradient_estimation_method == M2_GRADIENT_ESTIMATION_PLM) {
	set_grad_to_zero = 0;
	for (i=0; i<3; ++i) {
	  C[i] = m2->mesh.cells + m[i];
	  x[i] = C[i]->x[d];
	  m2sim_to_interpolated(m2, &C[i]->prim, &C[i]->aux, &y[i], 3);
	}
	for (q=0; q<8; ++q) {
	  G[q] = plm_gradient(x, &y[3*q], plm, &err);
	}
	if (err > 0) {
	  m2sim_error_cell(m2, C[1], M2_ERROR_GRADIENT_ESTIMATION);
	  set_grad_to_zero = 1;
	}
      }
      else {
	set_grad_to_zero = 1;
	m2sim_error_general(m2, "invalid gradient estimation method");
      }
      if (set_grad_to_zero) {
	for (q=0; q<8; ++q) {
	  G[q] = 0.0;
	}
      }
    }
  }
  return 0;
}



int m2sim_add_source_terms(m2sim *m2, double dt)
{
  struct mesh_cell *C;
  struct mesh_face *F[6];
  int n;
  double x0[4];
  double x1[4];
  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    C = &m2->mesh.cells[n];
    mesh_cell_faces(&m2->mesh, C, F);
    x0[0] = 0.0;
    x1[0] = dt;
    x0[1] = F[0]->x[1]; /* figure out the rectangular extent of the cell */
    x1[1] = F[1]->x[1];
    x0[2] = F[2]->x[2];
    x1[2] = F[3]->x[2];
    x0[3] = F[4]->x[3];
    x1[3] = F[5]->x[3];
    m2aux_add_geometrical_source_terms(&C->aux, x0, x1, C->consA);
  }
  /* TODO: put in source terms from stirring here */
  return 0;
}



int m2sim_cache_conserved(m2sim *m2)
{
  struct mesh_cell *C;
  struct mesh_face *F;
  int n, d, q;
  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    C = &m2->mesh.cells[n];
    for (q=0; q<5; ++q) {
      C->consB[q] = C->consA[q];
    }
  }
  for (d=1; d<=3; ++d) {
    for (n=0; n<m2->mesh.faces_shape[d][0]; ++n) {
      F = &m2->mesh.faces[d][n];
      F->BfluxB = F->BfluxA;
    }
  }
  return 0;
}



int m2sim_average_runge_kutta(m2sim *m2, double b)
{
  struct mesh_cell *C;
  struct mesh_face *F;
  int n, d, q;
  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    C = &m2->mesh.cells[n];
    for (q=0; q<5; ++q) {
      C->consA[q] = b * C->consA[q] + (1.0 - b) * C->consB[q];
    }
  }
  for (d=1; d<=3; ++d) {
    for (n=0; n<m2->mesh.faces_shape[d][0]; ++n) {
      F = &m2->mesh.faces[d][n];
      F->BfluxA = b * F->BfluxA + (1.0 - b) * F->BfluxB;
    }
  }
  return 0;
}



int m2sim_magnetic_flux_to_cell_center(m2sim *m2)
/*
 * Convert the face-centered magnetic flux to a pointwise magnetic field at the
 * cell center. Although the left-most cell along each axis is given a
 * cell-centered field value (using only its +1/2 face), it should never be used
 * since the field value at its (-1/2) face is unknown.
 */
{
  struct mesh_cell *C;
  struct mesh_face *F[6];
  int n;
  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    C = &m2->mesh.cells[n];
    mesh_cell_faces(&m2->mesh, C, F);
    C->prim.B1 = (F[0]->BfluxA + F[1]->BfluxA) / (F[0]->area + F[1]->area);
    C->prim.B2 = (F[2]->BfluxA + F[3]->BfluxA) / (F[2]->area + F[3]->area);
    C->prim.B3 = (F[4]->BfluxA + F[5]->BfluxA) / (F[4]->area + F[5]->area);
  }
  return 0;
}



int m2sim_from_conserved_all(m2sim *m2)
/* It is required before calling this function that magnetic fields have already
   been averaged from faces to the cell centers. */
{
  struct mesh_cell *C;
  double B[4];
  int error, num_bad=0;
  int n;
  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    C = &m2->mesh.cells[n];
    B[0] = 0.0;
    B[1] = C->prim.B1;
    B[2] = C->prim.B2;
    B[3] = C->prim.B3;
    error = m2sim_from_conserved(m2, C->consA, B, NULL, C->volume,
                                 &C->aux, &C->prim);
    num_bad += (error != 0);
  }
  global_sum(m2, &num_bad);
  return num_bad;
}



int m2sim_from_primitive_all(m2sim *m2)
{
  struct mesh_cell *C;
  int n;
  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    C = &m2->mesh.cells[n];
    m2sim_from_primitive(m2, &C->prim, NULL, NULL, C->volume,
                         C->consA, &C->aux);
  }
  return 0;
}



double m2sim_minimum_courant_time(m2sim *m2)
{
  /* This function assumes that cells at the MPI and periodic boundaries have
     valid aux data. */
  struct mesh_cell *C;
  struct mesh_face *faces[6];
  int error;
  int n;
  double dt=-1;
  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    C = &m2->mesh.cells[n];
    mesh_cell_faces(&m2->mesh, C, faces);
    double dx1 = faces[1]->x[1] - faces[0]->x[1];
    double dx2 = faces[3]->x[2] - faces[2]->x[2];
    double dx3 = faces[5]->x[3] - faces[4]->x[3];
    double L = M2_MIN3(dx1, dx2, dx3);
    double a = m2aux_maximum_wavespeed(&C->aux, &error);
    if (dt < 0 || L/a < dt) {
      dt = L/a;
    }
  }
  global_min(m2, &dt);
  return dt;
}



int m2sim_enforce_user_constraint(m2sim *m2)
{
  return 0;
}



int m2sim_runge_kutta_substep(m2sim *m2, double dt, double rkparam)
{
  int err = 0;
  err += m2sim_calculate_gradient(m2);
  /* Outermost 1 cells have invalid gradient data. */

  err += m2sim_synchronize_cells(m2);
  /* Communicate outermost 1 cells across MPI and periodic boundaries. All
     cells have valid gradient data. */

  err += m2sim_calculate_flux(m2);
  /* Outermost 1 faces at MPI and periodic boundaries have invalid Godunov
     fluxes. */
  /* TODO: set physical BC's on Godunov fluxes at domain and topological
     boundary */

  err += m2sim_synchronize_faces(m2);
  /* Communicate outermost 1 faces across MPI and periodic boundaries. All
     faces have valid Godunov and magnetic fluxes. */

  err += m2sim_calculate_emf(m2);
  /* TODO: set physical BC's on Godunov EMF's at domain and topological
     boundary */
  /* TODO: source terms from stirring the current go after the calculation of
     Godunov EMF's */
  /* EMF's on MPI and periodic boundary are still wrong, because they only
     account for 3 of 4 faces on a planar boundary, and 2 of 4 faces at the
     intersection of two planar boundaries. */

  err += m2sim_exchange_flux(m2, dt);
  /* TODO: communicate outermost 1 faces across MPI and periodic boundaries, to
     account for incorrect EMF's on boundary faces */
  /* Conserved quantites in cells on MPI and periodic boundaries have the
     correct values because Godunov fluxes on even their face was
     communicated. */

  err += m2sim_add_source_terms(m2, dt);
  /* Valid for all cells. */

  err += m2sim_enforce_user_constraint(m2);
  /* This is where additional user constraints are being enforced. These might
     include driving the momentum, or a "soft" boundary condition where the
     solution is driven toward something user prescribed. Such constraints may
     be enforced on the conserved quantities consA and BfluxA only.
     Modifications to cell-centered magnetic field, prim, or aux will be
     over-written in the last two steps. */

  err += m2sim_average_runge_kutta(m2, rkparam);
  /* Applies to consA and BfluxA only, all cells and faces have valid conserved
     quantities. */

  err += m2sim_magnetic_flux_to_cell_center(m2);
  /* All cells have valid cell-centered fields. */

  err += m2sim_from_conserved_all(m2);
  /* All cells and faces are valid. */
  return err;
}



void m2sim_drive(m2sim *m2)
{
  double dt;
  int err = 0, rk_order = m2->rk_order;
  double kzps; /* kilozones per second */
  clock_t start_cycle, stop_cycle;



  m2sim_cache_conserved(m2);
  dt = m2->cfl_parameter * m2sim_minimum_courant_time(m2);
  m2->status.time_step = dt;



  /* ======================================================================== */
  /* Start the clock                                                          */
  /* ======================================================================== */
  start_cycle = clock();



  /* ======================================================================== */
  /* Take Runge-Kutta substeps                                                */
  /* ======================================================================== */
  switch (rk_order) {
  case 1:
    err = m2sim_runge_kutta_substep(m2, dt, 1.0); if (err) {err=1; break;}
    break;
  case 2:
    err = m2sim_runge_kutta_substep(m2, dt, 1.0); if (err) {err=1; break;}
    err = m2sim_runge_kutta_substep(m2, dt, 0.5); if (err) {err=2; break;}
    break;
  case 3:
    err = m2sim_runge_kutta_substep(m2, dt, 1.0/1.0); if (err) {err=1; break;}
    err = m2sim_runge_kutta_substep(m2, dt, 1.0/4.0); if (err) {err=2; break;}
    err = m2sim_runge_kutta_substep(m2, dt, 2.0/3.0); if (err) {err=3; break;}
    break;
  }



  m2->status.iteration_number += 1;
  m2->status.time_simulation += dt;
  m2->status.error_code = 0;



  stop_cycle = clock();
  kzps = 1e-3 * m2->local_grid_size[0] / (stop_cycle - start_cycle) *
    CLOCKS_PER_SEC;



  /* ======================================================================== */
  /* Print status message                                                     */
  /* ======================================================================== */
  snprintf(m2->status.message, 1024,
	   "%08d(%d): t=%5.4f dt=%5.4e kzps=%4.3f [%d %d %d %d %d]",
	   m2->status.iteration_number,
	   m2->status.drive_attempt,
	   m2->status.time_simulation,
	   dt,
	   kzps,
	   m2->status.num_failures[0],
	   m2->status.num_failures[1],
	   m2->status.num_failures[2],
	   m2->status.num_failures[3],
	   m2->status.num_failures[4]);
}






















static void global_sum(m2sim *m2, int *N)
{
#if M2_HAVE_MPI
  if (m2->cart_comm) {
    MPI_Comm cart_comm = *((MPI_Comm *) m2->cart_comm);
    MPI_Allreduce(MPI_IN_PLACE, N, 1, MPI_INT, MPI_SUM, cart_comm);
  }
#endif
}



static void global_min(m2sim *m2, double *x)
{
#if M2_HAVE_MPI
  if (m2->cart_comm) {
    MPI_Comm cart_comm = *((MPI_Comm *) m2->cart_comm);
    MPI_Allreduce(MPI_IN_PLACE, x, 1, MPI_DOUBLE, MPI_MIN, cart_comm);
  }
#endif
}



static double plm_gradient(double *xs, double *fs, double plm, int *error)
{
#define sign M2_SIGN
#define min3 M2_MIN3
  double xL = xs[0];
  double x0 = xs[1];
  double xR = xs[2];
  double yL = fs[0];
  double y0 = fs[1];
  double yR = fs[2];
  double a = (yL - y0) / (xL - x0) * plm;
  double b = (yL - yR) / (xL - xR);
  double c = (y0 - yR) / (x0 - xR) * plm;
  double sa = sign(a), sb = sign(b), sc = sign(c);
  double fa = fabs(a), fb = fabs(b), fc = fabs(c);
  double g = 0.25 * fabs(sa + sb) * (sa + sc) * min3(fa, fb, fc);
  *error += g != g;
  return g;
#undef sign
#undef min3
}



int global_index_cell(m2sim *m2, struct mesh_cell *C, int d)
{
  switch (d) {
  case 1: return C->global_i;
  case 2: return C->global_j;
  case 3: return C->global_k;
  default: return 0;
  }
}



int global_index_face(m2sim *m2, struct mesh_face *F, int d)
{
  int I[4];
  mesh_linear_to_multi(F->id, m2->mesh.faces_shape[d], I);
  return I[d] + m2->local_grid_start[d];
}



int global_index_edge(m2sim *m2, struct mesh_edge *E, int d)
{
  int I[4];
  mesh_linear_to_multi(E->id, m2->mesh.edges_shape[d], I);
  return I[d] + m2->local_grid_start[d];
}
