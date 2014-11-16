#include <string.h>
#include <math.h>
#include "m2.h"
#if M2_HAVE_MPI
#include <mpi.h>
#endif



static void global_sum(m2sim *m2, int *N);
static void global_min(m2sim *m2, double *x);
static double plm_gradient(double *xs, double *fs, double plm, int *error);



int m2sim_calculate_flux(m2sim *m2)
{
  int n, d;
  for (d=1; d<=3; ++d) {
    for (n=0; n<m2->mesh.faces_shape[d][0]; ++n) {
      m2sim_riemann_solver(m2, &m2->mesh.faces[d][n]);
    }
  }
  return 0;
}



int m2sim_calculate_emf(m2sim *m2)
{
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



int m2sim_exchange_flux(m2sim *m2, double dt)
{
  struct mesh_cell *cells[2];
  struct mesh_face *faces[4], *F;
  struct mesh_edge *E;
  int n, d, q;
  for (d=1; d<=3; ++d) {
    for (n=0; n<m2->mesh.faces_shape[d][0]; ++n) {
      F = m2->mesh.faces[d] + n;
      mesh_face_cells(&m2->mesh, F, cells);
      for (q=0; q<5; ++q) {
        if (cells[0]) cells[0]->consA[q] -= dt * F->area * F->flux[q];
        if (cells[1]) cells[1]->consA[q] += dt * F->area * F->flux[q];
      }
    }
    for (n=0; n<m2->mesh.edges_shape[d][0]; ++n) {
      E = m2->mesh.edges[d] + n;
      mesh_edge_faces(&m2->mesh, E, faces);
      if (faces[0]) faces[0]->BfluxA -= E->emf * dt; /* TODO: double-check sign */
      if (faces[1]) faces[1]->BfluxA += E->emf * dt;
      if (faces[2]) faces[2]->BfluxA -= E->emf * dt;
      if (faces[3]) faces[3]->BfluxA += E->emf * dt;
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
    C = m2->mesh.cells + n;
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
    C = m2->mesh.cells + n;
    for (q=0; q<5; ++q) {
      C->consB[q] = C->consA[q];
    }
  }
  for (d=1; d<=3; ++d) {
    for (n=0; n<m2->mesh.faces_shape[d][0]; ++n) {
      F = m2->mesh.faces[d] + n;
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
    C = m2->mesh.cells + n;
    for (q=0; q<5; ++q) {
      C->consA[q] = b * C->consA[q] + (1.0 - b) * C->consB[q];
    }
  }
  for (d=1; d<=3; ++d) {
    for (n=0; n<m2->mesh.faces_shape[d][0]; ++n) {
      F = m2->mesh.faces[d] + n;
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
    C = m2->mesh.cells + n;
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
    C = m2->mesh.cells + n;
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
    C = m2->mesh.cells + n;
    m2sim_from_primitive(m2, &C->prim, NULL, NULL, C->volume,
                         C->consA, &C->aux);
  }
  return 0;
}



double m2sim_minimum_courant_time(m2sim *m2)
{
  struct mesh_cell *C;
  struct mesh_face *faces[6];
  int error;
  int n;
  double dt=-1;
  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    /* TODO: unless guards have already been synchronized, this function must
       skip non-interior cells. */
    C = m2->mesh.cells + n;
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



int m2sim_enforce_boundary_condition(m2sim *m2)
{
  if (m2->boundary_conditions_cell) {
    m2->boundary_conditions_cell(m2);
  }
  return 0;
}



int m2sim_runge_kutta_substep(m2sim *m2, double dt, double rkparam)
{
  int err = 0;
  err += m2sim_calculate_gradient(m2);
  err += m2sim_calculate_flux(m2); /* FAILURES: bad eigenvalues */
  err += m2sim_synchronize_guard(m2); /* so that Godunov fluxes are communicated */
  err += m2sim_calculate_emf(m2);
  err += m2sim_exchange_flux(m2, dt);
  err += m2sim_add_source_terms(m2, dt);
  err += m2sim_average_runge_kutta(m2, rkparam);
  err += m2sim_enforce_boundary_condition(m2);
  err += m2sim_magnetic_flux_to_cell_center(m2);
  err += m2sim_from_conserved_all(m2);
  err += m2sim_synchronize_guard(m2);
  return err;
}



void m2sim_drive(m2sim *m2)
{
  double dt;
  int err = 0, rk_order = m2->rk_order;
  double kzps; /* kilozones per second */
  char checkpoint_name[1024];

  m2sim_run_analysis(m2);
  m2sim_cache_conserved(m2);
  dt = m2->cfl_parameter * m2sim_minimum_courant_time(m2);
  m2->status.time_step = dt;

  double cad = m2->cadence_checkpoint_tpl;
  double now = m2->status.time_simulation;
  double las = m2->status.time_last_checkpoint_tpl;

  if (cad > 0.0 && now - las > cad) {
    m2->status.checkpoint_number_tpl += 1;
    m2->status.time_last_checkpoint_tpl += cad;
    snprintf(checkpoint_name, sizeof(checkpoint_name), "chkpt.%04d.m2",
             m2->status.checkpoint_number_tpl);
    m2sim_save_checkpoint_tpl(m2, checkpoint_name);
  }



  /* ======================================================================== */
  /* Start the clock                                                          */
  /* ======================================================================== */
  clock_t start_cycle, stop_cycle;
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
