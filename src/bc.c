#include <math.h>
#include "m2.h"


static inline int global_index_face(m2sim *m2, struct mesh_face *F, int d);


void bc_face_polar(m2sim *m2, int d, char which)
{
  struct mesh_face *F;
  int n, i, q;
  for (n=0; n<m2->mesh.faces_shape[d][0]; ++n) {
    F = &m2->mesh.faces[d][n];
    i = global_index_face(m2, F, d);
    if ((which == 'U' && i == m2->domain_resolution[d]) ||
	(which == 'L' && i == 0)) {
      for (q=0; q<8; ++q) {
	F->flux[q] = 0.0;
      }
    }
  }
}



void bc_face_reflecting1(m2sim *m2, int d, char which)
{
  struct mesh_face *F, *G;
  int B1, B2, B3;
  int S1, S2, S3;
  int n, i;

  switch (d) {
  case 1: B1 = B11; B2 = B22; B3 = B33; S1 = S11; S2 = S22; S3 = S33; break;
  case 2: B1 = B22; B2 = B33; B3 = B11; S1 = S22; S2 = S33; S3 = S11; break;
  case 3: B1 = B33; B2 = B11; B3 = B22; S1 = S33; S2 = S11; S3 = S22; break;
  }

  for (n=0; n<m2->mesh.faces_shape[d][0]; ++n) {
    F = &m2->mesh.faces[d][n];
    i = global_index_face(m2, F, d);
    if (which == 'U' && i == m2->domain_resolution[d]) {
      G = F - m2->mesh.faces_stride[d][d];
    }
    else if (which == 'L' && i == 0) {
      G = F + m2->mesh.faces_stride[d][d];
    }
    else {
      G = NULL;
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



void bc_face_reflecting2(m2sim *m2, int d, char which)
{
  struct mesh_face *F;
  struct mesh_cell *C, *cells[2], D;
  int B1, B2, B3;
  int S1, S2, S3;
  int n, i;

  switch (d) {
  case 1: B1 = B11; B2 = B22; B3 = B33; S1 = S11; S2 = S22; S3 = S33; break;
  case 2: B1 = B22; B2 = B33; B3 = B11; S1 = S22; S2 = S33; S3 = S11; break;
  case 3: B1 = B33; B2 = B11; B3 = B22; S1 = S33; S2 = S11; S3 = S22; break;
  }

  for (n=0; n<m2->mesh.faces_shape[d][0]; ++n) {
    F = &m2->mesh.faces[d][n];
    i = global_index_face(m2, F, d);
    if (which == 'U' && i == m2->domain_resolution[d]) {
      mesh_face_cells(&m2->mesh, F, cells);
      C = cells[0];
    }
    else if (which == 'L' && i == 0) {
      mesh_face_cells(&m2->mesh, F, cells);
      C = cells[1];
    }
    else {
      C = NULL;
    }
    if (C) {
      D = *C;
      D.prim.v1 *= -1.0;
      D.prim.B2 *= -1.0;
      D.prim.B3 *= -1.0;
      m2sim_from_primitive(m2, &D.prim, NULL, NULL, D.volume,
			   D.consA, &D.aux);
      D.grad1[1] *= -1; /* v-transverse */
      D.grad1[2] *= -1; /* v-transverse */
      D.grad1[3] *= -1; /* B-normal (doesn't matter) */
      D.grad1[6] *= -1; /* density */
      D.grad1[7] *= -1; /* pressure */
      D.x[1] = F->x[1] - (C->x[1] - F->x[1]);
      m2sim_riemann_solver(m2, F, &D);
    }
  }
}





int global_index_face(m2sim *m2, struct mesh_face *F, int d)
{
  int I[4];
  mesh_linear_to_multi(F->id, m2->mesh.faces_shape[F->axis], I);
  return I[d] + m2->local_grid_start[d];
}
