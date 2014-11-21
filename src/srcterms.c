#include <math.h>
#include "m2.h"

static inline int global_index_face(m2sim *m2, struct mesh_face *F, int d);

void srcterms_edge_magnetosphere(m2sim *m2)
{
  struct mesh_edge *E;
  struct mesh_face *F[4];
  int n;

  /* theta-edges */
  for (n=0; n<m2->mesh.edges_shape[2][0]; ++n) {
    E = &m2->mesh.edges[2][n];
    mesh_edge_faces(&m2->mesh, E, F);
    double r = E->x[1];
    double t = E->x[2];
    double B = 0.0;
    double A = 0.0;
    double R = r * sin(t);
    double vc = 1.0; /* velocity at light cylinder */
    double rc = 3.0;
    double dr = 0.5;
    double omega = vc/rc * 0.5 * (1 - tanh((r - rc) / dr));
    /* E-theta = -v-phi * B-r */
    if (F[2]) {B += F[2]->BfluxA; A += F[2]->area;} /* r-faces */
    if (F[3]) {B += F[3]->BfluxA; A += F[3]->area;}
    E->emf -= E->length * omega * R * B/A;
  }


  /* r-edges */
  for (n=0; n<m2->mesh.edges_shape[1][0]; ++n) {
    E = &m2->mesh.edges[1][n];
    mesh_edge_faces(&m2->mesh, E, F);
    double r = E->x[1];
    double t = E->x[2];
    double B = 0.0;
    double A = 0.0;
    double R = r * sin(t);
    double vc = 1.0; /* velocity at light cylinder */
    double rc = 3.0;
    double dr = 0.5;
    double omega = vc/rc * 0.5 * (1 - tanh((r - rc) / dr));
    /* E-r = v-phi * B-theta */
    if (F[0]) {B += F[0]->BfluxA; A += F[0]->area;} /* theta-faces */
    if (F[1]) {B += F[1]->BfluxA; A += F[1]->area;}
    E->emf += E->length * omega * R * B/A;
  }
}



void srcterms_cell_massive_rotator(m2sim *m2, double dt)
{
  struct mesh_cell *C;
  struct mesh_face *F;
  int n;
  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    C = &m2->mesh.cells[n];
    int q;
    double r = C->x[1];
    double t = C->x[2];
    double R = r * sin(t);
    double vc = 1.0; /* velocity at light cylinder */
    double rc = 3.0;
    double dr = 0.15;
    double omega = vc/rc * 0.5 * (1 - tanh((R - rc) / dr));
    double tau = 1e-3 / omega;
    double E = omega > 0.001 ? exp(-dt/tau) : 1.0;
    m2prim P = C->prim;
    double U0[5];

    P.d = 0.1 * C->aux.magnetic_pressure;
    P.p = 0.1 * P.d;
    P.v1 = 0.0;
    P.v2 = 0.0;
    P.v3 = R * omega;

    if (C->global_i == 0) {
      E = 0.0;
    }

    m2sim_from_primitive(m2, &P, NULL, NULL, C->volume, U0, NULL);
    for (q=0; q<5; ++q) {
      C->consA[q] = E * C->consA[q] + (1 - E) * U0[q];
    }
  }

  for (n=0; n<m2->mesh.faces_shape[3][0]; ++n) {
    F = &m2->mesh.faces[3][n];
    double r = F->x[1];
    double t = F->x[2];
    double R = r * sin(t);
    double vc = 1.0; /* velocity at light cylinder */
    double rc = 3.0;
    double dr = 0.15;
    double omega = vc/rc * 0.5 * (1 - tanh((R - rc) / dr));
    double tau = 1e-3 / omega;
    double E = omega > 0.001 ? exp(-dt/tau) : 1.0;

    if (global_index_face(m2, F, 1) == 0) {
      E = 0.0;
    }

    F->BfluxA = E * F->BfluxA + (1 - E) * 0.0;
  }
}



void srcterms_cell_magnetarwind(m2sim *m2, double dt)
{
  struct mesh_cell *C;
  struct mesh_face *F;
  int n;
  double K, r, t;
  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    C = &m2->mesh.cells[n];
    r = C->x[1];
    t = C->x[2];
    K = exp(-pow(r - 2, 2)) * pow(sin(t), 2);
    C->consA[DDD] += C->volume * 1e1 * K * dt;
    C->consA[TAU] += C->volume * 1e2 * K * dt;
  }
  for (n=0; n<m2->mesh.faces_shape[3][0]; ++n) {
    F = &m2->mesh.faces[3][n];
    r = F->x[1];
    t = F->x[2];
    K = exp(-pow(r - 2, 2)) * pow(sin(t), 2);
    F->BfluxA += F->area * 5 * K * dt;
  }
}



void srcterms_cell_magnetarwind_inner_bnd(m2sim *m2, double dt)
{
  struct mesh_cell *C;
  struct mesh_face *F[6];
  int n;
  double B, r, t, rho = 1.0, sigma=0.25, gamma=10.0;
  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    C = &m2->mesh.cells[n];

    if (C->global_i != 0) continue;
      
    mesh_cell_faces(&m2->mesh, C, F);

    rho = 1.0;
    r = C->x[1];
    t = C->x[2];
    B = sqrt(sigma * rho) * sin(t);

    C->prim.B3 = B;
    C->prim.v1 = sqrt(1.0 - pow(gamma, -2));
    C->prim.d = 0.01;
    C->prim.p = 0.1 * C->prim.d;
    F[4]->BfluxA = B * F[4]->area;
    F[5]->BfluxA = B * F[5]->area;

    m2sim_from_primitive(m2, &C->prim, NULL, NULL, C->volume,
                         C->consA, &C->aux);
  }
}


   /* local function wind_inner_boundary_flux(V0) */
   /*    local nhat = m2lib.dvec4(0,1,0,0) */
   /*    local T = m2.status.time_simulation */
   /*    local t = 0.5 * (V0.x0[2] + V0.x1[2]) */
   /*    local p = 0.5 * (V0.x0[3] + V0.x1[3]) */
   /*    if V0.global_index[1] == -1 then */
   /*       local d = mps.d_wind */
   /*       local B = mps.B_wind */
   /*       local g = mps.g_wind */
   /*       local h = mps.g_pert */
   /*       g = g * (1.0 + h * math.sin(t) * math.sin(6*p)^2) */
   /* 	 g = g * (1.0 - math.exp(-T)) + math.exp(-T) -- ramp up */
   /* 	 B = B * (1.0 - math.exp(-T)) -- ramp up */
   /*       V0.prim.v1 =(1.0 - g^-2)^0.5 */
   /*       V0.prim.v2 = 0.0 */
   /*       V0.prim.v3 = 0.0 */
   /*       V0.prim.B1 = 0.0 */
   /*       V0.prim.B2 = 0.0 */
   /*       V0.prim.B3 = B * math.sin(t) */
   /*       V0.prim.d = d */
   /*       V0.prim.p = d * 0.01 */
   /*       V0:from_primitive() */
   /*       V0.flux1 = V0.aux:fluxes(nhat) */
   /*    else */
   /*       V0.flux1 = V0.aux:fluxes(nhat) */
   /*    end */
   /* end */





int global_index_face(m2sim *m2, struct mesh_face *F, int d)
{
  int I[4];
  mesh_linear_to_multi(F->id, m2->mesh.faces_shape[F->axis], I);
  return I[d] + m2->local_grid_start[d];
}

