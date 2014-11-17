#include <stddef.h>
#include <memory.h>
#include <mpi.h>
#include "m2.h"


#define copy_cell(C,D)					\
  do {							\
    C->prim = D->prim;					\
    C->aux  = D->aux;					\
    memcpy(C->consA, D->consA, 5*sizeof(double));	\
    memcpy(C->consB, D->consB, 5*sizeof(double));	\
    memcpy(C->grad1, D->grad1, 8*sizeof(double));	\
    memcpy(C->grad2, D->grad2, 8*sizeof(double));	\
    memcpy(C->grad3, D->grad3, 8*sizeof(double));	\
  } while (0)						\


#define copy_face(F,G)						\
  do {								\
    F->BfluxA = G->BfluxA;					\
    F->BfluxB = G->BfluxB;					\
    memcpy(F->flux, G->flux, 8*sizeof(double));			\
  } while (0)							\



static inline int global_index_cell(m2sim *m2, struct mesh_cell *C, int axis) {
  switch (axis) {
  case 1: return C->global_i;
  case 2: return C->global_j;
  case 3: return C->global_k;
  default: return 0;
  }
}



/* static inline int global_index_face(m2sim *m2, struct mesh_face *F, int d) */
/* { */
/*   int I[4]; */
/*   mesh_linear_to_multi(F->id, m2->mesh.faces_shape[d], I); */
/*   return I[d] + m2->local_grid_start[d]; */
/* } */



struct mpi_types
{
  MPI_Datatype vol;
  MPI_Datatype send_R[3], recv_R[3];
  MPI_Datatype send_L[3], recv_L[3];
} ;

void m2sim_create_mpi_types(m2sim *m2)
{
  m2->mpi_types = malloc(sizeof(struct mpi_types));
  struct mpi_types *MT = (struct mpi_types *) m2->mpi_types;
  int B[20];
  MPI_Aint D[20];
  MPI_Datatype T[20];
  int *L = m2->local_grid_size;
  int *Ng0 = m2->number_guard_zones0;
  int *Ng1 = m2->number_guard_zones1;
  int o = MPI_ORDER_C;
  int n;

  /* TODO: use the MPI upper bound type to improve safety
     http://www.mpich.org/static/docs/v3.1/www3/MPI_Type_struct.html
  */
  B[ 0] = 4;   D[ 0] = offsetof(m2vol, global_index); T[ 0] = MPI_INT;
  B[ 1] = 1;   D[ 1] = offsetof(m2vol, Bflux1A);      T[ 1] = MPI_DOUBLE;
  B[ 2] = 1;   D[ 2] = offsetof(m2vol, Bflux2A);      T[ 2] = MPI_DOUBLE;
  B[ 3] = 1;   D[ 3] = offsetof(m2vol, Bflux3A);      T[ 3] = MPI_DOUBLE;
  B[ 4] = 1;   D[ 4] = offsetof(m2vol, Bflux1B);      T[ 4] = MPI_DOUBLE;
  B[ 5] = 1;   D[ 5] = offsetof(m2vol, Bflux2B);      T[ 5] = MPI_DOUBLE;
  B[ 6] = 1;   D[ 6] = offsetof(m2vol, Bflux3B);      T[ 6] = MPI_DOUBLE;
  B[ 7] = 5;   D[ 7] = offsetof(m2vol, consA);        T[ 7] = MPI_DOUBLE;
  B[ 8] = 5;   D[ 8] = offsetof(m2vol, consB);        T[ 8] = MPI_DOUBLE;
  B[ 9] = 8;   D[ 9] = offsetof(m2vol, flux1);        T[ 9] = MPI_DOUBLE;
  B[10] = 8;   D[10] = offsetof(m2vol, flux2);        T[10] = MPI_DOUBLE;
  B[11] = 8;   D[11] = offsetof(m2vol, flux3);        T[11] = MPI_DOUBLE;
  B[12] = 8;   D[12] = offsetof(m2vol, grad1);        T[12] = MPI_DOUBLE;
  B[13] = 8;   D[13] = offsetof(m2vol, grad2);        T[13] = MPI_DOUBLE;
  B[14] = 8;   D[14] = offsetof(m2vol, grad3);        T[14] = MPI_DOUBLE;
  B[15] = 1;   D[15] = offsetof(m2vol, emf1);         T[15] = MPI_DOUBLE;
  B[16] = 1;   D[16] = offsetof(m2vol, emf2);         T[16] = MPI_DOUBLE;
  B[17] = 1;   D[17] = offsetof(m2vol, emf3);         T[17] = MPI_DOUBLE;
  B[18] = sizeof(m2prim); D[18] = offsetof(m2vol, prim); T[18] = MPI_BYTE;
  B[19] = sizeof(m2aux);  D[19] = offsetof(m2vol, aux);  T[19] = MPI_BYTE;

  MPI_Type_create_struct(20, B, D, T, &MT->vol);
  MPI_Type_commit(&MT->vol);

  MPI_Datatype v = MT->vol;
  MPI_Datatype *sp = MT->send_R;
  MPI_Datatype *sm = MT->send_L;
  MPI_Datatype *rp = MT->recv_R;
  MPI_Datatype *rm = MT->recv_L;

  for (n=0; n<3; ++n) {
    int count_sp[3] = { L[1], L[2], L[3] };
    int count_rp[3] = { L[1], L[2], L[3] };
    int count_sm[3] = { L[1], L[2], L[3] };
    int count_rm[3] = { L[1], L[2], L[3] };
    int sizes_sp[3] = { L[1], L[2], L[3] };
    int sizes_rp[3] = { L[1], L[2], L[3] };
    int sizes_sm[3] = { L[1], L[2], L[3] };
    int sizes_rm[3] = { L[1], L[2], L[3] };
    int start_sp[3] = { 0, 0, 0 };
    int start_rp[3] = { 0, 0, 0 };
    int start_sm[3] = { 0, 0, 0 };
    int start_rm[3] = { 0, 0, 0 };

    sizes_sp[n] = Ng1[n+1];
    sizes_rp[n] = Ng1[n+1];
    sizes_sm[n] = Ng0[n+1];
    sizes_rm[n] = Ng0[n+1];

    start_sp[n] = L[n+1] - 2*Ng1[n+1];
    start_rp[n] = L[n+1] - 1*Ng1[n+1];
    start_sm[n] = Ng0[n+1];
    start_rm[n] = 0;

    if (Ng1[n+1]) {
      MPI_Type_create_subarray(3, count_sp, sizes_sp, start_sp, o, v, &sp[n]);
      MPI_Type_create_subarray(3, count_rp, sizes_rp, start_rp, o, v, &rp[n]);
      MPI_Type_commit(&sp[n]);
      MPI_Type_commit(&rp[n]);
    }
    else {
      sp[n] = MPI_BYTE; /* indicates empty; it's never used */
      rp[n] = MPI_BYTE;
    }
    if (Ng0[n+1]) {
      MPI_Type_create_subarray(3, count_sm, sizes_sm, start_sm, o, v, &sm[n]);
      MPI_Type_create_subarray(3, count_rm, sizes_rm, start_rm, o, v, &rm[n]);
      MPI_Type_commit(&sm[n]);
      MPI_Type_commit(&rm[n]);
    }
    else {
      sm[n] = MPI_BYTE;
      rm[n] = MPI_BYTE;
    }
  }
}

void m2sim_delete_mpi_types(m2sim *m2)
{
  struct mpi_types *MT = (struct mpi_types *) m2->mpi_types;
  MPI_Datatype *sp = MT->send_R;
  MPI_Datatype *sm = MT->send_L;
  MPI_Datatype *rp = MT->recv_R;
  MPI_Datatype *rm = MT->recv_L;
  int n;

  MPI_Type_free(&MT->vol);

  for (n=0; n<3; ++n) {
    if (sp[n] != MPI_BYTE) MPI_Type_free(&sp[n]);
    if (sm[n] != MPI_BYTE) MPI_Type_free(&sm[n]);
    if (rp[n] != MPI_BYTE) MPI_Type_free(&rp[n]);
    if (rm[n] != MPI_BYTE) MPI_Type_free(&rm[n]);
  }

  free(m2->mpi_types);
  m2->mpi_types = NULL;
}



int m2sim_synchronize_cells(m2sim *m2)
{
  /* TODO: implement the MPI version of this */
  struct mesh_cell *C, *D;
  int n, d, gind;
  for (d=1; d<=3; ++d) {
    for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
      C = m2->mesh.cells + n;
      gind = global_index_cell(m2, C, d);
      if (m2->boundary_condition_lower[d] == M2_BOUNDARY_PERIODIC) {
	if (gind < 0) {
	  D = C + m2->domain_resolution[d] * m2->mesh.cells_stride[d];
	  copy_cell(C, D);
	}
      }
      if (m2->boundary_condition_upper[d] == M2_BOUNDARY_PERIODIC) {
	if (gind >= m2->domain_resolution[d]) {
	  D = C - m2->domain_resolution[d] * m2->mesh.cells_stride[d];
	  copy_cell(C, D);
	}
      }
    }
  }
  return 0;
}



int m2sim_synchronize_faces(m2sim *m2)
{
  /* TODO: implement serial and MPI version of this */
  /* for (d=1; d<=3; ++d) { */
  /*   for (n=0; n<m2->mesh.faces_shape[d][0]; ++n) { */
  /*     F = m2->mesh.faces[d] + n; */

  /*     I[1] = global_index_face(m2, F, 1); */
  /*     I[2] = global_index_face(m2, F, 2); */
  /*     I[3] = global_index_face(m2, F, 3); */

  /*     if (I[d] < 0) { */
  /* 	I[d] += m2->mesh.cells_shape[d]; */
  /* 	m = mesh_multi_to_linear(m2->mesh.faces_shape[d], I); */
  /* 	G = m2->mesh.faces[d] + m; */
  /* 	copy_face(F, G); */
  /* 	I[d] -= m2->mesh.cells_shape[d]; */
  /*     } */

  /*     else if (I[d] >= m2->domain_resolution[d]) { */
  /* 	I[d] -= m2->mesh.cells_shape[d]; */
  /* 	m = mesh_multi_to_linear(m2->mesh.faces_shape[d], I); */
  /* 	G = m2->mesh.faces[d] + m; */
  /* 	copy_face(F, G); */
  /* 	I[d] += m2->mesh.cells_shape[d]; */
  /*     } */
  /*   } */
  /* } */
  return 0;
}



int m2sim_synchronize_guard__(m2sim *m2)
{
  return 0; /* TODO */
  if (m2->cart_comm == NULL) {
    MSG(FATAL, "guard zone synch without MPI is not yet implemented");
  }
  struct mpi_types *MT = (struct mpi_types *) m2->mpi_types;
  MPI_Comm cart_comm = *((MPI_Comm *) m2->cart_comm);
  MPI_Status status;
  int rankS;
  int rankD;
  int tagS = 0;
  int tagD = 0;
  int n;

  MPI_Datatype *sp = MT->send_R;
  MPI_Datatype *sm = MT->send_L;
  MPI_Datatype *rp = MT->recv_R;
  MPI_Datatype *rm = MT->recv_L;

  for (n=0; n<3; ++n) {
    if (m2->domain_resolution[n+1] == 1) continue;

    /* send to right, receive from left */
    MPI_Cart_shift(cart_comm, n, +1, &rankS, &rankD);
    MPI_Sendrecv(m2->volumes, 1, sp[n], rankD, tagD,
		 m2->volumes, 1, rm[n], rankS, tagS,
		 cart_comm, &status);

    /* send to left, receive from right */
    MPI_Cart_shift(cart_comm, n, -1, &rankS, &rankD);
    MPI_Sendrecv(m2->volumes, 1, sm[n], rankD, tagD,
		 m2->volumes, 1, rp[n], rankS, tagS,
		 cart_comm, &status);
  }

  /* because MPI overwrites aux.m2 */
  for (n=0; n<m2->local_grid_size[0]; ++n) {
    m2->volumes[n].aux.m2 = m2;
  }
  return 0;
}
