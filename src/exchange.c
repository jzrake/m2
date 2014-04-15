#include <mpi.h>
#include "m2.h"

struct mpi_struct_member
{
  int blocklength;
  MPI_Aint displacement;
  MPI_Datatype type;
} ;

void m2sim_synchronize_guard(m2sim *m2)
{
  MPI_Comm cart_comm = *((MPI_Comm *) m2->cart_comm);
  MPI_Status status;
  MPI_Datatype vol;
  MPI_Datatype sp[3], rp[3];
  MPI_Datatype sm[3], rm[3];
  int order = MPI_ORDER_C;
  int rankS;
  int rankD;
  int tagS = 0;
  int tagD = 0;
  int *Ng0 = m2->number_guard_zones0;
  int *Ng1 = m2->number_guard_zones1;
  int *G = m2->domain_resolution;
  int *L = m2->local_grid_size;
  int n;
  int blocklengths[20];
  MPI_Aint displacements[20];
  MPI_Datatype types[20];
  struct mpi_struct_member vol_members[20] = {
    {4,              offsetof(m2vol, global_index), MPI_INT },
    {1,              offsetof(m2vol, Bflux1A),      MPI_DOUBLE},
    {1,              offsetof(m2vol, Bflux2A),      MPI_DOUBLE},
    {1,              offsetof(m2vol, Bflux3A),      MPI_DOUBLE},
    {1,              offsetof(m2vol, Bflux1B),      MPI_DOUBLE},
    {1,              offsetof(m2vol, Bflux2B),      MPI_DOUBLE},
    {1,              offsetof(m2vol, Bflux3B),      MPI_DOUBLE},
    {5,              offsetof(m2vol, consA),        MPI_DOUBLE},
    {5,              offsetof(m2vol, consB),        MPI_DOUBLE},
    {8,              offsetof(m2vol, flux1),        MPI_DOUBLE},
    {8,              offsetof(m2vol, flux2),        MPI_DOUBLE},
    {8,              offsetof(m2vol, flux3),        MPI_DOUBLE},
    {8,              offsetof(m2vol, grad1),        MPI_DOUBLE},
    {8,              offsetof(m2vol, grad2),        MPI_DOUBLE},
    {8,              offsetof(m2vol, grad3),        MPI_DOUBLE},
    {1,              offsetof(m2vol, emf1),         MPI_DOUBLE},
    {1,              offsetof(m2vol, emf2),         MPI_DOUBLE},
    {1,              offsetof(m2vol, emf3),         MPI_DOUBLE},
    {sizeof(m2prim), offsetof(m2vol, prim),         MPI_BYTE},
    {sizeof(m2aux),  offsetof(m2vol, aux),          MPI_BYTE},
  };
  for (n=0; n<20; ++n) {
    blocklengths[n] = vol_members[n].blocklength;
    displacements[n] = vol_members[n].displacement;
    types[n] = vol_members[n].type;
  }

  MPI_Type_create_struct(20, blocklengths, displacements, types, &vol);
  MPI_Type_commit(&vol);

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

    sizes_sp[n] = G[n+1] > 1 ? Ng1[n+1] : 1;
    sizes_rp[n] = G[n+1] > 1 ? Ng1[n+1] : 1;
    sizes_sm[n] = G[n+1] > 1 ? Ng0[n+1] : 1;
    sizes_rm[n] = G[n+1] > 1 ? Ng0[n+1] : 1;

    start_sp[n] = G[n+1] > 1 ? L[n+1] - 2*Ng1[n+1] : 0;
    start_rp[n] = G[n+1] > 1 ? L[n+1] - 1*Ng1[n+1] : 0;
    start_sm[n] = G[n+1] > 1 ? Ng0[n+1] : 0;
    start_rm[n] = G[n+1] > 1 ? 0 : 0;

    MPI_Type_create_subarray(3, count_sp, sizes_sp, start_sp, order, vol, &sp[n]);
    MPI_Type_create_subarray(3, count_sm, sizes_sm, start_sm, order, vol, &sm[n]);
    MPI_Type_create_subarray(3, count_rp, sizes_rp, start_rp, order, vol, &rp[n]);
    MPI_Type_create_subarray(3, count_rm, sizes_rm, start_rm, order, vol, &rm[n]);
    MPI_Type_commit(&sp[n]);
    MPI_Type_commit(&sm[n]);
    MPI_Type_commit(&rp[n]);
    MPI_Type_commit(&rm[n]);
  }

  for (n=0; n<3; ++n) {
    if (G[n+1] == 1) continue;

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


  MPI_Type_free(&vol);
  for (n=0; n<3; ++n) {
    MPI_Type_free(&sp[n]);
    MPI_Type_free(&sm[n]);
    MPI_Type_free(&rp[n]);
    MPI_Type_free(&rm[n]);
  }

  /* because MPI overwrites aux.m2 */
  for (n=0; n<L[0]; ++n) {
    m2->volumes[n].aux.m2 = m2;
  }
}
