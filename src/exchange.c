#include <mpi.h>
#include "m2.h"

void m2sim_synchronize_guard(m2sim *m2)
{
  MPI_Comm cart_comm = *((MPI_Comm *) m2->cart_comm);
  MPI_Status status;
  MPI_Datatype vol;
  MPI_Datatype sp, rp;
  MPI_Datatype sm, rm;

  int order = MPI_ORDER_C;
  int rankS;
  int rankD;
  int tagS = 0;
  int tagD = 0;

  int *Ng0 = m2->number_guard_zones0;
  int *Ng1 = m2->number_guard_zones1;
  int *L = m2->local_grid_size;

  int count_sp[1] = { L[1] };
  int count_rp[1] = { L[1] };
  int count_sm[1] = { L[1] };
  int count_rm[1] = { L[1] };
  int sizes_sp[1] = { Ng1[1] };
  int sizes_rp[1] = { Ng1[1] };
  int sizes_sm[1] = { Ng0[1] };
  int sizes_rm[1] = { Ng0[1] };
  int start_sp[1] = { L[1] - 2*Ng1[1] };
  int start_rp[1] = { L[1] - 1*Ng1[1] };
  int start_sm[1] = { Ng0[1] };
  int start_rm[1] = { 0 };

  MPI_Type_contiguous(sizeof(m2vol), MPI_BYTE, &vol);
  MPI_Type_create_subarray(1, count_sp, sizes_sp, start_sp, order, vol, &sp);
  MPI_Type_create_subarray(1, count_sm, sizes_sm, start_sm, order, vol, &sm);
  MPI_Type_create_subarray(1, count_rp, sizes_rp, start_rp, order, vol, &rp);
  MPI_Type_create_subarray(1, count_rm, sizes_rm, start_rm, order, vol, &rm);
  MPI_Type_commit(&sp);
  MPI_Type_commit(&sm);
  MPI_Type_commit(&rp);
  MPI_Type_commit(&rm);

  /* send to right, receive from left */
  MPI_Cart_shift(cart_comm, 1, +1, &rankS, &rankD);
  MPI_Sendrecv(m2->volumes, 1, sp, rankD, tagD,
	       m2->volumes, 1, rm, rankS, tagS,
	       cart_comm, &status);

  /* send to left, receive from right */
  MPI_Cart_shift(cart_comm, 1, -1, &rankS, &rankD);
  MPI_Sendrecv(m2->volumes, 1, sm, rankD, tagD,
	       m2->volumes, 1, rp, rankS, tagS,
	       cart_comm, &status);

  MPI_Type_free(&vol);
  MPI_Type_free(&sp);
  MPI_Type_free(&sm);
  MPI_Type_free(&rp);
  MPI_Type_free(&rm);
}
