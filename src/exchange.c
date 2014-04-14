#include <mpi.h>
#include "m2.h"

void m2sim_synchronize_guard(m2sim *m2)
{
  MPI_Comm cart_comm = *((MPI_Comm *) m2->cart_comm);
  int rankS=0;
  int rankD=0;
  MPI_Cart_shift(cart_comm, 1, +1, &rankS, &rankD);
  //printf("S: %d D: %d\n", rankS, rankD);
}
