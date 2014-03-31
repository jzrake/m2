#include <math.h>
#include "m2.h"

#define NU_PARAMETER 0.75 /* nu parameter */

void jet_boundary_conditions_cell(m2sim *m2)
{
  int *L = m2->local_grid_size;
  int *G = m2->domain_resolution;
  int *I;
  m2vol *V0;
  int i, j, k;
  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {
	V0 = M2_VOL(i, j, k);
	I = V0->global_index;

        if (I[1] == 0 || I[2] == G[2] - 1) {
          double X[4], x[4];

          m2vol_coordinate_centroid_3d(V0, X);
          m2_to_cartesian(X, x, M2_SPHERICAL);

	  double t0 = 0.2;
	  double t = X[2];
          double n = NU_PARAMETER;
          double r = sqrt(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
          double cost = x[3]/r;
          double rdisk = r * pow(1.0 - cost, 1.0/n);
          double omega = rdisk < 1.0 ? 0.25 : 0.25 * pow(rdisk, -1.5);
          double vkepl = omega * rdisk;
          double vr = t < t0 ? 0.1 : 0.0;
	  double d = t < t0 ? 1.0 : 0.01;
	  double p = t < t0 ? 0.01 : 0.0001;

	  V0->prim.v1 = vr;
          V0->prim.v2 = 0.0;
          V0->prim.v3 = vkepl;
          V0->prim.p = p;
          V0->prim.d = d;
          m2sim_from_primitive(m2,
                               &V0->prim, NULL, NULL,
                               V0 ->volume,
                               V0 ->consA,
                               &V0->aux);
        }
	else if (I[1] == G[1] - 1) {
	  V0->prim.p = 0.0001;
	  V0->prim.d = 0.01;
	  V0->prim.v1 = 0.0;
	  V0->prim.v2 = 0.0;
	  V0->prim.v3 = 0.0;
	  V0->prim.B1 = 0.0;
	  V0->prim.B2 = 0.0;
	  V0->prim.B3 = 0.0;
	  V0->Bflux1A = 0.0;
	  V0->Bflux2A = 0.0;
	  V0->Bflux3A = 0.0;
	  m2sim_from_primitive(m2,
	  		       &V0->prim, NULL, NULL,
	  		       V0 ->volume,
	  		       V0 ->consA,
	  		       &V0->aux);
	}
	if (I[2] == G[2] - 1) {
	  /* initial_data(V0); */
	  /* m2sim_from_primitive(m2, */
	  /* 		       &V0->prim, NULL, NULL, */
	  /* 		       V0 ->volume, */
	  /* 		       V0 ->consA, */
	  /* 		       &V0->aux); */
	}
      }
    }
  }
}
