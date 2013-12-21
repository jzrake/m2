#include <math.h>
#include "m2.h"


double m2_volume_measure(double x0[4], double x1[4], int geometry)
{
  double dV;
  switch (geometry) {
  case M2_CARTESIAN: dV = (x1[1] - x0[1])*(x1[2] - x0[2])*(x1[3] - x0[3]);
    break;
  case M2_CYLINDRICAL: dV = (x1[1]*x1[1] - x0[1]*x0[1]) * (x1[2] - x0[2]) *
      (x1[3] - x0[3])/2;
    break;
  case M2_SPHERICAL: dV = -1.0/3.0 * (x1[1]*x1[1]*x1[1] - x0[1]*x0[1]*x0[1]) *
      (cos(x1[2]) - cos(x0[2])) * (x1[3] - x0[3]);
    break;
  default:
    MSG(FATAL, "internal error");
    break;
  }
  return dV;
}


double m2_area_measure(double x0[4], double x1[4], int geometry, int axis)
{
  switch (geometry) {
  case M2_CARTESIAN:
    switch (axis) {
    case 1: return (x1[2] - x0[2]) * (x1[3] - x0[3]);
    case 2: return (x1[3] - x0[3]) * (x1[1] - x0[1]);
    case 3: return (x1[1] - x0[1]) * (x1[2] - x0[2]);
    default: MSG(FATAL, "internal error"); return 0.0;
    }
    break;
  case M2_CYLINDRICAL:
    switch (axis) {
    case 1: return (x1[2] - x0[2]) * (x1[3] - x0[3]) * x0[1];
    case 2: return (x1[3] - x0[3]) * (x1[1] - x0[1]);
    case 3: return (x1[1]*x1[1] - x0[1]*x0[1]) * (x1[2] - x0[2])/2;
    default: MSG(FATAL, "internal error"); return 0.0;
    }
    break;
  case M2_SPHERICAL:
    switch (axis) {
    case 1: return -x0[1]*x0[1] * (cos(x1[2]) - cos(x0[2])) * (x1[3] - x0[3]);
    case 2: return (x1[1]*x1[1] - x0[1]*x0[1]) * (x1[3] - x0[3]) * sin(x0[2])/2;
    case 3: return (x1[1]*x1[1] - x0[1]*x0[1]) * (x1[2] - x0[2])/2;
    default: MSG(FATAL, "internal error"); return 0.0;
    }
    break;
  default:
    MSG(FATAL, "internal error");
    return 0.0;
  }
}


void m2sim_index_to_position(m2sim *m2, double index[4], double x[4])
/*
 * Map from the global index space to physical coordinates. Cell volumes are
 * centered at integer values of the index parameter, and faces at the
 * half-index.
 */
{
  double *x0 = m2->domain_extent_lower;
  double *x1 = m2->domain_extent_upper;
  int *N = m2->domain_resolution;
  x[0] = 0.0;
  x[1] = x0[1] + (x1[1] - x0[1]) / N[1] * (index[1] + 0.5);
  x[2] = x0[2] + (x1[2] - x0[2]) / N[2] * (index[2] + 0.5);
  x[3] = x0[3] + (x1[3] - x0[3]) / N[3] * (index[3] + 0.5);
}
