#include <math.h>
#include "m2.h"


void m2_to_cartesian(double x[4], double xcart[4], int geometry)
{
  switch (geometry) {
  case M2_CARTESIAN:
    xcart[1] = x[1];
    xcart[2] = x[2];
    xcart[3] = x[3];
    break;
  case M2_CYLINDRICAL:
    xcart[1] = x[1] * cos(x[2]);
    xcart[2] = x[1] * sin(x[2]);
    xcart[3] = x[3];
    break;
  case M2_SPHERICAL:
    xcart[1] = x[1] * sin(x[2]) * cos(x[3]);
    xcart[2] = x[1] * sin(x[2]) * sin(x[3]);
    xcart[3] = x[1] * cos(x[2]);
    break;
  default:
    MSG(FATAL, "internal error");
    break;
  }
}


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


double m2_line_measure(double x0[4], double x1[4], int geometry, int axis)
{
  switch (geometry) {
  case M2_CARTESIAN:
    switch (axis) {
    case 1: return (x1[1] - x0[1]);
    case 2: return (x1[2] - x0[2]);
    case 3: return (x1[3] - x0[3]);
    default: MSG(FATAL, "internal error"); return 0.0;
    }
    break;
  case M2_CYLINDRICAL:
    switch (axis) {
    case 1: return (x1[1] - x0[1]);
    case 2: return (x1[2] - x0[2]) * x1[1]; /* r df */
    case 3: return (x1[3] - x0[3]);
    default: MSG(FATAL, "internal error"); return 0.0;
    }
    break;
  case M2_SPHERICAL:
    switch (axis) {
    case 1: return (x1[1] - x0[1]); /* dr */
    case 2: return (x1[2] - x0[2]) * x0[1]; /* r dt */
    case 3: return (x1[3] - x0[3]) * x0[1] * sin(x0[2]); /* r sin(t) df */
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
  int d;
  int scale[4] = { 0,
		   m2->coordinate_scaling1,
		   m2->coordinate_scaling2,
		   m2->coordinate_scaling3 };
  x[0] = 0.0;
  for (d=1; d<=3; ++d) {
    switch (scale[d]) {
    case M2_LINEAR:
      x[d] = x0[d] + (x1[d] - x0[d]) / N[d] * (index[d] + 0.5);
      break;
    case M2_LOGARITHMIC:
      x[d] = x0[d] * exp(log(x1[d]/x0[d]) / N[d] * (index[d] + 0.5));
      break;
    default:
      MSG(FATAL, "unknown coordinate scaling");
      break;
    }
  }
}
