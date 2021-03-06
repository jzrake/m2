#include <math.h>
#include "m2.h"

#define NU_PARAMETER 0.75 /* nu parameter */

typedef struct jet_model_parameters
{
  double ambient_density;
  double ambient_pressure;
  double nu_parameter;
} jet_mps;


void jet_boundary_conditions_cell(m2sim *m2)
{
  int *L = m2->local_grid_size;
  int *G = m2->domain_resolution;
  int *I;
  m2vol *V0;
  int i, j, k;
  double rinner = 1.0;

  for (i=0; i<L[1]; ++i) {
    for (j=0; j<L[2]; ++j) {
      for (k=0; k<L[3]; ++k) {
	V0 = M2_VOL(i, j, k);
	I = V0->global_index;

        if (I[1] == 0) {
	   /* inner radial boundary: phi-velocity and small radial */
          double X[4], x[4];

          m2vol_coordinate_centroid_3d(V0, X);
          m2_to_cartesian(X, x, M2_SPHERICAL);

          double n = NU_PARAMETER;
	  double t = X[2];
          double r = sqrt(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
          double cost = x[3]/r;
          double rdisk = r * pow(1.0 - cost, 1.0/n);
          double omega = rdisk < rinner ? 0.25 : 0.25 * pow(rdisk/rinner, -1.5);
          double vkepl = omega * rdisk;
          double vz = 0.5 * vkepl;
	  double zh_r =  cos(t);
	  double zh_t = -sin(t);
	  double d = 1.0;
	  double p = 0.01;

          V0->prim.p = p;
          V0->prim.d = d;
	  V0->prim.v1 = vz * zh_r;
          V0->prim.v2 = vz * zh_t;
          V0->prim.v3 = vkepl;
          m2sim_from_primitive(m2,
                               &V0->prim, NULL, NULL,
                               V0 ->volume,
                               V0 ->consA,
                               &V0->aux);
        }
	else if (I[2] == G[2] - 1) {
	  /* equator: phi-velocity and small negative v-theta */
          double X[4], x[4];

          m2vol_coordinate_centroid_3d(V0, X);
          m2_to_cartesian(X, x, M2_SPHERICAL);

          double n = NU_PARAMETER;
	  double t = X[2];
          double r = sqrt(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
          double cost = x[3]/r;
          double rdisk = r * pow(1.0 - cost, 1.0/n);
          double omega = rdisk < rinner ? 0.25 : 0.25 * pow(rdisk/rinner, -1.5);
          double vkepl = omega * rdisk;
          double vz = 0.5 * vkepl;
	  double zh_r =  cos(t);
	  double zh_t = -sin(t);
	  double d = 1.0;
	  double p = 0.01;

          V0->prim.p = p;
          V0->prim.d = d;
	  V0->prim.v1 = vz * zh_r;
          V0->prim.v2 = vz * zh_t;
          V0->prim.v3 = vkepl;
          m2sim_from_primitive(m2,
                               &V0->prim, NULL, NULL,
                               V0 ->volume,
                               V0 ->consA,
                               &V0->aux);
	}
	else if (I[1] == G[1] - 1) {
	  V0->prim.p = 0.01;
	  V0->prim.d = 1.00;
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
      }
    }
  }
}


void jet_add_physical_source_terms(m2vol *V)
{
  double r = m2vol_coordinate_centroid(V, 1);
  double dt = V->x1[0] - V->x0[0];
  double dV = V->volume;
  double vr = V->prim.v1;
  double d0 = V->prim.d;
  double Fr = -1.0 / (r*r);
  V->consA[S11] += dt * dV * d0 * Fr;
  V->consA[TAU] += dt * dV * d0 * Fr * vr;
}



/* ========================================================================== */
/* Jet problem Lua interface                                                  */
/* ========================================================================== */
#include <stddef.h>
#include "struct.h"

static int L_model_parameters__init(lua_State *L)
{
  jet_mps *m = (jet_mps *) lua_struct_checkstruct(L, 1, "jet_model_parameters");
  m->ambient_density = 1.0;
  m->ambient_pressure = 0.01;
  m->nu_parameter = 0.75;
  return 0;
}

int luaopen_m2jet(lua_State *L)
{
#define MD(m) {#m, offsetof(jet_mps, m), LSTRUCT_DOUBLE}
  lua_struct_member_t members[] = {
    MD(ambient_density),
    MD(ambient_pressure),
    MD(nu_parameter),
    {NULL, 0, 0},
  };
#undef MD
  lua_newtable(L); /* jet module */
  lua_struct_t type = lua_struct_newtype(L);
  type.type_name = "jet_model_parameters";
  type.alloc_size = sizeof(jet_mps);
  type.members = members;
  type.__init = L_model_parameters__init;
  lua_struct_register(L, type);

  lua_pushlightuserdata(L, jet_boundary_conditions_cell);
  lua_setfield(L, -2, "jet_boundary_conditions_cell");
  lua_pushlightuserdata(L, jet_add_physical_source_terms);
  lua_setfield(L, -2, "jet_add_physical_source_terms");

  return 1;
}
