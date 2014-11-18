#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
#include "m2.h"
#include "zulu.h"
#include "struct.h"
#include "buffer.h"

#if M2_HAVE_MPI
#include <mpi.h>
#endif

#define DEBUG_GC 0



/* -----------------------------------------------------------------------------
 *
 * ivec4
 *
 * -------------------------------------------------------------------------- */
static int L_ivec4__tostring(lua_State *L)
{
  int *v = (int *) lua_struct_checkstruct(L, 1, "ivec4");
  lua_pushfstring(L, "[%d %d %d %d]", v[0], v[1], v[2], v[3]);
  return 1;
}
static int L_ivec4__init(lua_State *L)
{
  int n;
  int *v = (int*) lua_struct_checkstruct(L, 1, "ivec4");
  for (n=0; n<4; ++n) v[n] = luaL_optinteger(L, n+2, 0);
  return 0;
}
int register_ivec4(lua_State *L)
{
  lua_struct_member_t members[] = {
    {"0", 0*sizeof(int), LSTRUCT_INT},
    {"1", 1*sizeof(int), LSTRUCT_INT},
    {"2", 2*sizeof(int), LSTRUCT_INT},
    {"3", 3*sizeof(int), LSTRUCT_INT},
    {NULL, 0, 0},
  };
  lua_struct_t type = lua_struct_newtype(L);
  type.type_name = "ivec4";
  type.alloc_size = 4 * sizeof(int);
  type.members = members;
  type.__tostring = L_ivec4__tostring;
  type.__init = L_ivec4__init;
  lua_struct_register(L, type);
  return 0;
}



/* -----------------------------------------------------------------------------
 *
 * ivec8
 *
 * -------------------------------------------------------------------------- */
static int L_ivec8__tostring(lua_State *L)
{
  int *v = (int *) lua_struct_checkstruct(L, 1, "ivec8");
  lua_pushfstring(L, "[%d %d %d %d %d %d %d %d]",
		  v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);
  return 1;
}
static int L_ivec8__init(lua_State *L)
{
  int n;
  int *v = (int*) lua_struct_checkstruct(L, 1, "ivec8");
  for (n=0; n<8; ++n) v[n] = luaL_optinteger(L, n+2, 0);
  return 0;
}
int register_ivec8(lua_State *L)
{
  lua_struct_member_t members[] = {
    {"0", 0*sizeof(int), LSTRUCT_INT},
    {"1", 1*sizeof(int), LSTRUCT_INT},
    {"2", 2*sizeof(int), LSTRUCT_INT},
    {"3", 3*sizeof(int), LSTRUCT_INT},
    {"4", 4*sizeof(int), LSTRUCT_INT},
    {"5", 5*sizeof(int), LSTRUCT_INT},
    {"6", 6*sizeof(int), LSTRUCT_INT},
    {"7", 7*sizeof(int), LSTRUCT_INT},
    {NULL, 0, 0},
  };
  lua_struct_t type = lua_struct_newtype(L);
  type.type_name = "ivec8";
  type.alloc_size = 8 * sizeof(int);
  type.members = members;
  type.__tostring = L_ivec8__tostring;
  type.__init = L_ivec8__init;
  lua_struct_register(L, type);
  return 0;
}




/* -----------------------------------------------------------------------------
 *
 * dvec4
 *
 * -------------------------------------------------------------------------- */
static int L_dvec4__tostring(lua_State *L)
{
  double *v = (double*) lua_struct_checkstruct(L, 1, "dvec4");
  lua_pushfstring(L, "[%f %f %f %f]", v[0], v[1], v[2], v[3]);
  return 1;
}
static int L_dvec4__init(lua_State *L)
{
  int n;
  double *v = (double*) lua_struct_checkstruct(L, 1, "dvec4");
  for (n=0; n<4; ++n) v[n] = luaL_optnumber(L, n+2, 0.0);
  return 0;
}
int register_dvec4(lua_State *L)
{
  lua_struct_member_t members[] = {
    {"0", 0*sizeof(double), LSTRUCT_DOUBLE},
    {"1", 1*sizeof(double), LSTRUCT_DOUBLE},
    {"2", 2*sizeof(double), LSTRUCT_DOUBLE},
    {"3", 3*sizeof(double), LSTRUCT_DOUBLE},
    {NULL, 0, 0},
  };
  lua_struct_t type = lua_struct_newtype(L);
  type.type_name = "dvec4";
  type.alloc_size = 4 * sizeof(double);
  type.members = members;
  type.__tostring = L_dvec4__tostring;
  type.__init = L_dvec4__init;
  lua_struct_register(L, type);
  return 0;
}




/* -----------------------------------------------------------------------------
 *
 * dvec8
 *
 * -------------------------------------------------------------------------- */
static int L_dvec8__tostring(lua_State *L)
{
  double *v = (double*) lua_struct_checkstruct(L, 1, "dvec8");
  lua_pushfstring(L, "[%f %f %f %f %f %f %f %f]", 
		  v[0], v[1], v[2], v[3],
		  v[4], v[5], v[6], v[7]);
  return 1;
}
static int L_dvec8__init(lua_State *L)
{
  int n;
  double *v = (double*) lua_struct_checkstruct(L, 1, "dvec8");
  for (n=0; n<8; ++n) v[n] = luaL_optnumber(L, n+2, 0.0);
  return 0;
}
int register_dvec8(lua_State *L)
{
  lua_struct_member_t members[] = {
    {"0", 0*sizeof(double), LSTRUCT_DOUBLE},
    {"1", 1*sizeof(double), LSTRUCT_DOUBLE},
    {"2", 2*sizeof(double), LSTRUCT_DOUBLE},
    {"3", 3*sizeof(double), LSTRUCT_DOUBLE},
    {"4", 4*sizeof(double), LSTRUCT_DOUBLE},
    {"5", 5*sizeof(double), LSTRUCT_DOUBLE},
    {"6", 6*sizeof(double), LSTRUCT_DOUBLE},
    {"7", 7*sizeof(double), LSTRUCT_DOUBLE},
    {NULL, 0, 0},
  };
  lua_struct_t type = lua_struct_newtype(L);
  type.type_name = "dvec8";
  type.alloc_size = 8 * sizeof(double);
  type.members = members;
  type.__tostring = L_dvec8__tostring;
  type.__init = L_dvec8__init;
  lua_struct_register(L, type);
  return 0;
}



/* -----------------------------------------------------------------------------
 *
 * m2rand
 *
 * -------------------------------------------------------------------------- */
static int L_rand__init(lua_State *L)
{
  jsw_rand_t *R = (jsw_rand_t*) lua_struct_checkstruct(L, 1, "m2rand");
  unsigned long seed = luaL_optlong(L, 2, 0);
  jsw_seed(R, seed);
  return 0;
}
static int L_rand_seed(lua_State *L)
{
  jsw_rand_t *R = (jsw_rand_t*) lua_struct_checkstruct(L, 1, "m2rand");
  unsigned long seed = luaL_checklong(L, 2);
  jsw_seed(R, seed);
  return 0;
}
static int L_rand_rand(lua_State *L)
{
  jsw_rand_t *R = (jsw_rand_t*) lua_struct_checkstruct(L, 1, "m2rand");
  unsigned long res = jsw_rand(R);
  lua_pushnumber(L, res);
  return 1;
}
int register_rand(lua_State *L)
{
  int n;
  char member_names[JSW_STATE_SIZE][8];
  lua_struct_member_t members[JSW_STATE_SIZE+2];

  for (n=0; n<JSW_STATE_SIZE; ++n) {
    snprintf(member_names[n], 8, "%d", n);
  }

  for (n=0; n<JSW_STATE_SIZE; ++n) {
    lua_struct_member_t M = {member_names[n], offsetof(jsw_rand_t, x[n]),
			     LSTRUCT_ULONG};
    members[n] = M;
  }

  members[624].member_name = "next";
  members[624].offset = offsetof(jsw_rand_t, next);
  members[624].data_type = LSTRUCT_INT;
  members[625].member_name = NULL;

  lua_struct_method_t methods[] = {
    {"seed", L_rand_seed},
    {"rand", L_rand_rand},
    {NULL, NULL},
  } ;
  lua_struct_t type = lua_struct_newtype(L);
  type.type_name = "m2rand";
  type.alloc_size = sizeof(jsw_rand_t);
  type.members = members;
  type.methods = methods;
  type.__init = L_rand__init;
  lua_struct_register(L, type);
  return 0;
}



/* -----------------------------------------------------------------------------
 *
 * m2prim
 *
 * -------------------------------------------------------------------------- */
static int L_m2prim__gc(lua_State *L)
{
  if (DEBUG_GC) printf("killing m2prim\n");
  return 0;
}
static int L_m2prim__tostring(lua_State *L)
{
  m2prim *P = (m2prim *) lua_struct_checkstruct(L, 1, "m2prim");
  lua_pushfstring(L, "<m2prim> {"
		  "v1=%f v2=%f v3=%f "
		  "B1=%f B2=%f B3=%f "
		  "p=%f d=%f}",
		  P->v1, P->v2, P->v3,
		  P->B1, P->B2, P->B3,
		  P->p, P->d);
  return 1;
}
int register_m2prim(lua_State *L)
{
  lua_struct_member_t members[] = {
    {"v1", offsetof(m2prim, v1), LSTRUCT_DOUBLE},
    {"v2", offsetof(m2prim, v2), LSTRUCT_DOUBLE},
    {"v3", offsetof(m2prim, v3), LSTRUCT_DOUBLE},
    {"B1", offsetof(m2prim, B1), LSTRUCT_DOUBLE},
    {"B2", offsetof(m2prim, B2), LSTRUCT_DOUBLE},
    {"B3", offsetof(m2prim, B3), LSTRUCT_DOUBLE},
    {"p", offsetof(m2prim, p), LSTRUCT_DOUBLE},
    {"d", offsetof(m2prim, d), LSTRUCT_DOUBLE},
    {NULL, 0, 0},
  };
  lua_struct_t type = lua_struct_newtype(L);
  type.type_name = "m2prim";
  type.alloc_size = sizeof(m2prim);
  type.members = members;
  type.methods = NULL;
  type.__tostring = L_m2prim__tostring;
  type.__gc = L_m2prim__gc;
  lua_struct_register(L, type);
  return 0;
}




/* -----------------------------------------------------------------------------
 *
 * m2aux
 *
 * -------------------------------------------------------------------------- */
static int L_m2aux__gc(lua_State *L)
{
  if (DEBUG_GC) printf("killing m2aux\n");
  return 0;
}
static int L_m2aux_fluxes(lua_State *L)
{
  m2aux *aux = (m2aux *) lua_struct_checkstruct(L, 1, "m2aux");
  double *n = (double *) lua_struct_checkstruct(L, 2, "dvec4");
  double *F = (double *) lua_struct_new(L, "dvec8");
  m2aux_fluxes(aux, n, F);
  return 1;
}
static int L_m2aux_eigenvalues(lua_State *L)
{
  m2aux *aux = (m2aux *) lua_struct_checkstruct(L, 1, "m2aux");
  double *n = (double *) lua_struct_checkstruct(L, 2, "dvec4");
  double *E = (double *) lua_struct_new(L, "dvec8");
  m2aux_eigenvalues(aux, n, E);
  return 1;
}
int register_m2aux(lua_State *L)
{
#define MD(m) {#m, offsetof(m2aux, m), LSTRUCT_DOUBLE}
#define MS(m,n) {#m, offsetof(m2aux, m), LSTRUCT_STRUCT, n}
  lua_struct_member_t members[] = {
    MS(velocity_four_vector, "dvec4"),
    MS(magnetic_four_vector, "dvec4"),
    MS(momentum_density, "dvec4"),
    MD(comoving_mass_density),
    MD(gas_pressure),
    MD(magnetic_pressure),
    {NULL, 0, 0},
  };
  lua_struct_method_t methods[] = {
    {"fluxes", L_m2aux_fluxes},
    {"eigenvalues", L_m2aux_eigenvalues},
    {NULL, NULL},
  } ;
#undef MD
#undef MS
  lua_struct_t type = lua_struct_newtype(L);
  type.type_name = "m2aux";
  type.alloc_size = sizeof(m2aux);
  type.members = members;
  type.methods = methods;
  type.__gc = L_m2aux__gc;
  lua_struct_register(L, type);
  return 0;
}




/* -----------------------------------------------------------------------------
 *
 * m2reductions
 *
 * -------------------------------------------------------------------------- */
int register_m2reductions(lua_State *L)
{
#define MD(m) {#m, offsetof(m2reductions, m), LSTRUCT_DOUBLE}
  lua_struct_member_t members[] = {
    MD(total_volume),
    MD(total_mass),
    MD(total_energy),
    MD(total_kinetic_energy),
    MD(total_internal_energy),
    MD(total_magnetic_energy),
    {NULL, 0, 0},
  };
  lua_struct_method_t methods[] = {
    {NULL, NULL},
  } ;
#undef MD
#undef MS
  lua_struct_t type = lua_struct_newtype(L);
  type.type_name = "m2reductions";
  type.alloc_size = sizeof(m2reductions);
  type.members = members;
  type.methods = methods;
  lua_struct_register(L, type);
  return 0;
}




/* -----------------------------------------------------------------------------
 *
 * m2vol
 *
 * -------------------------------------------------------------------------- */
static int L_m2vol_next(lua_State *L)
{
  m2vol *V = (m2vol *) lua_struct_checkstruct(L, 1, "m2vol");
  int n = luaL_checknumber(L, 2);
  if (n < V->m2->local_grid_size[0] - 1) {
    lua_struct_pushparent(L, 1);
    lua_pushnumber(L, n + 1);
    lua_struct_newref(L, "m2vol", &V->m2->volumes[n+1], -2);
    lua_remove(L, -3);
    return 2;
  }
  else {
    return 0;
  }
}
static int L_m2vol_neighbor(lua_State *L)
{
  m2vol *V = (m2vol *) lua_struct_checkstruct(L, 1, "m2vol");
  int axis = luaL_checkinteger(L, 2);
  int dist = luaL_checkinteger(L, 3);
  m2vol *V1 = m2vol_neighbor(V, axis, dist);
  if (V1) {
    lua_struct_pushparent(L, 1);
    lua_struct_newref(L, "m2vol", V1, -1);
    lua_remove(L, -2);
    return 1;
  }
  else {
    return 0;
  }
  return 0;
}
static int L_m2vol_from_primitive(lua_State *L)
{
  m2vol *V = (m2vol *) lua_struct_checkstruct(L, 1, "m2vol");
  m2vol_from_primitive(V);
  return 0;
}
static int L_m2vol__gc(lua_State *L)
{
  if (DEBUG_GC) printf("killing m2vol\n");
  return 0;
}
int register_m2vol(lua_State *L)
{
  lua_struct_member_t members[] = {
    {"global_index", offsetof(m2vol, global_index), LSTRUCT_STRUCT, "ivec4"},
    {"x0", offsetof(m2vol, x0), LSTRUCT_STRUCT, "dvec4"},
    {"x1", offsetof(m2vol, x1), LSTRUCT_STRUCT, "dvec4"},
    {"prim", offsetof(m2vol, prim), LSTRUCT_STRUCT, "m2prim"},
    {"aux", offsetof(m2vol, aux), LSTRUCT_STRUCT, "m2aux"},
    {"flux1", offsetof(m2vol, flux1), LSTRUCT_STRUCT, "dvec8"},
    {"flux2", offsetof(m2vol, flux2), LSTRUCT_STRUCT, "dvec8"},
    {"flux3", offsetof(m2vol, flux3), LSTRUCT_STRUCT, "dvec8"},
    {"volume", offsetof(m2vol, volume), LSTRUCT_DOUBLE},
    {"area1", offsetof(m2vol, area1), LSTRUCT_DOUBLE},
    {"area2", offsetof(m2vol, area2), LSTRUCT_DOUBLE},
    {"area3", offsetof(m2vol, area3), LSTRUCT_DOUBLE},
    {"line1", offsetof(m2vol, line1), LSTRUCT_DOUBLE},
    {"line2", offsetof(m2vol, line2), LSTRUCT_DOUBLE},
    {"line3", offsetof(m2vol, line3), LSTRUCT_DOUBLE},
    {"emf1", offsetof(m2vol, emf1), LSTRUCT_DOUBLE},
    {"emf2", offsetof(m2vol, emf2), LSTRUCT_DOUBLE},
    {"emf3", offsetof(m2vol, emf3), LSTRUCT_DOUBLE},
    {"Bflux1A", offsetof(m2vol, Bflux1A), LSTRUCT_DOUBLE},
    {"Bflux2A", offsetof(m2vol, Bflux2A), LSTRUCT_DOUBLE},
    {"Bflux3A", offsetof(m2vol, Bflux3A), LSTRUCT_DOUBLE},
    {"zone_type", offsetof(m2vol, zone_type), LSTRUCT_CHAR},
    {NULL, 0, 0},
  };
  lua_struct_method_t methods[] = {
    {"next", L_m2vol_next},
    {"neighbor", L_m2vol_neighbor},
    {"from_primitive", L_m2vol_from_primitive},
    {NULL, NULL},
  };
  lua_struct_t type = lua_struct_newtype(L);
  type.type_name = "m2vol";
  type.alloc_size = sizeof(m2vol);
  type.members = members;
  type.methods = methods;
  type.__gc = L_m2vol__gc;
  lua_struct_register(L, type);
  return 0;
}



/* -----------------------------------------------------------------------------
 *
 * m2stirring
 *
 * -------------------------------------------------------------------------- */
static int L_m2stirring_next_realization(lua_State *L)
{
  m2stirring *S = (m2stirring *) lua_struct_checkstruct(L, 1, "m2stirring");
  m2stirring_next_realization(S);
  return 0;
}
int register_m2stirring(lua_State *L)
{
#define MI(m)   {#m, offsetof(m2stirring, m), LSTRUCT_INT}
#define MD(m)   {#m, offsetof(m2stirring, m), LSTRUCT_DOUBLE}
#define MS(m,n) {#m, offsetof(m2stirring, m), LSTRUCT_STRUCT, n}
  lua_struct_member_t members[] = {
    MI(stirring_type),
    MI(num_waves),
    MD(wavelength),
    MD(amplitude),
    MS(random, "m2rand"),
    {NULL, 0, 0},
  };
  lua_struct_method_t methods[] = {
    {"next_realization", L_m2stirring_next_realization},
    {NULL, NULL},
  };
#undef MI
#undef MD
#undef MS
  lua_struct_t type = lua_struct_newtype(L);
  type.type_name = "m2stirring";
  type.alloc_size = sizeof(m2stirring);
  type.members = members;
  type.methods = methods;
  lua_struct_register(L, type);
  return 0;
}



/* -----------------------------------------------------------------------------
 *
 * m2status
 *
 * -------------------------------------------------------------------------- */
static int L_m2sim_status__gc(lua_State *L)
{
  if (DEBUG_GC) printf("killing m2sim_status\n");
  return 0;
}
static int L_m2sim_status_get_message(lua_State *L)
{
  m2sim_status *status = (m2sim_status *)
    lua_struct_checkstruct(L, 1, "m2sim_status");
  lua_pushstring(L, status->message);
  return 1;
}
int register_m2status(lua_State *L)
{
#define MI(m)   {#m, offsetof(m2sim_status, m), LSTRUCT_INT}
#define MD(m)   {#m, offsetof(m2sim_status, m), LSTRUCT_DOUBLE}
#define MS(m,n) {#m, offsetof(m2sim_status, m), LSTRUCT_STRUCT, n}
  lua_struct_member_t members[] = {
    MI(iteration_number),
    MD(time_simulation),
    MD(time_step),
    MD(time_last_checkpoint_hdf5),
    MD(time_last_checkpoint_tpl),
    MD(time_last_analysis),
    MI(checkpoint_number_hdf5),
    MI(checkpoint_number_tpl),
    MI(error_code),
    MI(drive_attempt),
    MS(num_failures, "ivec8"),
    {"message", offsetof(m2sim_status, message), LSTRUCT_STRING, NULL, 1024},
    {NULL, 0, 0},
  };
  lua_struct_method_t methods[] = {
    {"get_message", L_m2sim_status_get_message},
    {NULL, NULL},
  };
#undef MI
#undef MD
#undef MS
  lua_struct_t type = lua_struct_newtype(L);
  type.type_name = "m2sim_status";
  type.alloc_size = sizeof(m2sim_status);
  type.members = members;
  type.methods = methods;
  type.__gc = L_m2sim_status__gc,
  lua_struct_register(L, type);
  return 0;
}




/* -----------------------------------------------------------------------------
 *
 * m2sim (unique to Lua API)
 *
 * -------------------------------------------------------------------------- */
static int L_m2sim__init(lua_State *L)
{
  m2sim *m2 = (m2sim *) lua_struct_checkstruct(L, 1, "m2sim");
  m2sim_new(m2);
  return 0;
}
static int L_m2sim__gc(lua_State *L)
{
  m2sim *m2 = (m2sim *) lua_struct_checkstruct(L, 1, "m2sim");
  if (DEBUG_GC) printf("killing m2sim\n");
  m2sim_del(m2);
  return 0;
}
static int L_m2sim_get_volume(lua_State *Ls)
{
  m2sim *m2 = (m2sim *) lua_struct_checkstruct(Ls, 1, "m2sim");
  int i = luaL_optinteger(Ls, 2, 0);
  int j = luaL_optinteger(Ls, 3, 0);
  int k = luaL_optinteger(Ls, 4, 0);
  int *L = m2->local_grid_size;
  m2vol *V = M2_VOL(i, j, k);
  if (V == NULL) {
    return 0;
  }
  else {
    lua_struct_newref(Ls, "m2vol", V, 1);
    return 1;
  }
}
static int L_m2sim_volumes(lua_State *L)
{
  m2sim *m2 = (m2sim *) lua_struct_checkstruct(L, 1, "m2sim");
  lua_pushcfunction(L, L_m2vol_next);
  lua_struct_newref(L, "m2vol", m2->volumes, 1);
  lua_pushnumber(L, -1);
  return 3;
}
static int L_m2sim_get_volume_data(lua_State *L)
{
  const char *const lst[] = {"v1", "v2", "v3",
			     "B1", "B2", "B3", "p", "d", NULL};
  m2sim *m2 = (m2sim *) lua_struct_checkstruct(L, 1, "m2sim");
  int offset = luaL_checkoption(L, 2, NULL, lst);
  int n, N = m2->local_grid_size[0];
  double *D = (double *) buf_new_buffer(L, NULL, N * sizeof(double));
  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    D[n] = *((double*) &m2->mesh.cells[n].prim + offset);
  }
  return 1;
}
static int L_m2sim_set_volume_data(lua_State *L)
{
  const char *const lst[] = {"v1", "v2", "v3",
			     "B1", "B2", "B3", "p", "d", NULL};
  m2sim *m2 = (m2sim *) lua_struct_checkstruct(L, 1, "m2sim");
  int offset = luaL_checkoption(L, 2, NULL, lst);
  int n, N = m2->local_grid_size[0];
  double *D = (double *) buf_check_buffer(L, 3, N * sizeof(double));
  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    *((double*) &m2->mesh.cells[n].prim + offset) = D[n];
  }
  return 0;
}
static int L_m2sim_get_face_data(lua_State *L)
{
  m2sim *m2 = (m2sim *) lua_struct_checkstruct(L, 1, "m2sim");
  int axis = luaL_checkunsigned(L, 2);
  int n, N = m2->local_grid_size[0];
  double *D = (double *) buf_new_buffer(L, NULL, N * sizeof(double));
  if (axis != 1 && axis != 2 && axis != 3) {
    luaL_error(L, "argument 2 must be 1, 2, or 3");
  }
  for (n=0; n<m2->mesh.faces_shape[axis][n]; ++n) {
    D[n] = *((double*) &m2->mesh.faces[axis][n].BfluxA);
  }
  return 1;
}
static int L_m2sim_set_face_data(lua_State *L)
{
  m2sim *m2 = (m2sim *) lua_struct_checkstruct(L, 1, "m2sim");
  int axis = luaL_checkunsigned(L, 2);
  int n, N = m2->local_grid_size[0];
  double *D = (double *) buf_new_buffer(L, NULL, N * sizeof(double));
  if (axis != 1 && axis != 2 && axis != 3) {
    luaL_error(L, "argument 2 must be 1, 2, or 3");
  }
  for (n=0; n<m2->mesh.faces_shape[axis][n]; ++n) {
    *((double*) &m2->mesh.faces[axis][n].BfluxA) = D[n];
  }
  return 1;
}
struct user_call
{
  int cell;
  int face;
  int edge;
  int bc_flux1;
  int bc_flux2;
  int bc_flux3;
  int bc_emf1;
  int bc_emf2;
  int bc_emf3;
  int bc_cell;
  int srcterm;
  lua_State *L;
} ;
#define ID_CALLBACK(index,Nq)						\
  static void id_##index(m2sim *m2, double X[4], double N[4], double *Y) \
  {									\
    struct user_call *u = (struct user_call*) m2->user_struct;		\
    lua_State *L = u->L;						\
    int q;								\
    lua_pushvalue(L, u->index);						\
    lua_struct_newref(L, "dvec4", X, LSTRUCT_ORPHAN);			\
    lua_struct_newref(L, "dvec4", N, LSTRUCT_ORPHAN);			\
    lua_call(L, 2, 1);							\
    if (!lua_istable(L, -1)) {						\
      luaL_error(L, "initial data callback '" #index "' must return a table"); \
    }									\
    for (q=0; q<Nq; ++q) {						\
      lua_rawgeti(L, -1, q+1);						\
      if (!lua_isnumber(L, -1)) {					\
	luaL_error(L, "initial data callback '"#index "' at index %d "	\
		   "must be a number", q);				\
      }									\
      Y[q] = lua_tonumber(L, -1);					\
      lua_pop(L, 1);							\
    }									\
    lua_pop(L, 1);							\
  }									\

#define BC_CALLBACK(m)							\
  static void bc_##m(m2vol *V)						\
  {									\
    struct user_call *u = (struct user_call*) V->m2->user_struct;	\
    lua_State *L = u->L;						\
    if (u->bc_##m != -1) {						\
      lua_pushvalue(L, u->bc_##m);					\
      lua_struct_newref(L, "m2vol", V, 1);				\
      lua_call(L, 1, 0);						\
    }									\
  }									\

static void bc_cell(m2sim *m2)
{
  struct user_call *u = (struct user_call*) m2->user_struct;
  lua_State *L = u->L;
  if (u->bc_cell != -1) {
    lua_pushvalue(L, u->bc_cell);
    lua_pushvalue(L, 1); /* m2sim on the stack */
    lua_call(L, 1, 0);
  }
}

static void cb_srcterm(m2vol *V)
{
  struct user_call *u = (struct user_call*) V->m2->user_struct;
  lua_State *L = u->L;
  if (u->srcterm != -1) {
    lua_pushvalue(L, u->srcterm);
    lua_struct_newref(L, "m2vol", V, 1);
    lua_call(L, 1, 0);
  }
}

ID_CALLBACK(cell, 5)
ID_CALLBACK(face, 1)
ID_CALLBACK(edge, 1)
BC_CALLBACK(flux1)
BC_CALLBACK(flux2)
BC_CALLBACK(flux3)
BC_CALLBACK(emf1)
BC_CALLBACK(emf2)
BC_CALLBACK(emf3)

#undef ID_CALLBACK
#undef BC_CALLBACK

static int L_m2sim_run_initial_data(lua_State *L)
{
  m2sim *m2 = (m2sim *) lua_struct_checkstruct(L, 1, "m2sim");
  m2crd_function orig_cell = m2->initial_data_cell;
  m2crd_function orig_face = m2->initial_data_edge;
  m2crd_function orig_edge = m2->initial_data_edge;

  struct user_call *u = (struct user_call*) malloc(sizeof(struct user_call));
  u->cell = lua_type(L, 2) == LUA_TFUNCTION ? 2 : -1;
  u->face = lua_type(L, 3) == LUA_TFUNCTION ? 3 : -1;
  u->edge = lua_type(L, 4) == LUA_TFUNCTION ? 4 : -1;
  u->L = L;
  m2->user_struct = u;

  m2->initial_data_cell = u->cell != -1 ? id_cell : NULL;
  m2->initial_data_face = u->face != -1 ? id_face : NULL;
  m2->initial_data_edge = u->edge != -1 ? id_edge : NULL;

  m2sim_run_initial_data(m2);
  m2->initial_data_cell = orig_cell;
  m2->initial_data_face = orig_face;
  m2->initial_data_edge = orig_edge;
  free(u);
  m2->user_struct = NULL;
  return 0;
}
static int _drive_or_visualize(lua_State *L, char mode)
{
  m2sim *m2 = (m2sim *) lua_struct_checkstruct(L, 1, "m2sim");
  m2vol_operator tbc_flux1 = m2->boundary_conditions_flux1;
  m2vol_operator tbc_flux2 = m2->boundary_conditions_flux2;
  m2vol_operator tbc_flux3 = m2->boundary_conditions_flux3;
  m2vol_operator tbc_emf1 = m2->boundary_conditions_emf1;
  m2vol_operator tbc_emf2 = m2->boundary_conditions_emf2;
  m2vol_operator tbc_emf3 = m2->boundary_conditions_emf3;
  m2sim_operator tbc_cell = m2->boundary_conditions_cell;
  m2vol_operator tsrctrms = m2->add_physical_source_terms;

  struct user_call *u = (struct user_call*) malloc(sizeof(struct user_call));
  u->bc_flux1 = -1;
  u->bc_flux2 = -1;
  u->bc_flux3 = -1;
  u->bc_emf1 = -1;
  u->bc_emf2 = -1;
  u->bc_emf3 = -1;
  u->bc_cell = -1;
  u->srcterm = -1;
  u->L = L;
  m2->user_struct = u;

#define CHECKF(m)						\
  do {								\
    lua_struct_pushmember(L, 1, "boundary_conditions_"#m);	\
    if (lua_isnil(L, -1)) {					\
								\
    }								\
    else if (lua_type(L, -1) == LUA_TFUNCTION) {		\
      u->bc_##m = lua_gettop(L);				\
	m2->boundary_conditions_##m = bc_##m;			\
    }								\
    else if (lua_type(L, -1) == LUA_TLIGHTUSERDATA) {		\
      u->bc_##m = -1;						\
	m2->boundary_conditions_##m = lua_touserdata(L, -1);	\
    }								\
    else {							\
      luaL_error(L,						\
		 "callback boundary_conditions_"#m" not set "	\
		 "to a function value");			\
    }								\
  } while (0)							\

  CHECKF(flux1);
  CHECKF(flux2);
  CHECKF(flux3);
  CHECKF(emf1);
  CHECKF(emf2);
  CHECKF(emf3);
  CHECKF(cell);
#undef CHECKF

  lua_struct_pushmember(L, 1, "add_physical_source_terms");
  if (lua_isnil(L, -1)) { }
  else if (lua_type(L, -1) == LUA_TFUNCTION) {
    u->srcterm = lua_gettop(L);
    m2->add_physical_source_terms = cb_srcterm;
  }
  else if (lua_type(L, -1) == LUA_TLIGHTUSERDATA) {
    u->srcterm = -1;
    m2->add_physical_source_terms = lua_touserdata(L, -1);
  }
  else {
    luaL_error(L,
	       "callback add_physical_source_terms not set "
	       "to a function value");
  }

  if (mode == 'd') {
    m2sim_drive(m2);
  }
  else if (mode == 'v') {
    m2sim_visualize(m2, 0, NULL);
  }

  free(u);
  m2->user_struct = NULL;
  m2->boundary_conditions_flux1 = tbc_flux1;
  m2->boundary_conditions_flux2 = tbc_flux2;
  m2->boundary_conditions_flux3 = tbc_flux3;
  m2->boundary_conditions_emf1 = tbc_emf1;
  m2->boundary_conditions_emf2 = tbc_emf2;
  m2->boundary_conditions_emf3 = tbc_emf3;
  m2->boundary_conditions_cell = tbc_cell;
  m2->add_physical_source_terms = tsrctrms;
  return 0;
}

static int L_m2sim_drive(lua_State *L)
{
  return _drive_or_visualize(L, 'd');
}

static int L_m2sim_visualize(lua_State *L)
{
  return _drive_or_visualize(L, 'v');
}
static int L_m2sim_get_reductions(lua_State *L)
{
  int n;
  m2sim *m2 = (m2sim *) lua_struct_checkstruct(L, 1, "m2sim");
  m2reductions *R = (m2reductions *) lua_struct_new(L, "m2reductions");
  R->total_volume          = 0.0;
  R->total_mass            = 0.0;
  R->total_energy          = 0.0;
  R->total_kinetic_energy  = 0.0;
  R->total_internal_energy = 0.0;
  R->total_magnetic_energy = 0.0;
  for (n=0; n<m2->mesh.cells_shape[0]; ++n) {
    struct mesh_cell *C = m2->mesh.cells + n;
    double v = C->volume;
    if (C->zone_type != M2_ZONE_TYPE_FULL) continue; /* TODO: ensure zone_type is set */
    R->total_volume          += v;
    R->total_mass            += v*m2aux_get(&C->aux, M2_OBSERVER_MASS_DENSITY);
    R->total_energy          += v*m2aux_get(&C->aux, M2_TOTAL_ENERGY_DENSITY);
    R->total_kinetic_energy  += v*m2aux_get(&C->aux, M2_KINETIC_ENERGY_DENSITY);
    R->total_internal_energy += v*m2aux_get(&C->aux, M2_INTERNAL_ENERGY_DENSITY);
    R->total_magnetic_energy += v*m2aux_get(&C->aux, M2_MAGNETIC_ENERGY_DENSITY);
  }
#if M2_HAVE_MPI
    MPI_Comm cart_comm = *((MPI_Comm *) m2->cart_comm);
    MPI_Allreduce(MPI_IN_PLACE, (double*)R, 6, MPI_DOUBLE, MPI_SUM, cart_comm);
#endif
  return 1;
}




/* -----------------------------------------------------------------------------
 *
 * m2sim (direct calls to C functions)
 *
 * -------------------------------------------------------------------------- */
#define WRAP0(nm)							\
  static int L_m2sim_##nm(lua_State *L)					\
  {									\
    m2sim *m2 = (m2sim *) lua_struct_checkstruct(L, 1, "m2sim");	\
    m2sim_##nm(m2);							\
      return 0;								\
  }									\

#define WRAP1(nm)							\
  static int L_m2sim_##nm(lua_State *L)					\
  {									\
    m2sim *m2 = (m2sim *) lua_struct_checkstruct(L, 1, "m2sim");	\
    double x = luaL_checknumber(L, 2);					\
    m2sim_##nm(m2, x);							\
      return 0;								\
  }

#define WRAP2(nm)							\
  static int L_m2sim_##nm(lua_State *L)					\
  {									\
    m2sim *m2 = (m2sim *) lua_struct_checkstruct(L, 1, "m2sim");	\
    lua_pushnumber(L, m2sim_##nm(m2));					\
    return 1;								\
  }									\

WRAP0(initialize)
WRAP0(calculate_gradient)
WRAP0(calculate_flux)
WRAP0(calculate_emf)
WRAP1(exchange_flux)
WRAP1(add_source_terms)
WRAP0(cache_conserved)
WRAP1(average_runge_kutta)
WRAP0(magnetic_flux_to_cell_center)
WRAP0(enforce_user_constraint)
WRAP0(print)
WRAP2(memory_usage)
WRAP2(from_conserved_all)
WRAP2(from_primitive_all)
WRAP2(minimum_courant_time)

#undef WRAP0
#undef WRAP1
#undef WRAP2

int register_m2sim(lua_State *L)
{
#define MD(m)   {#m, offsetof(m2sim, m), LSTRUCT_DOUBLE}
#define MI(m)   {#m, offsetof(m2sim, m), LSTRUCT_INT}
#define MS(m,n) {#m, offsetof(m2sim, m), LSTRUCT_STRUCT, n}
#define MO(m)   {#m, offsetof(m2sim, m), LSTRUCT_OBJECT}
#define MP(m)   {#m, offsetof(m2sim, m), LSTRUCT_POINTER}
  lua_struct_member_t members[] = {
    MS(domain_extent_lower, "dvec4"),
    MS(domain_extent_upper, "dvec4"),
    MS(domain_resolution, "ivec4"),
    MS(local_grid_size, "ivec4"),
    MS(local_grid_start, "ivec4"),
    MS(number_guard_zones0, "ivec4"),
    MS(number_guard_zones1, "ivec4"),
    MS(boundary_condition_lower, "ivec4"),
    MS(boundary_condition_upper, "ivec4"),
    MS(stirring, "m2stirring"),
    MS(status, "m2sim_status"),
    MI(coordinate_scaling1),
    MI(coordinate_scaling2),
    MI(coordinate_scaling3),
    MI(geometry),
    MI(relativistic),
    MI(magnetized),
    MI(rk_order),
    MI(gradient_fields),
    MI(gradient_estimation_method),
    MI(riemann_solver),
    MI(quartic_solver),
    MI(suppress_extrapolation_at_unhealthy_zones),
    MD(pressure_floor),
    MD(plm_parameter),
    MD(cfl_parameter),
    MD(gamma_law_index),
    MD(cadence_checkpoint_hdf5),
    MD(cadence_checkpoint_tpl),
    MD(cadence_analysis),
    MO(boundary_conditions_cell),
    MO(boundary_conditions_flux1),
    MO(boundary_conditions_flux2),
    MO(boundary_conditions_flux3),
    MO(boundary_conditions_emf1),
    MO(boundary_conditions_emf2),
    MO(boundary_conditions_emf3),
    MO(add_physical_source_terms),
    MP(cart_comm),
    {NULL, 0, 0},
  };
#undef MD
#undef MI
#undef MS

#define METHOD(m) {#m, L_m2sim_##m}
  lua_struct_method_t methods[] = {
    METHOD(volumes),
    METHOD(get_volume),
    METHOD(get_volume_data),
    METHOD(set_volume_data),
    METHOD(get_face_data),
    METHOD(set_face_data),
    METHOD(initialize),
    METHOD(calculate_gradient),
    METHOD(calculate_flux),
    METHOD(calculate_emf),
    METHOD(exchange_flux),
    METHOD(add_source_terms),
    METHOD(cache_conserved),
    METHOD(average_runge_kutta),
    METHOD(magnetic_flux_to_cell_center),
    METHOD(enforce_user_constraint),
    METHOD(print),
    METHOD(run_initial_data),
    METHOD(memory_usage),
    METHOD(from_conserved_all),
    METHOD(from_primitive_all),
    METHOD(minimum_courant_time),
    METHOD(drive),
    METHOD(visualize),
    METHOD(get_reductions),
    {NULL, NULL},
  };
#undef METHOD

  lua_struct_t type = lua_struct_newtype(L);
  type.type_name = "m2sim";
  type.alloc_size = sizeof(m2sim),
  type.members = members;
  type.methods = methods;
  type.__init = L_m2sim__init;
  type.__gc = L_m2sim__gc;
  lua_struct_register(L, type);
  return 0;
}



/* -----------------------------------------------------------------------------
 *
 * m2 utilities
 *
 * -------------------------------------------------------------------------- */
static int L_m2_self_test(lua_State *L)
{
  m2_self_test();
  return 0;
}
static int L_m2_force_free_vector_potential(lua_State *L)
{
  double **x = (double **) luaL_checkudata(L, 1, "dvec4");
  double **n = (double **) luaL_checkudata(L, 2, "dvec4");
  int M = luaL_checkint(L, 3);
  double A = m2_force_free_vector_potential(*x, *n, M);
  lua_pushnumber(L, A);
  return 1;
}
static int L_m2_force_free_magnetic_field(lua_State *L)
{
  double **x = (double **) luaL_checkudata(L, 1, "dvec4");
  double **n = (double **) luaL_checkudata(L, 2, "dvec4");
  int M = luaL_checkint(L, 3);
  double A = m2_force_free_magnetic_field(*x, *n, M);
  lua_pushnumber(L, A);
  return 1;
}
static int L_m2_magnetic_rope_vector_potential(lua_State *L)
{
  double **x = (double **) luaL_checkudata(L, 1, "dvec4");
  double **n = (double **) luaL_checkudata(L, 2, "dvec4");
  double A = m2_magnetic_rope_vector_potential(*x, *n);
  lua_pushnumber(L, A);
  return 1;
}




int luaopen_m2lib(lua_State *L)
{
  lua_newtable(L);
  register_m2sim(L);
  register_m2vol(L);
  register_m2aux(L);
  register_m2prim(L);
  register_m2stirring(L);
  register_m2status(L);
  register_m2reductions(L);
  register_ivec4(L);
  register_ivec8(L);
  register_dvec4(L);
  register_dvec8(L);
  register_rand(L);

  lua_pushstring(L, M2_GIT_SHA);
  lua_setfield(L, -2, "M2_GIT_SHA");
  lua_pushstring(L, __DATE__" "__TIME__);
  lua_setfield(L, -2, "M2_BUILD_DATE");
  lua_pushcfunction(L, L_m2_self_test);
  lua_setfield(L, -2, "self_test");
  lua_pushcfunction(L, L_m2_force_free_vector_potential);
  lua_setfield(L, -2, "force_free_vector_potential");
  lua_pushcfunction(L, L_m2_force_free_magnetic_field);
  lua_setfield(L, -2, "force_free_magnetic_field");
  lua_pushcfunction(L, L_m2_magnetic_rope_vector_potential);
  lua_setfield(L, -2, "magnetic_rope_vector_potential");

  lua_pushnumber(L, M2_PI);
  lua_setfield(L, -2, "M2_PI");

#define CONSTANT(m) lua_pushnumber(L, m); lua_setfield(L, -2, #m)
  CONSTANT(M2_CARTESIAN);
  CONSTANT(M2_CYLINDRICAL);
  CONSTANT(M2_SPHERICAL);
  CONSTANT(M2_PARABOLIC);
  CONSTANT(M2_LINEAR);
  CONSTANT(M2_LOGARITHMIC);
  CONSTANT(M2_CONSERVED);
  CONSTANT(M2_AUXILIARY);
  CONSTANT(M2_PRIMITIVE);
  CONSTANT(M2_PRIMITIVE_AND_FOUR_VELOCITY);
  CONSTANT(M2_VELOCITY1);
  CONSTANT(M2_VELOCITY2);
  CONSTANT(M2_VELOCITY3);
  CONSTANT(M2_MAGNETIC1);
  CONSTANT(M2_MAGNETIC2);
  CONSTANT(M2_MAGNETIC3);
  CONSTANT(M2_VELOCITY_FOUR_VECTOR0);
  CONSTANT(M2_VELOCITY_FOUR_VECTOR1);
  CONSTANT(M2_VELOCITY_FOUR_VECTOR2);
  CONSTANT(M2_VELOCITY_FOUR_VECTOR3);
  CONSTANT(M2_MAGNETIC_FOUR_VECTOR0);
  CONSTANT(M2_MAGNETIC_FOUR_VECTOR1);
  CONSTANT(M2_MAGNETIC_FOUR_VECTOR2);
  CONSTANT(M2_MAGNETIC_FOUR_VECTOR3);
  CONSTANT(M2_COMOVING_MASS_DENSITY);
  CONSTANT(M2_OBSERVER_MASS_DENSITY);
  CONSTANT(M2_TOTAL_ENERGY_DENSITY);
  CONSTANT(M2_KINETIC_ENERGY_DENSITY);
  CONSTANT(M2_INTERNAL_ENERGY_DENSITY);
  CONSTANT(M2_GAS_PRESSURE);
  CONSTANT(M2_MAGNETIC_PRESSURE);
  CONSTANT(M2_PLASMA_BETA);
  CONSTANT(M2_SIGMA);
  CONSTANT(M2_SOUND_SPEED);
  CONSTANT(M2_MACH_NUMBER);
  CONSTANT(M2_ENTROPY);
  CONSTANT(M2_MACH_ALFVEN);
  CONSTANT(M2_MACH_SLOW);
  CONSTANT(M2_MACH_FAST);
  CONSTANT(M2_CT_OUTOFPAGE3);
  CONSTANT(M2_CT_TRANSVERSE1B3);
  CONSTANT(M2_CT_FULL3D);
  CONSTANT(M2_RIEMANN_HLLE);
  CONSTANT(M2_RIEMANN_HLLC);
  CONSTANT(M2_GRADIENT_ESTIMATION_PCM);
  CONSTANT(M2_GRADIENT_ESTIMATION_PLM);
  CONSTANT(M2_QUARTIC_NONE);
  CONSTANT(M2_QUARTIC_FULL_COMPLEX);
  CONSTANT(M2_QUARTIC_ALGORITHMIC);
  CONSTANT(M2_ZONE_TYPE_FULL);
  CONSTANT(M2_ZONE_TYPE_GUARD);
  CONSTANT(M2_ZONE_TYPE_SHELL);
  CONSTANT(M2_STIRRING_NONE);
  CONSTANT(M2_STIRRING_KINETIC);
  CONSTANT(M2_STIRRING_CURRENT);
  CONSTANT(M2_BOUNDARY_PERIODIC);
  CONSTANT(M2_BOUNDARY_OPEN);
  CONSTANT(M2_BOUNDARY_REFLECTING);
  CONSTANT(M2_BOUNDARY_POLAR);
  CONSTANT(M2_BOUNDARY_ORIGIN);
  CONSTANT(M2_BOUNDARY_PROCESSOR);
#undef CONSTANT

  return 1;
}


/* problem-specific Lua modules */
int luaopen_m2jet(lua_State *L);


int main(int argc, char **argv)
{
#if M2_HAVE_MPI
  char *DISABLE_MPI = getenv("M2_DISABLE_MPI");
  if (DISABLE_MPI == NULL) {
    MPI_Init(&argc, &argv);
  }
  else if (strcmp(DISABLE_MPI, "1") != 0) {
    MPI_Init(&argc, &argv);
  }
#endif

  lua_State *L = luaL_newstate();
  char *script;

  if (argc > 1 &&
      strcmp(argv[1] + strlen(argv[1]) - 4, ".lua") == 0) {
    zulu_create_argtable(L, argc-1, argv+1);
    script = argv[1];
  }
  else {
    zulu_create_argtable(L, argc, argv);
    script = M2_INSTALL_PATH"/lib/main.lua";
  }

  luaL_openlibs(L);
  luaL_requiref(L, "buffer", luaopen_buffer, 0); lua_pop(L, 1);
  luaL_requiref(L, "H5", luaopen_H5, 0); lua_pop(L, 1);

  luaL_requiref(L, "30log", luaopen_30log, 0); lua_pop(L, 1);
  luaL_requiref(L, "ansicolors", luaopen_ansicolors, 0); lua_pop(L, 1);
  luaL_requiref(L, "argparse", luaopen_argparse, 0); lua_pop(L, 1);
  luaL_requiref(L, "array", luaopen_array, 0); lua_pop(L, 1);
  luaL_requiref(L, "class", luaopen_class, 0); lua_pop(L, 1);
  luaL_requiref(L, "hdf5", luaopen_hdf5, 0); lua_pop(L, 1);
  luaL_requiref(L, "MPI", luaopen_mpi, 0); lua_pop(L, 1);
  luaL_requiref(L, "json", luaopen_json, 0); lua_pop(L, 1);
  luaL_requiref(L, "serpent", luaopen_serpent, 0); lua_pop(L, 1);
  luaL_requiref(L, "lfs", luaopen_lfs, 0); lua_pop(L, 1);
  luaL_requiref(L, "struct", luaopen_struct, 0); lua_pop(L, 1);

  luaL_requiref(L, "m2lib", luaopen_m2lib, 0); lua_pop(L, 1);
  luaL_requiref(L, "m2jet", luaopen_m2jet, 0); lua_pop(L, 1);

  // Set the Lua path
  // ---------------------------------------------------------------------------
  lua_getglobal(L, "package");
  lua_pushstring(L, M2_INSTALL_PATH"/lib/?.lua;");
  lua_setfield(L, -2, "path");
  lua_pop(L, 1);

  zulu_runscript(L, script);
  lua_close(L);

#if M2_HAVE_MPI
  int mpi_started;
  MPI_Initialized(&mpi_started);
  if (mpi_started) {
    MPI_Finalize();
  }
#endif

  return 0;
}
