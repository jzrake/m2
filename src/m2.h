#ifndef M2_HEADER
#define M2_HEADER

#include <stdio.h>
#include <stdlib.h>
#include "m2-cfg.h"


/* ------------------------------------------------------------------
 * M2 PUBLIC API
 * --------------------------------------------------------------- */
enum m2flag {
  M2_CARTESIAN,
  M2_CYLINDRICAL,
  M2_SPHERICAL,
  M2_CONSERVED,
  M2_PRIMITIVE,
  M2_AUXILIARY,
  M2_NONRELATIVISTIC=0,
  M2_UNMAGNETIZED=0,
  M2_RELATIVISTIC=(1<<1),
  M2_MAGNETIZED=(1<<2),
} ;

enum { TAU, S11, S22, S33, DDD };

typedef struct m2sim m2sim;
typedef struct m2vol m2vol;
typedef struct m2aux m2aux;
typedef struct m2prim m2prim;
typedef void (*m2vol_operator)(m2vol *vol);

struct m2prim
{
  double v1, v2, v3; /* three-velocity (beta) */
  double p, d;
} ;

struct m2aux
{
  double velocity_four_vector[4];
  double magnetic_four_vector[4];
  double momentum_density[4]; /* T^{0,:} */
  double comoving_mass_density;
  double gas_pressure;
  double magnetic_pressure;
  m2sim *m2;
} ;

struct m2vol
{
  int global_index[4];
  double consA[5];
  double consB[5];
  double flux1[5];
  double flux2[5];
  double flux3[5];
  double area1;
  double area2;
  double area3;
  double volume;
  m2prim prim;
  m2sim *m2;
} ;

struct m2sim
{
  double domain_extent_lower[4];
  double domain_extent_upper[4];
  int domain_resolution[4];
  int local_grid_size[4];
  int number_guard_zones;
  int geometry;
  int physics;
  m2vol *volumes;
} ;


double m2_volume_measure(double x0[4], double x1[4], int geometry);
double m2_area_measure(double x0[4], double x1[4], int geometry, int axis);

double m2vol_minimum_dimension(m2vol *V);

m2sim *m2sim_new();
void m2sim_del(m2sim *m2);
void m2sim_set_resolution(m2sim *m2, int n1, int n2, int n3);
void m2sim_set_guard_zones(m2sim *m2, int ng);
void m2sim_set_extent0(m2sim *m2, double x1, double x2, double x3);
void m2sim_set_extent1(m2sim *m2, double x1, double x2, double x3);
void m2sim_set_geometry(m2sim *m2, int geometry);
void m2sim_set_physics(m2sim *m2, int modes);
void m2sim_initialize(m2sim *m2);
void m2sim_map(m2sim *m2, m2vol_operator f);
void m2sim_index_to_position(m2sim *m2, double index[4], double x[4]);
void m2sim_calculate_conserved(m2sim *m2);
void m2sim_calculate_flux(m2sim *m2);
void m2sim_exchange_flux(m2sim *m2, double dt);
void m2sim_cache_conserved(m2sim *m2);
void m2sim_enforce_boundary_condition(m2sim *m2);


int m2sim_from_conserved_all(m2sim *m2);
int m2sim_from_primitive(m2sim *m2, m2prim *P, double *B, double *X, double dV,
			 double *U, m2aux *aux);
int m2sim_from_conserved(m2sim *m2, double *U, double *B, double *X, double dV,
			 m2aux *aux, m2prim *P);
int m2sim_from_auxiliary(m2sim *m2, m2aux *aux, double *X, double dV,
			 m2prim *P, double *U);
double m2sim_minimum_courant_time(m2sim *m2);


double m2aux_maximum_wavespeed(m2aux *aux);
int m2aux_eigenvalues(m2aux *aux, double n[4], double *evals);
int m2aux_fluxes(m2aux *aux, double n[4], double *F);
int m2aux_add_geometrical_source_terms(m2aux *aux, double x0[4], double x1[4],
				       double *U);







/* ------------------------------------------------------------------
 * DEBUG MACROS AND UTILITY
 * --------------------------------------------------------------- */
#define M2_STRING_LEN 256
#define M2_STRING_SET(d,s) strncpy(d, s, 256); d[255] = '\0';

#define M2_IND(i,j,k) (m2->local_grid_size[2] * m2->local_grid_size[3] * i + \
		       m2->local_grid_size[3] * j +			\
		       k)						\

#define M2_MAX3(a,b,c) ((a)>(b))?(((a)>(c))?(a):(c)):(((b)>(c))?(b):(c))
#define M2_MIN3(a,b,c) ((a)<(b))?(((a)<(c))?(a):(c)):(((b)<(c))?(b):(c))

#define DEBUG 'D'
#define INFO 'I'
#define WARNING 'W'
#define ERROR 'E'
#define FATAL 'F'

#ifndef M2_DEBUG_FILE
#define DEBUG_FILE NULL
#else
#define DEBUG_FILE M2_DEBUG_FILE
#endif // M2_DEBUG_FILE
#define INFO_FILE stdout
#define WARNING_FILE stderr
#define ERROR_FILE stderr
#define FATAL_FILE stderr

#define MSGF(level, format, ...) do {				\
    FILE *out = NULL;						\
    char *pre = "";						\
    switch (level) {						\
    case DEBUG: out = DEBUG_FILE; pre = "DEBUG"; break;		\
    case INFO: out = INFO_FILE; pre = "INFO"; break;		\
    case WARNING: out = WARNING_FILE; pre = "WARNING"; break;	\
    case ERROR: out = ERROR_FILE; pre = "ERROR"; break; 	\
    case FATAL: out = FATAL_FILE; pre = "FATAL"; break; 	\
    }								\
    if (out) {							\
      fprintf(out, "[%s:%s] ", pre, __FUNCTION__);		\
      fprintf(out, format, __VA_ARGS__);			\
      fprintf(out, "\n");					\
    }								\
    if (level == FATAL) {					\
      exit(1);							\
    }								\
  } while (0)							\

#define MSG(level, message) MSGF(level, "%s", message)

#include <time.h>
#define TIME(cmd) do {                                  \
    clock_t start = clock();                            \
    cmd;                                                \
    printf("[%s] %s took %5.4f ms\n", __FUNCTION__,	\
           #cmd, 1e3*(clock() - start)/CLOCKS_PER_SEC); \
  } while (0)

/* assertion macro */
#include <assert.h>
#define NOHUP 0 // continue even if an assertion fails
#define ASSERT_L "[assertion:%s]$ %s == %s : %"PRIu64"\n"
#define ASSERT_I "[assertion:%s]$ %s == %s : %d\n"
#define ASSERT_F "[assertion:%s]$ %s == %s : %f\n"
#define ASSERTEQL(E,v)printf(ASSERT_L,__FUNCTION__,#E,#v,E);assert(E==v||NOHUP);
#define ASSERTEQI(E,v)printf(ASSERT_I,__FUNCTION__,#E,#v,E);assert(E==v||NOHUP);
#define ASSERTEQF(E,v)printf(ASSERT_F,__FUNCTION__,#E,#v,E);	\
  assert(fabs(E-v) < 1e-14||NOHUP);

#endif // M2_HEADER
