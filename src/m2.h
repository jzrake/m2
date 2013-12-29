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

  /* integer values for getting m2aux data members */
  M2_VELOCITY_FOUR_VECTOR0,
  M2_VELOCITY_FOUR_VECTOR1,
  M2_VELOCITY_FOUR_VECTOR2,
  M2_VELOCITY_FOUR_VECTOR3,
  M2_MAGNETIC_FOUR_VECTOR0,
  M2_MAGNETIC_FOUR_VECTOR1,
  M2_MAGNETIC_FOUR_VECTOR2,
  M2_MAGNETIC_FOUR_VECTOR3,
  M2_COMOVING_MASS_DENSITY,
  M2_GAS_PRESSURE,
  M2_MAGNETIC_PRESSURE,

  M2_NONRELATIVISTIC=0,
  M2_UNMAGNETIZED=0,
  M2_RELATIVISTIC=(1<<1),
  M2_MAGNETIZED=(1<<2),
} ;

/* indices into cons[AB] and flux[123] arrays */
enum { TAU, S11, S22, S33, DDD, B11, B22, B33 };

typedef struct m2sim m2sim;
typedef struct m2vol m2vol;
typedef struct m2aux m2aux;
typedef struct m2prim m2prim;
typedef struct m2sim_status m2sim_status;
typedef void (*m2vol_operator)(m2vol *vol);

struct m2prim
{
  double v1, v2, v3; /* three-velocity (beta) */
  double B1, B2, B3; /* cell-centered magnetic field */
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
  double x0[4];
  double x1[4];
  double consA[5];
  double consB[5];
  double flux1[8];
  double flux2[8];
  double flux3[8];
  double grad1[8];
  double grad2[8];
  double grad3[8];
  double area1;
  double area2;
  double area3;
  double line1;
  double line2;
  double line3;
  double volume;
  m2prim prim;
  double Bflux1A, Bflux1B; /* magnetic flux through (+1/2) faces */
  double Bflux2A, Bflux2B;
  double Bflux3A, Bflux3B;
  double emf1; /* line integral of electric field along (+1/2, +1/2) edges */
  double emf2;
  double emf3;
  m2aux aux;
  m2sim *m2;
} ;

struct m2sim_status
{
  double time_simulation;
  int iteration_number;
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
  m2sim_status status;
  m2vol *volumes;
} ;


double m2_volume_measure(double x0[4], double x1[4], int geometry);
double m2_area_measure(double x0[4], double x1[4], int geometry, int axis);
double m2_line_measure(double x0[4], double x1[4], int geometry, int axis);
void m2_to_cartesian(double x[4], double xcart[4], int geometry);

double m2vol_minimum_dimension(m2vol *V);
double m2vol_coordinate_centroid(m2vol *V, int axis);
void m2vol_coordinate_centroid_3d(m2vol *V, double x[4]);

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
void m2sim_calculate_emf(m2sim *m2);
void m2sim_exchange_flux(m2sim *m2, double dt);
void m2sim_add_source_terms(m2sim *m2, double dt);
void m2sim_cache_conserved(m2sim *m2);
void m2sim_average_runge_kutta(m2sim *m2, double b);
void m2sim_magnetic_flux_to_cell_center(m2sim *m2);
void m2sim_enforce_boundary_condition(m2sim *m2);
void m2sim_print(m2sim *m2);
void m2sim_save_checkpoint(m2sim *m2, char *fname);
void m2sim_load_checkpoint(m2sim *m2, char *fname);
void m2sim_write_ascii_1d(m2sim *m2, char *fname);
void m2sim_write_ascii_2d(m2sim *m2, char *fname);
int m2sim_memory_usage(m2sim *m2);

int m2sim_from_conserved_all(m2sim *m2);
int m2sim_from_primitive(m2sim *m2, m2prim *P, double *B, double *X, double dV,
			 double *U, m2aux *aux);
int m2sim_from_conserved(m2sim *m2, double *U, double *B, double *X, double dV,
			 m2aux *aux, m2prim *P);
int m2sim_from_auxiliary(m2sim *m2, m2aux *aux, double *X, double dV,
			 m2prim *P, double *U);
double m2sim_minimum_courant_time(m2sim *m2);

void m2sim_drive(m2sim *m2);
void m2sim_visualize(m2sim *m2, int argc, char **argv);

double m2aux_maximum_wavespeed(m2aux *aux);
int m2aux_eigenvalues(m2aux *aux, double n[4], double *evals);
int m2aux_fluxes(m2aux *aux, double n[4], double *F);
int m2aux_add_geometrical_source_terms(m2aux *aux, double x0[4], double x1[4],
				       double *U);
int m2aux_get(m2aux *aux, int member, double *value);






/* ------------------------------------------------------------------
 * DEBUG MACROS AND UTILITY
 * --------------------------------------------------------------- */
#define M2_PI 3.14159265358979
#define M2_STRING_LEN 256
#define M2_STRING_SET(d,s) strncpy(d, s, 256); d[255] = '\0';

#define M2_IND(i,j,k)(i)*L[2]*L[3]+(j)*L[3]+(k)
#define M2_VOL(i,j,k)(i)<L[1]&&(j)<L[2]&&(k)<L[3]?m2->volumes+M2_IND(i,j,k):NULL
#define M2_MAX3(a,b,c)((a)>(b))?(((a)>(c))?(a):(c)):(((b)>(c))?(b):(c))
#define M2_MIN3(a,b,c)((a)<(b))?(((a)<(c))?(a):(c)):(((b)<(c))?(b):(c))

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
