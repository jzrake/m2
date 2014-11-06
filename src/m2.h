#ifndef M2_HEADER
#define M2_HEADER

#include <stdio.h>
#include <stdlib.h>
#include "m2-cfg.h"
#include "jsw_rand.h"


/* ------------------------------------------------------------------
 * M2 PUBLIC API
 * --------------------------------------------------------------- */
enum m2flag {
  /* geometry */
  M2_CARTESIAN,
  M2_CYLINDRICAL,
  M2_SPHERICAL,
  M2_PARABOLIC,
  M2_LINEAR,
  M2_LOGARITHMIC,

  /* may have many uses including interpolation variables */
  M2_CONSERVED,
  M2_AUXILIARY,
  M2_PRIMITIVE,
  M2_PRIMITIVE_AND_FOUR_VELOCITY,

  /* m2aux data members and other measurements */
  M2_VELOCITY1,
  M2_VELOCITY2,
  M2_VELOCITY3,
  M2_MAGNETIC1,
  M2_MAGNETIC2,
  M2_MAGNETIC3,
  M2_VELOCITY_FOUR_VECTOR0,
  M2_VELOCITY_FOUR_VECTOR1,
  M2_VELOCITY_FOUR_VECTOR2,
  M2_VELOCITY_FOUR_VECTOR3,
  M2_MAGNETIC_FOUR_VECTOR0,
  M2_MAGNETIC_FOUR_VECTOR1,
  M2_MAGNETIC_FOUR_VECTOR2,
  M2_MAGNETIC_FOUR_VECTOR3,
  M2_COMOVING_MASS_DENSITY,
  M2_OBSERVER_MASS_DENSITY,
  M2_TOTAL_ENERGY_DENSITY,
  M2_KINETIC_ENERGY_DENSITY,
  M2_INTERNAL_ENERGY_DENSITY,
  M2_MAGNETIC_ENERGY_DENSITY,
  M2_GAS_PRESSURE,
  M2_MAGNETIC_PRESSURE,
  M2_PLASMA_BETA,
  M2_SIGMA,
  M2_SOUND_SPEED,
  M2_MACH_NUMBER,
  M2_ENTROPY,

  /* Lorentz-scalar directional Mach number measurements used with m2aux_mach */
  M2_MACH_ALFVEN,
  M2_MACH_SLOW,
  M2_MACH_FAST,

  /* CT schemes */
  M2_CT_OUTOFPAGE3,
  M2_CT_TRANSVERSE1B3,
  M2_CT_FULL3D,

  /* Riemann solvers */
  M2_RIEMANN_HLLE,
  M2_RIEMANN_HLLC,

  /* Quartic polynomial solvers */
  M2_QUARTIC_NONE,
  M2_QUARTIC_FULL_COMPLEX,
  M2_QUARTIC_ALGORITHMIC,

  /* Zone types */
  M2_ZONE_TYPE_FULL,
  M2_ZONE_TYPE_GUARD,
  M2_ZONE_TYPE_SHELL, /* indicates that only faces and edges are used */

  /* Zone status */
  M2_ZONE_HEALTH_GOOD,
  M2_ZONE_HEALTH_BAD,

  /* Stirring types */
  M2_STIRRING_NONE,
  M2_STIRRING_KINETIC,
  M2_STIRRING_CURRENT,
} ;

enum { M2_BIT_BAD_EIGENVALUES,
       M2_BIT_FAILED_CONSERVED_INVERSION,
       M2_BIT_NEARLY_FAILED_CONSERVED_INVERSION,
       M2_BIT_NEGATIVE_PRESSURE,
       M2_BIT_NEGATIVE_ENERGY,
       M2_BIT_NEGATIVE_DENSITY,
} ;

/* indices into cons[AB] and flux[123] arrays */
enum { TAU, S11, S22, S33, DDD, B11, B22, B33 };
/* indices into prim array if not m2prim struct */
enum { RHO, PRE, V11, V22, V33 };

typedef struct m2sim m2sim;
typedef struct m2vol m2vol;
typedef struct m2aux m2aux;
typedef struct m2prim m2prim;
typedef struct m2sim_status m2sim_status;
typedef struct m2reductions m2reductions;
typedef struct m2stirring m2stirring;

typedef void (*m2vol_operator)(m2vol *vol);
typedef void (*m2sim_operator)(m2sim *m2);
typedef void (*m2crd_function)(m2sim *m2, double X[4], double N[4], double *Y);

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
  m2prim prim;
  double Bflux1A; /* magnetic flux through (+1/2) faces */
  double Bflux2A;
  double Bflux3A;
  /* end serialized data members */
  double x0[4];
  double x1[4];
  double volume;
  double consA[5];
  double consB[5]; /* RK cache for cons */
  double Bflux1B; /* RK cache for Bflux */
  double Bflux2B;
  double Bflux3B;
  double flux1[8]; /* Godunov fluxes through (+1/2) faces */
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
  double emf1; /* line integral of electric field along (+1/2, +1/2) edges */
  double emf2;
  double emf3;
  char zone_type;
  char zone_health;
  m2sim *m2;
  m2aux aux;
} ;
#define M2_VOL_SERIALIZE(b) ("A(S(i#$(ffffffff)fff))",b,4)

struct m2stirring
{
  double *wavenumbers; /* num_waves x 4 */
  double *fourrieramp; /* num_waves x 4 */
  double correlation_time;
  double amplitude;
  double wavelength;
  int num_waves;
  int stirring_type;
  jsw_rand_t random;
} ;

struct m2sim_status
{
  double time_simulation;
  double time_step;
  double time_last_checkpoint_hdf5;
  double time_last_checkpoint_tpl;
  double time_last_analysis;
  int iteration_number;
  int checkpoint_number_tpl;
  int checkpoint_number_hdf5;
  int error_code;
  int drive_attempt;
  int num_failures[8];
  char message[1024];
} ;

struct m2reductions
{
  double total_volume;
  double total_energy;
  double total_mass;
  double total_kinetic_energy;
  double total_internal_energy;
  double total_magnetic_energy;
} ;

struct m2sim
{
  double domain_extent_lower[4];
  double domain_extent_upper[4];
  int domain_resolution[4];
  int local_grid_size[4];
  int local_grid_start[4];
  int number_guard_zones0[4];
  int number_guard_zones1[4];
  int periodic_dimension[4];
  int coordinate_scaling1; /* logarithmic, linear, etc */
  int coordinate_scaling2;
  int coordinate_scaling3;
  int geometry;
  int relativistic;
  int magnetized;
  int rk_order;
  int interpolation_fields;
  int riemann_solver;
  int quartic_solver;
  int suppress_extrapolation_at_unhealthy_zones;
  double pressure_floor;
  double plm_parameter;
  double cfl_parameter;
  double gamma_law_index;
  double cadence_checkpoint_hdf5;
  double cadence_checkpoint_tpl;
  double cadence_analysis;
  m2stirring stirring;
  m2sim_status status;
  m2sim_operator analysis;
  m2sim_operator boundary_conditions_gradient;
  m2sim_operator boundary_conditions_cell;
  m2vol_operator boundary_conditions_flux1;
  m2vol_operator boundary_conditions_flux2;
  m2vol_operator boundary_conditions_flux3;
  m2vol_operator boundary_conditions_emf1;
  m2vol_operator boundary_conditions_emf2;
  m2vol_operator boundary_conditions_emf3;
  m2vol_operator initial_data;
  m2crd_function initial_data_cell;
  m2crd_function initial_data_face;
  m2crd_function initial_data_edge;
  m2vol_operator add_physical_source_terms;
  m2vol *volumes;
  void *user_struct;
  void *cart_comm;
  void *mpi_types;
} ;
#define M2_SIM_SERIALIZE(m2) ("S(f#f#i#i#i#i#i#iiiiiiiiiiffffff$(fffffiii))", \
			      m2,4,4,4,4,4,4,4)


void m2_self_test();
void m2_print_state(m2prim *P, m2aux *aux, double *U);
double m2_volume_measure(double x0[4], double x1[4], int geometry);
double m2_area_measure(double x0[4], double x1[4], int geometry, int axis);
double m2_line_measure(double x0[4], double x1[4], int geometry, int axis);
void m2_to_cartesian(double x[4], double xcart[4], int geometry);
int m2_solve_quartic_equation1(double d4, double d3,
			       double d2, double d1, double d0,
			       double roots[4]);
int m2_solve_quartic_equation2(double d4, double d3,
			       double d2, double d1, double d0,
			       double roots[4]);
double m2_force_free_vector_potential(double x[4], double n[4], int model);
double m2_force_free_magnetic_field(double x[4], double n[4], int model);
double m2_magnetic_rope_vector_potential(double x[4], double n[4]);

double m2_stochastic_vector_potential(double x[4], double n[4]);
double m2_magnetic_rope_vector_potential(double x[4], double n[4]);


/* vol */
double m2vol_minimum_dimension(m2vol *V);
double m2vol_coordinate_centroid(m2vol *V, int axis);
void m2vol_coordinate_centroid_3d(m2vol *V, double x[4]);
void m2vol_to_interpolated(m2vol *V, double *y, int stride);
m2vol *m2vol_neighbor(m2vol *V, int axis, int dist);
int m2vol_from_primitive(m2vol *V);


/* sim */
void m2sim_new(m2sim *m2);
void m2sim_del(m2sim *m2);
void m2sim_initialize(m2sim *m2);
void m2sim_create_mpi_types(m2sim *m2);
void m2sim_delete_mpi_types(m2sim *m2);
void m2sim_index_to_position(m2sim *m2, double index[4], double x[4]);
void m2sim_volume_at_position(m2sim *m2, double x[4], m2vol **V, double dx[4]);
void m2sim_primitive_at_position(m2sim *m2, double x[4], m2prim *P);
void m2sim_index_global_to_local(m2sim *m2, int global_index[4], int I[4]);
int m2sim_calculate_gradient(m2sim *m2);
int m2sim_calculate_flux(m2sim *m2);
int m2sim_calculate_emf(m2sim *m2);
int m2sim_exchange_flux(m2sim *m2, double dt);
int m2sim_add_source_terms(m2sim *m2, double dt);
int m2sim_average_runge_kutta(m2sim *m2, double b);
int m2sim_magnetic_flux_to_cell_center(m2sim *m2);
int m2sim_cache_conserved(m2sim *m2);
int m2sim_enforce_boundary_condition(m2sim *m2);
int m2sim_synchronize_guard(m2sim *m2);
int m2sim_save_checkpoint_tpl(m2sim *m2, const char *fname);
int m2sim_load_checkpoint_tpl(m2sim *m2, const char *fname);
void m2sim_print(m2sim *m2);
void m2sim_write_ascii_1d(m2sim *m2, const char *fname);
void m2sim_write_ascii_2d(m2sim *m2, const char *fname);
void m2sim_run_analysis(m2sim *m2);
void m2sim_run_initial_data(m2sim *m2);
void m2sim_from_interpolated(m2sim *m2, double *y, m2prim *P);
int m2sim_memory_usage(m2sim *m2);
int m2sim_from_conserved_all(m2sim *m2);
int m2sim_from_primitive_all(m2sim *m2);
int m2sim_from_conserved(m2sim *m2, double *U, double *B, double *X, double dV,
			 m2aux *aux, m2prim *P);
int m2sim_from_primitive(m2sim *m2, m2prim *P, double *B, double *X, double dV,
			 double *U, m2aux *aux);
int m2sim_from_auxiliary(m2sim *m2, m2aux *aux, double *X, double dV,
			 m2prim *P, double *U);
double m2sim_minimum_courant_time(m2sim *m2);
double m2sim_volume_integral(m2sim *m2, int member, int (*cut_cb)(m2vol *V));
void m2sim_drive(m2sim *m2);
void m2sim_visualize(m2sim *m2, int argc, char **argv);


/* aux */
double m2aux_maximum_wavespeed(m2aux *aux, int *err);
int m2aux_eigenvalues(m2aux *aux, double n[4], double *evals);
int m2aux_fluxes(m2aux *aux, double n[4], double *F);
int m2aux_add_geometrical_source_terms(m2aux *aux, double x0[4], double x1[4],
				       double *U);
double m2aux_get(m2aux *aux, int member);
double m2aux_measure(m2aux *aux, int member);
double m2aux_mach(m2aux *aux, double n[4], int mode);


/* stirring */
void m2stirring_init(m2stirring *S);
void m2stirring_del(m2stirring *S);
void m2stirring_next_realization(m2stirring *S);
double m2stirring_get_field(m2stirring *S, double x[4], double n[4]);



/* ------------------------------------------------------------------
 * DEBUG MACROS AND UTILITY
 * --------------------------------------------------------------- */
#define M2_PI (4*atan(1))
#define M2_STRING_LEN 256
#define M2_STRING_SET(d,s) strncpy(d, s, 256); d[255] = '\0';

#define M2_IND(i,j,k)((i)*L[2]*L[3]+(j)*L[3]+(k))
#define M2_VOL(i,j,k)(((i)>=0)&&(i)<L[1]&&				\
		      ((j)>=0)&&(j)<L[2]&&				\
		      ((k)>=0)&&(k)<L[3]?m2->volumes+M2_IND(i,j,k):NULL)
#define M2_MAX3(a,b,c)(((a)>(b))?(((a)>(c))?(a):(c)):(((b)>(c))?(b):(c)))
#define M2_MIN3(a,b,c)(((a)<(b))?(((a)<(c))?(a):(c)):(((b)<(c))?(b):(c)))

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

#endif /* M2_HEADER */
