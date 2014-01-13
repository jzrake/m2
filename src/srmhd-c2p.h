#ifndef SRMHD_C2P_HEADER
#define SRMHD_C2P_HEADER

#ifdef __cplusplus
extern "C" {
#endif

  enum {
    SRMHD_C2P_SUCCESS,
    SRMHD_C2P_CONS_CONTAINS_NAN,
    SRMHD_C2P_CONS_NEGATIVE_DENSITY,
    SRMHD_C2P_CONS_NEGATIVE_ENERGY,
    SRMHD_C2P_PRIM_CONTAINS_NAN,
    SRMHD_C2P_PRIM_NEGATIVE_PRESSURE,
    SRMHD_C2P_PRIM_NEGATIVE_RESTMASS,
    SRMHD_C2P_PRIM_SUPERLUMINAL,
    SRMHD_C2P_MAXITER
  } ;

  void srmhd_c2p_set_gamma(double adiabatic_gamma);
  void srmhd_c2p_set_pressure_floor(double pf);
  void srmhd_c2p_new_state(const double *U);
  void srmhd_c2p_estimate_from_cons();
  void srmhd_c2p_set_starting_prim(const double *P);
  void srmhd_c2p_get_starting_prim(double *P);
  int srmhd_c2p_reconstruct_prim(double Z, double W, double *P);
  int srmhd_c2p_solve_anton2dzw(double *P);
  int srmhd_c2p_solve_noble1dw(double *P);
  int srmhd_c2p_get_iterations();
  int srmhd_c2p_check_cons(const double *U);
  int srmhd_c2p_check_prim(const double *P);
  char *srmhd_c2p_get_error(int error);

  void srmhd_c2p_eos_set_eos(const void *eos_);
  void srmhd_c2p_eos_new_state(const double *U);
  void srmhd_c2p_eos_estimate_from_cons();
  void srmhd_c2p_eos_set_starting_prim(const double *P);
  int srmhd_c2p_eos_solve_noble2dzt(double *P);
  int srmhd_c2p_eos_solve_duffell3d(double *P);
  int srmhd_c2p_eos_get_iterations();

#ifdef __cplusplus
}
#endif

#endif /* SRMHD_C2P_HEADER */
