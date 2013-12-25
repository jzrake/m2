#include "m2.h"
#include "tpl.h"

struct reduced_volume
{
  int global_index[4];
  m2prim prim;
  double Bflux1;
  double Bflux2;
  double Bflux3;
} ;

void m2sim_print(m2sim *m2)
{
  printf("domain_extent_lower: [%f %f %f]\n",
	 m2->domain_extent_lower[1],
	 m2->domain_extent_lower[2],
	 m2->domain_extent_lower[3]);
  printf("domain_extent_upper: [%f %f %f]\n",
	 m2->domain_extent_upper[1],
	 m2->domain_extent_upper[2],
	 m2->domain_extent_upper[3]);
  printf("domain_resolution: [%d %d %d]\n",
	 m2->domain_resolution[1],
	 m2->domain_resolution[2],
	 m2->domain_resolution[3]);
  printf("magnetized: %s\n", m2->physics & M2_MAGNETIZED ? "yes" : "no");
  printf("relativistic: %s\n", m2->physics & M2_RELATIVISTIC ? "yes" : "no");
}

void m2sim_save_header(m2sim *m2, FILE *of)
{
  printf("\t[save header]\n");
  tpl_node *tpl;
  tpl = tpl_map("S(f#f#i#i#iii)", m2, 4, 4, 4, 4);
  tpl_pack(tpl, 0);
  tpl_dump(tpl, TPL_FD, fileno(of));
  tpl_free(tpl);
}

void m2sim_load_header(m2sim *m2, FILE *of)
{
  printf("\t[load header]\n");
  tpl_node *tpl;
  tpl = tpl_map("S(f#f#i#i#iii)", m2, 4, 4, 4, 4);
  tpl_load(tpl, TPL_FD, fileno(of));
  tpl_unpack(tpl, 0);
  tpl_free(tpl);
}

void m2sim_save_volumes(m2sim *m2, FILE *of)
{
  printf("\t[save volumes]\n");
  int n, d;
  tpl_node *tpl;
  m2vol *V;
  struct reduced_volume v;
  tpl = tpl_map("A(S(i#$(ffffffff)fff))", &v, 4);
  for (n=0; n<m2->local_grid_size[0]; ++n) {
    V = m2->volumes + n;
    for (d=0; d<4; ++d) {
      v.global_index[d] = V->global_index[d];
    }
    v.prim = V->prim;
    v.Bflux1 = V->Bflux1A;
    v.Bflux2 = V->Bflux2A;
    v.Bflux3 = V->Bflux3A;
    tpl_pack(tpl, 1);
  }
  tpl_dump(tpl, TPL_FD, fileno(of));
  tpl_free(tpl);
}

void m2sim_load_volumes(m2sim *m2, FILE *of)
{
  printf("\t[load volumes]\n");
  int n, d;
  tpl_node *tpl;
  m2vol *V;
  struct reduced_volume v;
  tpl = tpl_map("A(S(i#$(ffffffff)fff))", &v, 4);
  tpl_load(tpl, TPL_FD, fileno(of));
  for (n=0; n<m2->local_grid_size[0]; ++n) {
    V = m2->volumes + n;
    tpl_unpack(tpl, 1);
    for (d=0; d<4; ++d) {
      if (V->global_index[d] != v.global_index[d]) {
	MSG(FATAL, "data format mismatch");
      }
    }
    V->prim = v.prim;
    V->Bflux1A = v.Bflux1;
    V->Bflux2A = v.Bflux2;
    V->Bflux3A = v.Bflux3;
  }
  tpl_free(tpl);
}

void m2sim_save_checkpoint(m2sim *m2, char *fname)
{
  printf("[save checkpoint: %s]\n", fname);
  FILE *of = fopen(fname, "w");
  m2sim_save_header(m2, of);
  m2sim_save_volumes(m2, of);
  fclose(of);
}

void m2sim_load_checkpoint(m2sim *m2, char *fname)
{
  printf("[load checkpoint: %s]\n", fname);
  FILE *of = fopen(fname, "r");
  m2sim_load_header(m2, of);
  m2sim_initialize(m2);
  m2sim_load_volumes(m2, of);
  fclose(of);
}

void m2sim_write_ascii_1d(m2sim *m2, char *fname)
{
  printf("[write ascii 1d: %s]\n", fname);
  FILE *outfile = fopen(fname, "w");
  int *L = m2->local_grid_size;
  int i;
  for (i=0; i<m2->local_grid_size[1]; ++i) {
    fprintf(outfile, "%8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e\n",
	    m2vol_coordinate_centroid(&m2->volumes[M2_IND(i,0,0)], 1),
	    m2->volumes[M2_IND(i,0,0)].prim.d,
	    m2->volumes[M2_IND(i,0,0)].prim.p,
	    m2->volumes[M2_IND(i,0,0)].prim.v1,
	    m2->volumes[M2_IND(i,0,0)].prim.v2,
	    m2->volumes[M2_IND(i,0,0)].prim.v3,
	    m2->volumes[M2_IND(i,0,0)].prim.B1,
	    m2->volumes[M2_IND(i,0,0)].prim.B2,
	    m2->volumes[M2_IND(i,0,0)].prim.B3);
  }
  fclose(outfile);
}

void m2sim_write_ascii_2d(m2sim *m2, char *fname)
{
  printf("[write ascii 2d: %s]\n", fname);
  FILE *outfile = fopen(fname, "w");
  int *L = m2->local_grid_size;
  int i, j;
  for (i=0; i<m2->local_grid_size[1]; ++i) {
    for (j=0; j<m2->local_grid_size[2]; ++j) {
      fprintf(outfile, "%8.6e %8.6e %8.6e\n",
	      m2vol_coordinate_centroid(&m2->volumes[M2_IND(i,j,0)], 1),
	      m2vol_coordinate_centroid(&m2->volumes[M2_IND(i,j,0)], 2),
	      m2->volumes[M2_IND(i,j,0)].prim.d);
    }
  }
  fclose(outfile);
}
