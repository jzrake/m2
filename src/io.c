#include "m2.h"
#include "tpl.h"

struct reduced_volume
{
  int global_index[4];
  m2prim prim;
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
  struct reduced_volume v;
  tpl = tpl_map("A(S(i#$(ffffffff)))", &v, 4);

  for (n=0; n<m2->local_grid_size[0]; ++n) {
    for (d=0; d<4; ++d) {
      v.global_index[d] = m2->volumes[n].global_index[d];
    }
    v.prim = m2->volumes[n].prim;
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
  struct reduced_volume v;
  tpl = tpl_map("A(S(i#$(ffffffff)))", &v, 4);

  tpl_load(tpl, TPL_FD, fileno(of));

  for (n=0; n<m2->local_grid_size[0]; ++n) {
    tpl_unpack(tpl, 1);
    for (d=0; d<4; ++d) {
      m2->volumes[n].global_index[d] = v.global_index[d];
    }
    m2->volumes[n].prim = v.prim;
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
