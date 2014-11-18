#include <stdlib.h>
#include "m2.h"



static inline void linear_to_multi(int n, int dims[4], int index[4])
{
  int S[4] = {0, dims[3] * dims[2], dims[3], 1};
  index[1] = n / S[1]; n -= index[1] * S[1];
  index[2] = n / S[2]; n -= index[2] * S[2];
  index[3] = n / S[3]; n -= index[3] * S[3]; 
}



static inline int multi_to_linear(int dims[4], int index[4])
{
  int S[4] = {0, dims[3] * dims[2], dims[3], 1};
  int n = S[1]*index[1] + S[2]*index[2] + S[3]*index[3];
  if (index[1] < 0 || dims[1] <= index[1] ||
      index[2] < 0 || dims[2] <= index[2] ||
      index[3] < 0 || dims[3] <= index[3]) {
    return -1;
  }
  else {
    return n;
  }
}



void mesh_linear_to_multi(int n, int dims[4], int index[4])
{
  linear_to_multi(n, dims, index);
}



int mesh_multi_to_linear(int dims[4], int index[4])
{
  return multi_to_linear(dims, index);
}



void mesh_new(struct mesh *M)
{
  int d, e;
  M->cells = NULL;
  for (d=0; d<=3; ++d) {
    M->cells_shape[d] = 0;
    for (e=0; e<=3; ++e) {
      M->edges_shape[d][e] = 0;
      M->faces_shape[d][e] = 0;
    }
    M->edges[d] = NULL;
    M->faces[d] = NULL;
  }
}



void mesh_allocate(struct mesh *M)
{
#define E M->edges_shape
#define F M->faces_shape
#define C M->cells_shape
  int d, n;
  for (d=1; d<=3; ++d) {
    F[d][1] = C[1] + (d == 1);
    F[d][2] = C[2] + (d == 2);
    F[d][3] = C[3] + (d == 3);
    E[d][1] = C[1] + (d != 1);
    E[d][2] = C[2] + (d != 2);
    E[d][3] = C[3] + (d != 3);
  }
  for (d=1; d<=3; ++d) {
    E[d][0] = E[d][1] * E[d][2] * E[d][3];
    F[d][0] = F[d][1] * F[d][2] * F[d][3];
  }
  C[0] = C[1] * C[2] * C[3];
  for (d=1; d<=3; ++d) {
    M->edges[d] = (struct mesh_edge*) malloc(E[d][0]*sizeof(struct mesh_edge));
    M->faces[d] = (struct mesh_face*) malloc(F[d][0]*sizeof(struct mesh_face));
  }
  M->cells = (struct mesh_cell*) malloc(C[0]*sizeof(struct mesh_cell));
  for (d=1; d<=3; ++d) {
    for (n=0; n<E[d][0]; ++n) {
      M->edges[d][n].id = n;
      M->edges[d][n].axis = d;
    }
    for (n=0; n<F[d][0]; ++n) {
      M->faces[d][n].id = n;
      M->faces[d][n].axis = d;
    }
  }
  for (n=0; n<C[0]; ++n) {
    M->cells[n].id = n;
  }
  M->cells_stride[1] = C[3] * C[2];
  M->cells_stride[2] = C[3];
  M->cells_stride[3] = 1;
  for (d=1; d<=3; ++d) {
    M->edges_stride[d][1] = E[d][3] * E[d][2];
    M->edges_stride[d][2] = E[d][3];
    M->edges_stride[d][3] = 1;
    M->faces_stride[d][1] = F[d][3] * F[d][2];
    M->faces_stride[d][2] = F[d][3];
    M->faces_stride[d][3] = 1;
  }
#undef E
#undef F
#undef C
}



void mesh_deallocate(struct mesh *M)
{
  int d;
  for (d=1; d<=3; ++d) {
    free(M->faces[d]);
    free(M->edges[d]);
  }
  free(M->cells);
}



int mesh_cell_faces1(struct mesh *M, int cell_id, int face_ids[6])
{
#define C M->cells_shape
#define F M->faces_shape
  int I[4], n, d;
  linear_to_multi(cell_id, C, I);
  for (d=1; d<=3; ++d) {
    n = 2*(d-1);
    face_ids[n+0] = multi_to_linear(F[d], I); I[d] += 1;
    face_ids[n+1] = multi_to_linear(F[d], I); I[d] -= 1;
  }
  return ((face_ids[0] != -1) +
	  (face_ids[1] != -1) +
	  (face_ids[2] != -1) +
	  (face_ids[3] != -1) +
	  (face_ids[4] != -1) +
	  (face_ids[5] != -1));
#undef C
#undef F
}



int mesh_face_edges1(struct mesh *M, int axis, int face_id, int edge_ids[4])
{
#define E M->edges_shape
#define F M->faces_shape
  int I[4];
  int a2 = (axis + 1 - 1) % 3 + 1;
  int a3 = (axis + 2 - 1) % 3 + 1;
  linear_to_multi(face_id, F[axis], I);
  edge_ids[0] = multi_to_linear(E[a2], I); I[a2] += 1; /* TODO: test (probably correct) */
  edge_ids[1] = multi_to_linear(E[a2], I); I[a2] -= 1;
  edge_ids[2] = multi_to_linear(E[a3], I); I[a3] += 1;
  edge_ids[3] = multi_to_linear(E[a3], I); I[a3] -= 1;
  return ((edge_ids[0] != -1) +
	  (edge_ids[1] != -1) +
	  (edge_ids[2] != -1) +
	  (edge_ids[3] != -1));
#undef E
#undef F
}



int mesh_face_cells1(struct mesh *M, int axis, int face_id, int cell_ids[2])
{
#define C M->cells_shape
#define F M->faces_shape
  int I[4];
  linear_to_multi(face_id, F[axis], I);
  I[axis] -= 1; cell_ids[0] = multi_to_linear(C, I);
  I[axis] += 1; cell_ids[1] = multi_to_linear(C, I);
  return ((cell_ids[0] != -1) +
	  (cell_ids[1] != -1));
#undef C
#undef F
}



int mesh_edge_faces1(struct mesh *M, int axis, int edge_id, int face_ids[4])
{
#define E M->edges_shape
#define F M->faces_shape
  int I[4];
  int a2 = (axis + 1 - 1) % 3 + 1;
  int a3 = (axis + 2 - 1) % 3 + 1;
  linear_to_multi(edge_id, E[axis], I);
  I[a3] -= 1; face_ids[0] = multi_to_linear(F[a2], I);
  I[a3] += 1; face_ids[1] = multi_to_linear(F[a2], I);
  I[a2] -= 1; face_ids[2] = multi_to_linear(F[a3], I);
  I[a2] += 1; face_ids[3] = multi_to_linear(F[a3], I);
  return ((face_ids[0] != -1) +
	  (face_ids[1] != -1) +
	  (face_ids[2] != -1) +
	  (face_ids[3] != -1));
#undef E
#undef F  
}



int mesh_edge_cells1(struct mesh *M, int axis, int edge_id, int cell_ids[4])
{
#define E M->edges_shape
#define F M->faces_shape
#define C M->cells_shape
  int I[4];
  int a2 = (axis + 1 - 1) % 3 + 1;
  int a3 = (axis + 2 - 1) % 3 + 1;
  linear_to_multi(edge_id, E[axis], I);
  cell_ids[0] = multi_to_linear(C, I); I[a2] += 1; /* TODO: test */
  cell_ids[1] = multi_to_linear(C, I); I[a2] -= 1;
  cell_ids[2] = multi_to_linear(C, I); I[a3] += 1;
  cell_ids[3] = multi_to_linear(C, I); I[a3] -= 1;
  return ((cell_ids[0] != -1) +
	  (cell_ids[1] != -1) +
	  (cell_ids[2] != -1) +
	  (cell_ids[3] != -1));
#undef E
#undef F
#undef C
}



int mesh_cell_faces(struct mesh *M, struct mesh_cell *cell, struct mesh_face *faces[6])
{
  int face_ids[6];
  int n, d;
  int num = mesh_cell_faces1(M, cell->id, face_ids);
  for (d=1; d<=3; ++d) {
    n = 2*(d-1);
    faces[n+0] = face_ids[n+0] == -1 ? NULL : M->faces[d] + face_ids[n+0];
    faces[n+1] = face_ids[n+1] == -1 ? NULL : M->faces[d] + face_ids[n+1];
  }
  return num;
}



int mesh_face_edges(struct mesh *M, struct mesh_face *face, struct mesh_edge *edges[4])
{
  int edge_ids[4];
  int a2 = (face->axis + 1 - 1) % 3 + 1;
  int a3 = (face->axis + 2 - 1) % 3 + 1;
  int num = mesh_face_edges1(M, face->axis, face->id, edge_ids);
  edges[0] = edge_ids[0] == -1 ? NULL : M->edges[a2] + edge_ids[0];
  edges[1] = edge_ids[1] == -1 ? NULL : M->edges[a2] + edge_ids[1];
  edges[2] = edge_ids[2] == -1 ? NULL : M->edges[a3] + edge_ids[2];
  edges[3] = edge_ids[3] == -1 ? NULL : M->edges[a3] + edge_ids[3];
  return num;
}



int mesh_face_cells(struct mesh *M, struct mesh_face *face, struct mesh_cell *cells[2])
{
  int cell_ids[4];
  int num = mesh_face_cells1(M, face->axis, face->id, cell_ids);
  cells[0] = cell_ids[0] == -1 ? NULL : M->cells + cell_ids[0];
  cells[1] = cell_ids[1] == -1 ? NULL : M->cells + cell_ids[1];
  return num;
}



int mesh_edge_faces(struct mesh *M, struct mesh_edge *edge, struct mesh_face *faces[4])
{
  int face_ids[4];
  int a2 = (edge->axis + 1 - 1) % 3 + 1;
  int a3 = (edge->axis + 2 - 1) % 3 + 1;
  int num = mesh_edge_faces1(M, edge->axis, edge->id, face_ids);
  faces[0] = face_ids[0] == -1 ? NULL : M->faces[a2] + face_ids[0];
  faces[1] = face_ids[1] == -1 ? NULL : M->faces[a2] + face_ids[1];
  faces[2] = face_ids[2] == -1 ? NULL : M->faces[a3] + face_ids[2];
  faces[3] = face_ids[3] == -1 ? NULL : M->faces[a3] + face_ids[3];
  return num;
}



int mesh_edge_cells(struct mesh *M, struct mesh_edge *edge, struct mesh_cell *cells[4])
{
  int cell_ids[4];
  int num = mesh_edge_cells1(M, edge->axis, edge->id, cell_ids);
  cells[0] = cell_ids[0] == -1 ? NULL : M->cells + cell_ids[0];
  cells[1] = cell_ids[1] == -1 ? NULL : M->cells + cell_ids[1];
  cells[2] = cell_ids[2] == -1 ? NULL : M->cells + cell_ids[2];
  cells[3] = cell_ids[3] == -1 ? NULL : M->cells + cell_ids[3];
  return num;
}




#ifdef __MAIN__
int main()
{
  struct mesh Mesh;

  Mesh.cells_shape[1] = 8;
  Mesh.cells_shape[2] = 8;
  Mesh.cells_shape[3] = 8;
  mesh_allocate(&Mesh);

  int index[4] = { 0, 4, 4, 4 };
  int cell_id;
  int face_id;
  int face_ids[6];
  int edge_ids[6];

  cell_id = multi_to_linear(Mesh.cells_shape, index);

  mesh_cell_faces1(&Mesh, cell_id, face_ids);

  linear_to_multi(face_ids[0], Mesh.faces_shape[1], index); printf("x0: %d %d %d\n", index[1], index[2], index[3]);
  linear_to_multi(face_ids[1], Mesh.faces_shape[1], index); printf("x1: %d %d %d\n", index[1], index[2], index[3]);
  linear_to_multi(face_ids[2], Mesh.faces_shape[2], index); printf("y0: %d %d %d\n", index[1], index[2], index[3]);
  linear_to_multi(face_ids[3], Mesh.faces_shape[2], index); printf("y1: %d %d %d\n", index[1], index[2], index[3]);
  linear_to_multi(face_ids[4], Mesh.faces_shape[3], index); printf("z0: %d %d %d\n", index[1], index[2], index[3]);
  linear_to_multi(face_ids[5], Mesh.faces_shape[3], index); printf("z1: %d %d %d\n", index[1], index[2], index[3]);

  face_id = multi_to_linear(Mesh.faces_shape[3], index);
  mesh_face_edges1(&Mesh, 3, face_id, edge_ids);

  linear_to_multi(edge_ids[0], Mesh.edges_shape[1], index); printf("%d %d %d\n", index[1], index[2], index[3]);
  linear_to_multi(edge_ids[1], Mesh.edges_shape[1], index); printf("%d %d %d\n", index[1], index[2], index[3]);
  linear_to_multi(edge_ids[2], Mesh.edges_shape[2], index); printf("%d %d %d\n", index[1], index[2], index[3]);
  linear_to_multi(edge_ids[3], Mesh.edges_shape[2], index); printf("%d %d %d\n", index[1], index[2], index[3]);

  mesh_deallocate(&Mesh);

  return 0;
}
#endif /* __MAIN__ */
