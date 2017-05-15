
#include "ec.h"

#ifndef W
# error "please define W"
#endif
#define NW (1 << W)   /* In other words, NW equals 2 to the w-th power */
#if W == 4
# define SIZEW(size) ((size)*2)
#elif W == 8
# define SIZEW(size) (size)
#elif W == 16
# define SIZEW(size) ((size)/2)
#endif

u_int prim_poly_4 = 023;
u_int prim_poly_8 = 0435;
u_int prim_poly_16 = 0210013;
unsigned short *gflog = NULL;
unsigned short *gfilog = NULL;

size_t sizew(size_t size)
{
  return SIZEW(size);
}

size_t freadw(void *ptr, FILE *stream)
{
#if W == 4
  return 0; //TBD
#elif W == 8
  return fread(ptr, 1, 1, stream);
#elif W == 16
  return fread(ptr, 2, 1, stream);
#endif
}

size_t fwritew(const void *ptr, FILE *stream)
{
#if W == 4
  return 0; //TBD
#elif W == 8
  return fwrite(ptr, 1, 1, stream);
#elif W == 16
  return fwrite(ptr, 2, 1, stream);
#endif
}

int check_w(int n)
{
  if (n > NW)
    return -1;
  return 0;
}

int setup_tables()
{
  u_int b, log, x_to_w, prim_poly;

  switch(W) {
  case 4:  prim_poly = prim_poly_4;  break;
  case 8:  prim_poly = prim_poly_8;  break;
  case 16: prim_poly = prim_poly_16; break;
  default: return -1;
  }
  x_to_w = 1 << W;
  gflog  = (unsigned short *) xmalloc (sizeof(unsigned short) * x_to_w);
  gfilog = (unsigned short *) xmalloc (sizeof(unsigned short) * x_to_w);
  b = 1;
  for (log = 0; log < x_to_w-1; log++) {
    gflog[b] = (unsigned short) log;
    gfilog[log] = (unsigned short) b;
    b = b << 1;
    if (b & x_to_w) b = b ^ prim_poly;
  }
  return 0; 
}

void dump_tables()
{
  u_int log, x_to_w;

  x_to_w = 1 << W;
  for (log = 0; log < x_to_w-1; log++) {
    printf("%d ", gflog[log]);
  }
  printf("\n");
  for (log = 0; log < x_to_w-1; log++) {
    printf("%d ", gfilog[log]);
  }
  printf("\n");
}

int gmul(int a, int b)
{
  int sum_log;
  if (a == 0 || b == 0) return 0;
  sum_log = gflog[a] + gflog[b];
  if (sum_log >= NW-1) sum_log -= NW-1;
  return gfilog[sum_log];
}

int gdiv(int a, int b)
{
  int diff_log;
  if (a == 0) return 0;
  if (b == 0) return -1;
  diff_log = gflog[a] - gflog[b];
  if (diff_log < 0) diff_log += NW-1;
  return gfilog[diff_log];
}

int gexp(int a, int b)
{
  int r, i;

  if (0 == b)
    return 1;

  if (1 == b)
    return a;

  r = a;
  for (i = 1; i < b;i++) {
    r = gmul(r, a);
  }

  return r;
}

void utest()
{
#if W == 4
  assert(gmul(3, 7) == 9);
  assert(gmul(13, 10) == 11);  
  assert(gdiv(13, 10) == 3);
  assert(gdiv(3, 7) == 10);
  /* non-MDS vandermonde matrix */
  t_mat *mat = mat_vandermonde(3, 3);
  t_vec *vec = vec_xcalloc(3);
  VEC_ITEM(vec, 0) = 3;
  VEC_ITEM(vec, 1) = 13;
  VEC_ITEM(vec, 2) = 9;
  t_vec *output = vec_xcalloc(3);
  mat_mult(output, mat, vec);
  assert(VEC_ITEM(output, 0) == 7);
  assert(VEC_ITEM(output, 1) == 2);
  assert(VEC_ITEM(output, 2) == 9);
  VEC_ITEM(vec, 0) = 3;
  VEC_ITEM(vec, 1) = 1;
  VEC_ITEM(vec, 2) = 9;
  mat_mult(output, mat, vec);
  assert(VEC_ITEM(output, 0) == 11);
  assert(VEC_ITEM(output, 1) == 9);
  assert(VEC_ITEM(output, 2) == 12);
  MAT_ITEM(mat, 0, 0) = 1;
  MAT_ITEM(mat, 0, 1) = 0;
  MAT_ITEM(mat, 0, 2) = 0;
  MAT_ITEM(mat, 1, 0) = 1;
  MAT_ITEM(mat, 1, 1) = 1;
  MAT_ITEM(mat, 1, 2) = 1;
  MAT_ITEM(mat, 2, 0) = 1;
  MAT_ITEM(mat, 2, 1) = 2;
  MAT_ITEM(mat, 2, 2) = 3;
  VEC_ITEM(vec, 0) = 3;
  VEC_ITEM(vec, 1) = 11;
  VEC_ITEM(vec, 2) = 9;
  mat_inv(mat);
  mat_mult(output, mat, vec);
  assert(VEC_ITEM(output, 1) == 1);
  assert(VEC_ITEM(output, 2) == 9);
  vec_free(output);
  vec_free(vec);
  mat_free(mat);
  mat = mat_vandermonde_correct(3, 3);
  mat_free(mat);
  mat = mat_cauchy(3, 3);
  mat_free(mat);
#elif W == 8
  assert(gmul(3, 7) == 9);
  assert(gmul(13, 10) == 114);
  assert(gdiv(13, 10) == 40);
  assert(gdiv(3, 7) == 211);
#else
  //TBD
#endif
}

