
#include "ec.h"

void vec_zero(t_vec *vec)
{
  memset(vec->mem, 0, sizeof (int) * vec->n);
}

t_vec *vec_xcalloc(u_int n)
{
  t_vec *vec;

  vec = xmalloc(sizeof (*vec));
  vec->n = n;
  vec->mem = xmalloc(sizeof (int) * n);
  vec_zero(vec);
  return vec;
}

void vec_free(t_vec *vec)
{
  if (vec) {
    free(vec->mem);
    free(vec);
  }
}

void vec_dump(t_vec *vec)
{
  int i;

  for (i = 0;i < vec->n;i++)
    printf("%d\n", VEC_ITEM(vec, i));
}

