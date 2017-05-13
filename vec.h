
typedef struct s_vec
{
  u_int n;
  int *mem;
#define VEC_ITEM(vec, i) ((vec)->mem[(i)])
} t_vec;

extern void vec_zero(t_vec *vec);
extern t_vec *vec_xcalloc(u_int n);
extern void vec_free(t_vec *vec);
extern void vec_dump(t_vec *vec);
