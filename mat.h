
/*
 * mat[n_rows][n_cols]:
 * mat[0][0] mat[0][1] ... mat[0][n_cols]
 * mat[1][0] ...
 * mat[n_rows-1] ...       mat[n_rows][n_cols]
 */
typedef struct s_mat
{
  u_int n_rows;
  u_int n_cols;
  int *mem;
#define MAT_ITEM(mat, i, j) ((mat)->mem[(i) * (mat)->n_cols + (j)])
} t_mat;

extern void mat_zero(t_mat *mat);
extern t_mat *mat_xcalloc(u_int n_rows, u_int n_cols);
extern void mat_free(t_mat *mat);
extern void mat_dump(t_mat *mat);
extern t_mat *mat_vandermonde(u_int n_rows, u_int n_cols);
extern t_mat *mat_cauchy(u_int n_rows, u_int n_cols);
extern int mat_check_row_is_identity(t_mat *mat, int row);
extern void mat_swap_cols(t_mat *mat, int c1, int c2);
extern void mat_transform1(t_mat *tmp, int i);
extern void mat_transform2(t_mat *tmp, int i, int j);
extern t_mat *mat_vandermonde_correct(u_int n_rows, u_int n_cols);
extern void mat_inv(t_mat *mat);
extern void mat_mult(t_vec *output, t_mat *a, t_vec *b);
