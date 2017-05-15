
#include "ec.h"

void mat_zero(t_mat *mat)
{
  memset(mat->mem, 0, sizeof (int) * mat->n_rows * mat->n_cols);
}

t_mat *mat_xcalloc(u_int n_rows, u_int n_cols)
{
  t_mat *mat;

  mat = xmalloc(sizeof (*mat));
  mat->n_rows = n_rows;
  mat->n_cols = n_cols;
  mat->mem = xmalloc(sizeof (int) * n_rows * n_cols);
  mat_zero(mat);
  return mat;
}

void mat_free(t_mat *mat)
{
  if (mat) {
    free(mat->mem);
    free(mat);
  }
}

void mat_dump(t_mat *mat)
{
  int i, j;
  
  for (i = 0;i < mat->n_rows;i++) {
    for (j = 0;j < mat->n_cols;j++) {
      printf("%d ", MAT_ITEM(mat, i, j));
    }
    printf("\n");
  }
}

t_mat *mat_vandermonde(u_int n_rows, u_int n_cols)
{
  t_mat *mat;
  int i, j;
  
  mat = mat_xcalloc(n_rows, n_cols);
  for (i = 0;i < n_rows;i++) {
    for (j = 0;j < n_cols;j++) {
      MAT_ITEM(mat, i, j) = gexp(j + 1, i); 
    }
  }
  return mat;
}

t_mat *mat_cauchy(u_int n_rows, u_int n_cols)
{
  t_mat *mat;
  int i, j;

  mat = mat_xcalloc(n_rows, n_cols);
  for (i = 0;i < n_rows;i++) {
    for (j = 0;j < n_cols;j++) {
      MAT_ITEM(mat, i, j) = gdiv(1, (i ^ (j + n_rows)));
    }
  }

  /* do optimise */
  // convert 1st row to all 1s
  for (j = 0;j < n_cols;j++) {
    for (i = 0;i < n_rows;i++) {
      MAT_ITEM(mat, i, j) = gdiv(MAT_ITEM(mat, i, j), MAT_ITEM(mat, 0, j));
    }
  }
  // convert 1st element of each row to 1
  for (i = 1;i < n_rows;i++) {
    for (j = 0;j < n_cols;j++) {
      MAT_ITEM(mat, i, j) = gdiv(MAT_ITEM(mat, i, j), MAT_ITEM(mat, i, 0));
    }
  }

  return mat;
}

int mat_check_row_is_identity(t_mat *mat, int row)
{
  int j;

  for (j = 0;j < mat->n_cols;j++) {
    if (MAT_ITEM(mat, row, j) != ((j == row) ? 1 : 0))
      return 0;
  }
  return 1;
}

void mat_swap_cols(t_mat *mat, int c1, int c2)
{
  int i;

  assert(c1 < mat->n_cols && c2 < mat->n_cols);
  for (i = 0;i < mat->n_rows;i++) {
    int val = MAT_ITEM(mat, i, c1);
    MAT_ITEM(mat, i, c1) = MAT_ITEM(mat, i, c2);
    MAT_ITEM(mat, i, c2) = val;
  }
}

/*
 * transform c_i into f_i_i_minus_1 * c_i 
 */
void mat_transform1(t_mat *tmp, int i)
{
  int k;
  int f_minus_1 = gdiv(1, MAT_ITEM(tmp, i, i));

  for (k = 0;k < tmp->n_rows;k++) {
    MAT_ITEM(tmp, k, i) = gmul(f_minus_1, MAT_ITEM(tmp, k, i));
  }
}

/*
 * transform c_j into c_j - f_i_j * c_i 
 */
void mat_transform2(t_mat *tmp, int i, int j)
{
  int k;
  int f_i_j = MAT_ITEM(tmp, i, j);

  for (k = 0;k < tmp->n_rows;k++) {
    MAT_ITEM(tmp, k, j) = MAT_ITEM(tmp, k, j) ^ gmul(f_i_j, MAT_ITEM(tmp, k, i));
  }
}

t_mat *mat_vandermonde_correct(u_int n_rows, u_int n_cols)
{
  t_mat *mat, *tmp;
  int i, j, dim;
  
  dim = n_rows + n_cols;
  tmp = mat_xcalloc(dim, n_cols);
  for (i = 0;i < dim;i++) {
    for (j = 0;j < n_cols;j++) {
      MAT_ITEM(tmp, i, j) = gexp(i, j); 
    }
  }

  /* perform transformations to get the identity matrix on the top rows */
  i = 0;
  while (i < n_cols) {
    
    if (mat_check_row_is_identity(tmp, i)) {
      i++;
      continue ;
    }
    
    //this case is mentionned in the paper but cannot happen
    /* if (0 == MAT_ITEM(tmp, i, i)) {
       for (j = i + 1;j < tmp->n_cols;j++) {
       if (0 != MAT_ITEM(tmp, i, j)) {
       mat_swap_cols(tmp, i, j);
       continue ;
       }
       } */
    
    //check if f_i_i == 1
    if (1 != MAT_ITEM(tmp, i, i)) {
      //check for inverse since f_i_i != 0
      mat_transform1(tmp, i);
    }
    
    //now f_i_i == 1
    for (j = 0;j < tmp->n_cols;j++) {
      if (i != j) {
        if (0 != MAT_ITEM(tmp, i, j)) {
          mat_transform2(tmp, i, j);
        }
      }
    }
    
    i++;
  }

  mat = mat_xcalloc(n_rows, n_cols);

  //copy last n_rows rows of tmp into mat
  for (i = 0;i < n_rows;i++) {
    for (j = 0;j < n_cols;j++) {
      MAT_ITEM(mat, i, j) = MAT_ITEM(tmp, n_cols + i, j);
    }
  }

  free(tmp);
  return mat;
}

void mat_inv(t_mat *mat)
{
  t_mat *aug;
  int dim, i, j, k, tpos, tval, r;

  assert(mat->n_rows == mat->n_cols);
  dim = mat->n_rows;
  aug = mat_xcalloc(dim, dim * 2);

  for (i = 0;i < dim;i++) {
    for (j = 0;j < dim;j++) {
      MAT_ITEM(aug, i, j) = MAT_ITEM(mat, i, j);
    }
  }

  for (i = 0;i < dim;i++) {
    for (j = dim;j < 2*dim;j++) {
      if (i == (j % dim))
        MAT_ITEM(aug, i, j) = 1;
      else
        MAT_ITEM(aug, i, j) = 0;
    }
  }

  /* using gauss-jordan elimination */
  for (j = 0;j < dim;j++) {
    tpos = j;
    
    /* finding maximum jth column element in last (dimension-j) rows */
    for (i = j+1;i < dim;i++) {
      if (MAT_ITEM(aug, i, j) > MAT_ITEM(aug, tpos, j))
        tpos = i;
    }
    
    /* swapping row which has maximum jth column element */
    if (tpos != j) {
      for (k=0;k < 2*dim;k++) {
        tval = MAT_ITEM(aug, j, k);
        MAT_ITEM(aug, j, k) = MAT_ITEM(aug, tpos, k);
        MAT_ITEM(aug, tpos, k) = tval;
      }
    }
    
    /* performing row operations to form required identity matrix out of the input matrix */
    for (i = 0;i < dim;i++) {
      if (i != j) {
        r = MAT_ITEM(aug, i, j);
        for (k = 0;k < 2*dim;k++) {
          MAT_ITEM(aug, i, k) ^= gmul(gdiv(MAT_ITEM(aug, j, k), MAT_ITEM(aug, j, j)), r);
        }
      } else {
        r = MAT_ITEM(aug, i, j);
        for (k = 0;k < 2*dim;k++) {
          MAT_ITEM(aug, i, k) = gdiv(MAT_ITEM(aug, i, k), r);
        }
      }
    }
  }

  for (i = 0;i < dim;i++) {
    for (j = 0;j < dim;j++) {
      MAT_ITEM(mat, i, j) = MAT_ITEM(aug, i, dim + j);
    }
  }

  mat_free(aug);
}

void mat_mult(t_vec *output, t_mat *a, t_vec *b)
{
  int i, j;

  assert(output->n == b->n);
  assert(output->n == a->n_cols);
  for (i = 0;i < a->n_rows;i++) {
    for (j = 0;j < a->n_cols;j++) {
      int x = gmul(MAT_ITEM(a, i, j), VEC_ITEM(b, j));
      if (0 == j)
        VEC_ITEM(output, i) = x;
      else
        VEC_ITEM(output, i) ^= x;
    }
  }
}

