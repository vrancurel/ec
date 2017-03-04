/**
 * @file   ec.c
 * @author vr
 * @date   Tue Feb 28 13:02:01 2017
 * 
 * @brief  Reed-Solomon for Fault-Tolerance in RAID-like Systems
 *         Straightforward implementation of James S. Plank technical report:
 *         http://www.cs.utk.edu/~plank/plank/papers/CS-96-332.html
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <getopt.h>
#include <unistd.h>

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

typedef struct s_vec
{
  u_int n;
  int *mem;
#define VEC_ITEM(vec, i) ((vec)->mem[(i)])
} t_vec;

u_int prim_poly_4 = 023;
u_int prim_poly_8 = 0435;
u_int prim_poly_16 = 0210013;
unsigned short *gflog = NULL;
unsigned short *gfilog = NULL;

int vflag = 0;

void xperror(char *str)
{
  perror(str);
  exit(1);
}

void xerrormsg(char *str1, char *str2)
{
  fprintf(stderr, "%s %s: %s\n", str1, str2, strerror(errno));
  exit(1);
}

void xmsg(char *str1, char *str2)
{
  fprintf(stderr, "%s %s\n", str1, str2);
  exit(1);
}

void *xmalloc(size_t size)
{
  void *p;

  if (NULL == (p = malloc(size)))
    xperror("malloc");
  return p;
}

char *xstrdup(char *str)
{
  char *n;

  if (NULL == (n = strdup(str)))
    xperror("malloc");
  return n;
}

int setup_tables(int w)
{
  u_int b, log, x_to_w, prim_poly;

  switch(w) {
  case 4:  prim_poly = prim_poly_4;  break;
  case 8:  prim_poly = prim_poly_8;  break;
  case 16: prim_poly = prim_poly_16; break;
  default: return -1;
  }
  x_to_w = 1 << w;
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

void dump_tables(int w)
{
  u_int log, x_to_w;

  x_to_w = 1 << w;
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
#elif W == 8
  assert(gmul(3, 7) == 9);
  assert(gmul(13, 10) == 114);
  assert(gdiv(13, 10) == 40);
  assert(gdiv(3, 7) == 211);
#else
  //TBD
#endif
}

/** 
 * (re-)create missing prefix.c1 ... cm files acc/to Vandermonde matrix
 * 
 * @param prefix prefix of files
 * @param mat Vandermonde matrix
 */
void create_coding_files(char *prefix, t_mat *mat)
{
  int i, j;
  FILE *d_files[mat->n_cols];
  FILE *c_files[mat->n_rows];
  char filename[1024];
  struct stat stbuf;
  size_t size = -1, wc;
  t_vec *words = NULL;
  t_vec *output = NULL;

  if (vflag) {
    fprintf(stderr, "vandermonde matrix:\n");
    mat_dump(mat);
  }
  
  for (i = 0;i < mat->n_cols;i++) {
    snprintf(filename, sizeof (filename), "%s.d%d", prefix, i);
    if (NULL == (d_files[i] = fopen(filename, "r")))
      xerrormsg("error opening", filename);
    if (-1 == fstat(fileno(d_files[i]), &stbuf))
      xerrormsg("error stating", filename);
    if (-1 == size)
      size = stbuf.st_size;
    else if (size != stbuf.st_size)
      xmsg("bad size", filename);
  }
  
  for (i = 0;i < mat->n_rows;i++) {
    snprintf(filename, sizeof (filename), "%s.c%d", prefix, i);
    if (NULL == (c_files[i] = fopen(filename, "w")))
      xerrormsg("error opening", filename);
  }
  
  words = vec_xcalloc(mat->n_cols);
  output = vec_xcalloc(mat->n_cols);
  
  for (i = 0;i < SIZEW(size);i++) {
    vec_zero(words);
    for (j = 0;j < mat->n_cols;j++) {
      wc = freadw(&VEC_ITEM(words, j), d_files[j]);
      if (1 != wc)
        xperror("short read data");
    }
    mat_mult(output, mat, words);
    for (j = 0;j < mat->n_rows;j++) {
      wc = fwritew(& VEC_ITEM(output, j), c_files[j]);
      if (1 != wc)
        xperror("short write coding");
    }
  } 
    
  for (i = 0;i < mat->n_cols;i++) {
    fclose(d_files[i]);
  }
  
  for (i = 0;i < mat->n_rows;i++) {
    fclose(c_files[i]);
  }

  vec_free(words);
  vec_free(output);
}

/** 
 * repair data files
 * 
 * @param prefix prefix of files 
 * @param mat 
 */
int repair_data_files(char *prefix, t_mat *mat)
{
  int i, j, k;
  FILE *d_files[mat->n_cols];
  FILE *r_files[mat->n_cols];
  FILE *c_files[mat->n_rows];
  char filename[1024];
  struct stat stbuf;
  size_t size = -1, wc;
  t_mat *a_prime = NULL;
  u_int n_data_ok = 0;
  u_int n_coding_ok = 0;
  u_int n_total_ok;
  t_vec *words = NULL;
  t_vec *output = NULL;
  int ret;
  
  for (i = 0;i < mat->n_cols;i++) {
    snprintf(filename, sizeof (filename), "%s.d%d", prefix, i);
    if (-1 == access(filename, F_OK)) {
      if (vflag)
        fprintf(stderr, "%s is missing\n", filename);
      d_files[i] = NULL;
      if (NULL == (r_files[i] = fopen(filename, "w")))
        xerrormsg("error opening", filename);
    } else {
      r_files[i] = NULL;
      if (NULL == (d_files[i] = fopen(filename, "r")))
        xerrormsg("error opening", filename);
      if (-1 == fstat(fileno(d_files[i]), &stbuf))
        xerrormsg("error stating", filename);
      if (-1 == size)
        size = stbuf.st_size;
      else if (size != stbuf.st_size)
        xmsg("bad size", filename);
      n_data_ok++;
    }
  }
  
  for (i = 0;i < mat->n_rows;i++) {
    snprintf(filename, sizeof (filename), "%s.c%d", prefix, i);
    if (access(filename, F_OK)) {
      if (vflag)
        fprintf(stderr, "%s is missing\n", filename);
      c_files[i] = NULL;
    } else {
      if (NULL == (c_files[i] = fopen(filename, "r")))
        xerrormsg("error opening", filename);
      n_coding_ok++;
    }
  }

  if (n_data_ok == mat->n_cols) {
    ret = 0;
    goto end;
  }

  if (n_coding_ok < (mat->n_cols-n_data_ok)) {
    fprintf(stderr, "too many losses\n");
    ret = -1;
    goto end;
  }

  n_total_ok = n_data_ok + n_coding_ok;
  if (vflag)
    fprintf(stderr, "n_data_ok=%d n_coding_ok=%d\n", n_data_ok, n_coding_ok);

  //generate a_prime
  a_prime = mat_xcalloc(n_total_ok, mat->n_cols);
  //for each data available generate the corresponding identity
  k = 0;
  for (i = 0;i < mat->n_cols;i++) {
    if (NULL != d_files[i]) {
      for (j = 0;j < mat->n_cols;j++) {
        if (i == j)
          MAT_ITEM(a_prime, k, j) = 1;
        else
          MAT_ITEM(a_prime, k, j) = 0;
      }
      k++;
    }
  }
  //finish the matrix with every coding available
  for (i = 0;i < mat->n_rows;i++) {
    if (NULL != c_files[i]) {
      //copy corresponding row in vandermonde matrix
      for (j = 0;j < mat->n_cols;j++) {
        MAT_ITEM(a_prime, k, j) = MAT_ITEM(mat, i, j);
      }
      k++;
      //stop when we have enough codings
      if (mat->n_cols == k)
        break ;
    }
  }

  if (vflag) {
    fprintf(stderr, "rebuild matrix:\n");
    mat_dump(a_prime);
  }

  mat_inv(a_prime);

  //read-and-repair
  words = vec_xcalloc(mat->n_cols);
  output = vec_xcalloc(mat->n_cols);
  
  for (i = 0;i < SIZEW(size);i++) {
    vec_zero(words);
    k = 0;
    for (j = 0;j < mat->n_cols;j++) {
      if (NULL != d_files[j]) {
        wc = freadw(&VEC_ITEM(words, k), d_files[j]);
        if (1 != wc)
          xperror("short read data");
        k++;
      }
    }
    for (j = 0;j < mat->n_rows;j++) {
      if (NULL != c_files[j]) {
        wc = freadw(&VEC_ITEM(words, k), c_files[j]);
        if (1 != wc)
          xperror("short read coding");
        k++;
        //stop when we have enough codings
        if (mat->n_cols == k)
          break ;
      }
    }

    mat_mult(output, a_prime, words);

    for (j = 0;j < mat->n_cols;j++) {
      if (NULL != r_files[j]) {
        wc = fwritew(& VEC_ITEM(output, j), r_files[j]);
        if (1 != wc)
          xperror("short write coding");
      }
    }
  } 
   
  ret = 0;
 end:
  for (i = 0;i < mat->n_cols;i++) {
    if (NULL != d_files[i])
      fclose(d_files[i]);
    if (NULL != r_files[i])
      fclose(r_files[i]);
  }
  
  for (i = 0;i < mat->n_rows;i++) {
    if (NULL != c_files[i])
      fclose(c_files[i]);
  }

  vec_free(words);
  vec_free(output);
  mat_free(a_prime);

  return ret;
}

void xusage()
{
  fprintf(stderr, 
          "Usage: erasure [-n n_data][-m n_coding][-p prefix][-v (verbose)] -c (encode) | -r (repair) | -u (utest)\n");
  exit(1);
}  

int main(int argc, char **argv)
{
  int n_data, n_coding, opt;
  t_mat *vander_mat;
  char *prefix = NULL;
  int cflag = 0;
  int rflag = 0;
  int uflag = 0;
  
  n_data = n_coding = -1;
  prefix = NULL;
  while ((opt = getopt(argc, argv, "n:m:p:cruv")) != -1) {
    switch (opt) {
    case 'v':
      vflag = 1;
      break ;
    case 'u':
      uflag = 1;
      break ;
    case 'c':
      cflag = 1;
      break ;
    case 'r':
      rflag = 1;
      break ;
    case 'n':
      n_data = atoi(optarg);
      break;
    case 'm':
      n_coding = atoi(optarg);
      break;
    case 'p':
      prefix = xstrdup(optarg);
      break;
    default: /* '?' */
      xusage();
    }
  }
  
  if (!(uflag || cflag || rflag))
    xusage();

  setup_tables(W);
  //dump_tables(W);

  if (uflag) {
    utest();
    goto end;
  }

  if (-1 == n_data || -1 == n_coding || NULL == prefix)
    xusage();

  vander_mat = mat_vandermonde_correct(n_coding, n_data);
  if (vflag)
    mat_dump(vander_mat);

  if (rflag) {
    if (0 != repair_data_files(prefix, vander_mat)) {
      exit(1);
    }
  }
  create_coding_files(prefix, vander_mat);

  mat_free(vander_mat);

 end:
  free(prefix);
  return 0;
}
