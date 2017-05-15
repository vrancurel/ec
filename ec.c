/**
 * @file   ec.c
 * @author vr
 * @date   Tue Feb 28 13:02:01 2017
 * 
 * @brief  Reed-Solomon for Fault-Tolerance in RAID-like Systems
 *         Straightforward implementation of James S. Plank technical report:
 *         http://www.cs.utk.edu/~plank/plank/papers/CS-96-332.html
 */

#include "ec.h"

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
    fprintf(stderr, "encoding matrix:\n");
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
  
  for (i = 0;i < sizew(size);i++) {
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
  
  for (i = 0;i < sizew(size);i++) {
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

