
#include "ec.h"

int vflag = 0;

void xusage()
{
  fprintf(stderr,
          "Usage: erasure [-n n_data][-m n_coding][-s (use cauchy instead of vandermonde)][-p prefix][-v (verbose)] -c (encode) | -r (repair) | -u (utest)\n");
  exit(1);
}

int main(int argc, char **argv)
{
  int n_data, n_coding, opt;
  t_mat *mat;
  char *prefix = NULL;
  int cflag = 0;
  int rflag = 0;
  int uflag = 0;
  int sflag = 0;

  n_data = n_coding = -1;
  prefix = NULL;
  while ((opt = getopt(argc, argv, "n:m:p:scruv")) != -1) {
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
    case 's':
      sflag = 1;
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

  if (0 != check_w(n_data + n_coding)) {
    fprintf(stderr, "Number of fragments is too big compared to Galois field size\n");
    exit(1);
  }

  setup_tables();
  //dump_tables();

  if (uflag) {
    utest();
    goto end;
  }

  if (-1 == n_data || -1 == n_coding || NULL == prefix)
    xusage();

  if (sflag) {
    mat = mat_cauchy(n_coding, n_data);
  } else {
    mat = mat_vandermonde_correct(n_coding, n_data);
  }
  if (vflag)
    mat_dump(mat);

  if (rflag) {
    if (0 != repair_data_files(prefix, mat)) {
      exit(1);
    }
  }
  create_coding_files(prefix, mat);

  mat_free(mat);

 end:
  free(prefix);
  return 0;
}
