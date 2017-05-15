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

#include "vec.h"
#include "mat.h"
#include "misc.h"
#include "gf.h"
#include "main.h"

extern void create_coding_files(char *prefix, t_mat *mat);
extern int repair_data_files(char *prefix, t_mat *mat);
