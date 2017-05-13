
#include "ec.h"

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

