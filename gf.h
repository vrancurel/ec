
extern size_t sizew(size_t size);
extern size_t freadw(void *ptr, FILE *stream);
extern size_t fwritew(const void *ptr, FILE *stream);
extern int check_w();
extern int setup_tables();
extern void dump_tables();
extern int gmul(int a, int b);
extern int gdiv(int a, int b);
extern int gexp(int a, int b);
extern void utest();
