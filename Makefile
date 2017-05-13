CFLAGS = -Werror -Wall -g
LDFLAGS =

PROGS = ecgf4 ecgf8 ecgf16

COMMON_OBJS = ec.o main.o mat.o misc.o vec.o

all: $(PROGS)

ecgf4: gf4.o $(COMMON_OBJS)
	cc -o ecgf4 gf4.o $(COMMON_OBJS) $(LDFLAGS)

gf4.o: gf.c
	cc -o gf4.o -c gf.c $(CFLAGS) -DW=4

ecgf8: gf8.o $(COMMON_OBJS)
	cc -o ecgf8 gf8.o $(COMMON_OBJS) $(LDFLAGS)

gf8.o: ec.c
	cc -o gf8.o -c gf.c $(CFLAGS) -DW=8

ecgf16: gf16.o $(COMMON_OBJS)
	cc -o ecgf16 gf16.o $(COMMON_OBJS) $(LDFLAGS)

gf16.o: ec.c
	cc -o gf16.o -c gf.c $(CFLAGS) -DW=16

clean:
	rm -f $(PROGS) *.o
