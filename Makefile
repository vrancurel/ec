CFLAGS = -Werror -Wall -g
LDFLAGS =

PROGS = ecgf4 ecgf8 ecgf16

all: $(PROGS)

ecgf4: ecgf4.o
	cc -o ecgf4 ecgf4.o $(LDFLAGS)

ecgf4.o: ec.c
	cc -o ecgf4.o -c ec.c $(CFLAGS) -DW=4

ecgf8: ecgf8.o
	cc -o ecgf8 ecgf8.o $(LDFLAGS)

ecgf8.o: ec.c
	cc -o ecgf8.o -c ec.c $(CFLAGS) -DW=8

ecgf16: ecgf16.o
	cc -o ecgf16 ecgf16.o $(LDFLAGS)

ecgf16.o: ec.c
	cc -o ecgf16.o -c ec.c $(CFLAGS) -DW=16

clean:
	rm -f $(PROGS) *.o
