CFLAGS =  -g0 -c
  
OBJFILES = control.o dslash_gn.o hamiltonian.o invert.o io_help.o \
	io_sigma.o kramer.o layout.o measurement.o piup.o setup.o

default:        
	echo TYPE make FILENAME.

gnkmc: $(OBJFILES)
	gcc -O3 $(OBJFILES) -lm -lgsl -lgslcblas -o gnkmc

control.o: control.c
	gcc $(CFLAGS) control.c

dslash_gn.o: dslash_gn.c
	gcc $(CFLAGS) dslash_gn.c

hamiltonian.o: hamiltonian.c
	gcc $(CFLAGS) hamiltonian.c

invert.o: invert.c
	gcc $(CFLAGS) invert.c

io_help.o: io_help.c
	gcc $(CFLAGS) io_help.c

io_sigma.o: io_sigma.c
	gcc $(CFLAGS) io_sigma.c

kramer.o: kramer.c
	gcc $(CFLAGS) kramer.c

layout.o: layout.c
	gcc $(CFLAGS) layout.c

measurement.o: measurement.c
	gcc $(CFLAGS) measurement.c
piup.o: piup.c
	gcc $(CFLAGS) piup.c

setup.o: setup.c
	gcc $(CFLAGS) setup.c

clean:
	rm -f *.o
