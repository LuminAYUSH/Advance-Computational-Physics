CFLAGS =  -g0 -c

OBJFILES = gn_hmc.o config.o setup_gn.o layout.o gather.o \
	hmc.o hamil.o piup.o matp2p.o matp2d.o matd2p.o matd2d.o \
	cg_md.o randomize.o ran2.o gasdev.o average_sigma.o\
	autocorel.o gauss.o propagator.o cg_prop.o

default:	
	echo TYPE make FILENAME.

gn_hmc31: $(OBJFILES)
	cc -O3 $(OBJFILES) -lm -o /user1/de/mcrun/gn_hmc31

gn_hmc.o: gn_hmc.c lattice_gn.h
	cc $(CFLAGS) gn_hmc.c

config.o: config.c lattice_gn.h
	cc $(CFLAGS) config.c

setup_gn.o: setup_gn.c lattice_gn.h
	cc $(CFLAGS) setup_gn.c

layout.o: layout.c lattice_gn.h
	cc $(CFLAGS) layout.c

gather.o: gather.c lattice_gn.h
	cc $(CFLAGS) gather.c

hmc.o: hmc.c lattice_gn.h
	cc $(CFLAGS) hmc.c

hamil.o: hamil.c lattice_gn.h
	cc $(CFLAGS) hamil.c

piup.o: piup.c lattice_gn.h
	cc $(CFLAGS) piup.c

matp2p.o: matp2p.c lattice_gn.h
	cc $(CFLAGS) matp2p.c

matp2d.o: matp2d.c lattice_gn.h
	cc $(CFLAGS) matp2d.c

matd2p.o: matd2p.c lattice_gn.h
	cc $(CFLAGS) matd2p.c

matd2d.o: matd2d.c lattice_gn.h
	cc $(CFLAGS) matd2d.c

cg_md.o: cg_md.c lattice_gn.h
	cc $(CFLAGS) cg_md.c

randomize.o: randomize.c 
	cc $(CFLAGS) randomize.c

ran2.o: ran2.c lattice_gn.h
	cc $(CFLAGS) ran2.c

gasdev.o: gasdev.c lattice_gn.h
	cc $(CFLAGS) gasdev.c

average_sigma.o: average_sigma.c lattice_gn.h
	cc $(CFLAGS) average_sigma.c

autocorel.o: autocorel.c lattice_gn.h
	cc $(CFLAGS) autocorel.c

gauss.o: gauss.c lattice_gn.h
	cc $(CFLAGS) gauss.c

propagator.o: propagator.c lattice_gn.h
	cc $(CFLAGS) propagator.c

cg_prop.o: cg_prop.c lattice_gn.h
	cc $(CFLAGS) cg_prop.c

clean:
	rm -f *.o

