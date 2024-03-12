#ifndef _LATTICE_H
#define _LATTICE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "defines.h"
#include "params.h"

/* definition of globals */
#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

typedef struct{ double f[4]; } psferm;
typedef struct {

	short x,t;	// coords for a site
	int index;	// coord index for a site
	int sign;	// staggered phase

	double mom,pi;// momenta
	double sigma;	// auxillary field = g^2(\psi_bar \psi)
	double phi;	// working storage for sigma
	double rho;	// gaussian noise

	psferm eta,zeta;// pseudofermionic fields
	psferm omega,xi;// pseudofermionic fields

	double r,p,mp,mmp;	// vectors for conjugate gradient
	
} site;

EXTERN site *lattice;

/* Global Lattice Parameters */
EXTERN int Nx,Nt,N;		// N=Nx=Nt
EXTERN int Volume;		// Nx*Nt=N^2
EXTERN int Nf;			// # of flavors
EXTERN double G;		// coupling
EXTERN double Step;		// MD time step
EXTERN int Warmup;		// thermalization loop
EXTERN int NMdtraj;		// # of molecular dynamics trajectory
EXTERN int Nmeas,Measlen;	// max number & interval of measurement
EXTERN int CgiterH,CgiterP;	// for Hamil and MD (piup)
EXTERN double ResidueH,ResidueP;// for Hamil and MD (piup)

EXTERN int sites_on_node;       // = volume

EXTERN int readflag,saveflag;	// fresh or reload / forget or save
EXTERN char read_binary_sigma[MAXFILENAME];
EXTERN char save_binary_sigma[MAXFILENAME];

EXTERN double *av_sigma;	// storage for average sigma

EXTERN int iseed;               // random number seed
EXTERN int junk_id;             // junk number, for checking

EXTERN const gsl_rng_type *T;
EXTERN gsl_rng *r;

/* pointers for accessing and gathering fields */
#define N_POINTERS 4
EXTERN char **gen_pt;
//EXTERN char **gen_pt[N_POINTERS];

// define macros, we will not be using MILC macros for 2-dim
typedef int field_offset;
#define F_OFFSET(a) ((field_offset)(((char *)&(lattice[0].a))- ((char *)&(lattice[0]))))
#define F_PT(site,fo) ((char *)(site) + (fo))

#endif /* _LATTICE_H */

