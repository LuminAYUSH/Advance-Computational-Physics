/* *********************lattice_gn_kr.h********************* */

/* This file contains the defined quantities, global variables -
   to be read from setup.c, and the lattice structure for lattice
   Gross_Neveu model.  						*/

#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

#define PI 3.14159265358979323846
#define COLDLAT 0
#define HOTLAT 1
#define PLUS 1
#define MINUS -1
#define EVEN 0x02
#define ODD 0x01
#define EVENANDODD 0x03
#define LSIZE 16
#define DELTAMAX 100.0
#define KMAX 4

/* Remove the comments when switch = 0 */
#define T_cut 400
#define D_cut 100
#define NOT_cut 1
#define MAXT_cut 400

/* Remove the comments when switch = 1 */
/*#define T_cut 400
#define D_cut 100
#define NOT_cut 1
#define MAXT_cut 400*/ 

/* With the "EVENFIRST" option the even sites are listed contiguously
   in the first part of lattice[], and the odd sites in the last part.
   In this version there must be an equal number of even and odd sites
   in each plane - in other words one of the two shortest dimensions must
   be even.
*/

/*#define EVENFIRST */


EXTERN int nx, nt;		/* lattice dimensions */
EXTERN unsigned int volume;	/* nx * nt */
EXTERN unsigned int nf;		/* # of flavors */
EXTERN double gama;		/* gamma is a tunable parameter */
EXTERN int cgiter1, cgiter2;	/* # of conjugate grad. iterations */
EXTERN long *iseed;		/* random number seed */
EXTERN double g;		/* Gross_Neveu coupling */
EXTERN double step;		/* leap frog step size */
EXTERN int mid;			/* middle of the lattice */
EXTERN double residue1, residue2; /* residue for cgiter1 and 2 */
EXTERN int no_even_sites;	/* # of even sites */
EXTERN int no_odd_sites;	/* # of odd sites */
EXTERN int no_garbage;          /* # of garbage loops */
EXTERN int bin_length;		/* length of each bin in garbage loops */
EXTERN int no_bin;		/* # of bins */
EXTERN int meas_loop;		/* # of measurement loops after garbage */
EXTERN int meas_length;		/* length after which each meas. is made */
EXTERN int prop_length;		/* length after which Acl. calculation 
				   on propagator is made */
EXTERN int no_meas;		/* no_meas = meas_loop/meas_length */
EXTERN int seg_length;		/* length of each segments */
EXTERN int no_a_seg;		/* no. of segments for autocoreln. meas. */
EXTERN int no_prop_seg;		/* no. of segments for Prop-Acl calcu */
EXTERN int kr_it;               /* # of KRAMER sweeps, counting both 
                                   accepted & rejected ones  */
EXTERN int counter;		/* counter for KRAMER loop */
EXTERN int sw_flag;		/* flag to switch action between garbage and
				   autocorln. loops and measurement loops.
				   sw_flag = 0 => garbage & autocorln. loop;
				   sw_flag = 1 => measurement loop; */

/* The lattice structure: each element of lattice (i.e. each site) 
stores all the informations of that site.   */

typedef struct{ double f[4]; }psferm;
typedef struct{
		int x, t;	/* co ordinates of this site */
		int sign;	/* phase for staggered fermions */
		int parity;	/* even or odd site */
		double sigma;	/* real scalar field */
		double phi;	/* temporary space for sigma */
		double mom;	/* canonical momenta of sigma */
		double mom1;	/* canonical momenta of sigma */
		double zeta;	/* store for gaussian values */
		psferm chi;	/* pseudo fermionic field, source vector */
		psferm eta;	/* Gaussian noise, from which chi is derived, 
		                   also used for dest vector etc */
		double p;	/* pseudofermionic vector for CG */
		double r;	/* pseudofermionic vector for CG, 
		                   used for residue */
		double mp;	/* pseudofermionic vector for CG */
		double mmp;	/* pseudofermionic vector for CG */	
	      } site;
EXTERN site *lattice;

EXTERN double *store;		/* storage space for sigma configuration while
			   	   the code is running. */
EXTERN double *conf;		/* storage space for sigma config. last
				   generated. */
EXTERN double *con;		/* storage space for sigma config. generated 
				   last but one. */
EXTERN double *garbage;		/* storage for <sigma> during garbage loops */
EXTERN double *ac_store;	/* storage for <sigma> during autocrln. loops*/
EXTERN double *ac_prop[LSIZE];	/* storage for propagator during autocrln. 
				   calculation loops*/
EXTERN double *bin_av;		/* storage space for bin averages */
EXTERN double *psi;		/* storage space for <psi_bar-psi> */
EXTERN double *psi_acl;		/* storage space for <psi_bar_psi> for Acl. 
				   calculation for each confign. */ 
EXTERN double *G_prop;		/* storage for propagator */
EXTERN double *G_store[LSIZE];	/* storage for propagator */
EXTERN double *G_temp[LSIZE];	/* storage for propagator for each confign. */ 
EXTERN double *prop;		/* working storage for propagator */
EXTERN double *tprop;		/* working storage for propagator */
EXTERN double *T_int[NOT_cut][MAXT_cut];   /* storage for t_int for various 
					      segments */
EXTERN double *T_int_prop[NOT_cut][MAXT_cut][LSIZE];      /* storage for t_int for
							  Prop-Acl */

EXTERN char ** gen_pt;

/* macros for "field offset" and "field pointer", used when fields
   are passed as arguments to subroutines.
   fo = F_OFFSET(field) converts the name of the field into an
        integer, gives the offset in the site structure of the
        named field.
   address = F_PT(sitepointer, field offset) converts this "fo" 
             integer back into an address at a given site.   */

typedef int field_offset;
#define F_OFFSET(a) ((field_offset)(((char *)&(lattice[0].a))- ((char *)&(lattice[0]))))
#define F_PT(site, fo) ((char *)(site) + (fo))

/* macros to loop over sites. 
   Usage:
	 int i; 	 --> index of the sites, used as a counter here 
	 site *s;	 --> pointer to the current site 
	 FORALLSITES(i,s){
			   Body of the for loop;
			  }                  */ 
			  
#define FORALLSITES(i,s) \
    for(i=0,s=lattice;i<volume;i++,s++)
#ifdef EVENFIRST
#define FOREVENSITES(i,s) \
    for(i=0,s=lattice;i<no_even_sites;i++,s++)
#define FORODDSITES(i,s) \
    for(i=no_even_sites,s= &(lattice[i]);i<volume;i++,s++)
#define FORSOMEPARITY(i,s,choice) \
    for( i=((choice)==ODD ? no_even_sites : 0 ),  \
    s= &(lattice[i]); \
    i< ( (choice)==EVEN ? no_even_sites : volume); \
    i++,s++)
/**#define FORSOMEPARITY(i,s,choice) \
    for( i=((choice)==ODD ? no_even_sites : 0 ),  \
    s= &(lattice[i]), \
    last_in_loop = ((choice)==EVEN ? no_odd_sites : volume); \
    i< last_in_loop; \
    i++,s++)**/
#else
#define FOREVENSITES(i,s) \
    for(i=0,s=lattice;i<volume;i++,s++)if(s->parity==EVEN)
#define FORODDSITES(i,s) \
    for(i=0,s=lattice;i<volume;i++,s++)if(s->parity==ODD)
#define FORSOMEPARITY(i,s,choice) \
    for(i=0,s=lattice;i<volume;i++,s++)if( (s->parity & (choice)) != 0)
#endif	/* end ifdef EVENFIRST */


	
	
/* Directions */	
/* They must go from 0 to 3 in 2 dimensions because they will be used
   to index arrays */
#define XUP 0
#define TUP 1
#define TDN 2
#define XDN 3
#define NDIRS 4                   /* no. of dirns. */
#define OPP_DIR(dir)  (7-(dir))   /* opp. dirn. */   

/* ************ end of lattice_gn_kr.h ************ */
