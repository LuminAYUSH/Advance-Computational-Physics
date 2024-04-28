/****************** gather.c **************************************/
/* Routines for accessing (gathering) neighbors */
/*  Alpha machine, version 2. */

/*   make_nn_gathers()  makes the necessary setup for communications
  among the neighboring sites.

  gather() gathers neighbors.
*/

#include <stdio.h>
#include <math.h>
#include "lattice_gn_kr.h"

    /* addresses of neighboring sites */
    /* neighbor[X][i] is a pointer to the neighbor of site lattice[i] in
       direction X */
       
site **neighbor[NDIRS]; 
    
/* FUNCTION MAKE_NN_GATHERS */    
void make_nn_gathers(void){
int i,j,dir;
int xpt,tpt,mul;
register site *s;

void neighbor_coords(int x,int t,int dir,int *xp,int *tp);

    /* Allocate space for lists of pointers to neighbor sites. */
    for(dir=0;dir<NDIRS;dir++){
    neighbor[dir] = (site **)malloc(volume*sizeof(site *)); 
    }
    if(neighbor == NULL) printf("no room for neighbor\n");
    
  /*  FORALLSITES(i,s) {*/
    	for(dir=XUP;dir<NDIRS;dir++) {
    	FORALLSITES(i,s){
    		neighbor_coords(s->x,s->t,dir,&xpt,&tpt);
    		j=site_index(xpt,tpt);
    		neighbor[dir][i] = &(lattice[j]);
    	}
    }	    	

}  /* end of make_nn_gathers */


/* FUNCTION NEIGHBOR_COORDS */
/* utility function for finding coordinates of neighbor */
/*int x,t,dir;	 coordinates of site, and direction (eg XUP) */
/*int *xp,*tp;	 pointers to coordinates of neighbor */
void neighbor_coords(int x,int t,int dir,int *xp,int *tp) 
{
    *xp = x; *tp = t;
    switch(dir){
	case XUP: *xp = (x+1)%nx; break;
	case XDN: *xp = (x+nx-1)%nx; break;
	case TUP: *tp = (t+1)%nt; break;
	case TDN: *tp = (t+nt-1)%nt; break;
	default: printf("BOTCH: bad direction\n"); exit(1);
    }
}


/* FUNCTION GATHER */
/* gather() returns a pointer to the neighbor of given field.
   usage: gather( source, direction, parity, dest )
   example:
	gather( F_OFFSET(phi), XUP, EVEN, gen_pt );
	   gen_pt[i] now contains the address of the phi
	   vector (or a copy therof) on the neighbor of site i in the
	   XUP direction for all even sites i.
	   Do whatever you want with it next.
*/


gather(field,index,parity,dest)
/* arguments */
field_offset field;	/* which field? Some member of structure "site" */
int index;		/* direction to gather from. eg XUP - index into
			   neighbor tables */
int parity;		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
char ** dest;		/* one of the vectors of pointers */
{
/* local variables */
register int i,j;	/* scratch */
register site *s;	/* scratch pointer to site */


    switch(parity){
	case EVEN:
	    FOREVENSITES(j,s){ 
                dest[j] = F_PT(neighbor[index][j],field);
	    }
	    break;
	case ODD:
	    FORODDSITES(j,s){ 
                dest[j] = F_PT(neighbor[index][j],field);
	    }
	    break;
	case EVENANDODD:
	    FORALLSITES(j,s){ 
                dest[j] = F_PT(neighbor[index][j],field);
	    }
	    break;
    }

} /* end of function gather */


/* FUNCTION DCLOCK */
/* Double precision time */
/* This one wraps around after 36 minutes!! It gives the cpu time,
   not the wall clock time */
double dclock(void){
long fine;
    fine = clock();
    return( ((double)fine)/1000000.0 );
} /* end of dclock */


/* FUNCTION TERMINATE */
terminate(status) int status; {
    printf("Termination: status = %d\n",status);
    exit(status);
} /* end of dclock */
