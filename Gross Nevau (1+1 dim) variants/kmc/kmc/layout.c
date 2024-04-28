/*****************layout.c**************************************/
/* ***** Alpha Machine, Version 2 ***** */

#include <stdio.h>
#include <math.h>
#include "lattice_gn_kr.h"

/* FUNCTION LAYOUT */
void layout(void)   {
#ifdef EVENFIRST
	printf("EVENFIRST,");

#endif
	printf("\n");
#ifdef EVENFIRST
    /* Need even number of lattice sites */
    if( volume%2 != 0){
	printf("SORRY, CAN'T LAY OUT THIS LATTICE WITH EVENFIRST\n");
	/*terminate(0);*/
	exit(1);
    }
#endif
    no_even_sites = no_odd_sites = volume/2;
    
}   /* end of layout */


/* FUNCTION SITE_INDEX */
int site_index(x,t) int x,t; {
register int i,xr,tr;
    xr = x%nx;
    tr = t%nt;
    i = xr + nx * tr;
#ifdef EVENFIRST
    if( (x+t)%2==0 ){	/* even site */
	return( i/2 );
    }
    else {
	return( (i + volume)/2 );
    }
#else
    return( i );
#endif

}   /* end of site_index */

