/* ************ propagator.c ************ */

/* Main procedure for calculation of propagator and <psi_bar
   psi>, to be called thro' gn_hmc.c */

#include <stdio.h>
#include <math.h>
#include "lattice_gn.h"    /* global variables for lattice fields */

double propagator(void)
{
int prompt, i, m, n, n0, t, x, xl, xu, source;
double lpsi;
site *s;

void cg_prop(field_offset src,field_offset dest,int cgiter,float residue,
             int cgflag);

lpsi = 0.0;
for(t=0;t<nt;t++) tprop[t] = 0.0;
FORALLSITES(i,s) s->r = 0.0;
 
for(t=0;t<nt;t++){

    source = (t * nx) + mid;
    lattice[source].r = 1.0;

    cg_prop(F_OFFSET(r),F_OFFSET(mmp),cgiter1,residue1,0);

    for(x=0;x<nx;x++){
        xl = x * nx;
    	xu = (((x+1) * nx) -1);
 	n0 = x - t;
	if(n0<0) n0 += nt;
	prop[n0] = 0.0;

        for(m=xl;m<xu;m++) prop[n0] += lattice[m].mmp;

	tprop[n0] += prop[n0];
        }

    lpsi += lattice[source].mmp;
    }

for(t=0;t<nt;t++){
    tprop[t] = (tprop[t]/nt); 
    /* G_prop[t] += tprop[t]; */ 
    }

lpsi = lpsi/(double)nt ;

return(lpsi);

} /* end of propagator.c */
