/* ******************** piup.c ******************** */
/* ***** Alpha Machine, Version 2 ***** */

/* This routine updates the momentum in the leap-frog
   process of the Hybrid Monte Carlo program.       */

#include<stdio.h>
#include<math.h>
#include "lattice_gn.h"


void piup(double t)
{

int i, n;
/*double g_p;*/
site *s;
void matp2d(field_offset src,field_offset dest,int isign,int parity,int flavor);
void cg_md(field_offset src,field_offset dest,int cgiter,double residue,
           int cgflag,int flavor);

/*g_p = (2.01 * 2.01);*/

/*printf("Piup Calcu.\n");*/
printf("Piup Calcu.\n");

/* Momentum update loop */
FORALLSITES(i,s) s->mom -= ((1/(g*g)) * (s->phi) * t);
/*FORALLSITES(i,s) s->mom -= ((1/(g_p)) * (s->phi) * t);*/

/* xi=(M^dagger*M)^(-1)*chi -- xi stored in eta; so eta is dest */
for(n=0;n<nf;n++){
    cg_md(F_OFFSET(chi),F_OFFSET(eta),cgiter2,residue2,1,n);

    matp2d(F_OFFSET(eta),F_OFFSET(p),PLUS,EVENANDODD,n);

    FORALLSITES(i,s)  s->mom += ( 2 * (s->p) * (s->eta.f[n]) * t );	
 
   }

} /* end of piup.c */
