/* ******************** matp2p.c ******************** */
/* Alpha machine, Version 2 */
/* 2d GN model with staggered fermions */

#include<stdio.h>
#include<math.h>
#include "lattice_gn.h"

/* FUNCTION MATP2P */

void matp2p(src,dest,isign,parity,flavor) 
field_offset src,dest; int isign,parity,flavor;
 
{

/* This program multiplies the fermion matrix or its adjoint with a
   a pseudovector. It accepts a "psferm" and returns a "psferm". */

int i, n;
site *s;
void gather(field_offset field,int index,int parity,char **dest);


gather(src,XUP,EVENANDODD,gen_pt);
FORALLSITES(i,s)
   { 
   ((psferm *)F_PT(s,dest))->f[flavor] = s->sign * 0.5 * 
                                       ((psferm *)gen_pt[i])->f[flavor];
   }
gather(src,XDN,EVENANDODD,gen_pt);
FORALLSITES(i,s)
   {
   ((psferm *)F_PT(s,dest))->f[flavor] -= s->sign * 0.5 * 
                                       ((psferm *)gen_pt[i])->f[flavor];
   }
gather(src,TUP,EVENANDODD,gen_pt);
FORALLSITES(i,s)
   {
   ((psferm *)F_PT(s,dest))->f[flavor] += 0.5 * 
                                       ((psferm *)gen_pt[i])->f[flavor];
   }
gather(src,TDN,EVENANDODD,gen_pt);
FORALLSITES(i,s)
   { 
   ((psferm *)F_PT(s,dest))->f[flavor] -= 0.5 * 
                                       ((psferm *)gen_pt[i])->f[flavor];
   }
FORALLSITES(i,s)
   {
   ((psferm *)F_PT(s,dest))->f[flavor] = 
                           (isign * ((psferm *)F_PT(s,dest))->f[flavor]) +
                           ((s->phi) * ((psferm *)F_PT(s,src))->f[flavor]);
               
   }

} /* end of routine matp2p.c */

