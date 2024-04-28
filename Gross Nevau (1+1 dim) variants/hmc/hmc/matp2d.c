/* ******************** matp2d.c ******************** */
/* Alpha machine, Version 2 */
/* 2d GN model with staggered fermions */

#include<stdio.h>
#include<math.h>
#include "lattice_gn.h"


/* FUNCTION MATP2D */

void matp2d(src,dest,isign,parity,flavor) 
field_offset src,dest; int isign,parity,flavor;
 
{

/* This program multiplies the fermion matrix or its adjoint with a
   a pseudovector. It accepts a "psferm" and returns a "double". */

int i, n;
site *s;
void gather(field_offset field,int index,int parity,char **dest);

gather(src, XUP, EVENANDODD, gen_pt);
FORALLSITES(i,s)
   { 
   *((double *)F_PT(s,dest)) = s->sign * 0.5 * 
                               ((psferm *)gen_pt[i])->f[flavor];
   }
gather(src, XDN, EVENANDODD, gen_pt);
FORALLSITES(i,s)
   {
   *((double *)F_PT(s,dest)) -= s->sign * 0.5 * 
                                ((psferm *)gen_pt[i])->f[flavor];
   }
gather(src, TUP, EVENANDODD, gen_pt);
FORALLSITES(i,s)
   {
   *((double *)F_PT(s,dest)) += 0.5 * ((psferm *)gen_pt[i])->f[flavor];
   }
gather(src, TDN, EVENANDODD, gen_pt);
FORALLSITES(i,s)
   { 
   *((double *)F_PT(s,dest)) -= 0.5 * ((psferm *)gen_pt[i])->f[flavor];
   }
FORALLSITES(i,s)
   {
   *((double *)F_PT(s,dest)) = (isign * (*((double *)F_PT(s,dest)))) + 
                              ((s->phi) * ((psferm *)F_PT(s,src))->f[flavor]);
   }

} /* end of routine matp2d.c */

