/* ******************** matd2d.c ******************** */
/* Alpha machine, Version 2 */
/* 2d GN model with staggered fermions */

#include<stdio.h>
#include<math.h>
#include "lattice_gn.h"


/* FUNCTION MATD2D */

void matd2d(src,dest,isign,parity) 
field_offset src,dest; int isign,parity;
 
{

/* This program multiplies the fermion matrix or its adjoint with a
   a pseudovector. It accepts a "double" and returns a "double". */

int i, n;
site *s;
void gather(field_offset field,int index,int parity,char **dest);

gather(src, XUP, EVENANDODD, gen_pt);
FORALLSITES(i,s)
   { 
   *((double *)F_PT(s,dest)) = s->sign * 0.5 * (*((double *)gen_pt[i]));
   }

gather(src, XDN, EVENANDODD, gen_pt);
FORALLSITES(i,s)
   {
   *((double *)F_PT(s,dest)) -= s->sign * 0.5 * (*((double *)gen_pt[i]));
   }

gather(src, TUP, EVENANDODD, gen_pt);
FORALLSITES(i,s)
   {
   *((double *)F_PT(s,dest)) += 0.5 * (*((double *)gen_pt[i]));
   }

gather(src, TDN, EVENANDODD, gen_pt);
FORALLSITES(i,s)
   { 
   *((double *)F_PT(s,dest)) -= 0.5 * (*((double *)gen_pt[i]));
   }

FORALLSITES(i,s)
   {
   *((double *)F_PT(s,dest)) = (isign * (*((double *)F_PT(s,dest)))) + 
                               ((s->phi) * (*((double *)F_PT(s,src))));
   }

} /* end of routine matd2d.c */

