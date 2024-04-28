/* ******************** matd2p.c ******************** */
/* Alpha machine, Version 2 */
/* 2d GN model with staggered fermions */

#include<stdio.h>
#include<math.h>
#include "lattice_gn_kr.h"


/* FUNCTION MATD2P */

void matd2p(src,dest,isign,parity,flavor) 
field_offset src,dest; int isign,parity,flavor;
 
{

/* This program multiplies the fermion matrix or its adjoint with a
   a pseudovector. It accepts a "double" and returns a "psferm". */

int i, n;
site *s;
void gather(field_offset field,int index,int parity,char **dest);

gather(src, XUP, EVENANDODD, gen_pt);
FORALLSITES(i,s)
   { 
   ((psferm *)F_PT(s,dest))->f[flavor] = s->sign * 0.5 * 
                                         (*((double *)gen_pt[i]));
   }
gather(src, XDN, EVENANDODD, gen_pt);
FORALLSITES(i,s)
   {
   ((psferm *)F_PT(s,dest))->f[flavor] -= s->sign * 0.5 * 
                                          (*((double *)gen_pt[i]));
   }
gather(src, TUP, EVENANDODD, gen_pt);
FORALLSITES(i,s)
   {
   ((psferm *)F_PT(s,dest))->f[flavor] += 0.5 * 
                                          (*((double *)gen_pt[i]));
   }
gather(src, TDN, EVENANDODD, gen_pt);
FORALLSITES(i,s)
   { 
   ((psferm *)F_PT(s,dest))->f[flavor] -= 0.5 * 
                                          (*((double *)gen_pt[i]));
   }
FORALLSITES(i,s)
   {
   ((psferm *)F_PT(s,dest))->f[flavor] = 
                           (isign * ((psferm *)F_PT(s,dest))->f[flavor]) +
                           ((s->phi) * (*((double *)F_PT(s,src))));
   }

} /* end of routine matd2p.c */

