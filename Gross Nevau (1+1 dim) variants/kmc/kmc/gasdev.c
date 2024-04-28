/* **********gasdev.c********** */
/* ***** Alpha Machine, Version 2 ***** */
/* This routine produces random # with Gaussian distribution --
   taken from Numerical Recipes in C, page:289  , to be used in
   generating configuration and/or Gaussian noise in lattice. */

#include <stdio.h>
#include <math.h>
#include "lattice_gn_kr.h"


double gasdev(void)
{
static int iset = 0;
static double gset;
double fac, rsq, v1, v2;
double ran2(void);

if(iset == 0){
   do{
      v1 = (2.0 * ran2()) - 1.0;
      v2 = (2.0 * ran2()) - 1.0;
      rsq = (v1 * v1) + (v2 * v2);
     }while(rsq >= 1.0 || rsq == 0.0);
   fac = sqrt(-2.0 * log(rsq)/rsq);
   gset = v1 * fac;
   iset = 1;
   return(v2*fac);
  }
else{
     iset = 0;
     return(gset);
    }

} /* end of gasdev.c */
 
