/* **********randomize.c********** */
/* ***** Alpha Machine, version 2 ***** */

long *iseed;

#include <stdio.h>
#include <math.h>
#include "lattice_gn_kr.h"

/* function to initialize the random # generator using the
   system clock.    */

void randomize(void)
{
iseed = (long *)malloc(sizeof(long));
*iseed = -849;

printf("\n SEED = %ld\n", *iseed);

} /* end of randomize.c */
 
