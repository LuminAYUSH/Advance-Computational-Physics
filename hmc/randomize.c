/* **********randomize.c********** */
/* ***** Alpha Machine, version 2 ***** */

long *iseed;

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

/* function to initialize the random # generator using the
   system clock.    */

void randomize(void)
{
iseed = (long *)malloc(sizeof(long));
*iseed = -772;

printf("\n SEED = %ld\n", *iseed);

} /* end of randomize.c */
 
