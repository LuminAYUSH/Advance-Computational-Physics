/* ************ average_sigma.c ************ */
/* Alpha Mechine, Version 2 */
/* This program calculates the expectation value of sigma-field
   by summing all the sigmas' in the sites and dividing by the 
   volume. This <sigma> is stored in store[] and it returns a
   double <sigma>. */

#include <stdio.h>
#include <math.h>
#include "lattice_gn.h"

double average_sigma(void)
{
int i;
double av_sigma, t_sigma;
site *s;

t_sigma = 0;

FORALLSITES(i,s) t_sigma += s->sigma;

av_sigma = t_sigma / volume;

return(av_sigma);

} /* end of average_sigma.c */

