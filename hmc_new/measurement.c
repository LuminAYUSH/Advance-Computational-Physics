// Average sigma value
// Average psi-bar -- psi

#include "include_gn.h"

double average_sigma(void)
{

  int i;
  double av_sigma;
  site *s;

  av_sigma=0;
  FORALLSITES(i,s) av_sigma+=s->sigma;

  return(av_sigma/Volume);

} // end of average_sigma()

double average_pbp(void)
{

  int i,j,t;
  double av_pbp,tprop[Nt];
  site *s;

  av_pbp=0.0;
  for(t=0;t<Nt;t++) tprop[t]=0.0;

  return(av_pbp/(double)Nt);

} // end of av_pbp()

