// I/O utilities for GN sigma field

#include "include_gnkr.h"

///////////////////////////
void save_bin_double(FILE *f_bin)
{

  int i;
  site *s;

  FORALLSITES(i,s) fwrite(&lattice[i].sigma,sizeof(double),1,f_bin);

} // end of save_bin_double()

///////////////////////////
void read_bin_double(FILE *f_bin)
{

  int i,ir;
  site *s;

  FORALLSITES(i,s) ir=fread(&lattice[i].sigma,sizeof(double),1,f_bin);

} // end of read_bin_double()

///////////////////////////
void coldlat(void)
{

  int i;
  site *s;
  FORALLSITES(i,s){ s->sigma=1.0; }

} // end of coldlat()

///////////////////////////
void hotlat(void)
{

  int i;
  site *s;
  //FORALLSITES(i,s){ s->sigma=gsl_ran_gaussian(r,1.0); }
  FORALLSITES(i,s){ s->sigma=gsl_ran_flat(r,-1.0,1.0); }

} // end of hotlat()

