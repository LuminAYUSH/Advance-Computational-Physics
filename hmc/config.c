/* ************ config.c ************ */
/* Version 2, Alpha Machine */
/* This program generates various initial scalar field
   configurations for lattice Gross_Neveu model in 2-d */

#include<stdio.h>
#include<math.h>
#include "lattice_gn.h"


/* routine ZEROLAT */
void zerolat(void)
{
int i;
site *s;

printf(" ZEROLAT: All zero initial config. of `sigma' field.\n");

FORALLSITES(i,s) s->sigma = 0.0;

} /* end of zerolat */

/* routine COLDLAT */
void coldlat(void)
{
int i;
site *s;

printf(" COLDLAT: Cold initial config. of `sigma' field.\n"); 

FORALLSITES(i,s) s->sigma = 1.0;

} /* end of coldlat */

/* routine COLDLAT2 */
void coldlat2(void)
{
int i;
site *s;

printf(" COLDLAT.4: Cold initial config. of `sigma' field.\n"); 

FORALLSITES(i,s) s->sigma = 0.4;

} /* end of coldlat */

/* routine HOTLAT */
void hotlat(void)
{
int i;
site *s;
double ran2(void);

printf(" HOTLAT: Hot initial config. of `sigma' field.\n"); 

FORALLSITES(i,s) s->sigma = (2 * ran2() - 1.0);

} /* end of hotlat */

/* routine FILELAT */
void filelat(FILE *pin)
{
int i;
site *s;

printf(" Configuration to be read from the file sigma.in \n");

FORALLSITES(i,s) fscanf(pin,"%lf", &(s->sigma));

} /* end of filelat */

/* routine FUNNYLAT */ 
void funnylat(void) 
{ 
int i;
site *s;

printf("\n FUNNYLAT: Funny initial config. of `sigma' field.\n");

FORALLSITES(i,s) s->sigma = i;

} /* end of funnylat */
