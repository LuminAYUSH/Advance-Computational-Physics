/* **********ran2.c********** */
/* From Numerical Recipes in C, page: 282 */
/* This is a random # generator, to be used to generate
   random configuration on lattice. */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "lattice_gn.h" 

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/* routine RAN2() */
double ran2()
{
int j;
long k;
static long iseed2=123456789;
static long iy=0;
static long iv[NTAB];
double temp;

if(*iseed <= 0){
   if(-(*iseed) < 1) *iseed = 1;
   else{
        *iseed = -(*iseed);
	iseed2 = (*iseed);
	for(j=NTAB+7;j>=0;j--){
	    k = (*iseed)/IQ1;
	    *iseed = IA1 * (*iseed-k*IQ1) -k * IR1;
	    if(*iseed < 0) *iseed += IM1;
	    if(j < NTAB) iv[j] = *iseed;
	   }
	iy = iv[0];
       }
  }
k = (*iseed)/IQ1;
*iseed = IA1*(*iseed-k*IQ1)-k*IR1;
if(*iseed < 0) *iseed += IM1;
k = iseed2/IQ2;
iseed2 = IA2*(iseed2-k*IQ2)-k*IR2;
if(iseed2 < 0) iseed2 += IM2;
j = iy/NDIV;
iy = iv[j] - iseed2;
iv[j] = *iseed;
if(iy < 1) iy += IMM1;
if((temp=AM*iy) > RNMX) return (RNMX);
else return(temp);

}


