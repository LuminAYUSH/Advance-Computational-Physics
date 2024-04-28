/* ******************** kramer.c ******************** */

/* This routine performs the Kramers' algorithm
   using the routines gasdev.c, piup.c and hamil.c.      */

#include<stdio.h>
#include<math.h>
#include "lattice_gn_kr.h"


int kramer(void) 
{

int i, j, m, n;
double hold, hnew, term1, term2; 
double xx=1.0, deltah, z=0;
site *s;

double ran2(void);
double gasdev(void);
double hamil(int hflag,int flag);
void piup(double t);

++counter;


   FORALLSITES(i,s){
      		   s->phi = s->sigma; 
      		   s->zeta = gasdev();
     		   }

   /* momentum calculation*/
   term1 = gama * step;
   term2 = sqrt( 1 - exp(-(2*term1)) ); 
   FORALLSITES(i,s) s->mom = ( ((s->mom1)*(exp(-term1))) + (term2*(s->zeta)) );

   /* initial value of Hamiltonian */
   hold = hamil(1,1);

   /* Leapfrog integration */
   piup(step/2); 	
  
   FORALLSITES(i,s) s->phi += (s->mom * step);

   piup(step/2); 
   /* Leapfrog ends */

    /* calculation of new Hamiltonian */ 
    hnew = hamil(1,1); 

    /* Accept/Reject */
    xx = ran2();
    deltah = hnew - hold;
    z = exp(- deltah);

    /*printf("xx= %lf \t deltah= %lf \n", xx, deltah);*/

    if(xx>z){
	    if(deltah>DELTAMAX){
		printf("KRAMER loop REJECTED in CALL #%d for LARGE DELTAH.\n",
     		    counter);
		FORALLSITES(i,s) printf("%lf\n", con[i]);
		printf("\n AVERAGE SIGMAS \n");
		for(m=0;m<meas_loop;m++) printf("%lf\n", ac_store[m]);
		terminate(1);
			   } 
	    /*printf("KRAMER loop REJECTED IN CALL #%d.\n", counter);*/
	    FORALLSITES(i,s) s->mom1 = -(s->mom);
	    return(0);
	    }
    else    {
	    /*printf("\n Accepted in call #%d\n\n", counter);*/
	    FORALLSITES(i,s) s->sigma = s->phi;
	    FORALLSITES(i,s) s->mom1 = s->mom;
	    return(1);
            }

} /* end of kramer.c */

