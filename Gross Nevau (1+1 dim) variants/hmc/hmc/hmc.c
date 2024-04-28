/* ******************** hmc.c ******************** */
/* ***** Alpha Machine, version 2 ***** */

/* This routine performs the Hybrid Monte Carlo calculations
   using the routines gasdev.c, gauss.c, piup.c and hamil.c.      */

#include<stdio.h>
#include<math.h>
#include "lattice_gn.h"


int hmc(void) 
{

int i, j, m, n, max_count=100;	/* max_count to be read from outside, later */
int count;
double hold, hnew, xx=1.0, deltah, z=0;
site *s;
double ran2(void);
double gasdev(void);
double gauss(void);
double hamil(int hflag,int flag);
void piup(double t);
void matp2p(field_offset src,field_offset dest,int isign,int parity,int flavor);

++counter;

for(count=0; count<max_count && xx>z; count++){

   FORALLSITES(i,s){
      s->phi = s->sigma; 
      s->mom = gasdev();		/* momentum assignment */
      for(n=0; n<nf; n++){
          s->eta.f[n] = gauss();	/* eta assignment */
	 }
     }

   /* chi = M^dagger *eta */
   for(n=0;n<nf;n++){
       matp2p(F_OFFSET(eta),F_OFFSET(chi),MINUS,EVENANDODD,n);
      }
  
   /* initial value of Hamiltonian */
   hold = hamil(0,1); 

   /* initial half-step */
   piup(step/2);
  
   /* leap-frog for n-1 steps */
   for(j=0; j<mdstep-1; j++){
       FORALLSITES(i,s){
          s->phi += (s->mom * step);
         }
       piup(step);
      }

    /* final n-th step */
    FORALLSITES(i,s){ 
       s->phi += (s->mom * step); 
      }

    /* final half-step */
    piup(step/2);

    /* calculation of new Hamiltonian */ 
    hnew = hamil(1,1); 

    /* Accept/Reject */
    xx = ran2();
    deltah = hnew - hold;

    if(deltah>DELTAMAX){ 
       printf("HMC loop #%d REJECTED in CALL #%d for LARGE DELTAH.\n", 
               count+1, counter);
       printf("\n Program is terminated. \n");
       FORALLSITES(i,s) printf("%lf\n", con[i]);
       printf("\n AVERAGE SIGMAS \n");
       for(m=0;m<meas_loop;m++) printf("%lf\n", ac_store[m]);

       terminate(1);
             } 
	else {
    	      z = exp(- deltah);
             }

   } /* for loop ends here */

/* Assign the phis' to sigma in the struct->lattice */
if(xx <= z){
   	   FORALLSITES(i,s) s->sigma = s->phi;
           }

return(count); 

} /* end of hmc.c */

