/* ********** autocorel.c ********** */
/* *** Alpha Machine, version 2 ***  */
/* This routine calculates the unnormalized autocorrelated functions and
   normalized autocorrelated function along with integrated autocorrela-
   tion time. It reads some T_cut and calculates var(int. ac time).  */

#include <stdio.h>
#include <math.h>
#include "lattice_gn.h"

void autocorel(double sigma_av, int lb, int a_index)
{

int j, k, t, N_t, tcut, u;
double ct, c0, rho[MAXT_cut];

c0=0.0; N_t=0; tcut=T_cut;

/* Calculating the unnormalized autocorrelation functions ct and c0
   and normalized autocorrelation function rho. */

for(j=lb;j<(lb+seg_length);j++){
    c0 += ((ac_store[j] - sigma_av) * (ac_store[j] - sigma_av))/seg_length;
    }	

for(u=0;u<NOT_cut;u++){

    for(t=1;t<=tcut;t++){

        N_t = seg_length - t;
        ct = 0.0;

        for(j=lb;j<(lb+N_t);j++){
	    ct += ((ac_store[j] - sigma_av) * (ac_store[j+t] - sigma_av))/N_t;
	                        }

        rho[t-1] = ct/c0;

        } /* t-loop ends */


    for(t=0;t<tcut;t++) T_int[u][t][a_index] = 0.5;

    /* Calculation of tau_int, the integrated autocorrelation time */
    for(t=0;t<tcut;t++){
    
        for(k=0;k<=t;k++) T_int[u][t][a_index] += rho[k];

                       }
         
    tcut += D_cut;
    
   } /* end of variable T_cut loop */      

} /* end of autocorel.c */

