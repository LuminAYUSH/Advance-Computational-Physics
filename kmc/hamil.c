/* ******************** hamil.c ******************** */
/* ***** Alpha Machine, version 2 ***** */

/* This routine calculates the hamiltonian, required by the
   Accept/reject test, and it accepts a matrix that contains 
   the old or updated sigma fields. It returns a double 
   precision value for hamiltonian.               */

#include<stdio.h>
#include<math.h>
#include "lattice_gn_kr.h"


double hamil(int hflag,int flag)
{

int i, n;
double h=0;
site *s;
void cg_md(field_offset src,field_offset dest,int cgiter,double residue,
           int cgflag,int flavor);

/*printf("\nHamiltonian Calculation\n");*/

FORALLSITES(i,s){   
	h += (((0.5/(g*g)) * (s->phi)* (s->phi)) + (0.5 * (s->mom) * (s->mom)));
     		}

if(hflag!=0){

     printf("\nHamiltonian Calculation\n");
     /* the inversion: eta = (M^dagger * M)^(-1) * chi; eta is sent as dest */
     for(n=0;n<nf;n++){
         cg_md(F_OFFSET(chi),F_OFFSET(eta),cgiter1,residue1,flag,n);

         FORALLSITES(i,s) h += ((s->chi.f[n]) * (s->eta.f[n]));
        } 
     } 
else {
     for(n=0;n<nf;n++){
 	 FORALLSITES(i,s) h += s->eta.f[n] * s->eta.f[n];
        }
     }

return(h);

} /* end of hamil() */ 

