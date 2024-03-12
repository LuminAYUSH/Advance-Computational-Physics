// Kramer equation algorithm --
// 1) starts with stochastic initialization of momentum
// 2) next is momentum initial half step
// 3) full single step for sigma and momentum
// 4) momentum final half step
// 5) Hamiltonian calculation before and after MD step
//    followed by Metropolis accept/reject

#include "include_gn.h"

int hmc(void)
{

  int i,j,nf,ctr;
  double term1,term2,xx,z;
  double hold,hnew,deltah;
  site *s;

  xx=1.0;z=0.0;
  for(ctr=0;ctr<(max_hmc_iter && xx>z);ctr++){	// loop till we get
						// an acceptance or max_iter

     // initialize phi=sigma and momentum=pi
     FORALLSITES(i,s){
	s->phi=s->sigma;
	s->pi=gsl_ran_gaussian(r,1.0);
	for(nf=0;nf<Nf;nf++) s->zeta.f[nf]=gsl_ran_gaussian(r,1.0/sqrt(2.0));
	}

     // initial hamiltonian
     hold=hamil(1,1);
     //printf("\nOld Hamiltonian = %e\n\n",hold);

     // momentum first half-step
     piup((double)Step/2.0);

     // leapfrog for n-1 steps
     for(j=0;j<(Mdstep-1);j++){
	FORALLSITES(i,s){ s->phi+=(s->pi*(double)Step);
	piup((double)Step);
	}

     // final n-th step
     FORALLSITES(i,s){ s->phi+=(s->pi*(double)Step); }

     // momentum final half-step
     piup((double)Step/2.0);

     // final hamiltonian
     hnew=hamil(1,1);
     //printf("\nNew Hamiltonian = %e\n\n",hnew);

     // Metropolis accept/reject
     xx=gsl_ran_flat(r,0.0,1.0);
     deltah=hnew-hold;
     z=exp(-1.0*deltah);
     //printf("\nDeltaH = %e\n",deltah);

    // keep looping till xx <= z

     } // ctr-loop ends

  if(xx<z){     // accepted
    //printf("Accepted!\n\n");
    FORALLSITES(i,s){ s->sigma=s->phi; }
    }

  return(ctr);

} // hmc()
