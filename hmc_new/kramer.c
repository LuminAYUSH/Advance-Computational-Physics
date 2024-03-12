// Kramer equation algorithm --
// 1) starts with stochastic initialization of momentum
// 2) next is momentum initial half step
// 3) full single step for sigma and momentum
// 4) momentum final half step
// 5) Hamiltonian calculation before and after MD step
//    followed by Metropolis accept/reject

#include "include_gn.h"

int kramer(void)
{

  int i,nf;
  double term1,term2,xx,z;
  double hold,hnew,deltah;
  site *s;

  // initialize phi=sigma and rho
  FORALLSITES(i,s){ s->phi=s->sigma; s->rho=gsl_ran_gaussian(r,1.0); }

  // momentum initialization
  term1=-1.0*Eps*Step; term2=sqrt(1-exp(2*term1));
  FORALLSITES(i,s) s->pi=exp(term1)*s->mom + term2*s->rho;

  // initial hamiltonian
  hold=hamil(1,1);
  //printf("\nOld Hamiltonian = %e\n\n",hold);

  // momentum first half-step
  piup((double)Step/2.0);

  // leapfrog
  FORALLSITES(i,s) s->phi+=((double)Step*s->pi);

  // momentum final half-step
  piup((double)Step/2.0);

  // final hamiltonian
  hnew=hamil(1,1);
  //printf("\nNew Hamiltonian = %e\n\n",hnew);

  // Metropolis accept/reject
  xx=gsl_ran_flat(r,0.0,1.0);
  deltah=hnew-hold;
  //printf("\nDeltaH = %e\n",deltah);
  z=exp(-1.0*deltah);
  if(xx<z){	// accepted
    //printf("Accepted!\n\n");
    FORALLSITES(i,s){ s->sigma=s->phi; s->mom=s->pi; }
    return(1);
    }
  else{		// rejected
    //printf("Rejected!\n\n");
    FORALLSITES(i,s){ s->mom=-1.0*s->pi; }
    return(0);
    }

} // kramer()
