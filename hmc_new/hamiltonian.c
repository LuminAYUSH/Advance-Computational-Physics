// Gross-Neveu Hamiltonian for KS fermion
// (eqn. 4.8, mythesis.pdf)
// H = sum_x [ sigma^2/2g^2 + p^2/2 + sum_f zeta_f^d (A^dA)^{-1} zeta_f ]
// hflag == skip inversion, flag = A^{-1} or (A^d A)^{-1}

#include "include_gn.h"

double hamil(int hflag,int flag)
{

  int i,nf;
  double hml=0,psf_hml;
  site *s;

  // sigma^2/2g^2 + p^2/2, phi=sigma
  FORALLSITES(i,s) hml+=( (0.5/(G*G))*(s->phi)*(s->phi) +
			  (0.5)*(s->pi)*(s->pi) );

  if(hflag!=0){
  
    // sum_f zeta_f^d (A^dA)^{-1} zeta_f
    for(nf=0;nf<Nf;nf++){

       // eta_f=(A^dagger A)^{-1} zeta_f
       congrad(F_OFFSET(zeta),F_OFFSET(eta),CgiterH,ResidueH,flag,nf);
       FORALLSITES(i,s) hml+=(s->zeta.f[nf])*(s->eta.f[nf]);
       }
    } // if-ends
  else{
    for(nf=0;nf<Nf;nf++){
       FORALLSITES(i,s) hml+=(s->zeta.f[nf])*(s->eta.f[nf]);
       }
    }

  return(hml);

} // end of hamil()

