// Momentum equation of motion -- half step in momentum
// involves calculation of force (eqn 4.10, mythesis.pdf)

#include "include_gnkr.h"

double piup(double tau)
{

  int i,nf;
  site *s;

  // initial momentum is already in s->pi
  // first add (-sigma/g^2)*tau
  FORALLSITES(i,s) s->pi-=(1/(G*G))*(s->phi)*tau;

  // now, sum_f [ omega^d_f xi_f + h.c. ]
  for(nf=0;nf<Nf;nf++){
     congrad(F_OFFSET(zeta),F_OFFSET(xi),CgiterP,ResidueP,1,nf);
     matp2p(F_OFFSET(xi),F_OFFSET(omega),PLUS,nf);
     FORALLSITES(i,s) s->pi+=(2*s->omega.f[nf]*s->xi.f[nf])*tau;
     } // nf-loop

} // piup()

