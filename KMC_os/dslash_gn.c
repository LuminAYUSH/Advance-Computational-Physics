// Dirac matrix multiplication for Gross-Neveu model
// Multiplication involves 4 different scenario :
// M * pseudofermion = pseudofermion	=> matp2p()
// M * pseudofermion = scalar (double)	=> matp2d()
// M * scalar (double) = pseudofermion	=> matd2p()
// M * scalar (double) = scalar (double)=> matd2d()

#include "include_gnkr.h"

////////////////////////////////////////////////
void matp2p(field_offset src,field_offset dest,int isign,int flav)
{

  int i;
  site *s;

  gather(src,Xup,gen_pt);
  FORALLSITES(i,s){
     ((psferm *)F_PT(s,dest))->f[flav]=(s->sign)*0.5*
				       ((psferm *)gen_pt[i])->f[flav];
     }
  gather(src,Xdn,gen_pt);
  FORALLSITES(i,s){
     ((psferm *)F_PT(s,dest))->f[flav]-=(s->sign)*0.5*
					((psferm *)gen_pt[i])->f[flav];
     }

  gather(src,Tup,gen_pt);
  FORALLSITES(i,s){
     ((psferm *)F_PT(s,dest))->f[flav]+=0.5*
				        ((psferm *)gen_pt[i])->f[flav];
     }
  gather(src,Tdn,gen_pt);
  FORALLSITES(i,s){
     ((psferm *)F_PT(s,dest))->f[flav]-=0.5*
					((psferm *)gen_pt[i])->f[flav];
     }

  FORALLSITES(i,s){
     ((psferm *)F_PT(s,dest))->f[flav]=(isign*
				((psferm *)F_PT(s,dest))->f[flav]) + 
                                ((s->phi)*((psferm *)F_PT(s,src))->f[flav]);
     }

} // matp2p()

////////////////////////////////////////////////
void matp2d(field_offset src,field_offset dest,int isign,int flav)
{

  int i;
  site *s;

  gather(src,Xup,gen_pt);
  FORALLSITES(i,s){
     *((double *)F_PT(s,dest))=(s->sign)*0.5*
				       ((psferm *)gen_pt[i])->f[flav];
     }
  gather(src,Xdn,gen_pt);
  FORALLSITES(i,s){
     *((double *)F_PT(s,dest))-=(s->sign)*0.5*
					((psferm *)gen_pt[i])->f[flav];
     }

  gather(src,Tup,gen_pt);
  FORALLSITES(i,s){
     *((double *)F_PT(s,dest))+=0.5*((psferm *)gen_pt[i])->f[flav];
     }
  gather(src,Tdn,gen_pt);
  FORALLSITES(i,s){
     *((double *)F_PT(s,dest))-=0.5*((psferm *)gen_pt[i])->f[flav];
     }

  FORALLSITES(i,s){
     *((double *)F_PT(s,dest))=(isign*(*((double *)F_PT(s,dest)))) + 
                               ((s->phi)*((psferm *)F_PT(s,src))->f[flav]);
     }

} // matp2d()

////////////////////////////////////////////////
void matd2p(field_offset src,field_offset dest,int isign,int flav)
{

  int i;
  site *s;

  gather(src,Xup,gen_pt);
  FORALLSITES(i,s){
     ((psferm *)F_PT(s,dest))->f[flav]=(s->sign)*0.5*
				       (*((double *)gen_pt[i]));
     }
  gather(src,Xdn,gen_pt);
  FORALLSITES(i,s){
     ((psferm *)F_PT(s,dest))->f[flav]-=(s->sign)*0.5*
					(*((double *)gen_pt[i]));
     }

  gather(src,Tup,gen_pt);
  FORALLSITES(i,s){
     ((psferm *)F_PT(s,dest))->f[flav]+=0.5*(*((double *)gen_pt[i]));
     }
  gather(src,Tdn,gen_pt);
  FORALLSITES(i,s){
     ((psferm *)F_PT(s,dest))->f[flav]-=0.5*(*((double *)gen_pt[i]));
     }

  FORALLSITES(i,s){
     ((psferm *)F_PT(s,dest))->f[flav]=(isign*
				((psferm *)F_PT(s,dest))->f[flav]) + 
                                ((s->phi)*(*((double *)F_PT(s,src))));
     }

} // matd2p()

////////////////////////////////////////////////
void matd2d(field_offset src,field_offset dest,int isign)
{

  int i;
  site *s;

  gather(src,Xup,gen_pt);
  FORALLSITES(i,s){
     *((double *)F_PT(s,dest))=(s->sign)*0.5*(*((double *)gen_pt[i]));
     }
  gather(src,Xdn,gen_pt);
  FORALLSITES(i,s){
     *((double *)F_PT(s,dest))-=(s->sign)*0.5*(*((double *)gen_pt[i]));
     }

  gather(src,Tup,gen_pt);
  FORALLSITES(i,s){
     *((double *)F_PT(s,dest))+=0.5*(*((double *)gen_pt[i]));
     }
  gather(src,Tdn,gen_pt);
  FORALLSITES(i,s){
     *((double *)F_PT(s,dest))-=0.5*(*((double *)gen_pt[i]));
     }

  FORALLSITES(i,s){
     *((double *)F_PT(s,dest))=(isign*(*((double *)F_PT(s,dest)))) + 
                               ((s->phi)*(*((double *)F_PT(s,src))));
     }

} // matd2d()

