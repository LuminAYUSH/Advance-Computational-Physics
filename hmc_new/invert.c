// Conjugate Gradient inversion : pg: 175 -- 177 of Carleton's book
// Source "src" and destination "dest" are of type "psferm"

#include "include_gn.h"

void congrad(field_offset src,field_offset dest,int cgiter,
	double residue,int cgflag,int flv)
{

  int i,niter;
  double a,b,c,cp,d;
  double size_src,dsize_src,size_r,dsize_r;
  site *s;

  // normalisation
  dsize_src=0.0;
  FORALLSITES(i,s){
     dsize_src+=( ((psferm *)F_PT(s,src))->f[flv] )*
                ( ((psferm *)F_PT(s,src))->f[flv] );                     
     }
  size_src=sqrt(dsize_src);
  // printf("beginning inversion--size_src=%e\n",size_src);

  // Initial guess r0=p0=src - M^dag M*guess_dest
  if(cgflag==0){
    FORALLSITES(i,s){
       ((psferm *)F_PT(s,dest))->f[flv]=0.0 ;
       s->r=((psferm *)F_PT(s,src))->f[flv];
       s->p=s->r;
       }
    dsize_r=1.00;
    size_r=dsize_r;
    } 
  else{
    matp2d(dest,F_OFFSET(mp),PLUS,flv);
    matd2d(F_OFFSET(mp),F_OFFSET(mmp),MINUS);
    dsize_r=0.0;
    FORALLSITES(i,s){
       s->r=((psferm *)F_PT(s,src))->f[flv]-s->mmp;
       dsize_r+=(s->r)*(s->r);
       s->p=s->r;                
       }
    size_r=sqrt(dsize_r)/size_src;
    }
  // printf("beginning inversion--size_r = %e\n",size_r);

  // cp=|r0|^2=|p0|^2
  cp=dsize_r;

  // start of CG iteration loop  */
  for(niter=0;niter<cgiter && size_r>residue;niter++){ 
     c=cp;
        
     // mmp=M^dag M*p
     matd2d(F_OFFSET(p),F_OFFSET(mp),PLUS);
     matd2d(F_OFFSET(mp),F_OFFSET(mmp),MINUS);

     //d=|mmp|^2=(p,M^dag M p)
     d=0.0;
     FORALLSITES(i,s) { d+=(s->p)*(s->mmp); }
     a=(c/d);

     // dest=dest+a*p, r=r-a*mmp
     cp=0.0;
     FORALLSITES(i,s){
        ((psferm *)F_PT(s,dest))->f[flv]+=a*(s->p);
        s->r-=a*(s->mmp);
        cp+=(s->r)*(s->r);        // cp=|r_(i+1)|^2
        }
     b=(cp/c);
        
     // p=r+b*p
     dsize_r=0.0;   
     FORALLSITES(i,s){
        s->p=s->r+b*(s->p);
        dsize_r+=(s->r)*(s->r);
        }
     size_r=sqrt(dsize_r)/size_src;
            
    }  // niter-loop

/*
  if((size_r)>residue){
    printf(" congrad() Not Converged!!\n");
    }
  else{
    printf("congrad() converged! Iteration=%d, Residue=%e\n",niter,size_r);
    }
*/

}  // end of congrad()

