/******* cg_prop.c - CG for staggered fermions to find propagator****/
/* version 1    Alpha machine */

/* The source vector is in "src", and the initial guess and answer
   in "dest".  "r" is the residual vector. "p","mp"and "mmp" are
   working vectors for the conjugate gradient. 
   cgiter = maximum number of iterations.
   residue = desired residual, quit when we reach it.
(Square root def for residue size_r = sqrt(r*r))
THIS PROGRAM SHOULD BE USED WITH src AS F_OFFSET(r)
*/

#include <stdio.h>
#include <math.h>
#include "lattice_gn_kr.h"

double dtime,dclock();

/* FUNCTION CG_PROP */
void cg_prop(src,dest,cgiter,residue,cgflag)
field_offset src,dest; /* src and dest are type double */
int cgiter,cgflag;
float residue;
{
    int N_iter;
    register int i;
    register site *s;
    float size_src;
    double cp,d,dsize_r,dsize_src,size_r;
    register float a, b, c;
    void matd2d(field_offset src,field_offset dest,int isign,int parity);
    

    /* Start Inversion */

dtime = -dclock();

    /* Normalisation  */
    dsize_src=0.0;
    
    FORALLSITES(i,s) {
        dsize_src += ( *((double *)F_PT(s,src)) ) *
                     ( *((double *)F_PT(s,src)) );                     
    }
    size_src = (float)sqrt(dsize_src);

    /*printf("beginning inversion--size_src=%e\n",size_src);*/
        
    /* Initial guess */
    
    /* code if you want to start with dest=0... */
    if(cgflag == 0) {
        /*printf("dest_0=0\n");*/
        FORALLSITES(i,s) {
             *((double *)F_PT(s,dest)) = 0.0 ;
             s->r = *((double *)F_PT(s,src));
        }
        dsize_r = 1.00;
        size_r = dsize_r;
        /**printf("size_r=%e\n",size_r);**/

    }
    /* code if you want to start dest with some particular starting value... */
    /*  r = src - M * dest  */
    if(cgflag != 0) {
        /*printf("dest_0  !=0\n");*/
        /* construct r , use mp as temp location */
        matd2d(dest,F_OFFSET(mp),PLUS,EVENANDODD);
        
        dsize_r =0.0;
        FORALLSITES(i,s) {
            s->r = *((double *)F_PT(s,src)) - s->mp;
                                  
            dsize_r += (s->r) * (s->r);                
        }
        size_r = (float)sqrt(dsize_r)/size_src;
        /*printf("beginning inversion--size_r=%e\n",size_r);*/

    }   /* end of if cgflag != 0 */
    
    /* construct p */
    matd2d(F_OFFSET(r),F_OFFSET(p),MINUS,EVENANDODD);
    

    /* cp = |p|^2 */
    cp=0.0;
    FORALLSITES(i,s) {
    cp += (s->p) * (s->p);
    }

/* start of CG iteration loop  */
    
    for( N_iter = 0; N_iter < cgiter && size_r > residue; N_iter++) { 

        c=cp;
        
        /*  mp = M * p */
        matd2d(F_OFFSET(p),F_OFFSET(mp),PLUS,EVENANDODD);

        /* d = |mp|^2  */
        d=0.0;

        FORALLSITES(i,s) {
            d += (s->mp) * (s->mp);
        }
    
        a = (c/d);

        /* dest = dest + a*p  */
        /* r = r - a*mp */
        
        FORALLSITES(i,s) {
            *((double *)F_PT(s,dest)) += a * (s->p);
            s->r -= a * (s->mp);
        }
        
        /* construct M^dagger r */
        matd2d(F_OFFSET(r),F_OFFSET(mp),MINUS,EVENANDODD);
        
        cp=0.0;
        FORALLSITES(i,s) {
             cp += (s->mp) * (s->mp);
        }     
        

        b = (cp/c);
        
        /*   p = M^dagger r + b*p  */
        dsize_r=0.0;
        
        FORALLSITES(i,s) {
            s->p = s->mp + b * (s->p);
            dsize_r += (s->r) * (s->r);
        }

        size_r = (float)sqrt(dsize_r)/size_src;
        /* printf("iteration= %d, residue= %e\n",N_iter,size_r); */
            
    }  /* end of N_iter loop */
    
dtime += dclock();
/*printf("CG_PROP: time = %e = %e/site-iter\n",
dtime,dtime/(N_iter*volume));*/

    if( (size_r) > residue ) {
        printf(" CG_PROP Not Converged\n");
    }
    
}  /* end of cg_prop */
