/******* cg_md.c - CG for staggered fermions ****/
/* version 2    Alpha machine */

/* The source vector is in "src", and the initial guess and answer
   in "dest".  "r" is the residual vector. "p" and "mp" are
   working vectors for the conjugate gradient. 
   MaxCG = maximum number of iterations.
   size_r = desired residual, quit when we reach it.
(Square root def for residue size_r = sqrt(r*r))
*/




#include <stdio.h>
#include <math.h>
#include "lattice_gn_kr.h"

double dtime,dclock();

/* FUNCTION CG_MD */
void cg_md(src,dest,cgiter,residue,cgflag,flavor)
field_offset src,dest; /* src and dest are type psferm */
int cgiter,cgflag,flavor;
double residue;
{
    int N_iter;
    register int i;
    register site *s;
    double size_src;
    double cp,d,dsize_r,dsize_src,size_r;
    register double a, b, c;
    void matp2d(field_offset src,field_offset dest,int isign,int parity,
                int flavor);
    void matd2d(field_offset src,field_offset dest,int isign,int parity);

    /* Start Inversion */

dtime = -dclock();

    /* Normalisation  */
    dsize_src=0.0;
    
    FORALLSITES(i,s) {
        dsize_src += ( ((psferm *)F_PT(s,src))->f[flavor] ) *
                     ( ((psferm *)F_PT(s,src))->f[flavor] );                     
    }
    size_src = sqrt(dsize_src);

    /*printf("beginning inversion--size_src=%e\n",size_src);*/
        
    /* Initial guess */
    
    /* code if you want to start with dest=0... */
    if(cgflag == 0) {
        /*printf("dest_0=0\n");*/
        FORALLSITES(i,s) {
             ((psferm *)F_PT(s,dest))->f[flavor] = 0.0 ;
             s->r = ((psferm *)F_PT(s,src))->f[flavor];
             s->p = s->r;
        }
        dsize_r = 1.00;
        size_r = dsize_r;

    } /* end of if cgflag == 0 */

    /* code if you want to start dest with some particular starting value... */
    /*  r = src - M^dagger M * dest  */
    if(cgflag != 0) {
        /*printf("dest_0  !=0\n");*/
        /* we use mp temporarily to construct r */
        matp2d(dest,F_OFFSET(p),PLUS,EVENANDODD,flavor);
        matd2d(F_OFFSET(p),F_OFFSET(mp),MINUS,EVENANDODD);
        
        dsize_r =0.0;
        FORALLSITES(i,s) {
            s->r = ((psferm *)F_PT(s,src))->f[flavor] - s->mp;
                                  
            dsize_r += (s->r) * (s->r);
            s->p = s->r;                
        }
        size_r = sqrt(dsize_r)/size_src;
        /*printf("beginning inversion--size_r=%e\n",size_r);*/

    }   /* end of if cgflag != 0 */

    /* cp = |r|^2 */
    cp=dsize_r;

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
        /* r = r - a*mmp */
        
        matd2d(F_OFFSET(mp),F_OFFSET(mmp),MINUS,EVENANDODD);
        
        cp = 0.0;
        FORALLSITES(i,s) {
            ((psferm *)F_PT(s,dest))->f[flavor] += a * (s->p);
            s->r -= a * (s->mmp);
            cp += (s->r) * (s->r);        /* cp = | r_(i+1) |^2 */
        }

        b = (cp/c);
        
        /*   p = r + b*p  */
        dsize_r=0.0;
        
        FORALLSITES(i,s) {
            s->p = s->r + b * (s->p);
            dsize_r += (s->r) * (s->r);
        }

        size_r = sqrt(dsize_r)/size_src;
            
    }  /* end of N_iter loop */

/*printf("%d\n", N_iter);*/    
printf("%d\n", N_iter);    

dtime += dclock();
/*if(N_iter != 0) printf("CG_MD: time = %e = %e/site-iter\n",
                        dtime,dtime/(N_iter*volume));*/

    if( (size_r) > residue ) {
        printf(" CG_MD Not Converged\n");
	printf("\n Program is terminated because of CG_MD Convergence. \n");
	terminate(1);
    }
    
}  /* end of cg_md */
