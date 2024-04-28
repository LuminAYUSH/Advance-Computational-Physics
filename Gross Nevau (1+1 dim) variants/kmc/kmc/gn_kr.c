/* ************ gn_kr.c ************ */

/* Main procedure for lattice gross_neveu model
   in 2-dimensions, using Kramer's algorithm. */

#define CONTROL
#include <stdio.h>
#include <math.h>
#include "lattice_gn_kr.h"    /* global variables for lattice fields */


FILE *ptlat;
FILE *ptacl;
FILE *pin;
FILE *ptout;
FILE *ptcon;
FILE *ptprop;
FILE *ptpropacl;

main()
{

/* ********************************************************** */
/*         Decleration of Variables and Functions             */
/* ********************************************************** */

int prompt;
int a, ac, acc, bin, g, i, j, k, l, kr_flag, m, meas; 
int n, no_kr, no_acc, no_auto, no_prop, u, v;
int A, Acc, J, lbd, M, Meas;
double acc_rate, av_sigma, t_ex_sigma, seg_av_sigma, av_psi, seg_av_prop;
double d_av_sigma, d_av_psi, sqdev, av_t_int, d_t_int;
site *s;

void randomize(void);
void readin(int prompt);
int setup_kr(void);
double gasdev(void);
double gauss(void);
int kramer(void);
double average_sigma(void);
void autocorel(double sigma_av, int lb, int a_index);
double propagator(void);
void matp2p(field_offset src,field_offset dest,int isign,int parity,int flavor);
/*void coldlat(void);*/
/*void zerolat(void);*/
/*void coldlat2(void);*/
/*void hotlat(void);*/
void filelat(FILE *pin);

/* ************************************************************* */
/* Calling setup_kr.c for reading input parameters, `layout'ing,
   calling make_nn_gather and generating starting configuration
   for scalar field. */
/* ************************************************************* */


pin=fopen("kr1.in", "r");
ptout = fopen("kr1.out", "w");
ptcon = fopen("kr1.con", "w");
ptacl = fopen("kr1.acl", "w");
ptlat = fopen("kr1.lat", "w");
ptprop = fopen("kr1.prop", "w");
ptpropacl = fopen("kr1.propacl", "w");

prompt = setup_kr();
readin(prompt);  

randomize();

/*coldlat();*/
/*zerolat();*/
/*coldlat2();*/
/*hotlat();*/
filelat(pin);


/* ***************************************************** */ 
/* Initializing various variables, counters and indices. */
/* ***************************************************** */

av_sigma = 0.0;		/* grand average of <sigma> */
av_psi = 0.0;		/* grand average of <psi_bar-psi> */
t_ex_sigma = 0.0;	/* grand sum over all <sigma>s */
seg_av_sigma = 0.0;	/* segment-average of <sigma> for tau_int */
seg_av_prop = 0.0;	/* segment-average of <propagator[m]> */

no_kr = 0;		/* # of kramer steps calculated */
no_acc = 0;		/* # of accepted configurations */
counter = 0;		/* counter used in kramer.c */
meas = 0;		/* meas is a switch for measurements */
no_auto = 0;		/* no_auto is the index for autocoreln. measurement */ 
no_prop = 0;		/* no_prop is the index for Prop-Acl calculation */

a = 0;		/* a index used in tau_int measurements */
ac = 0;		/* ac configuration index to ac_store[ac] */
acc = 0;	/* acc index used in G_temp[][acc] i.e. # of data points 
		   and also used as index in ac_store */
 
bin = 0;	/* bin is index to bin_average */
g = 0;		/* g configuration index to garbage i.e. garbage[g] */
k = 0;		/* k configuration index to store i.e. store[k] */
l = 0;		/* l configuration index to store i.e. store[l] */
j = 0;		/* j configuration index to store[j], during measurements */


 /* *********************************************************** */
 /* switches action between garbage loops and measurement loops */
 /* *********************************************************** */

switch(sw_flag){

  case 0:
	/* *************************************************** */
   	/* *************************************************** */
	/* GARBAGE LOOPS AND AUTOCORELATION MEASUREMENTS LOOPS */
	/* *************************************************** */
    	/* *************************************************** */

    /* ********** Loop for KRAMER for kr_it times ********** */
    	
    for(j=0;j<kr_it;j++){		/* Loop for Kramer for kr_it times */
   
	FORALLSITES(i,s){
			s->mom1 = gasdev();
			for(n=0;n<nf;++n) s->eta.f[n] = gauss();
			}

   	/* chi = M^dagger *eta */
   	for(n=0;n<nf;++n){
       		matp2p(F_OFFSET(eta),F_OFFSET(chi),MINUS,EVENANDODD,n);
      		}

	for(k=0;k<KMAX;++k){

	    kr_flag = kramer();

	    if(kr_flag==0){ ++no_kr; }
	    else{
             ++no_acc;
             ++no_kr;			/* only when accepted */

             if(no_acc<=no_garbage){
                garbage[g] = average_sigma();

	        ++g;
	        ++bin;

	     if((bin%bin_length)==0){
	         for(m=(g-bin);m<g;m++) bin_av[l] += (garbage[m]/bin_length);
	         ++l;
	         bin = 0;
	         }

	     } /* garbage loop ends */

	     /* ******* Autocorelation Calculation ******* */

             if(no_acc>no_garbage){
	        ac_store[acc] = average_sigma();	

	        ++acc;
	        ++no_auto;

	        if((no_auto%seg_length)==0){
		    lbd = acc - no_auto;
	            for(m=lbd;m<acc;m++){
		        seg_av_sigma += ac_store[m]/seg_length;
		        }
	            autocorel(seg_av_sigma, lbd, a);
	            ++a;
	            seg_av_sigma = 0.0;
	            no_auto = 0;
	            }

             } /* if (no_acc>no_garbage) loop ends */

		FORALLSITES(i,s) con[i] = conf[i];
		FORALLSITES(i,s) conf[i] = s->sigma;
	
	 } /* end of else loop */ 

	 } /* k-loop ends */

    } /* kr_it loop ends */

    Acc = acc;
    M = m;
    A = a;
      
    acc_rate = (double)no_acc/(double)no_kr;

    /* ******* Calculation of mean autocorelation time and s.d.s' ******* */

    for(n=0;n<MAXT_cut;n++){
        
        fprintf(ptacl, "%d", n+1);

        for(u=0;u<NOT_cut;u++){

        d_t_int =0.0;
         av_t_int = 0.0; 
 
         for(a=0;a<A;a++) av_t_int += T_int[u][n][a]/A;
         for(a=0;a<A;a++){
	     d_t_int += ((T_int[u][n][a]-av_t_int)*(T_int[u][n][a]-av_t_int)); 
	                 }
         d_t_int = sqrt( d_t_int/(A - 1) );

         fprintf(ptacl, "\t%lf\t%lf", av_t_int, d_t_int);
			} /* u-loop ends */
	
	fprintf(ptacl, "\n");
			
        } /* n-loop ends */ 			    

    /* ********** Print Statements ********** */

    printf("no_acc_traj= %d\t no_kr= %d\t acc_rate= %lf\n\n",
                  no_acc, no_kr, acc_rate);
    printf("Acc = %d\t M = %d\t A = %d\n\n", Acc, M, A);
    for(k=0;k<no_bin;k++) fprintf(ptout, "%d\t%lf\n", k+1, bin_av[k]);
    FORALLSITES(i,s) fprintf(ptlat, "%lf\n", s->sigma);
    for(m=0;m<meas_loop;m++) fprintf(ptcon, "%lf\n", ac_store[m]);

  break;

  case 1:
	/* ************************************************** */
	/* ************************************************** */
        /* <SIGMA> AND PROPAGATOR i.e. FERMIONIC MEASUREMENTS */
	/* ************************************************** */
	/* ************************************************** */

    /* ************ Loop for KRAMER for meas_loop times ************ */


    for(l=0;l<kr_it;l++){		/* Loop for KRAMER for kr_it times */
   
	FORALLSITES(i,s){
			s->mom1 = gasdev();
			for(n=0;n<nf;++n) s->eta.f[n] = gauss();
			}

   	/* chi = M^dagger * eta */
   	for(n=0;n<nf;++n){
       		matp2p(F_OFFSET(eta),F_OFFSET(chi),MINUS,EVENANDODD,n);
      		}

	for(k=0;k<KMAX;++k){

            kr_flag = kramer();

	    if(kr_flag==0){ ++no_kr; }
	    else{
 
	     psi_acl[no_acc] = propagator();	/* for Propagator_Acl calcu */
	     for(m=0;m<nt;m++) G_temp[m][meas] = tprop[m];

             ++no_acc;
             ++no_kr;			/* only when accepted */
 
	     ++meas;
	     ++no_prop;

	     /* ******* Measurement Segment ******* */

             if((meas%meas_length)==0){
	         store[j] = average_sigma();
	         t_ex_sigma += store[j];

	         psi[j] = psi_acl[no_acc-1];
	         av_psi += psi[j];

	         for(m=0;m<nt;m++){ 
		     G_store[m][j] = tprop[m];
		     G_prop[m] += tprop[m];
			          }

	         ++j;
		 Meas = meas;

	         } /* measurement loop ends */

             /* ******* Propagator_Acl calcu ******* */

	     if((meas%prop_length)==0){
	    
	         for(m=0;m<nt;m++){
	    	     for(acc=(meas-no_prop);acc<meas;acc++){ 
		         ac_store[ac] = G_temp[m][acc];
			 ac_prop[m][acc] = G_temp[m][acc];
		         seg_av_prop += ac_store[ac]/prop_length;
		         ++ac;
		        }

		     autocorel(seg_av_prop, 0, a);

		     for(u=0;u<NOT_cut;u++)for(v=0;v<MAXT_cut;v++)
		         T_int_prop[u][v][m][a] = T_int[u][v][a]; 
		         seg_av_prop = 0.0;
		         ac = 0;
			           } /* m-loop ends */

	         no_prop = 0;
	         ++a;

                }	/* Propagator_Acl calcu ends */


            } /* end of else loop */

	    } /* k-loop ends */

    } /* kr_it loop ends */

    Acc = acc;
    J = j;
    A = a;

    /* ******* Calculation of s.d. of <sigma> ******* */
    
    av_sigma = t_ex_sigma / J;
    sqdev = 0.0;
    for(k=0;k<J;k++){
        sqdev += ((store[k] - av_sigma)*(store[k] - av_sigma));
        } 
    d_av_sigma = sqrt(sqdev/(J - 1));

    /* ******* Calculation of s.d. of the average <psi_bar-psi> ******* */

    av_psi = (av_psi/J);
    sqdev = 0.0;
    for(j=0;j<J;j++) sqdev += ( (psi[j] - av_psi) * (psi[j] - av_psi) );
    d_av_psi = sqrt(sqdev/(J - 1));

    /* ******* Calculation of s.d. of the propagator ******* */

    for(m=0;m<nt;m++){
        sqdev = 0.0;
        G_prop[m] = (G_prop[m]/J);
        for(j=0;j<J;j++) sqdev += ( (G_store[m][j] - G_prop[m]) *
				          (G_store[m][j] - G_prop[m]) );	
        tprop[m] = sqrt(sqdev/(J - 1));
        }    

    acc_rate = (double)no_acc/(double)no_kr;

    /* ******* Calculation of mean and s.d. of the Propagator-Autocorelation ******* */

    for(n=0;n<MAXT_cut;n++){
	fprintf(ptpropacl, "%d", n+1);
	for(u=0;u<NOT_cut;u++){
	    for(m=0;m<nt;m++){
	    	d_t_int = 0.0;
	    	av_t_int = 0.0;
     	    	for(a=0;a<A;a++) av_t_int += T_int_prop[u][n][m][a]/A;
		for(a=0;a<A;a++) d_t_int += ( (T_int_prop[u][n][m][a] - av_t_int) * 
					      (T_int_prop[u][n][m][a] - av_t_int) );
		d_t_int = sqrt(d_t_int/(A - 1));
		fprintf(ptpropacl, "\t%lf\t%lf", av_t_int, d_t_int);
	       }
	   }
	fprintf(ptpropacl, "\n");
       } /* n-loop ends */
 
    /* ************ Print Statements ************ */

    printf("no_acc_traj= %d\t no_kr= %d\t acc_rate= %lf\n\n",
                  no_acc, no_kr, acc_rate);
    printf(" Acc = %d\t Meas = %d\t J = %d\t A = %d\n\n", Acc, Meas, J, A);
    printf("av_sigma = %lf\t av_psi= %lf \n\n", av_sigma, av_psi);
    printf("d_av_sigma= %lf\t d_av_psi = %lf\n\n", d_av_sigma, d_av_psi);

    for(m=0;m<nt;m++){
        printf("G_prop[%d] = %lf\n", m, G_prop[m]); 
        fprintf(ptprop, "%d\t%lf\t%lf\n", m, fabs(G_prop[m]), tprop[m]);
        }

    for(n=0;n<Acc;n++){
        for(m=0;m<nt;m++){
	    fprintf(ptcon, "%lf\t", ac_prop[m][n]);
	   }
	fprintf(ptcon, "\n");
       }

    for(k=0;k<J;k++) fprintf(ptout, "%d\t%lf\n", k+1, store[k]);

    FORALLSITES(i,s) fprintf(ptlat, "%lf\n", s->sigma);

  break;


  default:
	 /* ************************ */
	 /* ************************ */
	 /* DEFAULT CASE; DO NOTHING */
	 /* ************************ */
	 /* ************************ */

	 printf("\n Wrong Switch Value: It corresponds to nothing. \n");
	 printf(" The switch_flag value must be the following:\n");
	 printf("\t switch_flag 0 for garbage and autocorln. loops.\n");
	 printf("\t switch_flag 1 for measurement loops.\n");
 
  break;

  } /* switch terminates */



fclose(ptlat);
fclose(ptacl);
fclose(pin);
fclose(ptout);
fclose(ptcon);
fclose(ptprop);
fclose(ptpropacl);


} /* end of gn_kr.c */
