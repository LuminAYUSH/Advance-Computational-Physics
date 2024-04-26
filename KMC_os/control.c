// Front code for lattice Gross-Neveu model
// in 2-dim using Kramer's algorithm

#define CONTROL
#include "include_gnkr.h"

int main(int argc,char *argv[])
{

  int prompt;
  int i,j,k,l,m,n,nf;
  int krflag,numkr,numaccp;
  double avsigma,sdsigma;
  double dtime;
  site *s;

  FILE *f_sigma;
  FILE *fsig,*fpbp;

  // Initialize
  dtime=-dclock();
  prompt=setup();
  readin(prompt);

  // malloc storages
  Nmeas=(int)((KRtraj*KRmax)/Measlen);
  av_sigma=(double *)malloc(Nmeas*sizeof(double));

  /*****************************************/
  // read binary sigma fields, if requested
  // else populate it with flat / gaussian RN
/*
  if(readflag==FRESH) coldlat();
  else if(readflag==FRESHHOT) hotlat();
  else{
    f_sigma=fopen(read_binary_sigma,"r");
    read_bin_double(f_sigma);
    fclose(f_sigma);
   }
*/
  coldlat();
  l=numkr=numaccp=0;
  avsigma=0;

  for(j=0;j<KRtraj;j++){

     // Init momentum and pseudofermions
     FORALLSITES(i,s){
	s->mom=gsl_ran_gaussian(r,1.0);
	for(nf=0;nf<Nf;nf++) s->zeta.f[nf]=gsl_ran_gaussian(r,1/sqrt(2));
	}

     for(k=0;k<KRmax;){
	krflag=kramer();
	if(krflag==0){	// rejected
	  ++numkr;
	  }
	else{		// accepted
	  ++numkr; ++numaccp; ++k;
	  }
	} // krmax-loop

     // printf("Kriter= %d\t av_sigma= %e\n",j,fabs(average_sigma()));

     } // krtraj-loop
  printf("\nAcceptance rate = %.4lf\n",(double)numaccp/(double)numkr);

  // save sigma if requested
/*
  if(saveflag!=FORGET){
    f_sigma=fopen(save_binary_sigma,"w");
    save_bin_double(f_sigma);
    fclose(f_sigma);
    }
*/

  // Mark time
  dtime+=dclock();
  // printf("2-dim KS Gross-Neveu completed!\n");
  printf("Time = %e seconds\n",dtime);
  printf("The Value of coupling = %.2lf\n",G);
  printf("The Value of step_size = %.4lf\n",Step);
  fflush(stdout);

} // end of main

