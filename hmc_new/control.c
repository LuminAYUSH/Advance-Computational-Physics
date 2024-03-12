// Front code for lattice Gross-Neveu model
// in 2-dim using HMC / Kramer's algorithm

#define CONTROL
#include "include_gn.h"

int main(int argc,char *argv[])
{

  int prompt;
  int algoflagi,todo;
  int i,j,k,l,m,n,nf;
  int krflag,numtrial,numaccp;
  double avsigma,dtime;
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
  if(readflag==FRESH) coldlat();
  else if(readflag==FRESHHOT) hotlat();
  else{
    f_sigma=fopen(read_binary_sigma,"r");
    read_bin_double(f_sigma);
    fclose(f_sigma);
   }
  l=numtrial=numaccp=0;
  avsigma=0;

  for(todo=0;todo<NMdtraj;todo++){

     numtrial+=hmc();
     numaccp++;

     if((todo%Measlen)==0){
       printf("trajectory= %d\t av_sigma= %e\n",todo,fabs(average_sigma()));
       }

     } // todo-loop

  printf("\nAcceptance rate = %.4lf\n",(double)numaccp/(double)numtrial);

  // save sigma if requested
  if(saveflag!=FORGET){
    f_sigma=fopen(save_binary_sigma,"w");
    save_bin_double(f_sigma);
    fclose(f_sigma);
    }

  // Mark time
  dtime+=dclock();
  printf("2-dim KS Gross-Neveu completed!\n");
  printf("Time = %e seconds\n",dtime);
  fflush(stdout);

} // end of main

