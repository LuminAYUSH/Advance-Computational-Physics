/* setup lattice and read input parameters */

#include "include_gn.h"

#define IF_OK if(status==0)

int initial_set(void);

//////////////////////////
int setup(void)
{

  int prompt;

  prompt=initial_set();
  setup_layout();
  make_lattice();
  make_nn_gathers();

  gsl_rng_env_setup();
  T=gsl_rng_default;
  r=gsl_rng_alloc(T);

  return(prompt);

} /* end of setup() */

//////////////////////////
int initial_set(void)
{

  int prompt,status;

  printf("Gross-Neveu with Kramers algorithm for KS fermion\n");
  printf("MILC version 7.7.11\n");

  status=get_prompt(&prompt);
  if(status!=0){
    printf("error in input: initial prompt\n"); terminate(1);
    }
  Nx=get_i(prompt,"nx");
  Nt=get_i(prompt,"nt");
  printf("Lattice dimension -- %d %d\n",Nx,Nt);
  
  if(Nx==Nt) N=Nx;
  else{ printf("Nx must be equal to Nt.\n"); terminate(1); }
  Volume=Nx*Nt;
    
  return(prompt);

} /* End of initial_set() */

//////////////////////////
int readin(int prompt)
{

  int status=0;

  printf("\n");
  
  Nf=get_i(prompt,"num_flavors");
  G=get_lf(prompt,"coupling");
  Step=get_lf(prompt,"step_size");
  KRtraj=get_i(prompt,"kr_trajectory_length");
  KRmax=get_i(prompt,"kr_refreshing_steps");
  Measlen=get_i(prompt,"measurement_interval");
  CgiterH=get_i(prompt,"max_cg_iter_for_hamil");
  ResidueH=get_lf(prompt,"residue_cg_hamil");
  CgiterP=get_i(prompt,"max_cg_iter_for_piup");
  ResidueP=get_lf(prompt,"residue_cg_piup");

/*
  readflag=saveflag=99;
  status=get_s(prompt,"read_config",read_binary_sigma);
  if(strcmp(read_binary_sigma,"fresh")==0) readflag=FRESH;
  else if(strcmp(read_binary_sigma,"freshhot")==0) readflag=FRESHHOT;
  status=get_s(prompt,"save_config",save_binary_sigma);
  if(strcmp(read_binary_sigma,"forget")==0) saveflag=FORGET;
*/

  return(0);
    
} /* end of readin() */

