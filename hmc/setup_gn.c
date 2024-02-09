/******** setup_gn.c *********/
/*  set tabstop=2   for easy reading of this file */
/* alpha code, version 2 */

#include <stdio.h>
#include <math.h>
#include "lattice_gn.h"

int  setup_gn(void)   {

void make_lattice(void),make_nn_gathers(void),layout(void);
int prompt,initial_set(void);

	/* print banner, get lattice size and volume */
    prompt=initial_set();
    
        /* layout the lattice */
    layout();    

	/* allocate space for lattice, set up coordinate fields */
    make_lattice();
    
	/* set up nearest neighbor gathers */
    make_nn_gathers();

    return(prompt);
}  /* end of setup_gn */


/* FUNCTION INITIAL_SET */

int initial_set(void)   {

int prompt,status,get_i(int,char *),getprompt(int *);
float get_f(int,char *);

        /* print banner */
        printf("GN model with HMC algorithm\n");
        printf("Alpha machine, Version 1\n");

        printf(
	    "type 0 for no prompts, 1 for prompts or 2 for list of prompts\n");
        status=getprompt(&prompt);
        if (status != 0) 
	    {printf("error in input: initial prompt\n");return(-1);}

        nx=get_i(prompt,"nx");
        nt=get_i(prompt,"nt");
        if(nx%2 !=0 || nt%2 !=0) 
        {
            printf("nx,nt must be even!! \n"); 
            terminate(1);
	}

        printf("lattice dimensions = %d %d\n",nx,nt);

	/* Switch flag */
	sw_flag=get_i(prompt,"switch_flag");
	printf("Switch_Flag = %d\n", sw_flag);

	/* Number of measurements */
	no_garbage=get_i(prompt,"no_of_garbage_loops");
	printf("# of garbage = %d\n",no_garbage);

	/* the length of each bin */
	bin_length=get_i(prompt,"bin_length");
	printf("bin_length = %d\n", bin_length);

	/* Number of HMC iterations, only accepted ones count */
	hmc_it=get_i(prompt,"no_of_hmc_iterations");
	printf("# of hmc iterations = %d\n",hmc_it);

	/*  length, the length after which the fermionic
	   observables are measured. */
	meas_length=get_i(prompt,"meas_length");
	printf("meas_length = %d\n",meas_length);

	/*  length, the length after which the propagator
	   autocorelations are calculated. */
	prop_length=get_i(prompt,"prop_length");
	printf("prop_length = %d\n",prop_length);

	/* segment length, the length after which the auto
           corelations are measured. */
	seg_length=get_i(prompt,"seg_length");
	printf("seg_length = %d\n",seg_length);

    volume = nx * nt;
    mid = nt/2;
    meas_loop = hmc_it - no_garbage;
    if(sw_flag==1) seg_length = prop_length;

    return(prompt);

}  /* end of initial_set */

/* FUNCTION MAKE_LATTICE */

void make_lattice(void){
register int i,j,k;		/* dummy */
int x,t;		/* coordinates */

 no_bin = no_garbage/bin_length;
 no_meas = meas_loop/meas_length;
 no_a_seg = meas_loop/seg_length;
 no_prop_seg = meas_loop/prop_length;

    /* allocate space for lattice */
    lattice = (site *)malloc( volume * sizeof(site) );
    if(lattice==NULL){
	printf("no room for lattice\n");
	terminate(1);
    }

    /* allocate space for store, temporarily used to store sigma average
       per accepted trajectory */
    store = (double *)malloc( no_meas * sizeof(double) );
    if(store==NULL){
	printf("no room for store.\n");
	terminate(1);
    }

    /* allocate space for conf, used to store sigma config. generated in
       the last accepted loop. */
    conf = (double *)malloc( volume * sizeof(double) );
    if(conf==NULL){
	printf("no room for conf.\n");
	terminate(1);
    }  	

    /* allocate space for con, used to store sigma config. generated in
       the last but one accepted loop. */
    con = (double *)malloc( volume * sizeof(double) );
    if(con==NULL){
	printf("no room for con.\n");
	terminate(1);
    }  	

    /* allocate space for garbage, temporarily used to store sigma average
       per garbage trajectory */
    garbage = (double *)malloc( no_garbage * sizeof(double) );
    if(garbage==NULL){
	printf("no room for garbage.\n");
	terminate(1);
    }

    /* allocate space for ac_store, temporarily used to store sigma average
       per autocorelation trajectory */
    ac_store = (double *)malloc( meas_loop * sizeof(double) );
    if(ac_store==NULL){
	printf("no room for ac_store.\n");
	terminate(1);
    }

    /* allocate space for ac_prop, temporarily used to store sigma average
       per autocorelation trajectory */
    ac_prop = (double *)malloc( prop_length * sizeof(double) );
    if(ac_prop==NULL){
	printf("no room for ac_prop.\n");
	terminate(1);
    }

    /* allocate space for bin_av, used to store bin average
       per bin_length accepted trajectory */
    bin_av = (double *)malloc( no_bin * sizeof(double) );
    if(bin_av==NULL){
	printf("no room for bin_av.\n");
	terminate(1);
    }

    /* allocate space for psi, used to store psi_bar-psi */
    psi = (double *)malloc(no_meas * sizeof(double) );
    if(psi==NULL){
        printf("no room for psi. \n");
        terminate(1);
    } 

    /* allocate space for psi_acl, used to store psi_bar-psi */
    psi_acl = (double *)malloc(meas_loop * sizeof(double) );
    if(psi_acl==NULL){
        printf("no room for psi_acl. \n");
        terminate(1);
    } 

    /* allocate space for G_prop, the storage place for propagator. */
    G_prop = (double *)malloc(nt * sizeof(double) );
    if(G_prop==NULL){
       printf("no room for G_prop. \n");
       terminate(1);
    }

    /* allocate space for G_store, the storage place for propagator. */
    for(i=0;i<nt;i++){
        G_store[i] = (double *)malloc(no_meas * sizeof(double) );
        if(G_store[i]==NULL){
           printf("no room for G_store. \n");
           terminate(1);}
        }

    /* allocate space for G_temp, the storage place for propagator for
      Acl calculation. */
    for(i=0;i<nt;i++){
        G_temp[i] = (double *)malloc(meas_loop * sizeof(double) );
        if(G_temp[i]==NULL){
           printf("no room for G_temp. \n");
           terminate(1);}
        }

    /* allocate space for prop, the working place for the calculation
       of propagator */
    prop = (double *)malloc(nt * sizeof(double) );
    if(prop==NULL){
       printf("no room for prop. \n");
       terminate(1);
    }

    /* allocate space for prop, the working place for the calculation
       of propagator */
    tprop = (double *)malloc(nt * sizeof(double) );
    if(tprop==NULL){
       printf("no room for tprop. \n");
       terminate(1);
    }

    /* Allocate address vectors */
    gen_pt = (char **)malloc( volume * sizeof(char *) );
    if(gen_pt==NULL){
       printf("no room for pointer vector\n");
       terminate(1);
     }

    /* allocate space for T_int, the storage place for tau_int for
       various measurements one per segment. */
    for(i=0;i<NOT_cut;i++){
        for(j=0;j<MAXT_cut;j++){
            T_int[i][j] = (double *)malloc(no_a_seg * sizeof(double) );
            if(T_int[i][j]==NULL){
               printf("no room for T_int. \n");
               terminate(1);}
           }
       }

    /* allocate space for T_int_prop, the storage place for tau_int for
       various Prop-Acl measurements one per segment. */
    for(i=0;i<NOT_cut;i++){
        for(j=0;j<MAXT_cut;j++){
	    for(k=0;k<nt;k++){
                T_int_prop[i][j][k] = (double *)malloc(no_prop_seg * sizeof(double) );
                if(T_int_prop[i][j][k]==NULL){
                printf("no room for T_int_prop. \n");
                terminate(1);}
	       }
           }
       }


    for(t=0;t<nt;t++)for(x=0;x<nx;x++){
	    i=site_index(x,t);
	    lattice[i].x=x; lattice[i].t=t;
	    if ( t%2 == 0 )lattice[i].sign = 1; 
	    else lattice[i].sign = -1;
	    	    if( (x+t)%2 == 0)lattice[i].parity=EVEN;
	    else	         lattice[i].parity=ODD;
	    
    }

    for(i=0;i<meas_loop;i++) ac_store[i] = 0.0;
    for(i=0;i<no_bin;i++) bin_av[i] = 0.0;
    for(i=0;i<nt;i++) G_prop[i] = 0.0;
    for(i=0;i<NOT_cut;i++){
        for(j=0;j<MAXT_cut;j++){
	    for(k=0;k<no_a_seg;k++){
		T_int[i][j][k] = 0.0; }}}


}  /* end of make_lattice */


/* FUNCTIONS GET_F, GET_I & GETPROMPT */

/* get_f is used to get a floating point number.  If prompt is non-zero,
it will prompt for the input value with the variable_name_string.  If
prompt is zero, it will require that variable_name_string preceed the
input value.  get_i gets an integer.
get_i and get_f return the values, and exit on error */
/* getprompt gets the initial value of prompt */

float get_f(prompt,variable_name_string)
	int prompt;
	char *variable_name_string;
{
int s;
float x;
char checkname[80];
	if(prompt)  {
		printf("enter %s ",variable_name_string);
		s=scanf("%e",&x);
		if(s == 1) return(x);
	}
	else  {
		s=scanf("%s%e",checkname,&x);
		if (s == EOF) terminate(0);
		if(s == 2 && strcmp(checkname,variable_name_string) == 0)
		   return(x);
	}
	printf("error in input: %s\n",variable_name_string);
	terminate(1);
}
int get_i(prompt,variable_name_string)
	int prompt;
	char *variable_name_string;
{
int s,i;
char checkname[80];
	if(prompt)  {
		printf("enter %s ",variable_name_string);
		s=scanf("%d",&i);
		if (s == 1) return(i);
	}
	else  {
		s=scanf("%s%d",checkname,&i);
		if (s == EOF) terminate(0);
		if(s == 2 && strcmp(checkname,variable_name_string) == 0)
		   return(i);
	}
	printf("error in input: %s\n",variable_name_string);
	terminate(1);
}
int getprompt(prompt) int *prompt; {
char initial_prompt[80];
int stat;
void printprompts(void);

	scanf("%s",initial_prompt);
	if(strcmp(initial_prompt,"prompt") == 0)  {
	   stat=scanf("%d",prompt);
	   if (stat != 1) return(1);
	}
	else if(strcmp(initial_prompt,"0") == 0)
	   *prompt=0;
	else if(strcmp(initial_prompt,"1") == 0)
	   *prompt=1;
	else if(strcmp(initial_prompt,"2") == 0) {
	   *prompt=1; printprompts();
	}
	else return(1);

	return(0);
}

void printprompts(void) {}
/* printf("Here is the list of prompts in the order they are to appear.\n");
printf("A choice among keywords in denoted by {  }.\n");
printf("Optional filenames are only needed if reloading or saving a lattice.\n");
printf("\t prompt\n\t nx\n\t ny\n\t nz\n\t nt\n");
printf("\t beta\n\t kappa\n\t approx_kappa_c\n");
printf("\t max_iterations\n\t error_for_propagator\n"); 
printf("\t 1 for psi !=0, 0 for psi=0\n"); 
printf("\t number_hopping_steps\n"); 
printf("\t wallflag\n\t width^(-2)\n"); 
printf("\t prompt\n\t source_x\n\t source_y\n\t source_z\n\t source_t\n");
printf("\t {continue, fresh, reload, reload_binary} [for lattice]\n");
printf("\t [filename]\n");
printf("\t {continue, fresh, reload, reload_binary} [for propagator]\n");
printf("\t [filename]\n");
printf("\t {forget, save, save_binary} [for propagator]\n");
printf("\t [filename]\n");
printf("\t {forget, save, save_binary} [for mesons]\n");
printf("\t [filename]\n");
}
*/

/* read in parameters and coupling constants	*/

void readin(prompt) int prompt;  {
/* read in parameters for HMC on GN model	*/
/* argument "prompt" is 1 if prompts are to be given for input*/


int status;
float x,get_f();
int get_i();
/* char savebuf[20];
void save_ascii(),restore_ascii();
void save_binary(),restore_binary();
void coldlat(); */
double dtime;


        printf("\n\n");

	/* get coupling "g" and no. of flavor "nf"	*/
        g=(double)get_f(prompt,"g");
        nf=get_i(prompt,"nf");

	printf("g = %lf, nf = %d\n",g,nf);
	
	
	/* Number of HMC iterations, only accepted ones count */
	/*hmc_it=get_i(prompt,"no_of_hmc_iterations");
	printf("no of hmc iterations = %d\n",hmc_it);*/
	
	
	/* Molecular dynamics parameters */
	mdstep=get_i(prompt,"no_of_md_steps");
	step=(double)get_f(prompt,"step_size");
	printf("no of md steps = %d, step size = %lf\n",mdstep,step);
    
        /* maximum no. of conjugate gradient iterations */
        cgiter1=get_i(prompt,"max_cg_iterations_for_hamil");
        cgiter2=get_i(prompt,"max_cg_iterations_for_piup");
        printf("maximum no. of conj. grad. iterations = %d,%d\n",
        cgiter1,cgiter2);
    
    
        /* Max. error for conjugate gradients */
        residue1=(double)get_f(prompt,"residue_for_cg_hamil");
        residue2=(double)get_f(prompt,"residue_for_cg_piup");
        printf("residues for conjugate grad = %e,%e\n",residue1,residue2);

        
   
}/* end of readin */
