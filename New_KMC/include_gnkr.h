
#include "lattice.h"

// setup.c
int setup(void);
int initial_set(void);
int readin(int prompt);

// layout.c
void make_nn_gathers(void);
void gather(field_offset field,int dir,char **dest);
int site_index(int x,int t);
void nn_coords(int x,int t,int dir,int *xn,int *tn);
void make_lattice(void);
double dclock(void);
void terminate(int status);
void setup_layout(void);

// dslash_gn.c
void matp2p(field_offset src,field_offset dest,int isign,int flav);
void matp2d(field_offset src,field_offset dest,int isign,int flav);
void matd2p(field_offset src,field_offset dest,int isign,int flav);
void matd2d(field_offset src,field_offset dest,int isign);

// hamiltonian.c
double hamil(int hflag,int flag);

// piup.c
double piup(double tau);

// kramer.c
int kramer(void);

// measurement.c
double average_sigma(void);
double average_pbp(void);

// invert.c
void congrad(field_offset src,field_offset dest,int cgiter,
	   double residue,int cgflag,int flav);
// io_sigma.c
void read_bin_double(FILE *f_bin);
void save_bin_double(FILE *f_bin);
void coldlat(void);
void hotlat(void);

// io_help.c
int get_prompt(int *prompt);
double get_lf(int prompt,char *var_string);
int get_i(int prompt,char *var_string);
int get_s(int prompt,char *var_string,char *str);

