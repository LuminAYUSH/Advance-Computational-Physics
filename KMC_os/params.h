#ifndef _PARAMS_H
#define _PARAMS_H

typedef struct {

        int stopflag;
        int iseed;

        int nx,nt,nflv;
	double g;

	int step;
	double eps;
	int krmax,krtraj;

	int cgiter1,cgiter2;
	double residue1,residue2;

	int read_sigma,save_sigma;
        char read_binary_double[MAXFILENAME];
        char save_binary_double[MAXFILENAME];

} params;

#endif /* _PARAMS_H */
