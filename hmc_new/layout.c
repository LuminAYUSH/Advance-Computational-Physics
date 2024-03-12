/****** layout.c ******/

#include "include_gn.h"

site **neighbor[Ndirs];

/////////////////////////
void make_nn_gathers(void)
{

  int i,j,dir;
  int xpt,tpt;
  site *s;

  for(dir=0;dir<Ndirs;dir++){
     neighbor[dir]=(site **)malloc(Volume*sizeof(site *));
     if(neighbor[dir]==NULL){
	printf("No room for neighbor\n"); terminate(1);
	}
     }

  for(dir=Xup;dir<Ndirs;dir++){
     FORALLSITES(i,s){
	nn_coords(s->x,s->t,dir,&xpt,&tpt);
	j=site_index(xpt,tpt);
	neighbor[dir][i] = &(lattice[j]);
	}
     }

} // end of make_nn_gathers() 

/////////////////////////
void gather(field_offset field,int dir,char **dest)
{

  int i;
  site *s;
  FORALLSITES(i,s){ dest[i]=F_PT(neighbor[dir][i],field); }

} // end of gather()

/////////////////////////
int site_index(int x,int t)
{

  int xr, tr;
  xr= x%N; tr=t%N;
  return(xr+tr*N);

} // end of site_index ()

/////////////////////////
void nn_coords(int x,int t,int dir,int *xn,int *tn)
{

  *xn=x; *tn=t;
  switch(dir){
	case Xup: *xn=(x+1)%N; break;
	case Xdn: *xn=(x+N-1)%N; break;
	case Tup: *tn=(t+1)%N; break;
	case Tdn: *tn=(t+N-1)%N; break;
	default: printf("Bad direction\n"); terminate(1);
	}

} // end of nn_coords()

/////////////////////////
void make_lattice(void)
{

  int i,j,k;
  int x,t;

  lattice = (site *)malloc(Volume*sizeof(site));
  if(lattice==NULL){
    printf("No room for lattice.\n"); terminate(1);
    }

  gen_pt=(char **)malloc(sites_on_node*sizeof(char *));
  if(gen_pt==NULL){
    printf("no room for pointer vector\n"); terminate(1);
    }

  for(x=0;x<Nx;x++)for(t=0;t<Nt;t++){
     i=site_index(x,t);
     lattice[i].x=x; lattice[i].t=t;
     lattice[i].index=x+Nx*t;
     if(t%2==0) lattice[i].sign= 1;
     else	lattice[i].sign=-1;
     }

} // end of make_lattice()

/////////////////////////
double dclock(void)
{

  long fine;
  fine=clock();
  return(((double)fine)/1000000.0);

} // end of dclock()


/////////////////////////
void terminate(int status)
{
  printf("Termination: status = %d\n",status);
  exit(status);
} // end of dclock()

/////////////////////////
void setup_layout(void)
{
  sites_on_node = Volume;
} // end of setup_layout()

