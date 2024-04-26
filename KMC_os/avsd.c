#include <stdio.h>
#include <math.h>

#define N 1800
#define In 200

int main(int argc,char *arv[])
{

  int i,j;
  double sigma[N],avsigma,sdsigma;
  char trash[256];

  for(i=0;i<N;i++) scanf("%s %d %s %lf",trash,&j,trash,&sigma[i]);

  avsigma=sdsigma=0.0;
  for(i=In;i<N;i++) avsigma+=sigma[i];
  avsigma/=(double)(N-In);
  for(i=In;i<N;i++) sdsigma=+(sigma[i]-avsigma)*(sigma[i]-avsigma);
  sdsigma=sqrt(sdsigma/(double)(N-In));

  printf("%lf    %lf\n",avsigma,sdsigma);

}
