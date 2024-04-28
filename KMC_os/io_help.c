// ripped off from generic/io_helpers.c since compiling it needs bunch
// of sub-routine that are irrelavant for 2-dim
// get_i(), get_lf(), get_f(), get_sn(), check_read(),
// get_check_tag(), get_next_tag, get_prompt()

#include "include_gnkr.h"

//////////////////////////
int get_i(int prompt,char *var_string)
{

  int s,i;
  char checkname[80];
  char myname[]="get_i";

  if(prompt==1){
    printf("enter %s ",var_string);
    s=scanf("%d",&i);
    if(s==1) return(i);
    }
  else{
    s=scanf("%s %d",checkname,&i);
    if(s==EOF) terminate(0);
    if(s==2 && strcmp(checkname,var_string)==0){
      // printf("%s %d\n",var_string,i);
      return(i);
      }
    }

  printf("error in input: %s\n",var_string);
  terminate(1);

} // end of get_i()

//////////////////////////
double get_lf(int prompt,char *var_string)
{

  int s;
  double a;
  char checkname[80];
  char myname[]="get_lf";

  if(prompt==1){
    printf("enter %s ",var_string);
    s=scanf("%lf",&a);
    if(s==1) return(a);
    }
  else{
    s=scanf("%s %lf",checkname,&a);
    if(s==EOF) terminate(0);
    if(s==2 && strcmp(checkname,var_string)==0){
      // printf("%s %lf\n",var_string,a);
      return(a);
      }
    }

  printf("error in input: %s\n",var_string);
  terminate(1);

} // end of get_lf

//////////////////////////
int get_s(int prompt,char *var_string,char *str)
{

  int s;
  char checkname[80];

  if(prompt==1){
    printf("enter %s ",var_string);
    s=scanf("%s",str);
    if(s==1){
      printf("%s %s\n",var_string,str);
      return(1);
      }
    else{
      printf("Data format error.\n");
      return(0);
      }
    }
  else{
    s=scanf("%s %s",checkname,str);
    if(s==EOF) terminate(0);
    if(s==2 && strcmp(checkname,var_string)==0){
      // printf("%s %s\n",var_string,str);
      return(1);
      }
    }

} // end of get_s

//////////////////////////
int get_prompt(int *prompt)
{

  char initial_prompt[512];
  int status;
  char myname[]="get_prompt";

  // status=scanf("%s",initial_prompt);
  // if(strcmp(initial_prompt,"prompt")==0){
  //   status=scanf("%d",prompt);
  //   printf("%s %d\n",initial_prompt,*prompt);
  //   if(status!=1) terminate(1);
  //   }
  // else if(strcmp(initial_prompt,"0")==0)
  //   *prompt=0;
  // else if(strcmp(initial_prompt,"1")==0)
  //   *prompt=1;
  // else return(1);
  *prompt=0; 

  return(0);

} // end of get_prompt()

