#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "math.h"
#include <time.h>
#include "cg1.c"
void sustav(doublereal *B, doublereal *E, integer s, integer g, doublereal tol,
  int *info, char strana, integer granica){
  clock_t start, end;
  double cpu_time_used, cpu_time_used1;
  double otkucaji = 0.0, otkucaji0 = 0.0;


  integer i,j;
  doublereal *x0, *b;
  x0=calloc(s,sizeof(doublereal));
  b=calloc(s,sizeof(doublereal));
  if(strana=='D'){//s desne strane su ne-nul stupci
    for(j=granica;j<g;j++){
        //kopira j-ti stupac od E1 u b (desna strana)
        start = clock();
        for(i=0;i<s;i++) b[i]=E[j*s+i];
        end = clock();
        otkucaji = otkucaji + (end-start);
        start=clock();
        cg1(B,b,x0,tol,s,info);//rješavanje sustava
        //printf("\ninfo j=%d: %d\n", j, info);
        end=clock();
        otkucaji0 = otkucaji0 + (end-start);

        start = clock();
        for(i=0;i<s;i++){
          E[j*s+i]=x0[i];
          x0[i]=0.0;
          }
        end=clock();
        otkucaji = otkucaji + (end-start);
      }
  }
  else if(strana=='L'){//s lijeve strane su nam ne nul elemneti
    for(j=0;j<granica;j++){
      //kopira j-ti stupac od E1 u b (desna strana)
        start = clock();
        for(i=0;i<s;i++) b[i]=E[j*s+i];
        end = clock();
        otkucaji = otkucaji + (end-start);
        start=clock();
        cg1(B,b,x0,tol,s,info);//rješavanje sustava
        //printf("\ninfo j=%d: %d\n", j, info);
        end=clock();
        otkucaji0 = otkucaji0 + (end-start);

        start = clock();
        for(i=0;i<s;i++){
          E[j*s+i]=x0[i];
          x0[i]=0.0;
          }
          end=clock();
          otkucaji = otkucaji + (end-start);
        }
      }
  else if(granica==-1){//prolazimo cijelu maatricu(sve stupce)
    for(j=0;j<g;j++){
        //kopira j-ti stupac od E1 u b (desna strana)
        start = clock();
        for(i=0;i<s;i++) b[i]=E[j*s+i];
        end = clock();
        otkucaji = otkucaji + (end-start);
        start=clock();
        cg1(B,b,x0,tol,s,info);//rješavanje sustava
        //printf("\ninfo j=%d: %d\n", j, info);
        end=clock();
        otkucaji0 = otkucaji0 + (end-start);

        start = clock();
        for(i=0;i<s;i++){
          E[j*s+i]=x0[i];
          x0[i]=0.0;
        }
        end=clock();
        otkucaji = otkucaji + (end-start);
      }
  }

  cpu_time_used = ((double) (otkucaji)) / CLOCKS_PER_SEC;
  cpu_time_used1 = ((double) (otkucaji0)) / CLOCKS_PER_SEC;

  printf("\n\n\nvrijeme cg-a:%lf   \n", cpu_time_used);
  printf("\n\n\nvrijeme prekopiranja:%lf   \n", cpu_time_used1);

}
