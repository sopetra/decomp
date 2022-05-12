#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "math.h"
#include <time.h>
#include "cg.c"
void f(int*p){
  *p=10;
}
int main(){
  integer n=10,i,j;
  doublereal *A;
  A=calloc(n*n, sizeof(doublereal));
  // +1 i -1 dijagonala
  for(i=0;i<n-1;i++){
    A[(i+1)*n+i]=-1;
    A[i*n+(i+1)]=-2;
  }

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      printf("%.1lf   ", A[j*n+i]);
    }
    printf("\n");
  }
  int info = 0;
  printf("\ninfo prvi:%d", info);

  f(&info);
  printf("\ninfo:%d", info);

  i=100000000;
  clock_t start, end;
  double cpu_time_used;
  start = clock();
  while(i>0){
    i--;
  }
  end = clock();
  double otkucaji = end - start;
//trebalo bi ispast duplo
  j=100000000;
  end=clock();
  start = clock();
  while(j>0){
    j--;
  }
  end=clock();
  otkucaji = otkucaji + (end - start);
  cpu_time_used = ((double) (otkucaji) )/ CLOCKS_PER_SEC;
  printf("\n\n\nvrijeme sustava:%lf   \n", cpu_time_used);

  return 0;
}
