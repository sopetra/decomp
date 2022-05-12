#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "math.h"
//treba ispisivat aproks rješenje, potreban broj koraka i norme reziduala u svakom koraku
void cg(doublereal *A, doublereal *b, doublereal *x0, doublereal tol, integer m, int *info){
  doublereal *r,*d,*pom,alfa=1,beta=0,kriterij;
  integer i,k=0,inc=1, n=m;
  doublereal nb,nr,a,b1,c1;
  char trans='N';
  nb=dnrm2_(&n,b,&inc);
  if(nb==0.0) {(*info)=1;return;}
  else (*info)=0;
  //printf("||b||=%f\n", nb);
  r=malloc(n*sizeof(doublereal));
  pom=malloc(n*sizeof(doublereal));
  d=malloc(n*sizeof(doublereal));
  //pom=A*x0 tj. pom=alfa*A*x0 + beta*pom
  dgemv_(&trans,&n,&n,&alfa,A,&n,x0,&inc,&beta,pom,&inc);
    //for(i=0; i<n; i++) printf("%f", pom[i]);
  alfa=-1;
  dcopy_(&n,b,&inc,r,&inc); //r=b
  dcopy_(&n,b,&inc,d,&inc); //d=b
  //for(i=0; i<n; i++) printf("%f", pom[i]);
  daxpy_(&n,&alfa,pom,&inc,r,&inc); //r=r+alfa*pom=b-A*x0
  daxpy_(&n,&alfa,pom,&inc,d,&inc); //d=d+alfa*pom=b-A*x0
  doublereal dot1=ddot_(&n,r,&inc,r,&inc), dot2, dot;
  while(1){
    //ZA ALGORITAM
    alfa=1;
    dgemv_(&trans,&n,&n,&alfa,A,&n,d,&inc,&beta,pom,&inc); //pom=A*d tj pom=alfa*A*d
    dot2=ddot_(&n,d,&inc,pom,&inc);
    a=dot1/dot2;

    daxpy_(&n,&a,d,&inc,x0,&inc); //x(k+1)=x(k)+a*d

    a=a*(-1);
    daxpy_(&n,&a,pom,&inc,r,&inc); //r=r-a*pom
    dot=ddot_(&n,r,&inc,r,&inc);
    b1=dot/dot1;
    dcopy_(&n,r,&inc,pom,&inc); //pom=r
    daxpy_(&n,&b1,d,&inc,pom,&inc);
    dcopy_(&n,pom,&inc,d,&inc); //d=pom

    dot1=dot;

    //ZA KRITERIJ ZAUSTAVLJANJA
    nr=sqrt(dot);
  //  printf ("||r(%ld)||=%lg\n",k+1, nr);
    kriterij = nr/nb;
    if(kriterij<=tol) break;
    k++;

  }
//  printf("potreban broj iteracija je: %ld\n", k+1);
/*  for(i=0;i<n;i++){
    printf("%lg  ", x0[i]);
  }*/
}
