#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "math.h"
#include <time.h>
#include "cg.c"
#include "sustav.c"

void print(doublereal *A, integer m, integer n){
  printf("\n");
  integer i, j;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++) printf("%.1lf   ", A[j*m+i]);
    printf("\n");
  }
  printf("\n");
}
void Poisson_L(integer nx, integer ny){
  clock_t start1, end1;
  double cpu_time_used1;
  doublereal *A, *B1, *B2, *B3, *E1, *E2, *E3, *F, *C, *E;
  doublereal *f1, *f2, *f3, *gvec, *f;
  doublereal x, y, hx, hy;
  integer i, j, k, l, t, n;
  integer j1, j2, j3, i1, i2, i3;
  integer s1, s2, s3, k1, k2, g;

  s1=s2=s3=0;
  k1=0; k2=0; k=0; l=0;

  x=0.0; y=0.0;
  hx=2.0/(nx+1); hy=3.0/(ny+1);
  //printf("%lf, %lf", hx, hy);
  //poslije probat ovo efikasnije napravit bez prolaska kroz domenu

  for (i=1; i<=nx; i++){ //koliko je x-eva < 1
    x = i*hx;
    //printf ("%lf\n", x);
    if (x<1.0) k++;
    if (x==1.0){
      k1++; break;
    }
    if (x>1.0) break;
  }

  for (j=1; j<=ny; j++){
    y = j*hy;
    //printf ("%lf\n", y);
    if (y<1.0) l++;
    if (y==1.0){
      k2++; break;
    }
    if (y>1.0) break;
  }
  k1*=l;
  k2*=k;
  printf("\nk1 = %d, k2 = %d\n",k1,k2 );
start1=clock();
  i1=k;
  if(x==1.0) i2=nx-k-1;
  else i2=nx-k;

  j1=l;
  if(y==1.0) j3=ny-l-1;
  else j3=ny-l;
  j2=j1;
  i3=i1;

  s1=i1*j1; s2=i2*j2; s3=i3*j3;
  printf("s1 = %d, s2 = %d, s3 = %d",s1, s2, s3 );
  printf("\nj1 = %d, j2 = %d, j3 = %d",j1, j2, j3 );
  printf("\ni1 = %d, i2 = %d, i3 = %d",i1, i2, i3 );

  n = s1+s2+s3+k1+k2;
  g=k1+k2;
  A=calloc(n*n, sizeof(doublereal));

  B1=calloc(s1*s1, sizeof(doublereal));
  B2=calloc(s2*s2, sizeof(doublereal));
  B3=calloc(s3*s3, sizeof(doublereal));

  E1=calloc(s1*g, sizeof(doublereal));
  E2=calloc(s2*g, sizeof(doublereal));
  E3=calloc(s3*g, sizeof(doublereal));

  F=calloc((s1+s2+s3)*g, sizeof(doublereal));
  gvec=calloc(g,sizeof(doublereal));

  j=0; l=0; t=0;
  for(i=0;i<s1;i++){
    B1[i*s1+i]=2*(pow(nx+1,2)+pow(ny+1,2));
    //desni rub od B1 upada u gama1
    if(i%i1==0) B1[i*s1+i+1]=-pow(nx+1,2);
    else if((i+1)%i1==0 && k1!=0){
      B1[i*s1+i-1]=-pow(nx+1,2);
      E1[j*s1+i]=-pow(nx+1,2);//isti redak j-ti stupac
      F[i*g+j]=-pow(nx+1,2);//j-ti redak
      j++;
    }
    else B1[i*s1+i+1]=B1[i*s1+i-1]=-pow(nx+1,2);//ovde ne popunjavamo E

    if(i>=i1) B1[i*s1+i-i1]=-pow(ny+1,2);
    if(i<s1-i1) B1[i*s1+i+i1]=-pow(ny+1,2);

    //gornji rub od B1 ima za susjede gama2
    if(i>=s1-i1){//isti red j-ti stupac (j di smo stali od prije)
      if(t==0){l=j+1; t=1;}
      E1[l*s1+i]=-pow(ny+1,2);
      F[i*g+l]=-pow(ny+1,2);//l-ti redak
      l++;
    }
  }


  int gr2, gr3;
  j=0;
  for(i=0;i<s2;i++){
    B2[i*s2+i]=2*(pow(nx+1,2)+pow(ny+1,2));
    //lijevi rub upada u gama1
    if(i%i2==0 && k1!=0){
      B2[i*s2+i+1]=-pow(nx+1,2);
      E2[j*s2+i]=-pow(nx+1,2);
      F[(i+s1)*g+j]=-pow(nx+1,2);
      j++;
    }

    else if((i+1)%i2==0) B2[i*s2+i-1]=-pow(nx+1,2);
    else B2[i*s2+i+1]=B2[i*s2+i-1]=-pow(nx+1,2);
    if(i>=i2) B2[i*s2+i-i2]=-pow(ny+1,2);
    if(i<s2-i2) B2[i*s2+i+i2]=-pow(ny+1,2);
    gr2=j;
  }


  j=k1; t=0; //točke iz gama1 ne upadaju kao susjedi u ovu mtricu pa j ide od k1
  for(i=0;i<s3;i++){
    B3[i*s3+i]=2*(pow(nx+1,2)+pow(ny+1,2));

    if(i%i3==0) B3[i*s3+i+1]=-pow(nx+1,2);
    else if((i+1)%i3==0) B3[i*s3+i-1]=-pow(nx+1,2);
    else B3[i*s3+i+1]=B3[i*s3+i-1]=-pow(nx+1,2);

    if(i>=i3) B3[i*s3+i-i3]=-pow(ny+1,2);
    if(i<s3-i3) B3[i*s3+i+i3]=-pow(ny+1,2);
    //donji rub od B3 ima za susjede gama2
    if(i<i3 && k2!=0){//isti red j-ti stupac (j od 0 za gama2)
      if(t==0) {gr3=j; t=1;}
      E3[j*s3+i]=-pow(ny+1,2);
      F[(i+s2+s1)*g+j]=-pow(ny+1,2);
      j++;
    }
  }
  printf("\ngr2= %d, gr3= %d\n", gr2, gr3);

  C=calloc(g*g, sizeof(doublereal));
  for(i=0;i<g;i++){ //2 bloka tridijagonalnih matrica (dimenzija k1 i k2)
    C[i*g+i]=2*(pow(nx+1,2)+pow(ny+1,2));
    if(i==k1 || i==0) C[i*g+i+1]=-pow(ny+1,2);
    else if(i+1==k1 || i+1==g) C[i*g+i-1]=-pow(ny+1,2);
    else C[i*g+i+1]=C[i*g+i-1]=-pow(ny+1,2);

  }
printf("\nB1\n");print(B1,s1,s1);
printf("\nE1\n");print(E1,s2,g);


printf("\nB2\n");print(B2,s2,s2);
printf("\nE2\n");print(E2,s2,g);


printf("\nB3\n");print(B3,s3,s3);
printf("\nE3\n");print(E3,s3,g);

printf("\nF\n");print(F,g,s1+s2+s3);
  end1=clock();
  cpu_time_used1 = ((double) (end1-start1)) / CLOCKS_PER_SEC;

  printf("\n\n\nVRIJEME:%lf   \n", cpu_time_used1);

  j=s1+s2+s3;l=0; t=0;
  for(i=0;i<s1;i++){
    A[i*n+i]=2*(pow(nx+1,2)+pow(ny+1,2));

    if(i%i1==0) A[i*n+i+1]=-pow(nx+1,2);
    else if((i+1)%i1==0 && k1!=0){
      A[i*n+i-1]=-pow(nx+1,2);
      A[j*n+i]=-pow(nx+1,2);
      A[i*n+j]=-pow(nx+1,2);
      j++;
    }
    else A[i*n+i+1]=A[i*n+i-1]=-pow(nx+1,2);

    if(i>=i1) A[i*n+i-i1]=-pow(ny+1,2);
    if(i<s1-i1) A[i*n+i+i1]=-pow(ny+1,2);

    if(i>=s1-i1 && k2!=0){
      if(t==0){l=j+1; t=1;}
      A[l*n+i]=-pow(ny+1,2);
      A[i*n+l]=-pow(ny+1,2);
      l++;
    }

  }

  j=s1+s2+s3;
  for(i=s1;i<s1+s2;i++){
    A[i*n+i]=2*(pow(nx+1,2)+pow(ny+1,2));

    if((i-s1)%i2==0 && k1!=0){
      A[i*n+i+1]=-pow(nx+1,2);
      A[j*n+i]=-pow(nx+1,2);
      A[i*n+j]=-pow(nx+1,2);
      j++;
    }
    else if((i+1-s1)%i2==0) A[i*n+i-1]=-pow(nx+1,2);
    else A[i*n+i+1]=A[i*n+i-1]=-pow(nx+1,2);

    if(i-s1>=i2) A[i*n+i-i2]=-pow(ny+1,2);
    if(i-s1<s2-i2) A[i*n+i+i2]=-pow(ny+1,2);
  }

  j=s1+s2+s3+k1;
  for(i=s1+s2;i<s1+s1+s3;i++){
    A[i*n+i]=2*(pow(nx+1,2)+pow(ny+1,2));

    if((i-s1-s2)%i3==0) A[i*n+i+1]=-pow(nx+1,2);
    else if((i+1-s1-s2)%i3==0) A[i*n+i-1]=-pow(nx+1,2);//u zadnjem i uđe ovdje pa nema straha da će prekoračiti domenu
    else A[i*n+i+1]=A[i*n+i-1]=-pow(nx+1,2);

    if(i-s1-s2>=i3) A[i*n+i-i3]=-pow(ny+1,2);
    if(i-s2-s2<s3-i3) A[i*n+i+i3]=-pow(ny+1,2);

    if(i-s1-s2<i3 && k2!=0){
      A[j*n+i]=-pow(ny+1,2);
      A[i*n+j]=-pow(ny+1,2);
      j++;
    }
  }
  integer s=s1+s2+s3;
  for(i=s1+s1+s3;i<n;i++){
    A[i*n+i]=2*(pow(nx+1,2)+pow(ny+1,2));

    if((i-s)==k1 || (i-s)==0) A[i*n+i+1]=-pow(ny+1,2);
    else if((i-s+1)==k1 || (i-s+1)==n) A[i*n+i-1]=-pow(ny+1,2);
    else A[i*n+i+1]=A[i*n+i-1]=-pow(ny+1,2);
  }
/*
  printf("\nC iz A\n");
  for(i=s1+s1+s3;i<n;i++){
    for(j=s1+s2+s3;j<n;j++) printf("%.1lf   ", A[j*n+i]);
    printf("\n");
  }
  printf("\n\nC\n");
  print(C,g,g);

  printf("\n\nB1 iz A\n");
  for(i=0;i<s1;i++){
    for(j=0;j<s1;j++){
      printf("  %.1lf ", A[j*n+i]);
    }
    printf("\n");
  }
  printf("\n\nB1\n");
  print(B1,s1,s1);

  printf("\n\n");
  printf("\n\nE1 iz A\n");
  for(i=0;i<s1;i++){
    for(j=s1+s2+s3;j<n;j++){
      printf("  %.1lf ", A[j*n+i]);
    }
    printf("\n");
  }
  printf("\n\nE1\n");
  print(E1,s1,g);

  printf("\n\n");
  printf("\n\nF1 iz A\n");//provjera je li ispalo transopnirano
  for(j=s1+s2+s3;j<n;j++){
    for(i=0;i<s1;i++){
      printf("  %.1lf ", A[i*n+j]);
    }
    printf("\n");
  }

  printf("\nB2 iz A\n");
  for(i=s1;i<s1+s2;i++){
    for(j=s1;j<s1+s2;j++){
      printf("  %.1lf ", A[j*n+i]);
    }
    printf("\n");
  }
  printf("\n\nB2\n");
  print(B2,s2,s2);

  printf("\n\nE2 iz A\n");
  for(i=s1;i<s1+s2;i++){
    for(j=s1+s2+s3;j<n;j++){
      printf("  %.1lf ", A[j*n+i]);
    }
    printf("\n");
  }*//*
  printf("\n\nE2\n");
  print(E2,s2,g);

  printf("\n\nF2 iz A\n");//provjera je li ispalo transopnirano od E2 iz A
  for(j=s1+s2+s3;j<n;j++){
    for(i=s1;i<s1+s2;i++){
      printf("  %.1lf ", A[i*n+j]);
    }
    printf("\n");
  }


  printf("\nB3 iz A\n");
  for(i=s1+s2;i<s1+s2+s3;i++){
    for(j=s1+s2;j<s1+s2+s3;j++){
      printf("  %.1lf ", A[j*n+i]);
    }
    printf("\n");
  }
  printf("\n\nB3\n");
  print(B3,s3,s3);

  printf("\n\n");
  printf("\n\nE3 iz A\n");
  for(i=s1+s2;i<s1+s2+s3;i++){
    for(j=s1+s2+s3;j<n;j++){
      printf("  %.1lf ", A[j*n+i]);
    }
    printf("\n");
  }*//*
  printf("\n\nE3\n");
  print(E3,s3,g);

  printf("\n\n");
  printf("\n\nF3 iz A\n");//provjera je li ispalo transopnirano
  for(j=s1+s2+s3;j<n;j++){
    for(i=s1+s2;i<s1+s2+s3;i++){
      printf("  %.1lf ", A[i*n+j]);
    }
    printf("\n");
  }
  printf("\n\nF\n");


  print(F,g,s);*/
  char ch, strana;
  integer inc=1;
  int info;
  doublereal tol=1.0e-13;
  doublereal *b, *x1, *x2, *x3;

  //POČET S BROJANJEM VREMENA
  clock_t start_0, end_0;
  double cpu_time_used_0;
  start_0 = clock();

  f=calloc((s1+s1+s3),sizeof(doublereal));
  i=-1;
  printf("\nsustav E1\n\n");
  sustav (B1,E1,s1,g,tol,&info,strana,i);
  //print(E1,s1,g);

  //sustav Bf'=f
  x1=calloc(s1,sizeof(doublereal));
  f1=calloc(s1,sizeof(doublereal));//desna strana
  f1[s1-1]=10000.0;
  cg(B1,f1,x1,tol,s1,&info);
  for(i=0;i<s1;i++) f[i]=x1[i];

  printf("\n\n\nsustav E2\n\n");
  strana='L';
  sustav(B2,E2,s2,g,tol,&info,strana,gr2);
  //print(E2,s2,g);

  //sustav Bf'=f
  x2=calloc(s2,sizeof(doublereal));
  f2=calloc(s2,sizeof(doublereal));
  f2[s2-i2]=10000.0;
  cg(B2,f2,x2,tol,s2,&info);
  for(i=s1;i<s1+s2;i++) f[i]=x2[i-s1];

  printf("\n\n\nsustav E3\n\n");
  strana='D';
  sustav(B3,E3,s3,g,tol,&info,strana,gr3);
  //print(E3,s3,g);

  //sustav Bf'=f
  x3=calloc(s3,sizeof(doublereal));
  f3=calloc(s3,sizeof(doublereal));
  f3[i3-1]=10000.0;
  cg(B3,f3,x3,tol,s3,&info);//sustav Bf'=f
  for(i=s1+s2;i<s1+s2+s3;i++) f[i]=x3[i-s1-s2];

  // g' = g - Ff'
  doublereal alfa, beta;

  char N='N';
  alfa=-1.0; beta=1; i=1;
  dgemv_(&N,&g,&s,&alfa,F,&g,f,&i,&beta,gvec,&i);
        // g' je spremljen u g


  // S = C - FE'
  E=calloc(s*g,sizeof(doublereal));
  for(i=0;i<s1;i++){
    for(j=0;j<g;j++) E[j*s+i]=E1[j*s1+i];
  }
  for(i=s1;i<s1+s2;i++){
    for(j=0;j<g;j++) E[j*s+i]=E2[j*s2+i-s1];
  }
  for(i=s1+s2;i<s1+s2+s3;i++){
    for(j=0;j<g;j++) E[j*s+i]=E3[j*s3+i-s1-s2];
  }
printf("\nkonačni E\n");print(E,s,g);

  dgemm_(&N,&N,&g,&g,&s,&alfa,F,&g,E,&s,&beta,C,&g);
        // S je spremljen u C
  //printf("\nS:\n");
  //print(C,g,g);

  // cg na Sy=g'
  doublereal *xp;
  xp=calloc(g,sizeof(doublereal));
  cg(C,gvec,xp,tol,g,&info);

// DO OVDJE BI TREBALO BIT DOBRO JER JE S SIMETRIČAN
  // x = f' - E'y
  doublereal *z; i=1;
  z=calloc(s,sizeof(doublereal));

  dgemv_(&N,&s,&g,&alfa,E,&s,xp,&i,&beta,f,&i);

  end_0 = clock();

  cpu_time_used_0 = ((double) (end_0 - start_0)) / CLOCKS_PER_SEC;
  printf("\nvrijeme algoritma:%lf   \n", cpu_time_used_0);
  printf("\nn0=%d\n",s+g);
  printf("\naproksimacija algoritmom\n");
  for(i=0;i<s;i++) printf("%.4lf  ", f[i]);
  for(i=0;i<g;i++) printf("%.4lf  ", xp[i]);




  // na cijeloj matrici A
  clock_t start, end;
  double cpu_time_used;

  start = clock();

  doublereal *b_, *x0;
  b_=calloc(n,sizeof(doublereal));
  x0=calloc(n,sizeof(doublereal));
  b_[s1-1]=10000.0;
  b_[s1+s2-i2]=10000.0;
  b_[s1+s2-1+i1]=10000.0;
  cg(A, b_, x0, tol, n, &info);

  end = clock();

  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("\n\n\nvrijeme na cijeloj matrici:%lf   \n", cpu_time_used);
  printf("\nn=%d\n", n);
  printf("\n\naproksimacija na cijeloj matrici:\n");
  for(i=0;i<n;i++) printf("%.4lf  ", x0[i]);

// probat ne radit nul-stupce u sustavu BE'=E
// probat doc do podmatrice

/*
  printf("\nF\n");
  print(F,g,(s1+s2+s3));
  printf("\nE1\n");
  print(E1,s1,g);
  printf("\nE2\n");
  print(E2,s2,g);
  printf("\nE3\n");
  print(E3,s3,g);*/
}
