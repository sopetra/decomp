#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "math.h"
#include "Poisson_L.c"

int main(){
  integer nx, ny;
  doublereal hx, hy;

  nx=5; ny=11;

//  Poisson_L(nx,ny);
  nx=5; ny=8;

  Poisson_L(nx,ny);

  return 0;
}
