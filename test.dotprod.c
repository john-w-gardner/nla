#include <stdio.h>
#include <stdlib.h> // for rand()
#include <time.h> // for clock()
#include <math.h>
#include "mat.h"

// dotprod() seems to be acting up...

// testing
int main()
{
  int n = 3;
  double u[3] = {1, 1, 1};
  double v[3] = {1, 2, 1};
  double dd = 0.0;
  int i, j;

  printf("u:\n");
  printmat(u,n,1);
  printf("v:\n");
  printmat(v,n,1);

  dd = dotprod(u, v, n);
  printf("dd: %lf\n", dd);
}
