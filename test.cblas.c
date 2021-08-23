#include <stdio.h>
#include <math.h>
#include "cblas.h"
//#include <mkl.h>

int main()
{
  double x[] = {1,0,1};
  double A[] = {-1, 2, 0, 2, -1, 2, 0, 2, -1};
  double b[3];
  int i;
  
  cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, A, 3, x, 1, 1.0, b, 1);

  for (i=0;i<3;i++) printf("%lf ", b[i]);
  printf("\n");

}
