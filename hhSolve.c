#include <stdio.h>
#include "mat.h"

/* 
solve linear system via QR factorization 

UFD: how to read data 
for now taking a vector holding matrix in row-major order as function argument

Diagonalizes m \times n matrix A with m >= n

1. orthogonalize via gram-schmidt (A=QR)
2. apply Q -> y=Qb

Ax = b
QRx = b
Rx = Qb =: y

3. solve Rx = y via back substitution

*/


void hhSolve(double A[], int n, double (*x)[], double b[])
{
  int size = n*n;
  double Q[size];
  double Qt[size];
  double R[size];
  double y[n];
  int i;

  // QR
  qrHH(A, n, n, &Q, &R);

  // Q'b = y
  transpose1(Q,n,n,&Qt);
  matvec(Qt, n, n, b, n, &y);

  // solve via back substitution
  backsub(R, n, y, x);
}


  /* 
all these extract/update-s might be expensive...
How to implement if performance is the concern?

Pass in arrays of pointers to arrays rather than arrays of doubles. 
Each sub-array will then be a column. 
That way we can more easily access each column. 
  */
