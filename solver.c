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


void qrSolve(double A[], int n, double (*x)[], double b[])
{
  int size = n*n;
  double Q[size];
  double Qt[size];
  double R[size];

  // QR
  qrGS2(A, n, n, &Q, &R);

  // Q'b = y
  transpose1(Q, n, n, &Qt);

  // solve via back substitution
  backsub(R, n, y, x);
}


// given n \times n upper triangular matrix A and RHS b, solve Ax=b for x
void backsub(double A[], int n, double b[], double (*x)[])
{
  int i, j;
  double xaccum, rij;
  xaccum = rij = 0;

  for (i=n-1; i >= 0; i--)
    {
      xaccum = b[i];
      for (j=n-1; j > i; j--)
        {
          rij = extractEntry(A, n, n, i, j);
          xaccum -= rij*(*x)[j];
        }
      if ((rij = extractEntry(A, n, n, i, i)) != 0.0) // replace w ||<tol?
        {
          xaccum = xaccum/rij;
          (*x)[i] = xaccum;
        }
      else 
        printf("singular R: zero entry on diagonal\n");
    }
}

