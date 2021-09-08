#include <stdio.h>
#include <stdlib.h> // for rand()
#include <time.h> // for clock()
#include <math.h>
#include "mat.h"

/* 
testing method of solving linear system via QR decomposition and backward substitution
(see also solver.c and qrGS1.c, qrGS2.c)
Start with matrix A, vector x. 
Multiply Ax, call that b.
Solve Av=b.
v should well approximate x.

*/

int main()
{
  // ints & sizes
  int n = 10;
  int size = n*n;
  int i, j;

  //matrices 
  double A[size];
  double A2[size];
  double Q[size];
  double Qt[size];
  double R[size];
 
  for (i=0; i<size; i++)
    {
      A[i] = (double)rand()/RAND_MAX*10.0-5.0;
      Q[i] = Qt[i] = A2[i] = R[i] = 0.0;
    }

  //vectors 
  double x[n], b[n], Ax[n];
  for (i=0; i<n; i++)
    {
      b[i] = (double)rand()/RAND_MAX*10.0-5.0;
      x[i] = Ax[i] = 0.0;
    }

  
  printf("A:\n");
  printmat(A,n,n);
  printf("b:\n");
  printmat(b,n,1);
  
  // QR decomposition
  qrGS2(A, n, n, &Q, &R);
  printf("Q:\n");
  printmat(Q,n,n);
  printf("R:\n");
  printmat(R,n,n);

  // check
  printf("Check A=QR:\n");
  matmat(Q, n, n, R, n, n, &A2, n, n);
  printmat(A2, n, n);
  scaleVec(&A2, -1.0, size);  
  addVec(&A2, A, size);
  printf("Check A-QR=0:\n");
  printmat(A2,n,n);
  printf("\n");

  // solve 
  // NOTE: qrSolve computes Q and R separately from the rest of this exercise
  qrSolve(A, n, &x, b);

  // check
  printf("Resultant x:\n");
  printmat(x, n, 1);

  printf("Ax:\n");
  matvec(A, n, n, x, n, &Ax);
  printmat(Ax, n, 1);

  printf("Ax - b:\n");
  scaleVec(&b, -1.0, n);
  addVec(&Ax, b, n);
  printmat(Ax, n, 1);

  printf("printMax on Ax-b:\n");
  printMax(Ax, n);
  printf("\n");

  /*
  printf("printMax on Ax-Av:\n");
  matvec(A, n, n, v, n, &u);
  scaleVec(&u, -1.0, n);
  addVec(&u, b, n);
  printMax(u, n);
  */
}
