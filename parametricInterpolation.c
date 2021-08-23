#include <stdio.h>
#include <stdlib.h> // for rand()
#include <time.h> // for clock()
#include <math.h>
#include "mat.h"

/* 
Using solver for Vandermonde matrices in parametric polynomial interpolation. 
Motivated by the question 
"what if data we would like to interpolate is not unique in the dependent variable?"
i.e. some (x_i, y_i), (x_j, y_j) have x_i = x_j
Solution: give parametric function f(t) = (x(t), y(t)) where x(t) and y(t) interpolate 
on some new array of t_i-s. 
*/

int main()
{
  // ints & sizes
  int n = 4;
  int size = n*n;
  int i, j;


  double Vt[size]; // vandermonde matrix 
  double T[n]; // array of time points (somewhat arbitrary, user need not choose this)
  double X[n]; // coefficients of x-dim polynomial 
  double Y[n];
 
  // setup T, Vt
  double u[n], v[n], x[n], b[n], y[n];
  for (i=0; i<n; i++)
    {
      x[i] = (double)rand()/RAND_MAX*10.0-5.0;
      v[i] = b[i] = y[i] = u[i] = 0.0;
    }

  for (i=0; i<size; i++)
    {
      A[i] = (double)rand()/RAND_MAX*10.0-5.0;
      Q[i] = Qt[i] = A2[i] = R[i] = 0.0;
    }

  //vectors 
  double u[n], v[n], x[n], b[n], y[n];
  for (i=0; i<n; i++)
    {
      x[i] = (double)rand()/RAND_MAX*10.0-5.0;
      v[i] = b[i] = y[i] = u[i] = 0.0;
    }

  
  printf("A:\n");
  printmat(A,n,n);
  printf("x:\n");
  printmat(x,n,1);
  
  matvec(A, n, n, x, n, &b);
  printf("Ax:\n");
  printmat(b,n,1);

  // QR decomposition
  qrGS1(A, n, n, &Q, &R);
  printf("Q:\n");
  printmat(Q,n,n);
  printf("R:\n");
  printmat(R,n,n);

  //check
  printf("QR:\n");
  matmat(Q, n, n, R, n, n, &A2, n, n);
  printmat(A2, n, n);
  scaleVec(&A2, -1.0, size);  
  addVec(&A2, A, size);
  printf("A-QR:\n");
  printmat(A2,n,n);
  printf("printMax on A-QR:\n");
  printMax(A2, size);
  printf("\n");

  // Q'b = y
  transpose1(Q, n, n, &Qt);
  matvec(Qt, n, n, b, n, &y);

  // backsolve
  backsub(R, n, y, &v);

  // check
  printf("solution:\n");
  printmat(v, n, 1);

  scaleVec(&x, -1.0, n);
  addVec(&x, v, n);
  printf("printMax on v-x:\n");
  printMax(x, n);
  printf("\n");

  /*
  printf("printMax on Ax-Av:\n");
  matvec(A, n, n, v, n, &u);
  scaleVec(&u, -1.0, n);
  addVec(&u, b, n);
  printMax(u, n);
  */
}
