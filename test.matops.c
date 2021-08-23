#include <stdio.h>
#include <stdlib.h> // for rand()
#include <time.h> // for clock()
#include <math.h>
#include "mat.h"


// testing
int main()
{
  /*
  int m=5;
  int n=3;
  int k=2;

  double x[] = {1,2,3};
  double A[] = {-1, 2, 1, 2, -1, 2, 0, 2, -1, 1, 2, 3, 1, -1, 1};
  double B[] = {-1, 0, 2, 2, 1, 0};
  double C[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  double b[] = {0, 0, 0, 0, 0};
  int i, j;


  printf("test matmat:\n");
  matmat(A, m, n, B, n, k, &C, m, k);

  printf("A:\n");
  printmat(A,m,n);
  for (i=0; i<(m*n); i++) printf("%lf ", A[i]);
  printf("\n");

  printf("B:\n");
  printmat(B,n,k);

  printf("C:\n");
  printmat(C,m,k);
  printf("\n");

  printf("test matvec:\n");
  matvec(A, m, n, x, n, &b);
  for (i=0; i<m; i++) printf("%lf ", b[i]);
  printf("\n\n");

  printf("test dotprod:\n");
  double dd;
  double y[] = {1.3,5.77,-12.44};
  dd = dotprod(x, y, n);
  printf("%lf\n", dd);
  printf("\n");

  printf("test eucnorm:\n");
  dd = eucnorm(y, n);
  printf("%lf\n", dd);
  printf("\n");

  printf("test eucmetric:\n");
  printf("x: ");
  for (i=0; i<n; i++) printf("%lf ", x[i]);
  printf("\n");
  printf("y: ");
  for (i=0; i<n; i++) printf("%lf ", y[i]);
  printf("\n");
  dd = eucmetric(x, y, n);
  printf("||x-y||: %lf\n", dd);
  printf("\n");
  
  printf("test scalevec:\n");
  scaleVec(&y, 6.7, n);
  for (i=0; i<n; i++) printf("%lf ", y[i]);
  printf("\n");

  printf("\n");

  double A2[m*n];
  for (i=0; i<m*n; i++) A2[i] = A[i];
  double u[] = {9,9,9,9,9};

  printf("test updateCol:\n");
  printf("A: \n");
  printmat(A2,m,n);
  updateCol(&A2, m, n, 1, u);
  printf("A2: \n");
  printmat(A2,m,n);
  printf("\n");
  */

  /*
  double startTime, endTime, timeElapsed;
  startTime = (float)clock()/CLOCKS_PER_SEC;
  for (i=0; i<300; i++)
    transpose1(A,m,n,&At);
  endTime = (float)clock()/CLOCKS_PER_SEC;
  timeElapsed = endTime - startTime;
  printf("elapsed: %lf\n", timeElapsed);



  extractColumn(A, m, n, 3, &v);
  printf("A:\n");
  printmat(A,m,n);
  printf("\n");
  printf("v:\n");
  printmat(v,m,1);

  updateEntry(666, &A, m, n, 6, 2);
  printmat(A,m,n);


  printf("u:\n");
  printmat(u,m,1);
  printf("v:\n");
  printmat(v,m,1);

  printf("u+v:\n");
  addVec(&u, v, m);
  printmat(u,m,1);

  int m = 518;
  int n = m;
  int size = m*n;
  int squareSize = n*n;

  double A[size];
  double A2[size];
  double Q[size];
  double R[size];
  int i, j;
  
  for (i=0; i<size; i++)
    {
      A[i] = (double)rand()/RAND_MAX*20-10;
      Q[i] = 0.0;
      A2[i] = 0.0;
      if (i<squareSize)
        R[i] = 0;
    }

  double v[m], u[m];
  for (i=0; i<m; i++)
    {
      v[i] = 0.0;
      u[i] = 2.0;
    }

  
  printf("A:\n");
  //printmat(A,m,n);
  qrGS1(A, m, n, &Q, &R);
  printf("Q:\n");
  //printmat(Q,m,n);
  printf("R:\n");
  //printmat(R,n,n);

  printf("QR:\n");
  matmat(Q, m, n, R, n, n, &A2, m, n);
  //printmat(A2,m,n);

  //check
  scaleVec(&A2, -1.0, size);  
  addVec(&A2, A, size);
  //printf("A-QR:\n");
  //printmat(A2,m,n);
  printf("printMax on A-QR:\n");
  printMax(A2, m*n);



  int m = 3;
  int n = 3;
  int size = m*n;
  int squareSize = n*n;

  double A[] = {1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0};
  double A2[size];
  double Q[size];
  double R[size];
  int i, j;

  for (i=0; i<size; i++)
    {
      Q[i] = 0.0;
      A2[i] = 0.0;
      if (i<squareSize)
        R[i] = 0;
    }

  printf("A:\n");
  printmat(A,m,n);
  qrGS1(A, m, n, &Q, &R);
  printf("Q:\n");
  printmat(Q,m,n);
  printf("R:\n");
  printmat(R,n,n);

  printf("QR:\n");
  matmat(Q, m, n, R, n, n, &A2, m, n);
  printmat(A2,m,n);

  //check
  scaleVec(&A2, -1.0, size);  
  addVec(&A2, A, size);
  printf("A-QR:\n");
  printmat(A2,m,n);
*/


  int m = 5;
  int n = m;
  int size = m*n;
  int squareSize = n*n;

  double A[size];
  double A2[size];
  double Q[size];
  double R[size];
  int i, j;
  
  for (i=0; i<size; i++)
    {
      A[i] = (double)rand()/RAND_MAX*20-10;
      Q[i] = 0.0;
      A2[i] = 0.0;
      if (i<squareSize)
        R[i] = 0;
    }

  double v[m], u[m];
  for (i=0; i<m; i++)
    {
      v[i] = 0.0;
      u[i] = 2.0;
    }

  
  printf("A:\n");
  //printmat(A,m,n);
  qrGS2(A, m, n, &Q, &R);
  printf("Q:\n");
  //printmat(Q,m,n);
  printf("R:\n");
  //printmat(R,n,n);

  printf("QR:\n");
  matmat(Q, m, n, R, n, n, &A2, m, n);
  //printmat(A2,m,n);

  //check
  scaleVec(&A2, -1.0, size);  
  addVec(&A2, A, size);
  //printf("A-QR:\n");
  //printmat(A2,m,n);
  printf("printMax on A-QR:\n");
  printMax(A2, m*n);



  /*
  int m = 5;
  int n = 3;
  int size = m*n;
  int squareSize = n*n;

  double A[size];
  double A2[size];
  double Q[size];
  double R[size];
  int i, j;

  for (i=0; i<size; i++)
    {
      A[i] = (double)rand()/RAND_MAX*20-10;
      Q[i] = 0.0;
      A2[i] = 0.0;
      if (i<squareSize)
        R[i] = 0;
    }

  printf("A:\n");
  printmat(A,m,n);
  qrGS2(A, m, n, &Q, &R);
  printf("Q:\n");
  printmat(Q,m,n);
  printf("R:\n");
  printmat(R,n,n);

  printf("QR:\n");
  matmat(Q, m, n, R, n, n, &A2, m, n);
  printmat(A2,m,n);

  //check
  scaleVec(&A2, -1.0, size);  
  addVec(&A2, A, size);
  printf("A-QR:\n");
  printmat(A2,m,n);
  */
}
