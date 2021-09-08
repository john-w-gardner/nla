#include <stdio.h>
#include "mat.h"


/* Solve linear system via Cholesky Factorization
   see TB p. 175-176

   Input matrix A must be symmetric positive definite. 
   Procedure goes as follows (solving Ax=b for x):

   (1) Decompose A=R'R via Cholesky

   (2) Solve R'y = b for y
   (since R'Rx = Ax = b => Rx = (R')^-1 b, need to solve instead of just multiplying like QR via Gram-Schmidt )

   (3) Solve Rx = y for x

   UFD: make R an array with (n^2 - n)/2 entries since upper triangular
*/

/* TODO:
   write cholesky decomp
   write forward substitution solver
   write cholSolve
   move backsub to neutral location
*/

void cholSolve(double A[], int n, double (*x)[], double b[])
{
  int size = n*n;
  double Q[size];
  double Qt[size];
  double R[size];
  double y[n];
  int i;
  for (i=0; i<size; i++)
    {
      R[i] = 0;
    }

  // cholesky decomp

  // solve R'y = b via back substitution

  // solve Rx = y via forward substitution 

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



// given n \times n lower triangular matrix A and RHS b, solve Ax=b for x
void forwardSub(double A[], int n, double b[], double (*x)[])
/*
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
*/
}


// compute lower triangular R via Cholesky Factorization
void cholFac(double A[], int m, int n, double (*R)[]) 
{
  int i, j;
  double Rkk = 0.0;  
  double vi[m], vj[m], qi[m], aj[m], qtemp[m];
  for (i=0;i<m;i++)
    {
      vi[i] = qi[i] = vj[i] = qtemp[i] = 0;
    }

  // set Q = A
  for (i=0; i<m*n; i++)
    {
      (*Q)[i] = A[i];
    }

  for (i=0; i<n; i++)
    {
      // r_ii = ||v_i||
      extractColumn(*Q, m, n, i, &vi);
      rij = eucnorm(vi, m);
      updateEntry(rij, R, n, n, i, i);

      // q_i = v_i/r_ii
      updateCol(&qi, m, 1, 0, vi);
      scaleVec(&qi, 1/rij, m);
      updateCol(Q, m, n, i, qi);

      for (j=i+1; j<n; j++)
        {
          // rij = q_i'v_j
          extractColumn(*Q, m, n, j, &vj);
          rij = dotprod(qi, vj, m);
          updateEntry(rij, R, n, n, i, j);

          // v_j = v_j - r_{ij}q_i
          updateCol(&qtemp, m, 1, 0, qi); // need to keep q_i as is 
          scaleVec(&qtemp, -1.0*rij, m);
          addVec(&vj, qtemp, m);
          updateCol(Q, m, n, j, vj);
        }
    }

}



  /* 
all these extract/update-s might be expensive...
How to implement if performance is the concern?

Pass in arrays of pointers to arrays rather than arrays of doubles. 
Each sub-array will then be a column. 
That way we can more easily access each column. 
  */
