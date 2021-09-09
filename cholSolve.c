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
  double Rt[size];
  double y[n];
  int i;
  for (i=0; i<size; i++)
    {
      R[i] = Rt[i] = 0;
    }

  // cholesky decomp
  cholFac(A, m, &R);
  
  // solve R'y = b via back substitution
  // opportunity for speedup here, see note (2i)
  transpose1(R,m,m,Rt);
  backsub(R, m, b, &y);
  
  // solve Rx = y via forward substitution 
  forwardsub(R, m, y, &x);
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
        printf("backsub:\nsingular R: zero entry on diagonal\n");
    }
}



// given n \times n lower triangular matrix A and RHS b, solve Ax=b for x
void forwardsub(double A[], int n, double b[], double (*x)[])

{
  int i, j;
  double xaccum, rij;
  xaccum = rij = 0;

  for (i=0; i < n; i++)
    {
      xaccum = b[i];
      for (j=0; j<i; j++)
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
        printf("forwardsub:\nsingular input matrix: zero entry on diagonal\n");
    }
}


// compute lower triangular R via Cholesky Factorization
// Actually R will have garbage in the upper triangular part, 
// to be ignored by solvers downstream. see note 1i
void cholFac(double A[], int m, double (*R)[]) 
{
  int k, j;
  double Rkk, scale, tmp;  
  Rkk = scale = 0.0;
  double Rj[m], Rkm[m], Rkj[m];
  //  double vi[m], vj[m], qi[m], aj[m], qtemp[m];
  for (k=0;k<m;k++)
    {
      Rjm[k] = Rkm[k] = Rkj[k] = 0.0;
    }

  // set R = A
  for (k=0; k<m*m; k++)
    {
      (*R)[k] = A[k];
    }

  // go
  for (k=0; k<m; k++)
    {
      for (j=k+1; j<m; j++)
        {
          // R_j,j:m = R_j,j:m - R_k,j:m R'_kj/R_kk
          // see note (1i)
          extractRow(*R, m, m, j, &Rj);
          extractRow(*R, m, m, k, &Rkj);
          scale = -extractEntry(*R, m, m, k, j)/extractEntry(*R, m, m, k, k);
          scaleVec(&Rkj, scale, m);
          addVec(&Rj, Rkj, m);
          updateRow(R,m,m,j,Rj);
        }
      scale = 1/sqrt(extractEntry(*R,m,m,k,k));
      extractColumn(*R,m,m,k,Rkm);
      for (i=k; i<m; i++) // this column is not exempt from concern in (1i)
        {
          tmp = Rkm[i]*scale;
          Rkm[i] = tmp;
        }
      updateCol(R,m,m,k,Rkm);
    }
}



  /* 
all these extract/update-s might be expensive...
How to implement if performance is the concern?

Pass in arrays of pointers to arrays rather than arrays of doubles. 
Each sub-array will then be a column. 
That way we can more easily access each column. 
  */


/* notes
(1i) The text calls for using vectors from principal submatrices, 
but we can get away with using the full vectors, 
and ignoring the garbage in the entries outside of the principal submatrix, 
since we assume the matrix is upper triangular. 
If I was really diligent I would use a scheme that uses arrays that fill 
only half the matrix in question to avoid issues like this. 

(2i) Rather than transpose R here, 
we could write a "backwards back solver" that solves a lower triangular system 
with an upper triangular matrix
*/


