#include <stdio.h>
#include "mat.h"
#include <math.h>

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
  double R[size];
  double Rt[size];
  double y[n];
  int i;
  for (i=0; i<size; i++)
    {
      R[i] = Rt[i] = 0;
    }

  // cholesky decomp
  cholFac(A, n, &R);
  printf("R after cholfac:\n");
  printmat(R,n,n);

  // solve R'y = b via back substitution
  // opportunity for speedup here, see note (2i)
  transpose1(R, n, n, &Rt);
  forwardsub(Rt, n, b, &y);
  printf("y after backsub:\n");
  printmat(y,n,1);
  
  // solve Rx = y via forward substitution 
  backsub(R, n, y, x);
  printf("x after forwardsub:\n");
  printmat(*x,n,1);
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
  int i, k, j;
  double Rji, Rki, Rkk, scale, tmp;  
  Rkk = scale = 0.0;
  double Rj[m], Rkm[m], Rkj[m];
  //  double vi[m], vj[m], qi[m], aj[m], qtemp[m];
  for (k=0;k<m;k++)
    {
      Rj[k] = Rkm[k] = Rkj[k] = 0.0;
    }

  // set R = A
  for (k=0; k<m*m; k++)
    {
      (*R)[k] = A[k];
    }

  printf("matrix R before cholfac loop(s):\n");
  printmat(*R,m,m);
  // go
  for (k=0; k<m; k++)
    {
      for (j=k+1; j<m; j++)
        {
          // R_j,j:m = R_j,j:m - R_k,j:m R'_kj/R_kk
          // see note (1i)
          scale = -extractEntry(*R, m, m, k, j)/extractEntry(*R, m, m, k, k);
          for (i=j; i<m; i++)
            {
              Rji = extractEntry(*R,m,m,j,i);
              Rki = extractEntry(*R,m,m,k,i);
              Rji = Rji + scale*Rki;
              updateEntry(Rji, R, m, m, j, i);
            }
        }
      Rkk = extractEntry(*R,m,m,k,k);
      printf("Rkk at iter k=%d: %.15lf \n",k, Rkk);
      scale = 1.0/sqrt(Rkk);
      printf("scale at iter k=%d: %.15lf \n",k, scale);

      // R_k,k:m = R_k,k:m / sqrt(R_kk)     
      for (i=k; i<m; i++) 
        {
          Rki = extractEntry(*R,m,m,k,i);
          Rki = Rki*scale;
          updateEntry(Rki,R,m,m,k,i);
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


/* notes
(1i) The text uses matlab style slicing which I don't really have so I use an 
extra loop in i. 
The first version actually used extract/updateRow which was probably faster, 
better etc but does weird things to the subdiagonal entries that annoyed me. 

Also, the algorithm in text assumes A is upper triangular, 
I don't zero out sub-diagonal entries. 
Downstream uses of R may have to remember that it contains garbage.

If I was really diligent I would use a scheme that uses arrays that fill 
only half the matrix in question to avoid issues like this. 


(2i) Rather than transpose R here, 
we could write a "backwards back solver" that solves a lower triangular system 
with an upper triangular matrix
*/


