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
  double y[n];
  int i;
  for (i=0; i<size; i++)
    {
      Q[i] = R[i] = 0;
    }

  // QR
  qrGS2(A, n, n, &Q, &R);

  // Q'b = y
  transpose1(Q, n, n, &Qt);
  matvec(Qt, n, n, b, n, &y);

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

#include <stdio.h>
#include "mat.h"

/* "modified" Gram-Schmidt orthogonalization for reduced QR
 see TB p.58
 A and Q are m \times n with m >= n,
 R is n \times n
*/

/* overwriting v_i-s into A, since A is temporary anyway
   i.e. skip first two lines in Algorithm 8.1
   
*/   
void qrGS2(double A[], int m, int n, double (*Q)[], double (*R)[]) 
{
  int i, j;
  double rij = 0.0;  
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


// test Ax ~ b by returning norm |Ax-b|
double testSolve(double A[], int n, double x[], double b[])
{
  double Ax[n]; // product Ax
  double nAx; // norm of Ax
  int i;

  matvec(A, n, n, x, n, &Ax);
  scaleVec(&Ax, -1.0, n);
  for (i=0; i<n; i++) 
    {
      printf("i: %d -Ax: %.16lf b: %.16lf\n", i, Ax[i], b[i]);
    }
  addVec(&Ax, b, n);
  nAx = eucnorm(Ax, n);
  printf("|Ax-b|: %.16lf\n", nAx);
  return nAx;
}



  /* 
all these extract/update-s might be expensive...
How to implement if performance is the concern?

Pass in arrays of pointers to arrays rather than arrays of doubles. 
Each sub-array will then be a column. 
That way we can more easily access each column. 
  */
