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

/* "classical" Gram-Schmidt orthogonalization for reduced QR
 see TB p.51
 A and Q are m \times n with m >= n,
 R is n \times n
*/


/* WHYDOESTHISWORK:
there is a bug at rij = dotprod() where aj should be qi.
this sets rij = 0. 
the result is that R is diagonal, not upper triangular.
At any rate, the product QR still returns the same A we start with, 
so it appears to be a valid decomposition. 
Why??
 */

void qrGS1(double A[], int m, int n, double (*Q)[], double (*R)[]) 
{
  int i, j;
  double v[m], qi[m], aj[m];
  for (i=0;i<m;i++)
    {
      v[i] = qi[i] = aj[i] = 0;
    }
  double rij = 0.0;
  
  for (j=0; j<n; j++)
    {
      extractColumn(A, m, n, j, &v);
      for (i=0; i<j; i++)
        {
          // r_{ij} = q_i'a_j
          extractColumn((*Q), m, n, i, &qi);
          printf("qi and aj and v:\n");
          printmat(qi, m, 1);
          printf("\n");
          printmat(aj, m, 1);
          printf("\n");
          printmat(v, m, 1);
          printf("\n");
          rij = dotprod(qi, aj, m); // bug
          printf("r_ij at i=%d, j=%d: %lf\n", i, j, rij);
          updateEntry(rij, R, n, n, i, j);
          
          // v_j = v_j - r_{ij}q_i
          scaleVec(&qi, -rij, m);
          addVec(&v, qi, m);
        }
      // r_{jj} = ||v_j||
      rij = eucnorm(v, m);
      printf("r_ij at i=%d, j=%d: %lf\n", i, j, rij);
      updateEntry(rij, R, n, n, j, j);
      // q_j = v_j/r_jj
      if (rij > 0.0)
        scaleVec(&v, 1/rij, m);
      else
        printf("division by zero at rjj update\n");
      updateCol(Q, m, n, j, v);
    }

}


