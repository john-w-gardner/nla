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



  /* 
all these extract/update-s might be expensive...
How to implement if performance is the concern?

Pass in arrays of pointers to arrays rather than arrays of doubles. 
Each sub-array will then be a column. 
That way we can more easily access each column. 
  */
