#include <stdio.h>
#include "math.h"
#include "mat.h"

/* QR factorization via Householder reflectors
 see TB p.73
 this is a full QR factorization i.e.
 A is m \times n with m >= n,
 Q is m \times m 
 R is m \times n
*/

/* 
notes:
"easier but more expensive version".. 
attempting matlab style submatrix extraction at inner loop to hide innermost loop(s)
*/   

void qrHH(double A[], int m, int n, double (*Q)[], double (*R)[]) 
{
  // counters, accumulators
  int i, j, k, l;
  double r, xn, vn;
  r = xn = vn = 0.0;
  double tolvn = 0.00000000001; // potential division by zero..


  // store v_k vectors to construct Q after main loop
  double V[m*n];

  // set R = A
  insertSubmat(R,m,n,0,m,0,n,A);
  // set  Q = I
  for (i=0; i<m; i++)
    {
      for (j=0; j<m; j++)
        {
          if (i == j) 
            { 
              //printf("updating Q with 1 on diagonal.. i=%d j=%d\n",i,j);
              updateEntry(1.0, Q, m, m, i, j); 
              //printmat(*Q,m,m);
            }
          else 
            { 
              //printf("updating Q off diagonal.. i=%d j=%d\n",i,j);
              updateEntry(0.0, Q, m, m, i, j); 
              //printmat(*Q,m,m);
            }
          //printf("A:\n");
          //printmat(A,m,n);
        }
    }
  //printf("printmat Q then R then A in qrHH:\n");
  //printmat(*Q,m,m);
  //printf("R:\n");
  //printmat(*R,m,n);
  //printf("A:\n");
  //printmat(A,m,n);
  
  // construct R
  for (k=0; k<n; k++)
    {
      // inner arrays need to be able to change size
      double v[m-k];
      double Atmp[(m-k)*(n-k)];
      double Vtmp[(m-k)*(n-k)];

      // x = A_{k:m, k}
      //printf("inside qrHH construction of R.. k=%d\n",k);
      submat(*R, m,n,k,m,k,k+1,&v);
      /*printf("A:\n");
      //printmat(A,m,n);
      printf("v:\n");
      //printmat(v,m-k,1);
      printf("v[0]: %f\n",v[0]);
      */

      // v_k = sign(x_1) \norm{x} e_1 + x
      xn = (double)sign(v[0])*eucnorm(v, m-k);
      v[0] = v[0]+xn;

      // v_k = v_k/\norm{v_k}
      vn = eucnorm(v,m-k);
      if (vn < tolvn) {printf("WARNING: very small ||v_n||: %f\n", vn);}
      scaleVec(&v, 1.0/vn, m-k);
      //printf("v_k/norm(v_k):\n");
      //printmat(v,m-k,1);

      // A_{k:m, k:n} = A_{k:m, k:n} - 2v_k(v'_k A_{k:m, k:n})
      submat(*R,m,n,k,m,k,n,&Atmp);
      //printf("Atmp:\n");
      //printmat(Atmp,m-k,n-k);
      // v'A
      double vtmp[n-k];
      //printf("first matmat:\n");
      matmat(v,1,m-k,Atmp,m-k,n-k,&vtmp,1,n-k);
      //printf("vtmp after first matmat:\n");
      //printmat(vtmp,1,n-k);
      // vv'A
      
      matmat(v,m-k,1,vtmp,1,n-k,&Vtmp,m-k,n-k);
      //printf("Vtmp after second matmat:\n");
      //printmat(Vtmp,m-k,n-k);
      // A - 2vv'A
      scaleVec(&Vtmp,-2.0,(m-k)*(n-k));
      addVec(&Atmp, Vtmp,(m-k)*(n-k));

      // write into R
      insertSubmat(R,m,n,k,m,k,n,Atmp);
      // save V
      insertSubmat(&V,m,n,k,m,k,k+1,v);
      //printf("end of %d-th loop.. R:\n",k);
      //printmat(*R,m,n);
    }
  //printf("R completed, here's V:\n");
  //printmat(V,m,n);
  // construct Q
  for (k=n-1; k>=0; k--)
    {
      double v[m-k];
      double vtmp[0];
      
      for (i=0; i<m; i++)
        {
          submat(V,m,n,k,m,k,k+1,&v);
          //printf("v in construct Q k=%d, i=%d\n",k,i);
          //printmat(v,m-k,1);

          // e_i_k:m = e_i_k:m - 2v_k(v'_k e_i_k:m)
          double e[m-k];
          submat(*Q,m,m,k,m,i,i+1,&e);
          //printf("e:\n");
          //printmat(e,m-k,1);
          // v'e
          matmat(v,1,m-k,e,m-k,1,&vtmp,1,1);
          vn = vtmp[0]; //reusing vn its just a float
          // -2vv'e
          scaleVec(&v,-2.0*vn,m-k);
          //printf("v after all calculations:\n");
          //printmat(v,m-k,1);

          // e - 2vv'e
          addVec(&e, v, m-k);
          //printf("e after all calculations:\n");
          //printmat(e,m-k,1);

          // write into Q
          insertSubmat(Q,m,m,k,m,i,i+1,e);
        }
    }
  //printf("done.. Q:\n");
  //printmat(*Q,m,m);
  
}



/*
Q CONSTRUCTION: 
could extract submatrix of I rather than individual columns e_k

// e_k:m = e_k:m - 2v_k(v'_k e_k:m)
double I[(m-k)*(n-k)];
double v[m-k];
double vtmp[m-k];
double Vtmp[(m-k)*(n-k)];
submat(*Q,m,n,k,m,k,n,&I);
submat(V,m,n,k,m,k,k,&v);

// v'e
matmat(v,1,m-k,e,m-k,1,vtmp,1,1);
// -2vv'e
scaleVec(v,-2.0*vtmp,m-k);
// A - 2vv'A
addVec(&e, v, m-k);

// write into Q
insertSubmat(Q,m,n,m-k,m,n-k,n-k,e);
*/


/*
      for (i=0; i<n; i++)
        {
          e[0] = 1.0;
          for (i=1;i<m-k;i++) {e[i] = 0.0;}


      // vv' will be used twice
      double vvp[(m-k)*(n-k)];
      matmat(v,m-k,1,v,1,m-k,&vvp,m-k,n-k);

      
*/
