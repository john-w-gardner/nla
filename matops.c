#include <stdio.h>
#include <stdlib.h> // for rand()
#include <time.h> // for clock()
#include <math.h>
#include "mat.h"

/* matrix times matrix */
void matmat(double A[], int am, int an,
            double B[], int bm, int bn, 
            double (*C)[], int cm, int cn)
{
  int i,j,k;
  double accum;
  // check resultant matrix size
  if (cm != am || cn != bn) {
    printf("Error: invalid size of resultant matrix\n");
    return;
  } 

  // check compatibility of A and B
  if (an != bm) { 
    printf("Error: input matrices non-conformable\n");
    return;
  }

  // go
  for (i=0; i<am; i++) {
    for(j=0; j<bn; j++) {
      accum = 0;
      for(k=0; k<an; k++) {
        accum += A[(i*an)+k]*B[(k*bn)+j];
        if (isnan(accum)) 
          { 
            printf("nan at i=%d, j=%d\n", i,j);
          }
        //printf("%lf ", (*C)[(i*cn)+j]);
      }
      (*C)[(i*cn)+j] = accum;
      //printf("\n");
    }
    //printf("\n");
  }
}


/* matrix times vector */
void matvec(double A[], int am, int an,
            double x[], int xn,
            double (*b)[])
{
  int i, j;

  // check compatibility
  if (an != xn) {
    printf("Error: non-conformable arrays\n");
    return;
  }
  
  for (i=0; i<am; i++) {
    for(j=0; j<an; j++) {
      (*b)[i] += A[(i*an)+j]*x[j];
    }
  }
}

/* print matrix */
void printmat(double A[], int m, int n)
{
  int i,j;
  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
      printf("%.05lf ", A[i*n+j]);
    }
    printf("\n");
  }
  printf("\n");
}

/* dot product */
double dotprod(double x[], double y[], int n)
{
  int i;
  double result = 0;

  for (i=0; i<n; i++)
    {
      result += x[i]*y[i];
    }
  return result;
}

/* 2 norm */
double eucnorm(double x[], int n)
{
  double result = 0;
  int p = 2;
  int i;

  for (i=0; i<n; i++)
    {
      result += pow(x[i], p);
    }

  return sqrt(result);
}

/* euclidean metric, i.e. ||x-y||_2 */
double eucmetric(double x[], double y[], int n)
{
  int i;
  int p = 2;
  double result = 0;

  for (i=0; i<n; i++)
    {
      result += pow(x[i]-y[i], p);
    }
  return sqrt(result);
}

/* add vectors u and v in R^n 
   stores result in u
*/
void addVec(double (*u)[], double v[], int n)
{
  int i;
  for (i=0; i<n; i++)
    {
      (*u)[i] = (*u)[i] + v[i];
    }
}

/* scale vector in R^n by a */
void scaleVec(double (*x)[], double a, int n)
{
  int i;

  for (i=0; i<n; i++)
    {
      (*x)[i] = (*x)[i] * a;
    }
}

/* insert vector v into column i of matrix A stored in row-major order */
void updateCol(double (*A)[], int m, int n, int i, double v[])
{
  int k;
  for (k=0; k<m; k++)
    {
      (*A)[(k*n)+i] = v[k];
    }
}

/* insert vector v into row j of matrix A stored in row-major order */
void updateRow(double (*A)[], int m, int n, int j, double v[])
{
  int k;
  for (k=0; k<n; k++)
    {
      (*A)[(j*m)+k] = v[k];
    }
}

/* transpose: At must be pointer to m*n-length array */
/* result: pointer to matrix size n \times m in row-major order */
void transpose1(double A[], int m, int n, double (*At)[])
{
  int i, j;
  for (i=0; i<m; i++)
    {
      for (j=0; j<n; j++)
        {
          (*At)[m*j+i] = A[n*i+j];
        }
    }
 
}

// test pointer version for speed
void transpose2(double (*A)[], int m, int n, double (*At)[])
{
  int i, j;
  for (i=0; i<m; i++)
    {
      for (j=0; j<n; j++)
        {
          (*At)[m*j+i] = (*A)[n*i+j];
        }
    }
 
}

/* extract the row i, col j entry from matrix stored as vector 
 i, j indexed from zero 
*/
double extractEntry(double A[], int m, int n, int i, int j)
{
  if (j<n && i<m)
    return A[n*i+j];
  else
    printf("extractEntry: invalid entry i: %d  j: %d\n", i, j);
  return 0;
}

/* put column j of A into vector v */
void extractColumn(double A[], int m, int n, int j, double (*v)[])
{
  //double extractEntry(double A[], int m, int n, int i, int j);
  int i;

  for (i=0; i<m; i++)
    (*v)[i] = extractEntry(A, m, n, i, j);
}

/* put row j of A into vector v */
void extractRow(double A[], int m, int n, int j, double (*v)[])
{
  int i;

  for (i=0; i<n; i++)
    (*v)[i] = extractEntry(A, m, n, j, i);
}

// insert x into row i, col j entry in matrix A
void updateEntry(double x, double (*A)[], int m, int n, int i, int j)
{
  if (j<n && i<m)
    (*A)[n*i+j] = x;
  else
    printf("updateEntry: invalid entry i: %d  j: %d\n", i, j);
}

// print maximum entry of vector
double printMax(double v[], int n)
{
  int i=0;
  int loc=0;
  double max=0;

  while (i<n) 
    {
      if (fabs(v[i]) > fabs(max))
        {
          max = v[i];
          loc = i;
        }
      i++;
    }
     
  printf("max value: %.15lf location: %d\n", max, loc);
  return max;
}

// sign of double
int sign(double x)
{
  if (x < 0)
    {
      return -1;
    }
  else
    return 1;
}

// extract submatrix of A, store in As
// include row i1 up to but not including i2, 
// and column j1 up to but not including j2
void submat(double A[], int m, int n, int i1, int i2, int j1, int j2, 
            double (*As)[])
{
  int i, j, k, l;
  int Aslen = i2 - i1;
  int Aswid = j2 - j1;
  k = l = 0;

  //printf("%d %d\n", Aslen, Aswid);
  // WARNING: no size check for As..

  for (i=i1; i<i2; i++)
    {
      for (j=j1; j<j2; j++)
        {
          // printf("k,l inside submat loop: %d %d\n", k,l);
          updateEntry(extractEntry(A,m,n,i,j), As, Aslen,Aswid,k,l);
          l++;
        }
      k++;
      l = 0;
    }
}


// insert matrix As into bigger matrix A
void insertSubmat(double (*A)[], int m, int n, int i1, int i2, int j1, int j2, 
            double As[])
{
  int i, j, k, l;
  int Aslen = i2 - i1;
  int Aswid = j2 - j1;
  k = l = 0;

  //printf("Aslen: %d Aswid: %d\n", Aslen, Aswid);
  // WARNING: no size check for As..

  for (i=i1; i<i2; i++)
    {
      for (j=j1; j<j2; j++)
        {
          //printf("k,l inside submat loop: %d %d\n", k,l);
          updateEntry(extractEntry(As,Aslen,Aswid,k,l), A,m,n,i,j);
          l++;
        }
      k++;
      l = 0;
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

