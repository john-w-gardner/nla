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
      for(k=0; k<an; k++) {
        (*C)[(i*cn)+j] += A[(i*an)+k]*B[(k*bn)+j];
        //printf("%lf ", (*C)[(i*cn)+j]);
      }
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
    printf("extractEntry: invalid entry i: %d  j: %d", i, j);
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
    printf("updateEntry: invalid entry i: %d  j: %d", i, j);
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

