#include <stdio.h>
#include <stdlib.h> // for rand()
#include <time.h> // for clock()

int main()
{
  int m = 1043+434;
  int n = m*m;
  float v[n];
  int i;
  for (i=0; i<n; i++)
    {
      v[i] = (float)rand()/RAND_MAX*20-10;
    }
 
  printf("first entry of v: %.15f\n", v[0]);

}


/* biggest square matrix of doubles appears to be 1043*1043
 and the biggest square matrix of floats is 1043+a where a solves 2x^2 = (x+a)^2, 
 ie sqrt(2)-1
*/
