/* Monte Carlo integration of PI by sample mean */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
int main() {
  double x, pi, sum = 0.0, pi2, sum2 =0.0, stdv, fx;
  int try, ntry;
  printf("Input the number of MC trials\n");
  scanf("%d",&ntry);
  srand((unsigned)time((long *)0));
  for (try=0; try<ntry; try++) {
    x = rand()/(double)RAND_MAX;
    fx = 4.0/(1.0 + x*x);
    sum += fx;
    sum2 += fx*fx;
  }
  pi = sum/ntry;
  pi2 = sum2/ntry;
  stdv = sqrt((pi2-pi*pi)/(ntry-1));
  printf("MC estimate for PI = %f +- %e\n", pi, stdv);
  return 0;
}
