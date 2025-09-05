#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define MAX_DAY 365
#define DT (1.0/MAX_DAY)
#define MU 0.14
#define SIGMA 0.2
#define N_HIST 50
#define N_WALKER 100000
#define S_INIT 20.0

double rand_normal() {
  double r1,r2;
  r1 = rand()/(double)RAND_MAX;
  r2 = rand()/(double)RAND_MAX;
  return sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2);
}

int main() {
  int day,walker,k;
  double s; /* stock price */
  int hist[N_HIST];
  FILE *pst;
  
  pst = fopen("st.d","w");

  double mu_dt, sigma_rtdt;
  mu_dt = MU*DT;
  sigma_rtdt = SIGMA*sqrt(DT);
  
  for (k=0; k<N_HIST; k++)
    hist[k] = 0.0;

  srand((unsigned)time((long *)0)); /* Initialize the rondom-number sequence */

  for (walker=1; walker<=N_WALKER; walker++) {
    s = S_INIT;
    for (day=1; day<=MAX_DAY; day++) {
    	s += s*(mu_dt+sigma_rtdt*rand_normal());
    	if (walker == 1) fprintf(pst,"%d\t%f\n",day,s);
    	if (s <= 0.0) break;
    } /* Endfor step */
    s = s > 0.0 ? s : 0.0;
    k = (int)s;
    k = k < N_HIST ? k : N_HIST-1;
    ++hist[k];
  } /* Endfor walker */

  for (k=0; k<N_HIST; k++)
    printf("%d\t%d\n",k,hist[k]);
  fclose(pst);

  return 0;
}