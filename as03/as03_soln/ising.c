/* Monte Carlo integration of Ising model */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define L 20
int s[L][L];
double exp_dV[2][5];
double JdivT;
double HdivT;
int hist[2*L*L+1];

void table_set() {
  int k, l,s, s_new, s_neighbor;
  for (k=0; k<2; k++) {
	s_new = 2*k -1;
		for (l=0; l<5; l++) {
			s_neighbor = 2*l-4;
			exp_dV[k][l] = exp(2.0*s_new*(JdivT*s_neighbor+HdivT));
	}
  }
}

int main() {
//double x, pi, sum = 0.0;
//int try, ntry;
  double runM;
  double sumM = 0.0, sumM2 = 0.0;
  int Sta_step;
  double avgM, sigM, exp_val;
  int i,j,step,im,ip,jm,jp,s_new,k,s_neighbor,l;
  
  printf("Input the number of MC trials\n");
  scanf("%d",&Sta_step);
  printf("Input JdivT/n");
  scanf("%le",&JdivT);
  printf("Input HdivT/n");
  scanf("%le",&HdivT);  
  
  table_set();
  
  for (i=0; i<L; i++)
  	for (j=0; j<L; j++)
  		s[i][j] = 1; // Cold start
  runM = 1.0*L*L; 	
  for (i=0; i<2*L*L+1; i++) hist[i] = 0;
  		
  srand((unsigned)time((long *)0));
  
//for (try=0; try<ntry; try++) {
//  x = rand()/(double)RAND_MAX;
//  sum += 4.0/(1.0 + x*x);
//}
  for (step = 1; step<=Sta_step; step++) {
  	//random select a grid point, (i,j)
	i = rand()%L;
	j = rand()%L;
	im = (i + L - 1) % L;  
  	ip = (i + 1) % L;
  	jm = (j + L - 1) % L;
  	jp = (j + 1) % L;
	//compute the change in potential energy, DV
	//flip, si, j  >  si, j
	s_new = -s[i][j];
	k = (1+s_new)/2;
	s_neighbor = s[im][j]+s[ip][j] + s[i][jm] + s[i][jp];
	l = (4+s_neighbor)/2;
	exp_val = exp_dV[k][l];
	//if dV<0 accept the flip, si,j  > so,j
	if (exp_val > 1.0) {
	  s[i][j] = s_new;
	  runM += 2.0*s_new;	
	}
	//else if random() <=exp(-dV/kBT) then // 0 < random() < 1
	//accept the flip, si,j -> -si,j
	else if (rand()/(double)RAND_MAX < exp_val) {
	  s[i][j] = s_new;
	  runM += 2.0*s_new;
	}
	//endif
	//Sum_A = Sum_A + A(sN)  // Sample physical quantity A(sN)
	sumM += runM;
	sumM2 += runM*runM;
	++hist[(int)runM+L*L];
  }
  avgM = sumM/Sta_step;
  sigM = sqrt(sumM2/Sta_step - avgM*avgM);
  printf("Magnetization = %le +- %le\n", avgM, sigM);	
  for (i=0; i<2*L*L+1; i++)
  	printf("%d\t%d\n",i-L*L,hist[i]);
  return 0;
}
