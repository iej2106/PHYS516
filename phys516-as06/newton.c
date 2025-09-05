/*******************************************************************************
Newton root finding method for chemical potential of Fermi distribution.

USAGE

%cc newton.c -o newton -lm
%./newton < newton.in > newton.out
*******************************************************************************/
#include <stdio.h>
#include <math.h>
#define NMAX 256 /* # of eigen states */
#define M 256 /* # of valence electrons */
#define MAX_ITER 100 /* Max iterations in the Newton method */
#define TOL 1.0e-10 /* Error tolerance */
#define ARGMAX 20.0

int main(int argc, char **argv) {
	double e[NMAX]; /*Eigenergies*/
    double mu, dmu; /*Chemical potential*/
    double f, df, fn;
    int n, iter = 1, dummy;
    double T; /* Temperature */

    /* Read the eigenenergies */
    for(n = 0; n < NMAX; n++){
        scanf("%d %le", &dummy, &e[n]);
    }
    /* Read the temperature */
    scanf("%le",&T);

    /* Initial guess of chemical potential */
    mu = e[NMAX/2-1];
    while(iter<=MAX_ITER){
        f = df = 0.0;
        for(n = 0; n < NMAX; n++){
            fn = 2.0/(exp((e[n]-mu)/T)+1.0);
            f += fn;
            df += fn*(1.0-0.5*fn)/T;
        }
        f -= M;
        dmu = -f/df;

        printf("%d %le %le %le\n",iter,mu,dmu,f);

        if(fabs(dmu) < TOL){
            break;
        }
        else{
            mu += dmu;
        }
        ++iter;
    }

    if(iter > MAX_ITER){
        printf("Newton method did not converge!\n");
    }

    printf("MU = %le\n",mu);

    for(n = 0; n < NMAX; n++){
        fn = 2.0/(exp((e[n]-mu)/T)+1.0);
        printf("%d\t%le\t%le\n",n,e[n],fn);
    }

	return 0;
}