/******************************************************************************* 
Tight−binding (tb) model for silicon in eV and angstrom units .
USAGE
%cc tb.c eigen.c −o tb −lm
%./tb > tb.out 
*******************************************************************************/
#include <stdio.h> 
#include <math.h> 
#include "tb.h"

int main(int argc, char **argv) { 
		 int n4, i, j, k, nu, a, l;
		 double** h; 
		 double* d;
		 double* e;
		 double r0, n, Es, Ep, dummy, minEps, maxEps, sigma, Sum, rija , La, distance;
		 double Eps[Nbin];
		 double D[Nbin];
		 double dij[3];
		 double hlambda0[4], hlambda[4], nlambda[4] , rlambda[4];

/* Parameters for the silicon tb model in Kwon’s paper */
r0=2.360352; n=2.0; Es=-5.25; Ep=1.20;
hlambda0[0] = -2.038; 
hlambda0[1] = 1.745; 
hlambda0[2] = 2.75;
hlambda0[3] = -1.075;
nlambda[0] = 9.5;
nlambda[1] = 8.5;
nlambda[2] = 7.5;
nlambda[3] = 7.5;
rlambda[0] = 3.4;
rlambda[1] = 3.55;
rlambda[2] = 3.7;
rlambda[3] = 3.7;

/* min and max energies in eV gaussian function */
minEps = -15.0; maxEps =10.0; sigma =0.1;

InitConf (); /* Read input parameters */ \
n4 = 4*nAtom;
h = dmatrix(1,n4,1,n4);
d = dvector(1,n4);
e = dvector(1,n4);

/* Fill in first block of tb matrix */ 
for(i=1; i<=n4; i++){
	for(j=1; j<=n4; j++){
		h[i][j] = 0.0;
						}
}

/* Fill in block diagonals of tb matrix */ 
for(i=0; i<nAtom; i++){
	h[4*i+1][4*i+1] = Es; for (k=2; k<5; k++){
		h[4*i+k][4*i+k]=Ep;						
		}
	
	for ( j=i +1; j<nAtom; j++){ /* calculate atom distances and unit bond−direction vector */
		distance = 0.0;
		for(a=0; a<3; a++){
			rija =r[i][a] - r[j][a]; 
			La = LCNS*InitUcell[a];
			rija = rija-SignR(0.5*La, rija- 0.5*La) - SignR(0.5*La, rija + 0.5*La);
			dij[a]=rija;
			distance = distance + rija*rija;
			}
			distance = sqrt(distance);
		for(a=0; a<3; a++){ /* Normalize for unit bond−direction vector */
			dij[a] = dij[a]/distance;
			}
		for(k=0; k<4; k++){ /* calculate h {\lambda} */
			hlambda[k] = exp(n*(-pow(distance/rlambda[k], nlambda[k]) + pow(r0/rlambda[k], nlambda[k])));
			hlambda[k] = hlambda0[k] * pow(r0/distance,n) * hlambda [k];
			}
			/* Fill in off−diagonal block */
			h[4*i+1][4*j+1] = hlambda[0]; /* left−top−corner entry of off−diagonal block */
			/* Fill in the rest of off−diagonal block */
			for(a=0; a<3; a++){
				l = a+2;
				h[4*i+1][4*j+l] = dij[a]*hlambda[1]; h[4*i+l][4*j+1] =-h[4*i+1][4*j+l]; /* asymmetry in elements */
				for(k=l; k<5; k++){ h[4*i+l][4*j+k] =
				dij[a]*dij[k-2]*(hlambda[2] - hlambda[3]); h[4*i+k][4*j+l] = h[4*i+l][4*j+k]; /* symmetry in elements */
}
			h[4*i+l][4*j+l] = h[4*i+l][4*j+l] + hlambda[3]; 
			/* add to complete filling in diag of off−diag block */
}
			for(l=1; l<5; l++){
				for (k=1; k<5; k++){
					h[4*j+l][4*i+k]= h[4*i+k][4*j+l];
				}
}
}
}

tred2(h,n4,d,e); 
tqli(d,e,n4,h);

for(i=1; i<n4; i++){
	for(j=i+1; j<=n4; j++){
		if(d[i] > d[j]){
		dummy = d[i];
		d[i]=d[j]; 
		d[j]=dummy;
		}
	}
}

for(i=0; i<Nbin; i++){
	Eps[i] = minEps + i*(maxEps-minEps)/Nbin; 
	Sum = 0.0;
		for(nu=0; nu<n4; nu++){
			Sum = Sum+
			exp(-(Eps[i]-d[nu])*(Eps[i]-d[nu])/(sigma*sigma));
			}
		D[i] =Sum/ (sqrt(M_PI)*sigma);
		printf("%le\t%le\n",Eps[i],D[i]);
		}
		
//for(i=1; i<=n4; i++){
// printf(”%d\t%le\n”,i ,d[ i ]); //}

return 0;


}

/*----------------------------------------------------------------------------*/
void InitConf() {
/*------------------------------------------------------------------------------
	r are initialized to diamond lattice positions.  
------------------------------------------------------------------------------*/
	double gap[3];      /* Unit cell size */
	double c[3];
	int j,k,nX,nY,nZ;
	/* Atom positions in a unit diamond crystalline unit cell */
	double origAtom[NAUC][3] = {{0.0, 0.0, 0.0 }, {0.0, 0.5, 0.5 },
                                {0.5, 0.0, 0.5 }, {0.5, 0.5, 0.0 },
                                {0.25,0.25,0.25}, {0.25,0.75,0.75},
                                {0.75,0.25,0.75}, {0.75,0.75,0.25}};

	/* Read the # of unit cells in the x, y & z directions */
	scanf("%d%d%d",&InitUcell[0],&InitUcell[1],&InitUcell[2]);

	/* Sets up a diamond lattice */
	for (k=0; k<3; k++) gap[k] = LCNS;
	nAtom = 0;
	for (nZ=0; nZ<InitUcell[2]; nZ++) {
		c[2] = nZ*gap[2];
		for (nY=0; nY<InitUcell[1]; nY++) {
			c[1] = nY*gap[1];
			for (nX=0; nX<InitUcell[0]; nX++) {
				c[0] = nX*gap[0];
				for (j=0; j<NAUC; j++) {
					for (k=0; k<3; k++)
						r[nAtom][k] = c[k] + gap[k]*origAtom[j][k];
					++nAtom;
				}
			}
		}
	}
}


