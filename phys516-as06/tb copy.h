/* To be put in a header file *************************************************/
#define NMAX 100           /* Max # of atoms */
#define NAUC 8             /* # of atoms per unit cell */
#define LCNS (1*5.43)    /* Lattice constant of Si in angstrom */
#define Nbin 500		   /*Number of energy bins */
int nAtom;                 /* # of atoms */
double r[NMAX][3];         /* r[i][0|1|2] is the x|y|z coordinate of atom i */
int InitUcell[3];          /* # of unit cells */

/* Function pro t o t y p e s *************************************************/
void InitConf();
double SignR(double v, double x) {if (x>0) return v; else return -v;}
double **dmatrix (int nrl, int nrh, int ncl, int nch);
double *dvector (int nl, int nh);
void tred2(double **a, int n, double d[], double e[]);
void tqli(double d[], double e[], int n, double **z);