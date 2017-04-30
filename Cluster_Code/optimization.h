#ifndef OPTIMIZATION_H_
#define OPTIMIZATION_H_

#include "zone.h"
#include "inference.h"

#define PI 3.1415926535897932384626433832795028841971693993751


#define TINY 1.0e-2
#define NMAX 10000
#define GET_PSUM							\
for (j=0; j< ndim ; j++) {						\
for (sum=0.0, i=0; i < mpts; i++ ) sum += p[i][j];			\
psum[j]=sum; }
#define SWAP(a,b) {swap= (a); (a)=(b); (b)=swap;}

#ifndef isnan
	#define isnan(a) (a != a)
#endif

#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0

static double maxarg1,maxarg2;
#define FMAX(a,b) ( maxarg1 = (a), maxarg2 = (b) , (maxarg1) > (maxarg2) ? \
(maxarg1) : (maxarg2) )

static double sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

#define ITMAX 10000
#define EPS 3.e-11
#define TOLY (4.0*EPS)
#define STPMX 500.0

#define TOL 2.0e-4
#define ALF 1.0e-4
#define TOLX 1.0e-8

#define IA_bis 16807
#define IM_bis 2147483647
#define AM_bis (1.0/IM_bis)
#define IQ_bis 127773
#define IR_bis 2836
#define NTAB_bis 32
#define NDIV_bis (1+(IM_bis-1)/NTAB_bis)
#define EPS_bis 1.2e-7
#define RNMX_bis (1.0-EPS_bis)
#define IM1_bis 2147483563
#define IM2_bis 2147483399
#define AM2_bis (1.0/IM1_bis)
#define IMM1_bis (IM1_bis-1)
#define IA1_bis 40014
#define IA2_bis 40692
#define IQ1_bis 53668
#define IQ2_bis 52774
#define IR1_bis 12211
#define IR2_bis 3791
//#define NTAB 32     attention changements ici
#define NDIV2_bis (1+IMM1_bis/NTAB_bis)
#define EPS2_bis 1.e-5
#define RNMX2_bis (1.0-EPS2_bis)
#define JMAX_bis 30
#define isNaN(x) ((x) != (x))

#define NR_END 1
#define FREE_ARG char*

#define GTOL 1.e-16

extern         Zone *ZONES;
extern char    INPUT_FILE[FILENAME_MAX];
extern double  LOCALIZATION_ERROR;
extern char    INFERENCE_MODE[FILENAME_MAX];
extern double  D_PRIOR, V_PRIOR;
extern         FILE *INPUT;
extern int     NUMBER_OF_ZONES;
extern int     CURRENT_ZONE;
extern double  **HESSIAN;
extern int     DIMENSIONS;
/*******************************************************************************/
// finds minimum of function using Broyden-Fletcher-Goldfarb-Shanno quasi-Newton method
void dfpmin(double *p,int ndim,  double gtol, int *iter, double *fret, double (func)(double *), void (dfunc)(double (func)(double *),double *, double *) );
/*******************************************************************************/
// performs a line search
void lnsrch(int ndim,double xold[], double fold, double g[] , double p[], double xx[] , double f[], double stpmax, int check[],  double (func)(double [])  );
/*******************************************************************************/
// calculates gradient of a passed function
void dfunc(double (func)(double *), double *yy, double *ans);
/*******************************************************************************/
// polynomial interpolation function
void vander(double *xxx, double *cof, double *yyy, int nn);
/*******************************************************************************/
#endif /* OPTIMIZATION_H_ */
