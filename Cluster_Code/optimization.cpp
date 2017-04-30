#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "optimization.h"
/*******************************************************************************/
/*******************************************************************************/
/*** optimization VMMM ****/
/*******************************************************************************/
/*******************************************************************************/
void dfpmin(double *p,int ndim,  double gtol, int *iter, double *fret, double (func)(double *), void (dfunc)(double (func)(double *),double *, double *) ) {
	int check, i, ii, its, j;
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg, sumxi,temp,test;
	double* dg = new double[ndim];
	double* g = new double[ndim];
	double* hdg = new double[ndim];
	double* pnew = new double[ndim];
	double* xi = new double[ndim];

	char progress[100];

	for (ii=0; ii < ndim; ii++) {
		dg[ii] = 0.0;
		g[ii] = 0.0;
		hdg[ii] = 0.0;
		pnew[ii] = 0.0;
		xi[ii] = 0.0;
	}

	fp = (func)(p) ;
  ////fprintf(stderr, "%le\n", fp);
////fprintf(stderr, "%le\t %le\n", 10.53, log(10.53));
	(dfunc)((func),p,g);

	for (i = 0; i < ndim; i++){
		for (j=0; j<ndim; j++) {
			HESSIAN[i][j]=0.0;
		}
		HESSIAN[i][i]=1.0;
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}

	stpmax = STPMX*FMAX(sqrt(sum), (double)ndim);

	for ( its=0; its<ITMAX; its++) {
//	for ( its=0; its<0; its++) {
		iter[0]=its;

		////fprintf(stderr,"Iteration %i:\t%.3f\n",iter[0]+1,fp);

		lnsrch(ndim,p,fp,g,xi,pnew,fret,stpmax,&check,(func));

		fp = fret[0];
      //  //fprintf(stderr, "%le\n", fp);

		for (i=0; i<ndim; i++){
			xi[i] = pnew[i]-p[i];
			p[i] = pnew[i];
		}
		test = 0.0;
		for (i=0 ; i<ndim ; i++){
			temp = fabs(xi[i])/FMAX(fabs(p[i]), 1.0);
			if (temp > test) {
				test = temp;
			}
		}

		if (test<TOLY){
			delete [] dg;
			delete [] g;
			delete [] hdg;
			delete [] pnew;
			delete [] xi;

			return;
		}

		for (i=0 ; i<ndim ; i++){ dg[i]=g[i];}
		(dfunc)((func),p,g);
		test = 0.0;

		den = FMAX(fret[0],1.0);
		for (i=0; i<ndim ;  i++){
			temp = fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
			if (temp > test) {test =temp;}
		}

		if (test < gtol){
			delete [] dg;
			delete [] g;
			delete [] hdg;
			delete [] pnew;
			delete [] xi;

			return;
		}

		for (i=0; i<ndim; i++) {dg[i] = g[i]-dg[i];}

		for (i=0; i<ndim; i++) {
			hdg[i] = 0.0;
			for (j=0; j<ndim; j++){ hdg[i] += HESSIAN[i][j]*dg[j];}
		}

		fac =fae = sumdg = sumxi = 0.0;
		for (i=0; i<ndim ; i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}

		if (fac > sqrt(EPS*sumdg*sumxi)){
			fac = 1.0/fac;
			fad = 1.0/fae;
			for (i=0; i<ndim; i++ ){ dg[i]= fac*xi[i]-fad*hdg[i];}
			for (i=0; i<ndim ;i++){
				for (j=i ; j<ndim ; j++ ){
					HESSIAN[i][j] += fac*xi[i]*xi[j]
					- fad*hdg[i]*hdg[j] + fae*dg[i]*dg[j];
					HESSIAN[j][i] = HESSIAN[i][j];
				}
			}
		}
		for (i=0; i<ndim ; i++) {
			xi[i] =0.0;
			for (j=0; j<ndim; j++) xi[i] -= HESSIAN[i][j]*g[j];
		}
	}
	delete [] dg;
	delete [] g;
	delete [] hdg;
	delete [] pnew;
	delete [] xi;

	////fprintf(stderr,"dfpmin fail\n"); return; exit(-1);
}
/*******************************************************************************/
/*******************************************************************************/
/**** Line Search ****/
/*******************************************************************************/
/*******************************************************************************/

void lnsrch(int ndim,double xold[], double fold, double g[] , double p[], double xx[] , double f[], double stpmax, int check[],  double (func)(double [])  ){

	int i;
	double a,alam, alam2, alamin, b, disc, f2, rhs1, rhs2, slope, sum, temp, test, tmplam;

	check[0]=0;
	for (sum=0.0, i=0; i<ndim; i++){ sum += p[i]*p[i];}
	sum = sqrt(sum);

	if (sum > stpmax)
		for (i=0; i<ndim;i++) p[i] *= stpmax/sum;
	for (slope=0.0, i=0;i<ndim;i++ ) {
		slope += g[i]*p[i];
	}

	if (slope >= 0.0) {
		////fprintf(stderr,"Line Search Failure\n");
		return;
		exit(1);
	}
	test = 0.0;
	for (i=0; i<ndim; i++){
		temp = fabs(p[i])/FMAX(fabs(xold[i]),1.0);
		if (temp > test) test = temp;
	}
	alamin = TOLX/test;
	alam = 1.0;
	for (;;) {
		for (i=0;i<ndim;i++) xx[i]=xold[i]+alam*p[i];

		f[0] =   (func)( xx );


		if (alam < alamin) {
			for (i=0; i<ndim ; i++){ xx[i] = xold[i];}
			check[0]=1;
			return ;
		} else if (f[0] <= fold +ALF*alam*slope) {
			return;
		}
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(f[0] - fold -slope));
			else {
				rhs1 = f[0] -fold -alam*slope;
				rhs2 = f2 - fold -alam2*slope;
				a = (rhs1/(alam*alam) - rhs2/(alam2*alam2))/(alam - alam2);
				b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0){ tmplam = -slope/(2.0*b);}
				else {
					disc = b*b-3.0*a*slope;
					if (disc < 0.0){ tmplam = 0.5*alam;}
					else if (b <=0.0){ tmplam = (-b +sqrt(disc))/(3.0*a);}
					else {tmplam = -slope/(b+sqrt(disc));}
				}

				if (tmplam > 0.5 *alam)
					tmplam = 0.5 *alam;
			}
		}
		alam2 = alam;
		f2 =f[0];
		alam = FMAX(tmplam, 0.1*alam);
	}

}
/*******************************************************************************/
/*******************************************************************************/
/*** dereiv function multi dimentions ***/
/*******************************************************************************/
/*******************************************************************************/

void dfunc(double (func)(double *), double *yy, double *ans){

	double h=1.e-8, ERR;
	int i,j, indice, indice2;
	double errt,fac,hh;

	double* veclocal1 = new double[DIMENSIONS];
	double* veclocal2 = new double[DIMENSIONS];

	double a[NTAB][NTAB];

	for (int aa = 0; aa < NTAB; aa++) {
		for (int bb = 0; bb < NTAB; bb++) {
			a[aa][bb] = 0.0;
		}
	}

	for (indice=0; indice<DIMENSIONS ; indice++){

		hh=h;
		for (indice2=0; indice2<DIMENSIONS ; indice2++ ){
			if (indice2 == indice){
				veclocal1[indice2]= yy[indice2]+hh;
				veclocal2[indice2]= yy[indice2]-hh;
			} else {
				veclocal1[indice2] = yy[indice2];
				veclocal2[indice2] = yy[indice2];
			}
		}

		a[0][0]=( (func)(veclocal1 ) - (func)(veclocal2)  )/(2.0*hh);

		ERR=BIG;

		for(i=1 ; i<NTAB ;  i++ ){
			hh /= CON;

			for (indice2=0; indice2<DIMENSIONS ; indice2++ ){
				if (indice2 == indice){
					veclocal1[indice2]= yy[indice2]+hh;
					veclocal2[indice2]= yy[indice2]-hh;
				} else {
					veclocal1[indice2] = yy[indice2];
					veclocal2[indice2] = yy[indice2];
				}
			}

			a[0][i]=( (func)(veclocal1) - (func)(veclocal2)  )/(2.0*hh);

			fac = CON2;
			for (j=1 ; j<i; j++){
				a[j][i]= (a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
				fac = CON2*fac;
				errt=FMAX( fabs(a[j][i]-a[j-1][i]) , fabs(a[j][i]-a[j-1][i-1])  );
				if (errt <= ERR){
					ERR = errt;
					ans[indice] = a[j][i];
				}
			}
			if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(ERR)) break;
		}

	}

	delete [] veclocal1;
	delete [] veclocal2;

	return ;
}
/*******************************************************************************/
/*******************************************************************************/
/*** vandermone matrix *****/
/*******************************************************************************/
/*******************************************************************************/
void vander(double *xxx, double *cof, double *yyy, int nn){
	// resout x^{k}*w=V input x et V et cela donne w
	int k,j,i;
	double phi,ff,b;

	double* s = new double[nn+1];
	for (i=0;i<=nn;i++) s[i]=cof[i]=0.0;
	s[nn] = -xxx[0];
	for (i=1;i<=nn;i++) {
		for (j=nn-i;j<=nn-1;j++)
			s[j] -= xxx[i]*s[j+1];
		s[nn] -= xxx[i];
	}
	for (j=0;j<=nn;j++) {
		phi=nn+1;
		for (k=nn;k>=1;k--)
			phi=k*s[k]+xxx[j]*phi;
		ff=yyy[j]/phi;
		b=1.0;
		for (k=nn;k>=0;k--) {
			cof[k] += b*ff;
			b=s[k]+xxx[j]*b;
		}
	}
	delete [] s;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
