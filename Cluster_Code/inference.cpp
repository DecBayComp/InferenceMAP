#include <math.h>

#include "zone.h"
#include "inference.h"
#include "optimization.h"

#define PI 3.1415926535897932384626433832795028841971693993751

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
extern int     JEFFREYS_PRIOR;
extern         Zone* findZone(int id);
extern         Zone* findZone2(int id);

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
// zonal optimization
double dfPosterior(double *FxFyD) {

	double D_bruit;
	double result;

	result = 0.0;

	#pragma omp for
	for (int i = 0; i < ZONES[CURRENT_ZONE].translocations; i++) {
		// "next" index does not necessarily lie in the "current zone"
		const double dt = ZONES[CURRENT_ZONE].dt[i];
		const double dx = ZONES[CURRENT_ZONE].dx[i];
		const double dy = ZONES[CURRENT_ZONE].dy[i];

		D_bruit = LOCALIZATION_ERROR*LOCALIZATION_ERROR/dt;

		result += - log(4.0*PI*(FxFyD[2]+D_bruit)*dt ) - pow(fabs(dx-FxFyD[2]*FxFyD[0]*dt ),2.0)/(4.0*(FxFyD[2]+D_bruit)*dt) - pow(fabs(dy-FxFyD[2]*FxFyD[1]*dt ),2.0)/(4.0*(FxFyD[2]+D_bruit)*dt);
	}

	// Jeffrey's Prior
	if (JEFFREYS_PRIOR == 1) {
		result += 2.0*log(FxFyD[2]) - 2.0*log(FxFyD[2]*ZONES[CURRENT_ZONE].dtMean + LOCALIZATION_ERROR*LOCALIZATION_ERROR);
	}
	return -result;

}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
double dPosterior(double *D) {

	double D_bruit;

	double result = 0.0;

	#pragma omp for
	for (int i = 0; i < ZONES[CURRENT_ZONE].translocations; i++) {
		// "next" index does not necessarily lie in the "current zone"
		const double dt = ZONES[CURRENT_ZONE].dt[i];
		const double dx = ZONES[CURRENT_ZONE].dx[i];
		const double dy = ZONES[CURRENT_ZONE].dy[i];

		D_bruit = LOCALIZATION_ERROR*LOCALIZATION_ERROR/dt;

		result += - log(4.0*PI*(D[0]+D_bruit)*dt) - (dx*dx)/(4.0*(D[0]+D_bruit)*dt) - (dy*dy)/(4.0*(D[0]+D_bruit)*dt);
	}

	// Jeffrey's Prior
	if (JEFFREYS_PRIOR == 1) {
		result += - (D[0]*ZONES[CURRENT_ZONE].dtMean + LOCALIZATION_ERROR*LOCALIZATION_ERROR);
	}
	return -result;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
double dDDPosterior(double *DD) {

	double result;

	result = 0.0;

	#pragma omp for
	for (int a = 0; a < NUMBER_OF_ZONES; a++) {
		ZONES[a].gradDx = dvGradDx(DD,a);
		ZONES[a].gradDy = dvGradDy(DD,a);
		ZONES[a].priorActive = true;
	}


	#pragma omp for
	for (int z = 0; z < NUMBER_OF_ZONES; z++) {
		const double gradDx = ZONES[z].gradDx;
		const double gradDy = ZONES[z].gradDy;
		const double D = DD[z];

		for (int j = 0; j < ZONES[z].translocations; j++) {
			const double dt = ZONES[z].dt[j];
			const double dx = ZONES[z].dx[j];
			const double dy = ZONES[z].dy[j];
			const double Dnoise = LOCALIZATION_ERROR*LOCALIZATION_ERROR/dt;

			result += - log(4.0*PI*(D + Dnoise)*dt) - ( dx*dx + dy*dy)/(4.0*(D+Dnoise)*dt);
		//	////fprintf(stderr, "diffusion %lf\t diffusion noise %lf\t value result %lf\n", D,Dnoise,result );
			// Priors
		}



		if (ZONES[z].priorActive == true) {
		//	////fprintf(stderr, "here\n");
			result -= D_PRIOR*(gradDx*gradDx*ZONES[z].areaX + gradDy*gradDy*ZONES[z].areaY);
			if (JEFFREYS_PRIOR == 1) {
				// result += 2.0*log(D) - 2.0*log(D*dt + LOCALIZATION_ERROR*LOCALIZATION_ERROR);
				result += 2.0*log(D) - 2.0*log(D*ZONES[z].dtMean + LOCALIZATION_ERROR*LOCALIZATION_ERROR);
			}
		}


	}

	return -result;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
double dvPosterior(double *DV) {

	double result;

	result = 0.000;

	#pragma omp for
	for (int a = 0; a < NUMBER_OF_ZONES; a++) {
		ZONES[a].gradVx = dvGradVx(DV,a);
		ZONES[a].gradVy = dvGradVy(DV,a);
		ZONES[a].gradDx = dvGradDx(DV,a);
		ZONES[a].gradDy = dvGradDy(DV,a);
		ZONES[a].priorActive = true;
	}


	#pragma omp for
	for (int z = 0; z < NUMBER_OF_ZONES; z++) {
		const double gradVx = ZONES[z].gradVx;
		const double gradVy = ZONES[z].gradVy;
		const double gradDx = ZONES[z].gradDx;
		const double gradDy = ZONES[z].gradDy;

		const double D = DV[2*z];

	//	const double Dnoise = DV[2*z+1];
		for (int j = 0; j < ZONES[z].translocations; j++) {
			const double dt = ZONES[z].dt[j];
			const double dx = ZONES[z].dx[j];
			const double dy = ZONES[z].dy[j];
			const double  Dnoise = LOCALIZATION_ERROR*LOCALIZATION_ERROR/dt;

			result += - log(4.0*PI*(D + Dnoise)*dt) - ((dx-D*gradVx*dt)*(dx-D*gradVx*dt) + (dy-D*gradVy*dt)*(dy-D*gradVy*dt))/(4.0*(D+Dnoise)*dt);
		//	//fprintf(stderr, "%i\t %i\t %lf\t %lf\t %lf\t %lf\t %lf\t %le\n", z,j,dx,dy, D, gradVx,gradVy , result);
		//	////fprintf(stderr, "diffusion %lf\t diffusion noise %lf\t value result %lf\n", D,Dnoise,result );

			}

			if (ZONES[z].priorActive == true) {
				// Smoothing
			//	fprintf()
			//		////fprintf(stderr, "%lf\t %lf\t %lf\t %lf\n ", gradDx, gradDy, ZONES[z].areaX, ZONES[z].areaY);
			    result -= V_PRIOR*(gradVx*gradVx*ZONES[z].areaX + gradVy*gradVy*ZONES[z].areaY);
			    result -= D_PRIOR*(gradDx*gradDx*ZONES[z].areaX + gradDy*gradDy*ZONES[z].areaY);
				// result -= fabs(ZONES[z].potential);
				//////fprintf(stderr, "%lf\n", D_PRIOR*(gradDx*gradDx*ZONES[z].areaX + gradDy*gradDy*ZONES[z].areaY));
				//	////fprintf(stderr, "here\n");
				// Jeffreys
				  if (JEFFREYS_PRIOR == 1) {
					// result += 2.0*log(D*1.00) - 2.0*log(D*dt + LOCALIZATION_ERROR*LOCALIZATION_ERROR);
					    result += 2.0*log(D*1.00) - 2.0*log(D*ZONES[z].dtMean + LOCALIZATION_ERROR*LOCALIZATION_ERROR);
			   	}
		   }
	}

////fprintf(stderr,"in function %le\n", dt );
////fprintf(stderr,"in function %le\n", result );

	return -result;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/********* Gradients calculi ***************************************************/
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

double dvGradVx(double *DV, int i) {

	const double centre = ZONES[i].xCentre;
	const double potential = DV[2*i+1];
	double rCentre, lCentre;
	double rPotential, lPotential;
	int rNumber, lNumber;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rPotential = lPotential = 0.0;
	rNumber = lNumber = 0;


	for (int p = 0; p < ZONES[i].nLeftNeighbours; p++) {
		lCentre += (findZone(ZONES[i].leftNeighbours[p])->xCentre-centre);
		lPotential += DV[2*findZone(ZONES[i].leftNeighbours[p])->id+1];
		lNumber++;
	}
	for (int p = 0; p < ZONES[i].nRightNeighbours; p++) {
		rCentre += (findZone(ZONES[i].rightNeighbours[p])->xCentre-centre);
		rPotential += DV[2*findZone(ZONES[i].rightNeighbours[p])->id+1];
		rNumber++;
	}
////fprintf(stderr, "in grad %i\t %i\t %i\n", i, ZONES[i].nLeftNeighbours, ZONES[i].nRightNeighbours);
	// average left xCentres
	if (lNumber > 0) {
		lCentre /= lNumber*1.0;
		lPotential /= lNumber;
	}

	// average right xCentres
	if (rNumber > 0) {
		rCentre /= rNumber*1.0;
		rPotential /= rNumber;
	}

	lCentre += centre;
	rCentre += centre;

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lPotential;
	q[1] = potential;
	q[2] = rPotential;

	// for case where neighbours on both sides
	if (lNumber != 0 && rNumber != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (lNumber != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (rNumber != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
	//	//fprintf(stderr, "here error grad x\n");
		xForce = 0.0;
	}

	return xForce;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
double dvGradVy(double *DV, int i) {
	const double centre = ZONES[i].yCentre;
	const double potential = DV[2*i+1];
	double tCentre, bCentre;
	double tPotential, bPotential;
	int tNumber, bNumber;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	bCentre = tCentre = 0.0;
	bPotential = tPotential = 0.0;
	bNumber = tNumber = 0;

	for (int p = 0; p < ZONES[i].nTopNeighbours; p++) {
		tCentre += (findZone(ZONES[i].topNeighbours[p])->yCentre-centre);
		tPotential += DV[2*findZone(ZONES[i].topNeighbours[p])->id+1];
		tNumber++;
	}
	for (int p = 0; p < ZONES[i].nBottomNeighbours; p++) {
		bCentre += (findZone(ZONES[i].bottomNeighbours[p])->yCentre-centre);
		bPotential += DV[2*findZone(ZONES[i].bottomNeighbours[p])->id+1];
		bNumber++;
	}

	// average left xCentres
	if (tNumber > 0) {
		tCentre /= tNumber;
		tPotential /= tNumber;
	}

	// average right xCentres
	if (bNumber > 0) {
		bCentre /= bNumber;
		bPotential /= bNumber;
	}

	tCentre += centre;
	bCentre += centre;

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bPotential;
	q[1] = potential;
	q[2] = tPotential;

	// for case where neighbours on both sides
	if (tNumber != 0 && bNumber != 0) {
		vander(y,w,q,2);
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on bottom
	else if (bNumber != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on top
	else if (tNumber != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
double dvGradDx(double *DV, int i) {

	const double centre = ZONES[i].xCentre;
	const double diffusion = DV[2*i];
	double rCentre, lCentre;
	double rDiffusion, lDiffusion;
	int rNumber, lNumber;
	double xForce, x[3], w[3], q[3] ;

	// initialization
	rCentre = lCentre = 0.0;
	rDiffusion = lDiffusion = 0.0;
	rNumber = lNumber = 0;

	for (int p = 0; p < ZONES[i].nLeftNeighbours; p++) {
		lCentre += (findZone(ZONES[i].leftNeighbours[p])->xCentre-centre);
		lDiffusion += DV[2*findZone(ZONES[i].leftNeighbours[p])->id];
		lNumber++;
	}
	for (int p = 0; p < ZONES[i].nRightNeighbours; p++) {
		rCentre += (findZone(ZONES[i].rightNeighbours[p])->xCentre-centre);
		rDiffusion += DV[2*findZone(ZONES[i].rightNeighbours[p])->id];
		rNumber++;
	}

	// average left xCentres
	if (lNumber > 0) {
		lCentre /= lNumber;
		lDiffusion /= lNumber;
	}

	// average right xCentres
	if (rNumber > 0) {
		rCentre /= rNumber;
		rDiffusion /= rNumber;
	}

	lCentre += centre;
	rCentre += centre;

	// assign to arrays
	x[0] = lCentre;
	x[1] = centre;
	x[2] = rCentre;
	q[0] = lDiffusion;
	q[1] = diffusion;
	q[2] = rDiffusion;

	// for case where neighbours on both sides
	if (lNumber != 0 && rNumber != 0) {
		vander(x,w,q,2);
		// return derivative of polynomial
		xForce = - w[1] - 2.0*w[2]*x[1];
	}
	// for case where only neighbours on left side
	else if (lNumber != 0) {
		xForce = - (q[1]-q[0])/fabs(x[1]-x[0]);
	}
	// for case where only neighbours on right side
	else if (rNumber != 0) {
		xForce = - (q[2]-q[1])/fabs(x[2]-x[1]);
	}
	// erroneous case
	else {
		xForce = 0.0;
	}

	return xForce;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
double dvGradDy(double *DV, int i) {
	const double centre = ZONES[i].yCentre;
	const double diffusion = DV[2*i];
	double tCentre, bCentre;
	double tDiffusion, bDiffusion;
	int tNumber, bNumber;
	double yForce, y[3], w[3], q[3] ;

	// initialization
	bCentre = tCentre = 0.0;
	bDiffusion = tDiffusion = 0.0;
	bNumber = tNumber = 0;

	for (int p = 0; p < ZONES[i].nTopNeighbours; p++) {
		tCentre += (findZone(ZONES[i].topNeighbours[p])->yCentre-centre);
		tDiffusion += DV[2*findZone(ZONES[i].topNeighbours[p])->id];
		tNumber++;
	}
	for (int p = 0; p < ZONES[i].nBottomNeighbours; p++) {
		bCentre += (findZone(ZONES[i].bottomNeighbours[p])->yCentre-centre);
		bDiffusion += DV[2*findZone(ZONES[i].bottomNeighbours[p])->id];
		bNumber++;
	}

	// average left xCentres
	if (tNumber > 0) {
		tCentre /= tNumber;
		tDiffusion /= tNumber;
	}

	// average right xCentres
	if (bNumber > 0) {
		bCentre /= bNumber;
		bDiffusion /= bNumber;
	}

	tCentre += centre;
	bCentre += centre;

	// assign to arrays
	y[0] = bCentre;
	y[1] = centre;
	y[2] = tCentre;
	q[0] = bDiffusion;
	q[1] = diffusion;
	q[2] = tDiffusion;

	// for case where neighbours on both sides
	if (tNumber != 0 && bNumber != 0) {
		vander(y,w,q,2);
		yForce = - w[1] - 2.0*w[2]*y[1];
	}
	// for case where only neighbours on bottom
	else if (bNumber != 0) {
		yForce = - (q[1]-q[0])/fabs(y[1]-y[0]);
	}
	// for case where only neighbours on top
	else if (tNumber != 0) {
		yForce = - (q[2]-q[1])/fabs(y[2]-y[1]);
	}
	// erroneous case
	else {
		yForce = 0.0;
	}

	return yForce;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
