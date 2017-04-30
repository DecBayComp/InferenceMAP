/*
 * main.cpp
 *
 *  Created on: May 15, 2015
 *      Author: mohamed
 *  Modif JBM
 */




/***** EXTERNAL LIBRARIES *****/

#include <fstream>

#include <sys/types.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

/***** HEADER FILES *****/
#include "zone.h"
#include "inference.h"
#include "optimization.h"

/***** CONSTANTS *****/
#define PI 3.1415926535897932384626433832795028841971693993751
unsigned int microseconds=1000;

/***** FUNCTIONS *****/

/***** FUNCTION PROTOTYPES *****/
void parseCommandLine(int argc, char *argv[], char *inputFile, char *outputFile, double *localizationError, char *inferenceMode, double *dPrior, double *vPrior, int *jeffreysPrior);
void helpMessage();
void loadData(char *inputFile);
void saveData(char *outputFile);
void computeAreas();
void inferD(int zone);
void inferDF(int zone);
void inferDV();
void inferDD();

Zone* findZone(int id);
Zone* findZone2(int id);


/***** GLOBAL VARIABLES *****/
Zone   *ZONES = NULL;
char   INPUT_FILE[FILENAME_MAX];
char   OUTPUT_FILE[FILENAME_MAX];
double LOCALIZATION_ERROR;
char   INFERENCE_MODE[FILENAME_MAX];
double D_PRIOR, V_PRIOR;
FILE   *INPUT,*OUTPUT;
int    NUMBER_OF_ZONES;
int    CURRENT_ZONE;
double **HESSIAN;
int    DIMENSIONS;
int    JEFFREYS_PRIOR;


/*******************************************************************************/
/*******************************************************************************/
/***** FUNCTION DEFINITIONS *****/
/*******************************************************************************/
/*******************************************************************************/
void parseCommandLine(int argc, char *argv[], char *inputFile, char *outputFile, double *localizationError, char *inferenceMode, double *dPrior, double *vPrior, int *jeffreysPrior) {

	strcpy(inputFile, "NO_FILE");
	strcpy(inferenceMode, "NO_FILE");

	for (int k = 1; k < argc; k++) {
		if (*argv[k] == '-'){
			switch (*(argv[k] + 1))
			{
				case 'i':
					if (sscanf(argv[++k], "%s", inputFile) == 0) { exit(-1); }
				break;
				case 'o':
					if (sscanf(argv[++k], "%s", outputFile) == 0) { exit(-1); }
				break;
				case 'e':
					if (sscanf(argv[++k], "%lf", localizationError) == 0) { exit(-1); }
					break;
				case 'm':
					if (sscanf(argv[++k], "%s", inferenceMode) == 0) { exit(-1); }
					break;
				case 'd':
					if (sscanf(argv[++k], "%lf", dPrior) == 0) { *dPrior = 0.0; }
					break;
				case 'v':
					if (sscanf(argv[++k], "%lf", vPrior) == 0) { *vPrior = 0.0; }
					break;
				case 'j':
					if (sscanf(argv[++k], "%i", jeffreysPrior) == 0) { *jeffreysPrior = 0; }
					break;
				default:
					helpMessage();
					exit(-1);
			}
		}
	}
	if (!strcmp(inputFile, "NO_FILE")) {
		helpMessage();
		exit(-1);
	}
	if (!strcmp(inferenceMode, "NO_FILE")) {
		////fprintf(stderr, "No inference mode specified.\n");
		helpMessage();
		exit(-1);
	}

	*localizationError = *localizationError / 1000.0;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
void loadData(char *inputFile) {
	int id = -1, oldId = -1;
	int translocations;
	double xCentre;
	double yCentre;
	int nLeftNeighbours,nRightNeighbours,nTopNeighbours,nBottomNeighbours;
	double area,areaConvhull;

	char readString[100];

	int numberOfZones = 0;
	int ch;
	INPUT = fopen(inputFile,"r");

	// determine number of zones
	do {
		ch = fscanf(INPUT,"%s",&readString);
		if (!strcmp(readString,"ZONE:")) {
			numberOfZones++;
		}
	} while(ch != EOF);
//	//fprintf(stderr,"%i Zones in file.\n",numberOfZones);
	rewind(INPUT);

	// create zone array
	ZONES = new Zone[numberOfZones];

	int zoneCount = 0;

	// read data
	ch = 0;
	do {
		ch = fscanf(INPUT,"%s",&readString);
    //fprintf(stderr, "%s\n", readString);
		// read id
		if (!strcmp(readString,"ZONE:")) {
			fscanf(INPUT,"%i",&id);
			if (id != oldId) {
				ZONES[zoneCount].fileId = id;
				ZONES[zoneCount].id = zoneCount;
				zoneCount++;
			}
		}

		// read translocations
		if (!strcmp(readString,"NUMBER_OF_TRANSLOCATIONS:")) {
			fscanf(INPUT,"%i",&translocations);
			ZONES[zoneCount-1].translocations = translocations;
			ZONES[zoneCount-1].dt = new double[translocations];
			ZONES[zoneCount-1].dx = new double[translocations];
			ZONES[zoneCount-1].dy = new double[translocations];
//			//fprintf(stderr,"%i\ttranslocations = %i\n",zoneCount,ZONES[zoneCount-1].translocations);
		}

		// read centroid
		if (!strcmp(readString,"X-CENTRE:")) {
			fscanf(INPUT,"%lf",&xCentre);
			ZONES[zoneCount-1].xCentre = xCentre;
		}
		if (!strcmp(readString,"Y-CENTRE:")) {
			fscanf(INPUT,"%lf",&yCentre);
			ZONES[zoneCount-1].yCentre = yCentre;
		}
		if (!strcmp(readString,"AREA:")) {
			fscanf(INPUT,"%lf",&area);
			ZONES[zoneCount-1].area = area;
		}
		if (!strcmp(readString,"AREA_CONVHULL:")) {
			fscanf(INPUT,"%lf",&areaConvhull);
			ZONES[zoneCount-1].areaConvhull = areaConvhull;
			ZONES[zoneCount-1].areaX        = areaConvhull;
			ZONES[zoneCount-1].areaY        = areaConvhull;
		}
//usleep(microseconds);
		// read left neighbours
	//	if ( (!strcmp(readString,"NUMBER_OF_LEFT_NEIGHBOURS:"))||(!strcmp(readString,"MBER_OF_LEFT_NEIGHBOURS:")) ){
		//if (!strcmp(readString,"NUMBER_OF_LEFT_NEIGHBOURS:")) {
			if  (!strcmp(readString,"MBER_OF_LEFT_NEIGHBOURS:")) {
			fscanf(INPUT,"%i",&nLeftNeighbours);
			ZONES[zoneCount-1].nLeftNeighbours = nLeftNeighbours;
			ZONES[zoneCount-1].leftNeighbours = new int[nLeftNeighbours];
			//fprintf(stderr, "number left neighbors %i\n", ZONES[zoneCount-1].nLeftNeighbours );
		}
		if (!strcmp(readString,"LEFT_NEIGHBOURS:")) {
			const int N = ZONES[zoneCount-1].nLeftNeighbours;
			if (N > 0) {
				for (int n = 0; n < N; n++) {
					fscanf(INPUT,"%i",&ZONES[zoneCount-1].leftNeighbours[n]);
				}
			}
		}

		// read right neighbours
	//	if ( (!strcmp(readString,"NUMBER_OF_RIGHT_NEIGHBOURS:"))||(!strcmp(readString,"MBER_OF_RIGHT_NEIGHBOURS:")) ){
		//if (!strcmp(readString,"NUMBER_OF_RIGHT_NEIGHBOURS:")) {
		if (!strcmp(readString,"MBER_OF_RIGHT_NEIGHBOURS:")){
			fscanf(INPUT,"%i",&nRightNeighbours);
			ZONES[zoneCount-1].nRightNeighbours = nRightNeighbours;
			ZONES[zoneCount-1].rightNeighbours = new int[nRightNeighbours];
				//fprintf(stderr, "number left neighbors %i\n", ZONES[zoneCount-1].nRightNeighbours  );
		//	//fprintf(stderr, "number left neighbors %i\n", ZONES[zoneCount-1].nRightNeighbours );
		}
		if (!strcmp(readString,"RIGHT_NEIGHBOURS:")) {
			const int N = ZONES[zoneCount-1].nRightNeighbours;
			if (N > 0) {
				for (int n = 0; n < N; n++) {
					fscanf(INPUT,"%i",&ZONES[zoneCount-1].rightNeighbours[n]);
				}
			}
		}

		// read top neighbours
		//if ( (!strcmp(readString,"NUMBER_OF_TOP_NEIGHBOURS:"))||(!strcmp(readString,"MBER_OF_TOP_NEIGHBOURS:")) ){
		//if (!strcmp(readString,"NUMBER_OF_TOP_NEIGHBOURS:")) {
		if (!strcmp(readString,"MBER_OF_TOP_NEIGHBOURS:")){
			fscanf(INPUT,"%i",&nTopNeighbours);
			ZONES[zoneCount-1].nTopNeighbours = nTopNeighbours;
			ZONES[zoneCount-1].topNeighbours = new int[nTopNeighbours];
			//fprintf(stderr, "number left neighbors %i\n", ZONES[zoneCount-1].nTopNeighbours  );
		}
		if (!strcmp(readString,"TOP_NEIGHBOURS:")) {
			const int N = ZONES[zoneCount-1].nTopNeighbours;
			if (N > 0) {
				for (int n = 0; n < N; n++) {
					fscanf(INPUT,"%i",&ZONES[zoneCount-1].topNeighbours[n]);
				}
			}
		}
		// read bottom neighbours
	//	if ( (!strcmp(readString,"NUMBER_OF_BOTTOM_NEIGHBOURS:"))||(!strcmp(readString,"MBER_OF_BOTTOM_NEIGHBOURS:")) ){
			if ( !strcmp(readString,"MBER_OF_BOTTOM_NEIGHBOURS:")) {
			fscanf(INPUT,"%i",&nBottomNeighbours);
			ZONES[zoneCount-1].nBottomNeighbours = nBottomNeighbours;
			ZONES[zoneCount-1].bottomNeighbours = new int[nBottomNeighbours];
			//fprintf(stderr, "number left neighbors %i\n", ZONES[zoneCount-1].nBottomNeighbours );
		}
		if (!strcmp(readString,"BOTTOM_NEIGHBOURS:")) {
			const int N = ZONES[zoneCount-1].nBottomNeighbours;
			if (N > 0) {
				for (int n = 0; n < N; n++) {
					fscanf(INPUT,"%i",&ZONES[zoneCount-1].bottomNeighbours[n]);
				}
			}
		}

		// skip a line
		fscanf(INPUT, "%*[^\n]\n", NULL);

		// read dx dy dt
		for (int j = 0; j < ZONES[zoneCount-1].translocations; j++) {
			fscanf(INPUT,"%lf\t%lf\t%lf\t",&ZONES[zoneCount-1].dx[j],&ZONES[zoneCount-1].dy[j],&ZONES[zoneCount-1].dt[j]);
			////fprintf(stderr, "%lf\t%lf\t%lf\n",ZONES[zoneCount-1].dx[j],ZONES[zoneCount-1].dy[j],ZONES[zoneCount-1].dt[j]);
//			if (zoneCount-1 ==0){
	//			//fprintf(stderr, "%lf\t%lf\t%lf\n",ZONES[zoneCount-1].dx[j],ZONES[zoneCount-1].dy[j],ZONES[zoneCount-1].dt[j]);
		//	}
		}

		oldId = id;

	} while(ch != EOF);

	NUMBER_OF_ZONES = zoneCount;
//	for (int z = 0; z < numberOfZones; z++) { ZONES[z].Print(); }

/*int z = 0;
 ZONES[z].Print();
z++;
 ZONES[z].Print();*/

	fclose(INPUT);
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
void saveData(char *outputFile) {

	OUTPUT = fopen(outputFile,"wb");

	// File Information
	fprintf(OUTPUT,"Input File: %s\n",INPUT_FILE);
	fprintf(OUTPUT,"Inference Mode: %s\n",INFERENCE_MODE);
	fprintf(OUTPUT,"Localization Error: %lf [nm]\n",LOCALIZATION_ERROR*1000.0);
	fprintf(OUTPUT,"Number of Zones: %i\n",NUMBER_OF_ZONES);
	fprintf(OUTPUT,"Diffusion Prior: %lf\n",D_PRIOR);
	fprintf(OUTPUT,"Potential Prior: %lf\n",V_PRIOR);
	fprintf(OUTPUT,"Jeffreys' Prior: %i\n\n",JEFFREYS_PRIOR);

	fprintf(OUTPUT,"Zone ID\tTranslocations\tx-Centre\ty-Centre\tDiffusion\tx-Force\ty-Force\tPotential\n");
	for (int z = 0; z < NUMBER_OF_ZONES; z++) {
		fprintf(OUTPUT,"%i\t%i\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",ZONES[z].fileId,ZONES[z].translocations,ZONES[z].xCentre,ZONES[z].yCentre,ZONES[z].diffusion,ZONES[z].xForce,ZONES[z].yForce,ZONES[z].potential);
	}

	fclose(OUTPUT);
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
void helpMessage() {
	//fprintf(stderr, "\nInvalid Command.\n");
	//fprintf(stderr, "Input format:\n");
	//fprintf(stderr, "./InferenceMAP_Terminal -i INPUT_FILE -o OUTPUT_FILE -e LOCALIZATION_ERROR_IN_NM -m INFERENCE_MODE -d DIFFUSION_PRIOR -v POTENTIAL_PRIOR -j JEFFREYS_PRIOR\n\n");
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
void computeAreas() {
	for (int z = 0; z < NUMBER_OF_ZONES; z++) {
		ZONES[z].areaX = ZONES[z].areaConvhull;
		ZONES[z].areaY = ZONES[z].areaConvhull;
		for (int l = 0; l < ZONES[z].nLeftNeighbours; l++) { ZONES[z].areaX += findZone(ZONES[z].leftNeighbours[l])->areaConvhull; }
		for (int l = 0; l < ZONES[z].nRightNeighbours; l++) { ZONES[z].areaX += findZone(ZONES[z].rightNeighbours[l])->areaConvhull; }
		for (int l = 0; l < ZONES[z].nTopNeighbours; l++) { ZONES[z].areaY += findZone(ZONES[z].topNeighbours[l])->areaConvhull; }
		for (int l = 0; l < ZONES[z].nBottomNeighbours; l++) { ZONES[z].areaY += findZone(ZONES[z].bottomNeighbours[l])->areaConvhull; }
	}
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
void inferD(int zone) {
	CURRENT_ZONE = zone;

	// calculate relevant means
	double dxMean2 = 0.0;
	double dyMean2 = 0.0;
	double dtMean = 0.0;

	for (int j = 0; j < ZONES[zone].translocations; j++) {
		dxMean2 += pow(ZONES[CURRENT_ZONE].dx[j],2);
		dyMean2 += pow(ZONES[CURRENT_ZONE].dy[j],2);
		dtMean += ZONES[CURRENT_ZONE].dt[j];
	}
	dxMean2 /= ZONES[CURRENT_ZONE].translocations;
	dyMean2 /= ZONES[CURRENT_ZONE].translocations;
	dtMean /= ZONES[CURRENT_ZONE].translocations;

	// prepare optimization data structures
	DIMENSIONS = 1;

	HESSIAN = new double*[DIMENSIONS];
	for(int i = 0; i < DIMENSIONS; i++) {
		HESSIAN[i] = new double[DIMENSIONS];
	}
	for (int d = 0; d < DIMENSIONS; d++) {
		for (int e = 0; e < DIMENSIONS; e++) {
			HESSIAN[d][e] = 0.0;
		}
	}

	int iterations[1];
	double fret[1];

	double *optimizationArray = new double[DIMENSIONS];
	for (int i = 0; i < DIMENSIONS; i++) { optimizationArray[i] = 0.0; }

	const double D_eff_x = dxMean2/(2.*dtMean);
	const double D_eff_y = dyMean2/(2.*dtMean);

	optimizationArray[0] = 0.5*(D_eff_x + D_eff_y);
//	//fprintf(stderr, "initial values %lf\t,%lf\t %lf\t %lf\n ",optimizationArray[0], dxMean2, dyMean2,dtMean );

	dfpmin(optimizationArray, DIMENSIONS, GTOL, iterations, fret, dPosterior, (dfunc));

	ZONES[CURRENT_ZONE].diffusion = optimizationArray[0];

//	//fprintf(stderr,"Zone %i diffusion = %lf\t(%i Iterations)\n",CURRENT_ZONE,ZONES[CURRENT_ZONE].diffusion,iterations[0]);

	for (int d = 0; d < DIMENSIONS; d++) {
		delete [] HESSIAN[d];
	}
	delete [] HESSIAN;
	HESSIAN = NULL;
	delete [] optimizationArray;
	optimizationArray = NULL;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
void inferDF(int zone) {
	CURRENT_ZONE = zone;

	// calculate relevant means
	ZONES[CURRENT_ZONE].dxMean = 0.0;
	ZONES[CURRENT_ZONE].dyMean = 0.0;
	ZONES[CURRENT_ZONE].dtMean = 0.0;

	for (int j = 0; j < ZONES[zone].translocations; j++) {
		ZONES[CURRENT_ZONE].dxMean += pow(ZONES[CURRENT_ZONE].dx[j],2);
		ZONES[CURRENT_ZONE].dyMean += pow(ZONES[CURRENT_ZONE].dy[j],2);
		ZONES[CURRENT_ZONE].dtMean += ZONES[CURRENT_ZONE].dt[j];
	}
	ZONES[CURRENT_ZONE].dxMean /= ZONES[CURRENT_ZONE].translocations;
	ZONES[CURRENT_ZONE].dyMean /= ZONES[CURRENT_ZONE].translocations;
	ZONES[CURRENT_ZONE].dtMean /= ZONES[CURRENT_ZONE].translocations;

	// prepare optimization data structures
	DIMENSIONS = 3;

	HESSIAN = new double*[DIMENSIONS];
	for(int i = 0; i < DIMENSIONS; i++) {
		HESSIAN[i] = new double[DIMENSIONS];
	}
	for (int d = 0; d < DIMENSIONS; d++) {
		for (int e = 0; e < DIMENSIONS; e++) {
			HESSIAN[d][e] = 0.0;
		}
	}

	int iterations[1];
	double fret[1];

	double *optimizationArray = new double[DIMENSIONS];
	for (int i = 0; i < DIMENSIONS; i++) { optimizationArray[i] = 0.0; }

	const double D_eff_x = ZONES[CURRENT_ZONE].dxMean/(2.*ZONES[CURRENT_ZONE].dtMean);
	const double D_eff_y = ZONES[CURRENT_ZONE].dyMean/(2.*ZONES[CURRENT_ZONE].dtMean);

	optimizationArray[0] = 0.0;
	optimizationArray[1] = 0.0;
	optimizationArray[2] = 0.5*(D_eff_x + D_eff_y);

	dfpmin(optimizationArray, DIMENSIONS, GTOL, iterations, fret, dfPosterior, (dfunc));

	ZONES[CURRENT_ZONE].xForce = optimizationArray[0];
	ZONES[CURRENT_ZONE].yForce = optimizationArray[1];
	ZONES[CURRENT_ZONE].diffusion = optimizationArray[2];

//	//fprintf(stderr,"Zone %i : [%lf %lf %lf]\t(%i Iterations)\n",CURRENT_ZONE,ZONES[CURRENT_ZONE].diffusion,ZONES[CURRENT_ZONE].xForce,ZONES[CURRENT_ZONE].yForce,iterations[0]);

	for (int d = 0; d < DIMENSIONS; d++) {
		delete [] HESSIAN[d];
	}
	delete [] HESSIAN;
	HESSIAN = NULL;
	delete [] optimizationArray;
	optimizationArray = NULL;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
void inferDV() {

	int maxCount = -1;
	for (int z = 0; z < NUMBER_OF_ZONES; z++) {
		// calculate relevant means
		ZONES[z].dxMean = 0.0;
		ZONES[z].dyMean = 0.0;
		ZONES[z].dtMean = 0.0;
		for (int j = 0; j < ZONES[z].translocations; j++) {
			ZONES[z].dxMean += pow(ZONES[z].dx[j],2);
			ZONES[z].dyMean += pow(ZONES[z].dy[j],2);
			ZONES[z].dtMean += ZONES[z].dt[j];
			////fprintf(stderr, "%i\t %i\t %le\n",z, j, ZONES[z].dx[j]);
		}
		ZONES[z].dxMean /= ZONES[z].translocations;
		ZONES[z].dyMean /= ZONES[z].translocations;
		ZONES[z].dtMean /= ZONES[z].translocations;

		// find zone with maximum translocations
		if (maxCount < ZONES[z].translocations) { maxCount = ZONES[z].translocations; }
	}

	// prepare optimization data structures
	DIMENSIONS = 2*NUMBER_OF_ZONES;

	HESSIAN = new double*[DIMENSIONS];
	for(int i = 0; i < DIMENSIONS; i++) {
		HESSIAN[i] = new double[DIMENSIONS];
	}
	for (int d = 0; d < DIMENSIONS; d++) {
		for (int e = 0; e < DIMENSIONS; e++) {
			HESSIAN[d][e] = 0.0;
		}
	}

	int iterations[1];
	double fret[1];

	double *optimizationArray = new double[DIMENSIONS];
	for (int i = 0; i < DIMENSIONS; i++) { optimizationArray[i] = 0.0; }

	// initialize variables
	for (int z = 0; z < NUMBER_OF_ZONES; z++) {
		const double D_eff_x = ZONES[z].dxMean/(2.*ZONES[z].dtMean);
		const double D_eff_y = ZONES[z].dyMean/(2.*ZONES[z].dtMean);

		optimizationArray[2*z] = 0.5*(D_eff_x + D_eff_y);
		optimizationArray[2*z+1] = -1.0*log( (double)  ((double) ZONES[z].translocations *1.0)/((double) maxCount*1.0) );
	//	//fprintf(stderr, "initi %le\n", optimizationArray[2*z+1] );
	}
////fprintf(stderr, "%le\n",GTOL);
	dfpmin(optimizationArray, DIMENSIONS, GTOL, iterations, fret, dvPosterior, (dfunc));

	for (int z = 0; z < NUMBER_OF_ZONES; z++) {
		ZONES[z].diffusion = optimizationArray[2*z];
		ZONES[z].potential = optimizationArray[2*z+1];
	//	//fprintf(stderr,"Zone %i : [%lf %lf]\t(%i Iterations)\n",z,ZONES[z].diffusion,ZONES[z].potential,iterations[0]);
	}

	for (int d = 0; d < DIMENSIONS; d++) { delete [] HESSIAN[d]; }
	delete [] HESSIAN;
	HESSIAN = NULL;
	delete [] optimizationArray;
	optimizationArray = NULL;
}

Zone* findZone(int id) {
	for (int z = 0; z < NUMBER_OF_ZONES; z++) {
		if (ZONES[z].fileId == id) {
			return &ZONES[z];
			break;
		}
	}
	return NULL;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
void inferDD() {

	int maxCount = -1;
	for (int z = 0; z < NUMBER_OF_ZONES; z++) {
		// calculate relevant means
		ZONES[z].dxMean = 0.0;
		ZONES[z].dyMean = 0.0;
		ZONES[z].dtMean = 0.0;
		for (int j = 0; j < ZONES[z].translocations; j++) {
			ZONES[z].dxMean += pow(ZONES[z].dx[j],2);
			ZONES[z].dyMean += pow(ZONES[z].dy[j],2);
			ZONES[z].dtMean += ZONES[z].dt[j];
		}
		ZONES[z].dxMean /= ZONES[z].translocations;
		ZONES[z].dyMean /= ZONES[z].translocations;
		ZONES[z].dtMean /= ZONES[z].translocations;

		// find zone with maximum translocations
		if (maxCount < ZONES[z].translocations) { maxCount = ZONES[z].translocations; }
	}

	// prepare optimization data structures
	DIMENSIONS = NUMBER_OF_ZONES;

	HESSIAN = new double*[DIMENSIONS];
	for(int i = 0; i < DIMENSIONS; i++) {
		HESSIAN[i] = new double[DIMENSIONS];
	}
	for (int d = 0; d < DIMENSIONS; d++) {
		for (int e = 0; e < DIMENSIONS; e++) {
			HESSIAN[d][e] = 0.0;
		}
	}

	int iterations[1];
	double fret[1];

	double *optimizationArray = new double[DIMENSIONS];
	for (int i = 0; i < DIMENSIONS; i++) { optimizationArray[i] = 0.0; }

	// initialize variables
	for (int z = 0; z < NUMBER_OF_ZONES; z++) {
		const double D_eff_x = ZONES[z].dxMean/(2.*ZONES[z].dtMean);
		const double D_eff_y = ZONES[z].dyMean/(2.*ZONES[z].dtMean);

		optimizationArray[z] = 0.5*(D_eff_x + D_eff_y);
	}

	dfpmin(optimizationArray, DIMENSIONS, GTOL, iterations, fret, dDDPosterior, (dfunc));

	for (int z = 0; z < NUMBER_OF_ZONES; z++) {
		ZONES[z].diffusion = optimizationArray[z];
		}

	for (int d = 0; d < DIMENSIONS; d++) { delete [] HESSIAN[d]; }
	delete [] HESSIAN;
	HESSIAN = NULL;
	delete [] optimizationArray;
	optimizationArray = NULL;
}

Zone* findZone2(int id) {
	for (int z = 0; z < NUMBER_OF_ZONES; z++) {
		if (ZONES[z].fileId == id) {
			return &ZONES[z];
			break;
		}
	}
	return NULL;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

int main(int argc, char* argv[]) {

	// parse command line
	parseCommandLine(argc,argv,INPUT_FILE,OUTPUT_FILE,&LOCALIZATION_ERROR,INFERENCE_MODE,&D_PRIOR,&V_PRIOR,&JEFFREYS_PRIOR);

	// load data from input file
	loadData(INPUT_FILE);

	// perform inference
	if (strcmp(INFERENCE_MODE,"D")==0) { for (int z = 0; z < NUMBER_OF_ZONES; z++) { inferD(z); } }
	else if (strcmp(INFERENCE_MODE,"DF")==0) { for (int z = 0; z < NUMBER_OF_ZONES; z++) { inferDF(z); } }
	else if (strcmp(INFERENCE_MODE,"DD")==0) {  inferDD();  }
	else if (strcmp(INFERENCE_MODE,"DV")==0) {  inferDV();  }
	else { fprintf(stderr,"Invalid inference mode\n"); exit(-1); }

	// save data to output file
	saveData(OUTPUT_FILE);

}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
