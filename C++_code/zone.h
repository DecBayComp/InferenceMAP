#ifndef ZONE_H_
#define ZONE_H_

#include <fstream>
#include <sys/types.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

class Zone;
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/***** CLASSES *****/
class Zone {
	public:
	Zone() {
		this->fileId = -1;
		this->id = -1;
		this->translocations = 0;
		this->xCentre = 0.0;
		this->yCentre = 0.0;
		this->nLeftNeighbours = 0;
		this->nRightNeighbours = 0;
		this->nTopNeighbours = 0;
		this->nBottomNeighbours = 0;
		this->dx = NULL;
		this->dy = NULL;
		this->dt = NULL;
		this->leftNeighbours = NULL;
		this->rightNeighbours = NULL;
		this->topNeighbours = NULL;
		this->bottomNeighbours = NULL;
		this->areaX = 0.0;
		this->areaY = 0.0;
		this->area = 0.0;
		this->areaConvhull = 0.0;
		this->diffusion = 0.0;
		this->potential = 0.0;
		this->xForce = 0.0;
		this->yForce = 0.0;
		this->gradVx = 0.0;
		this->gradVy = 0.0;
		this->gradDx = 0.0;
		this->gradDy = 0.0;
		this->xDrift = 0.0;
		this->yDrift = 0.0;
		this->dtMean = 0.0;
		this->dxMean = 0.0;
		this->dyMean = 0.0;
		this->priorActive = false;
	}
	/*******************************************************************************/
	/*******************************************************************************/
	/*******************************************************************************/
	void Print() {
		fprintf(stderr,"***** ZONE %i *****\n",fileId);
		fprintf(stderr,"%i Translocations\n",translocations);
		fprintf(stderr,"Centred at (%f,%f)\n",xCentre,yCentre);
		fprintf(stderr,"Left Neighbours (%i) : ",nLeftNeighbours);
		for (int n = 0; n < nLeftNeighbours; n++) { fprintf(stderr,"%i ",leftNeighbours[n]); } fprintf(stderr,"\n");
		fprintf(stderr,"Right Neighbours (%i) : ",nRightNeighbours);
		for (int n = 0; n < nRightNeighbours; n++) { fprintf(stderr,"%i ",rightNeighbours[n]); } fprintf(stderr,"\n");
		fprintf(stderr,"Top Neighbours (%i) : ",nTopNeighbours);
		for (int n = 0; n < nTopNeighbours; n++) { fprintf(stderr,"%i ",topNeighbours[n]); } fprintf(stderr,"\n");
		fprintf(stderr,"Bottom Neighbours (%i) : ",nBottomNeighbours);
		for (int n = 0; n < nBottomNeighbours; n++) { fprintf(stderr,"%i ",bottomNeighbours[n]); } fprintf(stderr,"\n");
		fprintf(stderr,"dx [um]\t\tdy [um]\t\tdt [s]\n");
		for (int q = 0; q < translocations; q++) {
			fprintf(stderr,"%lf\t%lf\t%lf\n",dx[q],dy[q],dt[q]);
		}
	}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
	// loaded parameters
	int fileId;
	int id;
	int translocations;
	double *dx;
	double *dy;
	double *dt;
	double xCentre;
	double yCentre;
	int nLeftNeighbours,nRightNeighbours,nTopNeighbours,nBottomNeighbours;
	int *leftNeighbours;
	int *rightNeighbours;
	int *topNeighbours;
	int *bottomNeighbours;
	double dxMean;
	double dyMean;
	double dtMean;

	// inferred parameters
	double diffusion;
	double xForce;
	double yForce;
	double gradVx;
	double gradVy;
	double gradDx;
	double gradDy;
	double areaX;
	double areaY;
	double area;
	double areaConvhull;
	double potential;
	double xDrift;
	double yDrift;
	bool priorActive;
};

#endif /* ZONE_H_ */
