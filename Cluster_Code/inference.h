#ifndef INFERENCE_H_
#define INFERENCE_H_

#include <math.h>

#include "inference.h"
#include "optimization.h"

// OpenGL Libraries
#ifdef __APPLE__
	#include <omp.h>
#endif

// D Inference
double dPosterior(double *D);
double dPosteriorSmoothingSelection(double *D);
double dPosteriorSmoothing(double *D);
double dGradDx(double *D, int i);
double dGradDy(double *D, int i);

double dPosteriorSmoothingRandomizedOptimization(double *D);
double dPosteriorSmoothingRandomizedOptimizationSelection(double *D);
double dGradDxRandomizedOptimization(double *D, int i);
double dGradDyRandomizedOptimization(double *D, int i);

// DF inference functions
double dfPosterior(double *coeff_f);
double dfPosteriorSmoothing(double *DFxFy);
double dfPosteriorSmoothingSelection(double *DFxFy);
double dfGradDx(double *DFxFy, int i);
double dfGradDy(double *DFxFy, int i);
double dfFxValue(double *V_ij, int i);
double dfFyValue(double *V_ij, int i);
double squareDifference(double *V_ij);

// DV  inference functions
double dvPosterior(double *DV);
double dvPosteriorSmoothing(double *DV);
double dvPosteriorSelection(double *DV);
double dvPosteriorSmoothingSelection(double *DV);
double dvGradVx(double *DV, int i);
double dvGradVy(double *DV, int i);
double dvGradDx(double *DV, int i);
double dvGradDy(double *DV, int i);

// DD inference function
double dDDPosterior(double *DD);


#endif /* INFERENCE_H_ */
