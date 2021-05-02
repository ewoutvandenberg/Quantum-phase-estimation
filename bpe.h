/*
# Copyright 2021 IBM Inc. All rights reserved
# SPDX-License-Identifier: Apache2.0

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This file is part of the code to reproduce the results in the paper:
# E. van den Berg, "Efficient Bayesian phase estimation using mixed priors"
# arXiv:2007.11629.
*/

#ifndef __BPE_H__
#define __BPE_H__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef enum { BPE_MODE_MIXED, BPE_MODE_NORMAL_ONLY, BPE_MODE_FOURIER_ONLY } BPEMode;

typedef struct
{  double *gridpoints; /* The segment size must be less than pi. */
   double *mass;
   int     n;
   int     capacity;
} TDistribution;

typedef struct
{  double *coefCos;
   double *coefSin;
   double *buffer;
   int     kMax;
   int     kCapacity;
} TFourier;

typedef struct
{  double    *phi;
   double    *w;
   int        nEig;
} BPEPhases;

typedef struct
{  TFourier **fourier;
   double    *mu;
   double    *sigma;
   double    *weights;
   int        nEig;
} BPEDensity;


typedef struct __BPEInitDensity
{  void (*fptrInitDensity)(struct __BPEInitDensity *self, BPEDensity *);
   void (*fptrDelete)(struct __BPEInitDensity *self);
} BPEInitDensity;


typedef struct __BPEConfig
{  int       nIter;         /* Number of iterations                                         */
   int       nEig;          /* Number of eigenvalues in the model                           */
   BPEMode   mode;          /* Representation of the distributions                          */
   double    epsilon;       /* Threshold for determining critical sigma in mixed mode (1e-4)*/
   int       kMax;          /* Maximum number of coefficients in the Fourier representation */

   int       displaySolve;  /* Display the solver statistics or not                         */
   int       updateDouble;  /* Iteration at which to start exponential scheme for updates   */
   int       curvilinear;   /* Use of curvilinear line search in the evaluation             */

   BPEInitDensity *pluginInitDensity;
   void  (*fptrInitWeights)(struct __BPEConfig *config, double *);
} BPEConfig;


typedef struct __BPEMeasurementScheme
{  int       *k;
   double    *beta;
   double    *betaHat;
   int       *measurements;
   int        nRounds;
   void      *data;          /* Additional data */
   TFourier  *probability;
   BPEPhases *phases;        /* Ground truth */ 

   /* Commonly used parameters */
   int        kMax;          /* Maximum number of unitary applications */
   int        kErr;          /* Parameter for decoherence              */
   double     noise1;
   double     noise2;
   int       *state;   

   void  (*fptrSetK       )(struct __BPEMeasurementScheme *, BPEDensity *, int);
   void  (*fptrSetBeta    )(struct __BPEMeasurementScheme *, BPEDensity *, int);
   void  (*fptrMeasure    )(struct __BPEMeasurementScheme *, BPEDensity *, int);
   void  (*fptrFreeData)(void *);
} BPEMeasurementScheme;



/* ----------------------------------------------------- */
/* In this implementation we assume the weight vector to */
/* small (up to around, say, 10 elements). This allows   */
/* us to make the following implementation choices:      */
/*                                                       */
/* 1. We can copy the final result back into the input   */
/*    vector, rather than work with additional buffers   */
/*    that are logically swapped at the end, which would */
/*    result in an updated input pointer.                */
/*                                                       */
/* 2. Projection onto the simplex requires sorting. The  */
/*    small number of elements to sort allows us to use  */
/*    insertion sort.                                    */
/* ----------------------------------------------------- */

typedef struct
{  double *C;      /* Matrix C                      */
   double *x;      /* Current solution   (size = n) */
   double *v;      /* C*x                (size = m) */
   double *g;      /* C'*v               (size = n) */
   double *xbar;   /* x - alpha * g      (size = n) */
   double *buf;    /* Buffer for projection         */
   int    m;       /* Number of rows in C           */
   int    n;       /* Number of columns in C        */
   int    mMax;    /* Maximum number of rows        */
   int    xOwner;  /* Do we own the pointer to x?   */
} TSolveWeights;



/* ------------------------------------------------------------------------ */
/* Discretized distribution                                                 */
/* ------------------------------------------------------------------------ */

extern void initDistribution(TDistribution *distribution);
extern void freeDistribution(TDistribution *distribution);
extern int  resizeDistribution(TDistribution *distribution, int size);

/* Setup a linear grid */
extern int  linearGrid(TDistribution *distribution, int n);

/* Priors */
extern void uniformPrior(TDistribution *distribution);
extern int  uniformPriorAround(TDistribution *distribution, double phi, double delta, double gridspacing);

/* Normalize the probability mass to one */
extern void normalizeDistribution(TDistribution *distribution);

/* Integrate the probability over given segment */
/* x and y must be between [-2*pi,4*pi)         */
extern double integrateProbability(TDistribution *distribution, double x, double y);

/* Get mean and error (modulo 2*PI) */
extern double getMean        (TDistribution *distribution, double *resultMSE);
extern double getErrorAbs    (TDistribution *distribution, double ref);
extern double getErrorSquared(TDistribution *distribution, double ref);

/* Update the model */
extern void   updateNoiseless(TDistribution *distribution,
                              double phi, double s, double theta,
                              int measurement);

/* Refine the distribution */
extern int    refineDistribution(TDistribution *src, TDistribution *dst, double pMerge, double pSplit);

/* Sample from distributions */
extern double sampleUniform(void);
extern void   sampleAngles(double *ptr, int n);
extern double sampleDistribution(TDistribution *distribution);

/* Measure a circuit */
extern int measureNoiseless(double phi, double s, double theta);

/* Determine the error in angle */
extern double angleDiff(double alpha, double beta);

/* Miscellaneous */
extern void outputDistribution(FILE *fp, TDistribution *distribution);
extern void ensureDirectory(const char *str);



/* --------------------------------------------------------------------- */
/* Multiple eigenvalues                                                  */
/* --------------------------------------------------------------------- */
extern void sampleRounds(int nEig, const double *weights, const double *phi,
                         int nRounds, const int *k, const double *beta,
                         double *weightsBuffer, int *measurements);

extern void sampleFourierRounds(int nEig, const double *weights, const double *phi,
                                int nRounds, double *weightsBuffer, int *measurements);


/* --------------------------------------------------------------------- */
/* Fourier representation                                                */
/* --------------------------------------------------------------------- */
extern TFourier *newFourier        (void);
extern void      deleteFourier     (TFourier *dist);
extern void      freeFourier       (TFourier *dist);
extern void      initFourier       (TFourier *dist);
extern void      initFourierZero   (TFourier *dist);
extern void      initFourierOne    (TFourier *dist);
extern void      initFourierUniform(TFourier *dist);
extern void      initFourierCosine (TFourier *dist, int nCoef, int *k, double *beta, double *weights);
extern void      initFourierSine   (TFourier *dist, int nCoef, int *k, double *beta, double *weights);
extern void      initFourierProb   (TFourier *dist, int nCoef, int *k, double *beta, int *measurements);
extern void      initFourierNormal (TFourier *dist, int kMax, double mu, double sigma);
extern void      setFourierCapacity(TFourier *dist, int kCapacity);
extern void      resizeFourier     (TFourier *dist, int kMax, int clearMode);
extern void      clearFourier      (TFourier *dist); /* zero out coefficients */
extern void      swapFourier       (TFourier *distA, TFourier *distB);
extern void      copyFourier       (TFourier *distA, TFourier *distB); /* distA <- distB */
extern void      addFourier        (TFourier *distA, TFourier *distB); /* distA <- distA + distB */
extern void      subtractFourier   (TFourier *distA, TFourier *distB); /* distA <- distA - distB */
extern void      multiplyFourier   (TFourier *distA, TFourier *distB); /* distA <- distA * distB */
extern void      scaleFourier      (TFourier *dist, double scale);     /* distA <- distA * scale */
extern void      normalizeFourier  (TFourier *dist);
extern void      foldFourier       (TFourier *dist); /* Map negative indices to positive ones */
extern void      truncateFourier   (TFourier *dist, int kMax);
extern double    evaluateFourier   (TFourier *dist, double phi);
extern double    getFourierMean    (TFourier *dist);
extern double    getFourierVar     (TFourier *dist);
extern double    getFourierStd     (TFourier *dist);
extern void      printFourierCoef  (TFourier *dist);

extern double    fourierNormalError(int kMax, double sigma);
extern double    fourierCriticalSigma(int kMax, double epsilon);


/* ------------------------------------------------------------------------ */
/* Weight solves                                                            */
/* ------------------------------------------------------------------------ */
void initSolveWeights     (TSolveWeights *data, int mMax, int n, double *x);
void freeSolveWeights     (TSolveWeights *data);
void solveWeights         (TSolveWeights *data, int m, int verbose);
void solveWeightsBacktrack(TSolveWeights *data, int m, int verbose);
void solveWeightsProject  (double *vec, double *buf, int n);


/* ------------------------------------------------------------------------ */
/* Phases                                                                   */
/* ------------------------------------------------------------------------ */
void newPhases(BPEPhases **phases, int nEig);
void deletePhases(BPEPhases *phases);


/* ------------------------------------------------------------------------ */
/* Density                                                                  */
/* ------------------------------------------------------------------------ */
void initDensity(BPEDensity *density, BPEMode mode, int nEig);
void freeDensity(BPEDensity *density);


/* ------------------------------------------------------------------------ */
/* Measurement scheme                                                       */
/* ------------------------------------------------------------------------ */
void newMeasurementScheme   (BPEMeasurementScheme **scheme, int capacity);
void deleteMeasurementScheme(BPEMeasurementScheme  *scheme);
void resizeMeasurementScheme(BPEMeasurementScheme  *scheme, int capacity);


/* ------------------------------------------------------------------------ */
/* Bayesian phase estimation algorithm                                      */
/* ------------------------------------------------------------------------ */
void BayesianPhaseEstimation(BPEMeasurementScheme *scheme, BPEConfig *config);


/* ------------------------------------------------------------------------ */
/* Miscellaneous                                                            */
/* ------------------------------------------------------------------------ */

/* Malloc with exit on fail */
void *exitMalloc(size_t size);
#define nullFree(ptr) { if ((ptr) != NULL) { free(ptr); ptr = NULL; } }


/* Replace random numbers with file */
#ifdef DEBUG_RANDOM
extern void      initRandomBuffer(void);
extern void      freeRandomBuffer(void);
#endif

/* Debug settings */
#ifdef VALIDATE_MEAN
extern FILE *fp_validate_mean;
#endif

#endif
